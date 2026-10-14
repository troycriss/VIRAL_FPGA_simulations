//! Floating-point mean-field model: spins are unit vectors (x, y, z), and the
//! coupling J is either rank 1 (J_ij = a_i a_j, what VIRAL_n_spin does) or a
//! general dense matrix (what a bigger FPGA design would do).
//!
//! Per layer: h_i = sum_{j != i} J_ij z_j, rotate spin i about z by c * h_i,
//! then every spin about x by phi.

use crate::{cordic, fft, par, rng::Rng};

pub type Spin = [f64; 3];

/// Dense symmetric couplings J_ij = scale * q_ij, stored as 8-bit integers
/// (the precision an FPGA design would store): `q` is the full N x N matrix
/// with zero diagonal, and `tiles` the upper triangle again in TILE x TILE
/// tiles, each contiguous, for fast products.
pub struct DenseJ {
    pub n: usize,
    pub q: Vec<i8>,
    pub scale: f64,
    tiles: Vec<i8>,
    /// (block row, block column) of each tile in `tiles`
    tile_pos: Vec<(usize, usize)>,
}

impl DenseJ {
    /// Sherrington-Kirkpatrick couplings with variance 1/N: +-1/sqrt(N), or
    /// Gaussian rounded to 8 bits (steps of 4/127 standard deviations).
    pub fn sk(n: usize, gaussian: bool, rng: &mut Rng) -> DenseJ {
        let step = if gaussian { 4.0 / 127.0 } else { 1.0 };
        let mut q = vec![0i8; n * n];
        for i in 0..n {
            for v in &mut q[i * n + i + 1..(i + 1) * n] {
                *v = if gaussian {
                    (rng.normal() / step).round().clamp(-127.0, 127.0) as i8
                } else if rng.coin() {
                    1
                } else {
                    -1
                };
            }
        }
        // mirror the upper triangle, in blocks to stay in cache
        const B: usize = 64;
        for bi in (0..n).step_by(B) {
            for bj in (bi..n).step_by(B) {
                for i in bi..(bi + B).min(n) {
                    for j in bj.max(i + 1)..(bj + B).min(n) {
                        q[j * n + i] = q[i * n + j];
                    }
                }
            }
        }
        DenseJ::new(n, q, step / (n as f64).sqrt())
    }

    /// From a full symmetric matrix with zero diagonal.
    pub fn new(n: usize, q: Vec<i8>, scale: f64) -> DenseJ {
        let nb = n.div_ceil(TILE);
        let tile_pos: Vec<(usize, usize)> = (0..nb).flat_map(|bi| (bi..nb).map(move |bj| (bi, bj))).collect();
        // tiles keep only j > i, so a diagonal tile is used like any other
        let mut tiles = vec![0i8; tile_pos.len() * TILE * TILE];
        for (t, &(bi, bj)) in tile_pos.iter().enumerate() {
            let tile = &mut tiles[t * TILE * TILE..(t + 1) * TILE * TILE];
            for r in 0..TILE.min(n - bi * TILE) {
                let i = bi * TILE + r;
                for c in 0..TILE.min(n - bj * TILE) {
                    let j = bj * TILE + c;
                    if j > i {
                        tile[r * TILE + c] = q[i * n + j];
                    }
                }
            }
        }
        DenseJ { n, q, scale, tiles, tile_pos }
    }

    /// h = q z for z in Q1.15, exactly, reading the upper triangle once:
    /// tile (I, J) adds q_IJ z_J to h_I and q_IJ^T z_I to h_J.
    pub fn matvec(&self, z: &[i16], threads: usize) -> Vec<i64> {
        let n = self.n;
        let len = n.div_ceil(TILE) * TILE;
        let mut zp = vec![0i16; len];
        zp[..n].copy_from_slice(z);
        // a thread start costs about as much as 1-2 million couplings
        let threads = par::threads_for(n * n / 100, threads);
        let mut parts =
            par::for_items(self.tile_pos.len(), threads, || vec![0i64; len], |t, h| self.tile(t, &zp, h));
        let mut h = parts.swap_remove(0);
        for p in &parts {
            for (a, b) in h.iter_mut().zip(p) {
                *a += b;
            }
        }
        h.truncate(n);
        h
    }

    /// Adds tile t into h, both ways.
    fn tile(&self, t: usize, z: &[i16], h: &mut [i64]) {
        let (bi, bj) = self.tile_pos[t];
        let tile = &self.tiles[t * TILE * TILE..(t + 1) * TILE * TILE];
        let zc = &z[bj * TILE..(bj + 1) * TILE];
        let mut acc = [0i32; TILE];
        // four rows per pass, so acc is loaded and stored once per four rows
        for r in (0..TILE).step_by(4) {
            let i = bi * TILE + r;
            let rows: [&[i8]; 4] = std::array::from_fn(|k| &tile[(r + k) * TILE..(r + k + 1) * TILE]);
            let zr: [i16; 4] = std::array::from_fn(|k| z[i + k]);
            let d = dot_axpy4(&rows, zc, &mut acc, &zr);
            for k in 0..4 {
                h[i + k] += d[k] as i64;
            }
        }
        for (hc, &a) in h[bj * TILE..(bj + 1) * TILE].iter_mut().zip(acc.iter()) {
            *hc += a as i64;
        }
    }

    /// Exact integer fields sum_j q_ij s_j for spins s = +-1.
    pub fn fields_int(&self, s: &[i8], threads: usize) -> Vec<i32> {
        let s16: Vec<i16> = s.iter().map(|&v| v as i16).collect();
        self.matvec(&s16, threads).into_iter().map(|v| v as i32).collect()
    }
}

/// Tile size for DenseJ::matvec: a sum of up to 256 products q * z
/// (|q| <= 127, |z| <= 2^15) fits in an i32.
const TILE: usize = 256;

/// Returns sum_k q_k z_k and adds q_k * zr to acc_k, in one pass over q.
/// Integer lanes, so the compiler vectorizes both (widening multiply-adds).
#[cfg(not(target_arch = "aarch64"))]
#[inline]
fn dot_axpy(q: &[i8], z: &[i16], acc: &mut [i32], zr: i32) -> i32 {
    const L: usize = 16;
    let m = q.len();
    let split = m - m % L;
    let mut d = [0i32; L];
    for ((qc, zc), ac) in q[..split]
        .chunks_exact(L)
        .zip(z[..split].chunks_exact(L))
        .zip(acc[..split].chunks_exact_mut(L))
    {
        let qc: &[i8; L] = qc.try_into().unwrap();
        let zc: &[i16; L] = zc.try_into().unwrap();
        let ac: &mut [i32; L] = ac.try_into().unwrap();
        for k in 0..L {
            let qk = qc[k] as i32;
            d[k] += qk * zc[k] as i32;
            ac[k] += qk * zr;
        }
    }
    let mut sum: i32 = d.iter().sum();
    for k in split..m {
        let qk = q[k] as i32;
        sum += qk * z[k] as i32;
        acc[k] += qk * zr;
    }
    sum
}

/// dot_axpy for four rows of q at once, so acc is loaded and stored once per
/// four rows. NEON on ARM (Apple silicon), where the compiler doesn't
/// vectorize this well by itself.
#[cfg(target_arch = "aarch64")]
fn dot_axpy4(q: &[&[i8]; 4], z: &[i16], acc: &mut [i32], zr: &[i16; 4]) -> [i32; 4] {
    use std::arch::aarch64::*;
    let m = z.len();
    assert!(acc.len() == m && q.iter().all(|r| r.len() == m));
    let split = m - m % 16;
    let mut out = [0i32; 4];
    // SAFETY: NEON is always present on aarch64, and every load and store is
    // within the first `split` elements of slices of length m.
    unsafe {
        let mut d = [vdupq_n_s32(0); 4];
        let mut c = 0;
        while c < split {
            let z0 = vld1q_s16(z.as_ptr().add(c));
            let z1 = vld1q_s16(z.as_ptr().add(c + 8));
            let p = acc.as_mut_ptr().add(c);
            let (mut a0, mut a1, mut a2, mut a3) =
                (vld1q_s32(p), vld1q_s32(p.add(4)), vld1q_s32(p.add(8)), vld1q_s32(p.add(12)));
            for r in 0..4 {
                let qv = vld1q_s8(q[r].as_ptr().add(c));
                let q0 = vmovl_s8(vget_low_s8(qv));
                let q1 = vmovl_high_s8(qv);
                d[r] = vmlal_s16(d[r], vget_low_s16(q0), vget_low_s16(z0));
                d[r] = vmlal_high_s16(d[r], q0, z0);
                d[r] = vmlal_s16(d[r], vget_low_s16(q1), vget_low_s16(z1));
                d[r] = vmlal_high_s16(d[r], q1, z1);
                a0 = vmlal_n_s16(a0, vget_low_s16(q0), zr[r]);
                a1 = vmlal_high_n_s16(a1, q0, zr[r]);
                a2 = vmlal_n_s16(a2, vget_low_s16(q1), zr[r]);
                a3 = vmlal_high_n_s16(a3, q1, zr[r]);
            }
            vst1q_s32(p, a0);
            vst1q_s32(p.add(4), a1);
            vst1q_s32(p.add(8), a2);
            vst1q_s32(p.add(12), a3);
            c += 16;
        }
        for r in 0..4 {
            out[r] = vaddvq_s32(d[r]);
        }
    }
    for c in split..m {
        for r in 0..4 {
            let qv = q[r][c] as i32;
            out[r] += qv * z[c] as i32;
            acc[c] += qv * zr[r] as i32;
        }
    }
    out
}

#[cfg(not(target_arch = "aarch64"))]
fn dot_axpy4(q: &[&[i8]; 4], z: &[i16], acc: &mut [i32], zr: &[i16; 4]) -> [i32; 4] {
    std::array::from_fn(|r| dot_axpy(q[r], z, acc, zr[r] as i32))
}

pub enum Coupling {
    /// J_ij = a_i a_j
    Rank1 { a: Vec<f64> },
    Dense(DenseJ),
    /// J_ij = f[|i - j|]: a chain whose couplings depend only on distance
    Chain { f: Vec<f64>, conv: fft::Conv },
    /// J_ij = c sum_mu xi_i^mu xi_j^mu (i != j) with xi = +-1: rank P
    LowRank { xi: Vec<Vec<i8>>, c: f64 },
}

impl Coupling {
    pub fn chain(f: Vec<f64>) -> Coupling {
        let conv = fft::Conv::new(&f);
        Coupling::Chain { f, conv }
    }

    pub fn n(&self) -> usize {
        match self {
            Coupling::Rank1 { a } => a.len(),
            Coupling::Dense(d) => d.n,
            Coupling::Chain { f, .. } => f.len(),
            Coupling::LowRank { xi, .. } => xi.first().map_or(0, |x| x.len()),
        }
    }

    /// Local fields h_i = sum_{j != i} J_ij z_j.
    pub fn fields(&self, z: &[f64], threads: usize) -> Vec<f64> {
        match self {
            Coupling::Rank1 { a } => {
                let s: f64 = a.iter().zip(z).map(|(a, z)| a * z).sum();
                a.iter().zip(z).map(|(a, z)| a * (s - a * z)).collect()
            }
            Coupling::Dense(d) => {
                // z rounded to Q1.15 for the product, as an FPGA would store it
                let zq: Vec<i16> = z.iter().map(|&v| cordic::q15(v)).collect();
                let f = d.scale / 32768.0;
                d.matvec(&zq, threads).into_iter().map(|h| h as f64 * f).collect()
            }
            Coupling::Chain { conv, .. } => conv.apply(z),
            Coupling::LowRank { xi, c } => {
                // h_i = c (sum_mu xi_i^mu U_mu - P z_i), U_mu = sum_j xi_j^mu z_j
                let mut h: Vec<f64> = z.iter().map(|&z| -(xi.len() as f64) * z).collect();
                for x in xi {
                    let u: f64 = x.iter().zip(z).map(|(&x, &z)| x as f64 * z).sum();
                    for (hk, &xk) in h.iter_mut().zip(x) {
                        *hk += u * xk as f64;
                    }
                }
                h.iter_mut().for_each(|v| *v *= c);
                h
            }
        }
    }

    /// Adds factor * J_ki to h_k for every k != i (column i of J).
    pub fn add_column(&self, i: usize, factor: f64, h: &mut [f64]) {
        match self {
            Coupling::Rank1 { a } => {
                let fa = factor * a[i];
                for (hk, &ak) in h.iter_mut().zip(a) {
                    *hk += fa * ak;
                }
                h[i] -= fa * a[i];
            }
            Coupling::Dense(d) => {
                let fs = factor * d.scale;
                for (hk, &q) in h.iter_mut().zip(&d.q[i * d.n..(i + 1) * d.n]) {
                    *hk += fs * q as f64;
                }
            }
            Coupling::Chain { f, .. } => {
                for (k, hk) in h.iter_mut().enumerate() {
                    if k != i {
                        *hk += factor * f[k.abs_diff(i)];
                    }
                }
            }
            Coupling::LowRank { xi, c } => {
                for x in xi {
                    let fx = factor * c * x[i] as f64;
                    for (hk, &xk) in h.iter_mut().zip(x) {
                        *hk += fx * xk as f64;
                    }
                }
                h[i] -= factor * c * xi.len() as f64;
            }
        }
    }

    /// Size of a typical field, sqrt(sum_{i != j} J_ij^2 / N). Schedules
    /// divide by it so the same gamma works for any problem and size.
    pub fn field_scale(&self) -> f64 {
        let n = self.n() as f64;
        match self {
            Coupling::Rank1 { a } => {
                let a2: f64 = a.iter().map(|a| a * a).sum();
                let s: f64 = a.iter().map(|a| a * a * (a2 - a * a)).sum();
                (s / n).sqrt()
            }
            Coupling::Dense(d) => {
                let s: f64 = d.q.iter().map(|&q| (q as f64) * (q as f64)).sum();
                (s / n).sqrt() * d.scale
            }
            Coupling::Chain { f, .. } => {
                let s: f64 = (1..f.len()).map(|d| 2.0 * (n - d as f64) * f[d] * f[d]).sum();
                (s / n).sqrt()
            }
            Coupling::LowRank { xi, c } => {
                // sum_{i != j} (sum_mu xi_i xi_j)^2 = |O|^2 - N P^2, O = xi xi^T
                let p = xi.len();
                let s = if (p * p) as f64 * n > 2e9 {
                    n * (n - 1.0) * p as f64 // random patterns, on average
                } else {
                    let mut o2 = 0.0;
                    for a in xi {
                        for b in xi {
                            let o: i64 = a.iter().zip(b).map(|(&x, &y)| (x * y) as i64).sum();
                            o2 += (o * o) as f64;
                        }
                    }
                    o2 - n * (p * p) as f64
                };
                c.abs() * (s / n).sqrt()
            }
        }
    }
}

/// Rotates spin i about z by c * h_i, then every spin about x by phi.
pub fn rotate(spins: &mut [Spin], h: &[f64], c: f64, phi: f64, threads: usize) {
    let (sb, cb) = phi.sin_cos();
    let threads = par::threads_for(spins.len() * 40, threads);
    par::for_chunks(spins, threads, |off, chunk| {
        // in groups of L, so the sin/cos polynomials run in SIMD lanes
        const L: usize = 4;
        let h = &h[off..off + chunk.len()];
        let split = chunk.len() - chunk.len() % L;
        for (sp, hc) in chunk[..split].chunks_exact_mut(L).zip(h[..split].chunks_exact(L)) {
            let theta: [f64; L] = std::array::from_fn(|k| c * hc[k]);
            let (st, ct) = sincos_lanes(&theta);
            for k in 0..L {
                turn(&mut sp[k], st[k], ct[k], sb, cb);
            }
        }
        for (sp, &hi) in chunk[split..].iter_mut().zip(&h[split..]) {
            let (st, ct) = sincos(c * hi);
            turn(sp, st, ct, sb, cb);
        }
    });
}

#[inline(always)]
fn turn(sp: &mut Spin, st: f64, ct: f64, sb: f64, cb: f64) {
    let x1 = sp[0] * ct - sp[1] * st;
    let y1 = sp[0] * st + sp[1] * ct;
    sp[0] = x1;
    sp[1] = y1 * cb - sp[2] * sb;
    sp[2] = y1 * sb + sp[2] * cb;
}

#[inline(always)]
fn sincos_lanes<const L: usize>(x: &[f64; L]) -> ([f64; L], [f64; L]) {
    let mut s = [0.0; L];
    let mut c = [0.0; L];
    for k in 0..L {
        (s[k], c[k]) = sincos(x[k]);
    }
    (s, c)
}

/// sin and cos without branches, so loops over spins vectorize: reduce to
/// |r| <= pi/4 (two-part pi/2, fine for |x| up to about 10^6), then Taylor
/// series to r^14 (error below 1e-14).
#[inline(always)]
pub fn sincos(x: f64) -> (f64, f64) {
    const TWO_OVER_PI: f64 = 0.636_619_772_367_581_3;
    const PIO2_HI: f64 = 1.570_796_326_734_125_6; // first 33 bits of pi/2
    const PIO2_LO: f64 = 6.077_100_506_506_192e-11;
    let k = (x * TWO_OVER_PI).round();
    let r = (x - k * PIO2_HI) - k * PIO2_LO;
    let r2 = r * r;
    let s = r * (1.0
        + r2 * (-1.0 / 6.0
            + r2 * (1.0 / 120.0
                + r2 * (-1.0 / 5040.0
                    + r2 * (1.0 / 362_880.0 + r2 * (-1.0 / 39_916_800.0 + r2 * (1.0 / 6_227_020_800.0)))))));
    let c = 1.0
        + r2 * (-0.5
            + r2 * (1.0 / 24.0
                + r2 * (-1.0 / 720.0
                    + r2 * (1.0 / 40_320.0
                        + r2 * (-1.0 / 3_628_800.0
                            + r2 * (1.0 / 479_001_600.0 + r2 * (-1.0 / 87_178_291_200.0)))))));
    // quadrant k mod 4: swap sin and cos for odd k, then fix the signs
    let q = k as i64;
    let (s, c) = if q & 1 == 1 { (c, s) } else { (s, c) };
    let s = if q & 2 == 2 { -s } else { s };
    let c = if (q + 1) & 2 == 2 { -c } else { c };
    (s, c)
}

pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upper_triangle_matvec_matches_full() {
        let mut rng = Rng::new(7);
        for &(n, threads) in &[(1, 1), (5, 1), (37, 3), (300, 8), (700, 4)] {
            let d = DenseJ::sk(n, true, &mut rng);
            let z: Vec<i16> = (0..n).map(|_| rng.range(0, 65535) as u16 as i16).collect();
            let h = d.matvec(&z, threads);
            for i in 0..n {
                let want: i64 = (0..n).map(|j| d.q[i * n + j] as i64 * z[j] as i64).sum();
                assert_eq!(h[i], want, "n {n} row {i}");
                assert_eq!(d.q[i * n + i], 0);
            }
        }
    }

    #[test]
    fn fast_sincos() {
        let mut x = -2000.0;
        while x < 2000.0 {
            let (s, c) = sincos(x);
            assert!((s - x.sin()).abs() < 1e-12 && (c - x.cos()).abs() < 1e-12, "x = {x}");
            x += 0.0137;
        }
    }

    #[test]
    fn columns_match_fields() {
        let mut rng = Rng::new(12);
        let n = 40;
        let xi: Vec<Vec<i8>> = (0..5).map(|_| (0..n).map(|_| if rng.coin() { 1 } else { -1 }).collect()).collect();
        let couplings = [
            Coupling::Rank1 { a: (0..n).map(|_| rng.uniform()).collect() },
            Coupling::Dense(DenseJ::sk(n, true, &mut rng)),
            Coupling::chain((0..n).map(|d| if d == 0 { 0.0 } else { (d as f64).powf(-1.3) }).collect()),
            Coupling::LowRank { xi, c: -1.0 / n as f64 },
        ];
        for c in &couplings {
            // column i of J is the field from a unit vector at i
            for i in [0, 7, n - 1] {
                // 0.5, since 1.0 saturates in Q1.15 (Dense rounds z to 16 bits)
                let mut e = vec![0.0; n];
                e[i] = 0.5;
                let want = c.fields(&e, 1);
                let mut got = vec![0.0; n];
                c.add_column(i, 0.5, &mut got);
                for k in 0..n {
                    assert!((got[k] - want[k]).abs() < 1e-9, "i {i} k {k}: {} vs {}", got[k], want[k]);
                }
            }
            // field_scale against the definition
            let mut s = 0.0;
            for i in 0..n {
                let mut col = vec![0.0; n];
                c.add_column(i, 1.0, &mut col);
                s += col.iter().map(|v| v * v).sum::<f64>();
            }
            assert!(((s / n as f64).sqrt() - c.field_scale()).abs() < 1e-6 * c.field_scale());
        }
    }

    #[test]
    fn rank1_fields() {
        let a = vec![0.5, 1.0, 1.5, 0.25];
        let z = vec![0.1, -0.7, 0.3, 0.9];
        let h = Coupling::Rank1 { a: a.clone() }.fields(&z, 1);
        for i in 0..4 {
            let want: f64 = (0..4).filter(|&j| j != i).map(|j| a[i] * a[j] * z[j]).sum();
            assert!((h[i] - want).abs() < 1e-12);
        }
    }
}
