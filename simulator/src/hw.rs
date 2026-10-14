//! Bit-exact model of VIRAL_n_spin.v: same integers, same rounding, same
//! results, so its output can be checked against the Verilog simulation
//! (simulator/check_hw.sh).
//!
//! Per layer, with w_j = a_j * z_j rounded to Q1.15 and S = sum of all w_j:
//!   angle_i = gamma * a_i * (S - w_i), rounded to the 32-bit angle
//!   (x, y) rotated by angle_i, then (y, z) rotated by beta
//! Energy = (S^2 - sum w_i^2) >> 16, Q1.15 (= sum_{i<j} w_i w_j).

use crate::{cordic, par};

pub type HwSpin = [i16; 3];

/// angle_lut.v: gamma (Q1.31) and beta (32-bit angle) for layers 0..3.
pub const LUT_GAMMA: [i32; 4] = [0xA57E_0000u32 as i32, 0x2000_0000, 0, 0];
pub const LUT_BETA: [i32; 4] = [0x2000_0000, 0x4000_0000, 0, 0];

/// w = a * z, rounded to Q1.15 (17 bits, since a < 2).
#[inline(always)]
pub fn w_of(a: u16, z: i16) -> i64 {
    (a as i64 * z as i64 + 16384) >> 15
}

/// Bits [61:30] of gamma * a * h + 2^29: the z rotation angle. Only the
/// product mod 2^62 matters (as in the design), so 64-bit wrapping is exact.
#[inline(always)]
pub fn z_angle(gamma: i32, a: u16, h: i64) -> i32 {
    let p = (gamma as i64).wrapping_mul(a as i64).wrapping_mul(h).wrapping_add(1 << 29);
    ((p as u64) >> 30) as u32 as i32
}

#[inline]
pub fn spin_update(sp: HwSpin, a: u16, s: i64, gamma: i32, beta: i32) -> HwSpin {
    let h = s - w_of(a, sp[2]);
    let (x1, y1) = cordic::rotate(sp[0], sp[1], z_angle(gamma, a, h));
    let (y2, z2) = cordic::rotate(y1, sp[2], beta);
    [x1, y2, z2]
}

/// Verilog $clog2.
pub fn clog2(n: u64) -> u64 {
    if n <= 1 {
        0
    } else {
        64 - (n - 1).leading_zeros() as u64
    }
}

/// Clock count of VIRAL_n_spin from reset to finished.
pub fn clocks(n: u64, layers: u64) -> u64 {
    let w = 17 + clog2(n);
    (layers + 1) * (n + w) + 63 * layers + 17
}

pub struct Hw {
    pub spins: Vec<HwSpin>,
    /// weights, unsigned with 15 fractional bits (0x8000 = 1.0)
    pub a: Vec<u16>,
}

impl Hw {
    pub fn sum_w(&self) -> i64 {
        self.spins.iter().zip(&self.a).map(|(s, &a)| w_of(a, s[2])).sum()
    }

    /// S^2 - sum w^2 (the design's energy_diff).
    pub fn energy_diff(&self) -> i128 {
        // |w| < 2^16, so both sums fit 64 bits for up to 2^21 spins
        let (mut s, mut q) = (0i64, 0i64);
        for (sp, &a) in self.spins.iter().zip(&self.a) {
            let w = w_of(a, sp[2]);
            s += w;
            q += w * w;
        }
        s as i128 * s as i128 - q as i128
    }

    /// total_energy: sum_{i<j} w_i w_j in Q1.15.
    pub fn energy(&self) -> i64 {
        (self.energy_diff() >> 16) as i64
    }

    /// The 8 LEDs: energy >>> LED_LSB, saturated.
    pub fn led(&self, led_lsb: u32) -> u8 {
        let e = self.energy_diff() >> (16 + led_lsb);
        if e > 127 {
            0x7F
        } else if e < -128 {
            0x80
        } else {
            e as i8 as u8
        }
    }

    /// One layer. Every spin sees S from the start of the layer.
    pub fn layer(&mut self, gamma: i32, beta: i32, threads: usize) {
        let s = self.sum_w();
        let a = &self.a;
        let threads = par::threads_for(self.spins.len() * 60, threads);
        par::for_chunks(&mut self.spins, threads, |off, chunk| {
            update_spins(chunk, &a[off..off + chunk.len()], s, gamma, beta);
        });
    }
}

/// spin_update for a run of spins, LANES at a time.
fn update_spins(spins: &mut [HwSpin], a: &[u16], s: i64, gamma: i32, beta: i32) {
    const L: usize = cordic::LANES;
    let split = spins.len() - spins.len() % L;
    for (sp, aw) in spins[..split].chunks_exact_mut(L).zip(a[..split].chunks_exact(L)) {
        let ang: [i32; L] = std::array::from_fn(|k| z_angle(gamma, aw[k], s - w_of(aw[k], sp[k][2])));
        update_lanes(sp.try_into().unwrap(), &ang, beta);
    }
    for (sp, &aw) in spins[split..].iter_mut().zip(&a[split..]) {
        *sp = spin_update(*sp, aw, s, gamma, beta);
    }
}

/// Both rotations for LANES spins, given their z angles.
#[cfg(target_arch = "aarch64")]
#[inline(always)]
fn update_lanes(sp: &mut [HwSpin; cordic::LANES], ang: &[i32; cordic::LANES], beta: i32) {
    use std::arch::aarch64::*;
    // SAFETY: NEON is always present on aarch64. The 16 spins are 48
    // contiguous i16 ({x, y, z} each), read and written as two groups of 8
    // with the interleaving loads and stores.
    unsafe {
        let p = sp.as_mut_ptr() as *mut i16;
        let (v0, v1) = (vld3q_s16(p), vld3q_s16(p.add(24)));
        let av = std::array::from_fn(|g| {
            vcombine_s32(
                vcreate_s32(ang[4 * g] as u32 as u64 | (ang[4 * g + 1] as u32 as u64) << 32),
                vcreate_s32(ang[4 * g + 2] as u32 as u64 | (ang[4 * g + 3] as u32 as u64) << 32),
            )
        });
        let (x1, y1) = cordic::neon::rotate([v0.0, v1.0], [v0.1, v1.1], av);
        let (y2, z2) = cordic::neon::rotate(y1, [v0.2, v1.2], [vdupq_n_s32(beta); 4]);
        vst3q_s16(p, int16x8x3_t(x1[0], y2[0], z2[0]));
        vst3q_s16(p.add(24), int16x8x3_t(x1[1], y2[1], z2[1]));
    }
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(always)]
fn update_lanes(sp: &mut [HwSpin; cordic::LANES], ang: &[i32; cordic::LANES], beta: i32) {
    const L: usize = cordic::LANES;
    let x: [i16; L] = std::array::from_fn(|k| sp[k][0]);
    let y: [i16; L] = std::array::from_fn(|k| sp[k][1]);
    let z: [i16; L] = std::array::from_fn(|k| sp[k][2]);
    let (x1, y1) = cordic::rotate_lanes(&x, &y, ang);
    let (y2, z2) = cordic::rotate_lanes(&y1, &z, &[beta; L]);
    for k in 0..L {
        sp[k] = [x1[k], y2[k], z2[k]];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn energy_is_pairwise_sum() {
        let mut rng = Rng::new(3);
        let n = 37;
        let hw = Hw {
            spins: (0..n).map(|_| [0, 0, rng.range(0, 65535) as u16 as i16]).collect(),
            a: (0..n).map(|_| rng.range(1, 0xFFFF) as u16).collect(),
        };
        let w: Vec<i64> = hw.spins.iter().zip(&hw.a).map(|(s, &a)| w_of(a, s[2])).collect();
        let mut pairs = 0i64;
        for i in 0..n {
            for j in i + 1..n {
                pairs += w[i] * w[j];
            }
        }
        assert_eq!(hw.energy(), pairs >> 15);
    }

    #[test]
    fn z_angle_matches_wide_product() {
        let mut rng = Rng::new(4);
        for _ in 0..100_000 {
            let g = rng.next_u64() as i32;
            let a = rng.range(0, 65535) as u16;
            let h = (rng.next_u64() as i64) >> 25; // up to 39 bits: 4 million spins
            let want = ((g as i128 * a as i128 * h as i128 + (1i128 << 29)) >> 30) as i32;
            assert_eq!(z_angle(g, a, h), want);
        }
    }

    #[test]
    fn layer_matches_spin_update() {
        let mut rng = Rng::new(8);
        let n = 101;
        let mut m = Hw {
            spins: (0..n).map(|_| [0; 3].map(|_: i16| rng.range(0, 65535) as u16 as i16)).collect(),
            a: (0..n).map(|_| rng.range(1, 0xFFFF) as u16).collect(),
        };
        let s = m.sum_w();
        let want: Vec<HwSpin> =
            m.spins.iter().zip(&m.a).map(|(&sp, &a)| spin_update(sp, a, s, 0x1234_5678, -0x2345_6789)).collect();
        m.layer(0x1234_5678, -0x2345_6789, 3);
        assert_eq!(m.spins, want);
    }

    #[test]
    fn clock_formula() {
        // 49,152 spins, 2 layers: the Cyclone 10 GX case in the README
        assert_eq!(clocks(49_152, 2), 147_698);
    }
}
