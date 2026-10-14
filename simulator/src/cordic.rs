//! Bit-exact model of cordic_iterative.v (and so of cordic_pipelined.v, which
//! the testbenches check against it).
//!
//! Angles are 32 bits for a full turn (2^30 = 90 degrees). x and y are Q1.15.
//! The inputs are scaled by 1/K ~ 0.6074 first (shifts and adds), then 15
//! iterations run with 17-bit X and Y, and the outputs saturate to 16 bits.

/// atan(2^-i) in turns * 2^32, exactly as in cordic_iterative.v.
pub const ATAN: [i32; 15] = [
    0x2000_0000, 0x12E4_051D, 0x09FB_385B, 0x0511_11D4, 0x028B_0D43, 0x0145_D7E1, 0x00A2_F61E,
    0x0051_7C55, 0x0028_BE53, 0x0014_5F2E, 0x000A_2F98, 0x0005_17CC, 0x0002_8BE6, 0x0001_45F3,
    0x0000_A2F9,
];

/// Keep the low 17 bits, sign-extended (a 17-bit Verilog reg).
#[inline(always)]
fn wrap17(v: i32) -> i32 {
    (v << 15) >> 15
}

/// 1/K ~ 1/2 + 1/8 - 1/64 - 1/512, on the 17-bit sign-extended input.
#[inline(always)]
fn scale(v: i32) -> i32 {
    wrap17((v >> 1) + (v >> 3) - (v >> 6) - (v >> 9))
}

#[inline(always)]
fn saturate16(v: i32) -> i16 {
    v.clamp(-32768, 32767) as i16
}

/// Rotates (x, y) counterclockwise by `angle`.
#[inline]
pub fn rotate(x: i16, y: i16, angle: i32) -> (i16, i16) {
    let xs = scale(x as i32);
    let ys = scale(y as i32);

    // quadrant pre-rotation by +-90 degrees, so the iterations see |angle| <= 90
    let low = (angle as u32) & 0x3FFF_FFFF;
    let (mut cx, mut cy, mut z) = match (angle as u32) >> 30 {
        0 | 3 => (xs, ys, angle),
        1 => (wrap17(-ys), xs, low as i32),
        _ => (ys, wrap17(-xs), (low | 0xC000_0000) as i32),
    };

    for (i, &atan) in ATAN.iter().enumerate() {
        let xsh = cx >> i;
        let ysh = cy >> i;
        if z < 0 {
            cx = wrap17(cx + ysh);
            cy = wrap17(cy - xsh);
            z = z.wrapping_add(atan);
        } else {
            cx = wrap17(cx - ysh);
            cy = wrap17(cy + xsh);
            z = z.wrapping_sub(atan);
        }
    }
    (saturate16(cx), saturate16(cy))
}

/// Spins rotated together by rotate_lanes.
pub const LANES: usize = 16;

/// rotate() for LANES rotations at once, in NEON registers (the compiler
/// doesn't vectorize the plain version). Same results bit for bit (checked by
/// the tests below).
#[cfg(target_arch = "aarch64")]
#[inline]
pub fn rotate_lanes(
    x: &[i16; LANES],
    y: &[i16; LANES],
    angle: &[i32; LANES],
) -> ([i16; LANES], [i16; LANES]) {
    use std::arch::aarch64::*;
    let mut ox = [0i16; LANES];
    let mut oy = [0i16; LANES];
    // SAFETY: NEON is always present on aarch64; loads and stores cover
    // exactly the LANES = 16 elements of each array.
    unsafe {
        let xv = [vld1q_s16(x.as_ptr()), vld1q_s16(x.as_ptr().add(8))];
        let yv = [vld1q_s16(y.as_ptr()), vld1q_s16(y.as_ptr().add(8))];
        let av = std::array::from_fn(|g| vld1q_s32(angle.as_ptr().add(4 * g)));
        let (rx, ry) = neon::rotate(xv, yv, av);
        for g in 0..2 {
            vst1q_s16(ox.as_mut_ptr().add(8 * g), rx[g]);
            vst1q_s16(oy.as_mut_ptr().add(8 * g), ry[g]);
        }
    }
    (ox, oy)
}

/// The CORDIC on NEON registers: 16 lanes as two int16x8 for x and y and
/// four int32x4 for the angles. X and Y never leave 17 bits (|X|, |Y| <
/// 1.65 * 0.61 * 2^15 * sqrt 2), so the 17-bit wraparound of the Verilog
/// never happens and 32-bit lanes give the same values.
#[cfg(target_arch = "aarch64")]
pub mod neon {
    use super::ATAN;
    use std::arch::aarch64::*;

    #[inline(always)]
    unsafe fn scale(v: int32x4_t) -> int32x4_t {
        vsubq_s32(
            vsubq_s32(vaddq_s32(vshrq_n_s32::<1>(v), vshrq_n_s32::<3>(v)), vshrq_n_s32::<6>(v)),
            vshrq_n_s32::<9>(v),
        )
    }

    /// # Safety
    /// Needs NEON (always present on aarch64).
    #[inline(always)]
    pub unsafe fn rotate(
        x: [int16x8_t; 2],
        y: [int16x8_t; 2],
        angle: [int32x4_t; 4],
    ) -> ([int16x8_t; 2], [int16x8_t; 2]) {
        let wide = |v: [int16x8_t; 2]| {
            [vmovl_s16(vget_low_s16(v[0])), vmovl_high_s16(v[0]), vmovl_s16(vget_low_s16(v[1])), vmovl_high_s16(v[1])]
        };
        let (xw, yw) = (wide(x), wide(y));
        let mut cx = [vdupq_n_s32(0); 4];
        let mut cy = [vdupq_n_s32(0); 4];
        let mut z = [vdupq_n_s32(0); 4];
        for g in 0..4 {
            // quadrant pre-rotation
            let (xs, ys, a) = (scale(xw[g]), scale(yw[g]), angle[g]);
            let q = vshrq_n_u32::<30>(vreinterpretq_u32_s32(a));
            let is1 = vceqq_u32(q, vdupq_n_u32(1));
            let is2 = vceqq_u32(q, vdupq_n_u32(2));
            let low = vandq_s32(a, vdupq_n_s32(0x3FFF_FFFF));
            cx[g] = vbslq_s32(is1, vnegq_s32(ys), vbslq_s32(is2, ys, xs));
            cy[g] = vbslq_s32(is1, xs, vbslq_s32(is2, vnegq_s32(xs), ys));
            z[g] = vbslq_s32(is1, low, vbslq_s32(is2, vorrq_s32(low, vdupq_n_s32(0xC000_0000u32 as i32)), a));
        }
        for (i, &atan) in ATAN.iter().enumerate() {
            let sh = vdupq_n_s32(-(i as i32));
            let at = vdupq_n_s32(atan);
            for g in 0..4 {
                let neg = vcltzq_s32(z[g]);
                let xsh = vshlq_s32(cx[g], sh);
                let ysh = vshlq_s32(cy[g], sh);
                cx[g] = vbslq_s32(neg, vaddq_s32(cx[g], ysh), vsubq_s32(cx[g], ysh));
                cy[g] = vbslq_s32(neg, vsubq_s32(cy[g], xsh), vaddq_s32(cy[g], xsh));
                z[g] = vbslq_s32(neg, vaddq_s32(z[g], at), vsubq_s32(z[g], at));
            }
        }
        // saturating narrow to 16 bits, like resize_output
        let narrow = |v: [int32x4_t; 4]| {
            [vcombine_s16(vqmovn_s32(v[0]), vqmovn_s32(v[1])), vcombine_s16(vqmovn_s32(v[2]), vqmovn_s32(v[3]))]
        };
        (narrow(cx), narrow(cy))
    }
}

/// rotate() for LANES rotations at once, without branches. Same results
/// bit for bit (checked by the tests below).
#[cfg(not(target_arch = "aarch64"))]
#[inline]
pub fn rotate_lanes(
    x: &[i16; LANES],
    y: &[i16; LANES],
    angle: &[i32; LANES],
) -> ([i16; LANES], [i16; LANES]) {
    let (mut cx, mut cy, mut z) = ([0i32; LANES], [0i32; LANES], [0i32; LANES]);
    for k in 0..LANES {
        let xs = scale(x[k] as i32);
        let ys = scale(y[k] as i32);
        let q = (angle[k] as u32) >> 30;
        let low = (angle[k] as u32) & 0x3FFF_FFFF;
        cx[k] = wrap17(if q == 1 { -ys } else if q == 2 { ys } else { xs });
        cy[k] = wrap17(if q == 1 { xs } else if q == 2 { -xs } else { ys });
        z[k] = if q == 1 {
            low as i32
        } else if q == 2 {
            (low | 0xC000_0000) as i32
        } else {
            angle[k]
        };
    }
    for (i, &atan) in ATAN.iter().enumerate() {
        for k in 0..LANES {
            // m = 0 if z < 0, else -1; (v ^ m) - m is then v or -v
            let m = !(z[k] >> 31);
            let xsh = cx[k] >> i;
            let ysh = cy[k] >> i;
            cx[k] = wrap17(cx[k] + ((ysh ^ m) - m));
            cy[k] = wrap17(cy[k] - ((xsh ^ m) - m));
            z[k] = z[k].wrapping_add((atan ^ m).wrapping_sub(m));
        }
    }
    (cx.map(saturate16), cy.map(saturate16))
}

/// Real number -> Q1.15, rounded and saturated (like gen_spin_init.py).
pub fn q15(v: f64) -> i16 {
    (v * 32768.0).round().clamp(-32768.0, 32767.0) as i16
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn atan_table_matches_formula() {
        for (i, &a) in ATAN.iter().enumerate() {
            let want = (0.5f64.powi(i as i32)).atan() / (2.0 * PI) * 2f64.powi(32);
            assert!((a as f64 - want).abs() < 2.0, "entry {i}: {a:x} vs {want}");
        }
    }

    #[test]
    fn lanes_match_scalar() {
        let mut rng = crate::rng::Rng::new(9);
        for _ in 0..200_000 {
            let mut pick = || match rng.range(0, 9) {
                0 => -32768,
                1 => 32767,
                _ => rng.range(0, 65535) as u16 as i16,
            };
            let x: [i16; LANES] = std::array::from_fn(|_| pick());
            let y: [i16; LANES] = std::array::from_fn(|_| pick());
            let a: [i32; LANES] = std::array::from_fn(|_| rng.next_u64() as i32);
            let (rx, ry) = rotate_lanes(&x, &y, &a);
            for k in 0..LANES {
                assert_eq!((rx[k], ry[k]), rotate(x[k], y[k], a[k]));
            }
        }
    }

    #[test]
    fn rotation_is_close_to_exact() {
        let mut worst = 0.0f64;
        for k in 0..2000 {
            let angle = (k as u32).wrapping_mul(2_147_483_647).wrapping_add(12345) as i32;
            let theta = angle as f64 / 2f64.powi(32) * 2.0 * PI;
            let (x, y) = (20000i16, -9000i16);
            let (rx, ry) = rotate(x, y, angle);
            let ex = x as f64 * theta.cos() - y as f64 * theta.sin();
            let ey = x as f64 * theta.sin() + y as f64 * theta.cos();
            worst = worst.max((rx as f64 - ex).abs()).max((ry as f64 - ey).abs());
        }
        // 1/K is approximated to 3e-4, so errors are a few LSB per 10^4
        assert!(worst < 12.0, "worst error {worst} LSB");
    }
}
