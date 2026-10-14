//! Fast convolution for couplings that depend only on distance,
//! J_ij = f[|i - j|]: then h = J z is a convolution, O(N log N) by FFT
//! instead of O(N^2).

use std::f64::consts::PI;
use std::ops::{Add, Mul, Sub};

#[derive(Clone, Copy, Default)]
struct C {
    re: f64,
    im: f64,
}

impl Add for C {
    type Output = C;
    fn add(self, o: C) -> C {
        C { re: self.re + o.re, im: self.im + o.im }
    }
}

impl Sub for C {
    type Output = C;
    fn sub(self, o: C) -> C {
        C { re: self.re - o.re, im: self.im - o.im }
    }
}

impl Mul for C {
    type Output = C;
    fn mul(self, o: C) -> C {
        C { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
}

/// In-place radix-2 FFT; `tw[k]` = exp(-2 pi i k / len).
fn fft(a: &mut [C], tw: &[C], inverse: bool) {
    let n = a.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let (half, step) = (len / 2, n / len);
        for start in (0..n).step_by(len) {
            for k in 0..half {
                let mut w = tw[k * step];
                if inverse {
                    w.im = -w.im;
                }
                let u = a[start + k];
                let v = a[start + k + half] * w;
                a[start + k] = u + v;
                a[start + k + half] = u - v;
            }
        }
        len <<= 1;
    }
}

/// h_i = sum_{j != i} f[|i - j|] z_j for i, j < n.
pub struct Conv {
    n: usize,
    tw: Vec<C>,
    g_hat: Vec<C>,
}

impl Conv {
    /// f[d] for d = 0..n (f[0] is ignored: no self-coupling).
    pub fn new(f: &[f64]) -> Conv {
        let n = f.len();
        let len = (2 * n).next_power_of_two().max(2);
        let tw: Vec<C> = (0..len / 2)
            .map(|k| {
                let a = -2.0 * PI * k as f64 / len as f64;
                C { re: a.cos(), im: a.sin() }
            })
            .collect();
        // circular kernel: g[d] = g[len - d] = f[d]
        let mut g = vec![C::default(); len];
        for d in 1..n {
            g[d].re = f[d];
            g[len - d].re = f[d];
        }
        fft(&mut g, &tw, false);
        Conv { n, tw, g_hat: g }
    }

    pub fn apply(&self, z: &[f64]) -> Vec<f64> {
        let len = self.g_hat.len();
        let mut a = vec![C::default(); len];
        for (ai, &zi) in a.iter_mut().zip(z) {
            ai.re = zi;
        }
        fft(&mut a, &self.tw, false);
        for (ai, &g) in a.iter_mut().zip(&self.g_hat) {
            *ai = *ai * g;
        }
        fft(&mut a, &self.tw, true);
        a[..self.n].iter().map(|c| c.re / len as f64).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn matches_direct_sum() {
        let mut rng = Rng::new(2);
        for n in [1usize, 2, 3, 17, 100, 513] {
            let f: Vec<f64> = (0..n).map(|d| if d == 0 { 9.0 } else { rng.normal() }).collect();
            let z: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
            let h = Conv::new(&f).apply(&z);
            for i in 0..n {
                let want: f64 = (0..n).filter(|&j| j != i).map(|j| f[i.abs_diff(j)] * z[j]).sum();
                assert!((h[i] - want).abs() < 1e-9 * (1.0 + want.abs()), "n {n} i {i}");
            }
        }
    }
}
