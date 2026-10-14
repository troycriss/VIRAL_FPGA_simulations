//! Classical methods to judge the mean-field result against.

use crate::model::{Coupling, DenseJ};
use crate::rng::Rng;
use std::collections::BinaryHeap;

/// Greedy: largest number first, into the set with the smaller sum.
/// Returns the residual |sum A - sum B|.
pub fn greedy(weights: &[u16]) -> u64 {
    let mut w: Vec<u64> = weights.iter().map(|&w| w as u64).collect();
    w.sort_unstable_by(|a, b| b.cmp(a));
    let (mut a, mut b) = (0u64, 0u64);
    for x in w {
        if a <= b {
            a += x;
        } else {
            b += x;
        }
    }
    a.abs_diff(b)
}

/// Karmarkar-Karp differencing: keep replacing the two largest numbers by
/// their difference. Returns the residual.
pub fn karmarkar_karp(weights: &[u16]) -> u64 {
    let mut heap: BinaryHeap<u64> = weights.iter().map(|&w| w as u64).collect();
    while heap.len() > 1 {
        let a = heap.pop().unwrap();
        let b = heap.pop().unwrap();
        heap.push(a - b);
    }
    heap.pop().unwrap_or(0)
}

/// Smallest possible residual, by subset-sum dynamic programming on the
/// integer weights. None if it would take more than `max_work` word operations.
pub fn optimum(weights: &[u16], max_work: u64) -> Option<u64> {
    let total: u64 = weights.iter().map(|&w| w as u64).sum();
    let half = (total / 2) as usize;
    let words = half / 64 + 1;
    if weights.len() as u64 * words as u64 > max_work {
        return None;
    }
    let mut bits = vec![0u64; words];
    bits[0] = 1;
    for &w in weights {
        shift_or(&mut bits, w as usize);
    }
    (0..=half)
        .rev()
        .find(|&s| (bits[s / 64] >> (s % 64)) & 1 == 1)
        .map(|s| total - 2 * s as u64)
}

/// bits |= bits << shift
fn shift_or(bits: &mut [u64], shift: usize) {
    let (ws, bs) = (shift / 64, shift % 64);
    for i in (ws..bits.len()).rev() {
        let src = i - ws;
        let mut v = bits[src] << bs;
        if bs > 0 && src > 0 {
            v |= bits[src - 1] >> (64 - bs);
        }
        bits[i] |= v;
    }
}

/// Moves single numbers to the other set while that shrinks the residual.
/// Updates `signs` and returns the new residual.
pub fn partition_local(weights: &[u16], signs: &mut [i8]) -> u64 {
    let mut r: i64 = weights.iter().zip(signs.iter()).map(|(&w, &s)| w as i64 * s as i64).sum();
    loop {
        let mut best = (r.abs(), usize::MAX);
        for (i, (&w, &s)) in weights.iter().zip(signs.iter()).enumerate() {
            let nr = (r - 2 * w as i64 * s as i64).abs();
            if nr < best.0 {
                best = (nr, i);
            }
        }
        if best.1 == usize::MAX {
            return r.unsigned_abs();
        }
        let i = best.1;
        r -= 2 * weights[i] as i64 * signs[i] as i64;
        signs[i] = -signs[i];
    }
}

/// Steepest descent by single spin flips until no flip lowers the energy.
/// Updates `s` and returns sum_{i<j} q_ij s_i s_j (multiply by J's scale).
pub fn sk_descent(j: &DenseJ, s: &mut [i8], threads: usize) -> i64 {
    let n = j.n;
    let mut h = j.fields_int(s, threads);
    loop {
        // flipping spin i changes the energy by -2 s_i h_i
        let mut best = (0i64, usize::MAX);
        for i in 0..n {
            let g = s[i] as i64 * h[i] as i64;
            if g > best.0 {
                best = (g, i);
            }
        }
        if best.1 == usize::MAX {
            break;
        }
        let i = best.1;
        let si = s[i] as i32;
        s[i] = -s[i];
        for (hk, &q) in h.iter_mut().zip(&j.q[i * n..(i + 1) * n]) {
            *hk -= 2 * si * q as i32;
        }
    }
    h.iter().zip(s.iter()).map(|(&h, &s)| h as i64 * s as i64).sum::<i64>() / 2
}

/// Steepest descent by single spin flips for any coupling, in floating
/// point. Updates `s` and returns the energy sum_{i<j} J_ij s_i s_j.
pub fn descent(c: &Coupling, s: &mut [i8], threads: usize) -> f64 {
    let z: Vec<f64> = s.iter().map(|&v| v as f64).collect();
    let mut h = c.fields(&z, threads);
    let tol = 1e-9 * c.field_scale();
    for _ in 0..50 * s.len().max(1) {
        // flipping spin i changes the energy by -2 s_i h_i
        let mut best = (tol, usize::MAX);
        for (i, (&si, &hi)) in s.iter().zip(&h).enumerate() {
            let g = si as f64 * hi;
            if g > best.0 {
                best = (g, i);
            }
        }
        if best.1 == usize::MAX {
            break;
        }
        let i = best.1;
        let old = s[i] as f64;
        s[i] = -s[i];
        c.add_column(i, -2.0 * old, &mut h);
    }
    0.5 * s.iter().zip(&h).map(|(&s, &h)| s as f64 * h).sum::<f64>()
}

/// sk_descent from a random start.
pub fn sk_random_descent(j: &DenseJ, rng: &mut Rng, threads: usize) -> i64 {
    let mut s: Vec<i8> = (0..j.n).map(|_| if rng.coin() { 1 } else { -1 }).collect();
    sk_descent(j, &mut s, threads)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_force(w: &[u16]) -> u64 {
        let n = w.len();
        (0..1u32 << n)
            .map(|m| {
                let r: i64 = (0..n).map(|i| if m >> i & 1 == 1 { w[i] as i64 } else { -(w[i] as i64) }).sum();
                r.unsigned_abs()
            })
            .min()
            .unwrap()
    }

    #[test]
    fn optimum_matches_brute_force() {
        let mut rng = Rng::new(11);
        for n in 1..=14 {
            let w: Vec<u16> = (0..n).map(|_| rng.range(1, 0x8000) as u16).collect();
            let best = brute_force(&w);
            assert_eq!(optimum(&w, u64::MAX), Some(best), "n = {n}");
            assert!(karmarkar_karp(&w) >= best && greedy(&w) >= best);
            let mut s = vec![1i8; n];
            assert!(partition_local(&w, &mut s) >= best);
        }
    }

    #[test]
    fn descent_reaches_local_minimum() {
        let mut rng = Rng::new(5);
        let j = DenseJ::sk(60, false, &mut rng);
        let mut s: Vec<i8> = (0..60).map(|_| if rng.coin() { 1 } else { -1 }).collect();
        let e = sk_descent(&j, &mut s, 1);
        let h = j.fields_int(&s, 1);
        assert!(h.iter().zip(&s).all(|(&h, &s)| s as i32 * h <= 0));
        assert_eq!(e, h.iter().zip(&s).map(|(&h, &s)| h as i64 * s as i64).sum::<i64>() / 2);
    }
}
