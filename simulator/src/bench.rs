//! Times the layer update on this CPU and compares with the FPGA estimates.

use crate::model::{self, Coupling, DenseJ, Spin};
use crate::{cordic, hw, par, rng::Rng};
use std::f64::consts::PI;
use std::hint::black_box;
use std::time::Instant;

/// Fastest of `reps` runs, after one untimed run (page faults, clock ramp-up).
fn best_of(reps: usize, mut f: impl FnMut()) -> f64 {
    f();
    (0..reps)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64()
        })
        .fold(f64::INFINITY, f64::min)
}

fn fmt_time(s: f64) -> String {
    if s < 1e-3 {
        format!("{:.1} µs", s * 1e6)
    } else if s < 1.0 {
        format!("{:.2} ms", s * 1e3)
    } else {
        format!("{s:.2} s")
    }
}

fn row(case: &str, n: usize, cpu: f64, fpga: f64, extra: &str) {
    let ratio = if fpga < cpu {
        format!("FPGA {:.1}x faster", cpu / fpga)
    } else {
        format!("CPU {:.1}x faster", fpga / cpu)
    };
    println!("{case:<40} {n:>7} {:>10} {:>10}  {ratio:<18} {extra}", fmt_time(cpu), fmt_time(fpga));
}

fn random_unit(rng: &mut Rng) -> Spin {
    let z = 2.0 * rng.uniform() - 1.0;
    let a = 2.0 * PI * rng.uniform();
    let r = (1.0 - z * z).sqrt();
    [r * a.cos(), r * a.sin(), z]
}

fn plural(t: usize) -> String {
    format!("{t} thread{}", if t > 1 { "s" } else { "" })
}

pub fn run(threads: usize) {
    println!("CPU: {}, {}\n", par::cpu_name(), plural(threads));
    println!("{:<40} {:>7} {:>10} {:>10}", "case", "N", "this CPU", "FPGA est.");

    // rank 1: the Cyclone 10 GX case, 49,152 spins and the 2 angle_lut layers,
    // including the initial and per-layer energies, as the hardware does
    let n = 49_152;
    let mut rng = Rng::new(1);
    let a: Vec<u16> = (0..n).map(|_| rng.range(1, 0x8000) as u16).collect();
    let init: Vec<Spin> = (0..n).map(|_| random_unit(&mut rng)).collect();
    let init_q: Vec<hw::HwSpin> = init.iter().map(|s| s.map(cordic::q15)).collect();
    let fpga = hw::clocks(n as u64, 2) as f64 / 300e6;
    let threads_list = if threads > 1 { vec![1, threads] } else { vec![1] };
    for &t in &threads_list {
        let cpu = best_of(5, || {
            let mut m = hw::Hw { spins: init_q.clone(), a: a.clone() };
            black_box(m.energy());
            for k in 0..2 {
                m.layer(hw::LUT_GAMMA[k], hw::LUT_BETA[k], t);
                black_box(m.energy());
            }
            black_box(&m.spins);
        });
        row(&format!("rank 1, bit-exact, 2 layers, {}", plural(t)), n, cpu, fpga, "vs Cyclone 10 GX");
    }
    let coupling = Coupling::Rank1 { a: a.iter().map(|&w| w as f64 / 32768.0).collect() };
    let c: Vec<f64> = hw::LUT_GAMMA.iter().map(|&g| PI * g as f64 / 2f64.powi(31)).collect();
    let phi: Vec<f64> = hw::LUT_BETA.iter().map(|&b| 2.0 * PI * b as f64 / 2f64.powi(32)).collect();
    for &t in &threads_list {
        let cpu = best_of(5, || {
            let mut spins = init.clone();
            for k in 0..=2 {
                let z: Vec<f64> = spins.iter().map(|s| s[2]).collect();
                let h = coupling.fields(&z, t);
                black_box(model::dot(&z, &h));
                if k < 2 {
                    model::rotate(&mut spins, &h, c[k], phi[k], t);
                }
            }
            black_box(&spins);
        });
        row(&format!("rank 1, float, 2 layers, {}", plural(t)), n, cpu, fpga, "vs Cyclone 10 GX");
    }

    // dense 8-bit J, one layer: the proposed Agilex 7 design with J on chip
    for &n in &[2_000usize, 8_700] {
        let coupling = Coupling::Dense(DenseJ::sk(n, false, &mut rng));
        let mut spins: Vec<Spin> = (0..n)
            .map(|_| {
                let z = 0.3 * (2.0 * rng.uniform() - 1.0);
                [(1.0 - z * z).sqrt(), 0.0, z]
            })
            .collect();
        let nf = n as f64;
        let fpga = nf * nf / (12_000.0 * 400e6) + nf / (4.0 * 400e6) + 100.0 / 400e6;
        for &t in &threads_list {
            let cpu = best_of(10, || {
                let z: Vec<f64> = spins.iter().map(|s| s[2]).collect();
                let h = coupling.fields(&z, t);
                black_box(model::dot(&z, &h));
                model::rotate(&mut spins, &h, 0.01, -0.01, t);
            });
            let gbs = nf * (nf - 1.0) / 2.0 / cpu / 1e9;
            row(
                &format!("dense 8-bit J, 1 layer, {}", plural(t)),
                n,
                cpu,
                fpga,
                &format!("vs Agilex 7 (J read at {gbs:.0} GB/s)"),
            );
        }
    }
}
