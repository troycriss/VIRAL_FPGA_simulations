//! One run: build a problem, step the spins layer by layer, compare with
//! classical baselines, and describe it all as JSON for the web page.

use crate::json::Value;
use crate::model::{self, Coupling, Spin};
use crate::problem::{self, Kind, Problem, Spec};
use crate::{baseline, cordic, hw, par, rng::Rng};
use std::collections::HashSet;
use std::f64::consts::PI;
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Engine {
    /// bit-exact model of VIRAL_n_spin (rank-1 problems only)
    Hw,
    Float,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Schedule {
    /// linear ramp as in MF-AOA: gamma_k = s_k dt, beta_k = (1 - s_k) dt
    Ramp,
    /// the hardware's angle_lut (4 layers, a test pattern)
    Lut,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Init {
    /// along +x (the mixer's ground state), z tilted by up to +-tilt
    /// (for hopfield: tilted toward the cue)
    X,
    /// random unit vectors
    Random,
}

#[derive(Clone, Copy)]
pub struct Params {
    pub problem: Kind,
    pub n: usize,
    pub seed: u64,
    pub layers: usize,
    pub dt: f64,
    /// default: hw for rank-1 problems, float otherwise
    pub engine: Option<Engine>,
    pub schedule: Schedule,
    pub init: Init,
    pub tilt: f64,
    pub gamma_lo: Option<u16>,
    pub degree: usize,
    pub alpha: f64,
    pub patterns: usize,
    pub noise: f64,
    /// spins kept per layer for the web page's sphere
    pub points: usize,
    pub threads: usize,
}

impl Default for Params {
    fn default() -> Params {
        Params {
            problem: Kind::Partition,
            n: 24,
            seed: 1,
            layers: 100,
            dt: 0.5,
            engine: None,
            schedule: Schedule::Ramp,
            init: Init::X,
            tilt: 0.05,
            gamma_lo: None,
            degree: 3,
            alpha: 3.0,
            patterns: 10,
            noise: 0.2,
            points: 1000,
            threads: par::num_threads(),
        }
    }
}

/// Default (n, layers, dt) for each problem, used when they aren't given.
pub fn defaults(kind: Kind) -> (usize, usize, f64) {
    match kind {
        Kind::Partition => (24, 100, 0.5),
        Kind::Sk | Kind::SkGauss => (300, 200, 0.3),
        Kind::MaxCut => (100, 200, 0.3),
        Kind::Chain => (100, 1000, 0.3),
        Kind::Hopfield => (400, 100, 0.3),
    }
}

pub const KEYS: &str = "problem, n, seed, layers, dt, engine, schedule, init, tilt, gamma_lo, \
                        degree, alpha, patterns, noise, points, threads";

pub fn int(k: &str, v: &str) -> Result<usize, String> {
    v.replace('_', "").parse().map_err(|_| format!("{k}: not a whole number: {v}"))
}

pub fn float(k: &str, v: &str) -> Result<f64, String> {
    v.parse().map_err(|_| format!("{k}: not a number: {v}"))
}

impl Params {
    /// Parses key-value pairs; `extra` keys are skipped (for the studies).
    pub fn from_pairs_with(pairs: &[(String, String)], extra: &[&str]) -> Result<Params, String> {
        let mut p = Params::default();
        let mut given = HashSet::new();
        for (k, v) in pairs {
            let v = v.as_str();
            given.insert(k.as_str());
            match k.as_str() {
                "problem" => p.problem = Kind::parse(v)?,
                "n" => p.n = int(k, v)?,
                "seed" => p.seed = int(k, v)? as u64,
                "layers" => p.layers = int(k, v)?,
                "dt" => p.dt = float(k, v)?,
                "engine" => {
                    p.engine = Some(match v {
                        "hw" => Engine::Hw,
                        "float" => Engine::Float,
                        _ => return Err(format!("engine: {v:?} (hw, float)")),
                    })
                }
                "schedule" => {
                    p.schedule = match v {
                        "ramp" => Schedule::Ramp,
                        "lut" => Schedule::Lut,
                        _ => return Err(format!("schedule: {v:?} (ramp, lut)")),
                    }
                }
                "init" => {
                    p.init = match v {
                        "x" => Init::X,
                        "random" => Init::Random,
                        _ => return Err(format!("init: {v:?} (x, random)")),
                    }
                }
                "tilt" => p.tilt = float(k, v)?,
                "gamma_lo" => {
                    p.gamma_lo = Some(u16::try_from(int(k, v)?).map_err(|_| "gamma_lo: 0 to 65535")?)
                }
                "degree" => p.degree = int(k, v)?,
                "alpha" => p.alpha = float(k, v)?,
                "patterns" => p.patterns = int(k, v)?,
                "noise" => p.noise = float(k, v)?,
                "points" => p.points = int(k, v)?,
                "threads" => p.threads = int(k, v)?.max(1),
                k if extra.contains(&k) => {}
                _ => return Err(format!("unknown parameter {k:?} (use {KEYS})")),
            }
        }
        let (n, layers, dt) = defaults(p.problem);
        if !given.contains("n") {
            p.n = n;
        }
        if !given.contains("layers") {
            p.layers = layers;
        }
        if !given.contains("dt") {
            p.dt = dt;
        }
        p.check()?;
        Ok(p)
    }

    pub fn from_pairs(pairs: &[(String, String)]) -> Result<Params, String> {
        Params::from_pairs_with(pairs, &[])
    }

    pub fn check(&self) -> Result<(), String> {
        let rank1 = self.problem.rank1();
        let max_n = self.problem.max_n();
        if self.n < 1 || self.n > max_n {
            return Err(format!("n: 1 to {max_n} for this problem"));
        }
        if self.layers < 1 || self.layers > 5000 {
            return Err("layers: 1 to 5000".into());
        }
        if self.schedule == Schedule::Lut && (!rank1 || self.layers > 4) {
            return Err("schedule=lut is the hardware's 4-entry angle_lut: partition only, layers <= 4".into());
        }
        if self.engine == Some(Engine::Hw) && !rank1 {
            return Err("engine=hw (VIRAL_n_spin) only has rank-1 coupling: use problem=partition".into());
        }
        if !(self.dt > 0.0 && self.dt <= 100.0) {
            return Err("dt: above 0, at most 100".into());
        }
        if !(0.0..=1.0).contains(&self.tilt) {
            return Err("tilt: 0 to 1".into());
        }
        if self.points > 20_000 {
            return Err("points: at most 20000".into());
        }
        Ok(())
    }

    pub fn engine(&self) -> Engine {
        self.engine.unwrap_or(if self.problem.rank1() { Engine::Hw } else { Engine::Float })
    }

    pub fn spec(&self) -> Spec {
        Spec {
            kind: self.problem,
            n: self.n,
            degree: self.degree,
            alpha: self.alpha,
            patterns: self.patterns,
            noise: self.noise,
        }
    }

    pub fn to_json(&self) -> Value {
        Value::obj()
            .set("problem", self.problem.name())
            .set("n", self.n)
            .set("seed", self.seed)
            .set("layers", self.layers)
            .set("dt", self.dt)
            .set("engine", if self.engine() == Engine::Hw { "hw" } else { "float" })
            .set("schedule", if self.schedule == Schedule::Ramp { "ramp" } else { "lut" })
            .set("init", if self.init == Init::X { "x" } else { "random" })
            .set("tilt", self.tilt)
            .set("gamma_lo", self.gamma_lo.map(|v| v as i64))
            .set("degree", self.degree)
            .set("alpha", self.alpha)
            .set("patterns", self.patterns)
            .set("noise", self.noise)
            .set("points", self.points)
            .set("threads", self.threads)
    }
}

/// Per-layer record of a run.
#[derive(Default)]
struct Trace {
    energy: Vec<f64>,
    rounded: Vec<f64>,
    metric_mf: Vec<f64>,
    metric_rounded: Vec<f64>,
    snaps: Vec<Value>,
    hists: Vec<Value>,
    /// time in the simulation itself (fields, energy, rotations)
    seconds: f64,
}

const HIST_BINS: usize = 40;

impl Trace {
    fn record(&mut self, prob: &Problem, z: &[f64], shown: &[usize], spin: impl Fn(usize) -> Spin, threads: usize) {
        let e_round = prob.rounded_energy(&problem::signs(z), threads);
        let (m_mf, m_round) = prob.metric(z, *self.energy.last().unwrap(), e_round);
        self.rounded.push(e_round);
        self.metric_mf.push(m_mf);
        self.metric_rounded.push(m_round);
        self.snaps.push(Value::rounded(shown.iter().flat_map(|&i| spin(i)), 4));
        let mut h = vec![0u64; HIST_BINS];
        for &v in z {
            let b = ((v + 1.0) / 2.0 * HIST_BINS as f64).floor().clamp(0.0, (HIST_BINS - 1) as f64);
            h[b as usize] += 1;
        }
        self.hists.push(h.into());
    }
}

fn trace_float(
    prob: &Problem,
    init: &[Spin],
    c: &[f64],
    phi: &[f64],
    shown: Option<&[usize]>,
    threads: usize,
) -> (Vec<Spin>, Trace) {
    let mut spins = init.to_vec();
    let mut t = Trace::default();
    for k in 0..=c.len() {
        let z: Vec<f64> = spins.iter().map(|s| s[2]).collect();
        let start = Instant::now();
        let h = prob.coupling.fields(&z, threads);
        t.energy.push(0.5 * model::dot(&z, &h));
        t.seconds += start.elapsed().as_secs_f64();
        if let Some(shown) = shown {
            t.record(prob, &z, shown, |i| spins[i], threads);
        }
        if k < c.len() {
            let start = Instant::now();
            model::rotate(&mut spins, &h, c[k], phi[k], threads);
            t.seconds += start.elapsed().as_secs_f64();
        }
    }
    (spins, t)
}

fn to_float(s: hw::HwSpin) -> Spin {
    [s[0] as f64 / 32768.0, s[1] as f64 / 32768.0, s[2] as f64 / 32768.0]
}

fn trace_hw(
    prob: &Problem,
    init: &[hw::HwSpin],
    gamma: &[i32],
    beta: &[i32],
    shown: Option<&[usize]>,
    threads: usize,
) -> (hw::Hw, Trace) {
    let mut m = hw::Hw { spins: init.to_vec(), a: prob.weights.clone().unwrap() };
    let mut t = Trace::default();
    for k in 0..=gamma.len() {
        let start = Instant::now();
        t.energy.push(m.energy() as f64 / 32768.0);
        t.seconds += start.elapsed().as_secs_f64();
        if let Some(shown) = shown {
            let z: Vec<f64> = m.spins.iter().map(|s| s[2] as f64 / 32768.0).collect();
            t.record(prob, &z, shown, |i| to_float(m.spins[i]), threads);
        }
        if k < gamma.len() {
            let start = Instant::now();
            m.layer(gamma[k], beta[k], threads);
            t.seconds += start.elapsed().as_secs_f64();
        }
    }
    (m, t)
}

/// Time the FPGA would take for the same run: (seconds, design, note, clocks).
fn fpga_estimate(prob: &Problem, layers: usize) -> (f64, &'static str, &'static str, Option<u64>) {
    let n = prob.n() as f64;
    let products = (layers + 1) as f64; // one field computation per layer, plus the final energy
    let latency = 100.0 / 400e6; // pipeline latency per layer
    let rotate = n / (4.0 * 400e6); // 4 rotation pipelines
    match prob.kind {
        Kind::Partition => {
            let clocks = hw::clocks(prob.n() as u64, layers as u64);
            let note = if layers > 4 {
                "angle_lut holds 4 layers, so this many needs a longer table (a few M20K blocks)"
            } else {
                ""
            };
            (
                clocks as f64 / 300e6,
                "VIRAL_n_spin on a Cyclone 10 GX at 300 MHz (design simulated; clock speed estimated)",
                note,
                Some(clocks),
            )
        }
        Kind::Sk | Kind::SkGauss | Kind::MaxCut => {
            let note = if prob.kind == Kind::MaxCut {
                "J stored dense; a sparse-J design would need only N x degree multiply-adds per layer"
            } else {
                "12,000 multipliers and 4 rotation pipelines at 400 MHz, 100 clocks of latency per layer"
            };
            if n <= 8700.0 {
                let t = products * (n * n / (12_000.0 * 400e6) + rotate + latency);
                (t, "dense-J design on an Agilex 7 AGM 039, J in block RAM (not built; estimate)", note, None)
            } else {
                let t = products * (n * n / 2.0 / 7e11 + latency);
                (t, "dense-J design on an Agilex 7 AGM 039, J in HBM2e (not built; estimate)", "J read from HBM at 700 GB/s every layer", None)
            }
        }
        Kind::Chain => (
            products * (n * n / (12_000.0 * 400e6) + rotate + latency),
            "Agilex 7 design computing J = 1/r^α from a table by distance (not built; estimate)",
            "no J memory, N² multiply-adds per layer; the CPU uses an FFT instead (N log N)",
            None,
        ),
        Kind::Hopfield => {
            let p = prob.patterns().map_or(1, |x| x.len()) as f64;
            (
                products * (2.0 * n * p / (12_000.0 * 400e6) + rotate + latency),
                "rank-P version of VIRAL_n_spin on an Agilex 7 (not built; estimate)",
                "P sums per layer instead of 1: 2 N P multiply-adds per layer",
                None,
            )
        }
    }
}

/// Rough operation count of a run (setup, layers, local search), so studies
/// can report progress by work rather than by number of runs.
pub fn work_estimate(p: &Params) -> f64 {
    let n = p.n as f64;
    let l = (p.layers + 1) as f64;
    match p.problem {
        Kind::Partition => 50.0 * n * l,
        Kind::Sk | Kind::SkGauss | Kind::MaxCut => n * n * (l / 2.0 + 4.0),
        Kind::Chain => 40.0 * n * (2.0 * n).log2() * l + if p.n <= 20_000 { n * n } else { 0.0 },
        Kind::Hopfield => {
            let pats = p.patterns as f64;
            2.0 * n * pats * l + n * n * pats / 4.0
        }
    }
}

/// Everything a run needs before the layers start.
struct Prepared {
    prob: Problem,
    field_scale: f64,
    /// rotation per layer: about z by c[k] * h_i, then about x by phi[k]
    c: Vec<f64>,
    phi: Vec<f64>,
    hw_gamma: Vec<i32>,
    hw_beta: Vec<i32>,
    init: Vec<Spin>,
    init_q: Vec<hw::HwSpin>,
}

fn prepare(p: &Params) -> Result<Prepared, String> {
    let engine = p.engine();
    let mut rng = Rng::new(p.seed);
    let prob = Problem::generate(&p.spec(), &mut rng)?;
    let n = prob.n();
    let field_scale = prob.coupling.field_scale().max(1e-300);
    let layers = p.layers;

    let (mut c, mut phi) = (Vec::new(), Vec::new());
    let (mut hw_gamma, mut hw_beta) = (Vec::new(), Vec::new());
    match p.schedule {
        Schedule::Ramp => {
            // exp(-i gamma H) turns spin i by 2 gamma h_i about z; the mixer
            // exp(+i beta sum X) turns every spin by -2 beta about x
            for k in 0..layers {
                let s = (k as f64 + 0.5) / layers as f64;
                c.push(2.0 * s * p.dt / field_scale);
                phi.push(-2.0 * (1.0 - s) * p.dt);
            }
            if engine == Engine::Hw {
                for k in 0..layers {
                    let g = (c[k] / PI * 2f64.powi(31)).round();
                    if g.abs() >= 2f64.powi(31) {
                        return Err(format!(
                            "gamma in layer {} is too big for the hardware's Q1.31 gamma: use a smaller dt",
                            k + 1
                        ));
                    }
                    hw_gamma.push(g as i32);
                    hw_beta.push((phi[k] / (2.0 * PI) * 2f64.powi(32)).round() as i64 as u32 as i32);
                }
            }
        }
        Schedule::Lut => {
            for k in 0..layers {
                let mut g = hw::LUT_GAMMA[k];
                if let Some(lo) = p.gamma_lo {
                    g = (g & !0xFFFF) | lo as i32;
                }
                hw_gamma.push(g);
                hw_beta.push(hw::LUT_BETA[k]);
            }
        }
    }
    if !hw_gamma.is_empty() {
        // the float model then turns by exactly the hardware's angles
        c = hw_gamma.iter().map(|&g| PI * g as f64 / 2f64.powi(31)).collect();
        phi = hw_beta.iter().map(|&b| 2.0 * PI * b as f64 / 2f64.powi(32)).collect();
    }

    let init: Vec<Spin> = (0..n)
        .map(|i| match p.init {
            Init::X => {
                let z = match &prob.cue {
                    Some(cue) => p.tilt * cue[i] as f64,
                    None => p.tilt * (2.0 * rng.uniform() - 1.0),
                };
                [(1.0 - z * z).sqrt(), 0.0, z]
            }
            Init::Random => {
                let z = 2.0 * rng.uniform() - 1.0;
                let a = 2.0 * PI * rng.uniform();
                let r = (1.0 - z * z).sqrt();
                [r * a.cos(), r * a.sin(), z]
            }
        })
        .collect();
    let init_q = init.iter().map(|s| s.map(cordic::q15)).collect();
    Ok(Prepared { prob, field_scale, c, phi, hw_gamma, hw_beta, init, init_q })
}

/// Local search from spins s, if it's affordable for this problem's size.
/// Returns the new spins and their energy.
fn local_search(prob: &Problem, s: &[i8], threads: usize) -> Option<(Vec<i8>, f64)> {
    let n = prob.n();
    let mut s = s.to_vec();
    let e = match &prob.coupling {
        Coupling::Rank1 { .. } => {
            if n > 100_000 {
                return None;
            }
            baseline::partition_local(prob.weights.as_ref().unwrap(), &mut s);
            prob.rounded_energy(&s, threads)
        }
        Coupling::Dense(d) => baseline::sk_descent(d, &mut s, threads) as f64 * d.scale,
        Coupling::Chain { .. } => {
            if n > 20_000 {
                return None;
            }
            baseline::descent(&prob.coupling, &mut s, threads)
        }
        Coupling::LowRank { xi, .. } => {
            if (n * n * xi.len()) as f64 > 4e9 {
                return None;
            }
            baseline::descent(&prob.coupling, &mut s, threads)
        }
    };
    Some((s, e))
}

/// Local search from a random start (the classical method without mean field).
fn random_local_search(prob: &Problem, seed: u64, threads: usize) -> Option<(Vec<i8>, f64)> {
    let mut rng = Rng::new(seed ^ 0x5EED_5EED);
    let s: Vec<i8> = (0..prob.n()).map(|_| if rng.coin() { 1 } else { -1 }).collect();
    local_search(prob, &s, threads)
}

/// Result of a run without the per-layer record, for the studies.
pub struct Quick {
    pub energy: f64,
    pub rounded: f64,
    pub local: Option<f64>,
    pub metric: f64,
    pub metric_local: Option<f64>,
    pub seconds: f64,
    pub fpga_seconds: f64,
}

pub fn quick(p: &Params) -> Result<Quick, String> {
    let threads = p.threads.max(1);
    let pr = prepare(p)?;
    let (z, t) = match p.engine() {
        Engine::Hw => {
            let (m, t) = trace_hw(&pr.prob, &pr.init_q, &pr.hw_gamma, &pr.hw_beta, None, threads);
            (m.spins.iter().map(|s| s[2] as f64 / 32768.0).collect::<Vec<_>>(), t)
        }
        Engine::Float => {
            let (s, t) = trace_float(&pr.prob, &pr.init, &pr.c, &pr.phi, None, threads);
            (s.iter().map(|s| s[2]).collect(), t)
        }
    };
    let s = problem::signs(&z);
    let energy = *t.energy.last().unwrap();
    let rounded = pr.prob.rounded_energy(&s, threads);
    let (_, metric) = pr.prob.metric(&z, energy, rounded);
    let local = local_search(&pr.prob, &s, threads);
    let metric_local = local.as_ref().map(|(sl, e)| {
        let zl: Vec<f64> = sl.iter().map(|&v| v as f64).collect();
        pr.prob.metric(&zl, *e, *e).1
    });
    Ok(Quick {
        energy,
        rounded,
        local: local.map(|l| l.1),
        metric,
        metric_local,
        seconds: t.seconds,
        fpga_seconds: fpga_estimate(&pr.prob, p.layers).0,
    })
}

pub struct Output {
    pub json: Value,
    pub summary: String,
}

/// Show at most this many numbers (weights, assignments) in the JSON.
const MAX_LISTED: usize = 20_000;

/// One reference line: a named energy, and the same answer in the problem's
/// metric (when it has a meaning there).
struct Reference {
    name: String,
    energy: f64,
    metric: Option<f64>,
}

fn reference(name: &str, energy: f64, metric: Option<f64>) -> Reference {
    Reference { name: name.to_string(), energy, metric }
}

fn signs_json(s: &[i8]) -> Value {
    s.iter().map(|&v| v as i64).collect::<Vec<_>>().into()
}

pub fn run(p: &Params) -> Result<Output, String> {
    let threads = p.threads.max(1);
    let engine = p.engine();
    let pr = prepare(p)?;
    let prob = &pr.prob;
    let n = prob.n();
    let nf = n as f64;
    let layers = p.layers;
    let shown: Vec<usize> = if n <= p.points {
        (0..n).collect()
    } else {
        (0..p.points).map(|k| k * n / p.points).collect()
    };

    let mut compare = Value::Null;
    let mut hardware = Value::Null;
    let mut summary = Vec::new();
    let (final_z, trace) = match engine {
        Engine::Hw => {
            let (m, t) = trace_hw(prob, &pr.init_q, &pr.hw_gamma, &pr.hw_beta, Some(&shown), threads);
            // the same run in floating point, from the same rounded start
            let init_f: Vec<Spin> = pr.init_q.iter().map(|&s| to_float(s)).collect();
            let (fs, ft) = trace_float(prob, &init_f, &pr.c, &pr.phi, None, threads);
            let max_diff = m
                .spins
                .iter()
                .zip(&fs)
                .flat_map(|(a, b)| {
                    let a = to_float(*a);
                    (0..3).map(move |k| (a[k] - b[k]).abs())
                })
                .fold(0.0, f64::max);
            summary.push(format!("  bit-exact FPGA model vs floating point: largest spin difference {max_diff:.4}"));
            compare = Value::obj().set("engine", "float").set("energy", ft.energy).set("max_spin_diff", max_diff);
            hardware = Value::obj()
                .set("energy_q15", m.energy())
                .set("led", format!("{:08b}", m.led(8)))
                .set("gamma", pr.hw_gamma.clone())
                .set("beta", pr.hw_beta.clone());
            (m.spins.iter().map(|s| s[2] as f64 / 32768.0).collect::<Vec<_>>(), t)
        }
        Engine::Float => {
            let (s, t) = trace_float(prob, &pr.init, &pr.c, &pr.phi, Some(&shown), threads);
            (s.iter().map(|s| s[2]).collect(), t)
        }
    };

    let s = problem::signs(&final_z);
    let final_energy = *trace.energy.last().unwrap();
    let final_rounded = *trace.rounded.last().unwrap();
    let local = local_search(prob, &s, threads);
    let metric_of = |sv: &[i8], e: f64| {
        let z: Vec<f64> = sv.iter().map(|&v| v as f64).collect();
        prob.metric(&z, e, e).1
    };
    let mut refs: Vec<Reference> = Vec::new();
    let mut result = Value::obj().set("energy", final_energy).set("rounded_energy", final_rounded);
    let mut problem_json = Value::obj()
        .set("kind", prob.kind.name())
        .set("title", prob.kind.title())
        .set("n", n)
        .set("field_scale", pr.field_scale);
    let mut residuals = Value::Null;

    match prob.kind {
        Kind::Partition => {
            let Coupling::Rank1 { a } = &prob.coupling else { unreachable!() };
            let w = prob.weights.as_ref().unwrap();
            let a2: f64 = a.iter().map(|a| a * a).sum();
            let q = |v: u64| v as f64 / 32768.0;
            let e_of = |r: u64| (q(r) * q(r) - a2) / 2.0;
            let total: u64 = w.iter().map(|&w| w as u64).sum();
            let sum_a: u64 = w.iter().zip(&s).filter(|(_, &s)| s > 0).map(|(&w, _)| w as u64).sum();
            let residual = sum_a.abs_diff(total - sum_a);
            let kk = baseline::karmarkar_karp(w);
            let greedy = baseline::greedy(w);
            let opt = baseline::optimum(w, 300_000_000);
            if let Some(o) = opt {
                refs.push(reference("optimum (exact)", e_of(o), Some(q(o))));
            }
            refs.push(reference("Karmarkar–Karp", e_of(kk), Some(q(kk))));
            refs.push(reference("greedy", e_of(greedy), Some(q(greedy))));
            let local_res = local.as_ref().map(|(sl, _)| {
                w.iter().zip(sl).map(|(&w, &s)| w as i64 * s as i64).sum::<i64>().unsigned_abs()
            });
            if let Some(r) = local_res {
                refs.push(reference("mean field + local search", e_of(r), Some(q(r))));
            }
            let mut part = Value::obj()
                .set("total", q(total))
                .set("sum_a", q(sum_a))
                .set("sum_b", q(total - sum_a))
                .set("residual", q(residual))
                .set("residual_local", local_res.map(q))
                .set("karmarkar_karp", q(kk))
                .set("greedy", q(greedy))
                .set("optimum", opt.map(q));
            if n <= MAX_LISTED {
                part = part.set("assignment", signs_json(&s));
                if let Some((sl, _)) = &local {
                    part = part.set("assignment_local", signs_json(sl));
                }
                problem_json = problem_json.set("weights", Value::rounded(a.iter().copied(), 5));
            }
            result = result.set("partition", part);
            problem_json = problem_json.set("coupling", "rank 1: J_ij = a_i a_j");
            residuals = Value::obj()
                .set("rounded", trace.metric_rounded.clone())
                .set("mean_field", trace.metric_mf.clone());
            summary.push(format!(
                "  residual |sum A - sum B|: mean field {:.6}, + local search {}",
                q(residual),
                local_res.map_or("(too big)".into(), |r| format!("{:.6}", q(r)))
            ));
            summary.push(format!(
                "  for comparison: Karmarkar-Karp {:.6}, greedy {:.6}, optimum {}",
                q(kk),
                q(greedy),
                opt.map_or("(too big to compute)".into(), |o| format!("{:.6}", q(o)))
            ));
        }
        Kind::Sk | Kind::SkGauss => {
            let Coupling::Dense(d) = &prob.coupling else { unreachable!() };
            let parisi = -0.763166;
            refs.push(reference("ground state for N → ∞ (Parisi)", parisi * nf, Some(parisi)));
            let rl = random_local_search(prob, p.seed, threads);
            if let Some((_, e)) = &rl {
                refs.push(reference("random start + local search", *e, Some(e / nf)));
            }
            if let Some((_, e)) = &local {
                refs.push(reference("mean field + local search", *e, Some(e / nf)));
            }
            result = result.set(
                "sk",
                Value::obj()
                    .set("energy_per_spin", final_rounded / nf)
                    .set("local_per_spin", local.as_ref().map(|l| l.1 / nf))
                    .set("random_local_per_spin", rl.as_ref().map(|l| l.1 / nf))
                    .set("parisi_per_spin", parisi),
            );
            if n <= MAX_LISTED {
                result = result.set("assignment", signs_json(&s));
            }
            problem_json = problem_json.set("coupling", "dense, 8-bit");
            if n <= 200 {
                problem_json = problem_json.set(
                    "j",
                    Value::obj().set("scale", d.scale).set("q", d.q.iter().map(|&q| q as i64).collect::<Vec<_>>()),
                );
            }
            summary.push(format!(
                "  energy per spin: mean field {:.4}, + local search {:.4}; random start + local search {:.4}; N -> infinity ground state {parisi}",
                final_rounded / nf,
                local.as_ref().map_or(f64::NAN, |l| l.1 / nf),
                rl.as_ref().map_or(f64::NAN, |l| l.1 / nf)
            ));
        }
        Kind::MaxCut => {
            let edges = prob.edges.as_ref().unwrap();
            let m = edges.len() as f64;
            let degree = 2 * edges.len() / n;
            let frac = |e: f64| (m - e) / (2.0 * m);
            if degree == 3 {
                // QAOA on large random 3-regular graphs (locally tree-like):
                // fraction of edges cut at depth p = 1 (Farhi et al. 2014) and p = 2
                // (Wurtz and Love 2021)
                for (name, f) in [("QAOA p = 1: 69.24% of edges", 0.6924), ("QAOA p = 2: 75.59% of edges", 0.7559)] {
                    refs.push(reference(name, m * (1.0 - 2.0 * f), Some(f)));
                }
            }
            let rl = random_local_search(prob, p.seed, threads);
            if let Some((_, e)) = &rl {
                refs.push(reference("random start + local search", *e, Some(frac(*e))));
            }
            if let Some((_, e)) = &local {
                refs.push(reference("mean field + local search", *e, Some(frac(*e))));
            }
            result = result.set(
                "maxcut",
                Value::obj()
                    .set("edges", edges.len())
                    .set("degree", degree)
                    .set("cut", (m - final_rounded) / 2.0)
                    .set("fraction", frac(final_rounded))
                    .set("local_fraction", local.as_ref().map(|l| frac(l.1)))
                    .set("random_local_fraction", rl.as_ref().map(|l| frac(l.1))),
            );
            if n <= MAX_LISTED {
                result = result.set("assignment", signs_json(&s));
                if let Some((sl, _)) = &local {
                    result = result.set("assignment_local", signs_json(sl));
                }
            }
            problem_json = problem_json.set("coupling", "J = 1 on each edge (stored dense, 8-bit)");
            if n <= 500 {
                problem_json = problem_json.set(
                    "edges",
                    Value::Arr(edges.iter().map(|&(a, b)| Value::from(vec![a as i64, b as i64])).collect()),
                );
            }
            summary.push(format!(
                "  edges cut: mean field {:.4}, + local search {:.4}; random start + local search {:.4}{}",
                frac(final_rounded),
                local.as_ref().map_or(f64::NAN, |l| frac(l.1)),
                rl.as_ref().map_or(f64::NAN, |l| frac(l.1)),
                if degree == 3 { "; QAOA p=1 0.6924, p=2 0.7559" } else { "" }
            ));
        }
        Kind::Chain => {
            let Coupling::Chain { f, .. } = &prob.coupling else { unreachable!() };
            // alternating spins: sum_{d} (n - d) f[d] (-1)^d
            let neel: f64 = (1..n).map(|d| (n - d) as f64 * f[d] * if d % 2 == 1 { -1.0 } else { 1.0 }).sum();
            refs.push(reference("alternating ↑↓↑↓ (the ground state)", neel, Some(1.0)));
            // overall alternating order |sum (-1)^i s_i| / N: domains of the two
            // alternating patterns cancel, so a few defects lower it a lot
            let order = |sv: &[i8]| {
                sv.iter().enumerate().map(|(i, &v)| if i % 2 == 0 { v as f64 } else { -(v as f64) }).sum::<f64>().abs() / nf
            };
            let rl = random_local_search(prob, p.seed, threads);
            if let Some((sl, e)) = &rl {
                refs.push(reference("random start + local search", *e, Some(metric_of(sl, *e))));
            }
            if let Some((sl, e)) = &local {
                refs.push(reference("mean field + local search", *e, Some(metric_of(sl, *e))));
            }
            let walls = |sv: &[i8]| sv.windows(2).filter(|w| w[0] == w[1]).count();
            result = result.set(
                "chain",
                Value::obj()
                    .set("alpha", prob.alpha)
                    .set("alternating", metric_of(&s, final_rounded))
                    .set("order", order(&s))
                    .set("defects", walls(&s))
                    .set("order_local", local.as_ref().map(|(sl, _)| order(sl)))
                    .set("defects_local", local.as_ref().map(|(sl, _)| walls(sl)))
                    .set("ground_energy", neel),
            );
            if n <= MAX_LISTED {
                result = result.set("assignment", signs_json(&s));
                if let Some((sl, _)) = &local {
                    result = result.set("assignment_local", signs_json(sl));
                }
            }
            problem_json = problem_json
                .set("coupling", "J_ij = 1/|i − j|^α, computed by distance")
                .set("alpha", prob.alpha)
                .set("couplings", Value::rounded(f[..n.min(64)].iter().copied(), 6));
            summary.push(format!(
                "  defects (neighbours alike): mean field {}, + local search {}; energy {:.4} vs ground state {:.4}",
                walls(&s),
                local.as_ref().map_or("(too big)".into(), |(sl, _)| walls(sl).to_string()),
                final_rounded,
                neel
            ));
        }
        Kind::Hopfield => {
            let xi = prob.patterns().unwrap();
            let cue = prob.cue.as_ref().unwrap();
            let overlap = |x: &[i8], sv: &[i8]| x.iter().zip(sv).map(|(&a, &b)| (a * b) as f64).sum::<f64>() / nf;
            let overlaps: Vec<f64> = xi.iter().map(|x| overlap(x, &s)).collect();
            let best = (0..xi.len()).max_by(|&a, &b| overlaps[a].abs().total_cmp(&overlaps[b].abs())).unwrap();
            // the classic Hopfield recall: zero-temperature descent from the cue
            let recall = local_search(prob, cue, threads);
            refs.push(reference("stored pattern 1", prob.rounded_energy(&xi[0], threads), Some(1.0)));
            if let Some((sl, e)) = &recall {
                refs.push(reference("cue + local search (classic recall)", *e, Some(metric_of(sl, *e))));
            }
            if let Some((sl, e)) = &local {
                refs.push(reference("mean field + local search", *e, Some(metric_of(sl, *e))));
            }
            let pf = xi.len() as f64;
            result = result.set(
                "hopfield",
                Value::obj()
                    .set("patterns", xi.len())
                    .set("load", pf / nf)
                    .set("capacity", 0.138)
                    // below this load the stored patterns are the lowest-energy
                    // states (Amit, Gutfreund, Sompolinsky 1985), so annealing finds them
                    .set("global_limit", 0.051)
                    .set("overlap", overlaps[0])
                    .set("overlap_cue", overlap(&xi[0], cue))
                    .set("overlap_local", local.as_ref().map(|(sl, _)| overlap(&xi[0], sl)))
                    .set("overlap_recall", recall.as_ref().map(|(sl, _)| overlap(&xi[0], sl)))
                    .set("best_pattern", best + 1)
                    .set("best_overlap", overlaps[best])
                    .set("overlaps", Value::rounded(overlaps.iter().take(100).copied(), 4)),
            );
            if n <= MAX_LISTED {
                result = result.set("assignment", signs_json(&s));
                if let Some((sl, _)) = &local {
                    result = result.set("assignment_local", signs_json(sl));
                }
                problem_json = problem_json.set("cue", signs_json(cue));
            }
            let side = (nf.sqrt().round()) as usize;
            problem_json = problem_json
                .set("coupling", "J = −(1/N) Σ_μ ξ^μ ξ^μᵀ, rank P")
                .set("noise", p.noise)
                .set("shape", if side * side == n { Value::from(vec![side as i64, side as i64]) } else { Value::Null });
            if xi.len() * n <= 100_000 {
                problem_json = problem_json.set("patterns", Value::Arr(xi.iter().map(|x| signs_json(x)).collect()));
            }
            summary.push(format!(
                "  overlap with pattern 1: cue {:.3}, mean field {:.3}, + local search {}, classic recall {}; {} patterns, load {:.3} (capacity about 0.138)",
                overlap(&xi[0], cue),
                overlaps[0],
                local.as_ref().map_or("(too big)".into(), |(sl, _)| format!("{:.3}", overlap(&xi[0], sl))),
                recall.as_ref().map_or("(too big)".into(), |(sl, _)| format!("{:.3}", overlap(&xi[0], sl))),
                xi.len(),
                pf / nf
            ));
        }
    }

    let (fpga_seconds, design, note, clocks) = fpga_estimate(prob, layers);
    summary.insert(
        0,
        format!(
            "{}, N = {n}, seed {}, {layers} layers ({}), engine {}",
            prob.kind.title(),
            p.seed,
            if p.schedule == Schedule::Ramp { format!("ramp, dt {}", p.dt) } else { "angle_lut".into() },
            if engine == Engine::Hw { "hw (bit-exact VIRAL_n_spin model)" } else { "float" },
        ),
    );
    summary.insert(
        1,
        format!("  energy: start {:.6}, end {:.6}; rounded to +-1: {:.6}", trace.energy[0], final_energy, final_rounded),
    );
    summary.push(format!(
        "  time: {:.3} ms on {} ({threads} threads); FPGA estimate {:.3} ms",
        trace.seconds * 1e3,
        par::cpu_name(),
        fpga_seconds * 1e3
    ));

    let schedule = Value::obj()
        .set("gamma", pr.c.iter().map(|&c| c * pr.field_scale / 2.0).collect::<Vec<_>>())
        .set("beta", pr.phi.iter().map(|&p| -p / 2.0).collect::<Vec<_>>())
        .set("z_angle_per_field", pr.c.clone())
        .set("x_angle", pr.phi.clone());

    let info = prob.metric_info();
    let metrics = Value::obj()
        .set("title", info.title)
        .set("log", info.log)
        .set("mean_field", Value::obj().set("name", info.mean_field).set("values", trace.metric_mf.clone()))
        .set("rounded", Value::obj().set("name", info.rounded).set("values", trace.metric_rounded.clone()))
        .set(
            "refs",
            Value::Arr(
                refs.iter()
                    .filter_map(|r| r.metric.map(|m| Value::obj().set("name", r.name.as_str()).set("value", m)))
                    .collect(),
            ),
        );

    let layers_json: Vec<Value> = trace
        .snaps
        .into_iter()
        .zip(trace.hists)
        .enumerate()
        .map(|(k, (spins, hist))| {
            Value::obj()
                .set("energy", trace.energy[k])
                .set("rounded_energy", trace.rounded[k])
                .set("spins", spins)
                .set("z_hist", hist)
        })
        .collect();

    let summary = summary.join("\n");
    let json = Value::obj()
        .set("version", 2i64)
        .set("params", p.to_json())
        .set("problem", problem_json)
        .set("schedule", schedule)
        .set("shown", shown)
        .set("layers", Value::Arr(layers_json))
        .set("metrics", metrics)
        .set("residuals", residuals)
        .set("compare", compare)
        .set("hardware", hardware)
        .set("result", result)
        .set(
            "baselines",
            Value::Arr(
                refs.iter()
                    .map(|r| Value::obj().set("name", r.name.as_str()).set("energy", r.energy).set("metric", r.metric))
                    .collect(),
            ),
        )
        .set(
            "timing",
            Value::obj().set("cpu_seconds", trace.seconds).set("threads", threads).set("cpu", par::cpu_name()).set(
                "fpga",
                Value::obj().set("seconds", fpga_seconds).set("clocks", clocks).set("design", design).set("note", note),
            ),
        )
        .set("summary", summary.clone());
    Ok(Output { json, summary })
}
