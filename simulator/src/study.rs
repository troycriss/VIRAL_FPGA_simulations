//! Studies: many runs, with the results as curves.
//! - scaling: SK energy per spin against N, extrapolated to N -> infinity
//! - sweep: the final result over a grid of layers and dt (schedule tuning)
//! - squeezing: one-axis twisting by the discrete truncated Wigner method
//!   (many classical runs from random starts), against the exact quantum result

use crate::json::Value;
use crate::model::sincos;
use crate::problem::Kind;
use crate::rng::Rng;
use crate::run::{self, float, int, Params, Quick};
use crate::par;
use std::f64::consts::LN_2;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Study {
    Scaling,
    Sweep,
    Squeezing,
}

pub struct StudyParams {
    pub study: Study,
    pub base: Params,
    pub sizes: Vec<usize>,
    pub seeds: usize,
    pub layers_list: Vec<usize>,
    pub dt_list: Vec<f64>,
    pub trajectories: usize,
    pub times: usize,
    pub t_max: Option<f64>,
}

const STUDY_KEYS: [&str; 8] = ["study", "sizes", "seeds", "layers_list", "dt_list", "trajectories", "times", "t_max"];

fn list<T>(k: &str, v: &str, parse: impl Fn(&str, &str) -> Result<T, String>) -> Result<Vec<T>, String> {
    v.split(',').filter(|s| !s.trim().is_empty()).map(|s| parse(k, s.trim())).collect()
}

impl StudyParams {
    pub fn from_pairs(pairs: &[(String, String)]) -> Result<StudyParams, String> {
        let get = |k: &str| pairs.iter().find(|(key, _)| key == k).map(|(_, v)| v.as_str());
        let study = match get("study").unwrap_or("") {
            "scaling" => Study::Scaling,
            "sweep" => Study::Sweep,
            "squeezing" => Study::Squeezing,
            s => return Err(format!("study: {s:?} (scaling, sweep, squeezing)")),
        };
        // the scaling study is about the SK model; the others default to it too
        let mut pairs = pairs.to_vec();
        if get("problem").is_none() && study != Study::Squeezing {
            pairs.push(("problem".into(), "sk".into()));
        }
        if study == Study::Squeezing && get("n").is_none() {
            pairs.push(("n".into(), "1000".into()));
        }
        let base = Params::from_pairs_with(&pairs, &STUDY_KEYS)?;
        let mut sp = StudyParams {
            study,
            base,
            sizes: vec![250, 500, 1000, 2000, 4000, 8000],
            seeds: if study == Study::Scaling { 10 } else { 5 },
            layers_list: vec![25, 50, 100, 200, 400, 800],
            dt_list: vec![0.1, 0.2, 0.3, 0.5, 0.8, 1.2],
            trajectories: 2000,
            times: 60,
            t_max: None,
        };
        for (k, v) in &pairs {
            match k.as_str() {
                "sizes" => sp.sizes = list(k, v, int)?,
                "seeds" => sp.seeds = int(k, v)?,
                "layers_list" => sp.layers_list = list(k, v, int)?,
                "dt_list" => sp.dt_list = list(k, v, float)?,
                "trajectories" => sp.trajectories = int(k, v)?,
                "times" => sp.times = int(k, v)?,
                "t_max" => sp.t_max = Some(float(k, v)?),
                _ => {}
            }
        }
        sp.check()?;
        Ok(sp)
    }

    fn check(&self) -> Result<(), String> {
        let kind = self.base.problem;
        match self.study {
            Study::Scaling => {
                if kind != Kind::Sk && kind != Kind::SkGauss {
                    return Err("the scaling study is for problem=sk or sk-gauss".into());
                }
                if self.sizes.len() < 2 || self.sizes.len() > 20 {
                    return Err("sizes: 2 to 20 sizes, e.g. 250,500,1000".into());
                }
                if let Some(&n) = self.sizes.iter().find(|&&n| n < 10 || n > kind.max_n()) {
                    return Err(format!("sizes: {n} is outside 10 to {}", kind.max_n()));
                }
            }
            Study::Sweep => {
                if self.layers_list.is_empty() || self.dt_list.is_empty() || self.layers_list.len() * self.dt_list.len() > 400 {
                    return Err("layers_list and dt_list: 1 to 400 combinations".into());
                }
                for &l in &self.layers_list {
                    Params { layers: l, ..self.base }.check()?;
                }
                for &dt in &self.dt_list {
                    Params { dt, ..self.base }.check()?;
                }
            }
            Study::Squeezing => {
                if self.base.n < 2 || self.base.n > 100_000 {
                    return Err("n: 2 to 100000 for the squeezing study".into());
                }
                if self.trajectories < 10 || self.trajectories > 10_000_000 {
                    return Err("trajectories: 10 to 10000000".into());
                }
                if self.times < 2 || self.times > 1000 {
                    return Err("times: 2 to 1000".into());
                }
            }
        }
        if self.seeds < 1 || self.seeds > 1000 {
            return Err("seeds: 1 to 1000".into());
        }
        Ok(())
    }
}

/// Called with (fraction done, message, partial result).
pub type Progress<'a> = &'a (dyn Fn(f64, &str, Option<Value>) + Sync);

/// Sweep results per (layers, dt): mean metric, mean metric after local
/// search, mean energy per spin; None while still running.
type Grid = Vec<Vec<Option<(f64, f64, f64)>>>;

/// Runs a study; returns its JSON and a text summary.
pub fn run(sp: &StudyParams, progress: Progress, cancel: &AtomicBool) -> Result<(Value, String), String> {
    match sp.study {
        Study::Scaling => scaling(sp, progress, cancel),
        Study::Sweep => sweep(sp, progress, cancel),
        Study::Squeezing => squeezing(sp, progress, cancel),
    }
}

/// Runs `seeds` quick runs of p (seeds base.seed, base.seed + 1, ...).
/// Small problems run in parallel, one thread each; large ones one at a time
/// on all threads.
fn runs(p: &Params, seeds: usize, cancel: &AtomicBool, done: impl Fn() + Sync) -> Result<Vec<Quick>, String> {
    let threads = p.threads.max(1);
    let (outer, inner) = if p.n < 2000 { (threads, 1) } else { (1, threads) };
    let out: Mutex<Vec<(usize, Result<Quick, String>)>> = Mutex::new(Vec::new());
    par::for_items(seeds, outer, || (), |k, _| {
        if cancel.load(Ordering::Relaxed) {
            return;
        }
        let q = Params { seed: p.seed + k as u64, threads: inner, ..*p };
        let r = run::quick(&q);
        out.lock().unwrap().push((k, r));
        done();
    });
    if cancel.load(Ordering::Relaxed) {
        return Err("cancelled".into());
    }
    let mut out = out.into_inner().unwrap();
    out.sort_by_key(|(k, _)| *k);
    out.into_iter().map(|(_, r)| r).collect()
}

fn mean_std(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = if v.len() > 1 { v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (n - 1.0) } else { 0.0 };
    (m, var.sqrt())
}

/// Weighted least squares y = a + b x; returns (a, error of a, b).
fn fit_line(x: &[f64], y: &[f64], w: &[f64]) -> (f64, f64, f64) {
    let (mut s, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for ((&x, &y), &w) in x.iter().zip(y).zip(w) {
        s += w;
        sx += w * x;
        sy += w * y;
        sxx += w * x * x;
        sxy += w * x * y;
    }
    let det = s * sxx - sx * sx;
    let b = (s * sxy - sx * sy) / det;
    let a = (sy - b * sx) / s;
    (a, (sxx / det).sqrt(), b)
}

fn scaling(sp: &StudyParams, progress: Progress, cancel: &AtomicBool) -> Result<(Value, String), String> {
    let start = Instant::now();
    let total = sp.sizes.len() * sp.seeds;
    let done = AtomicUsize::new(0);
    // progress by estimated work, since big sizes take far longer
    let work = |n: usize| run::work_estimate(&Params { n, ..sp.base });
    let total_work: f64 = sp.sizes.iter().map(|&n| work(n) * sp.seeds as f64).sum();
    let done_work = Mutex::new(0.0);
    let parisi = -0.763166;
    let mut points = Vec::new();
    let (mut cpu, mut fpga) = (0.0, 0.0);
    let build = |points: &Vec<(usize, Vec<Quick>)>, cpu: f64, fpga: f64, wall: f64| {
        let mut rows = Vec::new();
        let (mut xs, mut ys, mut ws, mut yl, mut wl) = (vec![], vec![], vec![], vec![], vec![]);
        for (n, qs) in points {
            let nf = *n as f64;
            let e: Vec<f64> = qs.iter().map(|q| q.rounded / nf).collect();
            let l: Vec<f64> = qs.iter().filter_map(|q| q.local.map(|v| v / nf)).collect();
            let (m, s) = mean_std(&e);
            let (lm, ls) = mean_std(&l);
            // weights 1/(standard error)^2, with a floor so one lucky size can't dominate
            let x = nf.powf(-2.0 / 3.0);
            xs.push(x);
            ys.push(m);
            ws.push(qs.len() as f64 / (s * s).max(1e-8));
            yl.push(lm);
            wl.push(qs.len() as f64 / (ls * ls).max(1e-8));
            rows.push(
                Value::obj()
                    .set("n", *n)
                    .set("x", x)
                    .set("count", qs.len())
                    .set("mean", m)
                    .set("std", s)
                    .set("local_mean", lm)
                    .set("local_std", ls)
                    .set("values", e)
                    .set("seconds", qs.iter().map(|q| q.seconds).sum::<f64>()),
            );
        }
        let fit = |y: &[f64], w: &[f64]| {
            if xs.len() < 2 {
                return Value::Null;
            }
            let (a, err, b) = fit_line(&xs, y, w);
            Value::obj().set("e_inf", a).set("e_inf_err", err).set("slope", b)
        };
        Value::obj()
            .set("study", "scaling")
            .set("params", sp.base.to_json())
            .set("sizes", sp.sizes.clone())
            .set("seeds", sp.seeds)
            .set("points", Value::Arr(rows))
            .set("fit", fit(&ys, &ws))
            .set("fit_local", fit(&yl, &wl))
            .set("parisi", parisi)
            .set("x_label", "N^(−2/3)")
            .set("cpu_seconds", cpu)
            .set("wall_seconds", wall)
            .set("fpga_seconds", fpga)
            .set("cpu", par::cpu_name())
    };
    for &n in &sp.sizes {
        let p = Params { n, ..sp.base };
        let w = work(n);
        let qs = runs(&p, sp.seeds, cancel, || {
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            let mut dw = done_work.lock().unwrap();
            *dw += w;
            progress(*dw / total_work, &format!("N = {n}: {d} of {total} runs"), None);
        })?;
        cpu += qs.iter().map(|q| q.seconds).sum::<f64>();
        fpga += qs.iter().map(|q| q.fpga_seconds).sum::<f64>();
        points.push((n, qs));
        let dw = *done_work.lock().unwrap();
        progress(
            dw / total_work,
            &format!("N = {n} done"),
            Some(build(&points, cpu, fpga, start.elapsed().as_secs_f64())),
        );
    }
    let mut text = format!(
        "{} energy per spin, {} layers, dt {}, {} instances per size:\n",
        sp.base.problem.title(),
        sp.base.layers,
        sp.base.dt,
        sp.seeds
    );
    for (n, qs) in &points {
        let e: Vec<f64> = qs.iter().map(|q| q.rounded / *n as f64).collect();
        let l: Vec<f64> = qs.iter().filter_map(|q| q.local.map(|v| v / *n as f64)).collect();
        let ((m, s), (lm, ls)) = (mean_std(&e), mean_std(&l));
        text += &format!("  N = {n:>6}: mean field {m:.4} ± {s:.4}, + local search {lm:.4} ± {ls:.4}\n");
    }
    let xs: Vec<f64> = points.iter().map(|(n, _)| (*n as f64).powf(-2.0 / 3.0)).collect();
    for (name, local) in [("mean field", false), ("mean field + local search", true)] {
        let mut ys = Vec::new();
        let mut ws = Vec::new();
        for (n, qs) in &points {
            let v: Vec<f64> = qs.iter().filter_map(|q| if local { q.local } else { Some(q.rounded) }).map(|e| e / *n as f64).collect();
            let (m, s) = mean_std(&v);
            ys.push(m);
            ws.push(v.len() as f64 / (s * s).max(1e-8));
        }
        let (a, err, _) = fit_line(&xs, &ys, &ws);
        text += &format!("  {name}, extrapolated to N -> infinity: {a:.4} ± {err:.4}\n");
    }
    text += &format!(
        "  ground state for N -> infinity (Parisi): {parisi}\n  time: {:.1} s of simulation on {}; FPGA estimate {:.2} s",
        cpu,
        par::cpu_name(),
        fpga
    );
    Ok((build(&points, cpu, fpga, start.elapsed().as_secs_f64()), text))
}

fn sweep(sp: &StudyParams, progress: Progress, cancel: &AtomicBool) -> Result<(Value, String), String> {
    let start = Instant::now();
    let (nl, nd) = (sp.layers_list.len(), sp.dt_list.len());
    let total = nl * nd * sp.seeds;
    let done = AtomicUsize::new(0);
    let work = |layers: usize, dt: f64| run::work_estimate(&Params { layers, dt, ..sp.base });
    let total_work: f64 =
        sp.layers_list.iter().flat_map(|&l| sp.dt_list.iter().map(move |&d| (l, d))).map(|(l, d)| work(l, d)).sum::<f64>()
            * sp.seeds as f64;
    let done_work = Mutex::new(0.0);
    let info = sp.base.problem.metric_info();
    let metric = (info.title, info.log, info.higher_is_better);
    let mut grid: Grid = vec![vec![None; nd]; nl];
    let (mut cpu, mut fpga) = (0.0, 0.0);
    let n = sp.base.n as f64;
    let build = |grid: &Grid, cpu: f64, fpga: f64, wall: f64| {
        let pick = |f: fn(&(f64, f64, f64)) -> f64| {
            Value::Arr(
                grid.iter()
                    .map(|row| Value::Arr(row.iter().map(|c| c.as_ref().map_or(Value::Null, |c| Value::Num(f(c)))).collect()))
                    .collect(),
            )
        };
        let mut best: Option<(usize, usize, f64)> = None;
        for (i, row) in grid.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if let Some(c) = c {
                    let better = match best {
                        None => true,
                        Some((_, _, b)) => if metric.2 { c.0 > b } else { c.0 < b },
                    };
                    if better {
                        best = Some((i, j, c.0));
                    }
                }
            }
        }
        Value::obj()
            .set("study", "sweep")
            .set("params", sp.base.to_json())
            .set("layers_list", sp.layers_list.clone())
            .set("dt_list", sp.dt_list.clone())
            .set("seeds", sp.seeds)
            .set("metric", Value::obj().set("title", metric.0).set("log", metric.1).set("higher_is_better", metric.2))
            .set("grid", pick(|c| c.0))
            .set("grid_local", pick(|c| c.1))
            .set("grid_energy", pick(|c| c.2))
            .set(
                "best",
                best.map_or(Value::Null, |(i, j, v)| {
                    Value::obj().set("layers", sp.layers_list[i]).set("dt", sp.dt_list[j]).set("value", v)
                }),
            )
            .set("cpu_seconds", cpu)
            .set("wall_seconds", wall)
            .set("fpga_seconds", fpga)
            .set("cpu", par::cpu_name())
    };
    for (i, &layers) in sp.layers_list.iter().enumerate() {
        for (j, &dt) in sp.dt_list.iter().enumerate() {
            let p = Params { layers, dt, ..sp.base };
            let w = work(layers, dt);
            let qs = runs(&p, sp.seeds, cancel, || {
                let d = done.fetch_add(1, Ordering::Relaxed) + 1;
                let mut dw = done_work.lock().unwrap();
                *dw += w;
                progress(*dw / total_work, &format!("{layers} layers, dt {dt}: {d} of {total} runs"), None);
            })?;
            cpu += qs.iter().map(|q| q.seconds).sum::<f64>();
            fpga += qs.iter().map(|q| q.fpga_seconds).sum::<f64>();
            let m = mean_std(&qs.iter().map(|q| q.metric).collect::<Vec<_>>()).0;
            let l = mean_std(&qs.iter().map(|q| q.metric_local.unwrap_or(q.metric)).collect::<Vec<_>>()).0;
            let e = mean_std(&qs.iter().map(|q| q.rounded / n).collect::<Vec<_>>()).0;
            grid[i][j] = Some((m, l, e));
            let dw = *done_work.lock().unwrap();
            progress(dw / total_work, &format!("{layers} layers, dt {dt} done"), Some(build(&grid, cpu, fpga, start.elapsed().as_secs_f64())));
        }
    }
    let out = build(&grid, cpu, fpga, start.elapsed().as_secs_f64());
    let mut text = format!(
        "{}, N = {}, {} instances per setting; {} (rows: layers, columns: dt)\n{:>8}",
        sp.base.problem.title(),
        sp.base.n,
        sp.seeds,
        metric.0,
        ""
    );
    for dt in &sp.dt_list {
        text += &format!(" {dt:>8}");
    }
    for (row, layers) in grid.iter().zip(&sp.layers_list) {
        text += &format!("\n{layers:>8}");
        for c in row {
            text += &format!(" {:>8.4}", c.map_or(f64::NAN, |c| c.0));
        }
    }
    text += &format!("\n  time: {cpu:.1} s of simulation on {}; FPGA estimate {fpga:.2} s", par::cpu_name());
    Ok((out, text))
}

// ---------------------------------------------------------------
// One-axis twisting: H = chi S_z^2 from all spins along +x.
// Time is chi t throughout.
// ---------------------------------------------------------------

/// <S_x> / (N/2) and the squeezing parameter xi^2 = N V_min / <S_x>^2, where
/// V_min is the smallest variance perpendicular to the mean spin (Wineland).
#[derive(Clone, Copy)]
pub struct Squeeze {
    pub sx: f64,
    pub xi2: f64,
}

/// Exact quantum result in the symmetric (Dicke) states, O(N) per time.
pub fn oat_exact(n: usize, times: &[f64]) -> Vec<Squeeze> {
    let s = n as f64 / 2.0;
    let mut lnf = vec![0.0; n + 1];
    for k in 1..=n {
        lnf[k] = lnf[k - 1] + (k as f64).ln();
    }
    // coherent state along +x: amplitude sqrt(C(n, k) / 2^n) on m = k - n/2
    let c: Vec<f64> = (0..=n).map(|k| (0.5 * (lnf[n] - lnf[k] - lnf[n - k]) - 0.5 * n as f64 * LN_2).exp()).collect();
    let m = |k: usize| k as f64 - s;
    let up = |k: usize| (s * (s + 1.0) - m(k) * (m(k) + 1.0)).max(0.0).sqrt(); // <m+1|S+|m>
    let sz2: f64 = (0..=n).map(|k| c[k] * c[k] * m(k) * m(k)).sum();
    times
        .iter()
        .map(|&t| {
            // <S+>, Im <S+ Sz + Sz S+>, Re <S+^2>
            let (mut sp_re, mut sp_im, mut x_im, mut sp2_re) = (0.0, 0.0, 0.0, 0.0);
            for k in 0..n {
                let w = c[k + 1] * c[k] * up(k);
                let (sn, cs) = (t * (2.0 * m(k) + 1.0)).sin_cos();
                sp_re += w * cs;
                sp_im += w * sn;
                x_im += w * (2.0 * m(k) + 1.0) * sn;
            }
            for k in 0..n.saturating_sub(1) {
                let w = c[k + 2] * c[k] * up(k) * up(k + 1);
                sp2_re += w * (t * (4.0 * m(k) + 4.0)).cos();
            }
            let syy = 0.5 * (s * (s + 1.0) - sz2) - 0.5 * sp2_re - sp_im * sp_im;
            let vmin = min_variance(syy, sz2, 0.5 * x_im);
            Squeeze { sx: sp_re / s, xi2: n as f64 * vmin / (sp_re * sp_re) }
        })
        .collect()
}

/// Smallest eigenvalue of the covariance [[a, c], [c, b]].
fn min_variance(a: f64, b: f64, c: f64) -> f64 {
    (a + b) / 2.0 - (((a - b) / 2.0).powi(2) + c * c).sqrt()
}

/// Sums over DTWA trajectories, per time.
#[derive(Clone)]
struct Moments {
    count: f64,
    sum: Vec<[f64; 3]>,
    /// yy, zz, yz
    sq: Vec<[f64; 3]>,
}

impl Moments {
    fn new(times: usize) -> Moments {
        Moments { count: 0.0, sum: vec![[0.0; 3]; times], sq: vec![[0.0; 3]; times] }
    }

    fn add(&mut self, o: &Moments) {
        self.count += o.count;
        for (a, b) in self.sum.iter_mut().zip(&o.sum) {
            for k in 0..3 {
                a[k] += b[k];
            }
        }
        for (a, b) in self.sq.iter_mut().zip(&o.sq) {
            for k in 0..3 {
                a[k] += b[k];
            }
        }
    }

    fn squeeze(&self, n: usize) -> Vec<Squeeze> {
        let c = self.count.max(1.0);
        self.sum
            .iter()
            .zip(&self.sq)
            .map(|(s, q)| {
                let (mx, my, mz) = (s[0] / c, s[1] / c, s[2] / c);
                let vyy = q[0] / c - my * my;
                let vzz = q[1] / c - mz * mz;
                let vyz = q[2] / c - my * mz;
                Squeeze { sx: mx / (n as f64 / 2.0), xi2: n as f64 * min_variance(vyy, vzz, vyz) / (mx * mx) }
            })
            .collect()
    }
}

/// One DTWA trajectory: each spin starts at (1, +-1, +-1) (the discrete
/// Wigner sampling of +x), then turns about z by chi t sum_{j != i} s_j^z,
/// spin by spin, as the FPGA pipeline would.
fn trajectory(n: usize, times: &[f64], rng: &mut Rng, acc: &mut Moments) {
    let sy: Vec<f64> = (0..n).map(|_| if rng.coin() { 1.0 } else { -1.0 }).collect();
    let sz: Vec<f64> = (0..n).map(|_| if rng.coin() { 1.0 } else { -1.0 }).collect();
    let z_total: f64 = sz.iter().sum();
    for (k, &t) in times.iter().enumerate() {
        let (mut x_sum, mut y_sum) = (0.0, 0.0);
        for (&y, &z) in sy.iter().zip(&sz) {
            let (s, c) = sincos(t * (z_total - z));
            x_sum += c - y * s;
            y_sum += s + y * c;
        }
        let (sx, sy_, sz_) = (x_sum / 2.0, y_sum / 2.0, z_total / 2.0);
        acc.sum[k][0] += sx;
        acc.sum[k][1] += sy_;
        acc.sum[k][2] += sz_;
        acc.sq[k][0] += sy_ * sy_;
        acc.sq[k][1] += sz_ * sz_;
        acc.sq[k][2] += sy_ * sz_;
    }
    acc.count += 1.0;
}

fn db(xi2: f64) -> f64 {
    10.0 * xi2.log10()
}

fn squeezing(sp: &StudyParams, progress: Progress, cancel: &AtomicBool) -> Result<(Value, String), String> {
    let start = Instant::now();
    let n = sp.base.n;
    let t_max = sp.t_max.unwrap_or(2.0 / (n as f64).sqrt());
    let times: Vec<f64> = (0..sp.times).map(|k| t_max * k as f64 / (sp.times - 1) as f64).collect();
    let exact = oat_exact(n, &times);
    let threads = sp.base.threads.max(1);
    let batches = 20.min(sp.trajectories);
    let mut acc = Moments::new(times.len());
    let mut cpu = 0.0;
    // one spin per clock at 300 MHz, like VIRAL_n_spin, per trajectory and time
    let fpga = sp.trajectories as f64 * times.len() as f64 * (n as f64 + 80.0) / 300e6;
    let curve = |pts: &[Squeeze]| {
        Value::obj()
            .set("xi2_db", pts.iter().map(|p| db(p.xi2)).collect::<Vec<_>>())
            .set("sx", pts.iter().map(|p| p.sx).collect::<Vec<_>>())
    };
    let best = |pts: &[Squeeze]| {
        let (i, p) = pts.iter().enumerate().skip(1).min_by(|a, b| a.1.xi2.total_cmp(&b.1.xi2)).unwrap();
        Value::obj().set("t", times[i]).set("xi2_db", db(p.xi2))
    };
    let build = |acc: &Moments, cpu: f64, wall: f64| {
        let dtwa = acc.squeeze(n);
        Value::obj()
            .set("study", "squeezing")
            .set("n", n)
            .set("trajectories", sp.trajectories)
            .set("done", acc.count as u64)
            .set("times", times.clone())
            .set("exact", curve(&exact))
            .set("dtwa", curve(&dtwa))
            .set("best_exact", best(&exact))
            .set("best_dtwa", best(&dtwa))
            .set("cpu_seconds", cpu)
            .set("wall_seconds", wall)
            .set("fpga_seconds", fpga)
            .set("cpu", par::cpu_name())
    };
    let mut first = 0;
    for b in 0..batches {
        let last = sp.trajectories * (b + 1) / batches;
        let count = last - first;
        let t0 = Instant::now();
        let parts = par::for_items(count, par::threads_for(count * n * times.len(), threads), || Moments::new(times.len()), |k, m| {
            if cancel.load(Ordering::Relaxed) {
                return;
            }
            let mut rng = Rng::new(sp.base.seed.wrapping_mul(1_000_003).wrapping_add((first + k) as u64));
            trajectory(n, &times, &mut rng, m);
        });
        if cancel.load(Ordering::Relaxed) {
            return Err("cancelled".into());
        }
        for m in &parts {
            acc.add(m);
        }
        cpu += t0.elapsed().as_secs_f64();
        first = last;
        progress(
            last as f64 / sp.trajectories as f64,
            &format!("{last} of {} trajectories", sp.trajectories),
            Some(build(&acc, cpu, start.elapsed().as_secs_f64())),
        );
    }
    let dtwa = acc.squeeze(n);
    let (ie, pe) = exact.iter().enumerate().skip(1).min_by(|a, b| a.1.xi2.total_cmp(&b.1.xi2)).unwrap();
    let (id, pd) = dtwa.iter().enumerate().skip(1).min_by(|a, b| a.1.xi2.total_cmp(&b.1.xi2)).unwrap();
    let text = format!(
        "One-axis twisting, N = {n}, {} trajectories:\n  best squeezing, exact: {:.2} dB at chi t = {:.4}\n  best squeezing, DTWA:  {:.2} dB at chi t = {:.4}\n  mean field (one classical run): 0 dB, no squeezing\n  time: {cpu:.1} s on {}; FPGA estimate {fpga:.2} s",
        sp.trajectories,
        db(pe.xi2),
        times[ie],
        db(pd.xi2),
        times[id],
        par::cpu_name()
    );
    Ok((build(&acc, cpu, start.elapsed().as_secs_f64()), text))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Kitagawa and Ueda (1993), with mu = 2 chi t.
    fn kitagawa_ueda(n: usize, t: f64) -> f64 {
        let nf = n as f64;
        let mu = 2.0 * t;
        let a = 1.0 - mu.cos().powf(nf - 2.0);
        let b = 4.0 * (mu / 2.0).sin() * (mu / 2.0).cos().powf(nf - 2.0);
        let xs = 1.0 + (nf - 1.0) / 4.0 * (a - (a * a + b * b).sqrt());
        xs / (mu / 2.0).cos().powf(2.0 * (nf - 1.0))
    }

    #[test]
    fn exact_matches_kitagawa_ueda() {
        for n in [4usize, 50, 1000] {
            let times: Vec<f64> = (1..20).map(|k| k as f64 * 0.5 / (n as f64).sqrt() / 10.0).collect();
            let ex = oat_exact(n, &times);
            for (t, e) in times.iter().zip(&ex) {
                let want = kitagawa_ueda(n, *t);
                assert!((e.xi2 - want).abs() < 1e-8 * want.max(1.0), "n {n} t {t}: {} vs {want}", e.xi2);
                let sx = (*t).cos().powf(n as f64 - 1.0);
                assert!((e.sx - sx).abs() < 1e-9, "n {n} t {t}: sx {} vs {sx}", e.sx);
            }
        }
    }

    #[test]
    fn dtwa_is_close_to_exact() {
        let n = 100;
        let times = [0.0, 0.03, 0.06, 0.09];
        let mut acc = Moments::new(times.len());
        let mut rng = Rng::new(1);
        for _ in 0..20_000 {
            trajectory(n, &times, &mut rng, &mut acc);
        }
        let d = acc.squeeze(n);
        let e = oat_exact(n, &times);
        for k in 0..times.len() {
            assert!((db(d[k].xi2) - db(e[k].xi2)).abs() < 0.5, "t {}: {} vs {} dB", times[k], db(d[k].xi2), db(e[k].xi2));
        }
    }
}
