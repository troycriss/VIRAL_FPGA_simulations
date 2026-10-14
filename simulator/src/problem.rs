//! The problems:
//! - number partitioning: rank-1 J, what VIRAL_n_spin runs
//! - Sherrington-Kirkpatrick spin glass: dense random J
//! - Max-Cut on random regular graphs: sparse J, stored dense
//! - long-range antiferromagnetic chain, J_ij = 1/|i-j|^alpha (trapped ions)
//! - Hopfield associative memory: rank-P J from stored patterns

use crate::model::{Coupling, DenseJ};
use crate::rng::Rng;
use std::collections::HashMap;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Partition,
    Sk,
    SkGauss,
    MaxCut,
    Chain,
    Hopfield,
}

impl Kind {
    pub fn parse(s: &str) -> Result<Kind, String> {
        match s {
            "partition" => Ok(Kind::Partition),
            "sk" => Ok(Kind::Sk),
            "sk-gauss" => Ok(Kind::SkGauss),
            "maxcut" => Ok(Kind::MaxCut),
            "chain" => Ok(Kind::Chain),
            "hopfield" => Ok(Kind::Hopfield),
            _ => Err(format!("unknown problem {s:?} (partition, sk, sk-gauss, maxcut, chain, hopfield)")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Kind::Partition => "partition",
            Kind::Sk => "sk",
            Kind::SkGauss => "sk-gauss",
            Kind::MaxCut => "maxcut",
            Kind::Chain => "chain",
            Kind::Hopfield => "hopfield",
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Kind::Partition => "Number partitioning",
            Kind::Sk => "Sherrington-Kirkpatrick spin glass, couplings ±1/√N",
            Kind::SkGauss => "Sherrington-Kirkpatrick spin glass, Gaussian couplings",
            Kind::MaxCut => "Max-Cut on a random regular graph",
            Kind::Chain => "Long-range antiferromagnetic chain, J = 1/r^α",
            Kind::Hopfield => "Hopfield associative memory",
        }
    }

    pub fn rank1(self) -> bool {
        self == Kind::Partition
    }

    /// The problem's own measure of the answer, per layer: for the spin
    /// vectors (mean field) and for the spins rounded to +-1.
    pub fn metric_info(self) -> MetricInfo {
        match self {
            Kind::Partition => MetricInfo {
                title: "Imbalance |sum A − sum B|",
                log: true,
                mean_field: "mean field |Σ aᵢ zᵢ|",
                higher_is_better: false,
                rounded: "rounded |Σ aᵢ sᵢ|",
            },
            Kind::Sk | Kind::SkGauss => MetricInfo {
                title: "Energy per spin",
                log: false,
                mean_field: "mean field",
                higher_is_better: false,
                rounded: "rounded to ±1",
            },
            Kind::MaxCut => MetricInfo {
                title: "Fraction of edges cut",
                log: false,
                mean_field: "mean field (1 − zᵢzⱼ)/2",
                higher_is_better: true,
                rounded: "rounded to ±1",
            },
            Kind::Chain => MetricInfo {
                title: "Fraction of neighbours alternating",
                log: false,
                mean_field: "mean field (1 − zᵢzᵢ₊₁)/2",
                higher_is_better: true,
                rounded: "rounded to ±1",
            },
            Kind::Hopfield => MetricInfo {
                title: "Largest overlap with a stored pattern",
                log: false,
                mean_field: "mean field max |Σ ξᵢ zᵢ| / N",
                higher_is_better: true,
                rounded: "rounded max |Σ ξᵢ sᵢ| / N",
            },
        }
    }

    /// Largest N for this problem (memory and time).
    pub fn max_n(self) -> usize {
        match self {
            Kind::Partition => 2_000_000,
            Kind::Sk | Kind::SkGauss | Kind::MaxCut => 12_000,
            Kind::Chain => 1_000_000,
            Kind::Hopfield => 200_000,
        }
    }
}

/// What to generate.
#[derive(Clone, Copy)]
pub struct Spec {
    pub kind: Kind,
    pub n: usize,
    /// maxcut: degree of every vertex
    pub degree: usize,
    /// chain: J_ij = 1/|i-j|^alpha
    pub alpha: f64,
    /// hopfield: number of stored patterns
    pub patterns: usize,
    /// hopfield: fraction of pattern 1 flipped to make the cue
    pub noise: f64,
}

pub struct Problem {
    pub kind: Kind,
    pub coupling: Coupling,
    /// partition: the numbers, unsigned with 15 fractional bits
    pub weights: Option<Vec<u16>>,
    /// maxcut: the graph
    pub edges: Option<Vec<(u32, u32)>>,
    /// chain: alpha
    pub alpha: f64,
    /// hopfield: pattern 1 with a fraction of its bits flipped
    pub cue: Option<Vec<i8>>,
}

impl Problem {
    pub fn generate(spec: &Spec, rng: &mut Rng) -> Result<Problem, String> {
        let n = spec.n;
        let mut p = Problem {
            kind: spec.kind,
            coupling: Coupling::Rank1 { a: Vec::new() },
            weights: None,
            edges: None,
            alpha: spec.alpha,
            cue: None,
        };
        match spec.kind {
            // like gen_weights.py --pattern random: uniform in (0, 1]
            Kind::Partition => return Ok(Problem::partition((0..n).map(|_| rng.range(1, 0x8000) as u16).collect())),
            Kind::Sk | Kind::SkGauss => {
                p.coupling = Coupling::Dense(DenseJ::sk(n, spec.kind == Kind::SkGauss, rng));
            }
            Kind::MaxCut => {
                let edges = random_regular(n, spec.degree, rng)?;
                // J = 1 on every edge: the energy is (edges - 2 cut)
                let mut q = vec![0i8; n * n];
                for &(a, b) in &edges {
                    q[a as usize * n + b as usize] = 1;
                    q[b as usize * n + a as usize] = 1;
                }
                p.coupling = Coupling::Dense(DenseJ::new(n, q, 1.0));
                p.edges = Some(edges);
            }
            Kind::Chain => {
                if !(spec.alpha >= 0.0 && spec.alpha <= 10.0) {
                    return Err("alpha: 0 to 10".into());
                }
                let f = (0..n).map(|d| if d == 0 { 0.0 } else { (d as f64).powf(-spec.alpha) }).collect();
                p.coupling = Coupling::chain(f);
            }
            Kind::Hopfield => {
                let pats = spec.patterns;
                if pats < 1 || pats > n.max(1) || (pats as f64) * (n as f64) > 5e7 {
                    return Err("patterns: 1 to N, and N x patterns at most 50 million".into());
                }
                if !(0.0..=0.5).contains(&spec.noise) {
                    return Err("noise: 0 to 0.5".into());
                }
                let xi: Vec<Vec<i8>> =
                    (0..pats).map(|_| (0..n).map(|_| if rng.coin() { 1 } else { -1 }).collect()).collect();
                p.cue = Some(xi[0].iter().map(|&x| if rng.uniform() < spec.noise { -x } else { x }).collect());
                // J = -(1/N) sum_mu xi xi^T: each stored pattern is a low-energy state
                p.coupling = Coupling::LowRank { xi, c: -1.0 / n as f64 };
            }
        }
        Ok(p)
    }

    pub fn partition(weights: Vec<u16>) -> Problem {
        let a = weights.iter().map(|&w| w as f64 / 32768.0).collect();
        Problem {
            kind: Kind::Partition,
            coupling: Coupling::Rank1 { a },
            weights: Some(weights),
            edges: None,
            alpha: 0.0,
            cue: None,
        }
    }

    pub fn n(&self) -> usize {
        self.coupling.n()
    }

    pub fn patterns(&self) -> Option<&Vec<Vec<i8>>> {
        match &self.coupling {
            Coupling::LowRank { xi, .. } => Some(xi),
            _ => None,
        }
    }

    /// Energy sum_{i<j} J_ij s_i s_j of spins rounded to s = +-1.
    pub fn rounded_energy(&self, s: &[i8], threads: usize) -> f64 {
        match &self.coupling {
            Coupling::Rank1 { a } => {
                let r: f64 = a.iter().zip(s).map(|(a, &s)| a * s as f64).sum();
                let a2: f64 = a.iter().map(|a| a * a).sum();
                (r * r - a2) / 2.0
            }
            Coupling::Dense(d) => {
                let h = d.fields_int(s, threads);
                let e: i64 = h.iter().zip(s).map(|(&h, &s)| h as i64 * s as i64).sum();
                e as f64 / 2.0 * d.scale
            }
            c => {
                let z: Vec<f64> = s.iter().map(|&v| v as f64).collect();
                let h = c.fields(&z, threads);
                0.5 * z.iter().zip(&h).map(|(a, b)| a * b).sum::<f64>()
            }
        }
    }

    pub fn metric_info(&self) -> MetricInfo {
        self.kind.metric_info()
    }

    /// (mean field, rounded) values of the metric, from the spins' z and
    /// the two energies.
    pub fn metric(&self, z: &[f64], e_mf: f64, e_round: f64) -> (f64, f64) {
        let n = z.len() as f64;
        let sign = |v: f64| if v < 0.0 { -1.0 } else { 1.0 };
        match &self.coupling {
            Coupling::Rank1 { a } => {
                let a2: f64 = a.iter().map(|a| a * a).sum();
                let mf: f64 = a.iter().zip(z).map(|(a, z)| a * z).sum();
                (mf.abs(), (2.0 * e_round + a2).max(0.0).sqrt())
            }
            _ => match self.kind {
                Kind::MaxCut => {
                    let m = self.edges.as_ref().map_or(1, |e| e.len()) as f64;
                    ((m - e_mf) / (2.0 * m), (m - e_round) / (2.0 * m))
                }
                Kind::Chain => {
                    let pairs = (n - 1.0).max(1.0);
                    let mf: f64 = z.windows(2).map(|w| (1.0 - w[0] * w[1]) / 2.0).sum();
                    let r = z.windows(2).filter(|w| sign(w[0]) != sign(w[1])).count() as f64;
                    (mf / pairs, r / pairs)
                }
                Kind::Hopfield => {
                    let (mut mf, mut r) = (0.0f64, 0.0f64);
                    for x in self.patterns().unwrap() {
                        let a: f64 = x.iter().zip(z).map(|(&x, &z)| x as f64 * z).sum();
                        let b: f64 = x.iter().zip(z).map(|(&x, &z)| x as f64 * sign(z)).sum();
                        mf = mf.max(a.abs() / n);
                        r = r.max(b.abs() / n);
                    }
                    (mf, r)
                }
                _ => (e_mf / n, e_round / n),
            },
        }
    }
}

pub struct MetricInfo {
    pub title: &'static str,
    pub log: bool,
    pub higher_is_better: bool,
    pub mean_field: &'static str,
    pub rounded: &'static str,
}

/// Rounds spins to +-1 by the sign of z (z = 0 counts as +1).
pub fn signs(z: &[f64]) -> Vec<i8> {
    z.iter().map(|&z| if z < 0.0 { -1 } else { 1 }).collect()
}

/// Random d-regular graph: random pairing of the n*d edge ends, then random
/// double edge swaps until there are no self-loops or repeated edges.
pub fn random_regular(n: usize, d: usize, rng: &mut Rng) -> Result<Vec<(u32, u32)>, String> {
    if d < 1 || d >= n || !(n * d).is_multiple_of(2) || d > 50 {
        return Err("maxcut: degree 1 to 50, below N, and N x degree must be even".into());
    }
    let mut ends: Vec<u32> = (0..n as u32).flat_map(|v| std::iter::repeat_n(v, d)).collect();
    for i in (1..ends.len()).rev() {
        let j = rng.range(0, i as u64) as usize;
        ends.swap(i, j);
    }
    let key = |a: u32, b: u32| if a < b { (a, b) } else { (b, a) };
    let mut edges: Vec<(u32, u32)> = ends.chunks(2).map(|c| key(c[0], c[1])).collect();
    let mut count: HashMap<(u32, u32), u32> = HashMap::new();
    for &e in &edges {
        *count.entry(e).or_insert(0) += 1;
    }
    let bad = |e: (u32, u32), count: &HashMap<(u32, u32), u32>| e.0 == e.1 || count[&e] > 1;
    let m = edges.len();
    for _ in 0..200 * m {
        let Some(i) = (0..m).find(|&i| bad(edges[i], &count)) else {
            return Ok(edges);
        };
        // swap ends with a random other edge: (a,b),(c,d) -> (a,c),(b,d) or (a,d),(b,c)
        let j = rng.range(0, m as u64 - 1) as usize;
        let ((a, b), (c, dd)) = (edges[i], edges[j]);
        let (e1, e2) = if rng.coin() { (key(a, c), key(b, dd)) } else { (key(a, dd), key(b, c)) };
        if i == j || e1.0 == e1.1 || e2.0 == e2.1 || count.contains_key(&e1) || count.contains_key(&e2) || e1 == e2 {
            continue;
        }
        for e in [edges[i], edges[j]] {
            let c = count.get_mut(&e).unwrap();
            *c -= 1;
            if *c == 0 {
                count.remove(&e);
            }
        }
        edges[i] = e1;
        edges[j] = e2;
        count.insert(e1, 1);
        count.insert(e2, 1);
    }
    Err("could not make a simple random regular graph; try another seed".into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_graph_is_simple_and_regular() {
        let mut rng = Rng::new(3);
        for (n, d) in [(10, 3), (100, 3), (1000, 3), (200, 7)] {
            let e = random_regular(n, d, &mut rng).unwrap();
            assert_eq!(e.len(), n * d / 2);
            let mut deg = vec![0; n];
            let mut seen = std::collections::HashSet::new();
            for &(a, b) in &e {
                assert!(a < b && seen.insert((a, b)));
                deg[a as usize] += 1;
                deg[b as usize] += 1;
            }
            assert!(deg.iter().all(|&k| k == d));
        }
    }
}
