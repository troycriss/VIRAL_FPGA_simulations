# viral-sim: CPU simulator and web viewer

A Rust program that runs the same mean-field QAOA as the FPGA design, on a
CPU:

- **A bit-exact model of `VIRAL_n_spin`:** the same integers, rounding and
  CORDIC, so it gives the same final spins, energy, LEDs and clock count as
  the Verilog. `check_hw.sh` checks this against Icarus Verilog.
- **A floating-point model** for six problems (below): the rank-1 coupling
  of the hardware, dense J, couplings by distance, and rank-P couplings.
- **Classical baselines** to judge each result against: exact optimum,
  Karmarkar-Karp and greedy for partitioning; local search from the
  mean-field result and from random starts; known limits (the SK ground
  state, QAOA's cut fractions, the chain's ground state, a stored pattern).
- **Studies:** many runs with the results as curves (below).
- **A web page** (`../web`) that shows the problem, the spins layer by layer
  on the Bloch sphere, the energy and the problem's own measure per layer,
  the result against the baselines, and the studies.
- **A benchmark** of this CPU against the FPGA estimates.

## Problems

| Problem | Coupling J | What it shows |
|---|---|---|
| `partition` | rank 1, J<sub>ij</sub> = a<sub>i</sub>a<sub>j</sub> (what `VIRAL_n_spin` computes) | Split numbers into two sets with equal sums |
| `sk`, `sk-gauss` | dense random, ±1/√N or Gaussian | The Sherrington-Kirkpatrick spin glass |
| `maxcut` | 1 on each edge of a random regular graph (sparse, stored dense) | The standard QAOA benchmark |
| `chain` | J<sub>ij</sub> = 1/\|i − j\|<sup>α</sup>, antiferromagnetic | The Ising chain of trapped-ion and dipolar quantum simulators |
| `hopfield` | J = −(1/N) Σ<sub>μ</sub> ξ<sup>μ</sup>ξ<sup>μᵀ</sup>, rank P | Associative memory: P stored patterns |

## Studies

| Study | What it computes |
|---|---|
| `scaling` | SK energy per spin for several N and instances, fitted against N<sup>−2/3</sup> and extrapolated to N → ∞ |
| `sweep` | Any problem's result over a grid of layers and dt: the schedule tuning MF-AOA needs |
| `squeezing` | One-axis twisting (H = χS<sub>z</sub>², all spins along x at the start), averaged over many classical runs from random starts (the discrete truncated Wigner method), against the exact quantum result |

It needs only Rust (no other packages), so it builds offline.

## Usage

```sh
cd simulator
cargo build --release
./target/release/viral-sim serve          # then open http://127.0.0.1:8080/
./target/release/viral-sim run --problem sk --n 500 --layers 200 --dt 0.3
./target/release/viral-sim study --study scaling --sizes 500,1000,2000,4000,8000 --seeds 20
./target/release/viral-sim bench
sh check_hw.sh                            # Rust vs Verilog, bit for bit
cargo test --release
```

`serve` listens on 127.0.0.1 only, so nothing outside this computer can
reach it. The page runs simulations through `/api/run`, and studies through
`/api/study/start`, `/status` and `/cancel`, which report progress while
they run. It can also open a file from `viral-sim run ... --json out.json`
or `viral-sim study ... --json out.json`, which works without the server.

`run` options (also the `/api/run` parameters):

| Option | Meaning |
|---|---|
| `--problem` | `partition`, `sk`, `sk-gauss`, `maxcut`, `chain`, `hopfield` |
| `--n` | number of spins: up to 2,000,000 for partition, 1,000,000 for the chain, 200,000 for Hopfield, 12,000 for dense J (SK, Max-Cut) |
| `--layers`, `--dt` | linear ramp as in MF-AOA: γ<sub>k</sub> = s<sub>k</sub>·dt, β<sub>k</sub> = (1 − s<sub>k</sub>)·dt, s<sub>k</sub> = (k + ½)/layers |
| `--engine` | `hw` (bit-exact `VIRAL_n_spin`, partition only, the default there) or `float` |
| `--schedule` | `ramp`, or `lut` (the 4 entries of `angle_lut`, the hardware's test pattern) |
| `--init` | `x` (along +x, the mixer's ground state, with z tilted randomly by up to ±`--tilt`) or `random` |
| `--degree` | Max-Cut: edges per vertex (3) |
| `--alpha` | chain: J = 1/r<sup>α</sup> (3, dipolar) |
| `--patterns`, `--noise` | Hopfield: stored patterns (10); fraction of pattern 1 flipped to make the cue (0.2), which the start is tilted toward |
| `--seed`, `--tilt`, `--gamma-lo`, `--points`, `--threads` | random seed; tilt (default 0.05); low 16 bits of γ for `lut`; spins kept per layer for the sphere; threads |

Without `--n`, `--layers` or `--dt`, each problem has its own defaults.

`study` options: `--study scaling|sweep|squeezing`, the `run` options for the
problem, and `--sizes 250,500,...` and `--seeds` (scaling), `--layers-list
25,50,...`, `--dt-list 0.1,0.2,...` and `--seeds` (sweep), `--n`,
`--trajectories`, `--times` and `--t-max` (squeezing).

`viral-sim hw --init FILE [--weights FILE] [--layers L] [--gamma-lo X]
[--dump FILE]` runs the bit-exact model on the testbench's hex files with
the `angle_lut` schedule.

## The model

Per layer k, with h<sub>i</sub> = Σ<sub>j≠i</sub> J<sub>ij</sub> z<sub>j</sub>:

1. rotate spin i about z by 2γ<sub>k</sub> h<sub>i</sub> / ĥ
2. rotate every spin about x by −2β<sub>k</sub>

This is mean-field (product-state) evolution under e<sup>−iγH</sup> followed by
e<sup>+iβΣX</sup>, starting from the mixer's ground state along +x. ĥ =
√(Σ<sub>i≠j</sub> J<sub>ij</sub>² / N) is the size of a typical field, so the same dt works
for any problem and size. The energy is Σ<sub>i&lt;j</sub> J<sub>ij</sub> z<sub>i</sub> z<sub>j</sub>, and at the end
each spin is rounded to the sign of z.

All spins along +x is a fixed point (every field is zero), so the start is
tilted slightly to break the symmetry.

With `--engine hw` the angles are converted to the hardware's γ (Q1.31) and
β (32-bit angle). The floating-point model then runs from the same rounded
start with exactly those angles, for comparison.

## What it shows

- **SK spin glass:** mean field works well. At N = 500, 200 layers and
  dt = 0.3, the rounded result is −0.73 per spin, against −0.67 for a random
  start plus local search. The ground state for N → ∞ is −0.763.
- **Number partitioning** is hard for mean field. The rank-1 coupling has one
  very stiff direction (eigenvalue Σa<sub>i</sub>²) and a flat landscape
  everywhere else. At 24 spins it gets residuals around 1, against about 3 for
  a random split, but Karmarkar-Karp and the exact optimum do far better. At
  1,000 spins the large steps along the stiff direction heat the spins and
  the energy rises. Local search from the mean-field result helps a lot.
- **Max-Cut (random 3-regular graphs):** mean field cuts 90% of the edges at
  100 vertices and 89.8% at 2,000, against 84.6% for local search from random
  starts. QAOA with one layer cuts 69.24% on these graphs (Farhi et al.), and with two layers
  75.59% (Wurtz and Love).
- **Long-range chain (α = 3):** at 100 sites and 1,000 layers the anneal
  reaches the exact ground state, perfect alternation, in most instances. Faster anneals leave
  defects, neighbours that point the same way: on average 11 at 50 layers
  and 2 at 1,000 (5 instances each). Single spin flips can't remove a
  defect, so local search doesn't help.
- **Hopfield memory (400 spins):** annealing finds a stored pattern only
  when there are few: always at 10 patterns (load 0.025), sometimes at 20,
  never at 30 or more. That matches theory: below a load of about 0.05 the
  stored patterns are the lowest-energy states, and above it they are only
  local minima (Amit, Gutfreund and Sompolinsky, 1985). Annealing looks for
  low energy anywhere, so it also ignores the cue and may find another
  stored pattern, or its inverse. Classic recall, downhill from the cue,
  works up to a load of about 0.138.
- **SK scaling:** with 200 layers and dt = 0.3 (the quick preset), mean field
  extrapolates to −0.732 ± 0.001 per spin as N → ∞, 4% above the ground
  state. With 500 layers and dt = 0.5, over N = 500 to 12,000 with 20
  instances each (the big preset), it extrapolates to −0.7543 ± 0.0002, only
  1.2% above the ground state of −0.7632. That run takes 31 s on the M3 (20 s
  of it simulation), against about 1.3 s estimated for the Agilex 7 design.
- **Squeezing:** a single mean-field run gives no squeezing at all. Averaged
  over 2,000 runs from random starts it gives −19.24 dB at 1,000 spins,
  against −19.22 dB exactly (0.1 s on the M3). At 10,000 spins and 10,000
  runs it gives −26.33 dB against −26.13 dB (2.7 s).
- **Fixed point against float:** over many layers the bit-exact hardware model
  drifts from the floating-point model, because the dynamics amplify the
  rounding. Over 100 layers at 24 spins the largest difference is about 0.26.
  At the 2 to 4 layers of `angle_lut` it is below 0.01 (see `../sim`).

## Speed

`viral-sim bench` on an Apple M3 (4 performance and 4 efficiency cores):

| Case | N | M3, 1 thread | M3, 8 threads | FPGA estimate |
|---|---|---|---|---|
| Rank 1, bit-exact, 2 `angle_lut` layers and 3 energies | 49,152 | 0.80 ms | 0.32 ms | 0.49 ms (Cyclone 10 GX, 300 MHz) |
| Rank 1, floating point, the same | 49,152 | 0.81 ms | 0.54 ms | 0.49 ms |
| Dense 8-bit J, 1 layer | 2,000 | 117 µs | 111 µs | 2.3 µs (Agilex 7, J on chip) |
| Dense 8-bit J, 1 layer | 8,700 | 1.80 ms | 0.59 ms | 21.5 µs |

Timings vary by 10-30% from run to run while macOS background services are
busy; these are the fastest of several runs. At 2,000 spins a layer is too
short for 8 threads to pay off.

- **Rank 1:** the CPU matches the Cyclone 10 GX design, because each layer
  is O(N) work either way. The FPGA updates one spin per clock, 300 per
  microsecond; the M3 updates about 300 per microsecond on 8 threads too.
  (The floating-point path is slower than the bit-exact one because its
  sin/cos doesn't vectorize as well.)
- **Dense J:** this is where the FPGA design would pull ahead. Each layer
  reads all of J (half of it, since J is symmetric), so the CPU is limited by
  memory bandwidth: 64 GB/s here, of about 100 GB/s for this M3. An M6
  (170 GB/s) would take about 0.35 ms per layer at 8,700 spins, still about
  15 times the Agilex 7 estimate, which reads J from block RAM and is
  limited by its 12,000 multipliers instead.

**Is Rust as fast as C?** Yes. It compiles through LLVM like clang, and the
inner loops here are the same instructions C would produce. Two loops don't
vectorize automatically in either language, so on ARM they use NEON
intrinsics from Rust's standard library:

- the 16-lane CORDIC
- the 8-bit × 16-bit tiled matrix-vector product

Other CPUs use portable versions. The tests check that both give the same
results bit for bit.

How the speed was reached:

- **J layout:** J is stored as 256 × 256 tiles of the upper triangle, each
  contiguous, so it streams from memory.
- **Integer arithmetic:** spins are Q1.15 integers in the J product, as an
  FPGA would store them, which makes the product exact.
- **Work distribution:** work is handed out as threads finish, since
  efficiency cores run slower than performance cores.

## Files

| File | What it is |
|---|---|
| `src/cordic.rs` | Bit-exact `cordic_iterative`, scalar and 16 lanes |
| `src/hw.rs` | Bit-exact `VIRAL_n_spin` |
| `src/model.rs` | Floating-point model, rank-1 and dense J, fast sincos |
| `src/fft.rs` | Convolution by FFT, for couplings that depend on distance (the chain) |
| `src/problem.rs` | The six problems: instances, rounded energies, per-layer measures, random regular graphs |
| `src/baseline.rs` | Exact optimum (subset-sum DP), Karmarkar-Karp, greedy, local search |
| `src/run.rs` | One run, and the JSON the web page draws |
| `src/study.rs` | The studies, and the exact one-axis-twisting result (Dicke states) |
| `src/server.rs` | Local web server for `../web`, `/api/run` and the studies |
| `src/bench.rs` | CPU timings against the FPGA estimates |
| `check_hw.sh` | Rust against Icarus Verilog, bit for bit |
