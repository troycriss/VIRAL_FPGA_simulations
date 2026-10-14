# VIRAL FPGA simulations

Verilog for simulating interacting spins on an FPGA. Each spin is a classical
3D vector (x, y, z). In every layer, each spin is rotated about z by an angle
set by a parameter γ and the z components of the other spins, then rotated
about x by an angle β. After the last layer the design computes the energy
Σ<sub>i&lt;j</sub> z<sub>i</sub> z<sub>j</sub> and shows it on 8 LEDs. This is the structure of a
mean-field version of QAOA: the z rotation plays the role of the cost step
and the x rotation the mixer. `VIRAL_n_spin` can also give each spin a
weight a<sub>i</sub>, which turns the problem into number partitioning.

The repository has the original two- and three-spin designs, and
`VIRAL_n_spin`, a version that scales to tens of thousands of spins.

## What it computes

It computes a mean-field version of QAOA, the mean-field approximate
optimization algorithm (MF-AOA) of [Misra-Spieldenner et al., PRX Quantum
2023](https://arxiv.org/abs/2303.00329). Each qubit is replaced by a
classical spin vector, which captures no entanglement. Given γ and β for
each layer and a starting state, it evolves the spins and reports the energy.

With weights a<sub>i</sub> (`WEIGHT_FILE`), the coupling is J<sub>ij</sub> = a<sub>i</sub> a<sub>j</sub> and the
energy is Σ<sub>i&lt;j</sub> a<sub>i</sub> a<sub>j</sub> z<sub>i</sub> z<sub>j</sub>. Up to a constant this is the
number-partitioning cost (Σ<sub>i</sub> a<sub>i</sub> z<sub>i</sub>)² / 2: split the numbers a<sub>i</sub> into
two sets (z > 0 and z < 0) with sums as equal as possible. With random
weights this problem is NP-hard, and it is one of the two benchmarks in
the MF-AOA paper. Without weights every a<sub>i</sub> is 1, the uniform model of the
original designs.

The field on spin i is h<sub>i</sub> = Σ<sub>j≠i</sub> J<sub>ij</sub> z<sub>j</sub> = a<sub>i</sub>(S − a<sub>i</sub> z<sub>i</sub>), where
S = Σ<sub>j</sub> a<sub>j</sub> z<sub>j</sub>. Because J<sub>ij</sub> = a<sub>i</sub> a<sub>j</sub> has rank 1, each layer costs O(N)
instead of O(N²).

### Layers and spins

- **The number of spins, N,** is the problem size: how many spins
  (variables) there are, each a 3D vector. Memory grows with N, and so does
  the time per layer in `VIRAL_n_spin`, since one spin goes through the
  pipeline per clock.
- **Layers** are the steps in time, like the depth *p* in QAOA. Each layer k
  applies the cost step (rotate about z by γ<sub>k</sub>·h<sub>i</sub>) and then the mixer
  (rotate about x by β<sub>k</sub>), with its own γ<sub>k</sub> and β<sub>k</sub>. Layers run one after
  another and add time but no memory. There are at most 4, because
  `angle_lut` has 4 entries.

### Limitations

- **Without weights the problem is easy.** Every pair is coupled equally,
  so the energy is just (S² − Σ z<sub>i</sub>²)/2, where S = Σ z<sub>i</sub>. It depends only on
  the total z, and the minimum is any state with S = 0.
- **Rank-1 coupling is cheap on a CPU too.** Spins interact only through S,
  with or without weights. The bit-exact CPU model in `simulator/` runs
  49,152 spins over 2 layers in 0.32 ms on an Apple M3, against about
  0.49 ms estimated for the FPGA (see [CPU simulator](#cpu-simulator-and-web-page)).
- **Little comes out.** The only output is 8 LEDs of energy. Useful runs
  would need the final spins read out (the signs of z are a candidate
  solution), and γ and β loaded from a host computer so that an optimizer
  can search over them.

Number partitioning gives results that can be compared directly with the
MF-AOA paper and with [QAOA.jl](https://joss.theoj.org/papers/10.21105/joss.05364),
which implements MF-AOA in Julia. The FPGA itself becomes worthwhile with a
general coupling J<sub>ij</sub>, such as MaxCut on a sparse graph or the
Sherrington–Kirkpatrick model. Then
h<sub>i</sub> = Σ<sub>j</sub> J<sub>ij</sub> z<sub>j</sub> is the expensive part, which streaming hardware is good
at. That is the setting of published mean-field QAOA heuristics and of FPGA
Ising-solver hardware. On the Cyclone 10 GX below, J would have to be sparse
to fit in block RAM alongside the spins.

## Files

| File | What it is |
|---|---|
| `cordic_iterative.v` | Rotates a 2D vector by an angle (CORDIC), one iteration per clock |
| `VIRAL_rotation_2.v` | Rotates one spin about z, then about x, using one `cordic_iterative` twice |
| `angle_lut.v` | γ and β for each layer (4 entries). γ comes in 16 bits for the original designs and 32 bits for `VIRAL_n_spin` |
| `VIRAL_two_spin_2.v` | Original top level: 2 spins, 1 layer |
| `VIRAL_three_spin.v` | Original top level: 3 spins, 2 layers |
| `VIRAL_n_spin.v` | **Scalable top level**: N spins, same ports as `VIRAL_three_spin` |
| `VIRAL_rotation_pipe.v` | Streaming version of `VIRAL_rotation_2`: one spin per clock |
| `cordic_pipelined.v` | Pipelined version of `cordic_iterative` with bit-identical results |
| `spin_memory.v` | Initial and current spin states and the spin weights (block RAM), with registered inputs and outputs |
| `serial_mult.v` | Shift-and-add multiplier, one add per clock, for the products `VIRAL_n_spin` needs once per layer |
| `spin_init_3.hex` | Initial state of `VIRAL_three_spin`, for `VIRAL_n_spin` |
| `gen_spin_init.py` | Writes initial-state files for any number of spins |
| `gen_weights.py` | Writes spin weight files (`WEIGHT_FILE`) for `VIRAL_n_spin` |
| `sim/` | Icarus Verilog testbenches: `sh sim/run_sims.sh` |
| `synth/` | Yosys resource comparison: `sh synth/run_synth.sh` |
| `simulator/` | Rust CPU simulator: bit-exact model of `VIRAL_n_spin`, dense-J model, benchmark, web server ([README](simulator/README.md)) |
| `web/` | Web page that shows a run: the problem, the spins on the Bloch sphere, the energy per layer, the result |

## Number formats

- **Spin components** are 16-bit signed Q1.15, so `0x7FFF` ≈ +1.
  `0x5A82` is 1/√2 and `0xA57E` is −1/√2.
- **Angles** are 32 bits, where 2³² is a full turn: `0x20000000` = 45°,
  `0x40000000` = 90°. Overflow wraps around, which is still the correct angle.
- **γ** is Q1.31 in `VIRAL_n_spin` (`gamma_out` from `angle_lut`) and Q1.15
  in the original designs (`alpha_out`, the upper 16 bits of the same
  values). Spin i is rotated about z by γ·h<sub>i</sub>·180°, where h<sub>i</sub> is the field
  from the other spins (their z, times the weights if there are any).
- **Weights** are 16-bit unsigned with 15 fractional bits, so `0x8000` =
  1.0 and `0xFFFF` is just under 2. w<sub>i</sub> = a<sub>i</sub> z<sub>i</sub> is rounded to Q1.15; it
  can reach ±2, so it has 17 bits.

## How the original design works

`VIRAL_three_spin` runs each layer like this:

1. Read γ and β from `angle_lut`.
2. For each spin, add up the z of the other spins (h<sub>i</sub>) and multiply by γ
   to get its z-rotation angle.
3. Start all three `VIRAL_rotation_2` modules at once. Each uses its CORDIC
   twice: about z, then about x (about 40 clocks).
4. When all three finish, use the results as the next layer's input.

After the last layer it multiplies the z of every pair of spins, adds the
products, and puts bits [15:8] of the energy on the LEDs.

## Why the original design doesn't scale

Everything in it is copied once per spin or once per pair of spins:

- **Per spin:** a `VIRAL_rotation_2` with its own iterative CORDIC (about
  1,300 LUTs and 300 registers on Cyclone IV), about 150 registers of spin
  state and angles, and a multiplier for γ·h<sub>i</sub>.
- **Field:** each spin adds up the z of every other spin, so N(N−1) additions.
- **Energy:** one multiplier and one 32-bit register per pair, so N(N−1)/2.
  That's 4,950 multipliers at 100 spins.
- **One clock for every sum:** the field terms and the pair terms are each
  added in a single clock, so the logic between registers gets deeper as N
  grows.

The time per layer stays constant, because all spins rotate in parallel, but
the hardware grows as N². On Cyclone IV GX it runs out of multipliers at
around 25 spins and out of logic at around 50 (see the comparison below).

## The scalable design: `VIRAL_n_spin`

`VIRAL_n_spin` uses the same model, the same ports (`clk`, `reset`, `led`)
and the same `angle_lut`, and adds optional spin weights. It makes five
changes. Below, w<sub>i</sub> = a<sub>i</sub> z<sub>i</sub> (just z<sub>i</sub> without weights) and S = Σ w<sub>i</sub>.

1. **Energy in O(N).** Σ<sub>i&lt;j</sub> w<sub>i</sub> w<sub>j</sub> = (S² − Σ w<sub>i</sub>²) / 2. S and Σ w<sub>i</sub>² are
   added up as each spin is written back, so the energy costs the same for
   any N. It is also exact: the original rounds each pair term down, and
   those errors add up at large N.
2. **Field in O(N).** h<sub>i</sub> = a<sub>i</sub>(S − w<sub>i</sub>), using S from the start of the
   layer, instead of N − 1 additions per spin. The rotation angle is
   computed as γ·h<sub>i</sub> = a<sub>i</sub>(γ·S − γ·w<sub>i</sub>), with γ·S once per layer, so the
   per-spin multiplies don't grow with N.
3. **One shared, pipelined rotation unit.**
   - `cordic_pipelined` does the same arithmetic as `cordic_iterative`, but
     each iteration has its own pipeline stage, so it accepts a new vector
     every clock.
   - `VIRAL_rotation_pipe` chains two of them: rotate about z, then about x.
   - All spins stream through it, one per clock. The spin states live in
     block RAM (`spin_memory`) instead of registers.
4. **Generalized.**
   - The number of spins is a parameter (`N_SPINS`).
   - Field and energy widths grow with log₂N, so nothing overflows.
   - The initial state comes from a file and is reloaded on every reset.
   - Each spin can have a weight a<sub>i</sub>, from `WEIGHT_FILE`. With no file all
     weights are 1.0, results are bit for bit the same as without weights,
     and no memory is used for them.
   - The LED value saturates instead of wrapping, and `LED_LSB` sets its scale.
   - γ is 32 bits (Q1.31) instead of 16. With many spins, γ has to be about
     1/N to keep the rotation angles reasonable. In Q1.15 the smallest γ,
     2⁻¹⁵, would already turn fully aligned spins by about 330° at 60,000
     spins. The angle is γ·h<sub>i</sub> rounded to the nearest angle LSB.
5. **Pipelined for clock speed.** Between any two registers there is at most
   one adder or one multiply that fits a DSP multiplier, for any N:
   - Every per-spin product fits one DSP multiplier, with registers before
     and after: a<sub>i</sub>·z<sub>i</sub> (17 × 16 bits), γ·w<sub>i</sub> as two products (16 × 17
     and 17 × 17), and a<sub>i</sub> times the 62-bit γ·S − γ·w<sub>i</sub> as four (16 × 18
     or smaller). Nine multipliers in all, with the two for S and Σ w<sub>i</sub>².
   - S² for the energy and γ·S for the field, needed once per layer, come
     from shift-and-add multipliers (`serial_mult`): one add per clock,
     17 + log₂N clocks.
   - The CORDIC gain correction is split over two clocks.
   - The spin memories register their inputs and outputs. At large N each
     memory is spread over hundreds of blocks. The address and data
     fan-out to them, and the mux that picks a block's output, then get a
     clock each.

The datapath handles one spin per clock:

```
 spin_memory --> w = a*z --> angle = a*(gamma*S - gamma*w) --> rotate about z --> rotate about x --+
  (block RAM)     (DSP)      (DSP; gamma*S once per layer)    (pipelined CORDIC) (pipelined CORDIC)|
      ^                                                                                           |
      +------------------------------------- write back <-----------------------------------------+
                                                 |
                                                 +--> w = a*z, S += w, Q += w^2
                                                      (S for the next layer, energy = (S^2 - Q) / 2)
```

As in the original, all spins in a layer are updated from the z values at
the start of that layer.

### Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `N_SPINS` | 3 | Number of spins |
| `N_LAYERS` | 2 | Number of layers. At most 4 unless `angle_lut` is extended |
| `INIT_FILE` | `"spin_init_3.hex"` | Initial spin states |
| `WEIGHT_FILE` | `""` | Spin weights. Empty: all weights 1.0 |
| `LED_LSB` | 8 | `led` = energy >>> `LED_LSB`, saturated to 8 bits. 8 matches the original for energies between −1 and 1 |

### Initial state file

The file has one line per spin: 12 hex digits for x, y and z, each 16-bit
Q1.15. For example, `5A8200005A82` is (0.707, 0, 0.707).

```
python3 gen_spin_init.py 100 -o spin_init_100.hex                          # alternating ±45° tilts
python3 gen_spin_init.py 100 --pattern random --seed 1 -o spin_init_100.hex  # random directions
```

In Quartus, the path is relative to the project directory.

### Weight file

The file has one line per spin: 4 hex digits, unsigned with 15 fractional
bits. For example, `8000` is 1.0 and `4000` is 0.5.

```
python3 gen_weights.py 100 --pattern random --seed 1 -o weights_100.hex  # uniform in (0, 1]
```

### How its results compare with the originals

- `N_SPINS = 2`, `N_LAYERS = 1` reproduces `VIRAL_two_spin_2` bit for bit.
- `N_SPINS = 3`, `N_LAYERS = 2` with `spin_init_3.hex` runs the same case as
  `VIRAL_three_spin`, but one output value differs by 2 LSB. `VIRAL_three_spin`
  halves z<sub>j</sub> + z<sub>k</sub> before multiplying by γ, which drops the last bit when
  the sum is odd; `VIRAL_n_spin` keeps it. With the same rounding added, the
  results match exactly.

## Comparison: resources

These numbers are for Cyclone IV GX, from Yosys (`sh synth/run_synth.sh`).
"Original architecture" is `synth/VIRAL_n_spin_direct.v`: `VIRAL_three_spin`
generalized to N spins, with one `VIRAL_rotation_2` per spin, per-spin field
sums, and one multiplier per pair. Simulations confirm it gives the same
final spins as `VIRAL_n_spin`.

| Design | Spins | LUTs | Registers | Multipliers (18×18 blocks) | M9K | LUT levels |
|---|---|---|---|---|---|---|
| `VIRAL_two_spin_2` | 2 | 2,178 | 563 | 3 (3) | 0 | 17 |
| `VIRAL_three_spin` | 3 | 4,114 | 1,285 | 6 (6) | 0 | 17 |
| original architecture | 4 | 5,509 | 1,696 | 10 (10) | 0 | 17 |
| original architecture | 10 | 15,851 | 4,680 | 55 (65) | 0 | 21 |
| original architecture | 20 | 37,108 | 10,940 | 210 (230) | 0 | 32 |
| original architecture | 50 | 147,503 | 39,340 | 1,275 (1,325) | 0 | 63 |
| `VIRAL_n_spin` | 3 | 7,965 | 5,139 | 9 (9) | 4 | 17 |
| `VIRAL_n_spin` | 10 | 8,056 | 5,271 | 9 (9) | 4 | 17 |
| `VIRAL_n_spin` | 100 | 8,195 | 5,469 | 9 (9) | 4 | 17 |
| `VIRAL_n_spin` | 1,000 | 8,361 | 5,667 | 9 (9) | 12 | 16 |

How to read the table:

- **LUTs and registers:** a Cyclone IV logic element (LE) holds one 4-input
  LUT and one register, so the LE count lies between the larger of the two
  and their sum.
- **Multipliers and memory:** these go into the embedded 18×18 multipliers
  and M9K memory blocks. The script counts them from the bit widths, because
  that mapping is done by Quartus, not Yosys. Multipliers wider than 18 bits
  need more than one block.
- **LUT levels:** the longest chain of LUTs between two registers, not
  counting the multipliers and memories. It is a rough guide to how fast the
  design can be clocked (see below).
- **These are estimates.** Yosys's Cyclone IV flow is experimental and does
  not use the carry chains, so adders cost more LUTs and more levels than
  they will in Quartus. Use the numbers to compare designs, and Quartus for
  exact counts.

What the table shows:

- **Original architecture:** grows roughly as N².
  - LUTs: 37k at 20 spins, 148k at 50, and by extrapolation about 480k at
    100.
  - Multipliers: about N²/2, so 1,275 at 50 spins.
  - For scale, the largest Cyclone IV GX (EP4CGX150) has about 150k LEs and
    360 18×18 multipliers.
  - It runs out of multipliers at around 25 spins and out of logic at
    around 50.
- **`VIRAL_n_spin`:** stays at about 8–8.4k LUTs, 5–5.7k registers and 9
  multipliers (one 18×18 block each) from 3 to 1,000 spins. Only the
  memory grows: 2 × 48 bits per spin, 4 M9K blocks up to 256 spins and 12
  at 1,000. A weight file adds 16 bits per spin (not in the table).
- **Crossover:** below about 5 or 6 spins the original architecture is
  smaller, because the two pipelined CORDICs (30 stages) cost more than a
  few iterative ones.

## Comparison: timing

### Clock cycles

These are clock counts from simulation, for 2 layers including loading the
initial state and computing the energy:

| Spins | Original architecture | `VIRAL_n_spin` |
|---|---|---|
| 3 | 99 | 209 |
| 4 | 99 | 212 |
| 10 | 99 | 236 |
| 100 | 99 | 515 |
| 1,000 | 99 (not simulated) | 3,224 |

- **Original:** about 50 clocks per layer for any N, because every spin has
  its own rotation module and they all run at once.
- **`VIRAL_n_spin`:** N + W + 17 clocks to load the initial state and
  compute its energy, then N + W + 63 per layer, where W = 17 + log₂N.
  Spins go through the shared pipeline one per clock; 63 is the pipeline
  latency plus control, and W is the shift-and-add multiplies for S² and
  γ·S. The pipelining for clock speed adds about 50 clocks per layer,
  which matters only for small N. At 50 MHz, 100 spins and 2 layers take
  10 µs; 1,000 spins take 64 µs.

### Clock speed

Getting the maximum clock speed needs Quartus timing analysis (TimeQuest),
which this repository doesn't include. What sets it is the deepest logic
that has to finish in one clock:

| | Original architecture | `VIRAL_n_spin` |
|---|---|---|
| Rotation | Iterative CORDIC step: a variable shift by 0–14 bits, a 17-bit add or subtract, and a 32-bit angle add or subtract with a 15-entry table lookup | Pipelined CORDIC stage: fixed shift, a 17-bit add or subtract, a 32-bit add or subtract |
| Field | Sum of N − 1 values (adder tree log₂N deep) | a(γ·S − γ·w): one add or subtract of up to 62 bits per clock |
| Energy | Sum of N(N−1)/2 pair terms (adder tree 2·log₂N deep) | S² by shift-and-add, then S² − Q: one (34 + 2·log₂N)-bit add or subtract per clock |
| Multipliers | 16 × (16 + log₂N) for the field, 16 × 16 per pair, in one clock with the surrounding logic | Nine, 18 × 18 or smaller, each with its own input and output registers |
| Memory | None | Registered inputs and outputs |

The original also does the CORDIC gain correction (a sum of four shifted
copies of the input) in one clock; `VIRAL_n_spin` splits it over two.

So the original's clock gets slower as N grows: at 100 spins the energy is
a 13-level adder tree with 4,950 inputs, finished in one clock. Its
iterative CORDIC step is also deeper than a pipelined stage at any N. In
`VIRAL_n_spin` nothing scales with N except the widths, which grow with log₂N.

The "LUT levels" column in the resource table is Yosys's longest LUT chain
between registers. It fits this picture only roughly:

- **`VIRAL_n_spin`:** 16–17 levels from 3 to 1,000 spins, on the γ·S
  shift-and-add, a 62-bit add or subtract that Quartus builds as one carry
  chain.
- **Original architecture:** 17 levels at 3–4 spins, 21 at 10, 32 at 20
  and 63 at 50. The longest path is on the CORDIC angle adder, a field sum,
  or the energy sum, depending on N. At 50 spins it is a field sum: 49
  values added in one clock.

Because Yosys builds adders from LUTs rather than carry chains, these
numbers mostly measure adder length, and they vary with how its optimizer
happens to map the adders.

### Overall

- **Small N:** the original takes fewer clocks (99 vs 209 at 3 spins), but
  `VIRAL_n_spin`'s shallower logic per clock should let it run at a faster
  clock.
- **Large N:** `VIRAL_n_spin` takes about N + 100 clocks per layer at a
  clock speed that doesn't depend on N. The original keeps about 50 clocks per layer,
  but its clock slows down as N grows and it stops fitting on the chip
  beyond a few tens of spins.

## Estimate: Cyclone 10 GX development kit

The Intel Cyclone 10 GX FPGA Development Kit has a 10CX220YF780E5G, the
largest Cyclone 10 GX: 80,330 ALMs, 587 M20K memory blocks (11.7 Mb), 192
DSP blocks, speed grade −E5. These estimates are for `VIRAL_n_spin` as it
is, with the initial-state copy, and have not been checked in Quartus.

**Logic** is about the same for any N. Yosys's ALM flow targets Cyclone V,
whose logic cell is similar. At 65,536 spins it gives about 4,900 adder
bits, 2,100 other LUTs, 6,300 registers and 9 multipliers that each fit an
18 × 19 DSP multiplier. That is roughly 5,500 ALMs (7%) and 5 DSP blocks (3%).

**Memory** limits N. Each spin's state is 48 bits, stored twice (`init_rom`
and `spin_mem`). In the 2K × 10 configuration, 2,048 spins take 5 M20K
blocks per copy. A weight file adds 2 blocks per 2,048 spins (2K × 8, twice).
So 2,048 spins take 10 blocks without weights and 12 with:

| Spins | Weights | M20K blocks | Share of M20K |
|---|---|---|---|
| 57,344 | no | 280 | 48% |
| 65,536 | no | 320 | 55% |
| about 118,000 | no | 580 | 99% |
| 49,152 | yes | 288 | 49% |
| about 98,000 | yes | 576 | 98% |

Without the initial-state copy, each count would hold about twice as many
spins.

**Clock speed.** The datasheet limits for −E5 are 456 MHz for an 18 × 19
multiplier with its registers, 625 MHz for M20K blocks and 644 MHz for the
global clock network. The design's longest paths are single carry chains
(up to 34 + 2·log₂N = 66 bits, in the S² shift-and-add and the Σ w²
accumulator), and the routing from the memory blocks' output registers
through the mux that picks a block. A reasonable estimate is about 300 MHz.

**Final estimate: about 60,000 spins without weights, or about 49,000 with
weights (number partitioning), at about 300 MHz.** A run of L layers takes
(L + 1)(N + W) + 63L + 17 clocks, with W = 17 + log₂N:

| Spins | Weights | Layers | Clocks | Time at 300 MHz |
|---|---|---|---|---|
| 57,344 | no | 2 | 172,274 | 0.57 ms (191 µs per layer) |
| 65,536 | no | 2 | 196,850 | 0.66 ms (219 µs per layer) |
| 49,152 | yes | 2 | 147,698 | 0.49 ms (164 µs per layer) |

That is one spin update per clock, about 3 × 10⁸ per second. For
meaningful dynamics at this size, set γ to about 1/N (see the γ note under
the scalable design).

## Comparison with published hardware

How the Cyclone 10 GX estimate above compares with published FPGA and chip
Ising solvers. Figures for the other machines are from their papers and
press material; figures for this design are estimates.

| | Hardware | Coupling | Spins | Work per clock |
|---|---|---|---|---|
| **This design** (estimate) | One mid-range FPGA (Cyclone 10 GX) | Rank 1 only (J<sub>ij</sub> = a<sub>i</sub> a<sub>j</sub>) | About 49,000 (weights) or 60,000 | 1 spin update, about 3 × 10⁸ per second |
| **Toshiba simulated bifurcation machine** ([Kashimata et al.](https://arxiv.org/abs/2311.17370)) | One or more FPGAs | Dense, any J | Thousands per FPGA; 32,768 on 8 FPGAs | Thousands of multiply-adds (8,192 per chip in one account) |
| **Sparse p-bit Ising machine** ([Aadit et al., Nature Electronics 2022](https://arxiv.org/abs/2110.02481)) | FPGA | Sparse, any J | Thousands per FPGA; about a million on 18 FPGAs ([thesis](https://www.ece.ucsb.edu/events/all/2025/oct-23-thu-230pm-extreme-scale-ising-machines-distributed-p-computers-navid-anjum)) | Every p-bit of one colour group updated at once; about 10⁶ times faster than Gibbs sampling on a CPU |
| **Hitachi CMOS annealer** | Custom chip (an FPGA version also exists) | Sparse, up to 5 neighbours | About 20,000 | All spins in parallel |

- **The spin counts aren't comparable.** This design reaches 49,000 spins
  only because rank-1 coupling needs no stored J matrix. The other machines
  store a general J, and that storage is what limits them. With a dense
  16-bit J, the Cyclone 10 GX's 11.7 Mb of block RAM would hold only about
  900–1,200 spins.
- **They do far more work per clock.** This design pushes one spin through
  one shared pipeline per clock. The other machines run hundreds to
  thousands of units in parallel, so they do roughly 100 to 10,000 times
  more work per clock, usually on larger FPGAs.
- **Against a CPU, there is no gain for now.** Rank-1 coupling makes each
  layer O(N) on a CPU too. Measured with `simulator/`, an Apple M3 runs the
  same 49,152-spin, 2-layer case bit for bit in 0.80 ms on one core and
  0.32 ms on all 8, against about 0.49 ms for the FPGA. The other machines win because their problems
  need O(N²) work per step for dense J, or O(edges) for sparse J, and that
  is where parallel hardware pays off.
- **Solution quality is a separate question.** MF-AOA uses a few
  deterministic layers, while simulated bifurcation and p-bit machines
  iterate thousands of steps. We haven't found a published comparison of
  MF-AOA against either on solution quality.

What would close the gap:

1. **More pipelines.** Logic is only about 7% used, so about 8 copies of the
   pipeline would fit. That gives 8 times the throughput, or evaluates 8
   (γ, β) settings at once, which helps because MF-AOA needs γ and β
   optimized.
2. **General sparse J,** with many parallel multiply-add units. This is
   where FPGAs beat CPUs.
3. **A host link,** to load γ, β and the weights and read back the final
   spins.

For an overview of the field, see the review [Ising machines as hardware
solvers of combinatorial optimization problems](https://arxiv.org/abs/2204.00276).

## Simulation

`sh sim/run_sims.sh` (needs Icarus Verilog and Python 3) runs these checks:

- **`tb_compare_original.v`:** `VIRAL_n_spin` against `VIRAL_two_spin_2`
  (bit for bit) and `VIRAL_three_spin` (within 4 LSB, for the field rounding
  difference above).
- **`tb_VIRAL_n_spin.v`:** 1, 4, 10, 100 and 1,000 spins, with random and
  alternating initial states, and up to 4 layers. It checks:
  - final spins against N copies of the original `VIRAL_rotation_2`, bit
    for bit;
  - the energy against a pair-by-pair sum, exactly;
  - the LEDs;
  - spins against a floating-point model (max error 0.003 up to 100 spins,
    0.013 at 1,000, because the rotation angle is proportional to the field,
    which grows with N);
  - the full 32-bit γ: `GAMMA_LO` sets the low 16 bits of every γ, since
    the entries in `angle_lut` use only the upper 16;
  - spin weights (`WEIGHT_FILE`), random in (0, 1], at 10, 100 and 1,000
    spins, and a file of all 1.0, which must give the same results as no
    file;
  - that a reset reloads the initial state.
- **`tb_direct_vs_new.v`:** `synth/VIRAL_n_spin_direct.v` against
  `VIRAL_n_spin` at 3, 4, 10 and 100 spins: same final spins, plus clock
  counts.

## CPU simulator and web page

`simulator/` is a Rust program (no dependencies) that runs the same
computation on a CPU, and `web/` is a page that shows it:

```sh
cd simulator && cargo build --release
./target/release/viral-sim serve      # open http://127.0.0.1:8080/
./target/release/viral-sim bench      # this CPU against the FPGA estimates
sh check_hw.sh                        # Rust against the Verilog, bit for bit
```

- **A bit-exact model of `VIRAL_n_spin`.** It gives the same final spins,
  energy, LEDs and clock count as the Verilog simulation. `check_hw.sh`
  checks 11 cases, up to 3,000 spins with weights and the full 32-bit γ.
- **A floating-point model for six problems:**
  - number partitioning (the hardware's rank-1 coupling);
  - the Sherrington-Kirkpatrick spin glass;
  - Max-Cut on random 3-regular graphs;
  - the long-range Ising chain of trapped-ion simulators;
  - a Hopfield associative memory.
- **Classical baselines** to judge each result: exact optima, Karmarkar-Karp,
  local search, and known limits.
- **Studies** that take many runs:
  - SK energy extrapolated to N → ∞;
  - a sweep over layers and dt;
  - quantum spin squeezing from many classical runs (the discrete truncated
    Wigner method), checked against the exact result.
- **The page** shows the problem, the spins on the Bloch sphere layer by
  layer, the energy and the problem's own measure per layer, the result
  against the baselines, and the studies as they run.

What it shows:

- **SK spin glass:** with 500 layers, mean field extrapolates to −0.7543
  per spin as N → ∞, 1.2% above the ground state of −0.7632.
- **Max-Cut:** mean field cuts about 90% of the edges, against 69.24% for
  QAOA with one layer.
- **Long-range chain:** a slow anneal reaches perfect alternation; fast
  ones leave defects.
- **Hopfield memory:** annealing recalls a pattern only below a load of
  about 0.05, as theory predicts.
- **Squeezing:** averaging classical runs reproduces the quantum squeezing
  to within 0.2 dB.
- **Number partitioning:** hard for mean field; Karmarkar-Karp does far
  better.
- **Speed, rank 1:** an Apple M3 matches the Cyclone 10 GX estimate.
- **Speed, dense J:** the proposed Agilex 7 design would be about 25 times
  faster than the M3, at 8,700 spins.

Details are in [simulator/README.md](simulator/README.md).

## Synthesis comparison

`sh synth/run_synth.sh` (needs Yosys and Python 3) produces the resource
table above. Results go in `synth/build/results.md`. You can choose
configurations, e.g. `sh synth/run_synth.sh new:100 direct:10`. The original
architecture is slow to synthesize at large N: about 9 minutes at 50 spins,
and much longer at 100.

## Notes

- **Coupling:** `VIRAL_n_spin` supports rank-1 coupling J<sub>ij</sub> = a<sub>i</sub> a<sub>j</sub>
  (see [Limitations](#limitations)). For a general coupling matrix J<sub>ij</sub>,
  you could store J in block RAM and compute h<sub>i</sub> = Σ J<sub>ij</sub> z<sub>j</sub> with a few
  multiply-accumulate units.
- **Icarus Verilog and the original files:** Icarus rejects the original
  files as they are, because some regs are used before they are declared.
  Quartus and Yosys accept this. `sim/run_sims.sh` simulates copies with the
  declarations moved up (`sim/hoist_decls.py`).
- **Limitations of the original designs that `VIRAL_n_spin` doesn't have:**
  - `VIRAL_two_spin_2` doesn't feed results back between layers, and reads
    `angle_lut` one clock too early, so it only works for one layer.
  - In `VIRAL_three_spin`, spin 3 starts as the zero vector, so it never
    contributes.
  - Reset doesn't restore the initial spins.
  - The LEDs only show energies between −1 and 1 before wrapping.
