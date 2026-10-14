#!/bin/sh
# Checks the Rust model (viral-sim hw) against the Verilog simulation of
# VIRAL_n_spin bit for bit: the same initial states, weights and angle_lut
# layers through Icarus Verilog and through Rust must give identical final
# spins, energy, LEDs and clock count.
#
# Usage: sh simulator/check_hw.sh   (needs iverilog, python3 and cargo)

set -e
cd "$(dirname "$0")"
cargo build --release --quiet
SIM=../sim
OUT=target/check
mkdir -p "$OUT"

for f in VIRAL_rotation_2; do
	python3 $SIM/hoist_decls.py "../$f.v" "$OUT/$f.v"
done
SRC="../angle_lut.v ../cordic_iterative.v ../cordic_pipelined.v \
     $OUT/VIRAL_rotation_2.v ../VIRAL_rotation_pipe.v ../spin_memory.v \
     ../serial_mult.v ../VIRAL_n_spin.v"

fails=0
check() {
	n=$1; layers=$2; pattern=$3; gamma_lo=$4; weights=$5
	name="${n}_${layers}_${pattern}_${gamma_lo}_${weights:-none}"
	init="$OUT/init_$name.hex"
	python3 ../gen_spin_init.py "$n" --pattern "$pattern" -o "$init"
	wfile=""
	wopt=""
	if [ -n "$weights" ]; then
		wfile="$OUT/weights_$name.hex"
		python3 ../gen_weights.py "$n" --pattern "$weights" -o "$wfile"
		wopt="--weights $wfile"
	fi
	iverilog -g2005 -o "$OUT/tb_$name" \
		-Ptb_VIRAL_n_spin.N_SPINS="$n" \
		-Ptb_VIRAL_n_spin.N_LAYERS="$layers" \
		-Ptb_VIRAL_n_spin.INIT_FILE="\"$init\"" \
		-Ptb_VIRAL_n_spin.FLOAT_TOL=1.0 \
		-Ptb_VIRAL_n_spin.GAMMA_LO="$gamma_lo" \
		-Ptb_VIRAL_n_spin.WEIGHT_FILE="\"$wfile\"" \
		-Ptb_VIRAL_n_spin.DUMP_FILE="\"$OUT/verilog_$name.txt\"" \
		$SRC $SIM/tb_VIRAL_n_spin.v
	vvp -n "$OUT/tb_$name" > "$OUT/vvp_$name.log"
	./target/release/viral-sim hw --init "$init" $wopt --layers "$layers" \
		--gamma-lo "$gamma_lo" --dump "$OUT/rust_$name.txt" > /dev/null
	if cmp -s "$OUT/verilog_$name.txt" "$OUT/rust_$name.txt"; then
		echo "same: $n spins, $layers layers, $pattern, gamma low bits $gamma_lo, weights ${weights:-none}"
	else
		echo "DIFFERENT: $n spins, $layers layers, $pattern, gamma low bits $gamma_lo, weights ${weights:-none}"
		diff "$OUT/verilog_$name.txt" "$OUT/rust_$name.txt" | head -5
		fails=$((fails + 1))
	fi
}

check 1    2 random      0
check 4    2 random      0
check 10   4 alternating 0
check 100  2 random      0
check 1000 2 random      0
check 10   4 alternating 4660
check 100  4 random      43981
check 10   2 random      0     ones
check 100  2 random      0     random
check 1000 2 random      43981 random
check 3000 4 random      4660  random

if [ $fails -eq 0 ]; then
	echo "PASS: Rust and Verilog agree bit for bit"
else
	echo "FAIL: $fails cases differ"
	exit 1
fi
