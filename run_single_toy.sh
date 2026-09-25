#!/usr/bin/env bash
# Simulate one toy with simulation.cc and fit its total decay curve.
#
# Usage: ./run_single_toy.sh <parmsex_file> <simparmsex_file> <seed> [effparms|none]
#            [boundaryMarginSec=900] [ncpu=4] [out_prefix=toy] [startTime=0] [timeRange=10]
# Output: single_toy_results/<out_prefix>_seed<seed>.{root,log} and _fit.png/_fit.root
# Memory: limited to 50% of RAM (ulimit per process + watchdog on the
# total, which also covers NumCPU forks) -- see mem_guard.sh.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"
source "$SCRIPT_DIR/mem_guard.sh"

if [ $# -lt 3 ]; then
    echo "Usage: $0 <parmsex_file> <simparmsex_file> <seed> [effparms|none] [boundaryMarginSec=900] [ncpu=4] [out_prefix=toy] [startTime=0] [timeRange=10]" >&2
    exit 1
fi
abspath(){ echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
PARMSEX=$(abspath "$1")
SIMPARMS=$(abspath "$2")
SEED=$3
EFF=${4:-none}
[ "$EFF" != "none" ] && EFF=$(abspath "$EFF")
BOUNDARY=${5:-900}
NCPU=${6:-4}
PREFIX=${7:-toy}
START_TIME=${8:-0}
TIME_RANGE=${9:-10}

RESULTDIR="$SCRIPT_DIR/single_toy_results"
mkdir -p "$RESULTDIR" "$SCRIPT_DIR/builds"
BUILD_DIR=$(mktemp -d "$SCRIPT_DIR/builds/toy.XXXXXX")
trap 'rm -rf "$BUILD_DIR"' EXIT
"$SCRIPT_DIR/build.sh" "$PARMSEX" "$BUILD_DIR" withsim
cd "$BUILD_DIR"
OUT="$RESULTDIR/${PREFIX}_seed${SEED}"
BUDGET_KB=$(mem_budget_kb)
ulimit -v "$BUDGET_KB"
start_mem_watchdog "$BUDGET_KB"
./simulation "$PARMSEX" "$SIMPARMS" "$OUT.root" "$SEED" > sim.log 2>&1 || { tail sim.log >&2; exit 1; }
./fit_decay_curve "$OUT.root" "$PARMSEX" "$NCPU" "$EFF" "$OUT" "$BOUNDARY" "$START_TIME" "$TIME_RANGE" \
    2>&1 | tee "$OUT.log"
if ! stop_mem_watchdog; then echo "ABORTED: memory budget ($((BUDGET_KB/1024/1024)) GB) exceeded" >&2; exit 1; fi
echo "toy: $OUT.root  log: $OUT.log  plots: ${OUT}_fit.png"
