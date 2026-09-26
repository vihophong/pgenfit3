#!/usr/bin/env bash
# Simulate one toy of a MIXTURE of implanted species with simulation_mix and
# fit it with fit_decay_curve (parmsex list = mixture fit, see README).
#
# Usage: ./run_mix_toy.sh <parmsexA,parmsexB[,...]> <simparmsA,simparmsB[,...]> <seed>
#            [effparms|none] [boundaryMarginSec=900] [ncpu=4] [out_prefix=mix] [startTime=0] [timeRange=10]
# Set ALPHA_GATE=1 for the alpha-gated fit.
# Output: single_toy_results/<out_prefix>_seed<seed>.{root,log} and _fit.png/_fit.root
# Memory: limited to 50% of RAM (see mem_guard.sh).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"
source "$SCRIPT_DIR/mem_guard.sh"

if [ $# -lt 3 ]; then
    echo "Usage: $0 <parmsexA,parmsexB[,...]> <simparmsA,simparmsB[,...]> <seed> [effparms|none] [boundaryMarginSec=900] [ncpu=4] [out_prefix=mix] [startTime=0] [timeRange=10]" >&2
    exit 1
fi
abspath(){ echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
IFS=',' read -r -a PL <<< "$1"
IFS=',' read -r -a SL <<< "$2"
[ ${#PL[@]} -eq ${#SL[@]} ] || { echo "need one simparms per parmsex" >&2; exit 1; }
PARMS=(); SIMARGS=()
for i in "${!PL[@]}"; do
    PARMS+=("$(abspath "${PL[$i]}")")
    SIMARGS+=("$(abspath "${PL[$i]}")" "$(abspath "${SL[$i]}")")
done
PARMLIST=$(IFS=','; echo "${PARMS[*]}")
SEED=$3
EFF=${4:-none}
[ "$EFF" != "none" ] && EFF=$(abspath "$EFF")
BOUNDARY=${5:-900}
NCPU=${6:-4}
PREFIX=${7:-mix}
START_TIME=${8:-0}
TIME_RANGE=${9:-10}

RESULTDIR="$SCRIPT_DIR/single_toy_results"
mkdir -p "$RESULTDIR" "$SCRIPT_DIR/builds"
BUILD_DIR=$(mktemp -d "$SCRIPT_DIR/builds/mix.XXXXXX")
trap 'rm -rf "$BUILD_DIR"' EXIT
"$SCRIPT_DIR/build.sh" "$PARMLIST" "$BUILD_DIR" withsim
cd "$BUILD_DIR"
OUT="$RESULTDIR/${PREFIX}_seed${SEED}"
BUDGET_KB=$(mem_budget_kb)
ulimit -v "$BUDGET_KB"
start_mem_watchdog "$BUDGET_KB"
./simulation_mix "$OUT.root" "$SEED" "${SIMARGS[@]}" > sim.log 2>&1 || { tail sim.log >&2; exit 1; }
./fit_decay_curve "$OUT.root" "$PARMLIST" "$NCPU" "$EFF" "$OUT" "$BOUNDARY" "$START_TIME" "$TIME_RANGE" \
    2>&1 | tee "$OUT.log"
if ! stop_mem_watchdog; then echo "ABORTED: memory budget ($((BUDGET_KB/1024/1024)) GB) exceeded" >&2; exit 1; fi
echo "toy: $OUT.root  log: $OUT.log  plots: ${OUT}_fit.png"
