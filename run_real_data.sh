#!/usr/bin/env bash
# Fit the total decay curve of an existing ROOT file (real data, or any
# file with pgenfit simulation.cc's "tree": branches x, optional ionT).
#
# Usage: ./run_real_data.sh <root_file> <parmsex_file> [effparms|none] [ncpu=8]
#            [boundaryMarginSec=0] [out_prefix] [startTime=0] [timeRange=10]
#            [linBinFactor=4] [noplot]
#   out_prefix  default real_data_results/<root file basename>; the fit log
#               goes to <out_prefix>.log, plots to <out_prefix>_fit.png/.root
# Memory: limited to 50% of RAM (ulimit per process + watchdog on the
# total, which also covers NumCPU forks) -- see mem_guard.sh.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"
source "$SCRIPT_DIR/mem_guard.sh"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <root_file> <parmsex_file> [effparms|none] [ncpu=8] [boundaryMarginSec=0] [out_prefix] [startTime=0] [timeRange=10] [linBinFactor=4] [noplot]" >&2
    exit 1
fi
abspath(){ echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
ROOT_FILE=$(abspath "$1")
PARMSEX=$(abspath "$2")
EFF=${3:-none}
[ "$EFF" != "none" ] && EFF=$(abspath "$EFF")
NCPU=${4:-8}
BOUNDARY=${5:-0}
DEFAULT_PREFIX="real_data_results/$(basename "${ROOT_FILE%.root}")"
OUTPREFIX_IN=${6:-$DEFAULT_PREFIX}
START_TIME=${7:-0}
TIME_RANGE=${8:-10}
LIN_BIN_FACTOR=${9:-4}
NOPLOT=${10:-}

mkdir -p "$(dirname "$OUTPREFIX_IN")"
OUTPREFIX="$(cd "$(dirname "$OUTPREFIX_IN")" && pwd)/$(basename "$OUTPREFIX_IN")"
PLOTARG="$OUTPREFIX"
[ "$NOPLOT" = "noplot" ] && PLOTARG="noplot"

mkdir -p "$SCRIPT_DIR/builds"
BUILD_DIR=$(mktemp -d "$SCRIPT_DIR/builds/real.XXXXXX")
trap 'rm -rf "$BUILD_DIR"' EXIT
"$SCRIPT_DIR/build.sh" "$PARMSEX" "$BUILD_DIR"
cd "$BUILD_DIR"
BUDGET_KB=$(mem_budget_kb)
ulimit -v "$BUDGET_KB"
start_mem_watchdog "$BUDGET_KB"
./fit_decay_curve "$ROOT_FILE" "$PARMSEX" "$NCPU" "$EFF" "$PLOTARG" "$BOUNDARY" "$START_TIME" "$TIME_RANGE" "$LIN_BIN_FACTOR" \
    2>&1 | tee "$OUTPREFIX.log"
if ! stop_mem_watchdog; then echo "ABORTED: memory budget ($((BUDGET_KB/1024/1024)) GB) exceeded" >&2; exit 1; fi
echo "log: $OUTPREFIX.log"
[ "$PLOTARG" != "noplot" ] && echo "plots: ${OUTPREFIX}_fit.png ${OUTPREFIX}_fit.root"
exit 0
