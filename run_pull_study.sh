#!/usr/bin/env bash
# Toy-MC pull study: simulate nToys independent datasets with simulation.cc
# and fit each one with fit_decay_curve. One PULLDATA line per toy is
# collected in pull_study_results/<prefix>_all.txt; analyse it with
# ./plot_pulls.py.
#
# Usage: ./run_pull_study.sh <parmsex_file> <simparmsex_file> <startSeed> <nToys>
#            [effparms|none] [nWorkers=4] [boundaryMarginSec=900] [ncpu_per_fit=1]
#            [prefix=pull] [startTime=0] [timeRange=10]
#
# Memory: the whole study is limited to 50% of this machine's RAM (see
# mem_guard.sh). Each worker is capped at MEM_PER_WORKER_GB (default 8; a
# simulation with the As92 simparmsex.txt beam settings keeps millions of
# hits in memory). The worker count is reduced so that
# nWorkers*MEM_PER_WORKER_GB fits the budget, and a watchdog kills
# everything if the summed RSS still exceeds it. Measure one toy first
# (/usr/bin/time -v on ./simulation) and set MEM_PER_WORKER_GB from that.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"
source "$SCRIPT_DIR/mem_guard.sh"

if [ $# -lt 4 ]; then
    echo "Usage: $0 <parmsex_file> <simparmsex_file> <startSeed> <nToys> [effparms|none] [nWorkers=4] [boundaryMarginSec=900] [ncpu_per_fit=1] [prefix=pull] [startTime=0] [timeRange=10]" >&2
    exit 1
fi
abspath(){ echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
PARMSEX=$(abspath "$1")
SIMPARMS=$(abspath "$2")
START_SEED=$3
NTOYS=$4
EFF=${5:-none}
[ "$EFF" != "none" ] && EFF=$(abspath "$EFF")
NWORKERS=${6:-4}
BOUNDARY=${7:-900}
NCPU=${8:-1}
PREFIX=${9:-pull}
START_TIME=${10:-0}
TIME_RANGE=${11:-10}

MEM_PER_WORKER_GB=${MEM_PER_WORKER_GB:-8}
BUDGET_KB=$(mem_budget_kb)
PER_WORKER_KB=$(awk -v g="$MEM_PER_WORKER_GB" 'BEGIN{printf "%d", g*1024*1024}')
MAX_WORKERS=$(( BUDGET_KB / PER_WORKER_KB ))
if [ "$MAX_WORKERS" -lt 1 ]; then
    echo "MEM_PER_WORKER_GB=$MEM_PER_WORKER_GB exceeds the memory budget ($((BUDGET_KB/1024/1024)) GB)" >&2
    exit 1
fi
if [ "$NWORKERS" -gt "$MAX_WORKERS" ]; then
    echo "==> reducing workers $NWORKERS -> $MAX_WORKERS (budget $((BUDGET_KB/1024/1024)) GB, ${MEM_PER_WORKER_GB} GB/worker)"
    NWORKERS=$MAX_WORKERS
fi

RESULTDIR="$SCRIPT_DIR/pull_study_results"
mkdir -p "$RESULTDIR" "$SCRIPT_DIR/builds"
BUILD_DIR=$(mktemp -d "$SCRIPT_DIR/builds/pull.XXXXXX")
trap 'rm -rf "$BUILD_DIR"' EXIT
"$SCRIPT_DIR/build.sh" "$PARMSEX" "$BUILD_DIR" withsim

PER_WORKER=$(( (NTOYS + NWORKERS - 1) / NWORKERS ))
{
    echo "# parmsex=$PARMSEX simparmsex=$SIMPARMS effparms=$EFF"
    echo "# startSeed=$START_SEED nToys=$NTOYS boundaryMarginSec=$BOUNDARY startTime=$START_TIME timeRange=$TIME_RANGE"
} > "$RESULTDIR/${PREFIX}_config.txt"
echo "==> $NTOYS toys, $NWORKERS workers x $PER_WORKER, results in $RESULTDIR/${PREFIX}_all.txt"

run_worker() {
    local w=$1
    ulimit -v "$PER_WORKER_KB"   # hard cap for this worker's simulation/fit
    local wdir; wdir=$(mktemp -d "$BUILD_DIR/w${w}.XXXX")
    local out="$RESULTDIR/${PREFIX}_worker_${w}.txt"
    : > "$out"
    local i seed line
    for ((i=0; i<PER_WORKER; i++)); do
        local n=$(( w*PER_WORKER + i ))
        [ $n -ge "$NTOYS" ] && break
        seed=$(( START_SEED + n ))
        if ! "$BUILD_DIR/simulation" "$PARMSEX" "$SIMPARMS" "$wdir/toy.root" "$seed" > "$wdir/sim.log" 2>&1; then
            echo "worker $w toy $i: simulation FAILED, seed=$seed $(grep -m1 -o 'bad_alloc' "$wdir/sim.log" || true)" >> "$out"; continue
        fi
        line=$(cd "$BUILD_DIR" && ./fit_decay_curve "$wdir/toy.root" "$PARMSEX" "$NCPU" "$EFF" noplot \
                   "$BOUNDARY" "$START_TIME" "$TIME_RANGE" 2>/dev/null | grep "^PULLDATA" || true)
        if [ -z "$line" ]; then
            echo "worker $w toy $i: no PULLDATA, seed=$seed" >> "$out"
        else
            echo "worker $w toy $i seed=$seed $line" >> "$out"
        fi
        rm -f "$wdir/toy.root"
    done
}
echo "==> memory budget $((BUDGET_KB/1024/1024)) GB, $NWORKERS workers capped at ${MEM_PER_WORKER_GB} GB each"
start_mem_watchdog "$BUDGET_KB"
for ((w=0; w<NWORKERS; w++)); do run_worker "$w" & done
WORKER_PIDS=$(jobs -p | grep -v "^$MEM_WATCHDOG_PID\$" || true)
for p in $WORKER_PIDS; do wait "$p" || true; done
if ! stop_mem_watchdog; then
    echo "ABORTED: memory budget exceeded; partial results kept in $RESULTDIR/${PREFIX}_worker_*.txt" >&2
    exit 1
fi
cat "$RESULTDIR/${PREFIX}_worker_"*.txt > "$RESULTDIR/${PREFIX}_all.txt"
rm -f "$RESULTDIR/${PREFIX}_worker_"*.txt
NOK=$(grep -c "PULLDATA l0=" "$RESULTDIR/${PREFIX}_all.txt" || true)
echo "==> done: $NOK/$NTOYS toys fitted -> $RESULTDIR/${PREFIX}_all.txt"
[ "$NOK" -eq "$NTOYS" ] || echo "WARNING: $((NTOYS-NOK)) toy(s) failed, see lines without PULLDATA" >&2
