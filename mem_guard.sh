# Memory protection for pgenfit3 scripts (sourced, not executed).
#
# A single simulation.cc process can need several GB (it keeps every ion,
# beta and neutron hit in memory until the end of the run), and a fit of a
# large real dataset with NumCPU forks one copy per process. Running many at
# once exhausted this server's memory once, so every script:
#   1. computes a budget = 50% of MemTotal (PGENFIT3_MEM_BUDGET_GB may lower
#      it, never raise it);
#   2. caps each process with `ulimit -v` (a process over its cap fails with
#      bad_alloc instead of pushing the machine into swap/OOM);
#   3. runs a watchdog that sums the resident memory (RSS) of all its
#      descendant processes every 2 s and kills them all if the sum exceeds
#      the budget.

# budget in kB
mem_budget_kb() {
    local total_kb half_kb user_kb
    total_kb=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
    half_kb=$(( total_kb / 2 ))
    if [ -n "${PGENFIT3_MEM_BUDGET_GB:-}" ]; then
        user_kb=$(awk -v g="$PGENFIT3_MEM_BUDGET_GB" 'BEGIN{printf "%d", g*1024*1024}')
        if [ "$user_kb" -lt "$half_kb" ]; then echo "$user_kb"; return; fi
    fi
    echo "$half_kb"
}

# total RSS (kB) of all descendants of pid $1
_descendant_rss_kb() {
    ps -eo pid=,ppid=,rss= | awk -v root="$1" '
        { ppid[$1]=$2; rss[$1]=$3 }
        END {
            total=0
            for (p in ppid) {
                q=p; n=0
                while (q!="" && q!=root && q!=1 && q!=0 && n<64) { q=ppid[q]; n++ }
                if (q==root && p!=root) total+=rss[p]
            }
            print total
        }'
}

# all descendant pids of $1, deepest first
_descendant_pids() {
    ps -eo pid=,ppid= | awk -v root="$1" '
        { ppid[$1]=$2 }
        END {
            for (p in ppid) {
                q=p; d=0
                while (q!="" && q!=root && q!=1 && q!=0 && d<64) { q=ppid[q]; d++ }
                if (q==root && p!=root) print d, p
            }
        }' | sort -rn | awk '{print $2}'
}

# start_mem_watchdog <budget_kb>: watches the calling script's descendants.
# On overflow it prints an error, kills every descendant, and leaves a flag
# file whose path is in $MEM_WATCHDOG_FLAG. Stop it with stop_mem_watchdog.
start_mem_watchdog() {
    local budget_kb=$1 root=$$
    MEM_WATCHDOG_FLAG=$(mktemp "${TMPDIR:-/tmp}/pgenfit3_memflag.XXXXXX")
    rm -f "$MEM_WATCHDOG_FLAG"
    (
        while kill -0 "$root" 2>/dev/null; do
            used=$(_descendant_rss_kb "$root")
            if [ "$used" -gt "$budget_kb" ]; then
                echo "MEMORY WATCHDOG: descendants use $((used/1024)) MB > budget $((budget_kb/1024)) MB -- killing all jobs" >&2
                touch "$MEM_WATCHDOG_FLAG"
                for p in $(_descendant_pids "$root"); do
                    [ "$p" != "$BASHPID" ] && kill -TERM "$p" 2>/dev/null
                done
                sleep 2
                for p in $(_descendant_pids "$root"); do
                    [ "$p" != "$BASHPID" ] && kill -KILL "$p" 2>/dev/null
                done
                exit 0
            fi
            sleep 2
        done
    ) &
    MEM_WATCHDOG_PID=$!
}

stop_mem_watchdog() {
    [ -n "${MEM_WATCHDOG_PID:-}" ] && kill "$MEM_WATCHDOG_PID" 2>/dev/null
    wait "$MEM_WATCHDOG_PID" 2>/dev/null || true
    if [ -n "${MEM_WATCHDOG_FLAG:-}" ] && [ -f "$MEM_WATCHDOG_FLAG" ]; then
        rm -f "$MEM_WATCHDOG_FLAG"
        return 1
    fi
    return 0
}
