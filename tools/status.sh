#!/usr/bin/env bash
# What is running on this machine, in one screen (24 Nov 2026):
#   1. tmux sessions and the command running in each
#   2. who holds each GPU (process, age, and the tmux session it lives in, if any)
#   3. queue runners and the run queue of the PanArcticShenanigans clone (queued / running / last done / failed)
#   4. the newest results in les_under_ice/output
# Read-only: it changes nothing.
#
#     bash /mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-tools/tools/status.sh
#     bash .../status.sh /path/to/other/panarctic-clone        # another clone's queue

PAS=${1:-/mnt/autofs/bow.math/scratch/fpoulin/software/panarctic-061}
TOOLS=$(cd "$(dirname "$0")/.." && pwd)

# tmux pane process → session name, and a process's session by walking up its parents
declare -A pane_session
while read -r session pid; do
    pane_session[$pid]=$session
done < <(tmux list-panes -a -F '#{session_name} #{pane_pid}' 2>/dev/null)

session_of() {
    local p=$1
    while [ -n "$p" ] && [ "$p" -gt 1 ]; do
        if [ -n "${pane_session[$p]}" ]; then
            echo "${pane_session[$p]}"
            return
        fi
        p=$(ps -o ppid= -p "$p" 2>/dev/null | tr -d ' ')
    done
    echo "-"
}

echo "== tmux sessions =="
tmux list-panes -a -F '#{session_name}|#{pane_current_command}|#{pane_current_path}' 2>/dev/null |
    awk -F'|' '{ printf "  %-22s %-10s %s\n", $1, $2, $3 }'
[ ${#pane_session[@]} -eq 0 ] && echo "  (none)"

echo
echo "== GPUs (index: process, memory, age, tmux session) =="
declare -A gpu_index
while IFS=', ' read -r index uuid; do
    gpu_index[$uuid]=$index
done < <(nvidia-smi --query-gpu=index,uuid --format=csv,noheader)
busy=$(nvidia-smi --query-compute-apps=gpu_uuid,pid,used_memory --format=csv,noheader)
for index in $(nvidia-smi --query-gpu=index --format=csv,noheader); do
    line=""
    while IFS=', ' read -r uuid pid memory unit; do
        [ -z "$uuid" ] && continue
        [ "${gpu_index[$uuid]}" != "$index" ] && continue
        user=$(ps -o user= -p "$pid" 2>/dev/null | tr -d ' ')
        age=$(ps -o etime= -p "$pid" 2>/dev/null | tr -d ' ')
        line+=$(printf "pid %s (%s, %s %s, running %s, tmux: %s)  " "$pid" "$user" "$memory" "$unit" "$age" "$(session_of "$pid")")
    done <<< "$busy"
    echo "  GPU $index: ${line:-free}"
done

echo
echo "== Queue runners =="
pgrep -af "runner_loop" | grep -v pgrep | sed 's/^/  /' || true
pgrep -f runner_loop > /dev/null || echo "  (none running)"
for log in ~/logs/runner-*.log; do
    [ -f "$log" ] && echo "  $(basename "$log"): $(tail -n 1 "$log")"
done

echo
echo "== Run queue in $PAS/queue =="
for state in running queue failed; do
    files=$(ls "$PAS/queue/$state" 2>/dev/null | grep -v '\.claim$')     # .claim: the runner's lock file
    echo "  $state ($(echo -n "$files" | grep -c .)):"
    [ -n "$files" ] && echo "$files" | sed 's/^/    /'
done
echo "  done, last 5:"
ls -t "$PAS/queue/done" 2>/dev/null | head -5 | sed 's/^/    /'

echo
echo "== Newest results in les_under_ice/output =="
ls -lt --time-style=+'%a %H:%M' "$TOOLS/les_under_ice/output" 2>/dev/null | grep -E '\.(md|png)$' | head -8 |
    awk '{ printf "  %s %s  %s\n", $6, $7, $8 }'
