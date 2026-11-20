#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir"

size=${SIZE:-150}
read -r -a block_counts <<< "${BLOCK_COUNTS:-1 5 10 20 50 100}"
out=perf.nogit
summary="$out/summary-$size.csv"
mkdir -p "$out"
printf 'blocks,time_s,peak_rss_mib\n' > "$summary"

profiles=()
for blocks in "${block_counts[@]}"; do
    dat="$out/block-${size}-${blocks}.dat"
    log="$out/block-${size}-${blocks}.log"
    err="$out/block-${size}-${blocks}.mprof.err"
    profiles+=("$dat")

    echo "Running ${blocks} block(s)..."
    # mprof 0.61 samples before checking Popen.poll(). On macOS the child can
    # finish between those operations, producing a final psutil.NoSuchProcess
    # even though both the benchmark and profile completed successfully.
    # Truncate the append-only profile first so validation cannot accept data
    # left by an earlier run.
    : > "$dat"
    set +e
    PYENV_VERSION=3.9.19 MPLCONFIGDIR=/tmp/cfdtools-mpl \
        mprof run --python -T 0.05 -o "$dat" \
        python cube_connectivity.py \
        --size "$size" --implementation multi \
        --blocks "$blocks" --quiet-blocks > "$log" 2> "$err"
    mprof_status=$?
    set -e

    if ((mprof_status != 0)); then
        if grep -q 'psutil.NoSuchProcess' "$err" \
            && grep -q '^  multiblock (' "$log" \
            && grep -q '^MEM ' "$dat"; then
            echo "  mprof final PID race ignored; benchmark/profile are complete"
        else
            echo "mprof or benchmark failed for $blocks block(s)" >&2
            cat "$err" >&2
            exit "$mprof_status"
        fi
    fi
    if ! grep -q '^  multiblock (' "$log" || ! grep -q '^MEM ' "$dat"; then
        echo "Incomplete benchmark/profile for $blocks block(s)" >&2
        cat "$log" >&2
        cat "$err" >&2
        exit 1
    fi

    seconds=$(awk '/^  multiblock \(/{print $4; exit}' "$log")
    peak=$(PYENV_VERSION=3.9.19 mprof peak "$dat" | awk '{print $(NF-1)}')
    printf '%s,%s,%s\n' "$blocks" "$seconds" "$peak" |
        tee -a "$summary"
done

PYENV_VERSION=3.9.19 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/cfdtools-mpl \
    mprof plot "${profiles[@]}" \
    --flame \
    --title "Multiblock connectivity: RSS over time" \
    -o "$out/rss-profiles-$size.png"

PYENV_VERSION=3.9.19 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/cfdtools-mpl python - \
    "$summary" "$out/time-rss-summary-$size.png" <<'PY'
import sys
import numpy as np
import matplotlib.pyplot as plt

summary, output = sys.argv[1:]
data = np.genfromtxt(summary, delimiter=",", names=True)

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7, 7))
axes[0].plot(data["blocks"], data["time_s"], "o-")
axes[0].set_ylabel("Connectivity time [s]")
axes[0].grid(True)

axes[1].plot(data["blocks"], data["peak_rss_mib"], "o-", color="tab:red")
axes[1].set_xlabel("Number of blocks")
axes[1].set_ylabel("Peak RSS [MiB]")
axes[1].grid(True)

fig.tight_layout()
fig.savefig(output, dpi=160)
PY
