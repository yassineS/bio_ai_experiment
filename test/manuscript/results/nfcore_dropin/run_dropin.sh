#!/usr/bin/env bash
# Drop-in demo: run the UNMODIFIED nf-core samtools/flagstat module twice on the
# same input — once with upstream samtools on PATH, once with our Go samtools on
# PATH — and diff the results. The module file is never edited between runs; only
# PATH changes. Evidence (stdout logs + output .flagstat files) is written under
# evidence/.
set -euo pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$DEMO_DIR/../../../.." && pwd)"
UP_BIN="$REPO_ROOT/bin/upstream"
OURS_BIN="$REPO_ROOT/bin/ours"
EV="$DEMO_DIR/evidence"
mkdir -p "$EV"

command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
[ -x "$UP_BIN/samtools" ]   || { echo "missing $UP_BIN/samtools (run: make upstream)" >&2; exit 1; }
[ -x "$OURS_BIN/samtools" ] || { echo "missing $OURS_BIN/samtools (run: make ours)" >&2; exit 1; }

run_one() {
    local label="$1" binpath="$2"
    echo "=== [$label] samtools = $(PATH="$binpath:$PATH" command -v samtools) ===" | tee "$EV/${label}.stdout"
    echo "=== [$label] samtools version line: $(PATH="$binpath:$PATH" samtools version 2>&1 | head -1) ===" | tee -a "$EV/${label}.stdout"
    rm -rf "$DEMO_DIR/work_${label}" "$DEMO_DIR/results_${label}"
    ( cd "$DEMO_DIR" && PATH="$binpath:$PATH" nextflow -q run main.nf \
        -work-dir "work_${label}" \
        --outdir "results_${label}" 2>&1 ) | tee -a "$EV/${label}.stdout"
    # Collect the module's published/emitted output file.
    local out
    out="$(find "$DEMO_DIR/work_${label}" -name '*.flagstat' | head -1)"
    cp "$out" "$EV/${label}.flagstat"
    echo "--- [$label] test.flagstat ---" ; cat "$EV/${label}.flagstat"
}

run_one upstream "$UP_BIN"
echo
run_one ours "$OURS_BIN"

echo
echo "=== DIFF upstream vs ours (.flagstat) ==="
if diff -u "$EV/upstream.flagstat" "$EV/ours.flagstat"; then
    echo "RESULT: BYTE-IDENTICAL — our samtools is a drop-in for nf-core samtools/flagstat"
else
    echo "RESULT: DIFFERENCES FOUND" >&2
    exit 1
fi
