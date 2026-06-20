#!/usr/bin/env bash
# Recompute the manuscript's STATIC headline metrics from the actual tree in one
# deterministic pass, so docs/METRICS.md and the manuscript never drift (the
# internal red-team flagged 218k-vs-223k LOC and ×2.05/×2.17/×2.2 call-ratio
# inconsistencies). Performance ratios are NOT here — those come from the bench
# harness (`go run ./pipeline/cmd/full-validation` / `pipeline/bench`) on a
# pinned machine; this script covers only the counts that must be internally
# consistent across all docs.
#
# Usage:  scripts/recompute-metrics.sh            # print a markdown table
#         scripts/recompute-metrics.sh > out.md   # capture
#
# Run from the repo root. Uses `git ls-files`, so it only counts tracked files
# and ignores build artifacts.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

# --- helpers ---------------------------------------------------------------
# Sum line counts over a NUL-safe file list read on stdin (one path per line).
sum_lines() {
  local total=0 n
  while IFS= read -r f; do
    [ -f "$f" ] || continue
    n=$(wc -l < "$f")
    total=$((total + n))
  done
  echo "$total"
}
count() { wc -l | tr -d ' '; }

# --- Go code ---------------------------------------------------------------
go_nontest_files=$(git ls-files '*.go' | grep -v '_test\.go$' | grep -v '^reference_code/' || true)
go_test_files=$(git ls-files '*_test.go' | grep -v '^reference_code/' || true)

go_nontest_loc=$(printf '%s\n' "$go_nontest_files" | sed '/^$/d' | sum_lines)
go_test_loc=$(printf '%s\n' "$go_test_files" | sed '/^$/d' | sum_lines)
go_total_loc=$((go_nontest_loc + go_test_loc))

# --- Surface area ----------------------------------------------------------
clis=$(git ls-files 'tools/*/cmd/*/main.go' | count)
tool_dirs=$(git ls-files 'tools/*' | sed 's#^tools/\([^/]*\)/.*#\1#' | sort -u | grep -v '\.md$' | count)
test_files=$(printf '%s\n' "$go_test_files" | sed '/^$/d' | count)

# --- Test / benchmark / fuzz function counts -------------------------------
# Counted across the whole tracked tree (tools + pkg + pipeline), excluding submodules.
test_fns=$(git grep -hE '^func (Test|Example)[A-Z_]' -- '*_test.go' ':!reference_code/*' 2>/dev/null | count || echo 0)
bench_fns=$(git grep -hE '^func Benchmark[A-Z_]' -- '*_test.go' ':!reference_code/*' 2>/dev/null | count || echo 0)
fuzz_fns=$(git grep -hE '^func Fuzz[A-Z_]' -- '*_test.go' ':!reference_code/*' 2>/dev/null | count || echo 0)

# --- Upstream corpus (only if submodules are initialized) ------------------
upstream_loc="(submodules not initialized — run: git submodule update --init reference_code/<tool>)"
if [ -e reference_code/htslib/hts.c ] || [ -e reference_code/bcftools/main.c ]; then
  upstream_loc=$(find reference_code \
      \( -name '*.c' -o -name '*.cpp' -o -name '*.cc' -o -name '*.h' -o -name '*.hpp' -o -name '*.pl' -o -name '*.pm' \) \
      -not -path '*/test/*' 2>/dev/null \
    | sed '/^$/d' | sum_lines)
fi

# --- Dependency pins -------------------------------------------------------
deps=$(grep -E '^\t?(require|[[:space:]])?[a-z].*v[0-9]' go.mod | grep -oE '[a-z./0-9-]+ v[0-9][^ /]*' | sort -u | tr '\n' ';' | sed 's/;/; /g')

# --- Provenance ------------------------------------------------------------
head_commit=$(git rev-parse --short HEAD)
gen_date=$(date -u +%Y-%m-%dT%H:%M:%SZ)

# --- Output ----------------------------------------------------------------
cat <<EOF
# Recomputed static metrics

Generated $gen_date from \`$head_commit\` by \`scripts/recompute-metrics.sh\`.
These are the canonical counts; copy them verbatim into docs/METRICS.md and the
manuscript so the numbers never disagree. (Performance ratios are separate —
they come from the bench harness on a pinned machine.)

| metric | value | how |
|---|---|---|
| Go non-test LOC | **$go_nontest_loc** | \`git ls-files '*.go' \| grep -v _test\` |
| Go test LOC | **$go_test_loc** | \`git ls-files '*_test.go'\` |
| Go total LOC | **$go_total_loc** | sum |
| Drop-in CLI binaries | **$clis** | \`tools/*/cmd/*/main.go\` |
| Tool directories | **$tool_dirs** | distinct \`tools/<dir>/\` (≈ CLIs; the "13 families" is a conceptual grouping — QC/format, htslib-core, bedtools-surface — not directory-derived) |
| Test files | **$test_files** | \`*_test.go\` |
| Test/Example functions | **$test_fns** | \`func Test*\`/\`func Example*\` |
| Benchmark functions | **$bench_fns** | \`func Benchmark*\` |
| Fuzz functions | **$fuzz_fns** | \`func Fuzz*\` |
| Upstream C/C++/Perl LOC | **$upstream_loc** | \`find reference_code … *.c/*.cpp/*.h/*.pl\` (excl. test/) |
| Sanctioned deps | $deps | go.mod |

> Caveat (carry into the manuscript): the upstream LOC is *size-of-corpus*, NOT a
> 1:1 reimplementation ratio — it includes unexercised paths and the shared
> htslib library only partially reimplemented. Do not present LOC ratio as effort
> ratio.
EOF
