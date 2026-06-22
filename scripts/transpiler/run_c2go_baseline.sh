#!/usr/bin/env bash
#
# run_c2go_baseline.sh
#
# Manuscript claim C5 (counterfactual / desk-reject control): demonstrate that a
# NON-LLM, rule-based automated C->Go transpiler (c2go) does NOT clear the bar
# the LLM agents cleared (compiles + runs + byte-exact parity + idiomatic and
# memory-safe / cgo-free Go).
#
# This script is a reproducible harness. It installs c2go, attempts to transpile
# a trivial C program and a real vendored bioinformatics C source, and prints an
# evaluation checklist against the four bar criteria. It is written to be run on
# a local machine; it does not require an LLM and makes no network calls beyond
# the `go install` of c2go (and an optional clang install).
#
# Honesty note: on the environment where docs/manuscript/results/transpiler_baseline.md
# was produced, c2go installed but FAILED at the AST stage on every input because
# it does not understand the `BuiltinAttr` node emitted by modern clang (>= 10).
# This script will reproduce that, AND, if you have an old-enough clang (<= 9),
# will let you observe the generated Go directly. Either way the conclusion holds;
# see the results doc for the calibrated argument.
#
# Usage:
#   scripts/transpiler/run_c2go_baseline.sh
#
# Environment overrides:
#   C2GO_BIN   path to a prebuilt c2go (skip `go install`)
#   CLANG_BIN  clang to shadow on PATH (default: `clang`); set to clang-9 for best shot
#
set -uo pipefail

# Resolve repo root from this script's location (scripts/transpiler/..).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

PASS="PASS"
FAIL="FAIL"

log()  { printf '\n=== %s ===\n' "$*"; }
note() { printf '    %s\n' "$*"; }

# --- 0. Tooling -------------------------------------------------------------
log "0. Tooling"
command -v go    >/dev/null 2>&1 || { echo "go not found"; exit 127; }
go version
CLANG_BIN="${CLANG_BIN:-clang}"
if command -v "${CLANG_BIN}" >/dev/null 2>&1; then
  "${CLANG_BIN}" --version | head -1
  note "NOTE: c2go invokes the binary literally named 'clang' from PATH."
  note "      c2go supports clang up to ~9. clang >= 10 emits 'BuiltinAttr'"
  note "      AST nodes that c2go cannot parse (it panics). If your default"
  note "      clang is modern, the transpile step below is EXPECTED to fail."
else
  echo "clang not found -- c2go needs clang to dump the C AST"; exit 127
fi

# --- 1. Install c2go --------------------------------------------------------
log "1. Install c2go (non-LLM rule-based C->Go transpiler)"
if [ -n "${C2GO_BIN:-}" ] && [ -x "${C2GO_BIN}" ]; then
  note "Using preinstalled C2GO_BIN=${C2GO_BIN}"
else
  go install github.com/elliotchance/c2go@latest
  C2GO_BIN="$(go env GOPATH)/bin/c2go"
fi
[ -x "${C2GO_BIN}" ] || { echo "c2go not installed at ${C2GO_BIN}"; exit 1; }
note "c2go binary: ${C2GO_BIN}"

# --- 2. Trivial C program (best case) --------------------------------------
log "2. Trivial input: hello.c"
cat > "${WORK}/hello.c" <<'EOF'
#include <stdio.h>
int main(void) {
    printf("Hello, World!\n");
    return 0;
}
EOF
HELLO_TRANSPILE="${FAIL}"
if "${C2GO_BIN}" transpile -o "${WORK}/hello.go" "${WORK}/hello.c" 2>"${WORK}/hello.err"; then
  HELLO_TRANSPILE="${PASS}"
  note "Transpile produced hello.go:"
  sed 's/^/      /' "${WORK}/hello.go" | head -40
else
  note "Transpile FAILED. First lines of stderr:"
  sed 's/^/      /' "${WORK}/hello.err" | head -8
fi

# --- 3. Real vendored bioinformatics C source ------------------------------
# sickle is pure C (good case); htslib/bgzip.c needs a generated config.h
# (a build-system dependency the transpiler cannot satisfy without ./configure).
log "3. Real vendored C source: sickle/src/sickle.c"
SICKLE_SRC="${REPO_ROOT}/reference_code/sickle/src/sickle.c"
SICKLE_TRANSPILE="${FAIL}"
if [ -f "${SICKLE_SRC}" ]; then
  if "${C2GO_BIN}" transpile -o "${WORK}/sickle.go" "${SICKLE_SRC}" 2>"${WORK}/sickle.err"; then
    SICKLE_TRANSPILE="${PASS}"
    note "Transpile produced sickle.go (head):"
    sed 's/^/      /' "${WORK}/sickle.go" | head -30
  else
    note "Transpile FAILED. First lines of stderr:"
    sed 's/^/      /' "${WORK}/sickle.err" | head -8
  fi
else
  note "sickle.c not found -- run: git submodule update --init reference_code/sickle"
fi

# --- 4. Evaluation against the four bar criteria ---------------------------
# These checks only fire if a transpile actually produced Go. Each is the
# criterion the LLM-produced CLIs met; the transpiler must meet ALL four.
log "4. Evaluation checklist (the bar the LLM agents cleared)"
eval_go() {
  local label="$1" gofile="$2"
  [ -f "${gofile}" ] || { note "[${label}] no Go output -> COMPILES/RUNS/PARITY/IDIOMATIC all FAIL"; return; }

  # (a) compiles
  local c="${FAIL}"
  ( cd "${WORK}" && go vet "${gofile}" >/dev/null 2>&1 ) && c="${PASS}"
  note "[${label}] (a) compiles (go vet): ${c}"

  # (d.1) memory-safe / cgo-free: does the output import "unsafe" or call cgo?
  local unsafe_n cgo_n noarch_n
  unsafe_n=$(grep -c '"unsafe"\|unsafe\.Pointer' "${gofile}" 2>/dev/null || echo 0)
  cgo_n=$(grep -c '"C"\|import "C"' "${gofile}" 2>/dev/null || echo 0)
  noarch_n=$(grep -c 'noarch\.' "${gofile}" 2>/dev/null || echo 0)
  note "[${label}] (d) memory-safe: unsafe refs=${unsafe_n} cgo refs=${cgo_n} libc-shim(noarch) refs=${noarch_n}"
  note "      -> any unsafe>0 or noarch>0 means NON-idiomatic + NOT-from-stdlib (FAIL the bar)"

  # (b) runs / (c) byte-parity are left as manual steps because they require a
  # buildable program AND the original tool to diff against:
  note "[${label}] (b) runs:        MANUAL -- build the package and execute on a fixture"
  note "[${label}] (c) byte-parity: MANUAL -- diff stdout vs the original C binary on a trivial input"
}
eval_go "hello"  "${WORK}/hello.go"
eval_go "sickle" "${WORK}/sickle.go"

# --- 5. Scope reminder (cannot even be attempted) --------------------------
log "5. Out of scope for ANY rule-based transpiler"
note "c2go (and C2Rust-class tools) have NO C++ support. The bedtools/samtools/"
note "bcftools upstreams are C++ (reference_code/bedtools has 162 .cpp files),"
note "and prinseq is Perl (prinseq-lite.pl). No rule-based C->Go transpiler can"
note "even ingest these. They are excluded from the counterfactual by construction."

log "SUMMARY"
note "hello.c  transpile:  ${HELLO_TRANSPILE}"
note "sickle.c transpile:  ${SICKLE_TRANSPILE}"
note "See docs/manuscript/results/transpiler_baseline.md for the calibrated argument."
