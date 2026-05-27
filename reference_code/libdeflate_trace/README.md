# libdeflate_trace

Slice-0 oracle harness for the bio_ai_experiment libdeflate-parity port
scoped in `docs/HTSGO_ROADMAP.md` (appendix "libdeflate-parity port for
BGZF byte equality").

This directory contains a small C harness that wraps upstream libdeflate
with `#ifdef BTRACE` instrumentation at every internal decision point of
the deflate encoder.  Running the harness against a fixed corpus produces
trace text files that subsequent slices port the pure-Go encoder against.

## Layout

```
src/btrace.h                  BTRACE_*() macros that emit one event per
                              decision when -DBTRACE is set.
patches/apply_btrace.py       Copies upstream libdeflate sources into
                              build/src/ and patches deflate_compress.c
                              with BTRACE invocations via anchored
                              string replacement (fails loudly when an
                              anchor drifts upstream).
patches/make_corpus.py        Deterministic test-corpus generator
                              (empty / single byte / 100x'A' / 64K LCG /
                              SAM-shaped 65280 B BGZF payload).
cmd/dtrace/dtrace.c           CLI: compress one file at a given level
                              and dump the BTRACE events.
cmd/dverify/dverify.c         Round-trip sanity check (decode the .gz
                              and assert it matches the original).
Makefile                      build / oracle / test / clean.
```

## Build & test

The upstream source lives in the `reference_code/libdeflate` submodule.
Initialize it once (`git submodule update --init reference_code/libdeflate`)
then:

```
cd reference_code/libdeflate_trace
make build      # gcc / clang, no autotools, no cmake
make oracle     # regenerates corpus + traces under
                #   ../../pkg/htsgo/libdeflate/testdata/oracle/
make test       # build + oracle + round-trip verification
```

`make` honours `CC=` (defaults to `gcc`) and `PYTHON=` (defaults to
`python3`).

## Trace format

One event per line, whitespace-separated `KEY=VALUE` pairs.  Examples:

```
BTRACE HEADER version=1 source=<path> level=6 format=gzip
BTRACE BLOCK_BEGIN pos=0
BTRACE MIN_MATCH_LEN stage=initial min_len=3
BTRACE LIT pos=0 byte=0x41
BTRACE MATCH pos=2 len=98 offset=1
BTRACE HUFFMAN_LEN table=litlen sym=65 codelen=4
BTRACE SPLIT_CHECK block_length=5954 total_delta=71084 cutoff=901842 decision=0
BTRACE BLOCK_COSTS dynamic=124 static=42 uncompressed=840
BTRACE BLOCK_HEADER type=1 size=100 last=1
BTRACE FOOTER in_bytes=100 out_bytes=24
```

`type` in `BLOCK_HEADER` is `0=stored`, `1=static`, `2=dynamic`.

`BTRACE_BITS` events (one per bitstream write) are gated behind an extra
`-DBTRACE_BITS` so day-to-day oracle generation doesn't produce gigabyte
traces.  Slices 1 and 3 enable that flag for small fixtures only.

## Constraints

This is C code that links against the upstream libdeflate sources.  It
is a test-time tool only — no cgo bindings, no impact on the main Go
module's dependency policy.  The Go side consumes only the byte
artifacts written to `pkg/htsgo/libdeflate/testdata/oracle/`.

## Status

Slice 0 only.  Subsequent slices port the Go encoder against the traces
produced here, starting with constants + static-block emission (slice 1).
See the appendix in `docs/HTSGO_ROADMAP.md` for the full plan.
