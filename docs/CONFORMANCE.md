# Conformance and silent-corruption edge-case batteries

This document describes two test batteries that validate our re-implemented
binaries against the *originals' own* test corpora and against a ranked list of
silent-corruption failure classes:

- `pipeline/conformance/` — runs htslib's and htscodecs' own shipped test
  fixtures through our `samtools`/`bcftools` and our in-tree CRAM codec.
- `pipeline/edgecases/` — discrete, named Go tests, one per class of failure
  where a tool can produce *wrong* output without crashing or warning.

Both batteries SKIP gracefully (they do not fail) when a prerequisite — a
submodule corpus, a reference, or a vendored binary — is absent. That makes them
safe to run on a partially-provisioned machine while still being meaningful
where the data is present.

## Why the originals' own corpora

Passing the upstream projects' *own* adversarial fixtures is far more persuasive
than passing fixtures we authored ourselves: htslib's `test/` tree was built by
the htslib maintainers specifically to break SAM/BAM/CRAM parsers, and the
htscodecs corpus is the byte-exact compliance oracle the C reference
implementation is itself tested against.

## Prerequisites

### Submodule corpora

The conformance battery reads fixtures from two submodules. Initialise them
once:

```bash
git submodule update --init reference_code/htslib reference_code/htscodecs
```

- `reference_code/htslib/test/` — SAM/BAM/CRAM/VCF fixtures: empty files, CRLF
  (`index_dos.sam`), padded / bad-CIGAR (`c1#pad*.sam`, `c2#pad.sam`), no-SEQ
  (`c1#noseq.sam`), unknown-reference (`c1#unknown.sam`), reference-bounds
  (`c1#bounds.sam`), `bgzf_boundaries/` (records straddling BGZF blocks), and
  `longrefs/` (positions beyond 2^31).
- `reference_code/htscodecs/tests/dat/` — the rANS (`r4x8`, `r4x16`),
  adaptive-arithmetic (`arith`) and fqzcomp compliance vectors.

### Upstream binaries

The batteries compare our output against the upstream tools. Build them once:

```bash
git submodule update --init reference_code/samtools reference_code/bcftools reference_code/htslib
(cd reference_code/htslib   && autoreconf -i && ./configure && make)
(cd reference_code/samtools && autoreconf -i && ./configure && make)
(cd reference_code/bcftools && autoreconf -i && ./configure && make)
```

The resolver `pipeline/internal/upstream` finds these at
`reference_code/<tool>/<tool>` (and `reference_code/htslib/{bgzip,tabix}`). In an
isolated git worktree where the submodules are unpopulated, symlink the binaries
from your main checkout instead — see `pipeline/README.md`.

Our own binaries are built on demand by `upstream.OurBinary` (one `go build` per
tool into a temp cache), so no manual build of our tools is required.

## Running

### Via the driver command

```bash
# both batteries, verbose, human-readable:
go run ./pipeline/cmd/conformance

# filter by test-name regexp:
go run ./pipeline/cmd/conformance -run CRAM

# only one battery:
go run ./pipeline/cmd/conformance -pkgs conformance
go run ./pipeline/cmd/conformance -pkgs edgecases

# machine-readable:
go run ./pipeline/cmd/conformance -json
```

### Via `go test` directly

```bash
go test ./pipeline/conformance/... ./pipeline/edgecases/...
go test -v ./pipeline/edgecases/... -run CRAMReferenceHandling
```

## Reading the results

Each subtest reports one of:

- **PASS** — our binary/codec reproduced the upstream behaviour (round-trip
  parity, byte-identity, or decision-identity, as appropriate).
- **SKIP** — a prerequisite is absent (submodule not initialised, upstream
  binary missing, fixture not generated). The skip message names the missing
  item and the command that provides it.
- **FAIL** — a genuine divergence from upstream. Failures in these batteries are
  *findings*, not flakes: messages prefixed `PARITY GAP:` describe the exact
  divergence (e.g. an index that is functionally usable but not byte-identical,
  or a `norm` re-indexing mismatch).

## What each test covers

### Conformance (`pipeline/conformance/`)

| Test | Corpus | Assertion |
| --- | --- | --- |
| `TestHtslibSAM_RoundTripParity` | htslib `test/*.sam` | our SAM→BAM→SAM record set == upstream's |
| `TestHtslibSAM_CRAMRoundTrip` | htslib `test/*.sam` + `*.fa` | SEQ survives CRAM; our CRAM round-trip == upstream's |
| `TestHtslibBGZFBoundaries` | `test/bgzf_boundaries/*.bam` | records across BGZF block boundaries == upstream |
| `TestHtslibLongRefs` | `test/longrefs/longref.sam` | positions > 2^31 accepted like upstream |
| `TestHtslibEmptyFile` | `test/emptyfile` | empty/headerless input handled like upstream |
| `TestHtscodecs_RANS4x16` | `tests/dat/r4x16` | decode == raw and re-encode == reference vector (byte-exact) |
| `TestHtscodecs_RANS4x8` | `tests/dat/r4x8` | decode == raw and re-encode == reference vector (byte-exact) |
| `TestHtscodecs_Arith` | `tests/dat/arith` | decode == raw input (byte-exact) |

### Edge cases (`pipeline/edgecases/`)

Named after the manuscript's ranked silent-corruption list:

1. `TestCRAMReferenceHandling` — reference-relative CRAM base corruption
   (external `-T` reference and the multi-reference-slice / `RefSeqID=-2` case).
   Highest value: silent base corruption is the worst failure.
2. `TestBCFToolsNormReindex` — `norm -m-`/`-m+` re-indexing of `Number=A/R/G`
   INFO/FORMAT vectors (incl. AD, PL).
3. `TestIndexByteIdentity` — `.bai`/`.csi`/`.tbi` byte-identity vs upstream.
4. `TestSortStabilityStrnumCmp` — coordinate and `-n` queryname sort at `-@1`;
   tie-break order, `strnum_cmp` natural ordering, and unmapped/`*` placement.
   Multi-thread tie order is allowed to differ, so the test pins `-@1`.
5. `TestCalmdMDNMTags` — `samtools calmd` MD/NM recomputation across `=`/`X`
   CIGAR ops, `N` skips, and reference ambiguity codes.
6. `TestQualPLULPNonImpact` — any `bcftools call` QUAL last-ULP difference vs
   upstream must never flip GT or FILTER (lightweight version; the full
   statistical version lives in the GIAB harness).

## Running on an external machine

The batteries are self-contained Go tests with no third-party dependencies. On a
fresh checkout:

```bash
# 1. fetch the corpora and upstream sources
git submodule update --init \
  reference_code/htslib reference_code/htscodecs \
  reference_code/samtools reference_code/bcftools

# 2. build the upstream binaries (see "Upstream binaries" above)

# 3. run
go run ./pipeline/cmd/conformance
```

Tests that lack their corpus or an upstream binary SKIP with a message naming
exactly what to install, so a first run on a bare machine is informative rather
than a wall of failures.
