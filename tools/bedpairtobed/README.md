# bedpairtobed

Pure-Go reimplementation of `bedtools pairtobed`.

Reports overlaps between a **BEDPE** A file (paired-end intervals) and a
regular **BED** B file. For each pair in A, each of its two ends is queried
against B; the `-type` flag controls which pairs are emitted and in what
shape.

## Usage

```text
bedpairtobed -a <BEDPE> -b <BED> [OPTIONS]

I/O:
  -a FILE                BEDPE A input ('-' = stdin). Transparent gzip.
  -b FILE                BED B input.
  -o, --output FILE      Output (default stdout).

Reporting:
  -type {either|both|notboth|neither|xor|notxor}   (default either)

Filters:
  -f FRAC                Min fraction of A end covered by B (default 1e-9).
  -s, --same-strand      Require matching strands between A end and B hit.
  -S, --opposite-strand  Require opposite strands.
  -is, --ignore-strand   Ignore strand entirely (overrides -s/-S).

Standard:
  -h, --help             Show help.
  -v, --version          Show version.
```

### `-type` semantics

| value      | meaning                                                          |
| ---------- | ---------------------------------------------------------------- |
| `either`   | emit pair when at least one end overlaps any B record (default)  |
| `both`     | emit pair when both ends overlap some B record                   |
| `notboth`  | emit pair when not both ends overlap a B record                  |
| `neither`  | emit pair when no end overlaps any B record                      |
| `xor`      | emit pair when exactly one end overlaps                          |
| `notxor`   | emit pair when both ends or neither end overlaps (deviation: not in upstream) |

## Output format

Per upstream:

* For `either`, `both`, `xor`, the partial-overlap branch of `notboth`, and
  the both-hit branch of `notxor`: each emitted pair is followed by a tab
  and the full BED record of the hit.
* For `neither` and the no-hit branch of `notboth`/`notxor`: only the BEDPE
  line is emitted.

## Parity matrix

| feature                  | upstream `pairtobed`    | bedpairtobed |
| ------------------------ | ----------------------- | ------------ |
| `-type either/both/...`  | yes                     | yes (5 upstream + `notxor` extension) |
| `-f` fractional overlap  | yes                     | yes          |
| `-s` same-strand         | yes                     | yes          |
| `-S` opposite-strand     | yes                     | yes          |
| `-is` ignore-strand      | yes                     | yes          |
| BAM input                | yes                     | not implemented (out of scope: this PR is BEDPE-only) |
| BAM output `-bedpe/-ed`  | yes                     | not implemented |

Deviations:

* `notxor` is our addition — upstream only ships the 5 modes above.
* Score column treatment matches upstream: the score is preserved as a
  string and emitted verbatim.
* BAM input / output is not covered here — that path lives in `samtools`
  in this tree, and the upstream BAM-PE path was never the primary use case.

## Tests

* `bedpairtobed_test.go`: 14 unit tests (table-driven and end-to-end).
* `parity_test.go`: 5 hand-curated parity tests + 1 documented skip.
* Coverage: ≥90% on `pkg/bedpairtobed` (see PR body for exact number).

## Layout

```text
tools/bedpairtobed/
├── cmd/bedpairtobed/main.go
├── pkg/bedpairtobed/
│   ├── bedpairtobed.go
│   ├── bedpairtobed_test.go
│   └── parity_test.go
├── testdata/parity/
│   ├── a.bedpe
│   └── b.bed
└── README.md
```
