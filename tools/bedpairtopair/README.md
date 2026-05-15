# bedpairtopair

Pure-Go reimplementation of `bedtools pairtopair`.

Reports overlaps between **two BEDPE files**. For each A pair, the tool finds
B pairs where the chosen `-type` predicate over A's two ends is satisfied.

## Usage

```text
bedpairtopair -a <BEDPE> -b <BEDPE> [OPTIONS]

I/O:
  -a FILE                BEDPE A input ('-' = stdin). Transparent gzip.
  -b FILE                BEDPE B input.
  -o, --output FILE      Output (default stdout).

Reporting:
  -type {both|notboth|either|neither}   (default both)

Filters:
  -f FRAC                Min fraction of A end covered by B (default 1e-9).
  -slop N                Add N bp of slop to each end of A before search.
  -ss                    Make -slop strand-aware (extend in strand dir only).
  -is, --ignore-strand   Ignore strand entirely.
  -s, --same-strand      Require matching strands.
  -S, --opposite-strand  Require opposite strands.
  -rdn                   Require A and B pairs to have different names.

Standard:
  -h, --help             Show help.
  -v, --version          Show version.
```

### `-type` semantics

| value     | meaning                                                          |
| --------- | ---------------------------------------------------------------- |
| `both`    | A pair matches when one B pair has both its ends overlapped by A's two ends (in either orientation) — default |
| `notboth` | A pair is emitted bare when no B satisfies the `both` condition  |
| `either`  | A pair matches any B whose end overlaps either end of A          |
| `neither` | A pair is emitted bare when no end of A overlaps any end of any B |

`both` and `either` emit lines of the shape `<A pair>\t<B pair>`. `notboth`
and `neither` emit `<A pair>` alone.

## Algorithm

For each A pair we test the four end-vs-end combinations: A1×B1, A1×B2,
A2×B1, A2×B2. A B pair satisfies `both` when the same B row appears in
`{A1×B1 ∩ A2×B2} ∪ {A1×B2 ∩ A2×B1}` (i.e. each A end picks a different B end
on the same B record).

Each end-vs-end test requires the fraction of A-end length covered by the B
end to be ≥ `-f` (default `1e-9` ⇒ effectively 1 bp).

Strand defaults: matching strands are required on every end-vs-end pair
(matches upstream); `-is` disables this; `-s`/`-S` are equivalent flag-spelling
for the same enforcement / opposite enforcement.

`-slop N` extends each A end by N bp before testing overlaps; `-ss` makes
this directional based on the strand of the corresponding end.

## Parity matrix

| feature                  | upstream `pairtopair` | bedpairtopair |
| ------------------------ | --------------------- | ------------- |
| `-type both/notboth/either/neither` | yes      | yes |
| `-f` fractional overlap  | yes                   | yes |
| `-is` ignore-strand      | yes                   | yes |
| `-slop` / `-ss`          | yes                   | yes |
| `-rdn` require diff names | yes                  | yes |
| `-s`/`-S` same/opposite strand | not exposed (defaults to same-strand) | yes (extension, consistent with rest of toolkit) |

Deviations:

* The upstream tool does not expose `-s`/`-S` because same-strand is the
  default — we add the flags to match the rest of our bedtools-port flag
  vocabulary.
* Score field treated as opaque string (matches upstream which permits
  non-numeric).

## Tests

* `bedpairtopair_test.go`: 18 unit tests (basic, `-type` mode coverage,
  slop, strand options, `-rdn`, error paths).
* `parity_test.go`: 4 hand-curated parity tests + 1 with-slop test + 1
  documented skip for stranded slop.
* Coverage: ≥90% on `pkg/bedpairtopair` (see PR body for exact number).

## Layout

```text
tools/bedpairtopair/
├── cmd/bedpairtopair/main.go
├── pkg/bedpairtopair/
│   ├── bedpairtopair.go
│   ├── bedpairtopair_test.go
│   └── parity_test.go
├── testdata/parity/
│   ├── a.bedpe
│   └── b.bedpe
└── README.md
```
