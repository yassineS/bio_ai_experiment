# bedoverlap - Compute overlap/distance between two intervals on a line

A Go re-implementation of `bedtools overlap` (aka getOverlap). For each input
line it computes the overlap (positive) or distance (negative) between the two
intervals named by four 1-based columns, and appends the result as a new
trailing column on the same line.

The result is `min(end1, end2) - max(start1, start2)`: a positive value is the
number of overlapping bases; a non-positive value is the (negative) gap between
the intervals.

## Features

- Works on any tab-delimited line; the four interval columns are chosen with
  `-cols start1,end1,start2,end2`.
- Emits the original line verbatim, then a tab and the computed value, exactly
  like upstream.
- Skips lines that tokenize to one field or fewer.
- Reports a non-numeric column with the same message and exit status as
  upstream.
- Built-in gzip support (`.gz` inputs/outputs).

## Build

```bash
go build ./tools/bedoverlap/cmd/bedoverlap
```

## Usage

```bash
bedoverlap -i <file> -cols s1,e1,s2,e2
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input file (`-` or `stdin` for stdin, default stdin). |
| `-o, --output FILE` | Output file (`-` for stdout, default stdout). |
| `--cols COLS` | Comma-separated 1-based columns: `start1,end1,start2,end2`. Required. |
| `-h, --help` | Show help and exit. |
| `-v, --version` | Show version and exit. |

## Example

```bash
$ bedtools window -a A.bed -b B.bed -w 10 | bedoverlap -i stdin -cols 2,3,6,7
chr1	10	20	A	chr1	15	25	B	5
chr1	10	20	C	chr1	25	35	D	-5
```

## Parity notes

- The four columns are extracted by splitting the line on tabs (matching
  upstream's tab-delimited tokenizer); the original line text is preserved
  byte-for-byte in the output.
- A column is treated as numeric when it has a leading integer prefix
  (upstream's `strtol` semantics); otherwise the run aborts with
  `One of your columns appears to be non-numeric at line N. Exiting...`.
- Validated byte-for-byte against the upstream `bedtools overlap` binary for the
  `--help` example (`-cols 2,3,6,7`) plus touching, nested, and disjoint
  (negative-distance) interval pairs.

## Testing

```bash
go test ./tools/bedoverlap/...
go test -cover ./tools/bedoverlap/pkg/bedoverlap
```
