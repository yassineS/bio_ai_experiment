# bedunionbedg - Combine multiple BedGraph files

A Go re-implementation of `bedtools unionbedg` (aka unionBedGraphs). Combines
two or more BEDGRAPH files into a single matrix: for every interval boundary
across all inputs it emits one line of the form

```
chrom  start  end  val1  val2  ...  valN
```

where `valI` is the value of input file `I` over `[start, end)`.

## Features

- Merges any number (>= 2) of BEDGRAPH inputs by a coordinate sweep that
  reproduces upstream's priority-queue algorithm, including its handling of
  inputs that are not globally chromosome-sorted.
- `-header` prints a `chrom start end` header, extended with `-names`.
- `-empty` (with `-g`) reports the gaps with no coverage in any file.
- `-filler TEXT` controls the value printed for inputs with no coverage over an
  interval (default `0`).
- Values are carried as raw strings, so integer, float, and arbitrary tokens
  round-trip unchanged.
- Built-in gzip support (`.gz` inputs/outputs).

## Build

```bash
go build ./tools/bedunionbedg/cmd/bedunionbedg
```

## Usage

```bash
bedunionbedg [options] -i FILE1 FILE2 .. FILEn
```

Each BedGraph file is assumed sorted by chrom/start with non-overlapping
intervals.

## Options

| Option | Description |
|--------|-------------|
| `-i FILE1 FILE2 ..` | Input BedGraph files (two or more). Required. |
| `-o, --output FILE` | Output file (`-` for stdout, default stdout). |
| `-header` | Print a header line (`chrom/start/end` + names of each file). |
| `-names NAME ..` | A list of names (one per file) printed in the header. |
| `-g FILE` | Genome (chrom-sizes) file, used to calculate empty regions. |
| `-empty` | Report empty regions (requires `-g`). |
| `-filler TEXT` | Text for intervals with no value (default `0`). |
| `-h, --help` | Show help and exit. |
| `-v, --version` | Show version and exit. |

## Examples

```bash
bedunionbedg -i 1.bg 2.bg 3.bg
bedunionbedg -header -i 1.bg 2.bg 3.bg -names WT-1 WT-2 KO-1
bedunionbedg -header -empty -g sizes.txt -i 1.bg 2.bg 3.bg
bedunionbedg -empty -g sizes.txt -filler N/A -i 1.bg 2.bg 3.bg
```

## Parity notes

- The header for an unnamed run is exactly `chrom start end` (the `1 2 3`
  default shown in upstream's `-examples` help text does not actually appear in
  the binary's output; this port matches the binary).
- Each input line must have exactly four tab-separated fields; `track`,
  `browser`, and `#` lines are skipped, matching upstream's BedGraph parser.
- Validated byte-for-byte against the upstream `bedtools unionbedg` binary for:
  basic union, `-header`, `-names`, `-empty -g`, `-filler`, multi-chromosome
  unsorted inputs, and floating-point depth values.

## Testing

```bash
go test ./tools/bedunionbedg/...
go test -cover ./tools/bedunionbedg/pkg/bedunionbedg
```
