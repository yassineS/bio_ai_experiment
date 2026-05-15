# bednuc - Nucleotide content of BED intervals

A Go re-implementation of `bedtools nuc`. For each interval in a BED file,
fetches the sequence from an indexed FASTA reference and emits a
per-interval nucleotide composition profile.

## Features

- Per-interval `%AT`, `%GC`, A/C/G/T/N/other counts, and sequence length.
- Optional sequence emission (`-seq`).
- Optional substring counting per interval (`-pattern STR`) with optional
  case-insensitive mode (`-C`).
- Strand-aware mode (`-s`) reverse-complements `-`-strand intervals
  before counting.
- Best-effort `-fullHeader` support: contigs whose FASTA header carries
  whitespace can still be matched when the BED supplies the full string.
- Random-access FASTA via `pkg/bioformats/fasta` (uses the sibling
  `.fai` when present, or builds the index on the fly).
- Pure Go, no third-party dependencies.
- Transparent gzip/BGZF input on the BED side and `-` for stdin.

## Build

```bash
go build ./tools/bednuc/cmd/bednuc
```

## Usage

```bash
bednuc -fi <fasta> -bed <bed> [options]
```

### Options

- `-fi FASTA`        Indexed FASTA reference (required)
- `-bed BED`         BED/GFF/VCF intervals (`-` for stdin)
- `-s, --strand`     Reverse-complement `-` strand intervals before counting
- `-seq, --seq`      Emit the extracted sequence as an extra column
- `-pattern STR`     Count occurrences of substring STR per interval
- `-C, --ignorecase` Ignore case when matching `-pattern` (upstream default
  is case-sensitive)
- `-fullHeader`      Match contigs by their full FASTA header (including
  whitespace) via a best-effort fallback map
- `-o, --output FILE` Output file (default: stdout)
- `-h, --help`       Show help
- `-v, --version`    Show version

### Output

The first matching record emits a `#`-prefixed column-header line. Then,
for each record, the original BED columns are appended with:

```text
%AT  %GC  #A  #C  #G  #T  #N  #oth  seq_len  [seq]  [pattern_count]
```

The bracketed columns appear only when `-seq` / `-pattern` are set, in
the order `seq` then `pattern_count` (matching upstream).

Floating-point columns use printf `%f` (6 decimal digits).

## Deviations from upstream

- The FASTA index always keys on the first whitespace-delimited token of
  the header. `-fullHeader` resolves contigs by walking the FASTA once
  and building a fallback map from `full-header → first-token`; this
  covers the common case where the BED column is `"chr1 extra info"`
  but the index key is `chr1`. Cases where the contig name itself
  contains spaces *and* differs between the index and the lookup are
  not exercised here.
- Zero-length features and features beyond the contig end are skipped
  with a warning on stderr, matching upstream.

## Tests

```bash
go test ./tools/bednuc/...
```

The parity fixtures under `testdata/parity/` are hand-computed: upstream
ships no `nuc/` test directory, so the expected outputs were derived
against the upstream output format and counting rules in
`reference_code/bedtools/src/utils/sequenceUtilities/sequenceUtils.cpp`
and `src/nucBed/nucBed.cpp`.

## Performance

Linear scan over the BED with constant-time random access into the
FASTA via the geometry-based byte offsets in
`pkg/bioformats/fasta.RandomAccess`. Memory footprint is the FASTA
index (≈48 B/contig) plus a single per-record buffer.
