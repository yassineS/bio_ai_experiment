# bgzip — pure-Go Blocked GNU Zip Format

`bgzip` is a pure-Go reimplementation of htslib's `bgzip` command and the
underlying BGZF codec. BGZF is the foundational on-disk format used by:

- `.vcf.gz` (the random-access flavour of gzipped VCF) and tabix indices
- BAM and BCF
- any file that tabix expects to index

There is currently no actively-maintained Go BGZF library — `biogo/hts` has
been dormant since 2017. This package gives the repo a clean, well-tested,
stdlib-only BGZF building block so that downstream ports (`tabix`,
`samtools view`, `bcftools view`, and the `vcftools` random-access paths) can
share one implementation.

See `analysis/tool_ranking_2026.md` — `bgzip` is pick #1 of the 2026 next-up
list precisely because the rest of the SAM/BAM/VCF ecosystem rides on top of
it.

## What is BGZF?

A BGZF file is a sequence of independent gzip members ("blocks") concatenated
together. Each block:

1. Has an uncompressed payload of at most 64 KiB (htslib uses 65,280 bytes —
   a 256-byte safety margin so the worst-case deflate expansion still fits
   inside the 16-bit BSIZE field).
2. Carries the BGZF "BC" subfield in its gzip extra field: a 6-byte record
   `('B','C', SLEN=2, BSIZE)` where `BSIZE = total_compressed_block_size − 1`.
3. The terminal block of a well-formed BGZF stream is the canonical 28-byte
   gzip member representing an empty payload. Decoders use it to distinguish a
   complete stream from a silently truncated one.

Because every block is an independent gzip member, random access into the
stream only needs the compressed file offset of the start of a block plus the
in-block byte offset. That pair (compressed offset, uncompressed offset) is
exactly what `.gzi` index files record.

## Build

```bash
go build ./tools/bgzip/cmd/bgzip
```

The package has no third-party dependencies — it builds on `compress/flate`,
`encoding/binary`, and `hash/crc32` from the Go standard library.

## Usage

```text
bgzip [options] [file]
bgzip -d [options] file.gz
```

`-` as a filename means stdin (when used as input) or stdout (with `-c`).

### Flags

| Short | Long                  | Description |
|-------|-----------------------|-------------|
| `-c`  | `--stdout`            | Write output to stdout; keep input. |
| `-d`  | `--decompress`        | Decompress instead of compress. |
| `-f`  | `--force`             | Overwrite existing output file. |
| `-k`  | `--keep`              | Keep input file (do not delete it after success). |
| `-l N`| `--compress-level N`  | Compression level 0-9. Default 6. |
| `-t N`| `--threads N`         | Number of compression threads. Accepted; runs single-threaded in v1 (see "Deviations"). |
| `-b N`| `--offset N`          | Print the uncompressed byte offset that corresponds to compressed offset N. |
| `-s`  | `--size`              | Print the decompressed size of the file. |
| `-r`  | `--reindex`           | Write a `.gzi` index alongside `file.gz`. |
| `-h`  | `--help`              | Show help and exit. |
| `-v`  | `--version`           | Show version and exit. |

### Examples

Compress a VCF in place:

```bash
bgzip my.vcf            # produces my.vcf.gz, removes my.vcf
bgzip -k my.vcf         # keep my.vcf
bgzip -c my.vcf > my.vcf.gz   # stream to stdout
```

Decompress in place:

```bash
bgzip -d my.vcf.gz      # produces my.vcf, removes my.vcf.gz
bgzip -d -c my.vcf.gz   # write to stdout, keep my.vcf.gz
```

Query the index:

```bash
bgzip -s my.vcf.gz                # decompressed size (in bytes)
bgzip -b 65280 my.vcf.gz          # uncompressed offset at compressed offset
bgzip -r my.vcf.gz                # write my.vcf.gz.gzi
```

Pipe input via stdin (writes BGZF to stdout, never touches a file):

```bash
cat my.vcf | bgzip > my.vcf.gz
```

## Library

The CLI is a thin wrapper around `pkg/htsgo/bgzf`:

```go
// The package is named `bgzf` (htslib's term for the format); the
// convention across this repo is to import it under the `bgzip`
// alias so call sites read naturally.
import bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"

// Write
w := bgzip.NewWriter(out)
w.Write(payload)
w.Close()  // emits the BGZF EOF block

// Read
r, _ := bgzip.NewReader(in)
data, _ := io.ReadAll(r)
r.Close()

// Index / random-access helpers
offsets, _ := bgzip.Scan(in)        // every block's (offset, size)
size, _    := bgzip.DecompressedSize(in)
unc, _     := bgzip.UncompressedOffsetAt(in, compOff)
bgzip.WriteGZI(gziOut, offsets)
parsed, _ := bgzip.ReadGZI(gziIn)
```

## Deviations from upstream

This v1 implementation matches upstream `bgzip` for all documented behaviour
above with two known deviations:

1. **`-t` / `--threads` is single-threaded.** The flag is accepted (so wrapper
   scripts that pass `-t 4` keep working) but compression runs in a single
   goroutine. Block-level parallel compression is a natural follow-up — the
   block format is already independent — but is out of scope for this initial
   port. Decompression in upstream htslib is single-threaded too, so there is
   no behavioural difference there.
2. **`.gzi` semantics follow htslib's `bgzf_index_dump`:** the implicit
   leading `(0, 0)` block is *not* written, so the file is
   `8 + 16*(N−1)` bytes for an N-block input. This matches what tabix
   expects.

There are no other intentional differences. The Reader rejects gzip members
without the BC subfield (`ErrNoBCSubfield`), refuses streams missing the EOF
marker (`ErrTruncated`), and verifies per-block CRC32 and ISIZE.

## Layout

```text
tools/bgzip/
├── cmd/bgzip/main.go             # CLI entry point
├── pkg/bgzip/
│   ├── doc.go                    # package-level spec summary
│   ├── bgzip.go                  # Writer, Reader, block-header parser
│   ├── index.go                  # Scan, .gzi I/O, -b/-s helpers
│   ├── bgzip_test.go             # round-trip, EOF, header parsing tests
│   └── bgzip_extra_test.go       # error-path + edge-case tests
└── README.md
```

## What this unblocks

With a clean Go BGZF in hand, the natural follow-ups (NOT part of this PR) are:

- Wire `pkg/htsgo/bgzf` into `pkg/htsgo/iohelper` so any `.vcf.gz`
  with a BGZF extra field is decoded through this package — transparent
  support for the `vcftools` port and friends.
- Implement `tabix` on top of `BlockOffsets` and `ReadGZI`.
- Port `samtools view` / `bcftools view` whose entire I/O layer assumes BGZF.

## References

- SAM/BAM specification, section "The BGZF compression format".
- [htslib/bgzf.h](https://github.com/samtools/htslib/blob/develop/htslib/bgzf.h)
  and [htslib/bgzf.c](https://github.com/samtools/htslib/blob/develop/bgzf.c).
