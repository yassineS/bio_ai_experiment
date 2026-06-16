# bedtobam

A pure-Go, drop-in reimplementation of `bedtools bedtobam` (aka `bedToBam`).
It converts BED (or BED12) features into BAM alignments against a header built
from a genome (chrom-sizes) file.

## Usage

```
bedtobam [OPTIONS] -i <bed> -g <genome>
```

`-i` accepts a BED file (`-` = stdin, the default); `-g` is the required genome
file (two columns: `chrom<TAB>size`, one per line, order preserved). Output is a
BAM stream on stdout.

### Options

| Flag        | Meaning |
|-------------|---------|
| `-i FILE`   | BED/GFF/VCF input (`-` = stdin, the default). |
| `-g FILE`   | Genome (chrom-sizes) file. **Required.** |
| `-mapq INT` | MAPQ for the emitted records (0..255). Default `255`. |
| `-bed12`    | Treat the input as BED12; the CIGAR reflects the BED blocks. |
| `-ubam`     | Write uncompressed BAM (default writes compressed BAM). |
| `-h`, `--help` / `-v`, `--version` | Standard. |

BED input must be at least BED4 (a name field is required).

## Output

Each BED record becomes one BAM alignment: 0-based POS, the given MAPQ, FLAG `0`
(or `0x10` for a `-` strand), empty SEQ/QUAL, and a CIGAR of a single `<len>M`
(plain BED) or an N/M spliced CIGAR derived from the BED12 blocks. The BAM
header is `@HD VN:1.0 SO:unsorted`, `@PG ID:BEDTools_bedToBam VN:Vv2.31.1`, and
one `@SQ SN:<chrom> AS:<genome-file> LN:<size>` line per chromosome — matching
upstream byte-for-byte once decoded to SAM.

## Parity

`pkg/bedtobam/parity_test.go` builds the real upstream `bedtools` (and its
`htsutil` helper) from the vendored submodule, runs both implementations over
the same BED + genome, decodes both BAMs to SAM (BAM is BGZF-compressed, so the
decoded records and header are what must match), and diffs the SAM text
byte-for-byte across the `-mapq`, `-bed12`, and `-ubam` flags plus the upstream
`test/bedtobam` case. One upstream quirk is reproduced for parity and documented
in `docs/UPSTREAM_BUGS.md`: an out-of-genome chromosome is silently written
against the first `@SQ` reference.
