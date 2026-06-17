# bedbamtobed

A pure-Go, drop-in reimplementation of `bedtools bamtobed` (aka `bamToBed`).
It converts BAM/SAM alignments into BED6, blocked BED12, or BEDPE records.

## Usage

```
bedbamtobed [OPTIONS] -i <bam>
```

`-i` accepts a BAM file, a raw (decompressed) BAM body, or SAM text; `-` (the
default) reads from stdin. Output goes to stdout.

### Options

| Flag        | Meaning |
|-------------|---------|
| `-i FILE`   | BAM/SAM input (`-` = stdin, the default). |
| `-bedpe`    | Write BEDPE format. Requires the BAM to be grouped/sorted by query name. |
| `-mate1`    | With `-bedpe`, always report mate one as the first BEDPE block. |
| `-bed12`    | Write blocked BED12. The CIGAR's N/D-separated runs become BED blocks. |
| `-split`    | Report each split alignment block as a separate BED entry (splits on N). |
| `-splitD`   | Split on N **and** D CIGAR ops. Implies `-split`. |
| `-ed`       | Use BAM edit distance (NM tag) as the BED score. |
| `-tag TAG`  | Use another numeric BAM tag for the BED score. Disallowed with `-bedpe`. |
| `-color R,G,B` | itemRgb for BED12 (default `255,0,0`). |
| `-cigar`    | Append the CIGAR string as a trailing column (BED6 only). |
| `-h`, `--help` / `-v`, `--version` | Standard. |

## Output

- **Default (BED6):** `chrom start end name score strand`, where `name` is
  QNAME with a `/1` or `/2` mate suffix and `score` is MAPQ (or the `-tag`/`-ed`
  value).
- **`-bed12`:** the BED6 columns plus thickStart/thickEnd, itemRgb, blockCount,
  blockSizes and blockStarts derived from the CIGAR.
- **`-bedpe`:** `chrom1 start1 end1 chrom2 start2 end2 name score strand1
  strand2`. The two mates are ordered by `(chrom, start)` unless `-mate1` forces
  mate one first; the score is the minimum MAPQ of the pair, or (with `-ed`) the
  summed edit distance.

## Parity

`pkg/bedbamtobed/parity_test.go` builds the real upstream `bedtools` (and its
`htsutil` helper) from the vendored submodule and compares output byte-for-byte
across the upstream `test/bamtobed` fixtures and the documented flag set
(`-split`, `-splitD`, `-bed12`, `-cigar`, `-tag`, `-color`, `-bedpe`, `-mate1`,
`-ed`). Two upstream quirks are reproduced for parity and documented in
`docs/UPSTREAM_BUGS.md`: the spurious extra column emitted by `-tag` combined
with `-split`.

### Deliberate divergences (faithful upstream matches, not gaps)

- **`-cigar` together with `-splits`** is rejected with upstream's exact
  message ("Cannot use -cigar with -splits.  Not yet supported.").
  Upstream `bamToBed.cpp` itself errors on this combination — there is no
  per-block CIGAR to attach once a spliced alignment is split into blocks —
  so reproducing the error is byte-parity. Implementing it would emit output
  where upstream errors, i.e. a divergence, not a fix. (The guard mirrors
  upstream's `useEditDistance && useCigar` condition, message and all.)
