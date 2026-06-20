# GIAB variant-calling concordance harness

Manuscript experiment **C2 / P1**. This is the run recipe and data-acquisition
guide for `giab-concordance`, the turnkey harness that proves our `bcftools`
produces a variant call set that is **record-exact with upstream** (modulo a
QUAL/PL last-place wobble that never flips a genotype or a PASS/FAIL verdict)
and **biologically concordant with the GIAB truth set** (equal precision /
recall / F1 to upstream, per GA4GH stratification).

- Library: [`pipeline/giab/`](../pipeline/giab)
- CLI: [`pipeline/cmd/giab-concordance`](../pipeline/cmd/giab-concordance)

The harness degrades gracefully: every external prerequisite (GIAB VCF/BED,
reads BAM, `hap.py`/`vcfeval`, our and upstream `bcftools`) is checked, and a
missing one produces a clear **SKIP** with a pointer back to this document and
an exit code of `0`. It therefore runs unchanged in CI with no GIAB data; the
heavy run happens on an external machine that has the data.

## What it does

Given a JSON config, the harness runs three stages:

1. **Produce two call sets.** With OUR `bcftools` and the UPSTREAM `bcftools`,
   each:

   ```sh
   bcftools mpileup -f REF BAM | bcftools call -mv -O z -o ours.vcf.gz
   bcftools index -t ours.vcf.gz
   ```

   producing `ours.vcf.gz` and `up.vcf.gz`.

2. **Ours-vs-upstream concordance** (the core record-exact claim). The two call
   sets are compared record-by-record within the high-confidence BED. The
   report gives identical / differ counts and, crucially, a **ULP-flip
   detector**: where QUAL or PL differ only at the last unit-in-the-last-place
   (the Phred ULP — QUAL to the last printed decimal, PL by one integer point),
   it verifies the genotype (GT) and the FILTER (PASS/FAIL) are unchanged. The
   headline result is *"N sites differ in QUAL/PL, 0 of which flip a genotype or
   PASS/FAIL"* — the result that makes the libm floor a non-issue. Any genotype
   or PASS/FAIL flip makes the stage FAIL (exit 1).

3. **Biological concordance vs GIAB truth.** If `hap.py` or `vcfeval` is
   available, both call sets are scored against the GIAB truth VCF restricted to
   the high-confidence BED, stratified by the provided GA4GH stratification BEDs
   (whole-region plus one run per stratum). Precision / recall / F1 for SNVs and
   indels are reported **ours vs upstream side by side** — the Sentieon-style
   two-pronged template (record-exact *and* biologically equivalent).

It writes `giab_concordance.json` (machine-readable) and `giab_concordance.md`
(human-readable) into the output directory.

## Data acquisition

All paths below are for the GIAB v4.2.1 benchmark. Substitute the sample
(`HG001`..`HG007`) and reference build (`GRCh37`, `GRCh38`; T2T where
published).

### 1. GIAB truth VCF + high-confidence BED (NIST GIAB FTP)

The authoritative release tree is the NIST GIAB FTP:

```text
ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/
```

(HTTPS mirror: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/`)

Per-sample v4.2.1 benchmark directory layout:

```text
release/<SAMPLE>_<NA-ALIAS>/NISTv4.2.1/<BUILD>/
    <SAMPLE>_<BUILD>_1_22_v4.2.1_benchmark.vcf.gz
    <SAMPLE>_<BUILD>_1_22_v4.2.1_benchmark.vcf.gz.tbi
    <SAMPLE>_<BUILD>_1_22_v4.2.1_benchmark_noinconsistent.bed
```

Sample / alias map (v4.2.1):

| Sample | NA / alias | FTP subdirectory                       |
|--------|------------|-----------------------------------------|
| HG001  | NA12878    | `release/NA12878_HG001/NISTv4.2.1/`     |
| HG002  | NA24385    | `release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/` |
| HG003  | NA24149    | `release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/` |
| HG004  | NA24143    | `release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/` |
| HG005  | NA24631    | `release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/` |
| HG006  | NA24694    | `release/ChineseTrio/HG006_NA24694_father/NISTv4.2.1/` |
| HG007  | NA24695    | `release/ChineseTrio/HG007_NA24695_mother/NISTv4.2.1/` |

Example (HG002, GRCh38):

```sh
BASE=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
DIR=$BASE/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
curl -O $DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -O $DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
curl -O $DIR/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

For GRCh37 swap `GRCh38` → `GRCh37` in the path and filenames. The CMRG v1.00
benchmark (challenging medically-relevant genes) lives under each sample's
`*_CMRGv1.00/` directory on the same FTP.

### 2. GA4GH / GIAB genome stratifications (incl. CMRG & difficult regions)

The stratification BEDs come from the GIAB genome-stratifications project:

```text
https://github.com/genome-in-a-bottle/genome-stratifications
```

with the prebuilt BEDs mirrored on the same NIST FTP under
`release/genome-stratifications/` (e.g.
`.../genome-stratifications/v3.1/GRCh38/`). Relevant strata:

- **CMRG** — challenging medically-relevant genes:
  `GRCh38_CMRG_*.bed.gz` (and the matching CMRG benchmark VCF above).
- **All difficult regions** — the union difficulty mask:
  `union/GRCh38_alldifficultregions.bed.gz`.
- Lower-level difficulty: `LowComplexity/`, `SegmentalDuplications/`,
  `mappability/` (e.g. `GRCh38_lowmappabilityall.bed.gz`),
  `GCcontent/`, `FunctionalRegions/`.

`giab-concordance` consumes any number of these as `stratifications` entries
(name + BED path); `.gz` (gzip/bgzip) BEDs are read transparently.

### 3. Reference builds

- **GRCh37**: `human_g1k_v37` / `hs37d5` (1000 Genomes), or the GIAB-provided
  GRCh37 reference. Index with `samtools faidx REF.fa`.
- **GRCh38**: `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`
  (the GIAB / GA4GH standard "no-alt" analysis set).
- **T2T (CHM13v2.0)**: `chm13v2.0.fa` from the T2T consortium, where v4.2.1
  stratifications / benchmarks are published.

The reference FASTA must have a `.fai` (`samtools faidx`) and, for `vcfeval`,
an RTG **SDF** (`rtg format -o REF.sdf REF.fa`).

### 4. Reads (aligned BAM)

The call sets are produced from an aligned-reads BAM. Options:

- **GIAB-provided alignments.** GIAB publishes per-sample read alignments
  (Illumina, PacBio HiFi, ONT) under
  `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/<SAMPLE>/`
  (e.g. `.../data/AshkenazimTrio/HG002_NA24385_son/`). Pick a BAM aligned to
  the same build as the truth set.
- **Align from FASTQ.** Download GIAB FASTQs from the same `data/` tree and
  align with your aligner of choice (e.g. `bwa-mem2 mem REF reads | samtools
  sort -o reads.bam && samtools index reads.bam`).

A whole-genome BAM is large; for a fast smoke run restrict to a single
chromosome (`samtools view -b in.bam chr20 > chr20.bam`) and intersect the
high-confidence BED to that chromosome.

## hap.py and vcfeval

The harness shells out to whichever benchmarking engine is provided or on
`PATH`:

- **hap.py** (Illumina) is typically run via its container
  `pkrusche/hap.py`. To make it visible to the harness either install it on
  `PATH` as `hap.py`, or set `happy_bin` to a wrapper script that invokes the
  container, e.g.:

  ```sh
  cat > /usr/local/bin/hap.py <<'EOF'
  #!/bin/sh
  exec docker run --rm -v "$PWD:$PWD" -v /data:/data -w "$PWD" \
      pkrusche/hap.py /opt/hap.py/bin/hap.py "$@"
  EOF
  chmod +x /usr/local/bin/hap.py
  ```

  The harness invokes:
  `hap.py TRUTH QUERY -r REF -f HIGHCONF [-T STRATBED] -o PREFIX`
  and parses `PREFIX.summary.csv`.

- **vcfeval** ships with **RTG Tools** (`rtg vcfeval`). Provide it via
  `vcfeval_bin` (point at `rtg` or a standalone `vcfeval`) and set
  `sdf_template` to the RTG SDF reference (`rtg format -o REF.sdf REF.fa`). The
  harness invokes:
  `rtg vcfeval -b TRUTH -c QUERY -t SDF -e HIGHCONF [--bed-regions STRAT] -o OUT`
  and parses `OUT/summary.txt`.

If neither engine is found, the biological stage SKIPs; the ours-vs-upstream
record concordance (which needs no engine) still runs.

## Running it

1. Emit a template config and fill it in:

   ```sh
   go run ./pipeline/cmd/giab-concordance -print-config > run.json
   $EDITOR run.json
   ```

   Config fields (all paths absolute or relative to the working dir):

   | Field | Meaning |
   |-------|---------|
   | `sample` | GIAB sample label, e.g. `HG002` (informational). |
   | `build` | Reference build label, e.g. `GRCh38` (informational). |
   | `reference` | Indexed reference FASTA (`.fa` + `.fai`). |
   | `reads_bam` | Aligned-reads BAM the call sets are produced from. |
   | `truth_vcf` | GIAB benchmark truth VCF (v4.2.1). |
   | `high_conf_bed` | GIAB high-confidence region BED. |
   | `stratifications` | List of `{name, bed}` GA4GH stratification BEDs (CMRG, difficult regions, …). |
   | `our_bcftools` | Our `bcftools` binary. Empty → built from this repo. |
   | `upstream_bcftools` | Upstream `bcftools`. Empty → resolved via `pipeline/internal/upstream`. |
   | `happy_bin` / `vcfeval_bin` | Benchmarking engine. Empty → probe `PATH` (`hap.py`, then `rtg`/`vcfeval`). |
   | `sdf_template` | RTG SDF reference (vcfeval only). |
   | `qual_ulp` | QUAL ULP tolerance in Phred units (default `0.5`). |
   | `out_dir` | Output directory for the reports. |
   | `max_diffs` | Cap on differing sites embedded in the report (default 200; `-1` = unlimited). |

2. Run:

   ```sh
   go run ./pipeline/cmd/giab-concordance -config run.json -out ./reports -v
   # or build the binary:
   go build -o giab-concordance ./pipeline/cmd/giab-concordance
   ./giab-concordance -config run.json -out ./reports
   ```

## Outputs and how to read them

Two files land in `out_dir`:

- **`giab_concordance.json`** — the full structured result: per-stage status,
  the `concordance` block (`common`, `identical`, `differ`, `qual_ulp_only`,
  `genotype_or_filter_flips`, and a capped list of differing sites with their
  classification), and the per-stratum `biological` metrics for ours and
  upstream.
- **`giab_concordance.md`** — the human report: a stage-status table; the
  ours-vs-upstream concordance table with the **ULP-flip result** sentence; and
  a per-stratum **precision / recall / F1, ours vs upstream side by side**
  table for SNVs and indels.

How to read the key numbers:

- `genotype_or_filter_flips == 0` is the pass criterion for the record-exact
  claim: every difference between the two call sets is a QUAL/PL last-place
  wobble that leaves the genotype and the PASS/FAIL verdict intact. A non-zero
  value FAILs the run (exit 1) and the flipping sites are listed.
- In the biological table, ours and upstream P/R/F1 should match to the printed
  precision per variant type per stratum — equal accuracy against the GIAB
  truth, including in the hard strata (CMRG, difficult regions).

### Exit codes

| Code | Meaning |
|------|---------|
| `0` | All stages PASS or SKIP (no flips, no benchmarking error). Includes the no-data CI case. |
| `1` | A genotype/PASS-FAIL flip in the ours-vs-upstream comparison, or a benchmarking-engine error. |
| `2` | Usage / I/O error (bad config, unwritable output dir, …). |

## CI / no-data behaviour

The package unit tests (`go test ./pipeline/giab/...`) run with **no external
data**: the VCF record comparator, the ULP-flip detector, the BED region set,
the `hap.py`/`vcfeval` summary parsers (against embedded literal samples), the
config loader, and the SKIP paths are all covered with hand-built synthetic
inputs. Running the CLI with no `-config`, or with a config pointing at absent
data, SKIPs every stage and exits `0`.
