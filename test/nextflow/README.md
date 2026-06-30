# realbench — real-data parity + performance pipeline (Nextflow)

A DSL2 Nextflow pipeline that benchmarks **our** Go bioinformatics ports against
the **upstream** originals on a **real GIAB HG002 / GRCh38 WGS dataset**, and
reports both **parity** (byte-exact-after-provenance-stripping) and
**performance** (wall / CPU / RSS, ours-vs-upstream ratios).

There is **no synthetic data anywhere**. The pipeline stages a real GIAB
whole-genome dataset and derives smaller test tiers by **region-subsetting the
real inputs with the upstream tools** — so the per-tier inputs are independent
of our code (the thing under test never produces its own benchmark inputs).

It is launchable on the **Seqera Platform / Tower** (AWS Batch compute
environment) and locally with Docker.

## How it fits together

```
STAGE_REF  ┐
STAGE_ALIGN├─► DERIVE_<tier> ──► RUN_MATRIX(realbench) ──► AGGREGATE
STAGE_VCF  │      (upstream             (Go harness:           (jq merge →
STAGE_GFF  ┘       subsetting)           ours vs upstream)      summary.json/.md)
```

- **STAGE_REF / STAGE_ALIGN / STAGE_VCF / STAGE_GFF** localise the inputs
  (`s3://`, `https://`, `ftp://`, or local) and build any missing index
  (`samtools faidx`, `samtools index`, `tabix -p vcf`, `tabix -p gff`) using the
  **upstream** binaries in the container.
- **DERIVE_\<tier\>** region-subsets every input with **upstream** `samtools` /
  `bcftools` / `tabix` / `bgzip` and emits the 8-input file set the harness
  consumes: `(tier, [ref, bam, cram, vcf, fastq1, fastq2, bed, gff])`.
- **RUN_MATRIX** runs the Go matrix harness once per tier:

  ```sh
  realbench -tier <chr20|exome|wgs> -ref <fa> -bam <bam> -cram <cram> \
            -vcf <vcf.gz> -fastq1 <R1> -fastq2 <R2> -bed <bed> -gff <gff.gz> \
            -our-bin /opt/ours -upstream-bin /opt/upstream -reps <N> -out .
  ```

  publishing `realbench.<tier>.json` and `realbench.<tier>.md`.
- **AGGREGATE** merges the per-tier JSONs into `realbench.summary.json` and a
  combined `realbench.summary.md` table with `jq`.

> The Go `realbench` harness is built by a sibling component into
> `pipeline/cmd/realbench` and baked into the container at
> `/usr/local/bin/realbench`. This pipeline only orchestrates staging, tiering,
> and aggregation around it.

## Tiers

| Tier    | Size   | What it is | Default |
|---------|--------|------------|---------|
| `chr20` | small  | every input sliced to chromosome 20 | **on** |
| `exome` | medium | every input restricted to the exome targets — **merged CDS intervals derived from `params.gene_gff`** (no capture BED needed); reference stays whole-genome | **on** |
| `wgs`   | large  | the full whole-genome inputs, unchanged | **off** |

All three tiers come from the **same** WGS source — chr20 and exome are exact
region subsets of the whole-genome BAM/VCF/GFF. `wgs` is **defined but off by
default** (`params.tiers = ['chr20','exome']`) because the full run transcodes
and benchmarks against a ~300x whole-genome BAM and is expensive. Enable it
explicitly:

```sh
nextflow run test/nextflow/main.nf --tiers chr20,exome,wgs ...
```

## The container

One image carries both sides on the same architecture:

- `/opt/ours` — our Go binaries (`tools/<tool>/cmd/<tool>`).
- `/opt/upstream` — upstream oracles built from the `reference_code/`
  submodules: htslib (`bgzip`/`tabix`/`htsfile`), `samtools`, `bcftools`,
  `bedtools`, `seqtk`, `sickle`, `skewer`, `fastp`, `vcftools`, and `prinseq`
  (a Perl wrapper).
- `/usr/local/bin/realbench` — the Go matrix harness.
- `awscli`, `jq`, `bgzip`, `tabix` on `PATH`.

Build it from the **repo root** (the build context needs `tools/`, `pkg/`,
`pipeline/`, and the `reference_code/` submodules):

```sh
git submodule update --init \
    reference_code/htslib reference_code/samtools reference_code/bcftools \
    reference_code/bedtools reference_code/seqtk reference_code/sickle \
    reference_code/skewer reference_code/fastp reference_code/vcftools \
    reference_code/prinseq reference_code/patches

docker build -f test/nextflow/Dockerfile -t ghcr.io/yassines/bio-ai-realbench:latest .
docker push ghcr.io/yassines/bio-ai-realbench:latest   # for AWS Batch / Seqera
```

For AWS Batch always build/push **linux/amd64** so ours and upstream are
compared on the same hardware as the Batch nodes:

```sh
docker build --platform=linux/amd64 -f test/nextflow/Dockerfile \
    -t ghcr.io/yassines/bio-ai-realbench:latest .
```

> **mosdepth** is not built into the image (it is a Nim release-only binary). If
> you want the mosdepth cells, drop a `linux/amd64` `mosdepth` into
> `/opt/upstream/mosdepth` at runtime or set `MOSDEPTH_BIN`.

## Running it

### Local smoke test (Docker)

The `test` profile runs the **chr20 tier only**, 1 rep, against the smallest
real WGS alignment GIAB publishes (the **chr20-only 300x BAM**, ~11 GB):

```sh
nextflow run test/nextflow/main.nf -profile test,standard
```

This still pulls real GIAB data over the network. For a fully offline run, stage
the chr20 BAM + GRCh38 reference + benchmark VCF locally and override the params:

```sh
nextflow run test/nextflow/main.nf -profile standard \
    --ref_fasta     /data/GRCh38_no_alt.fa \
    --wgs_bam       /data/HG002.GRCh38.300x_chr20.bam \
    --benchmark_vcf /data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    --highconf_bed  /data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    --gene_gff      /data/GCF_000001405.40_GRCh38.p14_genomic.gff.gz \
    --tiers chr20 --reps 1
```

### Seqera Platform / Tower (AWS Batch)

1. **Compute environment.** In the Seqera Platform, create an **AWS Batch**
   compute environment whose job role has **read** access to `s3://giab` (the
   GIAB bucket is public — an unsigned/anonymous read works, but the role must
   not be blocked from it) and **read/write** on your own work/results bucket.
   Set the Nextflow **work directory** to `s3://<your-bucket>/work`.

2. **Container.** Point the pipeline at the pushed image. It is already the
   default (`params.container`); override per-launch if you publish elsewhere.

3. **Launch from the CLI** with the Tower CLI:

   ```sh
   tw launch https://github.com/yassineS/bio_ai_experiment \
       --compute-env  <your-aws-batch-ce> \
       --work-dir     s3://<your-bucket>/work \
       --revision     claude/realbench-harness \
       --main-script  test/nextflow/main.nf \
       --profile      awsbatch \
       --params-file  test/nextflow/conf/params.example.yaml \
       --config       test/nextflow/nextflow.config
   ```

   or **from the Platform UI**: *Launchpad → Add pipeline →* point at this repo,
   set `test/nextflow/main.nf` as the main script, choose the `awsbatch` profile, and
   paste/upload your params.

4. **Nothing Tower-specific is needed in-repo** — the pipeline is plain DSL2.
   The `awsbatch` profile sets `process.executor='awsbatch'`, `aws.region`, and
   the Batch CLI path; the Platform injects the queue/role from the compute
   environment. `s3://` inputs (e.g. `s3://giab/...`) are staged natively by
   Nextflow on AWS Batch.

A params file looks like:

```yaml
# test/nextflow/conf/params.example.yaml
tiers:        [chr20, exome]
reps:         3
ref_fasta:    s3://giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
wgs_bam:      s3://giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam
benchmark_vcf: s3://giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
highconf_bed:  s3://giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
gene_gff:      https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gff3.gz  # chr-named (drives exome CDS derivation)
outdir:        s3://<your-bucket>/results
```

## Data: GIAB HG002 / GRCh38

The default dataset is the GIAB Ashkenazim-trio **son** (HG002 / NA24385) on the
**GRCh38 no-alt analysis set**, with the **NIST v4.2.1** benchmark. The canonical
paths are in [`conf/hg002_grch38.config`](conf/hg002_grch38.config).

> **S3 keys drift — verify them.** Every default in `hg002_grch38.config` is
> tagged `verify this key`. Before a real run, confirm the objects exist:
>
> ```sh
> aws s3 ls --no-sign-request s3://giab/release/references/GRCh38/
> aws s3 ls --no-sign-request \
>     s3://giab/data/AshkenazimTrio/HG002_NA24385_son/
> aws s3 ls --no-sign-request \
>     s3://giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/
> ```
>
> **S3 vs HTTPS prefix gotcha:** the public bucket is `s3://giab/` with keys
> `release/...` and `data/...`. The NCBI **HTTPS** mirror inserts an extra `ftp/`
> segment (`https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/...`) and the newer
> `ReferenceSamples/giab/...` mirror uses yet another prefix — same files,
> different paths. Use `s3://` on AWS Batch; use an `https://` form for local
> runs.

One input you **must** verify yourself:

- **`gene_gff`** — a GRCh38 gene annotation GFF3. It drives **both** the
  `csq`/`subseq` cells **and** the `exome` tier (whose targets are the merged
  CDS intervals derived from this GFF — no separate capture BED is needed). It
  **must be `chr`-named** to match the GRCh38 no-alt reference: the default is a
  **Gencode** GFF3 (chr-named). A RefSeq genomic GFF uses `NC_` accessions and
  would derive zero exome intervals unless you rename its contigs first.

## Outputs

Published under `params.outdir`:

```
results/
├── realbench.chr20.json   realbench.chr20.md
├── realbench.exome.json   realbench.exome.md
├── realbench.summary.json combined machine-readable summary
├── realbench.summary.md   combined parity + ratio table
├── staged/                localised + indexed inputs
├── tiers/<tier>/          the derived per-tier subsets
└── pipeline_info/         timeline / report / trace / dag
```

Read `realbench.summary.md` first: it has one row per `(tier, cell)` with the
parity status and the ours/upstream wall/CPU/RSS ratios (`<1` = ours faster /
lighter), plus a per-tier verdict table.

## Resource model

Sized in [`nextflow.config`](nextflow.config) for a GIAB-class WGS BAM:

| Step | cpus | memory | time | notes |
|------|------|--------|------|-------|
| `STAGE_ALIGN` | 4 | 8 GB | 8 h | staging+indexing a ~300x WGS BAM is slow |
| `DERIVE_CHR20` / `DERIVE_EXOME` | 4 | 8 GB | 6 h | region slice + transcode (label `derive`) |
| `DERIVE_WGS` | 8 | 16 GB | 24 h | whole-genome BAM-CRAM transcode (label `derive_wgs`) |
| `RUN_MATRIX` (chr20/exome) | 4 | 8 GB | 8 h | the matrix harness (label `matrix`) |
| `RUN_MATRIX` (wgs) | 8 | 32 GB | 48 h | heaviest step; scaled in-process by a per-tier `cpus`/`memory`/`time` closure |

Tune these on the CLI/profile to match your Batch queue's instance types.

## Caveats

- **Performance numbers need an uncontended node.** Run the matrix on a Batch
  instance with no co-tenant compute; the parity verdict is robust to contention
  but the wall/CPU/RSS ratios are not.
- **Contig naming must be consistent** across the reference, BAM, VCF, and BEDs.
  GRCh38 no-alt is `chr`-prefixed (`chr20`); set `params.chrom` accordingly (use
  `20` for a GRCh37/b37 reference).
- This pipeline does **not** replace `pipeline/bench` or `pipeline/cmd/realparity`
  — it wraps the new `realbench` harness in a cloud-launchable, tiered,
  real-data form.
