#!/usr/bin/env nextflow

/*
 * realbench — real-data parity + performance benchmark for the Go ports.
 *
 * This DSL2 pipeline stages a real GIAB HG002 / GRCh38 WGS dataset, derives
 * test tiers by REGION-SUBSETTING the whole-genome inputs with UPSTREAM tools
 * (so the per-tier inputs are independent of our code), then runs the Go
 * `realbench` matrix harness once per tier. Each tier run benchmarks OUR
 * bioinformatics ports (/opt/ours) against the UPSTREAM binaries (/opt/upstream)
 * baked into the same container, and reports parity + performance.
 *
 * There is NO synthetic data anywhere: every input is a real GIAB file and the
 * tiers are exact region subsets of it.
 *
 * Tiers:
 *   chr20  — single chromosome (small / fast smoke).
 *   exome  — exome-capture target regions (medium).
 *   wgs    — the full whole-genome inputs, unchanged (large; off by default).
 *
 * See test/nextflow/README.md for the launch runbook (Seqera Platform / AWS Batch).
 */

nextflow.enable.dsl = 2

/* -------------------------------------------------------------------------- *
 *  Parameters (all overridable on the CLI or via -params-file / a profile).
 *  The GIAB HG002 / GRCh38 defaults live in conf/hg002_grch38.config.
 * -------------------------------------------------------------------------- */

// Reference + whole-genome inputs. Provide EITHER wgs_bam (+bai) OR
// wgs_cram (+crai); if both are given the BAM is preferred for subsetting.
params.ref_fasta     = null   // GRCh38 no-alt analysis set FASTA (.fa / .fa.gz)
params.wgs_bam       = null   // HG002 GRCh38 WGS BAM
params.wgs_bai       = null   // its .bai (optional; indexed if missing)
params.wgs_cram      = null   // HG002 GRCh38 WGS CRAM (alternative to BAM)
params.wgs_crai      = null   // its .crai (optional; indexed if missing)

params.benchmark_vcf = null   // NIST v4.2.1 benchmark VCF (.vcf.gz)
params.benchmark_tbi = null   // its .tbi (optional; indexed if missing)
params.highconf_bed  = null   // NIST v4.2.1 high-confidence BED

params.fastq_r1      = null   // HG002 raw Illumina R1 (.fastq.gz)
params.fastq_r2      = null   // HG002 raw Illumina R2 (.fastq.gz)

// exome targets are derived from the gene_gff CDS features inside DERIVE_EXOME;
// no external capture BED is required.
params.gene_gff      = null   // GRCh38 gene annotation GFF3 (.gff3.gz / .gff.gz)

// Run shape.
params.tiers         = ['chr20', 'exome']   // wgs is defined but off by default
params.reps          = 3
params.outdir        = 'results'

// Optional read subsampling for the derived tiers. When set to a fraction in
// (0,1), the chr20/exome BAM slice is downsampled with `samtools view -s` so
// the whole matrix stays fast on real data. The GIAB 300x source is far too
// deep for a quick smoke: 0.0167 -> ~5x, 0.0333 -> ~10x, 0.1 -> ~30x. Leave
// null (or 1) to keep full depth (the real benchmark). Applied in DERIVE_*.
params.subsample     = null

// Container + binary directories (defaults are the in-image install locations
// built by the Dockerfile; override only if you mount alternative binaries).
params.container     = 'ghcr.io/yassines/bio-ai-realbench:latest'
params.our_bin       = '/opt/ours'
params.upstream_bin  = '/opt/upstream'

// Contig naming. GRCh38 no-alt uses 'chr20'; switch to '20' for a GRCh37/b37
// reference. The VCF and BED contig names must match the BAM/reference.
params.chrom         = 'chr20'

/* -------------------------------------------------------------------------- *
 *  Helpers
 * -------------------------------------------------------------------------- */

// isSet reports whether a path param is a real value (not null, empty, or the
// literal string "null"/"false" that can leak in from CLI/params-file overrides).
def isSet(p) { p != null && p.toString().trim() && !(p.toString().trim().toLowerCase() in ['null', 'false']) }

// Resolve a path param to a Nextflow file() with checkIfExists, or null when the
// param is unset. s3://, https://, ftp:// and local paths are all handled by
// Nextflow's file().
def optFile(p) { isSet(p) ? file(p, checkIfExists: true) : null }

/* -------------------------------------------------------------------------- *
 *  STAGE_REF — localise the reference FASTA and ensure it has a .fai.
 *
 *  Indexing uses the UPSTREAM samtools in the container so the staged inputs
 *  never depend on our code. Handles a bgzipped reference transparently
 *  (samtools faidx indexes .fa.gz in place producing .fai + .gzi).
 * -------------------------------------------------------------------------- */
process STAGE_REF {
    tag 'reference'

    input:
    path ref

    output:
    tuple path(ref), path("${ref}.fai"), emit: ref

    script:
    """
    set -euo pipefail
    # samtools faidx works on plain .fa and bgzipped .fa.gz alike. If a .fai
    # already came along with the file Nextflow staged it; only build if absent.
    if [ ! -s "${ref}.fai" ]; then
        ${params.upstream_bin}/samtools faidx "${ref}"
    fi
    """
}

/* -------------------------------------------------------------------------- *
 *  STAGE_ALIGN — localise the WGS BAM/CRAM and ensure it is indexed.
 *
 *  Accepts a BAM or a CRAM. Emits a tuple (aln, index, kind) where kind is
 *  'bam' or 'cram'. Indexing uses the upstream samtools.
 * -------------------------------------------------------------------------- */
process STAGE_ALIGN {
    tag { aln.name }

    input:
    tuple path(aln), path(idx), val(kind)
    tuple path(ref), path(fai)

    output:
    tuple path(aln), path("${aln}.${kind == 'cram' ? 'crai' : 'bai'}"), val(kind), emit: aln

    script:
    def suffix = kind == 'cram' ? 'crai' : 'bai'
    """
    set -euo pipefail
    # A CRAM index needs the reference on PATH (REF_PATH=/dev/null is set for
    # reference-free decode; -T below is belt-and-braces for older htslib).
    if [ ! -s "${aln}.${suffix}" ]; then
        if [ "${kind}" = "cram" ]; then
            ${params.upstream_bin}/samtools index -@ ${task.cpus} "${aln}"
        else
            ${params.upstream_bin}/samtools index -@ ${task.cpus} "${aln}"
        fi
    fi
    """
}

/* -------------------------------------------------------------------------- *
 *  STAGE_VCF — localise the benchmark VCF and ensure it has a .tbi.
 * -------------------------------------------------------------------------- */
process STAGE_VCF {
    tag { vcf.name }

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path(vcf), path("${vcf}.tbi"), emit: vcf

    script:
    """
    set -euo pipefail
    if [ ! -s "${vcf}.tbi" ]; then
        ${params.upstream_bin}/tabix -p vcf "${vcf}"
    fi
    """
}

/* -------------------------------------------------------------------------- *
 *  STAGE_GFF — localise the gene GFF and ensure it is bgzipped + tabix-indexed
 *  (needed so DERIVE_TIER can tabix-subset it by region).
 * -------------------------------------------------------------------------- */
process STAGE_GFF {
    tag { gff.name }

    input:
    path gff

    output:
    tuple path("annot.sorted.gff.gz"), path("annot.sorted.gff.gz.tbi"), emit: gff

    script:
    """
    set -euo pipefail
    # Normalise to a sorted, bgzipped, tabix-indexed GFF3 regardless of how the
    # input arrived (plain, gzip, or already bgzipped). zcat -f reads all three.
    zcat -f "${gff}" \
        | grep -v '^#' \
        | sort -t\$'\\t' -k1,1 -k4,4n \
        | ${params.upstream_bin}/bgzip -@ ${task.cpus} > annot.sorted.gff.gz
    ${params.upstream_bin}/tabix -p gff annot.sorted.gff.gz
    """
}

/* -------------------------------------------------------------------------- *
 *  DERIVE_CHR20 — region-subset every input to a single chromosome with the
 *  UPSTREAM tools. Emits the 8-element file list realbench consumes.
 * -------------------------------------------------------------------------- */
process DERIVE_CHR20 {
    tag 'chr20'
    label 'derive'

    input:
    tuple path(ref), path(fai)
    tuple path(aln), path(alnidx), val(kind)
    tuple path(vcf), path(vcftbi)
    path  bed_in            // high-confidence (or exome) BED to subset
    tuple path(gff), path(gfftbi)

    output:
    tuple val('chr20'),
          path('chr20.ref.fa'), path('chr20.ref.fa.fai'),
          path('chr20.bam'),    path('chr20.bam.bai'),
          path('chr20.cram'),   path('chr20.cram.crai'),
          path('chr20.vcf.gz'), path('chr20.vcf.gz.tbi'),
          path('chr20.R1.fq.gz'), path('chr20.R2.fq.gz'),
          path('chr20.bed'),
          path('chr20.gff.gz'), path('chr20.gff.gz.tbi'),
          emit: tier

    script:
    def C  = params.chrom
    def ST = "${params.upstream_bin}/samtools"
    def BT = "${params.upstream_bin}/bcftools"
    def TB = "${params.upstream_bin}/tabix"
    def BG = "${params.upstream_bin}/bgzip"
    def SS = (params.subsample && "${params.subsample}".isNumber() && ("${params.subsample}" as float) > 0 && ("${params.subsample}" as float) < 1) ? "-s ${params.subsample}" : ""
    """
    set -euo pipefail

    # --- reference: just the one contig (+ .fai) -------------------------------
    ${ST} faidx "${ref}" "${C}" > chr20.ref.fa
    ${ST} faidx chr20.ref.fa

    # --- BAM: region-slice (optionally subsampled with ${SS}), then SORT + index.
    # The GIAB novoalign source is tagged @HD SO:unsorted even though its reads
    # are coordinate-ordered; samtools view copies that tag, which would force
    # downstream tools (ours and upstream) onto slow, memory-heavy buffered paths.
    # Sorting stamps a proper SO:coordinate header so every cell streams. --------
    if [ "${kind}" = "cram" ]; then
        ${ST} view -@ ${task.cpus} ${SS} -T "${ref}" -u "${aln}" "${C}" | ${ST} sort -@ ${task.cpus} -o chr20.bam
    else
        ${ST} view -@ ${task.cpus} ${SS} -u "${aln}" "${C}" | ${ST} sort -@ ${task.cpus} -o chr20.bam
    fi
    ${ST} index -@ ${task.cpus} chr20.bam

    # --- CRAM: transcode the chr20 BAM against the chr20 reference -------------
    ${ST} view -@ ${task.cpus} -C -T chr20.ref.fa chr20.bam -o chr20.cram
    ${ST} index -@ ${task.cpus} chr20.cram

    # --- VCF: region view + tabix ---------------------------------------------
    ${BT} view -r "${C}" "${vcf}" -O z -o chr20.vcf.gz
    ${TB} -p vcf chr20.vcf.gz

    # --- FASTQ: the chr20 reads, collated, as paired R1/R2 --------------------
    # collate gives name-grouped pairs so samtools fastq emits a proper pair.
    ${ST} collate -@ ${task.cpus} -u -O chr20.bam \
        | ${ST} fastq -@ ${task.cpus} -1 chr20.R1.fq.gz -2 chr20.R2.fq.gz \
              -0 /dev/null -s /dev/null -n

    # --- BED: keep only the chr20 records -------------------------------------
    awk -v c="${C}" 'BEGIN{OFS="\\t"} \$1==c' "${bed_in}" > chr20.bed

    # --- GFF: tabix-subset to chr20, re-bgzip + index -------------------------
    ${TB} "${gff}" "${C}" | ${BG} -@ ${task.cpus} > chr20.gff.gz
    ${TB} -p gff chr20.gff.gz
    """
}

/* -------------------------------------------------------------------------- *
 *  DERIVE_EXOME — region-subset every input to the exome capture targets.
 *  The reference stays whole-genome (exome calling needs the full reference).
 * -------------------------------------------------------------------------- */
process DERIVE_EXOME {
    tag 'exome'
    label 'derive'

    input:
    tuple path(ref), path(fai)
    tuple path(aln), path(alnidx), val(kind)
    tuple path(vcf), path(vcftbi)
    tuple path(gff), path(gfftbi)

    output:
    tuple val('exome'),
          path(ref), path(fai),
          path('exome.bam'),    path('exome.bam.bai'),
          path('exome.cram'),   path('exome.cram.crai'),
          path('exome.vcf.gz'), path('exome.vcf.gz.tbi'),
          path('exome.R1.fq.gz'), path('exome.R2.fq.gz'),
          path('exome.bed'),
          path('exome.gff.gz'), path('exome.gff.gz.tbi'),
          emit: tier

    script:
    def ST = "${params.upstream_bin}/samtools"
    def BT = "${params.upstream_bin}/bcftools"
    def TB = "${params.upstream_bin}/tabix"
    def BG = "${params.upstream_bin}/bgzip"
    def SS = (params.subsample && "${params.subsample}".isNumber() && ("${params.subsample}" as float) > 0 && ("${params.subsample}" as float) < 1) ? "-s ${params.subsample}" : ""
    """
    set -euo pipefail

    # Exome targets = merged CDS intervals from the gene-annotation GFF. No
    # external capture BED is needed; the targets are reproducible from the
    # staged GFF. NOTE: the GFF's contig names must match the reference (a
    # chr-named Gencode GRCh38 GFF3; a RefSeq NC_-accession GFF would need the
    # assembly-report rename first — see conf/hg002_grch38.config).
    zcat -f "${gff}" \
        | awk -F'\\t' 'BEGIN{OFS="\\t"} !/^#/ && \$3=="CDS" {print \$1, \$4-1, \$5}' \
        | sort -k1,1 -k2,2n \
        | ${params.upstream_bin}/bedtools merge -i - > exome.bed
    if [ ! -s exome.bed ]; then
        echo "DERIVE_EXOME: no CDS intervals derived from the GFF — check that the GFF contigs match the reference (chr-named Gencode, not RefSeq NC_)." >&2
        exit 1
    fi

    # --- BAM: region-slice over the exome targets (optionally subsampled), then
    # SORT + index. The GIAB source is tagged SO:unsorted despite being
    # coordinate-ordered; sorting stamps a proper SO:coordinate header so
    # downstream cells stream instead of buffering the whole slice in memory. ----
    if [ "${kind}" = "cram" ]; then
        ${ST} view -@ ${task.cpus} ${SS} -T "${ref}" -u -L exome.bed "${aln}" | ${ST} sort -@ ${task.cpus} -o exome.bam
    else
        ${ST} view -@ ${task.cpus} ${SS} -u -L exome.bed "${aln}" | ${ST} sort -@ ${task.cpus} -o exome.bam
    fi
    ${ST} index -@ ${task.cpus} exome.bam

    # --- CRAM: transcode the exome BAM against the whole-genome reference ------
    ${ST} view -@ ${task.cpus} -C -T "${ref}" exome.bam -o exome.cram
    ${ST} index -@ ${task.cpus} exome.cram

    # --- VCF: restrict to the exome regions -----------------------------------
    ${BT} view -R exome.bed "${vcf}" -O z -o exome.vcf.gz
    ${TB} -p vcf exome.vcf.gz

    # --- FASTQ: the exome reads as paired R1/R2 -------------------------------
    ${ST} collate -@ ${task.cpus} -u -O exome.bam \
        | ${ST} fastq -@ ${task.cpus} -1 exome.R1.fq.gz -2 exome.R2.fq.gz \
              -0 /dev/null -s /dev/null -n

    # --- GFF: subset to the exome regions, re-bgzip + index -------------------
    # tabix accepts a regions BED (-R) to slice the gff to exome targets.
    ${TB} -R exome.bed "${gff}" | ${BG} -@ ${task.cpus} > exome.gff.gz
    ${TB} -p gff exome.gff.gz
    """
}

/* -------------------------------------------------------------------------- *
 *  DERIVE_WGS — pass the full inputs through unchanged. Re-indexes are no-ops
 *  because STAGE_* already produced the indexes; we just rename to the tier
 *  layout and transcode a whole-genome CRAM (only here because wgs needs both
 *  a BAM and a CRAM cell, and only one alignment was staged).
 * -------------------------------------------------------------------------- */
process DERIVE_WGS {
    tag 'wgs'
    label 'derive_wgs'

    input:
    tuple path(ref), path(fai)
    tuple path(aln), path(alnidx), val(kind)
    tuple path(vcf), path(vcftbi)
    path  bed_in
    tuple path(gff), path(gfftbi)
    tuple path(fq1), path(fq2)

    output:
    tuple val('wgs'),
          path(ref), path(fai),
          path('wgs.bam'),    path('wgs.bam.bai'),
          path('wgs.cram'),   path('wgs.cram.crai'),
          path(vcf), path(vcftbi),
          path(fq1), path(fq2),
          path('wgs.bed'),
          path(gff), path(gfftbi),
          emit: tier

    script:
    def ST = "${params.upstream_bin}/samtools"
    """
    set -euo pipefail

    cp "${bed_in}" wgs.bed

    if [ "${kind}" = "cram" ]; then
        # Staged input is a CRAM: it becomes the cram cell; transcode to BAM.
        ln -sf "${aln}" wgs.cram
        cp "${alnidx}" wgs.cram.crai
        ${ST} view -@ ${task.cpus} -T "${ref}" -b "${aln}" -o wgs.bam
        ${ST} index -@ ${task.cpus} wgs.bam
    else
        # Staged input is a BAM: it becomes the bam cell; transcode to CRAM.
        ln -sf "${aln}" wgs.bam
        cp "${alnidx}" wgs.bam.bai
        ${ST} view -@ ${task.cpus} -C -T "${ref}" "${aln}" -o wgs.cram
        ${ST} index -@ ${task.cpus} wgs.cram
    fi
    """
}

/* -------------------------------------------------------------------------- *
 *  RUN_MATRIX — run the Go realbench harness for one tier. It benchmarks OUR
 *  ports vs UPSTREAM across the samtools/bcftools/bedtools/seqtk/... matrix and
 *  writes realbench.<tier>.json/.md.
 * -------------------------------------------------------------------------- */
process RUN_MATRIX {
    tag { tier }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
               pattern: 'realbench.*'
    label 'matrix'

    // Scale resources by tier: wgs runs every cell on the whole genome and is
    // the single heaviest step. cpus/memory/time accept closures over inputs.
    cpus   { tier == 'wgs' ? 8 : 4 }
    memory { tier == 'wgs' ? 32.GB : 8.GB }
    time   { tier == 'wgs' ? 48.h : 8.h }

    input:
    tuple val(tier),
          path(ref), path(fai),
          path(bam), path(bai),
          path(cram), path(crai),
          path(vcf), path(vcftbi),
          path(fq1), path(fq2),
          path(bed),
          path(gff), path(gfftbi)

    output:
    tuple val(tier), path("realbench.${tier}.json"), emit: json
    path "realbench.${tier}.md", emit: md

    script:
    """
    set -euo pipefail
    realbench \
        -tier "${tier}" \
        -ref "${ref}" \
        -bam "${bam}" \
        -cram "${cram}" \
        -vcf "${vcf}" \
        -fastq1 "${fq1}" \
        -fastq2 "${fq2}" \
        -bed "${bed}" \
        -gff "${gff}" \
        -our-bin "${params.our_bin}" \
        -upstream-bin "${params.upstream_bin}" \
        -reps ${params.reps} \
        -tmp . \
        -report-only \
        -out .
    """
}

/* -------------------------------------------------------------------------- *
 *  AGGREGATE — merge the per-tier realbench JSONs into one summary JSON and a
 *  combined Markdown table. Uses jq (in the container) so it is independent of
 *  the harness internals.
 * -------------------------------------------------------------------------- */
process AGGREGATE {
    tag 'aggregate'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path 'realbench.*.json'

    output:
    path 'realbench.summary.json'
    path 'realbench.summary.md'

    script:
    """
    set -euo pipefail

    # Merge every per-tier JSON into one array keyed by tier. The realbench JSON
    # is treated opaquely: we slurp them all into a {tiers:[...]} envelope so the
    # summary survives schema changes in the harness output.
    jq -s '{generated: (now | todateiso8601), tiers: .}' realbench.*.json \
        > realbench.summary.json

    # Combined Markdown. One row per (tier, cell) for the realbench schema
    # (fields: tier, cells[], pass/diff/skip/error; cell: name/tool/parity/
    # wall_x/cpu_x/rss_x). Field accessors fall back gracefully so the table
    # survives a casing/field tweak in the harness output.
    {
        echo "# realbench combined summary"
        echo
        echo "| Tier | Cell | Tool | Parity | wall× | CPU× | RSS× |"
        echo "|------|------|------|--------|-------|------|------|"
        jq -r '
            (.tier // .Tier // "?") as \$t
            | (.cells // .Cells // [])[]
            | [ \$t,
                (.name   // .Name   // .subcommand // "?"),
                (.tool   // .Tool   // "?"),
                (.parity // .Parity // .status // "?"),
                ((.wall_x // .WallX // .wall_ratio) | if . == null then "—" else (.*100|round/100|tostring) end),
                ((.cpu_x  // .CPUX  // .cpu_ratio)  | if . == null then "—" else (.*100|round/100|tostring) end),
                ((.rss_x  // .RSSX  // .rss_ratio)  | if . == null then "—" else (.*100|round/100|tostring) end)
              ]
            | "| " + join(" | ") + " |"
        ' realbench.*.json
        echo
        echo "## Per-tier verdicts"
        echo
        echo "| Tier | Verdict | PASS | DIFF | SKIP | ERROR |"
        echo "|------|---------|------|------|------|-------|"
        jq -r '
            (.pass    // .Pass    // 0)    as \$p
            | (.diff    // .Diff    // .diverge // 0) as \$d
            | (.skip    // .Skip    // 0)  as \$s
            | (.error   // .Errored // 0)  as \$e
            | (if (\$d + \$e) > 0 then "FAIL" elif \$p > 0 then "PASS" else "NO-DATA" end) as \$v
            | [ (.tier // .Tier // "?"), \$v, \$p, \$d, \$s, \$e ]
            | "| " + (map(tostring) | join(" | ")) + " |"
        ' realbench.*.json
    } > realbench.summary.md
    """
}

/* -------------------------------------------------------------------------- *
 *  Workflow wiring
 * -------------------------------------------------------------------------- */
workflow REALBENCH {

    main:

    // --- required inputs -----------------------------------------------------
    if (!isSet(params.ref_fasta)) { error "params.ref_fasta is required (the GRCh38 reference FASTA)" }
    if (!isSet(params.wgs_bam) && !isSet(params.wgs_cram)) {
        error "provide either params.wgs_bam (+bai) or params.wgs_cram (+crai)"
    }

    ref_ch = channel.value(optFile(params.ref_fasta))
    STAGE_REF(ref_ch)

    // Pick BAM or CRAM as the single staged alignment (BAM preferred).
    if (isSet(params.wgs_bam)) {
        aln_in = channel.value(
            tuple(optFile(params.wgs_bam),
                  isSet(params.wgs_bai) ? optFile(params.wgs_bai) : optFile(params.wgs_bam),
                  'bam')
        )
    } else {
        aln_in = channel.value(
            tuple(optFile(params.wgs_cram),
                  isSet(params.wgs_crai) ? optFile(params.wgs_crai) : optFile(params.wgs_cram),
                  'cram')
        )
    }
    // The reference is a second input so STAGE_ALIGN can index a CRAM in place.
    STAGE_ALIGN(aln_in, STAGE_REF.out.ref)

    if (!isSet(params.benchmark_vcf)) { error "params.benchmark_vcf is required (the NIST benchmark VCF)" }
    vcf_in = channel.value(
        tuple(optFile(params.benchmark_vcf),
              isSet(params.benchmark_tbi) ? optFile(params.benchmark_tbi) : optFile(params.benchmark_vcf))
    )
    STAGE_VCF(vcf_in)

    if (!isSet(params.gene_gff)) { error "params.gene_gff is required (the GRCh38 gene annotation GFF3)" }
    STAGE_GFF(channel.value(optFile(params.gene_gff)))

    // The high-confidence BED drives the chr20/wgs tiers; the exome tier derives
    // its targets from the gene_gff CDS features (see DERIVE_EXOME).
    if (!isSet(params.highconf_bed)) { error "params.highconf_bed is required (the NIST high-confidence BED)" }
    highconf_ch = channel.value(optFile(params.highconf_bed))

    // Convenience handles.
    ref_out   = STAGE_REF.out.ref
    aln_out   = STAGE_ALIGN.out.aln
    vcf_out   = STAGE_VCF.out.vcf
    gff_out   = STAGE_GFF.out.gff

    // The raw FASTQs are ONLY consumed by the wgs tier (chr20/exome derive their
    // FASTQs from the region-sliced BAM). Build + validate the FASTQ channel
    // lazily so a chr20/exome-only run never touches the (large, optional) raw
    // FASTQ inputs.
    if (params.tiers.contains('wgs')) {
        if (!isSet(params.fastq_r1) || !isSet(params.fastq_r2)) {
            error "tier 'wgs' requested but params.fastq_r1/params.fastq_r2 are not set"
        }
        fq_in = channel.value(tuple(optFile(params.fastq_r1), optFile(params.fastq_r2)))
    } else {
        fq_in = channel.empty()
    }

    // --- derive the requested tiers -----------------------------------------
    // Each DERIVE_<tier> emits one tuple(tier, [ref, bam, cram, vcf, fq1, fq2,
    // bed, gff, indexes]); mix the selected tiers into a single channel so a
    // single RUN_MATRIX instance fans out one job per tier.
    tier_ch = channel.empty()
    matched = false

    if (params.tiers.contains('chr20')) {
        DERIVE_CHR20(ref_out, aln_out, vcf_out, highconf_ch, gff_out)
        tier_ch = tier_ch.mix(DERIVE_CHR20.out.tier)
        matched = true
    }

    if (params.tiers.contains('exome')) {
        // Exome targets are derived from the gene-annotation GFF's CDS features
        // inside DERIVE_EXOME (no external capture BED needed). params.gene_gff
        // is already required above, and its contigs MUST match the reference
        // (use a chr-named Gencode GRCh38 GFF3, not a RefSeq NC_-accession GFF).
        DERIVE_EXOME(ref_out, aln_out, vcf_out, gff_out)
        tier_ch = tier_ch.mix(DERIVE_EXOME.out.tier)
        matched = true
    }

    if (params.tiers.contains('wgs')) {
        DERIVE_WGS(ref_out, aln_out, vcf_out, highconf_ch, gff_out, fq_in)
        tier_ch = tier_ch.mix(DERIVE_WGS.out.tier)
        matched = true
    }

    if (!matched) {
        error "no recognised tiers in params.tiers=${params.tiers} (expected chr20/exome/wgs)"
    }

    // --- run the matrix per tier --------------------------------------------
    RUN_MATRIX(tier_ch)

    // --- aggregate ----------------------------------------------------------
    AGGREGATE(RUN_MATRIX.out.json.map { _tier, j -> j }.collect())

    log.info "realbench tiers=${params.tiers} reps=${params.reps} -> ${params.outdir}"
}

// Anonymous entry point so `nextflow run test/nextflow/main.nf` works directly
// (the root ./main.nf wrapper imports REALBENCH for Seqera Platform, which
// requires a repo-root main.nf).
workflow {
    REALBENCH()
}
