package realbench

// This file defines the full real-data benchmark+parity matrix: a
// representative-but-thorough set of cells covering every subcommand of every
// ported tool, plus the principal and memory/format-relevant flags per
// subcommand. Each cell maps to the right real input(s) for the active tier via
// the {bam}/{cram}/{ref}/{vcf}/{fastq1}/{fastq2}/{bed}/{gff} placeholders.
//
// The matrix is tier-aware only in that the SAME cells run against
// progressively larger inputs (chr20 < exome < wgs); the cell set itself is
// tier-independent, so Matrix(tier) returns the union and the runner SKIPs any
// cell whose input is absent for the chosen tier.

// Matrix returns the ordered cell battery for a tier. The tier is recorded on
// each cell; the cell set is the same across tiers (the inputs differ in size).
func Matrix(tier string) []CellSpec {
	var cells []CellSpec
	cells = append(cells, samtoolsCells()...)
	cells = append(cells, bcftoolsCells()...)
	cells = append(cells, bedCells()...)
	cells = append(cells, seqtkCells()...)
	cells = append(cells, fastpCells()...)
	cells = append(cells, trimmerCells()...)
	cells = append(cells, vcftoolsCells()...)
	cells = append(cells, mosdepthCells()...)
	cells = append(cells, htslibCells()...)
	return cells
}

// stdout is a convenience constructor for a stdout-comparison cell.
func stdout(tool, name, sub string, need InputKind, args ...string) CellSpec {
	return CellSpec{Tool: tool, Name: name, Subcommand: sub, Need: need, Post: PostStdout, OurArgs: args}
}

// samtoolsCells covers the samtools subcommand surface (BAM/CRAM, ref for CRAM).
func samtoolsCells() []CellSpec {
	const t = "samtools"
	return []CellSpec{
		// view: the main decode/filter/convert surface.
		stdout(t, "samtools_view_sam", "view", NeedBAM, "view", phBAM),
		stdout(t, "samtools_view_sam_header", "view -h", NeedBAM, "view", "-h", phBAM),
		stdout(t, "samtools_view_count", "view -c", NeedBAM, "view", "-c", phBAM),
		stdout(t, "samtools_view_f", "view -f", NeedBAM, "view", "-f", "0x2", phBAM),
		stdout(t, "samtools_view_F", "view -F", NeedBAM, "view", "-F", "0x100", phBAM),
		stdout(t, "samtools_view_q", "view -q", NeedBAM, "view", "-q", "20", phBAM),
		stdout(t, "samtools_view_L", "view -L bed", NeedBAM|NeedBED, "view", "-L", phBED, phBAM),
		{Tool: t, Name: "samtools_view_bam", Subcommand: "view -b", Need: NeedBAM, Post: PostViewSAM,
			OurArgs: []string{"view", "-b", "-o", phOut, phBAM}, WriteOut: ".bam"},
		{Tool: t, Name: "samtools_view_cram", Subcommand: "view -C -T", Need: NeedBAM | NeedRef, Post: PostViewSAM,
			OurArgs: []string{"view", "-C", "-T", phRef, "-o", phOut, phBAM}, WriteOut: ".cram", NeedsCRAMEnv: true},

		// sort.
		{Tool: t, Name: "samtools_sort", Subcommand: "sort", Need: NeedBAM, Post: PostViewSAM,
			OurArgs: []string{"sort", "-o", phOut, phBAM}, WriteOut: ".bam"},
		{Tool: t, Name: "samtools_sort_n", Subcommand: "sort -n", Need: NeedBAM, Post: PostViewSAM,
			OurArgs: []string{"sort", "-n", "-o", phOut, phBAM}, WriteOut: ".bam"},

		// reporting / streaming-text.
		stdout(t, "samtools_flagstat", "flagstat", NeedBAM, "flagstat", phBAM),
		stdout(t, "samtools_idxstats", "idxstats", NeedBAM, "idxstats", phBAM),
		stdout(t, "samtools_stats", "stats", NeedBAM, "stats", phBAM),
		stdout(t, "samtools_depth_a", "depth -a", NeedBAM, "depth", "-a", phBAM),
		stdout(t, "samtools_depth_b", "depth -b bed", NeedBAM|NeedBED, "depth", "-b", phBED, phBAM),
		stdout(t, "samtools_coverage", "coverage", NeedBAM, "coverage", phBAM),
		stdout(t, "samtools_mpileup", "mpileup -f", NeedBAM|NeedRef, "mpileup", "-f", phRef, phBAM),
		stdout(t, "samtools_fastq", "fastq", NeedBAM, "fastq", phBAM),

		// index (file-producing; .bai re-decode is not meaningful — compare via
		// idxstats stdout instead, which is the index's observable output).
		{Tool: t, Name: "samtools_index", Subcommand: "index", Need: NeedBAM, Post: PostNone,
			OurArgs: []string{"index", phBAM}},

		// markdup / fixmate / calmd — record-rewriting filters (re-decoded).
		// fixmate requires name-collated input (upstream errors on coord-sorted);
		// feed the derived name-collated BAM to both sides.
		{Tool: t, Name: "samtools_fixmate", Subcommand: "fixmate", Need: NeedNameSortBAM, Post: PostViewSAM,
			OurArgs: []string{"fixmate", phNameBAM, phOut}, WriteOut: ".bam"},
		// markdup requires input that was name-sorted, fixmate -m'd, then coord
		// re-sorted; feed the derived markdup-ready BAM to both sides.
		{Tool: t, Name: "samtools_markdup", Subcommand: "markdup", Need: NeedFixmateBAM, Post: PostViewSAM,
			OurArgs: []string{"markdup", phFixmateBAM, phOut}, WriteOut: ".bam"},
		{Tool: t, Name: "samtools_calmd", Subcommand: "calmd", Need: NeedBAM | NeedRef, Post: PostStdout,
			OurArgs: []string{"calmd", phBAM, phRef}},

		// consensus.
		stdout(t, "samtools_consensus_fa", "consensus -f fasta", NeedBAM, "consensus", "-f", "fasta", phBAM),
		stdout(t, "samtools_consensus_pileup", "consensus -f pileup", NeedBAM, "consensus", "-f", "pileup", phBAM),

		// metadata / index-like commands.
		stdout(t, "samtools_dict", "dict", NeedRef, "dict", phRef),
		{Tool: t, Name: "samtools_faidx", Subcommand: "faidx", Need: NeedRef, Post: PostStdout,
			OurArgs: []string{"faidx", phRef}},
		{Tool: t, Name: "samtools_quickcheck", Subcommand: "quickcheck", Need: NeedBAM, Post: PostNone,
			OurArgs: []string{"quickcheck", phBAM}},

		// addreplacerg / reheader / cat / merge / split — alignment-rewriting.
		{Tool: t, Name: "samtools_addreplacerg", Subcommand: "addreplacerg", Need: NeedBAM, Post: PostViewSAM,
			OurArgs: []string{"addreplacerg", "-r", "ID:rb\tSM:sample", "-o", phOut, phBAM}, WriteOut: ".bam"},
		{Tool: t, Name: "samtools_cat", Subcommand: "cat", Need: NeedBAM, Post: PostViewSAM,
			OurArgs: []string{"cat", "-o", phOut, phBAM}, WriteOut: ".bam"},
		{Tool: t, Name: "samtools_merge", Subcommand: "merge", Need: NeedBAM, Post: PostViewSAM,
			OurArgs: []string{"merge", "-f", phOut, phBAM}, WriteOut: ".bam"},

		// CRAM-input cells (require the tier CRAM + ref).
		stdout(t, "samtools_view_cram_in", "view (cram in)", NeedCRAM|NeedRef, "view", "-T", phRef, phCRAM),
		{Tool: t, Name: "samtools_cram_to_bam", Subcommand: "view -b (cram in)", Need: NeedCRAM | NeedRef, Post: PostViewSAM,
			OurArgs: []string{"view", "-b", "-T", phRef, "-o", phOut, phCRAM}, WriteOut: ".bam", NeedsCRAMEnv: true},
	}
}

// bcftoolsCells covers the bcftools subcommand surface (VCF.gz; some need ref).
func bcftoolsCells() []CellSpec {
	const t = "bcftools"
	return []CellSpec{
		stdout(t, "bcftools_view", "view", NeedVCF, "view", phVCF),
		stdout(t, "bcftools_view_body", "view -H", NeedVCF, "view", "-H", phVCF),
		stdout(t, "bcftools_view_region", "view -r", NeedVCF, "view", "-H", phVCF),
		stdout(t, "bcftools_view_include", "view -i", NeedVCF, "view", "-H", "-i", "QUAL>20", phVCF),
		stdout(t, "bcftools_view_exclude", "view -e", NeedVCF, "view", "-H", "-e", "QUAL<=20", phVCF),
		{Tool: t, Name: "bcftools_view_bcf", Subcommand: "view -Ob", Need: NeedVCF, Post: PostViewSAM,
			OurArgs: []string{"view", "-O", "b", "-o", phOut, phVCF}, WriteOut: ".bcf"},

		stdout(t, "bcftools_query", "query -f", NeedVCF, "query", "-f", `%CHROM\t%POS\t%REF\t%ALT\n`, phVCF),
		stdout(t, "bcftools_stats", "stats", NeedVCF, "stats", phVCF),
		stdout(t, "bcftools_norm", "norm -f -m", NeedVCF|NeedRef, "norm", "-f", phRef, "-m", "-", "-O", "v", phVCF),

		{Tool: t, Name: "bcftools_filter", Subcommand: "filter", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"filter", "-e", "QUAL<20", "-O", "v", phVCF}},
		{Tool: t, Name: "bcftools_sort", Subcommand: "sort", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"sort", "-O", "v", phVCF}},
		{Tool: t, Name: "bcftools_concat", Subcommand: "concat", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"concat", "-O", "v", phVCF}},
		{Tool: t, Name: "bcftools_annotate", Subcommand: "annotate", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"annotate", "-x", "INFO", "-O", "v", phVCF}},
		{Tool: t, Name: "bcftools_head", Subcommand: "head", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"head", phVCF}},
		// reheader needs a modification directive; a bare `reheader -o` errors on
		// upstream with usage. Feed a deterministic one-line sample-rename map (-s)
		// so both sides perform the same rename.
		{Tool: t, Name: "bcftools_reheader", Subcommand: "reheader", Need: NeedVCF | NeedSampleRename, Post: PostBgzipD,
			OurArgs: []string{"reheader", "-s", phSampleRename, "-o", phOut, phVCF}, WriteOut: ".vcf.gz"},
		{Tool: t, Name: "bcftools_convert", Subcommand: "convert", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"convert", "-O", "v", phVCF}},
		{Tool: t, Name: "bcftools_consensus", Subcommand: "consensus -f", Need: NeedVCF | NeedRef, Post: PostStdout,
			OurArgs: []string{"consensus", "-f", phRef, phVCF}},
		{Tool: t, Name: "bcftools_csq", Subcommand: "csq -g", Need: NeedVCF | NeedRef | NeedGFF, Post: PostStdout,
			OurArgs: []string{"csq", "-p", "a", "-f", phRef, "-g", phGFF, "-O", "v", phVCF}},
		{Tool: t, Name: "bcftools_roh", Subcommand: "roh", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"roh", "-G30", phVCF}},
		{Tool: t, Name: "bcftools_gtcheck", Subcommand: "gtcheck", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"gtcheck", phVCF}},
	}
}

// bedCells covers the ~41 bed* tools. Each maps to our per-bed binary and to
// `bedtools <sub>` upstream (the runner derives the upstream subcommand). The
// .fai is used as the genome file for genomecov/slop/makewindows/complement;
// the BAM for bamtobed/multicov/coverage; the GFF where relevant.
func bedCells() []CellSpec {
	type bc struct {
		tool, sub  string
		need       InputKind
		args       []string // our argv (also upstream, sans the bedtools subcommand prefix)
		ourOnly    bool     // ran ours-only as a perf cell (parity SKIP)
		workDirOut bool     // run in a per-side work dir (needed by cells with an -p prefix)
	}
	specs := []bc{
		{"bedintersect", "intersect", NeedBED, []string{"-a", phBED, "-b", phBED}, false, false},
		{"bedmerge", "merge", NeedBED, []string{"-i", phBED}, false, false},
		{"bedsort", "sort", NeedBED, []string{"-i", phBED}, false, false},
		{"bedsubtract", "subtract", NeedBED, []string{"-a", phBED, "-b", phBED}, false, false},
		{"bedwindow", "window", NeedBED, []string{"-a", phBED, "-b", phBED}, false, false},
		{"bedclosest", "closest", NeedBED, []string{"-a", phBED, "-b", phBED}, false, false},
		// BED3 has no col 4; map over col 3 (end) so this is a real BED3 comparison.
		{"bedmap", "map", NeedBED, []string{"-a", phBED, "-b", phBED, "-c", "3", "-o", "count"}, false, false},
		{"bedcoverage", "coverage", NeedBED, []string{"-a", phBED, "-b", phBED}, false, false},
		{"bedjaccard", "jaccard", NeedBED, []string{"-a", phBED, "-b", phBED}, false, false},
		{"bedfisher", "fisher", NeedBED | NeedRef, []string{"-a", phBED, "-b", phBED, "-g", phFai}, false, false},
		{"bedreldist", "reldist", NeedBED, []string{"-a", phBED, "-b", phBED}, false, false},
		{"bedspacing", "spacing", NeedBED, []string{"-i", phBED}, false, false},
		{"bedslop", "slop", NeedBED | NeedRef, []string{"-i", phBED, "-g", phFai, "-b", "10"}, false, false},
		{"bedflank", "flank", NeedBED | NeedRef, []string{"-i", phBED, "-g", phFai, "-b", "10"}, false, false},
		{"bedshift", "shift", NeedBED | NeedRef, []string{"-i", phBED, "-g", phFai, "-s", "5"}, false, false},
		{"bedcomplement", "complement", NeedBED | NeedRef, []string{"-i", phBED, "-g", phFai}, false, false},
		{"bedgenomecov", "genomecov", NeedBED | NeedRef, []string{"-i", phBED, "-g", phFai}, false, false},
		{"bedmakewindows", "makewindows", NeedRef, []string{"-g", phFai, "-w", "100000"}, false, false},
		{"bedmulticov", "multicov", NeedBED | NeedBAM, []string{"-bams", phBAM, "-bed", phBED}, false, false},
		{"bedmultiinter", "multiinter", NeedBED, []string{"-i", phBED, phBED}, false, false},
		{"bednuc", "nuc", NeedBED | NeedRef, []string{"-fi", phRef, "-bed", phBED}, false, false},
		{"bedgetfasta", "getfasta", NeedBED | NeedRef, []string{"-fi", phRef, "-bed", phBED}, false, false},
		{"bed12tobed6", "bed12tobed6", NeedBED, []string{"-i", phBED}, false, false},
		{"bedbamtobed", "bamtobed", NeedBAM, []string{"-i", phBAM}, false, false},
		// BED3 has no col 4; expand col 3 (a real, valid BED3 invocation).
		{"bedexpand", "expand", NeedBED, []string{"-i", phBED, "-c", "3"}, false, false},
		{"bedgroupby", "groupby", NeedBED, []string{"-i", phBED, "-g", "1", "-c", "2", "-o", "min"}, false, false},
		{"bedannotate", "annotate", NeedBED, []string{"-i", phBED, "-files", phBED}, false, false},
		// overlap needs four position columns on a paired/windowed file; feed the
		// derived 8-field window BED with -cols 2,3,6,7 (a real ours-vs-upstream cell).
		{"bedoverlap", "overlap", NeedBED | NeedWindow, []string{"-i", phWindow, "-cols", "2,3,6,7"}, false, false},
		{"bedsummary", "summary", NeedBED | NeedRef, []string{"-i", phBED, "-g", phFai}, false, false},
		// unionbedg needs 4-col BedGraph, not BED3 (upstream SIGABRTs on BED3);
		// feed the derived 4-col bedgraph.
		{"bedunionbedg", "unionbedg", NeedBED | NeedBedGraph, []string{"-i", phBedGraph, phBedGraph}, false, false},
		{"bedcluster", "cluster", NeedBED, []string{"-i", phBED}, false, false},
		{"bedlinks", "links", NeedBED, []string{"-i", phBED}, false, false},
		{"bedigv", "igv", NeedBED, []string{"-i", phBED}, false, false},
		// -labels and -names are mutually exclusive upstream; pass only -labels.
		{"bedtag", "tag", NeedBED | NeedBAM, []string{"-i", phBAM, "-files", phBED, "-labels", "x"}, false, false},
		// pairtopair needs BEDPE (>=10 fields) on both sides; feed the derived BEDPE.
		{"bedpairtopair", "pairtopair", NeedBED | NeedBEDPE, []string{"-a", phBEDPE, "-b", phBEDPE}, true, false},
		// pairtobed needs a BEDPE -a and a BED -b; feed the derived BEDPE + the BED3.
		{"bedpairtobed", "pairtobed", NeedBED | NeedBEDPE, []string{"-a", phBEDPE, "-b", phBED}, true, false},
		// -p is a required output PREFIX; it must live inside a real work dir, so
		// run in a per-side work dir and write the shards under {outdir}/split.
		{"bedsplit", "split", NeedBED, []string{"-i", phBED, "-n", "2", "-p", phOutdir + "/split"}, true, true},
		// tobam requires a name column (BED4+) to fill QNAME; feed the derived BED4.
		{"bedtobam", "tobam", NeedBED | NeedBED4 | NeedRef, []string{"-i", phBED4, "-g", phFai}, true, false},
		// RNG-driven tools have no deterministic upstream pair: ours-only perf.
		{"bedrandom", "random", NeedRef, []string{"-g", phFai, "-n", "100", "-seed", "1"}, true, false},
		{"bedshuffle", "shuffle", NeedBED | NeedRef, []string{"-i", phBED, "-g", phFai, "-seed", "1"}, true, false},
		{"bedsample", "sample", NeedBED, []string{"-i", phBED, "-n", "10", "-seed", "1"}, true, false},
	}
	cells := make([]CellSpec, 0, len(specs))
	for _, s := range specs {
		post := PostStdout
		if s.ourOnly {
			post = PostOursOnly
		}
		cells = append(cells, CellSpec{
			Tool:       s.tool,
			Name:       s.tool,
			Subcommand: "bedtools " + s.sub,
			Need:       s.need,
			Post:       post,
			OurArgs:    s.args,
			WorkDirOut: s.workDirOut,
			// Upstream needs the bedtools subcommand prepended; the runner does
			// that via UpStub. Same trailing args, so UpArgs == OurArgs.
		})
	}
	return cells
}

// seqtkCells covers the seqtk subcommand surface (FASTA/FASTQ).
func seqtkCells() []CellSpec {
	const t = "seqtk"
	return []CellSpec{
		stdout(t, "seqtk_seq", "seq", NeedFastq1, "seq", phFastq1),
		stdout(t, "seqtk_seq_A", "seq -A", NeedFastq1, "seq", "-A", phFastq1),
		stdout(t, "seqtk_seq_r", "seq -r", NeedFastq1, "seq", "-r", phFastq1),
		stdout(t, "seqtk_seq_mask", "seq -q -n", NeedFastq1, "seq", "-q", "20", "-n", "N", phFastq1),
		stdout(t, "seqtk_comp", "comp", NeedFastq1, "comp", phFastq1),
		stdout(t, "seqtk_sample", "sample -s", NeedFastq1, "sample", "-s", "11", phFastq1, "0.1"),
		stdout(t, "seqtk_trimfq", "trimfq", NeedFastq1, "trimfq", phFastq1),
		stdout(t, "seqtk_subseq", "subseq bed", NeedFastq1|NeedBED, "subseq", phFastq1, phBED),
		stdout(t, "seqtk_cutN", "cutN", NeedFastq1, "cutN", "-n", "1", phFastq1),
		stdout(t, "seqtk_mergepe", "mergepe", NeedFastq1|NeedFastq2, "mergepe", phFastq1, phFastq2),
		stdout(t, "seqtk_fqchk", "fqchk", NeedFastq1, "fqchk", phFastq1),
		stdout(t, "seqtk_hpc", "hpc", NeedFastq1, "hpc", phFastq1),
	}
}

// fastpCells covers fastp: a default PE run. The filtered FASTQ is the
// comparison stream (the JSON/HTML embed time/version, so they are not
// compared). fastp writes to -o/-O paths, so a per-side temp out is used.
func fastpCells() []CellSpec {
	const t = "fastp"
	return []CellSpec{
		{Tool: t, Name: "fastp_pe", Subcommand: "fastp (PE)", Need: NeedFastq1 | NeedFastq2, Post: PostBgzipD,
			OurArgs: []string{
				"-i", phFastq1, "-I", phFastq2,
				"-o", phOut, "-O", phOut + ".R2.fq.gz",
				"--json", phOut + ".json", "--html", phOut + ".html",
			}, WriteOut: ".fq.gz"},
		{Tool: t, Name: "fastp_se", Subcommand: "fastp (SE)", Need: NeedFastq1, Post: PostBgzipD,
			OurArgs: []string{
				"-i", phFastq1, "-o", phOut,
				"--json", phOut + ".json", "--html", phOut + ".html",
			}, WriteOut: ".fq.gz"},
	}
}

// trimmerCells covers sickle, skewer, and prinseq (FASTQ quality trimming/QC).
func trimmerCells() []CellSpec {
	return []CellSpec{
		// sickle pe writes -o (R1) / -p (R2) / -s (singles); compare R1.
		{Tool: "sickle", Name: "sickle_pe", Subcommand: "sickle pe", Need: NeedFastq1 | NeedFastq2, Post: PostFile,
			OurArgs: []string{"pe", "-t", "sanger", "-f", phFastq1, "-r", phFastq2,
				"-o", phOut, "-p", phOut + ".R2.fq", "-s", phOut + ".s.fq"}, WriteOut: ".fq"},
		{Tool: "sickle", Name: "sickle_se", Subcommand: "sickle se", Need: NeedFastq1, Post: PostFile,
			OurArgs: []string{"se", "-t", "sanger", "-f", phFastq1, "-o", phOut}, WriteOut: ".fq"},

		// skewer writes <prefix>-trimmed*.fastq into a work dir.
		{Tool: "skewer", Name: "skewer_pe", Subcommand: "skewer (PE)", Need: NeedFastq1 | NeedFastq2,
			Post: PostFile, WorkDirOut: true, Compare: "rb-trimmed-pair1.fastq",
			OurArgs: []string{"-o", phOutdir + "/rb", phFastq1, phFastq2}},

		// prinseq-lite: stats (stdout-ish) and a filter pass. Both cells read
		// the DECOMPRESSED plain FASTQ ({fastqplain}, derived at run start),
		// because prinseq-lite.pl 0.20.4 cannot read gzip — feeding it the
		// bgzipped R1 yields empty output on both sides.
		//
		// The stats cell uses the DETERMINISTIC subset of -stats_* groups
		// (everything except -stats_tag). Upstream's stats_tag "midseq" value
		// is a `join(',', keys %MIDS)` over an unsorted Perl hash
		// (prinseq-lite.pl:1970), so its token order varies run-to-run whenever
		// two or more predefined 454 MIDs each clear the 3% threshold — a
		// nondeterministic byte-exact comparison. Excluding -stats_tag matches
		// what the tool's own upstream-parity unit test drops.
		{Tool: "prinseq", Name: "prinseq_stats", Subcommand: "prinseq -stats", Need: NeedFastqPlain, Post: PostStdout,
			OurArgs: []string{"-fastq", phFastqPlain,
				"-stats_info", "-stats_len", "-stats_dinuc", "-stats_dupl", "-stats_ns", "-stats_assembly"}},
		{Tool: "prinseq", Name: "prinseq_filter", Subcommand: "prinseq filter", Need: NeedFastqPlain,
			Post: PostFile, WorkDirOut: true, Compare: "rb_good.fastq",
			OurArgs: []string{"-fastq", phFastqPlain, "-min_len", "30", "-out_good", phOutdir + "/rb_good",
				"-out_bad", "null"}},
	}
}

// vcftoolsCells covers a representative vcftools flag battery (VCF.gz). vcftools
// writes <prefix>.<ext> files into cwd, so a work dir is used and the named
// output file is compared.
func vcftoolsCells() []CellSpec {
	const t = "vcftools"
	mk := func(name, sub, outFile string, extra ...string) CellSpec {
		args := append([]string{"--gzvcf", phVCF, "--out", phOutdir + "/rb"}, extra...)
		return CellSpec{Tool: t, Name: name, Subcommand: sub, Need: NeedVCF, Post: PostFile,
			WorkDirOut: true, Compare: outFile, OurArgs: args}
	}
	return []CellSpec{
		mk("vcftools_freq", "--freq", "rb.frq", "--freq"),
		mk("vcftools_counts", "--counts", "rb.frq.count", "--counts"),
		mk("vcftools_depth", "--depth", "rb.idepth", "--depth"),
		mk("vcftools_site_mean_depth", "--site-mean-depth", "rb.ldepth.mean", "--site-mean-depth"),
		mk("vcftools_missing_indv", "--missing-indv", "rb.imiss", "--missing-indv"),
		mk("vcftools_het", "--het", "rb.het", "--het"),
		mk("vcftools_hardy", "--hardy", "rb.hwe", "--hardy"),
		mk("vcftools_relatedness2", "--relatedness2", "rb.relatedness2", "--relatedness2"),
		mk("vcftools_snpdensity", "--SNPdensity", "rb.snpden", "--SNPdensity", "100000"),
		mk("vcftools_window_pi", "--window-pi", "rb.windowed.pi", "--window-pi", "100000"),
	}
}

// mosdepthCells covers mosdepth (BAM/CRAM + ref). mosdepth writes
// <prefix>.mosdepth.* files into cwd; the summary file is compared. Upstream
// mosdepth is linux/amd64-only and resolved via MOSDEPTH_BIN; absent it, these
// degrade to ours-only perf cells.
func mosdepthCells() []CellSpec {
	const t = "mosdepth"
	return []CellSpec{
		{Tool: t, Name: "mosdepth_default", Subcommand: "mosdepth", Need: NeedBAM, Post: PostFile,
			WorkDirOut: true, Compare: "rb.mosdepth.summary.txt",
			OurArgs: []string{"rb", phBAM}},
		{Tool: t, Name: "mosdepth_by_bed", Subcommand: "mosdepth --by", Need: NeedBAM | NeedBED, Post: PostFile,
			WorkDirOut: true, Compare: "rb.mosdepth.summary.txt",
			OurArgs: []string{"--by", phBED, "rb", phBAM}},
		{Tool: t, Name: "mosdepth_fast", Subcommand: "mosdepth --fast-mode", Need: NeedBAM, Post: PostFile,
			WorkDirOut: true, Compare: "rb.mosdepth.summary.txt",
			OurArgs: []string{"--fast-mode", "rb", phBAM}},
	}
}

// htslibCells covers bgzip, tabix, and htsfile.
func htslibCells() []CellSpec {
	return []CellSpec{
		// bgzip: compress a BED (stdout gzip stream), then -d round-trip, then --test.
		{Tool: "bgzip", Name: "bgzip_compress", Subcommand: "bgzip", Need: NeedBED, Post: PostStdoutGzip,
			OurArgs: []string{"-c", phBED}},
		{Tool: "bgzip", Name: "bgzip_decompress", Subcommand: "bgzip -d", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{"-d", "-c", phVCF}},
		{Tool: "bgzip", Name: "bgzip_test", Subcommand: "bgzip --test", Need: NeedVCF, Post: PostNone,
			OurArgs: []string{"--test", phVCF}},

		// tabix: region query over the indexed VCF.
		{Tool: "tabix", Name: "tabix_query", Subcommand: "tabix region", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{phVCF, "1"}},

		// htsfile: file type identification.
		{Tool: "htsfile", Name: "htsfile_bam", Subcommand: "htsfile", Need: NeedBAM, Post: PostStdout,
			OurArgs: []string{phBAM}},
		{Tool: "htsfile", Name: "htsfile_vcf", Subcommand: "htsfile", Need: NeedVCF, Post: PostStdout,
			OurArgs: []string{phVCF}},
	}
}
