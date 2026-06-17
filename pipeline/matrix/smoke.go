package matrix

// This file registers the SMOKE matrix: a small, representative set of entries
// spanning the major formats (BAM, CRAM, VCF, BED) that proves the ours-vs-
// upstream loop end to end. Later agents add the full per-tool matrices in
// their own files (e.g. samtools.go, bcftools.go) calling Register from init,
// or by appending to the spec slices here.
//
// Each block demonstrates the curated combinatorics expander (ExpandSpec) plus
// a couple of hand-written entries for cases the expander does not cover.

func init() {
	Register(smokeSamtoolsView()...)
	Register(smokeBcftools()...)
	Register(smokeBedIntersect()...)
}

// smokeSamtoolsView exercises `samtools view` over BAM and CRAM. The BAM cases
// use the expander (baseline + single flags + a curated combo); the CRAM cases
// are written out because they need the -T reference argument.
func smokeSamtoolsView() []Entry {
	bam := ExpandSpec{
		Tool:           "samtools",
		Subcommand:     "view",
		UsesSubcommand: true,
		Input:          InputBAM,
		Compare:        ByteExact,
		BaseArgs:       []string{"{bam}"},
		Flags: []Flag{
			{Name: "-H", Bool: true},               // header only
			{Name: "-c", Bool: true},               // count
			{Name: "-q", Values: []string{"30"}},   // MAPQ filter
			{Name: "-f", Values: []string{"0x10"}}, // require reverse-strand flag
			{Name: "-F", Values: []string{"0x10"}}, // exclude reverse-strand flag
		},
		Combos: []Combo{
			{Name: "count_q30", Flags: []string{"-c", "-q", "30"}},
			{Name: "header_and_reads", Flags: []string{"-h"}}, // -h: header + reads (lowercase)
		},
	}.Expand()
	for i := range bam {
		bam[i].Name = "view_bam_" + bam[i].Name
	}

	cram := []Entry{
		{
			Tool: "samtools", Subcommand: "view", UsesSubcommand: true,
			Name: "view_cram_body", Input: InputCRAM, Compare: ByteExact,
			Args: []string{"-T", "{fasta}", "{cram}"},
		},
		{
			Tool: "samtools", Subcommand: "view", UsesSubcommand: true,
			Name: "view_cram_count", Input: InputCRAM, Compare: ByteExact,
			Args: []string{"-c", "-T", "{fasta}", "{cram}"},
		},
		{
			Tool: "samtools", Subcommand: "view", UsesSubcommand: true,
			Name: "view_cram_decode_sam_heavy", Input: InputCRAM, Compare: ByteExact, Heavy: true,
			Args: []string{"-T", "{fasta}", "{cram}"},
			// Full CRAM->SAM decode of every record; the SAM text is byte-exact
			// (no provenance in the body). Marked heavy so the report surfaces
			// the decode timing ratio. (Binary BGZF BAM output is deliberately
			// NOT compared byte-exact here: our klauspost deflate backend frames
			// blocks differently from htslib though both decode identically;
			// that path is covered by the bench harness, not byte parity.)
		},
	}
	return append(bam, cram...)
}

// smokeBcftools exercises `bcftools view` and `bcftools query` over the plain
// VCF. view uses -H (body only) so VCF header line-ordering differences (which
// upstream normalises) do not create spurious divergence; query is body-only by
// construction.
func smokeBcftools() []Entry {
	view := ExpandSpec{
		Tool:           "bcftools",
		Subcommand:     "view",
		UsesSubcommand: true,
		Input:          InputVCFPlain,
		Compare:        ByteExact,
		BaseArgs:       []string{"-H", "{vcf_plain}"},
		Flags: []Flag{
			{Name: "-v", Values: []string{"snps", "indels"}}, // variant-type select
			{Name: "-i", Values: []string{"QUAL>30"}},        // include expression
			{Name: "-e", Values: []string{"QUAL<30"}},        // exclude expression
		},
		Combos: []Combo{
			{Name: "snps_qual", Flags: []string{"-v", "snps", "-i", "QUAL>30"}},
		},
	}.Expand()
	for i := range view {
		view[i].Name = "view_" + view[i].Name
	}

	query := []Entry{
		{
			Tool: "bcftools", Subcommand: "query", UsesSubcommand: true,
			Name: "query_chrom_pos_ref_alt", Input: InputVCFPlain, Compare: ByteExact,
			Args: []string{"-f", `%CHROM\t%POS\t%REF\t%ALT\n`, "{vcf_plain}"},
		},
		{
			Tool: "bcftools", Subcommand: "query", UsesSubcommand: true,
			Name: "query_gt", Input: InputVCFPlain, Compare: ByteExact,
			Args: []string{"-f", `%CHROM:%POS[\t%GT]\n`, "{vcf_plain}"},
		},
		{
			Tool: "bcftools", Subcommand: "query", UsesSubcommand: true,
			Name: "query_info_dp", Input: InputVCFPlain, Compare: ByteExact,
			Args: []string{"-f", `%CHROM\t%POS\t%INFO/DP\n`, "{vcf_plain}"},
		},
	}
	return append(view, query...)
}

// smokeBedIntersect exercises `bedtools intersect`. Our binary IS the
// subcommand (bedintersect), so UsesSubcommand is false: only the upstream
// invocation receives the "intersect" token, and UpstreamTool maps to the
// bedtools binary. The two inputs are the plain BED (-a) and the BED12 (-b).
func smokeBedIntersect() []Entry {
	spec := ExpandSpec{
		Tool:           "bedintersect",
		Subcommand:     "intersect",
		UpstreamTool:   "bedtools",
		UsesSubcommand: false,
		Input:          InputBED,
		Compare:        ByteExact,
		BaseArgs:       []string{"-a", "{bed}", "-b", "{bed12}"},
		Flags: []Flag{
			{Name: "-c", Bool: true}, // count overlaps per A
			{Name: "-v", Bool: true}, // A with no overlap
			{Name: "-u", Bool: true}, // unique A reported once
			{Name: "-wa", Bool: true},
			{Name: "-wb", Bool: true},
			{Name: "-s", Bool: true}, // same-strand
		},
		Combos: []Combo{
			{Name: "wa_wb", Flags: []string{"-wa", "-wb"}},
			{Name: "u_s", Flags: []string{"-u", "-s"}},
			{Name: "c_s", Flags: []string{"-c", "-s"}},
		},
	}
	entries := spec.Expand()
	// Mark the full write-A-write-B variant heavy to demonstrate timing on the
	// largest output of the set.
	for i := range entries {
		entries[i].Name = "intersect_" + entries[i].Name
		if entries[i].Name == "intersect_combo_wa_wb" {
			entries[i].Heavy = true
		}
	}
	return entries
}
