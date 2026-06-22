package main

// This file defines the real-world differential cell battery: the
// representative samtools/bcftools commands run, on whole-genome (multi-contig)
// input by default, against BOTH our port and the upstream binary, comparing
// provenance-stripped stdout and timing each side.
//
// A cell SKIPs (rather than fails) when its required input is absent, so the
// command runs partially with whatever data is present. Cross-contig behaviour
// (BAI/CSI multi-ref bins, RNEXT '=' vs mate-on-other-contig, coordinate sort
// across contigs, per-contig idxstats counts) is exactly what the multi-contig
// input exercises, so `view` of the whole file and `idxstats` are in the
// battery on purpose.

// inputKind names the required input for a cell. A cell whose required input is
// absent SKIPs.
type inputKind int

const (
	needBAM    inputKind = iota // requires the multi-contig BAM
	needBAMRef                  // requires the BAM and the reference FASTA (CRAM paths)
	needVCF                     // requires the bgzipped+indexed VCF
	needVCFRef                  // requires the VCF and the reference FASTA (norm -f)
)

// postKind names how a cell's primary output is turned into a comparable text
// stream. Streaming text cells (SAM/flagstat/stats/...) compare stdout
// directly; file-producing cells (BAM/CRAM/index) re-decode the written file
// through a `view` step so two outputs with different BGZF/CRAM framing are
// compared by their decoded records (the repo-documented caveat).
type postKind int

const (
	postStdout  postKind = iota // primary output is stdout text; compare directly
	postViewSAM                 // primary output is a written BAM/CRAM; `view -h` it to SAM
	postNone                    // command's success is the signal (quickcheck); compare exit-status text
)

// cellSpec is the declarative description of one battery cell. The same spec is
// instantiated against our binary and the upstream binary; argv is built by
// substituting the resolved input paths into args via the {bam}/{vcf}/{ref}
// placeholders. tool is "samtools" or "bcftools".
type cellSpec struct {
	Name  string
	Tool  string
	Need  inputKind
	Post  postKind
	Args  []string // argv after the subcommand; placeholders substituted at run time
	Multi bool     // true => exercises explicit cross-contig behaviour (for the report note)
	// WriteOut, when set, names the {out} placeholder file extension the command
	// writes to (e.g. ".bam", ".cram"). The runner allocates a temp path,
	// substitutes it for {out}, then (postViewSAM) re-decodes it.
	WriteOut string
}

// battery is the full ordered cell list. It is whole-genome / multi-contig by
// default: no cell pins a single chromosome unless the operator passes -region,
// which the runner appends where the command accepts a region argument.
func battery() []cellSpec {
	return []cellSpec{
		// ---- samtools: streaming text ----
		{Name: "samtools_view_sam", Tool: "samtools", Need: needBAM, Post: postStdout,
			Args: []string{"view", "{bam}"}, Multi: true},
		{Name: "samtools_view_sam_header", Tool: "samtools", Need: needBAM, Post: postStdout,
			Args: []string{"view", "-h", "{bam}"}, Multi: true},
		{Name: "samtools_flagstat", Tool: "samtools", Need: needBAM, Post: postStdout,
			Args: []string{"flagstat", "{bam}"}},
		{Name: "samtools_idxstats", Tool: "samtools", Need: needBAM, Post: postStdout,
			Args: []string{"idxstats", "{bam}"}, Multi: true},
		{Name: "samtools_stats", Tool: "samtools", Need: needBAM, Post: postStdout,
			Args: []string{"stats", "{bam}"}},
		{Name: "samtools_depth_a", Tool: "samtools", Need: needBAM, Post: postStdout,
			Args: []string{"depth", "-a", "{bam}"}},
		{Name: "samtools_quickcheck", Tool: "samtools", Need: needBAM, Post: postNone,
			Args: []string{"quickcheck", "{bam}"}},

		// ---- samtools: file-producing (re-decoded) ----
		{Name: "samtools_view_bam", Tool: "samtools", Need: needBAM, Post: postViewSAM,
			Args: []string{"view", "-b", "-o", "{out}", "{bam}"}, WriteOut: ".bam", Multi: true},
		{Name: "samtools_sort", Tool: "samtools", Need: needBAM, Post: postViewSAM,
			Args: []string{"sort", "-o", "{out}", "{bam}"}, WriteOut: ".bam", Multi: true},
		{Name: "samtools_view_cram", Tool: "samtools", Need: needBAMRef, Post: postViewSAM,
			Args: []string{"view", "-C", "-T", "{ref}", "-o", "{out}", "{bam}"}, WriteOut: ".cram", Multi: true},

		// ---- bcftools ----
		{Name: "bcftools_view", Tool: "bcftools", Need: needVCF, Post: postStdout,
			Args: []string{"view", "{vcf}"}, Multi: true},
		{Name: "bcftools_view_body", Tool: "bcftools", Need: needVCF, Post: postStdout,
			Args: []string{"view", "-H", "{vcf}"}, Multi: true},
		{Name: "bcftools_norm", Tool: "bcftools", Need: needVCFRef, Post: postStdout,
			Args: []string{"norm", "-f", "{ref}", "-O", "v", "{vcf}"}},
		{Name: "bcftools_stats", Tool: "bcftools", Need: needVCF, Post: postStdout,
			Args: []string{"stats", "{vcf}"}},
		{Name: "bcftools_query", Tool: "bcftools", Need: needVCF, Post: postStdout,
			Args: []string{"query", "-f", `%CHROM\t%POS\t%REF\t%ALT\n`, "{vcf}"}, Multi: true},
	}
}

// cramEnv is the child environment used for CRAM operations so a missing MD5 in
// the local cache never triggers a network fetch from the EBI reference server.
// Pointing REF_PATH at a non-existent cache plus passing -T <ref> on the command
// line makes CRAM fully local and offline.
func cramEnv(base []string) []string {
	out := make([]string, 0, len(base)+2)
	for _, e := range base {
		if len(e) >= 9 && e[:9] == "REF_PATH=" {
			continue
		}
		if len(e) >= 10 && e[:10] == "REF_CACHE=" {
			continue
		}
		out = append(out, e)
	}
	return append(out, "REF_PATH=/nonexistent/realparity/cram/cache/%2s/%2s/%s", "REF_CACHE=")
}
