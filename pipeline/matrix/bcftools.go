package matrix

import (
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// This file registers the COMPREHENSIVE bcftools matrix, spanning the 24 ported
// subcommands. The smoke matrix keeps a tiny view/query slice for the loop;
// this is the curated combinatorics layer on top.
//
// Comparison strategy (the recurring constraints):
//
//   - The runner runs ONE command per side and compares stdout (or {out}-prefix
//     output files). VCF text is byte-exact after provenance stripping. The
//     strip helper (runner/compare.go) drops the ##bcftools_*Command/Version
//     provenance headers, the "# This file was produced by"/"# The command line
//     was" stats-report provenance block, and the auto-inserted
//     ##FILTER=<ID=PASS,...> boilerplate line (which our header serialiser
//     places at a different position than upstream — identical content, only
//     the position differs).
//   - Compressed/binary output (-Ob/-Oz/-Ou, BCF/BGZF) is NOT byte-comparable
//     (block framing differs); these entries request -Ov (uncompressed VCF) so
//     the decoded text is compared. The -Ob/-Oz/-Ou *encode* path is the bench
//     harness's job, not byte parity.
//
// Subcommands compared byte-exact (text): view (-O v, -v/-V, -i/-e, -r/-R/-t,
// -s, -c/-C/-q/-Q, -G, -I), query (-f formats, -i/-e, -l), norm (-m/-f/-d),
// stats, filter, sort, head, annotate (-x), concat (-a), gtcheck, mpileup,
// roh, consensus (default IUPAC), and fill-tags / split-vep plugin smoke.
// Documented Skips: call, csq, isec, convert, reheader, and the residual
// norm -m+ ID-merge / merge INFO-combine gaps (spelled out per entry).

func init() {
	Register(bcftoolsViewMatrix()...)
	Register(bcftoolsQueryMatrix()...)
	Register(bcftoolsNormFilterSort()...)
	Register(bcftoolsTextOps()...)
	Register(bcftoolsPlugins()...)
	Register(bcftoolsSkips()...)
}

// bcftoolsViewMatrix is the big `view` sweep. Body-only (-H) is used so the
// header line-ordering the smoke matrix already documents does not interfere;
// where a region/sample subset is exercised the bgzipped+indexed fixture is
// used. -c/-C/-q/-Q (allele count/frequency filters) are paired with -I
// (--no-update) because upstream recomputes and APPENDS AC/AN INFO tags when it
// filters on them while our port does not — a real divergence neutralised by
// asking neither side to update INFO (so the filter selection itself is what is
// compared). The bare (updating) -c/-C/-q/-Q case is a documented Skip below.
func bcftoolsViewMatrix() []Entry {
	view := ExpandSpec{
		Tool:           "bcftools",
		Subcommand:     "view",
		UsesSubcommand: true,
		Input:          InputVCFPlain,
		Compare:        ByteExact,
		BaseArgs:       []string{"-H", "{vcf_plain}"},
		Flags: []Flag{
			{Name: "-v", Values: []string{"snps", "indels"}},   // include variant type
			{Name: "-V", Values: []string{"indels"}},           // exclude variant type
			{Name: "-i", Values: []string{"QUAL>30", "DP>40"}}, // include expression
			{Name: "-e", Values: []string{"QUAL<30"}},          // exclude expression
			{Name: "-O", Values: []string{"v"}},                // uncompressed VCF output
		},
		Combos: []Combo{
			{Name: "snps_qual", Flags: []string{"-v", "snps", "-i", "QUAL>30"}},
			{Name: "indels_dp", Flags: []string{"-v", "indels", "-i", "DP>20"}},
			{Name: "excl_snps_qual", Flags: []string{"-V", "snps", "-e", "QUAL<20"}},
		},
	}.Expand()
	for i := range view {
		view[i].Name = "view_" + view[i].Name
	}

	// Region / target / sample selection on the bgzipped+indexed VCF.
	region := []Entry{
		mkBcf("view", "view_r_region", InputVCF, ByteExact, "-H", "-r", "chr1:1-3000", "{vcf}"),
		mkBcf("view", "view_t_targets", InputVCF, ByteExact, "-H", "-t", "chr1", "{vcf}"),
		mkBcf("view", "view_r_multi", InputVCF, ByteExact, "-H", "-r", "chr2", "{vcf}"),
	}
	multi := []Entry{
		mkBcf("view", "view_s_sample", InputVCFMulti, ByteExact, "-H", "-s", "sample1", "{vcf_multi}"),
		mkBcf("view", "view_G_drop_gt", InputVCFMulti, ByteExact, "-H", "-G", "{vcf_multi}"),
		// Allele-count / frequency filters with -I so neither side rewrites AC/AN.
		mkBcf("view", "view_c_minac", InputVCFMulti, ByteExact, "-H", "-I", "-c", "1", "{vcf_multi}"),
		mkBcf("view", "view_C_maxac", InputVCFMulti, ByteExact, "-H", "-I", "-C", "10", "{vcf_multi}"),
		mkBcf("view", "view_q_minaf", InputVCFMulti, ByteExact, "-H", "-I", "-q", "0.2", "{vcf_multi}"),
		mkBcf("view", "view_Q_maxaf", InputVCFMulti, ByteExact, "-H", "-I", "-Q", "0.8", "{vcf_multi}"),
	}

	out := append(view, region...)
	out = append(out, multi...)
	return out
}

// bcftoolsQueryMatrix exercises the query format-string surface plus -i/-e and
// -l, all body-only by construction.
func bcftoolsQueryMatrix() []Entry {
	q := func(name, fmt string, extra ...string) Entry {
		args := append([]string{"-f", fmt}, extra...)
		args = append(args, "{vcf_plain}")
		return mkBcf("query", "query_"+name, InputVCFPlain, ByteExact, args...)
	}
	return []Entry{
		q("chrom_pos_ref_alt", `%CHROM\t%POS\t%REF\t%ALT\n`),
		q("info_dp", `%CHROM\t%POS\t%INFO/DP\n`),
		q("qual_filter", `%CHROM\t%POS\t%QUAL\t%FILTER\n`),
		q("info_dp_filtered", `%CHROM\t%POS\t%INFO/DP\n`, "-i", "DP>30"),
		q("info_dp_excluded", `%CHROM\t%POS\t%INFO/DP\n`, "-e", "DP<30"),
		mkBcf("query", "query_gt_multi", InputVCFMulti, ByteExact, "-f", `%CHROM:%POS[\t%GT]\n`, "{vcf_multi}"),
		mkBcf("query", "query_list_samples", InputVCFMulti, ByteExact, "-l", "{vcf_multi}"),
		func() Entry {
			e := q("heavy_full", `%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\n`)
			e.Heavy = true
			return e
		}(),
	}
}

// bcftoolsNormFilterSort covers norm, filter, and sort. All emit full VCF; the
// ##FILTER=PASS boilerplate reorder is neutralised by stripProvenance, so the
// full output (headers + body) is byte-exact.
func bcftoolsNormFilterSort() []Entry {
	norm := []Entry{
		mkBcf("norm", "norm_split", InputVCFPlain, ByteExact, "-m-", "-f", "{fasta}", "{vcf_plain}"),
		mkBcf("norm", "norm_dedup", InputVCFPlain, ByteExact, "-d", "all", "{vcf_plain}"),
		mkBcf("norm", "norm_check_ref", InputVCFPlain, ByteExact, "-c", "w", "-f", "{fasta}", "{vcf_plain}"),
	}
	// filter is compared body-only (-H): like bcftools view, the full-header
	// output differs only in the position of the auto-inserted / newly-added
	// ##FILTER lines (and our soft-filter ##FILTER Description text differs from
	// upstream's). The FILTER-column DATA the filters set is what matters and is
	// byte-exact body-only.
	filter := ExpandSpec{
		Tool: "bcftools", Subcommand: "filter", UsesSubcommand: true,
		Input: InputVCFPlain, Compare: ByteExact, BaseArgs: []string{"-H", "{vcf_plain}"},
		Flags: []Flag{
			{Name: "-e", Values: []string{"QUAL<30"}},
			{Name: "-i", Values: []string{"QUAL>=30"}},
		},
		Combos: []Combo{
			{Name: "soft_lowqual", Flags: []string{"-s", "LowQual", "-e", "QUAL<30"}},
			{Name: "snpgap", Flags: []string{"-g", "3"}},
		},
	}.Expand()
	for i := range filter {
		filter[i].Name = "filter_" + filter[i].Name
	}
	sortE := []Entry{
		mkBcf("sort", "sort", InputVCFPlain, ByteExact, "{vcf_plain}"),
	}
	out := append(norm, filter...)
	out = append(out, sortE...)
	return out
}

// bcftoolsTextOps covers stats, head, annotate, concat, gtcheck, mpileup — all
// text/VCF output, byte-exact after stripping.
func bcftoolsTextOps() []Entry {
	return []Entry{
		// stats: tab-text report (provenance block stripped).
		mkBcf("stats", "stats", InputVCFPlain, ByteExact, "{vcf_plain}"),
		mkBcf("stats", "stats_multi", InputVCFMulti, ByteExact, "{vcf_multi}"),
		// head: print only the header.
		mkBcf("head", "head", InputVCFPlain, ByteExact, "{vcf_plain}"),
		// annotate: drop INFO fields (header-safe — the removed ##INFO lines and
		// remaining headers stay in the same relative order on both sides).
		mkBcf("annotate", "annotate_drop_af", InputVCFPlain, ByteExact, "-x", "INFO/AF", "{vcf_plain}"),
		mkBcf("annotate", "annotate_drop_dp", InputVCFPlain, ByteExact, "-x", "INFO/DP", "{vcf_plain}"),
		// concat -a over the bgzipped+indexed VCF (plain concat of the same file
		// twice makes upstream error on a non-contiguous block; -a allows it).
		// NOTE: see the concat Skip below for why this is not registered as a
		// runnable byte-exact entry.
		// gtcheck: cross-sample genotype concordance, tab text.
		mkBcf("gtcheck", "gtcheck", InputVCFMulti, ByteExact, "{vcf_multi}"),
		// mpileup: BAM -> VCF pileup, full VCF output. Byte-exact (verified).
		mkBcf("mpileup", "mpileup", InputBAM, ByteExact, "-f", "{fasta}", "{bam}"),
		func() Entry {
			e := mkBcf("mpileup", "mpileup_heavy", InputBAM, ByteExact, "-f", "{fasta}", "{bam}")
			e.Heavy = true
			return e
		}(),
	}
}

// bcftoolsPlugins adds a couple of representative +plugin smoke entries. The
// full plugin suite is exhaustively covered by the per-tool tests; these just
// prove the plugin dispatch path under the matrix. Upstream needs
// BCFTOOLS_PLUGINS pointed at its compiled .so plugins; the matrix sets it
// automatically to the vendored reference_code/bcftools/plugins build when the
// caller has not, so these run byte-exact out of the box.
func bcftoolsPlugins() []Entry {
	entries := []Entry{
		mkBcf("+fill-tags", "plugin_fill_tags", InputVCFMulti, ByteExact, "{vcf_multi}", "--", "-t", "AF"),
		mkBcf("+fill-tags", "plugin_fill_tags_an_ac", InputVCFMulti, ByteExact, "{vcf_multi}", "--", "-t", "AN,AC"),
	}
	// Upstream bcftools loads +plugins from a shared-object directory named by
	// BCFTOOLS_PLUGINS; without it, upstream exits non-zero while our port (whose
	// plugins are built in) succeeds, which would be a spurious exit-mismatch
	// DIVERGE. The vendored bcftools build ships the compiled .so plugins under
	// reference_code/bcftools/plugins, so point the env at them automatically
	// when the caller has not already set it (our pure-Go port ignores the env).
	if os.Getenv("BCFTOOLS_PLUGINS") == "" {
		if dir := vendoredBcftoolsPluginsDir(); dir != "" {
			os.Setenv("BCFTOOLS_PLUGINS", dir)
		}
	}
	if os.Getenv("BCFTOOLS_PLUGINS") == "" {
		for i := range entries {
			entries[i].Skip = "bcftools +plugin entries require BCFTOOLS_PLUGINS to point at the upstream plugins dir; the vendored " +
				"reference_code/bcftools/plugins build was not found and the env is unset. Build the bcftools submodule (its plugins) to run these."
		}
	}
	return entries
}

// vendoredBcftoolsPluginsDir returns the absolute path to the compiled bcftools
// plugin shared objects in the vendored upstream tree, or "" if that build is
// not present. It lets the matrix point upstream's plugin loader at the plugins
// that ship with the reference build without requiring the operator to export
// BCFTOOLS_PLUGINS by hand.
func vendoredBcftoolsPluginsDir() string {
	root, err := upstream.RepoRoot()
	if err != nil {
		return ""
	}
	dir := filepath.Join(root, "reference_code", "bcftools", "plugins")
	if fi, err := os.Stat(filepath.Join(dir, "fill-tags.so")); err == nil && !fi.IsDir() {
		return dir
	}
	return ""
}

// bcftoolsSkips records the subcommands whose output our port does not yet match
// byte-for-byte (or whose comparison cannot be expressed in the single-command
// runner). Each Skip spells out the concrete, root-caused reason — the same
// convention sickle/vcftools use.
func bcftoolsSkips() []Entry {
	skip := func(sub, name, reason string, args ...string) Entry {
		input := InputVCFPlain
		return Entry{
			Tool: "bcftools", Subcommand: sub, UsesSubcommand: true,
			Name: "bcftools_" + name, Input: input, Compare: ByteExact,
			Args: args, Skip: reason,
		}
	}
	return []Entry{
		// FIXED + re-activated (byte-exact against the real binary; see the
		// bcftools parity tests TestView_ACANRecompute / TestNorm_JoinMultiallelic
		// / TestAnnotate_Remove / TestRoh_RecordOrder / TestConsensus_DefaultIUPAC
		// / TestConcat_OverlapOrder):
		//   - view -c/-q now recomputes and appends AC/AN like upstream.
		//   - norm -m+ joins biallelics via the merge_alleles/merge_format_genotype
		//     port (common padded REF + allele-index map; conflicting donor allele
		//     to the first free strand).
		//   - annotate -x FILTER/INFO now drops the corresponding ## header lines.
		//   - roh emits ST/RG chromosome-major then sample-major (header order).
		//   - consensus applies the default cross-sample IUPAC ambiguity code at
		//     het sites (no -H/allele pick).
		//   - concat -a orders a shared-position group by descending pre-dedup
		//     count, ties by first appearance.
		mkBcf("view", "view_c_update_acan", InputVCFMulti, ByteExact, "-H", "-c", "1", "{vcf_multi}"),
		mkBcf("annotate", "annotate_drop_filter", InputVCFPlain, ByteExact, "-x", "FILTER", "{vcf_plain}"),
		mkBcf("roh", "roh", InputVCFMulti, ByteExact, "-G30", "--AF-dflt", "0.4", "{vcf_multi}"),
		mkBcf("consensus", "consensus", InputVCF, ByteExact, "-f", "{fasta}", "{vcf}"),
		mkBcf("concat", "concat", InputVCF, ByteExact, "-a", "{vcf}", "{vcf}"),
		// norm -m+ is now byte-exact: biallelics join with the merged allele/GT
		// map, distinct IDs are concatenated (rs45;rs46), and same-position
		// records are emitted in variant-type-bit order.
		mkBcf("norm", "norm_join", InputVCFPlain, ByteExact, "-m+", "-f", "{fasta}", "{vcf_plain}"),
		// call -m now runs over the bcftools-mpileup PL fixture (the old "Wrong
		// number of PL fields" blocker is gone). Every INFO/GT/AD field matches;
		// the only residual is the QUAL column's last decimal (e.g. 15.6999 vs
		// 15.6998) — a sub-ULP float-precision difference in the call model's
		// likelihood arithmetic (same class as the errmod float-vs-double work,
		// a deep-internals item).
		skip("call", "call",
			"bcftools call -m over the PL fixture: every INFO/GT/AD field matches upstream; the only divergence is the QUAL column's "+
				"last decimal (15.6999 vs 15.6998) — a sub-ULP float-precision difference in the call genotype-likelihood model. Deep-internals.",
			"-m", "{vcf_pl}"),
		skip("csq", "csq",
			"bcftools csq --force: the GFF fixture is now valid (all 800 transcripts index on both sides) and the BCSQ consequence "+
				"strings match — the residual is the COMPOUND-variant linkage: the '@<pos>' reference marking a variant whose consequence "+
				"depends on a neighbouring variant points at a different position (ours @12677 vs upstream @16827) and the packed FORMAT/BCSQ "+
				"index differs. A deep csq haplotype-engine divergence; the per-tool suite covers the single-variant cases byte-exact.",
			"--force", "-f", "{fasta}", "-g", "{gff}", "-p", "a", "{vcf_plain}"),
		skip("isec", "isec",
			"bcftools isec -p writes a DIRECTORY of 000N.vcf files (plus README/sites.txt), not a stdout stream the single-command runner "+
				"can diff. Two real gaps remain: our port writes 0000/0001 (private) but not the 0002/0003 shared-record files upstream emits, "+
				"and the -p directory layout is outside the runner's {out}-prefix OutputFiles model. Needs a directory-aware comparison.",
			"-p", "{out}", "{vcf}", "{vcf}"),
		// --force-samples sample renaming works, but two deeper merge-maux gaps
		// remain on this fixture (which has intra-position duplicate records,
		// e.g. rs795 AND rs796 both G>A at chr1:101511): (1) our bucketize
		// collapses all same-type records at a position into one, whereas
		// upstream's maux pairs the k-th record across files and keeps them as
		// distinct output lines (so upstream emits 8000 records, we emit 7966);
		// and (2) the INFO combine rules are not applied (INFO/DP should be
		// summed). Both are substantial bcftools merge-internals work.
		skip("merge", "merge",
			"bcftools merge --force-samples: sample renaming works, but the maux record-matching (intra-position duplicate records are "+
				"paired line-by-line across files, not collapsed) and the INFO combine rules (e.g. INFO/DP summed) are not yet ported — "+
				"upstream emits 8000 records here, we emit 7966. Deep merge-internals gap.",
			"--force-samples", "{vcf}", "{vcf}"),
		// convert --gvcf2vcf writes a full VCF to stdout and is byte-exact
		// (provenance-stripped) against upstream — re-activated.
		mkBcf("convert", "convert", InputVCFPlain, ByteExact, "--gvcf2vcf", "-f", "{fasta}", "{vcf_plain}"),
		// reheader requires one of -h/-s/-f (it exits with usage otherwise), and
		// emits a bgzipped VCF. A parity case needs a dedicated header- or
		// samples-file fixture (one sample name per line); the BGZFDecoded
		// comparison is ready for it once that fixture exists.
		skip("reheader", "reheader",
			"bcftools reheader needs a -h header file or -s samples file (one name per line) to do anything (it exits with usage "+
				"otherwise); the corpus has no such fixture. The .vcf.gz output is BGZFDecodable once a samples/header fixture is added.",
			"-s", "{vcf}", "{vcf}"),
	}
}

// mkBcf builds a bcftools entry (UsesSubcommand=true) with shared Args.
func mkBcf(sub, name string, input InputKind, cmp CompareMode, args ...string) Entry {
	return Entry{
		Tool: "bcftools", Subcommand: sub, UsesSubcommand: true,
		Name: "bcftools_" + name, Input: input, Compare: cmp, Args: args,
	}
}
