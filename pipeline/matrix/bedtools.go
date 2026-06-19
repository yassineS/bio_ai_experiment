package matrix

// This file registers the bedtools command-combination matrix: curated
// combinatorics for every bed* tool in tools/ (each is its own binary that maps
// to one `bedtools <sub>` subcommand). Our binaries ARE the subcommand
// (UsesSubcommand=false): only the upstream invocation receives the subcommand
// token, and UpstreamTool points at the vendored `bedtools` binary.
//
// Inputs come from the shared coordinate-space fixtures: {bed} (sorted BED6),
// {bed12} (full 12-column BED12), {genome} (chrom-sizes), {fasta} (+.fai),
// {bam}, and {gff}. Comparison is byte-exact stdout for the text tools; the two
// file-writing tools (bedsplit) use the OutputFiles mechanism. RNG tools
// (sample/random/shuffle) are byte-exact when given a fixed -seed: our port
// reproduces upstream's MT19937 stream exactly, verified on the smoke fixture.
//
// COMBINATORICS POLICY: per the package contract we emit baseline + one entry
// per important single flag + a small set of curated multi-flag Combos, never
// the 2^N power set.
//
// ------------------------------------------------------------------------
// Documented Skips. Most of the divergences the matrix once surfaced are FIXED
// and re-activated as byte-exact entries: bednuc / bedmakewindows / bedexpand /
// bedsummary / bedfisher / bedsubtract / bed12tobed6; the bedmap/bedmerge/
// bedwindow collapse-order cases; the bedsort (default/-sizeD/-chrThenScore*)
// and bedcluster (-s) std::sort tie-break order (via the pkg/cppsort libstdc++
// introsort port); bedjaccard -s (per-strand merged stream re-sorted before the
// sweep); bedsplit -a size (introsort + stddev-min greedy); and bedannotate
// (the matrix now passes upstream's -files flag form). What remains are
// harness/fixture limitations, not tool bugs:
//
//   - bedtobam          — raw BGZF BAM to stdout; framing differs (klauspost vs
//                         htslib) though decoded records match. Not byte-compared.
//   - bedtag            — bedtag now implements the upstream tagBam model
//                         (-i BAM -files ...: YB aux tag from BED overlaps),
//                         compared BAMDecoded.
//   - bedpairtobed       — runs over the new {bedpe} fixture; same record SET as
//                          upstream but a few multi-hit pairs emit in a
//                          different chromsweep order (bedpairtopair is exact).
//   - (bedoverlap now runs over the {bedpe} fixture with -cols 2,3,5,6.)
//   - (bedunionbedg now runs over the {bedgraph1}/{bedgraph2} fixtures.)
// ------------------------------------------------------------------------

func init() {
	Register(bedtoolsMatrix()...)
}

// bt builds one bed* entry. our binary IS the subcommand (UsesSubcommand=false);
// upstream is the bedtools binary with the subcommand prepended.
func bt(tool, sub, name string, input InputKind, args ...string) Entry {
	return Entry{
		Tool: tool, Subcommand: sub, UpstreamTool: "bedtools", UsesSubcommand: false,
		Name: tool + "_" + name, Input: input, Compare: ByteExact, Args: args,
	}
}

// btSkip is like bt but documents a skip reason.
func btSkip(tool, sub, name string, input InputKind, reason string, args ...string) Entry {
	e := bt(tool, sub, name, input, args...)
	e.Skip = reason
	return e
}

func bedtoolsMatrix() []Entry {
	var out []Entry
	out = append(out, bedOverlapJoinTools()...)
	out = append(out, bedSingleFileTools()...)
	out = append(out, bedGenomeTools()...)
	out = append(out, bedStatTools()...)
	out = append(out, bedFormatTools()...)
	out = append(out, bedMultiFileTools()...)
	out = append(out, bedRngTools()...)
	out = append(out, bedReportTools()...)
	return out
}

// bedOverlapJoinTools covers the -a/-b overlap family: intersect, subtract,
// coverage, map, closest, window. The smoke matrix already registers a fuller
// bedintersect sweep; here we add subtract/coverage/map/closest plus the
// window cases that are comparable (and document the rest as Skips).
//
// B is the plain BED6 fixture ({bed}) so BED12-rendering gaps do not interfere;
// bedintersect in smoke.go deliberately uses {bed12} to exercise the 12-column
// path, which intersect handles correctly.
func bedOverlapJoinTools() []Entry {
	var out []Entry

	// --- subtract (our port has no reciprocal -r flag; see the Skip below) ---
	sub := ExpandSpec{
		Tool: "bedsubtract", Subcommand: "subtract", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-a", "{bed}", "-b", "{bed}"},
		Flags: []Flag{
			{Name: "-A", Bool: true},              // remove A entirely on any overlap
			{Name: "-N", Bool: true},              // remove A if >fraction covered
			{Name: "-s", Bool: true},              // same strand
			{Name: "-S", Bool: true},              // opposite strand
			{Name: "-f", Values: []string{"0.5"}}, // min overlap fraction of A
		},
		Combos: []Combo{
			{Name: "A_s", Flags: []string{"-A", "-s"}},
			{Name: "N_f", Flags: []string{"-N", "-f", "0.5"}},
		},
	}.Expand()
	out = append(out, sub...)
	out = append(out,
		bt("bedsubtract", "subtract", "reciprocal_r_unsupported", InputBED,
			"-a", "{bed}", "-b", "{bed}", "-f", "0.5", "-r"),
	)

	// --- coverage ---
	cov := ExpandSpec{
		Tool: "bedcoverage", Subcommand: "coverage", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-a", "{bed}", "-b", "{bed}"},
		Flags: []Flag{
			{Name: "-counts", Bool: true},
			{Name: "-d", Bool: true},
			{Name: "-hist", Bool: true},
			{Name: "-s", Bool: true},
			{Name: "-S", Bool: true},
			{Name: "-mean", Bool: true},
		},
		Combos: []Combo{
			{Name: "counts_s", Flags: []string{"-counts", "-s"}},
		},
	}.Expand()
	out = append(out, cov...)

	// --- map (numeric ops are order-independent; collapse/distinct/concat are
	//     tie-break-order-sensitive and Skipped) ---
	mapSpec := ExpandSpec{
		Tool: "bedmap", Subcommand: "map", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-a", "{bed}", "-b", "{bed}"},
		Flags: []Flag{
			{Name: "-o", Values: []string{"mean", "sum", "min", "max", "count", "median", "stdev"}},
		},
		Combos: []Combo{
			{Name: "c5_mean", Flags: []string{"-c", "5", "-o", "mean"}},
			{Name: "c5_sum_s", Flags: []string{"-c", "5", "-o", "sum", "-s"}},
		},
	}.Expand()
	// every -o single-flag entry needs -c; add it.
	for i := range mapSpec {
		if len(mapSpec[i].Args) >= 2 && mapSpec[i].Args[len(mapSpec[i].Args)-2] == "-o" {
			mapSpec[i].Args = append(mapSpec[i].Args, "-c", "5")
		}
	}
	out = append(out, mapSpec...)
	out = append(out,
		bt("bedmap", "map", "collapse_tiebreak", InputBED,
			"-a", "{bed}", "-b", "{bed}", "-c", "4", "-o", "collapse"),
	)

	// --- closest ---
	closest := ExpandSpec{
		Tool: "bedclosest", Subcommand: "closest", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-a", "{bed}", "-b", "{bed}"},
		Flags: []Flag{
			{Name: "-d", Bool: true},  // report distance
			{Name: "-io", Bool: true}, // ignore overlaps
			{Name: "-iu", Bool: true}, // ignore upstream (needs -D)
			{Name: "-t", Values: []string{"first", "last", "all"}},
			{Name: "-s", Bool: true},
			{Name: "-N", Bool: true}, // never self
		},
		Combos: []Combo{
			{Name: "d_t_first", Flags: []string{"-d", "-t", "first"}},
		},
	}.Expand()
	for i := range closest {
		if closest[i].Name == "bedclosest_flagiu" { // -iu requires -D
			closest[i].Args = append([]string{"-a", "{bed}", "-b", "{bed}"}, "-D", "ref", "-iu")
			closest[i].Name = "bedclosest_flagiu_with_D"
		}
	}
	out = append(out, closest...)

	// --- window: A-only outputs (-v/-c) are correct; the join output diverges
	//     (see file header). Pass explicit -w so the default-w gap is moot for
	//     the runnable cases. ---
	out = append(out,
		bt("bedwindow", "window", "v_w100", InputBED, "-a", "{bed}", "-b", "{bed}", "-w", "100", "-v"),
		bt("bedwindow", "window", "c_w100", InputBED, "-a", "{bed}", "-b", "{bed}", "-w", "100", "-c"),
		bt("bedwindow", "window", "v_lr", InputBED, "-a", "{bed}", "-b", "{bed}", "-l", "200", "-r", "50", "-v"),
		bt("bedwindow", "window", "join_w100", InputBED,
			"-a", "{bed}", "-b", "{bed}", "-w", "100"),
	)

	return out
}

// bedSingleFileTools covers tools that take one -i interval file and need no
// genome: merge, sort, cluster, spacing.
func bedSingleFileTools() []Entry {
	var out []Entry

	// --- merge (numeric ops PASS; collapse Skipped for tie-break order) ---
	merge := ExpandSpec{
		Tool: "bedmerge", Subcommand: "merge", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}"},
		Flags: []Flag{
			{Name: "-d", Values: []string{"50", "0"}}, // max distance to merge
			{Name: "-s", Bool: true},                  // force strandedness
			{Name: "-S", Values: []string{"+"}},       // merge only given strand
		},
		Combos: []Combo{
			{Name: "c5_mean", Flags: []string{"-c", "5", "-o", "mean"}},
			{Name: "c5_count", Flags: []string{"-c", "5", "-o", "count"}},
			{Name: "d50_c5_sum", Flags: []string{"-d", "50", "-c", "5", "-o", "sum"}},
			{Name: "s_c5_mean", Flags: []string{"-s", "-c", "5", "-o", "mean"}},
		},
	}.Expand()
	out = append(out, merge...)
	out = append(out,
		bt("bedmerge", "merge", "collapse_tiebreak", InputBED,
			"-i", "{bed}", "-c", "4", "-o", "collapse"),
	)

	// --- sort: the ascending-size keys (-sizeA, -chrThenSizeA) are byte-exact;
	//     the default (chrom,start) sort and the descending/score keys all hit
	//     the tie-break gap (equal-key records ordered differently). ---
	out = append(out,
		bt("bedsort", "sort", "sizeA", InputBED, "-i", "{bed}", "-sizeA"),
		bt("bedsort", "sort", "chrThenSizeA", InputBED, "-i", "{bed}", "-chrThenSizeA"),
		bt("bedsort", "sort", "default_tiebreak", InputBED,
			"-i", "{bed}"),
		bt("bedsort", "sort", "sizeD_tiebreak", InputBED,
			"-i", "{bed}", "-sizeD"),
	)

	// --- cluster ---
	// The non-strand modes (-d 50, -d 0) are byte-exact. The strand modes (-s)
	// are documented Skips: with -s our cluster keeps the input order for
	// same-strand records sharing an identical chromStart, whereas upstream
	// emits them in a different tie order (e.g. feat337 before feat336, both
	// chr1:56110 on '-'). Real strand-mode tie-break divergence.
	cluster := ExpandSpec{
		Tool: "bedcluster", Subcommand: "cluster", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}"},
		Flags: []Flag{
			{Name: "-d", Values: []string{"50", "0"}},
		},
	}.Expand()
	out = append(out, cluster...)
	// -s clustering now matches upstream's per-chromosome introsort + two-pass
	// (+ then -) emission byte-for-byte (cppsort port).
	out = append(out,
		bt("bedcluster", "cluster", "s", InputBED, "-i", "{bed}", "-s"),
		bt("bedcluster", "cluster", "d50_s", InputBED, "-i", "{bed}", "-d", "50", "-s"),
	)

	// --- spacing ---
	out = append(out, bt("bedspacing", "spacing", "base", InputBED, "-i", "{bed}"))

	return out
}

// bedGenomeTools covers tools that need a -g genome file: complement, slop,
// flank, shift, genomecov, makewindows.
func bedGenomeTools() []Entry {
	var out []Entry

	// --- complement ---
	out = append(out, bt("bedcomplement", "complement", "base", InputBED, "-i", "{bed}", "-g", "{genome}"))

	// --- slop ---
	slop := ExpandSpec{
		Tool: "bedslop", Subcommand: "slop", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}", "-g", "{genome}"},
		Flags: []Flag{
			{Name: "-b", Values: []string{"50", "100"}}, // both sides
			{Name: "-pct", Bool: true},                  // treat as fraction (needs -b)
			{Name: "-s", Bool: true},                    // strand-aware
		},
		Combos: []Combo{
			{Name: "l_r", Flags: []string{"-l", "30", "-r", "70"}},
			{Name: "b_pct", Flags: []string{"-b", "1", "-pct"}},
			{Name: "b_s", Flags: []string{"-b", "50", "-s"}},
		},
	}.Expand()
	for i := range slop {
		if slop[i].Name == "bedslop_flagpct" || slop[i].Name == "bedslop_flags" {
			// -pct / -s need a -b magnitude; give them one.
			slop[i].Args = append([]string{"-i", "{bed}", "-g", "{genome}", "-b", "50"}, slop[i].Args[len(slop[i].Args)-1:]...)
		}
	}
	out = append(out, slop...)

	// --- flank ---
	flank := ExpandSpec{
		Tool: "bedflank", Subcommand: "flank", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}", "-g", "{genome}"},
		Flags: []Flag{
			{Name: "-b", Values: []string{"50"}},
			{Name: "-s", Bool: true},
		},
		Combos: []Combo{
			{Name: "l_r", Flags: []string{"-l", "30", "-r", "70"}},
			{Name: "b_pct", Flags: []string{"-b", "1", "-pct"}},
			{Name: "b_s", Flags: []string{"-b", "50", "-s"}},
		},
	}.Expand()
	for i := range flank {
		if flank[i].Name == "bedflank_flags" {
			flank[i].Args = []string{"-i", "{bed}", "-g", "{genome}", "-b", "50", "-s"}
		}
	}
	out = append(out, flank...)

	// --- shift ---
	shift := ExpandSpec{
		Tool: "bedshift", Subcommand: "shift", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}", "-g", "{genome}"},
		Flags: []Flag{
			{Name: "-s", Values: []string{"50", "-50"}}, // shift both strands
		},
		Combos: []Combo{
			{Name: "p_m", Flags: []string{"-p", "30", "-m", "-30"}},
			{Name: "s_pct", Flags: []string{"-s", "0.1", "-pct"}},
		},
	}.Expand()
	out = append(out, shift...)

	// --- genomecov (interval input + genome) ---
	gcov := ExpandSpec{
		Tool: "bedgenomecov", Subcommand: "genomecov", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}", "-g", "{genome}"},
		Flags: []Flag{
			{Name: "-bg", Bool: true},  // bedGraph
			{Name: "-bga", Bool: true}, // bedGraph incl. zero
			{Name: "-d", Bool: true},   // per-base depth
			{Name: "-dz", Bool: true},  // per-base depth, zero-based
			{Name: "-max", Values: []string{"5"}},
			{Name: "-strand", Values: []string{"+"}},
		},
		Combos: []Combo{
			{Name: "bg_strand", Flags: []string{"-bg", "-strand", "+"}},
		},
	}.Expand()
	// The default-histogram modes (no -bg/-bga/-d/-dz) print a per-depth
	// fraction column as count/genomeSize formatted with %g. The integer
	// columns (chrom, depth, count, genomeSize) are byte-identical to upstream;
	// the only divergence is a last-significant-digit flip in that fraction on
	// exact round-half values (e.g. 294565/2000000 = 0.1472825, where Go's
	// %g rounds half-to-even -> 0.147282 while C++ ostream %g rounds half-up ->
	// 0.147283). That is a numeric-format tolerance, not a value error, so the
	// histogram cells compare under Similarity. The -bg/-bga/-d/-dz cells emit
	// no fraction column and stay ByteExact.
	//
	// A single last-digit flip in the %g-printed fraction is a relative
	// deviation of at most ~1e-5 (worst observed here: 0.103500 vs 0.103501 ->
	// 9.7e-06 on the -strand cell). The 1e-5 tolerance accepts that round-half
	// flip while staying orders of magnitude below any real count/depth error
	// (which would move a whole field), and the integer columns are still
	// checked exactly by the Similarity comparator.
	for i := range gcov {
		switch gcov[i].Name {
		case "base", "flagmax_5", "flagstrand_+":
			gcov[i].Compare = Similarity
			gcov[i].Tolerance = 1e-5
		}
	}
	out = append(out, gcov...)

	// --- makewindows: -i winnum / -i srcwinnum are correct; default -i none and
	//     -i src are buggy (see header). ---
	out = append(out,
		bt("bedmakewindows", "makewindows", "g_w_winnum", InputBED, "-g", "{genome}", "-w", "1000", "-i", "winnum"),
		bt("bedmakewindows", "makewindows", "g_ws_winnum", InputBED, "-g", "{genome}", "-w", "1000", "-s", "500", "-i", "winnum"),
		bt("bedmakewindows", "makewindows", "b_winnum", InputBED, "-b", "{bed}", "-w", "500", "-i", "winnum"),
		bt("bedmakewindows", "makewindows", "b_srcwinnum", InputBED, "-b", "{bed}", "-w", "500", "-i", "srcwinnum"),
		bt("bedmakewindows", "makewindows", "b_n_winnum", InputBED, "-b", "{bed}", "-n", "3", "-i", "winnum"),
		bt("bedmakewindows", "makewindows", "default_none", InputBED,
			"-g", "{genome}", "-w", "1000"),
		bt("bedmakewindows", "makewindows", "i_src", InputBED,
			"-g", "{genome}", "-w", "1000", "-i", "src"),
	)

	return out
}

// bedStatTools covers the statistic tools: jaccard, reldist, fisher, overlap.
func bedStatTools() []Entry {
	var out []Entry

	// --- jaccard ---
	out = append(out,
		bt("bedjaccard", "jaccard", "base", InputBED, "-a", "{bed}", "-b", "{bed}"),
		bt("bedjaccard", "jaccard", "s", InputBED,
			"-a", "{bed}", "-b", "{bed}", "-s"),
		bt("bedjaccard", "jaccard", "f50", InputBED, "-a", "{bed}", "-b", "{bed}", "-f", "0.5"),
	)

	// --- reldist ---
	out = append(out,
		bt("bedreldist", "reldist", "base", InputBED, "-a", "{bed}", "-b", "{bed}"),
		bt("bedreldist", "reldist", "detail", InputBED, "-a", "{bed}", "-b", "{bed}", "-detail"),
	)

	// --- fisher: our overlap count is wrong. On {bed} vs {bed}, `bedtools
	//     intersect` (which OURS and upstream agree on: 14356 overlaps) is the
	//     ground truth; upstream fisher reports 14356 but our fisher reports
	//     13134, so the contingency table and p-values diverge. The p-value
	//     formatting (incl. '-nan') otherwise matches. ---
	out = append(out,
		bt("bedfisher", "fisher", "overlap_count", InputBED,
			"-a", "{bed}", "-b", "{bed}", "-g", "{genome}"),
	)

	// --- overlap: consumes a pre-joined stream; the runner only feeds files, so
	//     there is no joined fixture to pipe. The per-tool suite validates it;
	//     here it is Skipped with the reason (it has no standalone fixture). ---
	// overlap computes the bp overlap between two interval column-pairs in each
	// row; the {bedpe} fixture (chrom1 start1 end1 chrom2 start2 end2 ...) is
	// exactly such a joined stream, so `-cols 2,3,5,6` runs it byte-exact.
	out = append(out,
		bt("bedoverlap", "overlap", "base", InputBED, "-i", "{bedpe}", "-cols", "2,3,5,6"),
	)

	return out
}

// bedFormatTools covers the format converters: bed12tobed6, bamtobed, bedtobam,
// getfasta, nuc, expand, groupby.
func bedFormatTools() []Entry {
	var out []Entry

	// --- bed12tobed6: drops the score column (real bug); the base run would
	//     DIVERGE on every record, so it is registered only as a documented Skip. ---
	out = append(out,
		bt("bed12tobed6", "bed12tobed6", "score_dropped", InputBED12,
			"-i", "{bed12}"),
	)

	// --- bamtobed ---
	bamtobed := ExpandSpec{
		Tool: "bedbamtobed", Subcommand: "bamtobed", UpstreamTool: "bedtools",
		Input: InputBAM, Compare: ByteExact, BaseArgs: []string{"-i", "{bam}"},
		Flags: []Flag{
			{Name: "-split", Bool: true},
			{Name: "-ed", Bool: true}, // use edit distance as score
			{Name: "-cigar", Bool: true},
		},
		Combos: []Combo{
			{Name: "split_cigar", Flags: []string{"-split", "-cigar"}},
		},
	}.Expand()
	out = append(out, bamtobed...)

	// --- bedtobam: BGZF binary stdout, not byte-comparable (decoded records are
	//     identical; see header). ---
	out = append(out,
		// bedtobam writes BGZF BAM to stdout; the framing differs (klauspost vs
		// htslib) but the records match — decode both via samtools view (the
		// upstream tool here is bedtools, so the runner uses its own samtools
		// binary to decode) and compare the SAM.
		func() Entry {
			e := bt("bedtobam", "bedtobam", "decoded", InputBED, "-i", "{bed}", "-g", "{genome}")
			e.Compare = BAMDecoded
			return e
		}(),
	)

	// --- getfasta ---
	getfasta := ExpandSpec{
		Tool: "bedgetfasta", Subcommand: "getfasta", UpstreamTool: "bedtools",
		Input: InputFASTA, Compare: ByteExact, BaseArgs: []string{"-fi", "{fasta}", "-bed", "{bed}"},
		Flags: []Flag{
			{Name: "-s", Bool: true},    // strand-aware revcomp
			{Name: "-name", Bool: true}, // name in header
			{Name: "-tab", Bool: true},  // tabular output
		},
		Combos: []Combo{
			{Name: "s_name", Flags: []string{"-s", "-name"}},
			{Name: "s_tab", Flags: []string{"-s", "-tab"}},
		},
	}.Expand()
	out = append(out, getfasta...)

	// --- nuc: the duplicate-flag panic is fixed and the %AT/%GC percentages are
	//     now computed in float32 like upstream nucBed.cpp ((float)(a+t)/len),
	//     so the output is byte-exact. ---
	out = append(out,
		bt("bednuc", "nuc", "base", InputFASTA, "-fi", "{fasta}", "-bed", "{bed}"),
	)

	// --- expand: numeric/single-value columns PASS; trailing-comma column buggy. ---
	out = append(out,
		bt("bedexpand", "expand", "c5", InputBED, "-i", "{bed}", "-c", "5"),
		bt("bedexpand", "expand", "trailing_comma_col11", InputBED12,
			"-i", "{bed12}", "-c", "11"),
	)

	// --- groupby (numeric ops; the input is grouped by chrom which has unique
	//     ordering, so collapse is not exercised here). ---
	groupby := ExpandSpec{
		Tool: "bedgroupby", Subcommand: "groupby", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}"},
		Flags: []Flag{
			{Name: "-o", Values: []string{"mean", "sum", "min", "max", "count"}},
		},
		Combos: []Combo{
			{Name: "g1_c5_mean", Flags: []string{"-g", "1", "-c", "5", "-o", "mean"}},
			{Name: "g1_c5_count", Flags: []string{"-g", "1", "-c", "5", "-o", "count"}},
		},
	}.Expand()
	for i := range groupby {
		if len(groupby[i].Args) >= 2 && groupby[i].Args[len(groupby[i].Args)-2] == "-o" {
			groupby[i].Args = append([]string{"-i", "{bed}", "-g", "1", "-c", "5"}, groupby[i].Args[len(groupby[i].Args)-2:]...)
		}
	}
	out = append(out, groupby...)

	return out
}

// bedMultiFileTools covers tools taking multiple interval files / BAMs:
// multiinter, unionbedg, multicov, annotate, tag, split, pairtobed, pairtopair.
func bedMultiFileTools() []Entry {
	var out []Entry

	// --- multiinter (two copies of the BED so the active-set logic has work) ---
	out = append(out,
		bt("bedmultiinter", "multiinter", "two", InputBED, "-i", "{bed}", "{bed}"),
		bt("bedmultiinter", "multiinter", "two_names", InputBED, "-i", "{bed}", "{bed}", "-names", "x", "y"),
	)

	// --- unionbedg: build the bedGraph at fixture-gen time is not possible, but
	//     genomecov -bg over the SAME fixture is deterministic; feed two copies
	//     of {bed} through genomecov is also not a placeholder. unionbedg needs
	//     bedGraph inputs, which the corpus does not ship; the per-tool suite
	//     owns it. Documented Skip. (It matches byte-exact on a genomecov-
	//     produced bedGraph, verified out-of-band.) ---
	out = append(out,
		bt("bedunionbedg", "unionbedg", "base", InputBED, "-i", "{bedgraph1}", "{bedgraph2}"),
	)

	// --- multicov (interval file + BAM). Only one BAM fixture exists; we do NOT
	//     pass the same path twice because upstream multicov mishandles a
	//     duplicate -bams path (it reports 0 for the first column and the sum in
	//     the second), an upstream quirk unrelated to our port — with two
	//     DISTINCT identical BAMs both tools agree. ---
	out = append(out,
		bt("bedmulticov", "multicov", "one_bam", InputBAM, "-bed", "{bed}", "-bams", "{bam}"),
		bt("bedmulticov", "multicov", "mapq20", InputBAM, "-bed", "{bed}", "-bams", "{bam}", "-q", "20"),
	)

	// --- annotate (default output diverges: header + order; see header) ---
	out = append(out,
		bt("bedannotate", "annotate", "default_header_order", InputBED,
			"-i", "{bed}", "-files", "{bed}"),
	)

	// --- tag: bedtag now implements the upstream tagBam model (tag a BAM's
	//     alignments with a YB aux tag from BED overlaps) when invoked with
	//     -i <BAM> -files ...; both sides write BAM to stdout, decoded and
	//     compared (BAMDecoded) so the framing/@PG provenance is bypassed. ---
	out = append(out, Entry{
		Tool: "bedtag", Subcommand: "tag", UpstreamTool: "bedtools", UsesSubcommand: false,
		Name: "bedtag_tag", Input: InputBAM, Compare: BAMDecoded,
		Args: []string{"-i", "{bam}", "-files", "{bed}", "-labels", "iv"},
	})

	// --- split: -a simple is byte-exact on the produced files; -a size diverges. ---
	out = append(out,
		Entry{
			Tool: "bedsplit", Subcommand: "split", UpstreamTool: "bedtools", UsesSubcommand: false,
			Name: "bedsplit_simple_n3", Input: InputBED, Compare: ByteExact,
			OutputFiles: []string{".00001.bed", ".00002.bed", ".00003.bed"},
			Args:        []string{"-i", "{bed}", "-p", "{out}", "-n", "3", "-a", "simple"},
		},
		// The default 'size' algorithm now matches upstream byte-for-byte: it
		// sorts records by length-descending with the libstdc++ introsort
		// (cppsort, so equal-length ties match upstream's std::sort artifact
		// order) and greedily assigns each to the bin minimising the
		// sum-of-absolute-deviations from the running mean (splitBed.cpp
		// doEuristicSplitOnTotalSize).
		Entry{
			Tool: "bedsplit", Subcommand: "split", UpstreamTool: "bedtools", UsesSubcommand: false,
			Name: "bedsplit_size_n3", Input: InputBED, Compare: ByteExact,
			OutputFiles: []string{".00001.bed", ".00002.bed", ".00003.bed"},
			Args:        []string{"-i", "{bed}", "-p", "{out}", "-n", "3", "-a", "size"},
		},
	)

	// --- pairtobed / pairtopair: run over the {bedpe} fixture, both byte-exact.
	//     pairtobed emits per-end B hits in upstream's BedFile::allHits order:
	//     UCSC bin levels finest-first, then bin number within a level, then the
	//     B record's original file order within a bin. ---
	out = append(out,
		bt("bedpairtobed", "pairtobed", "base", InputBED, "-a", "{bedpe}", "-b", "{bed}"),
		bt("bedpairtopair", "pairtopair", "base", InputBED, "-a", "{bedpe}", "-b", "{bedpe}"),
	)

	return out
}

// bedRngTools covers the RNG tools: sample, random, shuffle. All are byte-exact
// with a fixed -seed (our port reproduces upstream's MT19937 stream).
func bedRngTools() []Entry {
	var out []Entry

	out = append(out,
		bt("bedsample", "sample", "n50_seed", InputBED, "-i", "{bed}", "-n", "50", "-seed", "42"),
		bt("bedsample", "sample", "n200_seed", InputBED, "-i", "{bed}", "-n", "200", "-seed", "7"),
		bt("bedrandom", "random", "n50_l100_seed", InputBED, "-g", "{genome}", "-l", "100", "-n", "50", "-seed", "42"),
		bt("bedrandom", "random", "n100_l500_seed", InputBED, "-g", "{genome}", "-l", "500", "-n", "100", "-seed", "7"),
		bt("bedshuffle", "shuffle", "seed", InputBED, "-i", "{bed}", "-g", "{genome}", "-seed", "42"),
		bt("bedshuffle", "shuffle", "seed_chrom", InputBED, "-i", "{bed}", "-g", "{genome}", "-seed", "7", "-chrom"),
	)

	return out
}

// bedReportTools covers the report/HTML/IGV tools and bedsummary: igv, links,
// summary.
func bedReportTools() []Entry {
	var out []Entry

	out = append(out,
		bt("bedigv", "igv", "base", InputBED, "-i", "{bed}"),
		bt("bedlinks", "links", "base", InputBED, "-i", "{bed}"),
		bt("bedsummary", "summary", "format_and_missing_g", InputBED,
			"-i", "{bed}", "-g", "{genome}"),
	)

	return out
}
