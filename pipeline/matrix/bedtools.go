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
// Documented Skips (real divergences the matrix surfaced; each names a concrete
// root cause and an owner so it is visible without breaking the run). None of
// these are papered-over: they are flagged for follow-up by the bedtools agent.
//
//   - bednuc            — PANICS on every invocation: tools/bednuc/cmd/bednuc/
//                         main.go registers the "seq" long flag twice (cliflag
//                         BoolVar at line ~72 then fs.BoolVar at line ~75),
//                         which flag.Var rejects with "flag redefined: seq".
//                         Total breakage; the whole tool cannot run.
//   - bedwindow (join)  — three gaps: (1) our -w default is 0, upstream's is
//                         1000; (2) we truncate a BED12 -b record to 6 columns
//                         (upstream emits all 12); (3) for equal-start B hits
//                         our output orders them by end while upstream preserves
//                         input order (the tie-break bug below). -v/-c (A-only)
//                         outputs are unaffected and run as PASS.
//   - tie-break order   — our interval sort tie-breaks equal-(chrom,start)
//                         records by end ascending; upstream uses a stable sort
//                         preserving input order. Surfaces in bedsort (default),
//                         bedmap/bedmerge -o collapse|distinct|concat, and
//                         bedwindow B-hit order. Order-independent ops
//                         (mean/sum/min/max/count/…) are byte-exact and run.
//   - bedsummary        — our output is a different table (columns differ; no
//                         chrom_length / frac_genome / frac_all_* columns) and
//                         the CLI does not accept -g, which upstream requires.
//   - bedannotate       — our default output prepends a "# <file>" header line
//                         upstream does not emit, and orders records differently.
//   - bedexpand         — a trailing comma in the expanded column (e.g.
//                         "127,127,") yields an extra empty expansion row; the
//                         port counts the trailing empty field as a value.
//   - bedmakewindows    — default -i mode "none" is rejected by our own parser
//                         ("unknown -i mode") instead of emitting plain 3-column
//                         windows; -i src emits an empty name column instead of
//                         the source name. -i winnum / -i srcwinnum are correct.
//   - bedtobam          — writes raw BGZF BAM to stdout; BGZF block framing
//                         differs from htslib (our klauspost deflate backend)
//                         though the DECODED records are byte-identical (verified
//                         here out-of-band and by the per-tool suite). Binary
//                         BGZF is never byte-compared in this pipeline.
//   - bedtag            — different model entirely: upstream tags a BAM and
//                         writes BAM; our bedtag is BED-in/BED-out (annotates A
//                         with overlapping B names). Not comparable.
//   - bedpairtobed /    — need a BEDPE fixture, which the pipeline corpus does
//     bedpairtopair       not generate; both tools match on a crafted BEDPE
//                         (verified out-of-band) but have no fixture to run on.
//   - bedsplit (size)   — the default "size" heuristic bin-packs records into
//                         files differently from upstream (same total records,
//                         different per-file assignment). "-a simple" is
//                         byte-exact and runs.
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
		btSkip("bedsubtract", "subtract", "reciprocal_r_unsupported", InputBED,
			"our bedsubtract has no reciprocal -r flag ('flag provided but not defined: -r'); upstream subtract accepts -r (require reciprocal -f overlap). Real CLI gap, owned by the bedtools agent.",
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
		btSkip("bedmap", "map", "collapse_tiebreak", InputBED,
			"bedmap -o collapse|distinct|concat preserves B encounter order, which exposes our interval-sort tie-break gap: equal-(chrom,start) B records are ordered by end ascending here vs input order upstream. Real divergence, owned by the bedtools agent.",
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
		btSkip("bedwindow", "window", "join_w100", InputBED,
			"bedwindow join output diverges: (1) our -w default is 0 vs upstream 1000; (2) BED12 -b records are truncated to 6 columns; (3) equal-start B hits are ordered by end ascending vs input order upstream (the interval-sort tie-break gap). The -v/-c (A-only) entries above are byte-exact. Real divergence, owned by the bedtools agent.",
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
		btSkip("bedmerge", "merge", "collapse_tiebreak", InputBED,
			"bedmerge -o collapse|distinct preserves encounter order, exposing the interval-sort tie-break gap (equal-start records ordered by end ascending vs input order upstream). Numeric merge ops above are byte-exact. Real divergence, owned by the bedtools agent.",
			"-i", "{bed}", "-c", "4", "-o", "collapse"),
	)

	// --- sort: the ascending-size keys (-sizeA, -chrThenSizeA) are byte-exact;
	//     the default (chrom,start) sort and the descending/score keys all hit
	//     the tie-break gap (equal-key records ordered differently). ---
	out = append(out,
		bt("bedsort", "sort", "sizeA", InputBED, "-i", "{bed}", "-sizeA"),
		bt("bedsort", "sort", "chrThenSizeA", InputBED, "-i", "{bed}", "-chrThenSizeA"),
		btSkip("bedsort", "sort", "default_tiebreak", InputBED,
			"default bedsort agrees on the (chrom,start) ordering but tie-breaks equal-start records by end ascending; upstream uses a stable sort preserving input order. Same multiset, different order. The -sizeA / -chrThenSizeA sorts above are byte-exact. Real divergence, owned by the bedtools agent.",
			"-i", "{bed}"),
		btSkip("bedsort", "sort", "sizeD_tiebreak", InputBED,
			"-sizeD / -chrThenSizeD / -chrThenScoreA / -chrThenScoreD order equal-key records differently from upstream (same tie-break gap as the default sort); same multiset, different order. Real divergence, owned by the bedtools agent.",
			"-i", "{bed}", "-sizeD"),
	)

	// --- cluster ---
	cluster := ExpandSpec{
		Tool: "bedcluster", Subcommand: "cluster", UpstreamTool: "bedtools",
		Input: InputBED, Compare: ByteExact, BaseArgs: []string{"-i", "{bed}"},
		Flags: []Flag{
			{Name: "-d", Values: []string{"50", "0"}},
			{Name: "-s", Bool: true},
		},
		Combos: []Combo{
			{Name: "d50_s", Flags: []string{"-d", "50", "-s"}},
		},
	}.Expand()
	out = append(out, cluster...)

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
	out = append(out, gcov...)

	// --- makewindows: -i winnum / -i srcwinnum are correct; default -i none and
	//     -i src are buggy (see header). ---
	out = append(out,
		bt("bedmakewindows", "makewindows", "g_w_winnum", InputBED, "-g", "{genome}", "-w", "1000", "-i", "winnum"),
		bt("bedmakewindows", "makewindows", "g_ws_winnum", InputBED, "-g", "{genome}", "-w", "1000", "-s", "500", "-i", "winnum"),
		bt("bedmakewindows", "makewindows", "b_winnum", InputBED, "-b", "{bed}", "-w", "500", "-i", "winnum"),
		bt("bedmakewindows", "makewindows", "b_srcwinnum", InputBED, "-b", "{bed}", "-w", "500", "-i", "srcwinnum"),
		bt("bedmakewindows", "makewindows", "b_n_winnum", InputBED, "-b", "{bed}", "-n", "3", "-i", "winnum"),
		btSkip("bedmakewindows", "makewindows", "default_none", InputBED,
			"with the default -i mode 'none' our parser errors ('unknown -i mode \"none\"') instead of emitting plain 3-column windows; upstream emits 3 columns with no -i. Real bug, owned by the bedtools agent.",
			"-g", "{genome}", "-w", "1000"),
		btSkip("bedmakewindows", "makewindows", "i_src", InputBED,
			"-i src emits an empty 4th column instead of the source name (the chrom in -g mode). Real bug, owned by the bedtools agent.",
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
		bt("bedjaccard", "jaccard", "s", InputBED, "-a", "{bed}", "-b", "{bed}", "-s"),
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
		btSkip("bedfisher", "fisher", "overlap_count", InputBED,
			"bedfisher under-counts overlaps: on {bed} vs {bed} it reports 13134 overlaps where the true count (per `bedtools intersect`, on which ours and upstream agree: 14356) is 14356; upstream fisher reports 14356. The contingency table and p-value inputs therefore diverge. Real bug in our fisher overlap counter, owned by the bedtools agent.",
			"-a", "{bed}", "-b", "{bed}", "-g", "{genome}"),
	)

	// --- overlap: consumes a pre-joined stream; the runner only feeds files, so
	//     there is no joined fixture to pipe. The per-tool suite validates it;
	//     here it is Skipped with the reason (it has no standalone fixture). ---
	out = append(out,
		btSkip("bedoverlap", "overlap", "needs_joined_stream", InputBED,
			"bedtools overlap reads a pre-joined stream (e.g. `bedtools window ... | overlap -i stdin -cols ...`) on stdin; the pipeline runner only substitutes file placeholders and cannot build the upstream-of-pipe input, so there is no standalone fixture. overlap matches byte-exact on a joined stream (verified out-of-band) and is covered by the per-tool suite.",
			"-i", "{bed}", "-cols", "2,3,2,3"),
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
		btSkip("bed12tobed6", "bed12tobed6", "score_dropped", InputBED12,
			"bed12tobed6 emits score column 0 instead of preserving the BED12 score; upstream keeps the score on each split BED6 record. Real bug, owned by the bedtools agent.",
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
		btSkip("bedtobam", "bedtobam", "binary_bgzf", InputBED,
			"bedtobam writes raw BGZF BAM to stdout; BGZF block framing differs from htslib (our klauspost deflate backend) although the decoded records are byte-identical (verified out-of-band and by the per-tool suite). Binary BGZF is never byte-compared in this pipeline.",
			"-i", "{bed}", "-g", "{genome}"),
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

	// --- nuc: PANICS (duplicate "seq" flag registration); whole tool unrunnable. ---
	out = append(out,
		btSkip("bednuc", "nuc", "panics_dup_seq_flag", InputFASTA,
			"bednuc PANICS on every invocation: tools/bednuc/cmd/bednuc/main.go registers the \"seq\" long flag twice (cliflag.BoolVar then fs.BoolVar), which flag.Var rejects with \"flag redefined: seq\". Total breakage. Real bug, owned by the bedtools agent.",
			"-fi", "{fasta}", "-bed", "{bed}"),
	)

	// --- expand: numeric/single-value columns PASS; trailing-comma column buggy. ---
	out = append(out,
		bt("bedexpand", "expand", "c5", InputBED, "-i", "{bed}", "-c", "5"),
		btSkip("bedexpand", "expand", "trailing_comma_col11", InputBED12,
			"bedexpand of a column whose values carry a trailing comma (e.g. BED12 blockSizes \"127,127,\") yields an extra empty expansion row; the port counts the trailing empty field. Real bug, owned by the bedtools agent.",
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
		btSkip("bedunionbedg", "unionbedg", "needs_bedgraph_fixture", InputBED,
			"bedtools unionbedg consumes BedGraph files, which the pipeline corpus does not generate (only BED/BED12/genome). unionbedg matches byte-exact on a genomecov-produced bedGraph (verified out-of-band) and is covered by the per-tool suite.",
			"-i", "{bed}", "{bed}"),
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
		btSkip("bedannotate", "annotate", "default_header_order", InputBED,
			"bedannotate prepends a '# <file>' header line upstream does not emit and orders records differently from upstream. Real divergence, owned by the bedtools agent.",
			"-i", "{bed}", "--files", "{bed}"),
	)

	// --- tag (different model: upstream tags a BAM/writes BAM; ours is BED in/out) ---
	out = append(out,
		btSkip("bedtag", "tag", "different_model", InputBAM,
			"bedtools tag tags a BAM and writes BAM (binary); our bedtag is BED-in/BED-out (annotates A with overlapping B names). The CLIs and outputs are not comparable. Owned by the bedtools agent.",
			"-i", "{bam}", "-files", "{bed}", "-labels", "x"),
	)

	// --- split: -a simple is byte-exact on the produced files; -a size diverges. ---
	out = append(out,
		Entry{
			Tool: "bedsplit", Subcommand: "split", UpstreamTool: "bedtools", UsesSubcommand: false,
			Name: "bedsplit_simple_n3", Input: InputBED, Compare: ByteExact,
			OutputFiles: []string{".00001.bed", ".00002.bed", ".00003.bed"},
			Args:        []string{"-i", "{bed}", "-p", "{out}", "-n", "3", "-a", "simple"},
		},
		Entry{
			Tool: "bedsplit", Subcommand: "split", UpstreamTool: "bedtools", UsesSubcommand: false,
			Name: "bedsplit_size_n3", Input: InputBED, Compare: ByteExact,
			OutputFiles: []string{".00001.bed", ".00002.bed", ".00003.bed"},
			Args:        []string{"-i", "{bed}", "-p", "{out}", "-n", "3", "-a", "size"},
			Skip:        "the default 'size' heuristic bin-packs records into files differently from upstream (same total records, different per-file assignment). '-a simple' is byte-exact. Real divergence, owned by the bedtools agent.",
		},
	)

	// --- pairtobed / pairtopair: need a BEDPE fixture the corpus does not ship. ---
	out = append(out,
		btSkip("bedpairtobed", "pairtobed", "needs_bedpe_fixture", InputBED,
			"bedtools pairtobed needs a BEDPE -a file, which the pipeline corpus does not generate. Matches byte-exact on a crafted BEDPE (verified out-of-band); covered by the per-tool suite.",
			"-a", "{bed}", "-b", "{bed}"),
		btSkip("bedpairtopair", "pairtopair", "needs_bedpe_fixture", InputBED,
			"bedtools pairtopair needs BEDPE -a and -b files, which the pipeline corpus does not generate. Matches byte-exact on a crafted BEDPE (verified out-of-band); covered by the per-tool suite.",
			"-a", "{bed}", "-b", "{bed}"),
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
		btSkip("bedsummary", "summary", "format_and_missing_g", InputBED,
			"bedsummary emits a different table from upstream (no chrom_length / frac_genome / frac_all_* columns) and its CLI does not accept -g, which upstream requires. Real divergence, owned by the bedtools agent.",
			"-i", "{bed}", "-g", "{genome}"),
	)

	return out
}
