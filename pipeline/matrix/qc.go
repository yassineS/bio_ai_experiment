package matrix

// This file registers the QC / format tool matrices: seqtk, prinseq, sickle,
// skewer, and fastp. They consume the FASTA / FASTQ (single-end and paired-end)
// fixtures the generator produces.
//
// Comparison strategy per tool (see each builder's doc comment for detail):
//
//   - seqtk:   byte-exact on stdout. Our seqtk is subcommand-based exactly like
//              upstream, so a shared Args works for both sides.
//   - prinseq: byte-exact on the produced output file. Our prinseq is
//              subcommand-based (filter) with -i/-o flags while upstream is the
//              flat prinseq-lite.pl; OurArgs / UpstreamArgs express both sides.
//   - sickle:  byte-exact on the produced output file. Driven with -w 0 so our
//              CLI uses upstream's dynamic int(0.1*len) window (see the note on
//              the sickle default-window divergence below).
//   - skewer:  byte-exact on stdout. Our skewer is subcommand-based with -i/-o;
//              upstream is flat with positionals + -1/--stdout, so the two sides
//              need OurArgs / UpstreamArgs.
//   - fastp:   mostly documented Skips — fastp's heuristic / default-filter /
//              cut-window-boundary behaviour and several missing
//              disable-filter flags make generic byte-exact comparison
//              impractical here; the per-tool suite owns its byte-exact
//              validation. A couple of entries run as Similarity.

func init() {
	Register(seqtkMatrix()...)
	Register(prinseqMatrix()...)
	Register(sickleMatrix()...)
	Register(skewerMatrix()...)
	Register(fastpMatrix()...)
}

// seqtkMatrix exercises seqtk subcommands over the FASTA and FASTQ fixtures.
// All entries are byte-exact on stdout. Our seqtk takes the subcommand as its
// first argument exactly like upstream (UsesSubcommand=true), so the shared
// Args serves both sides. Subcommands whose output is a documented formatting
// divergence (cutN's unwrapped sequence) or that have no upstream counterpart
// (our fq2fa extension) are excluded; seq -A covers FASTQ->FASTA.
func seqtkMatrix() []Entry {
	fa := "{fasta}"
	fq := "{fastq}"

	var entries []Entry
	// Build the seq sweep with the expander (FASTQ input).
	seq := ExpandSpec{
		Tool: "seqtk", Subcommand: "seq", UsesSubcommand: true,
		Input: InputFASTQ, Compare: ByteExact, BaseArgs: []string{fq},
		Flags: []Flag{
			{Name: "-A", Bool: true},             // FASTQ -> FASTA
			{Name: "-r", Bool: true},             // reverse complement
			{Name: "-L", Values: []string{"95"}}, // drop reads shorter than 95
			{Name: "-q", Values: []string{"20"}}, // mask low-quality bases
			{Name: "-l", Values: []string{"60"}}, // line wrap
		},
		Combos: []Combo{
			{Name: "rev_upper", Flags: []string{"-r", "-U"}},
			{Name: "qmask_n", Flags: []string{"-q", "20", "-n", "N"}},
		},
	}.Expand()
	for i := range seq {
		seq[i].Name = "seqtk_seq_fq_" + seq[i].Name
	}
	entries = append(entries, seq...)

	// --- comp, fqchk, size, trimfq, sample, hpc, gap, subseq, mergepe,
	//     dropse, randbase, telo, listhet, hety ---
	add := func(name, sub string, input InputKind, cmp CompareMode, args ...string) {
		entries = append(entries, Entry{
			Tool: "seqtk", Subcommand: sub, UsesSubcommand: true,
			Name: "seqtk_" + name, Input: input, Compare: cmp, Args: args,
		})
	}
	add("comp_fa", "comp", InputFASTA, ByteExact, fa)
	add("comp_fq", "comp", InputFASTQ, ByteExact, fq)
	add("fqchk", "fqchk", InputFASTQ, ByteExact, fq)
	add("fqchk_q20", "fqchk", InputFASTQ, ByteExact, "-q", "20", fq)
	add("size_fq", "size", InputFASTQ, ByteExact, fq)
	add("size_fa", "size", InputFASTA, ByteExact, fa)
	add("trimfq", "trimfq", InputFASTQ, ByteExact, fq)
	add("trimfq_q", "trimfq", InputFASTQ, ByteExact, "-q", "0.05", fq)
	add("trimfq_be", "trimfq", InputFASTQ, ByteExact, "-b", "5", "-e", "5", fq)
	add("sample_count", "sample", InputFASTQ, ByteExact, "-s11", fq, "200")
	add("sample_frac", "sample", InputFASTQ, ByteExact, "-s11", fq, "0.1")
	add("hpc_fa", "hpc", InputFASTA, ByteExact, fa)
	add("hpc_fq", "hpc", InputFASTQ, ByteExact, fq)
	add("gap_fa", "gap", InputFASTA, ByteExact, fa)
	add("subseq_bed", "subseq", InputFASTA, ByteExact, fa, "{bed}")
	add("mergepe", "mergepe", InputFASTQPaired, ByteExact, "{fastq1}", "{fastq2}")
	add("dropse", "dropse", InputFASTQ, ByteExact, fq)
	add("randbase", "randbase", InputFASTA, ByteExact, fa)
	add("telo", "telo", InputFASTA, ByteExact, fa)
	add("listhet", "listhet", InputFASTA, ByteExact, fa)
	add("hety", "hety", InputFASTA, ByteExact, fa)

	// seq over the FASTA reference (baseline transform) and a heavy full
	// FASTQ->FASTA pass for the timing ratio.
	add("seq_fa", "seq", InputFASTA, ByteExact, fa)
	entries = append(entries, Entry{
		Tool: "seqtk", Subcommand: "seq", UsesSubcommand: true,
		Name: "seqtk_seq_fq_to_fa_heavy", Input: InputFASTQ, Compare: ByteExact,
		Args: []string{"-A", fq}, Heavy: true,
	})

	// cutN: byte-exact. Our cutN now wraps each emitted fragment at 60 bases,
	// matching upstream print_seq's `(i-begin)%60` rule (previously our writer
	// emitted the fragment unwrapped on a single line).
	add("cutN", "cutN", InputFASTA, ByteExact, "-n", "100", fa)

	return entries
}

// prinseqMatrix exercises our prinseq `filter` subcommand against the upstream
// prinseq-lite.pl Perl script. Both write a good-reads FASTQ; the comparison is
// byte-exact on that output file. The two CLIs differ in shape (ours is
// subcommand-based with -i/-o; upstream is flat with -fastq/-out_good), so
// OurArgs / UpstreamArgs express each side. The {out} placeholder resolves to a
// per-invocation prefix; ours writes "<prefix>.fastq" via -o while upstream's
// -out_good "<prefix>" yields "<prefix>.fastq" too, so OutputFiles is
// [".fastq"]. (-n 0 is avoided: our filter treats max-ns 0 as "disabled" where
// this prinseq build treats it as "drop any N", a divergence owned by the
// prinseq agent.)
func prinseqMatrix() []Entry {
	fq := "{fastq}"
	mk := func(name string, up, ours []string) Entry {
		return Entry{
			Tool: "prinseq", UpstreamTool: "prinseq", Name: "prinseq_" + name,
			Input: InputFASTQ, Compare: ByteExact, OutputFiles: []string{".fastq"},
			// Upstream prinseq-lite.pl: flat flags writing -out_good <prefix>.
			UpstreamArgs: append([]string{"-fastq", fq, "-out_bad", "null", "-out_good", "{out}"}, up...),
			// Ours: filter subcommand, -o <prefix>.fastq, FASTQ out_format 3.
			OurArgs: append([]string{"filter", "-i", fq, "--fastq", "-o", "{out}.fastq", "--out_format", "3"}, ours...),
		}
	}
	// No "base" (no-flag) entry: with no filter/trim, upstream prinseq-lite.pl
	// writes no out_good file at all while our filter passes every read
	// through, so the no-op case is not a meaningful comparison.
	return []Entry{
		mk("min_len", []string{"-min_len", "95"}, []string{"-l", "95"}),
		// max_len 160 keeps a real subset (~33k of 40k; reads are 135-165 bp).
		// A lower threshold like 105 would drop EVERY read, and upstream then
		// writes no out_good file at all while our filter writes an empty one —
		// a genuine but narrow empty-output edge case, avoided here by choosing a
		// threshold that keeps reads so the real length filtering is exercised.
		mk("max_len", []string{"-max_len", "160"}, []string{"-L", "160"}),
		mk("trim_left", []string{"-trim_left", "5"}, []string{"--trim-left", "5"}),
		mk("trim_right", []string{"-trim_right", "5"}, []string{"--trim-right", "5"}),
		mk("trim_qual_right", []string{"-trim_qual_right", "20"}, []string{"--trim-qual-right", "20"}),
		mk("trim_qual_left", []string{"-trim_qual_left", "20"}, []string{"--trim-qual-left", "20"}),
		mk("min_qual_mean", []string{"-min_qual_mean", "25"}, []string{"-q", "25"}),
		mk("max_ns", []string{"-ns_max_n", "2"}, []string{"-n", "2"}),
	}
}

// sickleMatrix exercises sickle se/pe against the upstream C binary, comparing
// the produced trimmed FASTQ byte-for-byte.
//
// IMPORTANT — sickle default-window divergence (a real, root-caused bug):
// upstream sickle's sliding window defaults to int(0.1*read_length); our CLI
// hardcodes a default window of 10 (tools/sickle/cmd/sickle/main.go), which
// defeats the library's dynamic fallback and trims ~1% of reads one window
// short. We therefore drive every comparison entry with `-w 0`, which selects
// our library's upstream-faithful dynamic window and yields byte-exact output.
// A single baseline entry runs at the CLI default and is Skipped with this
// explanation so the divergence is recorded (and visible to the sickle agent)
// without failing the run.
func sickleMatrix() []Entry {
	fq := "{fastq}"
	r1, r2 := "{fastq1}", "{fastq2}"
	se := func(name string, extra ...string) Entry {
		args := append([]string{"se", "-f", fq, "-t", "sanger", "-w", "0", "-o", "{out}.fastq"}, extra...)
		up := append([]string{"se", "-f", fq, "-t", "sanger", "-o", "{out}.fastq"}, extra...)
		return Entry{
			Tool: "sickle", UpstreamTool: "sickle", Name: "sickle_se_" + name,
			Input: InputFASTQ, Compare: ByteExact, OutputFiles: []string{".fastq"},
			OurArgs: args, UpstreamArgs: up,
		}
	}
	entries := []Entry{
		se("base"),
		se("q30", "-q", "30"),
		se("l30", "-l", "30"),
		se("no5prime", "-x"),
		se("truncn", "-n"),
		se("q30_l40", "-q", "30", "-l", "40"),
	}

	// PE: our pe and upstream pe both write -o/-p/-s; compare all three files.
	peOut := []string{".r1.fastq", ".r2.fastq", ".s.fastq"}
	peOurs := []string{"pe", "-f", r1, "-r", r2, "-t", "sanger", "-w", "0",
		"-o", "{out}.r1.fastq", "-p", "{out}.r2.fastq", "-s", "{out}.s.fastq"}
	peUp := []string{"pe", "-f", r1, "-r", r2, "-t", "sanger",
		"-o", "{out}.r1.fastq", "-p", "{out}.r2.fastq", "-s", "{out}.s.fastq"}
	entries = append(entries, Entry{
		Tool: "sickle", UpstreamTool: "sickle", Name: "sickle_pe_base",
		Input: InputFASTQPaired, Compare: ByteExact, OutputFiles: peOut,
		OurArgs: peOurs, UpstreamArgs: peUp, Heavy: true,
	})

	// CLI-default window (no -w): our sickle now defaults the window size to 0,
	// which selects the library's upstream-faithful dynamic int(0.1*len) window
	// (was hardcoded to 10, trimming ~1% of reads one window short). The
	// no-flag invocation is therefore byte-exact with upstream now (fixed by the
	// sickle agent; previously a documented Skip).
	entries = append(entries, Entry{
		Tool: "sickle", UpstreamTool: "sickle", Name: "sickle_se_cli_default_window",
		Input: InputFASTQ, Compare: ByteExact, OutputFiles: []string{".fastq"},
		OurArgs:      []string{"se", "-f", fq, "-t", "sanger", "-o", "{out}.fastq"},
		UpstreamArgs: []string{"se", "-f", fq, "-t", "sanger", "-o", "{out}.fastq"},
	})
	return entries
}

// skewerMatrix exercises skewer SE adapter trimming. Our skewer is
// subcommand-based (`skewer se -i IN -o OUT -x ADAPTER`); upstream is flat
// (`skewer -x ADAPTER -1 IN` with -1/--stdout writing to stdout). Because the
// two CLIs differ in shape we supply OurArgs / UpstreamArgs and compare stdout
// byte-for-byte (the trimmed FASTQ). The adapter is given explicitly so no
// heuristic auto-detection is involved and the output is deterministic.
func skewerMatrix() []Entry {
	fq := "{fastq}"
	adapter := "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	se := func(name string, ourExtra, upExtra []string) Entry {
		ours := append([]string{"se", "-i", fq, "-o", "-", "-x", adapter, "-t", "sanger"}, ourExtra...)
		up := append([]string{"-x", adapter, "-1", "-f", "sanger", fq}, upExtra...)
		return Entry{
			Tool: "skewer", UpstreamTool: "skewer", Name: "skewer_se_" + name,
			Input: InputFASTQ, Compare: ByteExact,
			OurArgs: ours, UpstreamArgs: up,
		}
	}
	base := se("base", nil, nil)
	minlen := se("minlen30", []string{"-l", "30"}, []string{"-l", "30"})
	heavy := se("full_heavy", nil, nil)
	heavy.Heavy = true
	return []Entry{base, minlen, heavy}
}

// fastpMatrix covers fastp ours-vs-upstream. The deterministic trimming/
// filtering paths are now byte-exact and run as ByteExact entries:
//
//   - default quality/length filtering (cut_tail and the no-flag default) match
//     upstream byte-for-byte after a spurious standalone end-quality-trim — one
//     that ran by default because qualified_quality_phred (-q) was mistakenly
//     treated as a trim threshold rather than a filter threshold — was removed.
//
// Only genuinely non-deterministic / non-comparable paths remain Skipped:
//
//   - adapter auto-detection (--detect_adapter_for_pe) is a sampling heuristic;
//     the per-tool suite validates it with a documented similarity bound.
//   - the --json/--html reports carry a version stamp and wall-clock time.
//
// Nothing here DIVERGEs.
func fastpMatrix() []Entry {
	fq := "{fastq}"
	r1, r2 := "{fastq1}", "{fastq2}"
	return []Entry{
		{
			// cut_tail is the cut_front/cut_tail/cut_right sliding-window trim.
			// It is now byte-exact against upstream: a spurious standalone
			// end-quality-trim that used to run by default (qualified_quality_phred
			// is a filter, not a trim threshold) was removed, so the sliding
			// window now sees the full read exactly as upstream does.
			Tool: "fastp", UpstreamTool: "fastp", Name: "fastp_cut_tail",
			Input: InputFASTQ, Compare: ByteExact, OutputFiles: []string{".fastq"},
			OurArgs:      []string{"-i", fq, "-o", "{out}.fastq", "-A", "--cut-tail", "--json", "{out}.json", "--html", "{out}.html"},
			UpstreamArgs: []string{"-i", fq, "-o", "{out}.fastq", "-A", "--cut_tail", "--json", "{out}.json", "--html", "{out}.html"},
		},
		{
			// Default quality/length filtering. Now byte-exact: with the spurious
			// default end-quality-trim removed, too-many-N reads keep their
			// N-laden tails and are dropped by the N filter exactly as upstream
			// does (previously we trimmed those tails so the reads slipped
			// through). --disable_quality_filtering/-Q and
			// --disable_length_filtering/-L are also now available.
			Tool: "fastp", UpstreamTool: "fastp", Name: "fastp_default_filter",
			Input: InputFASTQ, Compare: ByteExact, OutputFiles: []string{".fastq"},
			OurArgs:      []string{"-i", fq, "-o", "{out}.fastq", "-A", "--json", "{out}.json", "--html", "{out}.html"},
			UpstreamArgs: []string{"-i", fq, "-o", "{out}.fastq", "-A", "--json", "{out}.json", "--html", "{out}.html"},
		},
		{
			Tool: "fastp", UpstreamTool: "fastp", Name: "fastp_detect_adapter_pe_heavy",
			Input: InputFASTQPaired, Compare: Similarity, Heavy: true,
			OutputFiles:  []string{".r1.fastq", ".r2.fastq"},
			OurArgs:      []string{"-I", r1, "--in2", r2, "-O", "{out}.r1.fastq", "--out2", "{out}.r2.fastq", "--detect_adapter_for_pe", "--json", "{out}.json", "--html", "{out}.html"},
			UpstreamArgs: []string{"-i", r1, "-I", r2, "-o", "{out}.r1.fastq", "-O", "{out}.r2.fastq", "--detect_adapter_for_pe", "--json", "{out}.json", "--html", "{out}.html"},
			Skip: "fastp adapter auto-detection is a sampling heuristic; the per-tool suite validates it with a documented " +
				"similarity bound (TestUnitDetectAdapterSE). The two CLIs also differ in PE input flags. Owned by the fastp agent.",
		},
	}
}
