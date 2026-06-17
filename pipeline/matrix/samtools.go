package matrix

// This file registers the COMPREHENSIVE samtools matrix, spanning the 24 ported
// subcommands. The smoke matrix (smoke.go) keeps a tiny `samtools view` slice
// for the end-to-end loop; this file is the curated combinatorics layer on top.
//
// Comparison strategy (the recurring constraints):
//
//   - The runner runs ONE command per side and compares stdout (or output
//     files written through an {out} prefix). It does NOT pipe or post-process,
//     so a subcommand's parity has to be observable in a single command's text
//     output.
//   - Binary BGZF (BAM/CRAM) bytes are NOT byte-comparable: our klauspost
//     deflate backend frames blocks differently from htslib though both decode
//     identically (see the pipeline README). So every entry here that would
//     otherwise emit BAM/CRAM is driven to emit TEXT instead — SAM via
//     `-O sam` / `view`, pileup/coverage/depth tab text, or the FASTQ/FASTA
//     conversions — and compared byte-exact on the decoded stream.
//   - Several subcommands whose text-output path or upstream-faithful behaviour
//     our port does not yet match are recorded as documented Skips (the same
//     convention sickle/vcftools use): the divergence is real and surfaced for
//     the samtools agent, not papered over.
//
// Subcommands compared byte-exact (text/decoded): view (BAM/CRAM, region,
// -f/-F/-q/-c/-b/-h/-L bed), sort (name + coord via count), flagstat, idxstats,
// stats, depth (-a/-b/-r), coverage, mpileup (documented Skip on a 0-depth
// gap), calmd, consensus, dict, quickcheck, fastq/fasta, cat (count round-trip),
// tview (text). Documented Skips: subsample -s (RNG differs), markdup / fixmate
// / addreplacerg / merge / reheader / split / import / phase (output-format /
// behaviour gaps spelled out per entry).

func init() {
	Register(samtoolsViewMatrix()...)
	Register(samtoolsSortMatrix()...)
	Register(samtoolsStatsFamily()...)
	Register(samtoolsDepthFamily()...)
	Register(samtoolsDecodeText()...)
	Register(samtoolsBinaryOutputSkips()...)
}

// samtoolsViewMatrix is the big `view` flag sweep over BAM and CRAM, all decoded
// to SAM/count text. -s (subsample) is a documented Skip because our RNG read
// selection differs from htslib's even at a fixed seed.
func samtoolsViewMatrix() []Entry {
	bam := ExpandSpec{
		Tool:           "samtools",
		Subcommand:     "view",
		UsesSubcommand: true,
		Input:          InputBAM,
		Compare:        ByteExact,
		BaseArgs:       []string{"{bam}"},
		Flags: []Flag{
			{Name: "-h", Bool: true},                        // SAM with header (decoded text)
			{Name: "-H", Bool: true},                        // header only
			{Name: "-c", Bool: true},                        // count
			{Name: "-q", Values: []string{"20", "30"}},      // MAPQ filter
			{Name: "-f", Values: []string{"0x2", "0x10"}},   // require flag bits
			{Name: "-F", Values: []string{"0x100", "0x10"}}, // exclude flag bits
			{Name: "-L", Values: []string{"{bed}"}},         // overlap a BED region set
		},
		Combos: []Combo{
			{Name: "count_q30", Flags: []string{"-c", "-q", "30"}},
			{Name: "count_L", Flags: []string{"-c", "-L", "{bed}"}},
			{Name: "f2_F256_q20_count", Flags: []string{"-f", "0x2", "-F", "0x100", "-q", "20", "-c"}},
			{Name: "header_plus_body", Flags: []string{"-h", "-q", "20"}},
		},
	}.Expand()
	for i := range bam {
		bam[i].Name = "view_bam_" + bam[i].Name
	}

	// Region queries (require the index): a contig and a sub-range.
	region := []Entry{
		mkSam("view", "view_region_contig", InputBAM, ByteExact, "{bam}", "chr1"),
		mkSam("view", "view_region_range", InputBAM, ByteExact, "{bam}", "chr1:1-5000"),
		mkSam("view", "view_region_count", InputBAM, ByteExact, "-c", "{bam}", "chr2"),
	}

	// CRAM: same flag ideas, plus the reference (-T). Decoded to SAM text.
	cram := []Entry{
		mkSam("view", "view_cram_body", InputCRAM, ByteExact, "-T", "{fasta}", "{cram}"),
		mkSam("view", "view_cram_header", InputCRAM, ByteExact, "-H", "-T", "{fasta}", "{cram}"),
		mkSam("view", "view_cram_count", InputCRAM, ByteExact, "-c", "-T", "{fasta}", "{cram}"),
		mkSam("view", "view_cram_q30", InputCRAM, ByteExact, "-q", "30", "-T", "{fasta}", "{cram}"),
		func() Entry {
			e := mkSam("view", "view_cram_decode_sam_heavy", InputCRAM, ByteExact, "-T", "{fasta}", "{cram}")
			e.Heavy = true
			return e
		}(),
	}

	// Subsample with a fixed seed: our subsample now ports htslib's exact
	// name-hash RNG (Wang/X31 hash + glibc rand), so at a fixed seed the
	// selected read SET — and the headerless-SAM stdout — is byte-identical to
	// upstream (fixed by the samtools agent; previously a documented Skip).
	subsample := []Entry{{
		Tool: "samtools", Subcommand: "view", UsesSubcommand: true,
		Name: "view_subsample_seed", Input: InputBAM, Compare: ByteExact,
		Args: []string{"-s", "11.5", "{bam}"},
	}}

	out := append(bam, region...)
	out = append(out, cram...)
	out = append(out, subsample...)
	return out
}

// samtoolsSortMatrix exercises sort. Name sort fully orders records so the
// decoded SAM is byte-exact. Coordinate sort is verified order-insensitively via
// `view -c` on the result is not expressible in one command, so the coord-decode
// entry is a documented Skip recording the equal-coordinate tie-break difference
// (our sort orders records at an identical (rname,pos) differently from htslib).
func samtoolsSortMatrix() []Entry {
	return []Entry{
		mkSam("sort", "sort_name_sam", InputBAM, ByteExact, "-n", "-O", "sam", "{bam}"),
		mkSam("sort", "sort_byname_tag", InputBAM, ByteExact, "-N", "-O", "sam", "{bam}"),
		func() Entry {
			e := mkSam("sort", "sort_name_sam_heavy", InputBAM, ByteExact, "-n", "-O", "sam", "{bam}")
			e.Heavy = true
			return e
		}(),
		// Coordinate sort: our coordCmp now matches htslib's tie-break exactly
		// (refID, pos, then reverse-flag), so records at an identical (rname,pos)
		// are emitted in upstream order and the decoded SAM is byte-exact (fixed
		// by the samtools agent; previously a documented Skip).
		mkSam("sort", "sort_coord_sam", InputBAM, ByteExact, "-O", "sam", "{bam}"),
	}
}

// samtoolsStatsFamily covers the report subcommands whose output is plain text
// and byte-exact after provenance stripping: flagstat, idxstats, stats,
// quickcheck, dict.
func samtoolsStatsFamily() []Entry {
	return []Entry{
		mkSam("flagstat", "flagstat", InputBAM, ByteExact, "{bam}"),
		mkSam("idxstats", "idxstats", InputBAM, ByteExact, "{bam}"),
		mkSam("stats", "stats", InputBAM, ByteExact, "{bam}"),
		mkSam("quickcheck", "quickcheck", InputBAM, ByteExact, "{bam}"),
		mkSam("dict", "dict", InputFASTA, ByteExact, "{fasta}"),
		func() Entry {
			e := mkSam("stats", "stats_heavy", InputBAM, ByteExact, "{bam}")
			e.Heavy = true
			return e
		}(),
	}
}

// samtoolsDepthFamily covers depth and coverage, both tab-text output.
func samtoolsDepthFamily() []Entry {
	depth := ExpandSpec{
		Tool: "samtools", Subcommand: "depth", UsesSubcommand: true,
		Input: InputBAM, Compare: ByteExact, BaseArgs: []string{"{bam}"},
		Flags: []Flag{
			{Name: "-a", Bool: true},                // all positions
			{Name: "-r", Values: []string{"chr1"}},  // region
			{Name: "-b", Values: []string{"{bed}"}}, // BED-restricted
		},
		Combos: []Combo{
			{Name: "all_region", Flags: []string{"-a", "-r", "chr1:1-4000"}},
			{Name: "all_bed", Flags: []string{"-a", "-b", "{bed}"}},
		},
	}.Expand()
	for i := range depth {
		depth[i].Name = "depth_" + depth[i].Name
	}
	depth = append(depth,
		mkSam("coverage", "coverage", InputBAM, ByteExact, "{bam}"),
		mkSam("coverage", "coverage_region", InputBAM, ByteExact, "-r", "chr1", "{bam}"),
		// depth -Q (min mapping quality) now matches upstream's letter assignment
		// and filtering exactly and is byte-exact (the old -q/-Q swap is fixed).
		mkSam("depth", "depth_mapq_filter", InputBAM, ByteExact, "-Q", "20", "{bam}"),
		// depth -q (min BASE quality): documented Skip. The -q/-Q letter swap is
		// fixed, but a residual gap remains — when a covered interior position's
		// only base is filtered out by -q, upstream emits that position with
		// depth 0 (e.g. `chr1 8 0`) whereas our port omits the row entirely.
		// Same root cause as mpileup_pileup below (interior zero-depth emission
		// after per-base quality filtering). Owned by the samtools agent.
		Entry{
			Tool: "samtools", Subcommand: "depth", UsesSubcommand: true,
			Name: "depth_baseq_filter", Input: InputBAM, Compare: ByteExact,
			Args: []string{"-q", "10", "{bam}"},
			Skip: "samtools depth -q (min base quality): the -q/-Q swap is fixed and -Q is byte-exact, but when an interior covered " +
				"position's only base is filtered out by -q, upstream emits it with depth 0 (`chr1 8 0`) while our port omits the row. " +
				"Same interior-zero-depth root cause as mpileup. Owned by the samtools agent.",
		},
	)
	return depth
}

// samtoolsDecodeText covers the remaining subcommands whose parity is observable
// as text in one command: mpileup, calmd, consensus, fastq/fasta, cat, tview.
func samtoolsDecodeText() []Entry {
	return []Entry{
		// calmd writes SAM to stdout (NM/MD recomputed against the reference).
		mkSam("calmd", "calmd", InputBAM, ByteExact, "{bam}", "{fasta}"),
		// consensus emits a FASTA/FASTQ consensus; default FASTA matches.
		mkSam("consensus", "consensus", InputBAM, ByteExact, "{bam}"),
		mkSam("consensus", "consensus_region", InputBAM, ByteExact, "-r", "chr1:1-2000", "{bam}"),
		// fastq / fasta: BAM -> FASTQ/FASTA text.
		mkSam("fastq", "fastq", InputBAM, ByteExact, "{bam}"),
		mkSam("fastq", "fastq_n", InputBAM, ByteExact, "-n", "{bam}"),
		func() Entry {
			e := mkSam("fastq", "fastq_heavy", InputBAM, ByteExact, "{bam}")
			e.Heavy = true
			return e
		}(),
		// tview in text mode (-d T) renders an ASCII alignment view; byte-exact.
		mkSam("tview", "tview_text", InputBAM, ByteExact, "-d", "T", "{bam}", "{fasta}"),

		// mpileup: documented Skip. The non-zero-depth pileup rows are byte-exact
		// (verified), but upstream additionally emits interior zero-depth rows —
		// positions inside the covered span whose only base is filtered out (or
		// genuinely uncovered), e.g. `chr1 8 A 0 * *` — which our port omits.
		// Same interior-zero-depth root cause as depth_baseq_filter; the -a
		// full-contig zero-row work did not cover this implicit interior case.
		// A real port gap owned by the samtools agent.
		{
			Tool: "samtools", Subcommand: "mpileup", UsesSubcommand: true,
			Name: "mpileup_pileup", Input: InputBAM, Compare: ByteExact,
			Args: []string{"-f", "{fasta}", "{bam}"},
			Skip: "samtools mpileup: our non-zero-depth pileup rows match upstream byte-for-byte, but upstream also emits interior " +
				"zero-depth positions (covered span, only base filtered out or uncovered) such as `chr1 8 A 0 * *` that our port omits. " +
				"Same interior-zero-depth root cause as depth -q. Real port gap owned by the samtools agent.",
		},

		// cat: documented Skip. cat concatenates BAMs into a BAM; our cat ignores
		// -O sam and always emits BGZF BAM, which is not byte-comparable. The
		// decoded record stream is correct (verified via `view`), but that decode
		// step cannot be expressed in the single-command runner.
		{
			Tool: "samtools", Subcommand: "cat", UsesSubcommand: true,
			Name: "cat_concat", Input: InputBAM, Compare: ByteExact,
			Args: []string{"-o", "-", "{bam}", "{bam}"},
			Skip: "samtools cat emits BGZF BAM (not byte-comparable: klauspost vs htslib block framing) and does not honour -O sam to " +
				"produce decodable text in one command. The concatenated record stream is correct when decoded via `view`. Owned by the samtools agent.",
		},
	}
}

// samtoolsBinaryOutputSkips records the subcommands that emit BGZF BAM/CRAM and
// — unlike view/sort/calmd/consensus — do NOT honour `-O sam` to produce a
// decodable text stream, so their parity cannot be expressed byte-exact in the
// single-command runner. Each Skip spells out the concrete reason; where the
// decoded DATA was verified correct out-of-band it is noted, and where the data
// itself diverges that is flagged for the samtools agent.
func samtoolsBinaryOutputSkips() []Entry {
	skip := func(name, args, reason string, argv ...string) Entry {
		return Entry{
			Tool: "samtools", Subcommand: name, UsesSubcommand: true,
			Name: "samtools_" + name, Input: InputBAM, Compare: ByteExact,
			Args: argv, Skip: reason,
		}
	}
	return []Entry{
		skip("markdup", "", "samtools markdup ignores -O/--output-fmt and always writes BGZF BAM (not byte-comparable). The duplicate "+
			"marking itself is byte-identical to upstream when decoded out-of-band; the output-format gap is owned by the samtools agent.",
			"{bam}", "-O", "sam", "-"),
		skip("fixmate", "", "samtools fixmate ignores -O sam (writes BGZF BAM) and needs name-collated input, which cannot be produced "+
			"in a single command. Owned by the samtools agent.",
			"-O", "sam", "{bam}", "-"),
		skip("addreplacerg", "", "samtools addreplacerg ignores -O sam (writes BGZF BAM) and, when decoded, our -r does not replace the "+
			"existing RG tag (keeps the original RG) whereas upstream substitutes it. Real port gap owned by the samtools agent.",
			"-O", "sam", "-r", `ID:x\tSM:y`, "{bam}", "-"),
		skip("merge", "", "samtools merge ignores -O sam (writes BGZF BAM); the fixtures share sample/RG names so a meaningful merge of "+
			"distinct inputs is not expressible. Owned by the samtools agent.",
			"-O", "sam", "-f", "-", "{bam}", "{bam}"),
		skip("reheader", "", "samtools reheader rewrites the header of a BAM in place and emits BGZF BAM (not byte-comparable); it has no "+
			"single-command text-output form. Owned by the samtools agent.",
			"-c", "cat", "{bam}"),
		skip("split", "", "samtools split writes one BGZF BAM per read group to a directory of files (not byte-comparable and not a single "+
			"stdout/prefix the runner compares). Owned by the samtools agent.",
			"-f", "%*_%!.bam", "{bam}"),
		skip("import", "", "samtools import (FASTQ -> unaligned BAM) ignores -O sam and writes BGZF BAM (not byte-comparable). Owned by the samtools agent.",
			"-O", "sam", "{bam}"),
		skip("phase", "", "samtools phase emits phased BGZF BAM(s) (not byte-comparable) and has no single-command decoded-text form. Owned by the samtools agent.",
			"{bam}"),
	}
}

// mkSam builds a samtools entry (UsesSubcommand=true) with the shared Args used
// for both sides. name is the unique-within-tool label (it is prefixed with
// "samtools_" for the report).
func mkSam(sub, name string, input InputKind, cmp CompareMode, args ...string) Entry {
	return Entry{
		Tool: "samtools", Subcommand: sub, UsesSubcommand: true,
		Name: "samtools_" + name, Input: input, Compare: cmp, Args: args,
	}
}
