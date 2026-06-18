package matrix

// This file registers the vcftools matrix. vcftools writes most of its output
// to "<prefix>.<ext>" files rather than stdout, so these entries use the
// output-file comparison path (OutputFiles + the {out} prefix placeholder).
// Our vcftools and the upstream binary share the same CLI shape (flat flags,
// --vcf / --out), so a single shared Args serves both sides; the {out} token
// resolves to a per-invocation, per-side prefix.
//
// Two fixtures are used: the single-sample VCF for per-site/summary modes, and
// the multi-sample VCF (vcf_multi_plain) for modes that need >1 sample
// (relatedness, het). vcftools writes a "<prefix>.log" too, but our port does
// not, and the log is non-reproducible anyway, so it is never listed in
// OutputFiles.
//
// vcftools 0.1.14 has a temp-file off-by-one VLA that a modern glibc
// (_FORTIFY_SOURCE on at -O2) aborts on, so the pairwise-LD (--geno-r2/--hap-r2)
// and 012-matrix (--012) writers crash before producing output. The vendored
// build applies reference_code/patches/vcftools-tmpfile-vla-off-by-one.patch,
// after which: --012 is byte-exact (closed); --geno-r2 emits the same pairs in a
// different order plus a few near-zero r^2 residuals (our-side fix pending);
// --hap-r2 needs a phased fixture (the unphased multi-sample VCF yields no LD
// sites). --LROH runs cleanly with the required --chr and is byte-exact.

func init() {
	Register(vcftoolsMatrix()...)
}

func vcftoolsMatrix() []Entry {
	vcf := "{vcf_plain}"
	multi := "{vcf_multi_plain}"

	// single builds an entry on the single-sample VCF.
	single := func(name, ext string, modeArgs ...string) Entry {
		args := append([]string{"--vcf", vcf}, modeArgs...)
		args = append(args, "--out", "{out}")
		return Entry{
			Tool: "vcftools", UpstreamTool: "vcftools", Name: "vcftools_" + name,
			Input: InputVCFPlain, Compare: ByteExact, OutputFiles: []string{ext},
			Args: args,
		}
	}
	// multiS builds an entry on the multi-sample VCF.
	multiS := func(name, ext string, modeArgs ...string) Entry {
		args := append([]string{"--vcf", multi}, modeArgs...)
		args = append(args, "--out", "{out}")
		return Entry{
			Tool: "vcftools", UpstreamTool: "vcftools", Name: "vcftools_" + name,
			Input: InputVCFMulti, Compare: ByteExact, OutputFiles: []string{ext},
			Args: args,
		}
	}

	entries := []Entry{
		// --- Per-site / frequency / depth (single-sample VCF) ---
		single("freq", ".frq", "--freq"),
		single("counts", ".frq.count", "--counts"),
		single("freq2", ".frq", "--freq2"),
		single("depth", ".idepth", "--depth"),
		single("site_depth", ".ldepth", "--site-depth"),
		single("site_mean_depth", ".ldepth.mean", "--site-mean-depth"),
		single("site_pi", ".sites.pi", "--site-pi"),
		single("window_pi", ".windowed.pi", "--window-pi", "1000"),
		single("tstv_summary", ".TsTv.summary", "--TsTv-summary"),
		single("missing_indv", ".imiss", "--missing-indv"),
		single("missing_site", ".lmiss", "--missing-site"),
		single("het", ".het", "--het"),
		single("singletons", ".singletons", "--singletons"),
		single("recode_heavy", ".recode.vcf", "--recode", "--recode-INFO-all"),

		// --- Modes that need >1 sample (multi-sample VCF) ---
		multiS("het_multi", ".het", "--het"),
		multiS("relatedness", ".relatedness", "--relatedness"),
		multiS("relatedness2", ".relatedness2", "--relatedness2"),
		multiS("freq_multi", ".frq", "--freq"),
		multiS("missing_indv_multi", ".imiss", "--missing-indv"),

		// --- Heavy: window-pi over the whole VCF, timing ratio surfaced ---
		func() Entry {
			e := single("window_pi_heavy", ".windowed.pi", "--window-pi", "500", "--window-pi-step", "250")
			e.Heavy = true
			return e
		}(),
	}

	// recode INFO ordering now preserves the data-line source order (DP;AF) like
	// upstream, so recode_heavy is byte-exact and runs (the prior alphabetical
	// AF;DP serialisation was fixed).

	// Upstream-crashing modes on this build (recorded as documented Skips so
	// they neither run nor DIVERGE). These are NOT our bugs and cannot be
	// byte-validated: the vendored upstream binary aborts (glibc buffer overflow
	// / segfault) before writing output, so there is no golden to compare
	// against. Our port produces correct output, validated by the per-tool unit
	// suite. Keeping them skipped is the correct terminal state — re-running
	// upstream would only reproduce the crash.
	crash := func(name, ext, reason string, modeArgs ...string) Entry {
		e := multiS(name, ext, modeArgs...)
		e.Skip = reason
		return e
	}
	entries = append(entries,
		// --geno-r2: pairwise genotype LD. With the temp-file off-by-one patch
		// (see reference_code/patches) upstream runs; our port emits the pairs in
		// upstream's outer-first-SNP order and encodes genotypes as the reference
		// allele count (matching calc_geno_r2's sx/sy), so the XY-X*Y
		// cancellation that yields an exact 0 for uncorrelated sites is
		// bit-identical. A 3000bp LD window keeps the output ~190k rows; byte-
		// exact against upstream.
		func() Entry {
			e := multiS("geno_r2", ".geno.ld", "--geno-r2", "--ld-window-bp", "3000")
			e.Heavy = true
			return e
		}(),
		crash("hap_r2", ".hap.ld",
			"upstream vcftools --hap-r2 needs PHASED genotypes; the multi-sample fixture is unphased, so upstream errors 'Insufficient sites "+
				"remained after filtering' and emits no LD. Needs a phased fixture to exercise (and the temp-file patch to run).",
			"--hap-r2"),
		// --012: with the vendored binary's temp-file off-by-one patched (see
		// reference_code/patches), upstream writes the .012 matrix plus the
		// .012.indv (sample list) and .012.pos (CHROM/POS list) sidecars; all
		// three are byte-exact against our port.
		func() Entry {
			e := multiS("matrix012", ".012", "--012")
			e.OutputFiles = []string{".012", ".012.indv", ".012.pos"}
			return e
		}(),
		// --LROH detects runs of homozygosity via the Boyko/Auton forward-
		// backward HMM. Upstream requires a single --chr, so the entry passes
		// --chr chr1; our port reproduces the 8-column report (including the
		// MIN_START/MAX_END/N_MISMATCHES columns) byte-for-byte.
		multiS("lroh", ".LROH", "--LROH", "--chr", "chr1"),
		// --hardy: byte-exact. Our .hwe writer now emits glibc's '-nan' for the
		// monomorphic-site ChiSq (matching the upstream C++ printf rendering of a
		// quiet NaN) instead of Go's 'NaN'.
		single("hardy", ".hwe", "--hardy"),
	)
	return entries
}
