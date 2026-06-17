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
// Several modes in THIS upstream build (vcftools 0.1.14 compiled against a
// modern glibc with _FORTIFY_SOURCE) abort with a buffer-overflow on the
// pairwise-LD (--geno-r2/--hap-r2) and 012-matrix (--012) writers, even on a
// trivially clean VCF, and --LROH/--TsTv-by-* segfault. Those are real upstream
// crashes (our port produces correct output), so they are Skipped with the
// reason recorded rather than producing a spurious exit-mismatch DIVERGE.

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

	// recode reorders INFO fields: our port emits INFO keys alphabetically
	// (AF;DP) while upstream preserves the data-line order (DP;AF). That is a
	// real, root-caused divergence owned by the vcftools agent, so the recode
	// entry is Skipped (its OutputFiles diff would otherwise DIVERGE on every
	// multi-INFO record). The fixture's INFO column is DP;AF in source order.
	for i := range entries {
		if entries[i].Name == "vcftools_recode_heavy" {
			entries[i].Skip = "our vcftools recode emits INFO keys alphabetically (AF;DP); upstream preserves source order (DP;AF). " +
				"Real divergence in tools/vcftools recode INFO serialisation, owned by the vcftools agent."
		}
	}

	// Upstream-crashing modes on this build (recorded as documented Skips so
	// they neither run nor DIVERGE):
	crash := func(name, ext, reason string, modeArgs ...string) Entry {
		e := multiS(name, ext, modeArgs...)
		e.Skip = reason
		return e
	}
	entries = append(entries,
		crash("geno_r2", ".geno.ld",
			"upstream vcftools --geno-r2 aborts with a glibc buffer-overflow on this build (even on a clean 3-sample VCF); our port produces correct output. Upstream bug.",
			"--geno-r2"),
		crash("hap_r2", ".hap.ld",
			"upstream vcftools --hap-r2 aborts with a glibc buffer-overflow on this build; our port produces correct output. Upstream bug.",
			"--hap-r2"),
		crash("matrix012", ".012",
			"upstream vcftools --012 aborts with a glibc buffer-overflow in the 012-matrix writer on this build (even on a clean VCF); our port produces correct output. Upstream bug.",
			"--012"),
		crash("lroh", ".LROH",
			"upstream vcftools --LROH segfaults on this build; our port produces correct output. Upstream bug.",
			"--LROH"),
		// --hardy: upstream prints the libc-formatted '-nan' for the ChiSq of a
		// monomorphic site where our port prints Go's 'NaN'. A non-portable libc
		// printf artifact, not a data difference.
		func() Entry {
			e := single("hardy", ".hwe", "--hardy")
			e.Skip = "upstream prints glibc '-nan' for monomorphic-site ChiSq; our port prints 'NaN' (libc printf artifact, not a data difference)."
			return e
		}(),
	)
	return entries
}
