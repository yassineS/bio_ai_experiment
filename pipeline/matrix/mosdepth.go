package matrix

// This file registers the mosdepth matrix. mosdepth writes several
// "<prefix>.mosdepth.*" / "<prefix>.*.bed.gz" output files (no stdout), so
// these entries use the output-file comparison path. mosdepth's CLI is
// "mosdepth [options] <prefix> <in.bam>"; <prefix> is positional, so the {out}
// placeholder resolves it directly and a shared Args serves both sides.
//
// Upstream mosdepth ships ONLY as a linux/amd64 GitHub release asset (it is a
// Nim project, not built from source here). The runner resolves it from
// $MOSDEPTH_BIN or the temp-dir release cache the per-tool parity tests
// populate; on any other platform, or when neither is present, the entry
// reports ERROR with an actionable hint. The matrix marks the whole set with a
// platform Skip on non-linux/amd64 so it mirrors the existing per-tool gate.
//
// Two real parity gaps in our mosdepth --by path are documented as Skips rather
// than papered over: (1) the summary file omits upstream's per-region
// "<chrom>_region" / "total_region" rows, and (2) we do not emit the
// "<prefix>.mosdepth.region.dist.txt" file. The substantive --by output (the
// regions.bed.gz depths) IS byte-exact and is exercised directly.

import "github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"

func init() {
	Register(mosdepthMatrix()...)
}

func mosdepthMatrix() []Entry {
	bam := "{bam}"
	bed := "{bed}"

	// platformSkip is set on every entry when this platform has no upstream
	// mosdepth release binary, mirroring the per-tool parity gate.
	platformSkip := ""
	if !upstream.MosdepthSupported() {
		platformSkip = "upstream mosdepth release binary is only published for linux/amd64; skipping on this platform (mirrors the per-tool mosdepth parity gate)"
	}

	mk := func(name string, out []string, args ...string) Entry {
		full := append(args, "{out}", bam)
		return Entry{
			Tool: "mosdepth", UpstreamTool: "mosdepth", Name: "mosdepth_" + name,
			Input: InputBAM, Compare: ByteExact, OutputFiles: out,
			Args: full, Skip: platformSkip,
		}
	}

	// Fully byte-exact modes: per-base depth + summary + global dist.
	perBaseOut := []string{".mosdepth.summary.txt", ".mosdepth.global.dist.txt", ".per-base.bed.gz"}
	entries := []Entry{
		mk("default", perBaseOut),
		mk("fast_mode", perBaseOut, "--fast-mode"),
		mk("mapq20", perBaseOut, "--mapq", "20"),
		mk("flag", perBaseOut, "--flag", "1796"),
	}

	// --by: the regions.bed.gz per-region mean is byte-exact for both fixed
	// windows and BED-defined regions. The per-region mean is now computed as
	// upstream mosdepth's imean does — Σ(depth_i / L) accumulated per base, NOT
	// (Σ depth_i)/L — so the float rounding matches and the %.2f column agrees
	// (the prior ±0.01 boundary divergence on ~4/1240 regions is closed).
	regionsOnly := []string{".regions.bed.gz", ".mosdepth.global.dist.txt"}
	byBedThresh := func() Entry {
		e := mk("by_bed_thresholds", []string{".thresholds.bed.gz", ".regions.bed.gz"}, "--by", bed, "--thresholds", "1,5,10")
		e.Heavy = true
		return e
	}()
	entries = append(entries,
		mk("by_bed_regions", regionsOnly, "--by", bed),
		mk("by_window_regions", regionsOnly, "--by", "500"),
		byBedThresh,
	)

	// Heavy default run for the timing ratio.
	heavy := mk("default_heavy", perBaseOut)
	heavy.Heavy = true
	entries = append(entries, heavy)

	// Previously-documented --by parity gaps, now FIXED and re-activated: our
	// mosdepth --by summary now emits upstream's per-region "<chrom>_region" /
	// "total_region" rows, and we now write "<prefix>.mosdepth.region.dist.txt"
	// (region depth-distribution). Both are byte-exact against upstream and are
	// exercised directly here (subject only to the platform gate).
	entries = append(entries,
		mk("by_summary_region_rows", []string{".mosdepth.summary.txt"}, "--by", bed),
	)
	rd := Entry{
		Tool: "mosdepth", UpstreamTool: "mosdepth", Name: "mosdepth_by_region_dist",
		Input: InputBAM, Compare: ByteExact, OutputFiles: []string{".mosdepth.region.dist.txt"},
		Args: []string{"--by", bed, "{out}", bam}, Skip: platformSkip,
	}
	entries = append(entries, rd)

	return entries
}
