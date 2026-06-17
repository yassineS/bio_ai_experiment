package fastp

// Live-parity tests for the JSON report sub-fields added to close the second
// fastp residual: the per-read quality_curves / content_curves / kmer_count /
// q40_bases blocks (for BOTH the before- and after-filtering streams), the
// real per-read q20_bases / q30_bases totals, the summary.sequencing string,
// and the paired-end insert_size block.
//
// These assert against the upstream OpenGene/fastp binary (built/located by the
// shared ensureUpstream helper). Per the project's deterministic-vs-heuristic
// principle:
//
//   - DETERMINISTIC counts (q20/q30/q40 bases, total_cycles, kmer_count, the
//     insert_size histogram/peak/unknown, summary.sequencing) are compared
//     EXACTLY.
//   - The quality_curves / content_curves are deterministic but upstream emits
//     them through C++ ostream's default 6-significant-digit float formatting,
//     so an exact float64 round-trip is impossible. They are validated with a
//     tight absolute tolerance (1e-4) — a structural/numeric-equality check, not
//     a heuristic similarity bound. Key sets and lengths must match exactly.
//
// Intentionally EXCLUDED from comparison (non-reproducible, as the existing
// JSON parity already excludes): summary.fastp_version (upstream's version
// string vs our ToolVersion), the top-level command string, and the Go-only
// tool.time wall-clock field. These are documented in tools/fastp/README.md.

import (
	"math"
	"testing"
)

// assertCurvesClose compares two curve maps (metric name -> per-cycle values)
// for identical key sets/lengths and per-value agreement within absTol.
func assertCurvesClose(t *testing.T, label string, upJSON, goJSON map[string]any, path string, absTol float64) {
	t.Helper()
	up, upOK := lookupJSON(upJSON, path)
	got, goOK := lookupJSON(goJSON, path)
	if !upOK || !goOK {
		t.Fatalf("%s: path %s missing (upstream=%v go=%v)", label, path, upOK, goOK)
	}
	um, ok1 := up.(map[string]any)
	gm, ok2 := got.(map[string]any)
	if !ok1 || !ok2 {
		t.Fatalf("%s: %s not an object (upstream=%T go=%T)", label, path, up, got)
	}
	if len(um) != len(gm) {
		t.Fatalf("%s: %s key count differs: upstream=%d go=%d", label, path, len(um), len(gm))
	}
	for k, uv := range um {
		gv, ok := gm[k]
		if !ok {
			t.Fatalf("%s: %s missing curve %q in go", label, path, k)
		}
		ul, ok1 := uv.([]any)
		gl, ok2 := gv.([]any)
		if !ok1 || !ok2 {
			t.Fatalf("%s: %s.%s not an array", label, path, k)
		}
		if len(ul) != len(gl) {
			t.Fatalf("%s: %s.%s length differs: upstream=%d go=%d", label, path, k, len(ul), len(gl))
		}
		for i := range ul {
			uf, _ := toFloat(ul[i])
			gf, _ := toFloat(gl[i])
			if math.Abs(uf-gf) > absTol {
				t.Fatalf("%s: %s.%s[%d] differs beyond %g: upstream=%v go=%v", label, path, k, i, absTol, uf, gf)
			}
		}
	}
}

// assertKmerExact compares the 1024-entry kmer_count histograms for exact
// equality (deterministic integer counts).
func assertKmerExact(t *testing.T, label string, upJSON, goJSON map[string]any, path string) {
	t.Helper()
	up, upOK := lookupJSON(upJSON, path)
	got, goOK := lookupJSON(goJSON, path)
	if !upOK || !goOK {
		t.Fatalf("%s: path %s missing (upstream=%v go=%v)", label, path, upOK, goOK)
	}
	um, _ := up.(map[string]any)
	gm, _ := got.(map[string]any)
	if len(um) != 1024 || len(gm) != 1024 {
		t.Fatalf("%s: %s expected 1024 kmers, upstream=%d go=%d", label, path, len(um), len(gm))
	}
	for k, uv := range um {
		gv, ok := gm[k]
		if !ok {
			t.Fatalf("%s: %s missing kmer %q in go", label, path, k)
		}
		uf, _ := toFloat(uv)
		gf, _ := toFloat(gv)
		if uf != gf {
			t.Fatalf("%s: %s[%q] differs: upstream=%v go=%v", label, path, k, uf, gf)
		}
	}
}

// assertStringEqual compares a string-valued JSON field exactly.
func assertStringEqual(t *testing.T, label string, upJSON, goJSON map[string]any, path string) {
	t.Helper()
	up, upOK := lookupJSON(upJSON, path)
	got, goOK := lookupJSON(goJSON, path)
	if !upOK || !goOK {
		t.Fatalf("%s: path %s missing (upstream=%v go=%v)", label, path, upOK, goOK)
	}
	if up != got {
		t.Fatalf("%s: %s differs: upstream=%q go=%q", label, path, up, got)
	}
}

// readStatsSubfields is the per-read sub-field set asserted exact / close.
func assertReadStatsSubfields(t *testing.T, label string, upJSON, goJSON map[string]any, section string) {
	t.Helper()
	assertCounters(t, label, upJSON, goJSON,
		section+".total_reads",
		section+".total_bases",
		section+".q20_bases",
		section+".q30_bases",
		section+".q40_bases",
		section+".total_cycles",
	)
	assertCurvesClose(t, label, upJSON, goJSON, section+".quality_curves", 1e-4)
	assertCurvesClose(t, label, upJSON, goJSON, section+".content_curves", 1e-4)
	assertKmerExact(t, label, upJSON, goJSON, section+".kmer_count")
}

// TestParity_Fastp_JSONFields_SE validates the new SE JSON sub-fields against
// upstream for both the before- and after-filtering streams. The adapter +
// length filter make the after stream genuinely differ from the before stream,
// exercising the after-curve/kmer/q-count tracking.
func TestParity_Fastp_JSONFields_SE(t *testing.T) {
	bin := ensureUpstream(t)
	in := parityInput(t, "se_adapter.fq")
	dir := t.TempDir()

	_, upJSON := runUpstreamSE(t, bin, in, dir, []string{"-a", "AGATCGGAAGAGC", "-l", "50"})

	opts := DefaultProcessOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	opts.MinLength = 50
	opts.LengthRequired = 50
	_, goStats := runGoFastpSE(t, in, opts)
	goJSON := jsonFromStats(t, goStats)

	assertStringEqual(t, "SE sequencing", upJSON, goJSON, "summary.sequencing")
	assertReadStatsSubfields(t, "SE before", upJSON, goJSON, "read1_before_filtering")
	assertReadStatsSubfields(t, "SE after", upJSON, goJSON, "read1_after_filtering")
}

// TestParity_Fastp_JSONFields_PE validates the PE JSON sub-fields plus the
// paired-end-only insert_size block (peak, unknown, full histogram) against
// upstream.
func TestParity_Fastp_JSONFields_PE(t *testing.T) {
	bin := ensureUpstream(t)
	r1 := parityInput(t, "pe_r1.fq")
	r2 := parityInput(t, "pe_r2.fq")
	dir := t.TempDir()

	_, _, upJSON := runUpstreamPE(t, bin, r1, r2, dir, nil)

	opts := DefaultProcessOptions()
	_, _, goStats := runGoFastpPE(t, r1, r2, opts)
	goJSON := jsonFromStats(t, goStats)

	assertStringEqual(t, "PE sequencing", upJSON, goJSON, "summary.sequencing")
	for _, sect := range []string{
		"read1_before_filtering", "read2_before_filtering",
		"read1_after_filtering", "read2_after_filtering",
	} {
		assertReadStatsSubfields(t, "PE "+sect, upJSON, goJSON, sect)
	}

	// insert_size: deterministic peak + unknown bucket + full histogram.
	assertCounters(t, "PE insert_size scalars", upJSON, goJSON,
		"insert_size.peak", "insert_size.unknown")
	upHist, ok1 := lookupJSON(upJSON, "insert_size.histogram")
	goHist, ok2 := lookupJSON(goJSON, "insert_size.histogram")
	if !ok1 || !ok2 {
		t.Fatalf("PE insert_size.histogram missing (upstream=%v go=%v)", ok1, ok2)
	}
	ul := upHist.([]any)
	gl := goHist.([]any)
	if len(ul) != len(gl) {
		t.Fatalf("PE insert_size.histogram length: upstream=%d go=%d", len(ul), len(gl))
	}
	for i := range ul {
		uf, _ := toFloat(ul[i])
		gf, _ := toFloat(gl[i])
		if uf != gf {
			t.Fatalf("PE insert_size.histogram[%d]: upstream=%v go=%v", i, uf, gf)
		}
	}
}
