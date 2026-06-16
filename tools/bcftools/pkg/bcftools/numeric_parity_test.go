package bcftools

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Field-aware, tolerance-aware ("proximity") parity comparison for text output
// whose numeric fields are computed through libm transcendentals.
//
// Rationale: byte parity is the right bar for integer/combinatorial output
// (e.g. the trio-dnm3 NAIVE model, which is a pure Mendelian table lookup), but
// it is the WRONG bar for the trio-dnm3 float models (DMM/ALM/DNG). Those
// scores are log/exp/pow/lgamma reductions over hundreds of terms. The C build
// links the platform libm; our Go port uses math.Log/Exp/Pow and the in-tree
// AS245 kf_lgamma. Even when both are correctly rounded, a transcendental such
// as log(x) is only guaranteed to the last ULP, and a long sum of such terms
// accumulates a handful of ULPs of disagreement. After htslib narrows the score
// to a 32-bit float and prints it with %g (6 significant figures), that shows up
// as a difference in the last printed digit (e.g. -46.0521 vs -46.0522). That is
// not a bug; it is the floating-point reproducibility boundary. Insisting on
// byte parity there would be testing the libm implementation, not our port.
//
// compareProximity therefore splits each line into fields, compares STRING
// fields exactly, and compares NUMERIC fields (including %e/%f/%g forms and the
// nan/-nan/inf spellings) as equal when they agree after rounding to a
// configurable number of significant figures OR fall within a small combined
// relative+absolute epsilon. nan and -nan compare equal to each other (and to
// "nan"); +inf/-inf must match in sign. The reported diff names the line, the
// field index, both raw values, and the absolute/relative delta so a genuine
// model divergence is still caught and localised.

// proximityTolerance configures the numeric closeness test used by
// compareProximity. SigFigs is the number of leading significant decimal digits
// that must agree; RelEps and AbsEps are the fallback relative and absolute
// epsilons (a field passes if EITHER the rounded-sig-fig forms match OR the
// combined epsilon test passes). FieldSep, when non-empty, lists the byte
// runes that separate fields; the default splits on any run of spaces and tabs.
type proximityTolerance struct {
	SigFigs int
	RelEps  float64
	AbsEps  float64
}

// defaultProximityTolerance is the tolerance used by compareProximityDefault:
// agreement to ~6 significant figures, or within a relative 1e-5 / absolute
// 1e-6 epsilon. Six sig-figs matches htslib's %g default precision, so the
// common "last printed digit differs" case is absorbed, while a wrong model
// (off by whole points of phred score) still fails.
var defaultProximityTolerance = proximityTolerance{SigFigs: 6, RelEps: 1e-5, AbsEps: 1e-6}

// proximityDiff describes a single field-level mismatch found by
// compareProximity.
type proximityDiff struct {
	Line     int     // 1-based line number
	Field    int     // 0-based field index within the line
	Want     string  // upstream (reference) raw field text
	Got      string  // our (candidate) raw field text
	AbsDelta float64 // |want-got| when both numeric, else NaN
	RelDelta float64 // |want-got|/max(|want|,|got|) when both numeric, else NaN
	Reason   string  // human-readable cause
}

// String renders a proximityDiff as a single diagnostic line.
func (d proximityDiff) String() string {
	if math.IsNaN(d.AbsDelta) {
		return fmt.Sprintf("line %d field %d: %s (want %q, got %q)", d.Line, d.Field, d.Reason, d.Want, d.Got)
	}
	return fmt.Sprintf("line %d field %d: %s (want %q, got %q, |Δ|=%g rel=%g)",
		d.Line, d.Field, d.Reason, d.Want, d.Got, d.AbsDelta, d.RelDelta)
}

// compareProximityDefault is compareProximity with defaultProximityTolerance.
func compareProximityDefault(want, got string) []proximityDiff {
	return compareProximity(want, got, defaultProximityTolerance)
}

// compareProximity compares two multi-line texts field by field with the given
// tolerance. String fields must match exactly; numeric fields must be close per
// numericFieldsClose. It returns one proximityDiff per mismatch (capped, see
// below) or nil when the texts are proximity-equal. A differing line/field count
// is itself reported as a diff.
func compareProximity(want, got string, tol proximityTolerance) []proximityDiff {
	wantLines := splitProximityLines(want)
	gotLines := splitProximityLines(got)
	var diffs []proximityDiff

	maxLines := len(wantLines)
	if len(gotLines) > maxLines {
		maxLines = len(gotLines)
	}
	for li := 0; li < maxLines; li++ {
		if li >= len(wantLines) {
			diffs = append(diffs, proximityDiff{Line: li + 1, Field: -1, Want: "", Got: gotLines[li],
				AbsDelta: math.NaN(), RelDelta: math.NaN(), Reason: "extra line in candidate output"})
			continue
		}
		if li >= len(gotLines) {
			diffs = append(diffs, proximityDiff{Line: li + 1, Field: -1, Want: wantLines[li], Got: "",
				AbsDelta: math.NaN(), RelDelta: math.NaN(), Reason: "missing line in candidate output"})
			continue
		}
		diffs = append(diffs, compareLine(li+1, wantLines[li], gotLines[li], tol)...)
		if len(diffs) > 50 {
			break // cap diagnostics; a flood means the outputs are structurally different
		}
	}
	return diffs
}

// splitProximityLines splits text into lines, dropping a single trailing empty line so a
// terminal newline does not register as an extra line.
func splitProximityLines(s string) []string {
	lines := strings.Split(s, "\n")
	if n := len(lines); n > 0 && lines[n-1] == "" {
		lines = lines[:n-1]
	}
	return lines
}

// compareLine splits a single line into whitespace-separated fields and compares
// them pairwise.
func compareLine(lineNo int, want, got string, tol proximityTolerance) []proximityDiff {
	wf := strings.Fields(want)
	gf := strings.Fields(got)
	var diffs []proximityDiff
	maxF := len(wf)
	if len(gf) > maxF {
		maxF = len(gf)
	}
	for fi := 0; fi < maxF; fi++ {
		if fi >= len(wf) {
			diffs = append(diffs, proximityDiff{Line: lineNo, Field: fi, Want: "", Got: gf[fi],
				AbsDelta: math.NaN(), RelDelta: math.NaN(), Reason: "extra field in candidate line"})
			continue
		}
		if fi >= len(gf) {
			diffs = append(diffs, proximityDiff{Line: lineNo, Field: fi, Want: wf[fi], Got: "",
				AbsDelta: math.NaN(), RelDelta: math.NaN(), Reason: "missing field in candidate line"})
			continue
		}
		if d, ok := compareField(lineNo, fi, wf[fi], gf[fi], tol); !ok {
			diffs = append(diffs, d)
		}
	}
	return diffs
}

// compareField compares one field. A field is a single token that may itself be
// a colon-joined composite (e.g. a VCF sample column "0|1:10,10:-46.0521:0:50")
// or a comma-joined list. To stay field-aware without re-implementing the VCF
// grammar, it recurses on ':' and ',' sub-tokens so that a numeric sub-value is
// compared numerically even when embedded in a structured field. Atomic tokens
// fall through to the exact/numeric leaf comparison.
func compareField(lineNo, fi int, want, got string, tol proximityTolerance) (proximityDiff, bool) {
	for _, sep := range []byte{':', ','} {
		if strings.IndexByte(want, sep) >= 0 || strings.IndexByte(got, sep) >= 0 {
			ws := strings.Split(want, string(sep))
			gs := strings.Split(got, string(sep))
			if len(ws) != len(gs) {
				return proximityDiff{Line: lineNo, Field: fi, Want: want, Got: got,
					AbsDelta: math.NaN(), RelDelta: math.NaN(),
					Reason: fmt.Sprintf("sub-token count differs on %q", string(sep))}, false
			}
			for i := range ws {
				if d, ok := compareField(lineNo, fi, ws[i], gs[i], tol); !ok {
					d.Want, d.Got = want, got // report the whole field for context
					return d, false
				}
			}
			return proximityDiff{}, true
		}
	}
	return compareLeaf(lineNo, fi, want, got, tol)
}

// compareLeaf compares two atomic tokens: numbers numerically (with tolerance),
// everything else byte-for-byte.
func compareLeaf(lineNo, fi int, want, got string, tol proximityTolerance) (proximityDiff, bool) {
	wn, wIsNum := parseNumericField(want)
	gn, gIsNum := parseNumericField(got)
	if wIsNum != gIsNum {
		return proximityDiff{Line: lineNo, Field: fi, Want: want, Got: got,
			AbsDelta: math.NaN(), RelDelta: math.NaN(),
			Reason: "one field numeric, the other not"}, false
	}
	if !wIsNum {
		if want == got {
			return proximityDiff{}, true
		}
		return proximityDiff{Line: lineNo, Field: fi, Want: want, Got: got,
			AbsDelta: math.NaN(), RelDelta: math.NaN(), Reason: "string fields differ"}, false
	}
	if numericClose(wn, gn, tol) {
		return proximityDiff{}, true
	}
	abs := math.Abs(wn - gn)
	rel := abs / math.Max(math.Abs(wn), math.Abs(gn))
	return proximityDiff{Line: lineNo, Field: fi, Want: want, Got: got,
		AbsDelta: abs, RelDelta: rel, Reason: "numeric fields outside tolerance"}, false
}

// parseNumericField reports whether s is a numeric token (decimal, scientific
// %e, or the nan/-nan/inf/-inf spellings, optionally signed) and returns its
// float64 value. NaN spellings yield a NaN; infinities yield ±Inf.
func parseNumericField(s string) (float64, bool) {
	if s == "" || s == "." {
		return 0, false // VCF "missing" is a structural marker, compared as a string
	}
	low := strings.ToLower(s)
	switch low {
	case "nan", "+nan", "-nan":
		return math.NaN(), true
	case "inf", "+inf", "infinity", "+infinity":
		return math.Inf(1), true
	case "-inf", "-infinity":
		return math.Inf(-1), true
	}
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, false
	}
	return f, true
}

// numericClose reports whether two numeric values agree to within the tolerance:
// NaNs match any NaN; infinities must match in sign; otherwise the values pass
// if their significant-figure roundings match OR the combined relative/absolute
// epsilon test passes.
func numericClose(a, b float64, tol proximityTolerance) bool {
	aNaN, bNaN := math.IsNaN(a), math.IsNaN(b)
	if aNaN || bNaN {
		return aNaN && bNaN // -nan == nan; a NaN never equals a finite number
	}
	aInf, bInf := math.IsInf(a, 0), math.IsInf(b, 0)
	if aInf || bInf {
		return aInf && bInf && math.Signbit(a) == math.Signbit(b)
	}
	if a == b {
		return true
	}
	abs := math.Abs(a - b)
	if abs <= tol.AbsEps {
		return true
	}
	if abs <= tol.RelEps*math.Max(math.Abs(a), math.Abs(b)) {
		return true
	}
	return roundSig(a, tol.SigFigs) == roundSig(b, tol.SigFigs)
}

// roundSig rounds x to n significant decimal figures, the same notion of
// "agree to N digits" htslib's %g printing uses. Zero and non-finite inputs are
// returned unchanged.
func roundSig(x float64, n int) float64 {
	if x == 0 || n <= 0 || math.IsInf(x, 0) || math.IsNaN(x) {
		return x
	}
	d := math.Ceil(math.Log10(math.Abs(x)))
	power := float64(n) - d
	mag := math.Pow(10, power)
	return math.Round(x*mag) / mag
}
