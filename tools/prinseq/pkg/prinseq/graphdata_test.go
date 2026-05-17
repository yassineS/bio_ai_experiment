package prinseq

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"sort"
	"strconv"
	"strings"
	"testing"
)

// readGDBody strips the two leading `#`-comment header lines (which
// carry the upstream timestamp + argv) and returns the JSON body
// only. The upstream emitter always writes exactly two comment
// lines before the JSON, but we tolerate any number of leading
// `#` lines for robustness.
func readGDBody(t *testing.T, path string) []byte {
	t.Helper()
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	// Skip comment lines.
	idx := 0
	for idx < len(raw) && raw[idx] == '#' {
		// advance to newline
		for idx < len(raw) && raw[idx] != '\n' {
			idx++
		}
		if idx < len(raw) {
			idx++
		}
	}
	return raw[idx:]
}

// normaliseGD parses a `.gd` body (JSON) and recursively coerces
// numeric strings to numbers. Upstream's mean/std fields are
// JSON-quoted strings like "19.58"; we treat them as numbers for
// semantic diff. complvals.{minseq,maxseq} stay strings.
func normaliseGD(t *testing.T, body []byte) any {
	t.Helper()
	var v any
	if err := json.Unmarshal(body, &v); err != nil {
		t.Fatalf("unmarshal .gd: %v\nbody: %s", err, truncForLog(body, 400))
	}
	return coerceNumericStrings(v, "")
}

func truncForLog(b []byte, n int) string {
	if len(b) <= n {
		return string(b)
	}
	return string(b[:n]) + "..."
}

// coerceNumericStrings walks the JSON tree, parsing any string
// that matches a float/int pattern into a float64. The `path`
// argument is used to suppress coercion for the dust/entropy
// `minseq`/`maxseq` slots, which are nucleotide strings even
// though some look like they could be numeric.
func coerceNumericStrings(v any, path string) any {
	switch x := v.(type) {
	case map[string]any:
		out := make(map[string]any, len(x))
		for k, vv := range x {
			child := path + "/" + k
			out[k] = coerceNumericStrings(vv, child)
		}
		return out
	case []any:
		out := make([]any, len(x))
		for i, e := range x {
			out[i] = coerceNumericStrings(e, path)
		}
		return out
	case string:
		// Leave the complvals minseq/maxseq + filename1/filename2 fields as strings.
		if strings.HasSuffix(path, "/minseq") || strings.HasSuffix(path, "/maxseq") ||
			strings.HasSuffix(path, "/filename1") || strings.HasSuffix(path, "/filename2") ||
			strings.HasSuffix(path, "/format1") || strings.HasSuffix(path, "/format2") ||
			strings.HasSuffix(path, "/tagmidseq") {
			return x
		}
		if f, err := strconv.ParseFloat(x, 64); err == nil {
			return f
		}
		return x
	case float64:
		return x
	case bool:
		return x
	case nil:
		return nil
	default:
		return x
	}
}

// gdDiff returns a list of structural differences between two
// parsed/coerced JSON trees. Floats are compared with a small
// absolute tolerance (1e-6) to absorb the formatting round-trip
// of mean/std/dinucodds.
func gdDiff(a, b any, path string, tol float64) []string {
	switch av := a.(type) {
	case map[string]any:
		bv, ok := b.(map[string]any)
		if !ok {
			return []string{fmt.Sprintf("%s: type mismatch (got %T vs %T)", path, a, b)}
		}
		var diffs []string
		for k, v := range av {
			cp := path + "/" + k
			vv, ok := bv[k]
			if !ok {
				diffs = append(diffs, fmt.Sprintf("%s: missing in right tree", cp))
				continue
			}
			diffs = append(diffs, gdDiff(v, vv, cp, tol)...)
		}
		for k := range bv {
			if _, ok := av[k]; !ok {
				diffs = append(diffs, fmt.Sprintf("%s/%s: missing in left tree", path, k))
			}
		}
		return diffs
	case []any:
		bv, ok := b.([]any)
		if !ok {
			return []string{fmt.Sprintf("%s: type mismatch (got %T vs %T)", path, a, b)}
		}
		if len(av) != len(bv) {
			return []string{fmt.Sprintf("%s: array length %d vs %d", path, len(av), len(bv))}
		}
		var diffs []string
		for i := range av {
			diffs = append(diffs, gdDiff(av[i], bv[i], fmt.Sprintf("%s[%d]", path, i), tol)...)
		}
		return diffs
	case float64:
		bv, ok := b.(float64)
		if !ok {
			return []string{fmt.Sprintf("%s: type mismatch (float vs %T)", path, b)}
		}
		if math.Abs(av-bv) > tol {
			return []string{fmt.Sprintf("%s: %g vs %g (Δ=%g)", path, av, bv, av-bv)}
		}
		return nil
	case string:
		bv, ok := b.(string)
		if !ok {
			return []string{fmt.Sprintf("%s: type mismatch (string vs %T)", path, b)}
		}
		if av != bv {
			return []string{fmt.Sprintf("%s: %q vs %q", path, av, bv)}
		}
		return nil
	case nil:
		if b == nil {
			return nil
		}
		return []string{fmt.Sprintf("%s: nil vs %v", path, b)}
	default:
		if !reflect.DeepEqual(a, b) {
			return []string{fmt.Sprintf("%s: %v vs %v", path, a, b)}
		}
		return nil
	}
}

// TestGraphDataParityExample1 is the headline parity test. It
// drives our collector over the upstream-vendored `example1.fastq`
// fixture, parses both our emit and the upstream `.gd`, and asserts
// the parsed trees are equal up to a numerical tolerance. The
// upstream `.gd` carries Perl-hash random key order; comparing
// parsed maps removes that source of non-determinism.
func TestGraphDataParityExample1(t *testing.T) {
	fixturesDir := filepath.Join("..", "..", "testdata", "parity")
	fastqPath := filepath.Join(fixturesDir, "graphdata_example1.fastq")
	upstreamGDPath := filepath.Join(fixturesDir, "graphdata_example1.gd")

	f, err := os.Open(fastqPath)
	if err != nil {
		t.Fatalf("open fastq: %v", err)
	}
	defer f.Close()

	opts := DefaultGraphDataOptions()
	// Upstream embeds filename1 as the hex-encoded basename
	// ("6578616d706c65312e6661737471" = "example1.fastq"). Replicate
	// that so the filename1 comparison passes.
	opts.Filename1 = "6578616d706c65312e6661737471"
	g, err := CollectGraphData(f, true, opts)
	if err != nil {
		t.Fatalf("collect: %v", err)
	}

	var buf bytes.Buffer
	if err := g.EmitGD(&buf, GDHeader{}); err != nil {
		t.Fatalf("emit: %v", err)
	}

	upstreamBody := readGDBody(t, upstreamGDPath)
	ourBody := buf.Bytes()

	upstream := normaliseGD(t, upstreamBody)
	ours := normaliseGD(t, ourBody)

	diffs := gdDiff(upstream, ours, "", 1e-3)
	if len(diffs) > 0 {
		sort.Strings(diffs)
		// Truncate to a manageable number of diffs for the log.
		max := len(diffs)
		if max > 60 {
			max = 60
		}
		t.Fatalf("graph-data divergence (%d entries):\n%s", len(diffs),
			strings.Join(diffs[:max], "\n"))
	}
}

// TestGraphDataDeterminism asserts that running the collector
// twice on the same input produces byte-identical output. This
// is the contract that distinguishes us from upstream (Perl-hash
// randomness) and the strongest reason the parity test above is
// safe to commit.
func TestGraphDataDeterminism(t *testing.T) {
	fastqPath := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.fastq")
	run := func() []byte {
		f, err := os.Open(fastqPath)
		if err != nil {
			t.Fatalf("open: %v", err)
		}
		defer f.Close()
		opts := DefaultGraphDataOptions()
		opts.Filename1 = "6578616d706c65312e6661737471"
		g, err := CollectGraphData(f, true, opts)
		if err != nil {
			t.Fatalf("collect: %v", err)
		}
		var buf bytes.Buffer
		if err := g.EmitGD(&buf, GDHeader{}); err != nil {
			t.Fatalf("emit: %v", err)
		}
		return buf.Bytes()
	}
	a := run()
	b := run()
	if !bytes.Equal(a, b) {
		t.Fatalf("non-deterministic output: lengths %d vs %d", len(a), len(b))
	}
	// Also assert valid JSON.
	var v any
	if err := json.Unmarshal(a, &v); err != nil {
		t.Fatalf("emitted body is not valid JSON: %v\n%s", err, truncForLog(a, 300))
	}
}

// TestGenerateStatsType exercises the quantile/std/mode roll-up
// against hand-computed expectations.
func TestGenerateStatsType(t *testing.T) {
	cases := []struct {
		name   string
		counts map[int]int
		want   statsType
	}{
		{
			name:   "single value",
			counts: map[int]int{42: 1},
			want: statsType{
				Min: 42, Max: 42, Range: 1, Modeval: 1, Mode: 42,
				Median: 42, P25: 42, P75: 42,
				Mean: "42.00", Std: "0.00",
			},
		},
		{
			name:   "two same values",
			counts: map[int]int{10: 2},
			want: statsType{
				Min: 10, Max: 10, Range: 1, Modeval: 2, Mode: 10,
				Median: 10, P25: 10, P75: 10,
				Mean: "10.00", Std: "0.00",
			},
		},
		{
			name:   "two different values",
			counts: map[int]int{1: 1, 3: 1},
			want: statsType{
				Min: 1, Max: 3, Range: 3, Modeval: 1, Mode: 1,
				Median: 2, P25: 1, P75: 3,
				Mean: "2.00", Std: "1.00",
			},
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := generateStatsType(c.counts)
			if got != c.want {
				t.Fatalf("generateStatsType(%v) = %+v, want %+v", c.counts, got, c.want)
			}
		})
	}
}

// TestGetBinVal walks the four upstream branches of getBinVal.
func TestGetBinVal(t *testing.T) {
	cases := []struct {
		in, want int
	}{
		{0, 1},
		{50, 1},
		{100, 1},
		{101, 2},
		{200, 2},
		{250, 3},
		{9999, 100},
		{10000, 1000},
		{99999, 1000},
		{100000, 10000}, // val=100000, step=1000000, val%step=100000 -> xmax=1000000, /100=10000
	}
	for _, c := range cases {
		if got := getBinVal(c.in); got != c.want {
			t.Errorf("getBinVal(%d)=%d, want %d", c.in, got, c.want)
		}
	}
}

// TestParseGraphStatsCSV exercises both the happy path and the
// unknown-option error.
func TestParseGraphStatsCSV(t *testing.T) {
	opts := DefaultGraphDataOptions()
	if err := ParseGraphStatsCSV("gc,qd,ns", &opts); err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if opts.LD || !opts.GC || !opts.QD || !opts.NS || opts.PT {
		t.Fatalf("unexpected option mask: %+v", opts)
	}
	if err := ParseGraphStatsCSV("nope", &opts); err == nil {
		t.Fatalf("expected error for unknown option")
	}
}

// TestComplexityScores sanity-checks the dust/entropy scoring on
// canonical low- and high-complexity sequences.
func TestComplexityScores(t *testing.T) {
	// All-A — maximum dust, zero entropy.
	dust, entropy := complexityScores(strings.Repeat("A", 100))
	// Map through the upstream scaling for sanity.
	dustVal := int(dust * 100 / 31)
	entropyVal := int(entropy * 100)
	if dustVal < 50 {
		t.Errorf("low-complexity dust scaled to %d, expected >= 50", dustVal)
	}
	if entropyVal > 5 {
		t.Errorf("low-complexity entropy scaled to %d, expected near 0", entropyVal)
	}
	// Random-ish — much lower dust score.
	dust, entropy = complexityScores("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
	dustVal = int(dust * 100 / 31)
	entropyVal = int(entropy * 100)
	if dustVal > 50 {
		t.Errorf("high-complexity dust scaled to %d, expected lower", dustVal)
	}
	_ = entropyVal
}

// TestSplitNRuns confirms our split mimics Perl's split(/N+/).
func TestSplitNRuns(t *testing.T) {
	cases := []struct {
		in   string
		want []string
	}{
		{"ACGT", []string{"ACGT"}},
		{"AANCGNNTA", []string{"AA", "CG", "TA"}},
		{"NNN", nil},
		{"NACGTN", []string{"ACGT"}},
		{"", nil},
	}
	for _, c := range cases {
		got := splitNRuns(c.in)
		if !reflect.DeepEqual(got, c.want) {
			t.Errorf("splitNRuns(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

// TestEmitGDHeader confirms the header block is omitted when no
// fields are set and emitted (with placeholders) when at least
// one is set.
func TestEmitGDHeader(t *testing.T) {
	g := NewGraphData(DefaultGraphDataOptions())
	g.AddSeq("ACGT")
	g.AddQual([]int{30, 30, 30, 30})
	var buf bytes.Buffer
	if err := g.EmitGD(&buf, GDHeader{}); err != nil {
		t.Fatalf("emit: %v", err)
	}
	if bytes.HasPrefix(buf.Bytes(), []byte("#")) {
		t.Fatalf("expected no #-header when GDHeader is zero-valued")
	}

	buf.Reset()
	if err := g.EmitGD(&buf, GDHeader{Version: "0.20.4", Timestamp: "01/01/2026 00:00:00", Command: "prinseq -graph_data"}); err != nil {
		t.Fatalf("emit: %v", err)
	}
	if !bytes.HasPrefix(buf.Bytes(), []byte("#Graph data\n#[prinseq-lite-0.20.4]")) {
		t.Fatalf("expected upstream-shaped #-header prefix, got: %q", truncForLog(buf.Bytes(), 100))
	}
}

// TestResolveGraphDataPath checks the upstream default convention.
func TestResolveGraphDataPath(t *testing.T) {
	if p := ResolveGraphDataPath("custom.gd", "x.fastq"); p != "custom.gd" {
		t.Errorf("explicit path overridden: %q", p)
	}
	if p := ResolveGraphDataPath("", "x.fastq"); p != "x.fastq__.gd" {
		t.Errorf("default path = %q, want x.fastq__.gd", p)
	}
	if p := ResolveGraphDataPath("", ""); p != "nonamegiven__.gd" {
		t.Errorf("stdin default = %q, want nonamegiven__.gd", p)
	}
}
