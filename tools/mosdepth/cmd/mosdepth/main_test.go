package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// makeBAM is a slim duplicate of the same helper inside the pkg/mosdepth
// tests; kept local so this test package doesn't depend on test-only
// exports.
func makeBAM(t *testing.T, refs []sam.Reference, recs []*sam.Record) []byte {
	t.Helper()
	hdr := &sam.Header{Refs: refs}
	for _, r := range refs {
		hdr.Lines = append(hdr.Lines, sam.HeaderLine{
			Tag: "SQ",
			Fields: []sam.HeaderField{
				{Tag: "SN", Value: r.Name},
				{Tag: "LN", Value: itoa(int(r.Length))},
			},
		})
	}
	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, rec := range recs {
		if err := bw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

func itoa(v int) string {
	if v == 0 {
		return "0"
	}
	neg := v < 0
	if neg {
		v = -v
	}
	d := []byte{}
	for v > 0 {
		d = append([]byte{byte('0' + v%10)}, d...)
		v /= 10
	}
	if neg {
		return "-" + string(d)
	}
	return string(d)
}

func TestRun_HelpAndVersion(t *testing.T) {
	if rc := run([]string{"-h"}); rc != 0 {
		t.Errorf("help rc: %d", rc)
	}
	if rc := run([]string{"--version"}); rc != 0 {
		t.Errorf("version rc: %d", rc)
	}
}

func TestRun_BadPositional(t *testing.T) {
	if rc := run([]string{}); rc != 2 {
		t.Errorf("expected rc=2 for missing args, got %d", rc)
	}
}

func TestRun_BadFlag(t *testing.T) {
	if rc := run([]string{"--bogus-flag"}); rc != 2 {
		t.Errorf("expected rc=2 for bad flag, got %d", rc)
	}
}

func TestRun_Successful(t *testing.T) {
	dir := t.TempDir()
	cigar, err := sam.ParseCigar("5M")
	if err != nil {
		t.Fatal(err)
	}
	bam := makeBAM(t, []sam.Reference{{Name: "chr1", Length: 20}}, []*sam.Record{
		{QName: "r1", RName: "chr1", Pos: 1, Cigar: cigar, MapQ: 60, Seq: "AAAAA"},
	})
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bam, 0644); err != nil {
		t.Fatal(err)
	}
	prefix := filepath.Join(dir, "out")
	rc := run([]string{prefix, bamPath})
	if rc != 0 {
		t.Errorf("run rc: %d", rc)
	}
	// per-base file produced.
	if _, err := os.Stat(prefix + ".per-base.bed.gz"); err != nil {
		t.Errorf("per-base missing: %v", err)
	}
}

// TestRun_D4Output checks that `-d` writes a per-base D4 file (and no BED).
func TestRun_D4Output(t *testing.T) {
	dir := t.TempDir()
	cigar, _ := sam.ParseCigar("5M")
	bam := makeBAM(t, []sam.Reference{{Name: "chr1", Length: 20}}, []*sam.Record{
		{QName: "r1", RName: "chr1", Pos: 1, Cigar: cigar, MapQ: 60, Seq: "AAAAA"},
	})
	bamPath := filepath.Join(dir, "in.bam")
	os.WriteFile(bamPath, bam, 0644)
	prefix := filepath.Join(dir, "out")
	if rc := run([]string{"-d", prefix, bamPath}); rc != 0 {
		t.Fatalf("D4 rc: got %d, want 0", rc)
	}
	if _, err := os.Stat(prefix + ".per-base.d4"); err != nil {
		t.Errorf("per-base.d4 missing: %v", err)
	}
	if _, err := os.Stat(prefix + ".per-base.bed.gz"); err == nil {
		t.Errorf("per-base.bed.gz should not be written with -d")
	}
}

func TestRun_BadThresholds(t *testing.T) {
	if rc := run([]string{"-T", "abc", "/tmp/x", "/tmp/x.bam"}); rc != 1 {
		t.Errorf("bad thresholds rc: %d", rc)
	}
}

func TestSplitCSV(t *testing.T) {
	got := splitCSV(",a,,b,c,")
	want := []string{"a", "b", "c"}
	if strings.Join(got, ",") != strings.Join(want, ",") {
		t.Errorf("splitCSV: got %v, want %v", got, want)
	}
}

func TestUint8Clamp(t *testing.T) {
	if uint8Clamp(-1) != 0 || uint8Clamp(255) != 255 || uint8Clamp(300) != 255 || uint8Clamp(42) != 42 {
		t.Errorf("uint8Clamp behaves unexpectedly")
	}
}

// TestParseFlags_DocoptBundling proves the cliflag.Parse routing gives
// mosdepth the same clustered short-flag and value-concatenation parsing
// that upstream's docopt parser provides: "-nx" == "-n -x" and "-Q20" ==
// "-Q 20". We compare the parsed runOptions of the bundled and canonical
// argv forms field-by-field (positionals included) so any divergence
// fails the test.
func TestParseFlags_DocoptBundling(t *testing.T) {
	cases := []struct {
		name      string
		bundled   []string
		canonical []string
	}{
		{"two-bools", []string{"-nx", "p", "b.bam"}, []string{"-n", "-x", "p", "b.bam"}},
		{"bool-then-value-concat", []string{"-xQ20", "p", "b.bam"}, []string{"-x", "-Q", "20", "p", "b.bam"}},
		{"value-concat", []string{"-Q30", "p", "b.bam"}, []string{"-Q", "30", "p", "b.bam"}},
		{"value-next-arg", []string{"-Q", "30", "p", "b.bam"}, []string{"-Q", "30", "p", "b.bam"}},
		{"bool-cluster-then-value-next", []string{"-nxQ", "30", "p", "b.bam"}, []string{"-n", "-x", "-Q", "30", "p", "b.bam"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			gotOpts, gotPos, err := parseFlags(tc.bundled)
			if err != nil {
				t.Fatalf("bundled parse %v: %v", tc.bundled, err)
			}
			wantOpts, wantPos, err := parseFlags(tc.canonical)
			if err != nil {
				t.Fatalf("canonical parse %v: %v", tc.canonical, err)
			}
			if *gotOpts != *wantOpts {
				t.Errorf("bundled %v parsed to %+v, canonical %v parsed to %+v",
					tc.bundled, *gotOpts, tc.canonical, *wantOpts)
			}
			if strings.Join(gotPos, ",") != strings.Join(wantPos, ",") {
				t.Errorf("positionals: bundled %v, canonical %v", gotPos, wantPos)
			}
		})
	}
}

// TestParseFlags_ReadGroupsUpstreamShort confirms -R is the upstream short
// for --read-groups and that the lowercase -r port alias resolves to the
// same value, with -R winning when both are supplied.
func TestParseFlags_ReadGroupsUpstreamShort(t *testing.T) {
	for _, tc := range []struct {
		args []string
		want string
	}{
		{[]string{"-R", "rg1,rg2", "p", "b.bam"}, "rg1,rg2"},
		{[]string{"--read-groups", "rg3", "p", "b.bam"}, "rg3"},
		{[]string{"-r", "rg4", "p", "b.bam"}, "rg4"},
		{[]string{"-r", "lower", "-R", "upper", "p", "b.bam"}, "upper"},
	} {
		opts, _, err := parseFlags(tc.args)
		if err != nil {
			t.Fatalf("parse %v: %v", tc.args, err)
		}
		if opts.readGroups != tc.want {
			t.Errorf("parse %v: readGroups=%q, want %q", tc.args, opts.readGroups, tc.want)
		}
	}
}

// TestRun_ConflictingFlagsRejected confirms mutually exclusive flag
// combinations are rejected with exit code 2 rather than silently producing
// divergent output. (-m/--use-median is now implemented, so it is no longer
// rejected — see TestParseFlags_UseMedianAccepted.)
func TestRun_ConflictingFlagsRejected(t *testing.T) {
	for _, args := range [][]string{
		{"-a", "-x", "p", "b.bam"}, // fragment-mode and fast-mode conflict
	} {
		if rc := run(args); rc != 2 {
			t.Errorf("run %v: rc=%d, want 2", args, rc)
		}
	}
}

// TestParseFlags_UseMedianAccepted confirms -m/--use-median now parses and
// feeds opts.useMedian instead of being rejected.
func TestParseFlags_UseMedianAccepted(t *testing.T) {
	for _, args := range [][]string{
		{"-m", "p", "b.bam"},
		{"--use-median", "p", "b.bam"},
	} {
		opts, _, err := parseFlags(args)
		if err != nil {
			t.Fatalf("parse %v: %v", args, err)
		}
		if !opts.useMedian {
			t.Errorf("parse %v: useMedian not set", args)
		}
	}
}

// TestParseFlags_FragmentAndQuantizeAccepted confirms -a/--fragment-mode and
// -q/--quantize now parse and feed the Options without being rejected.
func TestParseFlags_FragmentAndQuantizeAccepted(t *testing.T) {
	opts, _, err := parseFlags([]string{"-a", "--quantize", "0:1:4", "p", "b.bam"})
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if !opts.fragmentLen {
		t.Errorf("fragmentLen not set")
	}
	if opts.quantize != "0:1:4" {
		t.Errorf("quantize=%q, want 0:1:4", opts.quantize)
	}
}

// TestParseFlags_FastaAccepted confirms -f/--fasta parses (CRAM reference)
// without error even though CRAM input is not yet supported.
func TestParseFlags_FastaAccepted(t *testing.T) {
	opts, _, err := parseFlags([]string{"-f", "ref.fa", "p", "b.bam"})
	if err != nil {
		t.Fatalf("parse -f: %v", err)
	}
	if opts.fasta != "ref.fa" {
		t.Errorf("fasta=%q, want %q", opts.fasta, "ref.fa")
	}
}
