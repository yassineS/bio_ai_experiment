package main

import (
	"os"
	"path/filepath"
	"testing"
)

// TestCheckGtcheckDeferred locks in the upstream-flag-name surface
// that runGtcheck hard-rejects rather than silently accepting. Per
// the project parity rule (docs/PARITY_ROADMAP.md "Definition of 1:1")
// every documented upstream flag must be recognised — either
// implemented or gracefully rejected with a roadmap pointer. A
// regression that drops any of these from the rejection set without
// implementing the underlying behaviour is a parity bug.
func TestCheckGtcheckDeferred(t *testing.T) {
	if got := checkGtcheckDeferred(checkGtcheckDeferredInputs{outputType: "t"}); got != "" {
		t.Fatalf("default outputType=t: got deferred=%q, want \"\"", got)
	}
	// PL scoring, --n-matches and --distinctive-sites are now
	// implemented; only the -c/--cluster dendrogram (which upstream
	// itself rejects) and -O z compressed output remain deferred.
	cases := []struct {
		name string
		in   checkGtcheckDeferredInputs
		want string
	}{
		// -c/--cluster is now implemented, so it is no longer deferred.
		{"cluster", checkGtcheckDeferredInputs{cluster: "2,4"}, ""},
		{"output-type-z", checkGtcheckDeferredInputs{outputType: "z"}, "-O z (compressed output)"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkGtcheckDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestRohOzOutputsBGZF verifies that roh -O z now emits BGZF-compressed
// output (the gzip magic bytes) instead of being rejected.
func TestRohOzOutputsBGZF(t *testing.T) {
	dir := t.TempDir()
	vcf := writeRohTempFile(t, dir, "roh.vcf", "##fileformat=VCFv4.2\n"+
		`##INFO=<ID=AF,Number=A,Type=Float,Description="af">`+"\n"+
		`##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">`+"\n"+
		"##contig=<ID=1>\n"+
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"+
		"1\t100\t.\tA\tT\t50\tPASS\tAF=0.2\tGT\t0/0\n"+
		"1\t200\t.\tG\tC\t50\tPASS\tAF=0.3\tGT\t1/1\n")
	out := filepath.Join(dir, "out.roh.gz")
	if rc := runRoh([]string{"-G30", "--AF-tag", "AF", "-Osrz", "-o", out, vcf}); rc != 0 {
		t.Fatalf("roh -Oz: rc=%d, want 0", rc)
	}
	b, err := os.ReadFile(out)
	if err != nil {
		t.Fatalf("read output: %v", err)
	}
	if len(b) < 2 || b[0] != 0x1f || b[1] != 0x8b {
		t.Fatalf("roh -Oz output is not gzip/BGZF framed: % x", b[:min(4, len(b))])
	}
}

// writeRohTempFile writes content to dir/name and returns its path.
func writeRohTempFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", p, err)
	}
	return p
}

// TestSplitSamplesArg validates the upstream "qry:" / "gt:" prefix
// handling on `-s`. Un-prefixed lists are rejected to match upstream
// vcfgtcheck.c's "Which one? Query samples (qry:...) or genotype
// samples (gt:...)?" diagnostic.
func TestSplitSamplesArg(t *testing.T) {
	cases := []struct {
		in      string
		wantQry []string
		wantGT  []string
		wantErr bool
	}{
		{"", nil, nil, false},
		{"qry:a,b", []string{"a", "b"}, nil, false},
		{"gt:x,y,z", nil, []string{"x", "y", "z"}, false},
		{"a,b,c", nil, nil, true},
	}
	for _, c := range cases {
		gotQ, gotG, gotErr := splitSamplesArg(c.in)
		if (gotErr != nil) != c.wantErr {
			t.Errorf("splitSamplesArg(%q) err = %v, wantErr = %v", c.in, gotErr, c.wantErr)
		}
		if c.wantErr {
			continue
		}
		if !stringsEqual(gotQ, c.wantQry) {
			t.Errorf("splitSamplesArg(%q) qry = %v, want %v", c.in, gotQ, c.wantQry)
		}
		if !stringsEqual(gotG, c.wantGT) {
			t.Errorf("splitSamplesArg(%q) gt = %v, want %v", c.in, gotG, c.wantGT)
		}
	}
}

func stringsEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
