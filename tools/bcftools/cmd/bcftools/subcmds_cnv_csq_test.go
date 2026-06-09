package main

import "testing"

// TestCheckCSQDeferred locks in the upstream-flag surface that runCSQ
// hard-rejects rather than silently accepting. Per the project parity
// rule (docs/PARITY_ROADMAP.md "Definition of 1:1") every documented
// upstream flag must be recognised — either implemented or gracefully
// rejected with a roadmap pointer. A regression that drops any of
// these from the rejection set without implementing the underlying
// behaviour is a parity bug.
func TestCheckCSQDeferred(t *testing.T) {
	if got := checkCSQDeferred(checkCSQDeferredInputs{}); got != "" {
		t.Fatalf("empty input: got deferred=%q, want \"\"", got)
	}
	if got := checkCSQDeferred(checkCSQDeferredInputs{outputType: "v"}); got != "" {
		t.Fatalf("-O v: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkCSQDeferredInputs
		want string
	}{
		// Unknown -O type must still be rejected with the format hint.
		{"-O t", checkCSQDeferredInputs{outputType: "t"}, "-O t (expect v|z|b|u)"},
		// Now-implemented flags must NOT be rejected:
		{"-O z", checkCSQDeferredInputs{outputType: "z"}, ""},
		{"-O b", checkCSQDeferredInputs{outputType: "b"}, ""},
		{"-O u", checkCSQDeferredInputs{outputType: "u"}, ""},
		{"unify-chr", checkCSQDeferredInputs{unifyChrNames: "chr,Chr,-"}, ""},
		{"dump-gff", checkCSQDeferredInputs{dumpGFF: "out.gff.gz"}, ""},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkCSQDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestParseCNVPair exercises the FLOAT[,FLOAT] parser used for
// -a/--aberrant, -d/--BAF-dev, -k/--LRR-dev.
func TestParseCNVPair(t *testing.T) {
	a, b, err := parseCNVPair("0.05", "-d/--BAF-dev")
	if err != nil || a != 0.05 || b != 0.05 {
		t.Errorf("single value: got (%v,%v,%v)", a, b, err)
	}
	a, b, err = parseCNVPair("0.05,0.1", "-d/--BAF-dev")
	if err != nil || a != 0.05 || b != 0.1 {
		t.Errorf("pair: got (%v,%v,%v)", a, b, err)
	}
	if _, _, err := parseCNVPair("bad", "-d/--BAF-dev"); err == nil {
		t.Errorf("expected error on bad input")
	}
	if _, _, err := parseCNVPair("0.05,bad", "-d/--BAF-dev"); err == nil {
		t.Errorf("expected error on bad second value")
	}
	a, b, err = parseCNVPair("", "-d/--BAF-dev")
	if err != nil || a != 0 || b != 0 {
		t.Errorf("empty: got (%v,%v,%v)", a, b, err)
	}
}

// TestCNVRunInputs sanity-checks the CLI binding (missing input,
// missing -o, help).
func TestCNVRunInputs(t *testing.T) {
	if rc := runCNV([]string{"--help"}); rc != 0 {
		t.Errorf("--help rc=%d want 0", rc)
	}
	if rc := runCNV([]string{}); rc != 2 {
		t.Errorf("no args rc=%d want 2", rc)
	}
	if rc := runCNV([]string{"some.vcf"}); rc != 2 {
		t.Errorf("no -o rc=%d want 2", rc)
	}
}

// TestCSQRunInputs sanity-checks the CLI binding.
func TestCSQRunInputs(t *testing.T) {
	if rc := runCSQ([]string{"--help"}); rc != 0 {
		t.Errorf("--help rc=%d want 0", rc)
	}
	if rc := runCSQ([]string{}); rc != 2 {
		t.Errorf("no args rc=%d want 2", rc)
	}
	if rc := runCSQ([]string{"some.vcf"}); rc != 2 {
		t.Errorf("missing -f and -g rc=%d want 2", rc)
	}
	if rc := runCSQ([]string{"-f", "ref.fa", "some.vcf"}); rc != 2 {
		t.Errorf("missing -g rc=%d want 2", rc)
	}
	// An unknown genetic-code table is rejected up front (rc=2) before
	// any file I/O.
	if rc := runCSQ([]string{"-C", "99", "-f", "ref.fa", "-g", "anno.gff", "some.vcf"}); rc != 2 {
		t.Errorf("-C 99 must be rejected, got rc=%d", rc)
	}
	// A non-numeric genetic code is rejected.
	if rc := runCSQ([]string{"-C", "xyz", "-f", "ref.fa", "-g", "anno.gff", "some.vcf"}); rc != 2 {
		t.Errorf("-C xyz must be rejected, got rc=%d", rc)
	}
	// -C l lists the tables and exits 0.
	if rc := runCSQ([]string{"-C", "l"}); rc != 0 {
		t.Errorf("-C l should list tables and exit 0, got rc=%d", rc)
	}
	// -b/--brief-predictions is now implemented: it passes the deferred
	// gate and fails later only because the input files are absent (rc=1),
	// not rc=2.
	if rc := runCSQ([]string{"-b", "-f", "ref.fa", "-g", "anno.gff", "some.vcf"}); rc != 1 {
		t.Errorf("-b should be accepted (file-open failure rc=1), got rc=%d", rc)
	}
	// -B/--trim-protein-seq with an explicit length is accepted (rc=1 from
	// the absent input files, not the rc=2 argument-rejection path).
	if rc := runCSQ([]string{"-B", "2", "-f", "ref.fa", "-g", "anno.gff", "some.vcf"}); rc != 1 {
		t.Errorf("-B 2 should be accepted (file-open failure rc=1), got rc=%d", rc)
	}
	// Upstream rejects -B < 1; we mirror that with an up-front rc=2 before
	// any file I/O.
	if rc := runCSQ([]string{"-B", "-1", "-f", "ref.fa", "-g", "anno.gff", "some.vcf"}); rc != 2 {
		t.Errorf("-B -1 must be rejected (rc=2), got rc=%d", rc)
	}
}
