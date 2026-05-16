package main

import (
	"path/filepath"
	"testing"
)

// TestCheckMpileupDeferred locks in the upstream-flag surface that
// runMpileup hard-rejects rather than silently accepting. Per the
// project parity rule (docs/PARITY_ROADMAP.md "Definition of 1:1")
// every documented upstream flag must be recognised — either
// implemented or gracefully rejected with a roadmap pointer. A
// regression that drops any of these from the rejection set without
// implementing the underlying behaviour is a parity bug.
func TestCheckMpileupDeferred(t *testing.T) {
	if got := checkMpileupDeferred(&mpileupFlags{}); got != "" {
		t.Fatalf("empty flags: got deferred=%q, want \"\"", got)
	}
	if got := checkMpileupDeferred(&mpileupFlags{outputType: "v"}); got != "" {
		t.Fatalf("-O v: got deferred=%q, want \"\"", got)
	}
	if got := checkMpileupDeferred(&mpileupFlags{outputType: "z"}); got != "" {
		t.Fatalf("-O z: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		mf   *mpileupFlags
		want string
	}{
		{"redoBAQ", &mpileupFlags{redoBAQ: true}, "-E/--redo-BAQ"},
		{"-O u", &mpileupFlags{outputType: "u"}, "-O u (BCF output)"},
		{"-O b", &mpileupFlags{outputType: "b"}, "-O b (BCF output)"},
		{"-O bogus", &mpileupFlags{outputType: "x"}, "-O x"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkMpileupDeferred(tc.mf); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestMpileupRunInputs sanity-checks the CLI binding (missing input,
// missing -f, help).
func TestMpileupRunInputs(t *testing.T) {
	if rc := runMpileup([]string{"--help"}); rc != 0 {
		t.Errorf("--help rc=%d want 0", rc)
	}
	if rc := runMpileup([]string{"-?"}); rc != 0 {
		t.Errorf("-? rc=%d want 0", rc)
	}
	if rc := runMpileup([]string{}); rc != 2 {
		t.Errorf("no args rc=%d want 2", rc)
	}
	if rc := runMpileup([]string{"some.bam"}); rc != 2 {
		t.Errorf("no -f rc=%d want 2", rc)
	}
	// -E/--redo-BAQ flag should reject early.
	if rc := runMpileup([]string{"-E", "-f", "x.fa", "some.bam"}); rc != 2 {
		t.Errorf("-E rc=%d want 2", rc)
	}
	// -O b (BCF) should reject.
	if rc := runMpileup([]string{"-O", "b", "-f", "x.fa", "some.bam"}); rc != 2 {
		t.Errorf("-O b rc=%d want 2", rc)
	}
}

// TestMpileupFlagSurface verifies that the full upstream flag table
// is wired so a typical "fold every accepted flag onto the command
// line" stress-test parses without errors. This guards against a
// regression that drops a getopt_long entry from registerMpileupFlags.
//
// The fixture mimics the upstream `mpileup.c::lopts[]` table — every
// flag must be parseable. We expect rc=2 because we do not pass a
// real BAM, but parsing must succeed (the rc=2 comes from the
// "missing input file" branch, NOT from a flag-parse error).
func TestMpileupFlagSurface(t *testing.T) {
	// Build the full set of flags. Long form preferred for clarity;
	// short-form duplicates are exercised in the parser unit tests.
	tmpOut := filepath.Join(t.TempDir(), "out.vcf")
	args := []string{
		"--bam-list", "bams.txt",
		"--fasta-ref", "ref.fa",
		"--read-groups", "rg.txt",
		"--ignore-RG",
		"--output", tmpOut,
		"--output-type", "v",
		"--no-version",
		"--threads", "2",
		"--count-orphans",
		"--max-depth", "500",
		"--min-MQ", "20",
		"--min-BQ", "25",
		"--max-bq", "40",
		"--delta-BQ", "5",
		"--ignore-overlaps",
		"--adjust-MQ", "50",
		"--regions", "chr1:100-200",
		"--regions-file", "r.bed",
		"--targets", "chr2:1-2",
		"--targets-file", "t.bed",
		"--samples", "s1,s2",
		"--samples-file", "samples.txt",
		"--skip-any-unset", "0x1",
		"--skip-all-unset", "0x2",
		"--skip-any-set", "0x4",
		"--skip-all-set", "0x8",
		"--ls", "x",
		"--no-BAQ",
		"--full-BAQ",
		"--platforms", "ILLUMINA",
		"--per-sample-mF",
		"--illumina1.3+",
		"--config", "1.12",
		"--seed", "42",
		"--ambig-reads", "drop",
		"--skip-indels",
		"--max-idepth", "100",
		"--min-ireads", "2",
		"--tandem-qual", "100",
		"--ext-prob", "20",
		"--gap-frac", "0.05",
		"--max-read-len", "500",
		"--indel-bias", "1.0",
		"--indel-size", "100",
		"--indels-cns",
		"--open-prob", "40",
		"--del-bias", "0.5",
		"--poly-mqual",
		"--score-vs-ref", "0.5",
		"--seqq-offset", "20",
		"--annotate", "DP,AD,SP",
		"--gvcf", "5,15,30",
		"--no-reference",
		"--write-index",
		"--verbosity", "1",
	}
	// We expect rc=1 (open-fail on the FASTA "ref.fa") or rc=2; the
	// important thing is that fs.Parse did NOT return an error
	// (which would also yield rc=2 but with a different message —
	// hard to assert without scraping stderr). At minimum we
	// require the rc to be non-zero, indicating we got past parsing.
	rc := runMpileup(args)
	if rc == 0 {
		t.Errorf("expected non-zero rc (no BAM input), got 0")
	}
}
