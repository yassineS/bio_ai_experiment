package difffuzz

import "testing"

func TestDiffClassifier(t *testing.T) {
	tests := []struct {
		name string
		ours RunOutcome
		up   RunOutcome
		want DivergenceClass
	}{
		{
			name: "identical clean",
			ours: RunOutcome{Stdout: []byte("hello\n"), ExitCode: 0},
			up:   RunOutcome{Stdout: []byte("hello\n"), ExitCode: 0},
			want: ClassNone,
		},
		{
			name: "stdout differs",
			ours: RunOutcome{Stdout: []byte("hello\n"), ExitCode: 0},
			up:   RunOutcome{Stdout: []byte("world\n"), ExitCode: 0},
			want: ClassStdoutDiffers,
		},
		{
			name: "stderr differs only",
			ours: RunOutcome{Stdout: []byte("x\n"), Stderr: []byte("warn A\n"), ExitCode: 0},
			up:   RunOutcome{Stdout: []byte("x\n"), Stderr: []byte("warn B\n"), ExitCode: 0},
			want: ClassStderrDiffers,
		},
		{
			name: "exit code differs",
			ours: RunOutcome{Stdout: []byte("x\n"), ExitCode: 0},
			up:   RunOutcome{Stdout: []byte("x\n"), ExitCode: 1},
			want: ClassExitDiffers,
		},
		{
			name: "exit code outranks stdout",
			ours: RunOutcome{Stdout: []byte("a\n"), ExitCode: 0},
			up:   RunOutcome{Stdout: []byte("b\n"), ExitCode: 2},
			want: ClassExitDiffers,
		},
		{
			name: "one crashed",
			ours: RunOutcome{Crashed: true, ExitCode: -1},
			up:   RunOutcome{Stdout: []byte("x\n"), ExitCode: 1},
			want: ClassOneCrashed,
		},
		{
			name: "one timed out is one-crashed",
			ours: RunOutcome{TimedOut: true},
			up:   RunOutcome{ExitCode: 0},
			want: ClassOneCrashed,
		},
		{
			name: "both crashed",
			ours: RunOutcome{Crashed: true, ExitCode: -1},
			up:   RunOutcome{Crashed: true, ExitCode: -1},
			want: ClassBothCrashed,
		},
		{
			name: "stdout outranks stderr",
			ours: RunOutcome{Stdout: []byte("a\n"), Stderr: []byte("p\n"), ExitCode: 0},
			up:   RunOutcome{Stdout: []byte("b\n"), Stderr: []byte("q\n"), ExitCode: 0},
			want: ClassStdoutDiffers,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, _ := Diff(tc.ours, tc.up)
			if got != tc.want {
				t.Fatalf("Diff()=%q want %q", got, tc.want)
			}
		})
	}
}

// TestDiffReusesProvenanceNormalization verifies that a pure provenance
// difference (a @PG line, a ##bcftools_*Command line) is NOT a divergence,
// because Diff applies the same StripProvenance the parity harness uses.
func TestDiffReusesProvenanceNormalization(t *testing.T) {
	ours := RunOutcome{Stdout: []byte(
		"@HD\tVN:1.6\n@PG\tID:ours\tPN:ourtool\tVN:9.9\n@SQ\tSN:chr1\tLN:100\n"), ExitCode: 0}
	up := RunOutcome{Stdout: []byte(
		"@HD\tVN:1.6\n@PG\tID:up\tPN:samtools\tVN:1.20\n@SQ\tSN:chr1\tLN:100\n"), ExitCode: 0}
	if got, detail := Diff(ours, up); got != ClassNone {
		t.Fatalf("provenance-only difference classified as %q (%s); want none", got, detail)
	}

	// A VCF command-header line is also benign.
	ours = RunOutcome{Stdout: []byte("##fileformat=VCFv4.2\n##bcftools_viewCommand=view a.vcf\n#CHROM\n"), ExitCode: 0}
	up = RunOutcome{Stdout: []byte("##fileformat=VCFv4.2\n##bcftools_viewCommand=view b.vcf; Date=...\n#CHROM\n"), ExitCode: 0}
	if got, _ := Diff(ours, up); got != ClassNone {
		t.Fatalf("VCF command-header difference classified as %q; want none", got)
	}

	// But a genuine DATA difference still diverges.
	ours = RunOutcome{Stdout: []byte("@SQ\tSN:chr1\tLN:100\n"), ExitCode: 0}
	up = RunOutcome{Stdout: []byte("@SQ\tSN:chr1\tLN:200\n"), ExitCode: 0}
	if got, _ := Diff(ours, up); got != ClassStdoutDiffers {
		t.Fatalf("genuine data difference classified as %q; want stdout-differs", got)
	}
}

func TestIsDivergence(t *testing.T) {
	if IsDivergence(ClassNone) {
		t.Error("ClassNone should not be a divergence")
	}
	for _, c := range []DivergenceClass{ClassStdoutDiffers, ClassStderrDiffers, ClassExitDiffers, ClassOneCrashed, ClassBothCrashed} {
		if !IsDivergence(c) {
			t.Errorf("%q should be a divergence", c)
		}
	}
}
