package bedbamtobed

import (
	"strings"
	"testing"
)

// samHeader is a minimal two-ref SAM header reused by the unit tests below.
const samHeader = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n"

func runSAM(t *testing.T, body string, opts Options) string {
	t.Helper()
	var out strings.Builder
	if _, err := Run(strings.NewReader(samHeader+body), &out, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return out.String()
}

func TestBED6_Default(t *testing.T) {
	body := "r1\t0\tchr1\t11\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	got := runSAM(t, body, Options{})
	want := "chr1\t10\t20\tr1\t60\t+\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestBED6_ReverseAndMateSuffix(t *testing.T) {
	// flag 0x10 (reverse) | 0x1 (paired) | 0x80 (read2) = 0x91 = 145.
	body := "r1\t145\tchr1\t11\t30\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	got := runSAM(t, body, Options{})
	want := "chr1\t10\t20\tr1/2\t30\t-\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSplit_OnN(t *testing.T) {
	body := "r1\t0\tchr1\t1\t40\t5M10N5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	got := runSAM(t, body, Options{ObeySplits: true})
	want := "chr1\t0\t5\tr1\t40\t+\nchr1\t15\t20\tr1\t40\t+\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSplitD_BreaksOnDeletion(t *testing.T) {
	body := "r1\t0\tchr1\t1\t40\t5M2D5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	withN := runSAM(t, body, Options{ObeySplits: true})
	if withN != "chr1\t0\t12\tr1\t40\t+\n" {
		t.Errorf("-split should keep D inside a single block, got %q", withN)
	}
	withD := runSAM(t, body, Options{ObeySplits: true, SplitOnDeletions: true})
	want := "chr1\t0\t5\tr1\t40\t+\nchr1\t7\t12\tr1\t40\t+\n"
	if withD != want {
		t.Errorf("-splitD got %q want %q", withD, want)
	}
}

func TestCigarColumn(t *testing.T) {
	body := "r1\t0\tchr1\t1\t40\t5M2I3M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	got := runSAM(t, body, Options{UseCigar: true})
	want := "chr1\t0\t8\tr1\t40\t+\t5M2I3M\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestTagScore(t *testing.T) {
	body := "r1\t0\tchr1\t1\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:3\n"
	got := runSAM(t, body, Options{Tag: "NM"})
	want := "chr1\t0\t10\tr1\t3\t+\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestTagMissing_IsError(t *testing.T) {
	body := "r1\t0\tchr1\t1\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	var out strings.Builder
	_, err := Run(strings.NewReader(samHeader+body), &out, Options{Tag: "ZZ"})
	if err == nil {
		t.Fatalf("expected error for missing tag")
	}
	if !strings.Contains(err.Error(), "ZZ") {
		t.Errorf("error should name the tag: %v", err)
	}
}

func TestUnmapped_Skipped(t *testing.T) {
	body := "u\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\nr1\t0\tchr1\t1\t40\t4M\t*\t0\t0\tACGT\tIIII\n"
	got := runSAM(t, body, Options{})
	if got != "chr1\t0\t4\tr1\t40\t+\n" {
		t.Errorf("unmapped read should be skipped, got %q", got)
	}
}

func TestValidate_Rejections(t *testing.T) {
	tests := []struct {
		name string
		o    Options
		want string
	}{
		{"ed+split", Options{UseEditDistance: true, ObeySplits: true}, "-ed with -splits"},
		{"ed+cigar", Options{UseEditDistance: true, UseCigar: true}, "-cigar"},
		{"tag+bedpe", Options{WriteBedPE: true, Tag: "NM"}, "-tag with -bedpe"},
		{"mate1 without bedpe", Options{Mate1First: true}, "-mate1 with -bedpe"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := tt.o.Validate()
			if err == nil || !strings.Contains(err.Error(), tt.want) {
				t.Errorf("Validate()=%v, want substring %q", err, tt.want)
			}
		})
	}
}

func TestBedPE_MinMapqAndOrder(t *testing.T) {
	body := "" +
		"pA\t99\tchr1\t100\t40\t30M\t=\t300\t230\t*\t*\tNM:i:1\n" +
		"pA\t147\tchr1\t300\t30\t30M\t=\t100\t-230\t*\t*\tNM:i:2\n"
	got := runSAM(t, body, Options{WriteBedPE: true})
	want := "chr1\t99\t129\tchr1\t299\t329\tpA\t30\t+\t-\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
	gotED := runSAM(t, body, Options{WriteBedPE: true, UseEditDistance: true})
	wantED := "chr1\t99\t129\tchr1\t299\t329\tpA\t3\t+\t-\n"
	if gotED != wantED {
		t.Errorf("ed: got %q want %q", gotED, wantED)
	}
}
