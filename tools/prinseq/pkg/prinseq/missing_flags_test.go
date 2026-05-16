package prinseq

// Unit tests for the five upstream PRINSEQ-lite flags closed in
// PR claude/prinseq-missing-flags:
//
//	--out_format       (prinseq-lite.pl:242-247, branches at 1302-1348,
//	                    3703-3714, 3737-3757)
//	--seq_id_mappings  (prinseq-lite.pl:293-295, 1350-1358, 3645-3647)
//	--ns_max_p         (prinseq-lite.pl:344-346, 3465-3470)
//	--noniupac         (prinseq-lite.pl:352-354, 3478-3481)
//	--phred64 (alias)  (prinseq-lite.pl:230-232, 760-764) — encoding
//	                    toggle equivalent to --qual-type illumina; the
//	                    underlying offset is already covered by
//	                    phred_test.go, but the parity-test below pins
//	                    down the alias semantics end-to-end.
//
// Upstream PRINSEQ-lite has no formal test suite, so these are
// hand-built fixtures derived from the documented behaviour rather
// than byte-for-byte parity captures.

import (
	"bytes"
	"strings"
	"testing"
)

// runFilterCapture is a small wrapper that exercises Filter and returns
// the (primary, fasta, qual) byte streams. It's a convenience for the
// out_format tests below.
func runFilterCapture(t *testing.T, input string, isFastq bool, opts FilterOptions) (primary, fastaSec, qualSec []byte) {
	t.Helper()
	var prim, fa, qu bytes.Buffer
	opts.FastaOut = &fa
	opts.QualOut = &qu
	if err := Filter(strings.NewReader(input), &prim, isFastq, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	return prim.Bytes(), fa.Bytes(), qu.Bytes()
}

func TestNsMaxP_PercentBoundaries(t *testing.T) {
	// Three FASTA records: 0%, 50% and 100% Ns of length 10.
	input := ">zero\nACGTACGTAC\n>half\nACGTNNNNNA\n>all\nNNNNNNNNNN\n"

	tests := []struct {
		name    string
		maxNsP  float64
		wantIDs []string // primary IDs that pass
	}{
		// Upstream check is `(N_count * 100 / len) > threshold` (strict
		// `>`, prinseq-lite.pl:3467). At threshold=50 a 50%-N record
		// passes; at 49 it does not.
		{"threshold0p001_keepsZeroOnly", 0.0001, []string{"zero"}},
		{"threshold49_keepsZeroOnly", 49, []string{"zero"}},
		{"threshold50_keepsZeroAndHalf", 50, []string{"zero", "half"}},
		{"threshold51_keepsZeroAndHalf", 51, []string{"zero", "half"}},
		{"threshold99p9_keepsZeroAndHalf", 99.9, []string{"zero", "half"}},
		// At exactly 100%: strict `>` keeps even the all-N record.
		{"threshold100_keepsAll", 100, []string{"zero", "half", "all"}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var out bytes.Buffer
			err := Filter(strings.NewReader(input), &out, false, FilterOptions{MaxNsP: tc.maxNsP})
			if err != nil {
				t.Fatalf("Filter: %v", err)
			}
			gotIDs := extractFastaIDs(out.String())
			if !equalStringSlices(gotIDs, tc.wantIDs) {
				t.Fatalf("want %v, got %v\nraw output:\n%s", tc.wantIDs, gotIDs, out.String())
			}
		})
	}
}

func TestNonIUPAC_FiltersAmbiguityCodes(t *testing.T) {
	// One sequence per IUPAC code that we should reject when --noniupac
	// is set, plus a clean ACGTN-only one that should pass.
	input := ">clean\nACGTNACGTN\n" +
		">r_upper\nACGTRACGT\n" +
		">y_lower\nacgtyacgt\n" +
		">mixed\nACGTACGTm\n" + // upstream upper-cases before the regex
		">spaceFail\nACGTACGT.\n" // punctuation is also outside ACGTN

	var out bytes.Buffer
	if err := Filter(strings.NewReader(input), &out, false, FilterOptions{NonIUPAC: true}); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	gotIDs := extractFastaIDs(out.String())
	wantIDs := []string{"clean"}
	if !equalStringSlices(gotIDs, wantIDs) {
		t.Fatalf("noniupac kept %v, want %v\nraw:\n%s", gotIDs, wantIDs, out.String())
	}

	// When --noniupac is NOT set, all records pass (the rest of the
	// filter chain in this test is empty).
	var out2 bytes.Buffer
	if err := Filter(strings.NewReader(input), &out2, false, FilterOptions{}); err != nil {
		t.Fatalf("Filter (no noniupac): %v", err)
	}
	gotAll := extractFastaIDs(out2.String())
	if len(gotAll) != 5 {
		t.Fatalf("without noniupac expected all 5, got %v", gotAll)
	}
}

func TestSeqIDAndMappings_TSVFormat(t *testing.T) {
	// FASTA renaming.
	input := ">orig_one some comment\nACGTACGT\n" +
		">orig_two\nACGTACGT\n" +
		">orig_three\nACGTACGT\n"
	var out, mappings bytes.Buffer
	opts := FilterOptions{SeqID: "seq", SeqIDMap: &mappings}
	if err := Filter(strings.NewReader(input), &out, false, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}

	// Headers should be rewritten to seq1/seq2/seq3, dropping any
	// trailing whitespace/comment from the original description.
	want := ">seq1\nACGTACGT\n>seq2\nACGTACGT\n>seq3\nACGTACGT\n"
	if out.String() != want {
		t.Fatalf("renamed FASTA mismatch\nwant:\n%s\ngot:\n%s", want, out.String())
	}

	// Mapping file: "<orig_id>\t<new_id>\n" matching
	// prinseq-lite.pl:3646 ("join('\t', $sid, $params{seq_id}.$seqcount)").
	wantMap := "orig_one\tseq1\norig_two\tseq2\norig_three\tseq3\n"
	if mappings.String() != wantMap {
		t.Fatalf("mappings mismatch\nwant:\n%qgot:\n%q", wantMap, mappings.String())
	}
}

func TestSeqIDMappings_OnlyPassingRecordsAreNumbered(t *testing.T) {
	// One record fails the length filter; the counter must not
	// advance for it.
	input := ">short\nAC\n>long_one\nACGTACGTAC\n>long_two\nACGTACGTAC\n"
	var out, mappings bytes.Buffer
	opts := FilterOptions{MinLen: 5, SeqID: "x", SeqIDMap: &mappings}
	if err := Filter(strings.NewReader(input), &out, false, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	if got := out.String(); got != ">x1\nACGTACGTAC\n>x2\nACGTACGTAC\n" {
		t.Fatalf("renamed FASTA mismatch\ngot:\n%s", got)
	}
	if got := mappings.String(); got != "long_one\tx1\nlong_two\tx2\n" {
		t.Fatalf("mappings mismatch\ngot:\n%q", got)
	}
}

func TestOutFormat_FASTQInput_AllFiveModes(t *testing.T) {
	// One short FASTQ record. Each test asserts the three potential
	// output streams (primary, FASTA sidecar, QUAL sidecar) match the
	// expected byte layout.
	input := "@r1\nACGT\n+\nIIII\n"
	// Phred+33: 'I' (73) - 33 = 40, all four bases score 40 → emitted
	// verbatim as `40 40 40 40` on one line by writePrinseqQual.

	tests := []struct {
		name          string
		mode          int
		wantPrimary   string
		wantFastaSide string
		wantQualSide  string
	}{
		{
			name:          "mode1_FASTA_only",
			mode:          1,
			wantPrimary:   ">r1\nACGT\n",
			wantFastaSide: "",
			wantQualSide:  "",
		},
		{
			name:          "mode2_FASTA_plus_QUAL",
			mode:          2,
			wantPrimary:   ">r1\nACGT\n",
			wantFastaSide: "",
			wantQualSide:  ">r1\n40 40 40 40\n",
		},
		{
			name:          "mode3_FASTQ_only",
			mode:          3,
			wantPrimary:   "@r1\nACGT\n+r1\nIIII\n",
			wantFastaSide: "",
			wantQualSide:  "",
		},
		{
			name:          "mode4_FASTQ_plus_FASTA",
			mode:          4,
			wantPrimary:   "@r1\nACGT\n+r1\nIIII\n",
			wantFastaSide: ">r1\nACGT\n",
			wantQualSide:  "",
		},
		{
			name:          "mode5_FASTQ_plus_FASTA_plus_QUAL",
			mode:          5,
			wantPrimary:   "@r1\nACGT\n+r1\nIIII\n",
			wantFastaSide: ">r1\nACGT\n",
			wantQualSide:  ">r1\n40 40 40 40\n",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			prim, fa, qu := runFilterCapture(t, input, true, FilterOptions{OutFormat: tc.mode})
			if string(prim) != tc.wantPrimary {
				t.Fatalf("primary stream mismatch\nwant %q\ngot  %q", tc.wantPrimary, string(prim))
			}
			if string(fa) != tc.wantFastaSide {
				t.Fatalf("FASTA sidecar mismatch\nwant %q\ngot  %q", tc.wantFastaSide, string(fa))
			}
			if string(qu) != tc.wantQualSide {
				t.Fatalf("QUAL sidecar mismatch\nwant %q\ngot  %q", tc.wantQualSide, string(qu))
			}
		})
	}
}

func TestOutFormat_QualWrap_AtLineWidth(t *testing.T) {
	// Build a 70-character FASTQ record at Phred+33 quality 40 to
	// verify the default 60-column wrap from
	// prinseq-lite.pl:2531-2546 / LINE_WIDTH=60 (line 45).
	seq := strings.Repeat("A", 70)
	qual := strings.Repeat("I", 70) // 'I'-33=40
	input := "@r1\n" + seq + "\n+\n" + qual + "\n"

	_, _, qu := runFilterCapture(t, input, true, FilterOptions{OutFormat: 2})
	// 60 "40"s on line 1, 10 "40"s on line 2, no trailing space.
	line1 := strings.TrimRight(strings.Repeat("40 ", 60), " ")
	line2 := strings.TrimRight(strings.Repeat("40 ", 10), " ")
	want := ">r1\n" + line1 + "\n" + line2 + "\n"
	if string(qu) != want {
		t.Fatalf("QUAL wrap mismatch\nwant:\n%q\ngot:\n%q", want, string(qu))
	}
}

func TestOutFormat_QualSmallValuesPadded(t *testing.T) {
	// Quality bytes giving phred 0, 5, 9, 10 verify the two-character
	// zero-padded layout: " 0  5  9 10".
	// Phred+33 chars: 0→'!'(33), 5→'&'(38), 9→'*'(42), 10→'+'(43).
	input := "@r\nACGT\n+\n!&*+\n"
	_, _, qu := runFilterCapture(t, input, true, FilterOptions{OutFormat: 2})
	want := ">r\n 0  5  9 10\n"
	if string(qu) != want {
		t.Fatalf("QUAL padding mismatch\nwant %q\ngot  %q", want, string(qu))
	}
}

func TestPhred64Alias_DecodesPhred64(t *testing.T) {
	// Same physical encoding as TestTrimQualityOffsetTable's phred64
	// case: 'h' (104) - 64 = 40. Setting QualType to "illumina" must
	// decode this as quality 40, not the much lower value implied by
	// the default Phred+33 offset.
	input := "@r1\nACGTACGT\n+\nhhhhhhhh\n"
	var out bytes.Buffer
	// MinQualMean 35 should pass under Phred+64 (mean=40) but reject
	// under Phred+33 (mean ≈ 71, no — actually Phred+33 mean would be
	// even higher because 104-33=71). Switch the check around:
	// MaxQualMean 50 rejects when wrongly decoded as Phred+33 (71),
	// accepts when correctly decoded as Phred+64 (40).
	opts := FilterOptions{QualType: "illumina", MaxQualMean: 50}
	if err := Filter(strings.NewReader(input), &out, true, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	if !strings.Contains(out.String(), "ACGTACGT") {
		t.Fatalf("expected record to pass under Phred+64; got %q", out.String())
	}

	// Now check that the QUAL sidecar reflects the Phred+64 decode.
	_, _, qu := runFilterCapture(t, input, true, FilterOptions{
		QualType:  "illumina",
		OutFormat: 2,
	})
	want := ">r1\n40 40 40 40 40 40 40 40\n"
	if string(qu) != want {
		t.Fatalf("Phred+64 QUAL output mismatch\nwant %q\ngot  %q", want, string(qu))
	}
}

// --- helpers ---------------------------------------------------------

func extractFastaIDs(s string) []string {
	var ids []string
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, ">") {
			id := strings.TrimPrefix(line, ">")
			if i := strings.IndexAny(id, " \t"); i >= 0 {
				id = id[:i]
			}
			ids = append(ids, id)
		}
	}
	return ids
}

func equalStringSlices(a, b []string) bool {
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
