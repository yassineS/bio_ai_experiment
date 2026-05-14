package prinseq

import (
	"bytes"
	"strings"
	"testing"
)

// TestTrimQualityOffsetTable exercises trimQualityLeft / trimQualityRight with
// both Phred+33 and Phred+64 encodings to verify the ASCII offset is honoured.
//
// Historical bug: both helpers subtracted 33 unconditionally, so when callers
// passed --qual-type illumina (Phred+64) the decoded score was off by 31 and
// near-threshold bases were either kept when they should have been trimmed or
// vice versa. The Phred+64 cases below are constructed so that the buggy
// implementation would NOT trim (decoded scores were inflated by 31) while the
// fixed implementation trims as expected.
func TestTrimQualityOffsetTable(t *testing.T) {
	tests := []struct {
		name      string
		seq       string
		qual      string
		threshold int
		offset    int
		left      bool // true => trimQualityLeft; false => trimQualityRight
		wantSeq   string
		wantQual  string
	}{
		{
			// Phred+33: '!' = 0, '5' = 20. Threshold 20 trims the leading '!' run.
			name:      "phred33_left_threshold20",
			seq:       "AAAACCCC",
			qual:      "!!!!5555",
			threshold: 20,
			offset:    33,
			left:      true,
			wantSeq:   "CCCC",
			wantQual:  "5555",
		},
		{
			// Phred+33: trailing '!' run trimmed from the right at threshold 20.
			name:      "phred33_right_threshold20",
			seq:       "CCCCAAAA",
			qual:      "5555!!!!",
			threshold: 20,
			offset:    33,
			left:      false,
			wantSeq:   "CCCC",
			wantQual:  "5555",
		},
		{
			// Phred+64: 'E' decodes to 5 (64+5=69), '^' decodes to 30 (64+30=94).
			// At threshold 20 the leading 'E' run should be trimmed.
			// Under the bug (offset 33), 'E' would decode to 36 — well above
			// threshold — and nothing would be trimmed, so this case fails today.
			name:      "phred64_left_threshold20",
			seq:       "AAAACCCC",
			qual:      "EEEE^^^^",
			threshold: 20,
			offset:    64,
			left:      true,
			wantSeq:   "CCCC",
			wantQual:  "^^^^",
		},
		{
			// Phred+64: same characters, reversed orientation — the trailing
			// 'E' run is trimmed from the right at threshold 20. Buggy code
			// would leave the string untouched.
			name:      "phred64_right_threshold20",
			seq:       "CCCCAAAA",
			qual:      "^^^^EEEE",
			threshold: 20,
			offset:    64,
			left:      false,
			wantSeq:   "CCCC",
			wantQual:  "^^^^",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			var gotSeq, gotQual []byte
			if tt.left {
				gotSeq, gotQual = trimQualityLeft([]byte(tt.seq), []byte(tt.qual), tt.threshold, tt.offset)
			} else {
				gotSeq, gotQual = trimQualityRight([]byte(tt.seq), []byte(tt.qual), tt.threshold, tt.offset)
			}
			if string(gotSeq) != tt.wantSeq {
				t.Errorf("seq: got %q, want %q", string(gotSeq), tt.wantSeq)
			}
			if string(gotQual) != tt.wantQual {
				t.Errorf("qual: got %q, want %q", string(gotQual), tt.wantQual)
			}
		})
	}
}

// TestPhredOffsetHelper documents the qual-type -> offset mapping.
func TestPhredOffsetHelper(t *testing.T) {
	cases := map[string]int{
		"":         33, // default
		"sanger":   33,
		"illumina": 64,
		"unknown":  33, // unknown strings fall back to Sanger
	}
	for in, want := range cases {
		if got := phredOffset(in); got != want {
			t.Errorf("phredOffset(%q) = %d, want %d", in, got, want)
		}
	}
}

// TestFilterFastqRespectsQualType verifies the end-to-end fix: running Filter
// over the same FASTQ payload with --qual-type illumina vs sanger produces
// different trimming results, as it should. Before the fix both invocations
// produced identical output because trimQuality* hard-coded offset 33.
func TestFilterFastqRespectsQualType(t *testing.T) {
	// 'E' (ASCII 69) is Phred+64 score 5 — should be trimmed at threshold 20
	// under illumina, but kept under sanger (decoded score 36 there).
	input := "@r1\nAAAACCCC\n+\nEEEE^^^^\n"

	runFilter := func(qualType string) string {
		var out bytes.Buffer
		err := Filter(strings.NewReader(input), &out, true, FilterOptions{
			TrimQualL: 20,
			QualType:  qualType,
		})
		if err != nil {
			t.Fatalf("Filter(%q) error: %v", qualType, err)
		}
		return out.String()
	}

	sanger := runFilter("sanger")
	illumina := runFilter("illumina")

	if sanger == illumina {
		t.Fatalf("expected different trimming for sanger vs illumina;\n sanger  = %q\n illumina= %q", sanger, illumina)
	}

	// Under illumina, the leading AAAA should be trimmed away.
	if !strings.Contains(illumina, "\nCCCC\n") {
		t.Errorf("illumina output should retain trimmed CCCC sequence, got %q", illumina)
	}
	// Under sanger, the sequence is left alone (all chars decode well above 20).
	if !strings.Contains(sanger, "\nAAAACCCC\n") {
		t.Errorf("sanger output should retain full AAAACCCC sequence, got %q", sanger)
	}
}
