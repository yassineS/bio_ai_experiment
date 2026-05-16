package seqtk

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// TestFamask_BasicRules exercises the three masking rules (X = keep,
// x = soft-mask, anything else = overwrite) on a hand-built fixture
// short enough to read at a glance.
func TestFamask_BasicRules(t *testing.T) {
	cases := []struct {
		name     string
		src      string
		mask     string
		want     string
		wantWarn string // substring expected on the warn stream
	}{
		{
			name: "all keep",
			src:  ">a\nACGT\n",
			mask: ">a\nXXXX\n",
			want: ">a\nACGT\n",
		},
		{
			name: "all softmask",
			src:  ">a\nACGT\n",
			mask: ">a\nxxxx\n",
			want: ">a\nacgt\n",
		},
		{
			name: "all overwrite",
			src:  ">a\nACGT\n",
			mask: ">a\nNNNN\n",
			want: ">a\nNNNN\n",
		},
		{
			name: "mixed",
			src:  ">a\nACGTACGT\n",
			mask: ">a\nXXxxNNXX\n",
			want: ">a\nACgtNNGT\n",
		},
		{
			name:     "name mismatch warns",
			src:      ">a\nACGT\n",
			mask:     ">b\nXXXX\n",
			want:     ">a\nACGT\n",
			wantWarn: "Different sequence names",
		},
		{
			name:     "length mismatch warns and uses shorter",
			src:      ">a\nACGT\n",
			mask:     ">a\nXX\n",
			want:     ">a\nAC\n",
			wantWarn: "Unequal sequence length",
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var out, warn bytes.Buffer
			if err := famaskImpl(strings.NewReader(tc.src), strings.NewReader(tc.mask), &out, &warn); err != nil {
				t.Fatalf("famaskImpl: %v", err)
			}
			if got := out.String(); got != tc.want {
				t.Errorf("output mismatch\nwant: %q\ngot:  %q", tc.want, got)
			}
			if tc.wantWarn != "" && !strings.Contains(warn.String(), tc.wantWarn) {
				t.Errorf("expected warn substring %q in: %q", tc.wantWarn, warn.String())
			}
		})
	}
}

// TestFamask_LineWrapAt60 verifies that long sequences are wrapped at
// 60 bases per line, matching upstream's `l%60==0` rule.
func TestFamask_LineWrapAt60(t *testing.T) {
	src := ">x\n" + strings.Repeat("A", 125) + "\n"
	mask := ">x\n" + strings.Repeat("X", 125) + "\n"
	var out bytes.Buffer
	if err := famaskImpl(strings.NewReader(src), strings.NewReader(mask), &out, io.Discard); err != nil {
		t.Fatalf("famaskImpl: %v", err)
	}
	wantLines := []string{
		">x",
		strings.Repeat("A", 60),
		strings.Repeat("A", 60),
		strings.Repeat("A", 5),
	}
	want := strings.Join(wantLines, "\n") + "\n"
	if got := out.String(); got != want {
		t.Errorf("wrap mismatch\nwant: %q\ngot:  %q", want, got)
	}
}

// TestFamask_MaskNonASCIIPasses confirms the toLowerByte helper only
// touches A-Z (matches C's tolower for the ASCII subset).
func TestFamask_PassThroughLowercaseSrc(t *testing.T) {
	src := ">a\nacgt\n"
	mask := ">a\nxxxx\n"
	var out bytes.Buffer
	if err := famaskImpl(strings.NewReader(src), strings.NewReader(mask), &out, io.Discard); err != nil {
		t.Fatalf("famaskImpl: %v", err)
	}
	// 'a' is already lowercase; toLowerByte is a no-op.
	if got, want := out.String(), ">a\nacgt\n"; got != want {
		t.Errorf("want %q, got %q", want, got)
	}
}
