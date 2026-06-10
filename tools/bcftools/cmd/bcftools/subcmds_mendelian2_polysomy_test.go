package main

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

// TestParseWriteIndexFormat checks that the optional -W/--write-index
// argument maps to the right index flavour. The bare flag (and an
// explicit "csi") yields a forced .csi index; "tbi" selects the tabix
// flavour. These were previously deferred; they are now wired so the
// mapping is the parity surface worth locking in.
func TestParseWriteIndexFormat(t *testing.T) {
	cases := []struct {
		name string
		arg  string
		want bcftools.IndexFormat
	}{
		{"default", "", bcftools.IndexCSI},
		{"csi", "csi", bcftools.IndexCSI},
		{"csi-equals", "=csi", bcftools.IndexCSI},
		{"tbi", "tbi", bcftools.IndexTBI},
		{"tbi-equals", "=tbi", bcftools.IndexTBI},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := parseWriteIndexFormat(tc.arg, "out.vcf.gz")
			if got.Format != tc.want {
				t.Errorf("parseWriteIndexFormat(%q): got %v, want %v", tc.arg, got.Format, tc.want)
			}
			if !got.Force {
				t.Errorf("parseWriteIndexFormat(%q): Force should be true", tc.arg)
			}
		})
	}
}
