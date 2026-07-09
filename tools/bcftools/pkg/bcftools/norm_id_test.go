package bcftools

import "testing"

// TestBcfAddIDJoin exercises bcfAddID and the norm -m+ ID-join sequence
// (seed with records[0].ID verbatim, then fold in each later record's whole ID
// via bcfAddID). Expected results mirror upstream htslib bcf_add_id
// (vcf.c:6004) + vcfnorm.c merge_lines: no per-token splitting, no dedup of the
// first record's own tokens, only skip a later ID when the whole string is
// already present as a ';'-delimited token.
func TestBcfAddIDJoin(t *testing.T) {
	// joinIDs replicates mergeBiallelicsToMultiallelic's ID logic for a list of
	// per-record ID strings.
	joinIDs := func(ids []string) string {
		dst := ids[0]
		for i := 1; i < len(ids); i++ {
			if ids[i] == "" || ids[i] == "." {
				continue
			}
			bcfAddID(&dst, ids[i])
		}
		return dst
	}

	cases := []struct {
		name string
		ids  []string
		want string
	}{
		{"dot-then-two", []string{".", "rs45", "rs46"}, "rs45;rs46"},
		{"overlap-token", []string{"rs1;rs2", "rs2;rs3"}, "rs1;rs2;rs2;rs3"},
		{"whole-not-token", []string{"rs1", "rs1;rs4"}, "rs1;rs1;rs4"},
		{"already-present", []string{"rs1;rs2", "rs2"}, "rs1;rs2"},
		{"internal-dup-kept", []string{"rs1;rs1", "rs2"}, "rs1;rs1;rs2"},
		{"substring-not-token", []string{"foo", "foobar"}, "foo;foobar"},
		{"all-dot", []string{".", ".", "."}, "."},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := joinIDs(tc.ids); got != tc.want {
				t.Errorf("joinIDs(%v) = %q, want %q", tc.ids, got, tc.want)
			}
		})
	}
}
