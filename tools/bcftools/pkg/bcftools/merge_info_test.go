package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// runMergeErr runs Merge over in-memory inputs and returns the error (if any)
// without failing the test, for asserting rule-validation errors.
func runMergeErr(t *testing.T, inputs []string, opts MergeOptions) (string, error) {
	t.Helper()
	dir := t.TempDir()
	paths := make([]string, len(inputs))
	for i, s := range inputs {
		p := filepath.Join(dir, "in"+string(rune('0'+i))+".vcf")
		if err := os.WriteFile(p, []byte(s), 0644); err != nil {
			t.Fatal(err)
		}
		paths[i] = p
	}
	var out bytes.Buffer
	_, err := MergeFiles(paths, &out, opts)
	return out.String(), err
}

// mergeInfoVCF builds a single-sample VCF with a rich INFO header (DP scalar,
// DP4 fixed vector, AF Number=A, MQ scalar float, SOMATIC flag, STR string)
// used by the INFO-combine parity tests.
func mergeInfoVCF(sample string, records ...string) string {
	hdr := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="dp">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="dp4">
##INFO=<ID=AF,Number=A,Type=Float,Description="af">
##INFO=<ID=MQ,Number=1,Type=Float,Description="mq">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="somatic">
##INFO=<ID=STR,Number=1,Type=String,Description="str">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	` + sample + "\n"
	return hdr + strings.Join(records, "\n") + "\n"
}

// mergeInfoField extracts the INFO column of the first data record from a merged
// VCF body.
func mergeInfoField(t *testing.T, body string) string {
	t.Helper()
	for _, ln := range strings.Split(body, "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		f := strings.Split(ln, "\t")
		if len(f) < 8 {
			t.Fatalf("record has %d columns: %q", len(f), ln)
		}
		return f[7]
	}
	t.Fatalf("no data record in body:\n%s", body)
	return ""
}

// TestMergeInfoByteExact drives a table of crafted merges and asserts the
// merged INFO column matches the byte-exact upstream result. The expected
// strings were validated against bcftools 1.23.1 `merge` (default and
// -i ...) on bgzipped/tabixed copies of the same records.
func TestMergeInfoByteExact(t *testing.T) {
	// Files sharing chr1:100 A>G with differing INFO values.
	a := mergeInfoVCF("SA", "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10;DP4=1,2,3,4;AF=0.25;MQ=60;SOMATIC;STR=foo\tGT\t0/1")
	b := mergeInfoVCF("SB", "chr1\t100\t.\tA\tG\t60\tPASS\tDP=20;DP4=5,6,7,8;AF=0.50;MQ=40;STR=bar\tGT\t1/1")
	// File c has a swapped/extended ALT set (T,G) exercising AGR allele remap.
	c := mergeInfoVCF("SC", "chr1\t100\t.\tA\tT,G\t70\tPASS\tDP=30;DP4=9,10,11,12;AF=0.10,0.90;MQ=30\tGT\t1/2")

	tests := []struct {
		name   string
		inputs []string
		opts   MergeOptions
		want   string // expected INFO column
	}{
		{
			// 3-class ordering: scalars (MQ,SOMATIC,STR first-seen) then rule
			// tags (DP,DP4 summed) then AGR no-rule (AF) last, with AF taking
			// the LAST non-missing value (0.5 from b, not 0.25 from a).
			name:   "default a b (ordering + DP:sum + AF last-wins)",
			inputs: []string{a, b},
			want:   "MQ=60;SOMATIC;STR=foo;DP=30;DP4=6,8,10,12;AF=0.5",
		},
		{
			// File-order swap: scalars re-order to first-seen of b, and AF
			// last-wins now yields 0.25 (from a, the later file).
			name:   "swap b a (AF last-wins order flip)",
			inputs: []string{b, a},
			want:   "MQ=40;STR=bar;SOMATIC;DP=30;DP4=6,8,10,12;AF=0.25",
		},
		{
			// AGR remap: c's ALT (T,G) merges onto union (G,T); AF remaps and
			// last-non-missing-wins per output allele slot.
			name:   "3-way a b c (AGR remap last-wins)",
			inputs: []string{a, b, c},
			want:   "MQ=60;SOMATIC;STR=foo;DP=60;DP4=15,18,21,24;AF=0.9,0.1",
		},
		{
			name:   "3-way c a b (swap, AGR remap)",
			inputs: []string{c, a, b},
			want:   "MQ=30;SOMATIC;STR=foo;DP=60;DP4=15,18,21,24;AF=0.1,0.5",
		},
		{
			// -i replaces defaults: DP,DP4 lose their default sum rule and
			// become plain scalars (first-wins); AF:max combines per-allele.
			name:   "-i DP:max,AF:max",
			inputs: []string{a, b},
			opts:   MergeOptions{InfoRules: "DP:max,AF:max"},
			want:   "DP4=1,2,3,4;MQ=60;SOMATIC;STR=foo;AF=0.5;DP=20",
		},
		{
			name:   "-i DP:avg (integer avg, float render)",
			inputs: []string{a, b},
			opts:   MergeOptions{InfoRules: "DP:avg"},
			want:   "DP4=1,2,3,4;MQ=60;SOMATIC;STR=foo;DP=15;AF=0.5",
		},
		{
			name:   "-i DP:min",
			inputs: []string{a, b},
			opts:   MergeOptions{InfoRules: "DP:min"},
			want:   "DP4=1,2,3,4;MQ=60;SOMATIC;STR=foo;DP=10;AF=0.5",
		},
		{
			name:   "-i STR:join (String tag join)",
			inputs: []string{a, b},
			opts:   MergeOptions{InfoRules: "STR:join"},
			want:   "DP=10;DP4=1,2,3,4;MQ=60;SOMATIC;STR=foo,bar;AF=0.5",
		},
		{
			// Rule tags emit in alphabetical order (DP before MQ) after the
			// non-rule scalars, regardless of the -i spec order.
			name:   "-i MQ:max,DP:sum (rule alpha order)",
			inputs: []string{a, b},
			opts:   MergeOptions{InfoRules: "MQ:max,DP:sum"},
			want:   "DP4=1,2,3,4;SOMATIC;STR=foo;DP=30;MQ=60;AF=0.5",
		},
		{
			// A ruled-AGR tag (AF:sum) is placed among the rule block, not the
			// trailing AGR block, and remaps+sums per allele.
			name:   "-i AF:sum 3-way (ruled AGR)",
			inputs: []string{a, b, c},
			opts:   MergeOptions{InfoRules: "AF:sum"},
			want:   "DP=10;DP4=1,2,3,4;MQ=60;SOMATIC;STR=foo;AF=1.65,0.1",
		},
		{
			name:   "-i - disables all rules (first-wins, AF last)",
			inputs: []string{a, b},
			opts:   MergeOptions{InfoRules: "-"},
			want:   "DP=10;DP4=1,2,3,4;MQ=60;SOMATIC;STR=foo;AF=0.5",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			body, _ := readMerged(t, tc.inputs, tc.opts)
			if got := mergeInfoField(t, body); got != tc.want {
				t.Errorf("INFO mismatch\n got: %s\nwant: %s", got, tc.want)
			}
		})
	}
}

// TestMergeInfoFlagOR verifies a Flag INFO tag is present in the merged record
// when any input carries it (logical OR): a has SOMATIC, b does not.
func TestMergeInfoFlagOR(t *testing.T) {
	a := mergeInfoVCF("SA", "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10;SOMATIC\tGT\t0/1")
	b := mergeInfoVCF("SB", "chr1\t100\t.\tA\tG\t60\tPASS\tDP=20\tGT\t1/1")
	body, _ := readMerged(t, []string{a, b}, MergeOptions{})
	info := mergeInfoField(t, body)
	if !strings.Contains(info, "SOMATIC") {
		t.Errorf("Flag SOMATIC should survive the merge (OR), got INFO=%s", info)
	}
}

// TestMergeInfoMissingTag verifies that when one input lacks a ruled tag, the
// rule combines only the records that carried it (avg/min divide/skip over the
// carrier set, not the full bucket).
func TestMergeInfoMissingTag(t *testing.T) {
	// b has no DP.
	a := mergeInfoVCF("SA", "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10;AF=0.2\tGT\t0/1")
	b := mergeInfoVCF("SB", "chr1\t100\t.\tA\tG\t60\tPASS\tAF=0.4\tGT\t1/1")
	for _, tc := range []struct {
		rule string
		want string // expected DP value
	}{
		{"", "10"},       // default DP:sum over carriers -> 10
		{"DP:avg", "10"}, // avg over 1 carrier -> 10
		{"DP:min", "10"}, // min over 1 carrier -> 10
	} {
		body, _ := readMerged(t, []string{a, b}, MergeOptions{InfoRules: tc.rule})
		info := mergeInfoField(t, body)
		var dp string
		for _, kv := range strings.Split(info, ";") {
			if strings.HasPrefix(kv, "DP=") {
				dp = strings.TrimPrefix(kv, "DP=")
			}
		}
		if dp != tc.want {
			t.Errorf("rule %q: DP=%q, want %q (INFO=%s)", tc.rule, dp, tc.want, info)
		}
	}
}

// TestMergeInfoRuleValidation checks the -i rule validation errors match
// upstream: numeric method on a String tag and an undeclared tag both fail.
func TestMergeInfoRuleValidation(t *testing.T) {
	a := mergeInfoVCF("SA", "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10;STR=foo\tGT\t0/1")
	b := mergeInfoVCF("SB", "chr1\t100\t.\tA\tG\t60\tPASS\tDP=20;STR=bar\tGT\t1/1")

	if _, err := runMergeErr(t, []string{a, b}, MergeOptions{InfoRules: "STR:max"}); err == nil {
		t.Errorf("STR:max should error (numeric op on String)")
	} else if !strings.Contains(err.Error(), "Numeric operation") {
		t.Errorf("unexpected error for STR:max: %v", err)
	}

	if _, err := runMergeErr(t, []string{a, b}, MergeOptions{InfoRules: "NOPE:sum"}); err == nil {
		t.Errorf("NOPE:sum should error (undeclared tag)")
	} else if !strings.Contains(err.Error(), "not defined in the header") {
		t.Errorf("unexpected error for NOPE:sum: %v", err)
	}
}
