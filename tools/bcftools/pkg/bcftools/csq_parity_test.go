package bcftools

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// csqTestdataDir is the vendored upstream csq fixture directory. Only the
// input fixtures (VCF/FASTA/GFF) are committed; the expected BCSQ values
// are produced live by the upstream `bcftools csq` binary at test time.
const csqTestdataDir = "../../testdata/csq"

// sortCSQField replicates the upstream test/csq/sort-csq behaviour:
// comma-separated consequence values within one BCSQ tag are sorted
// lexicographically before comparison.
func sortCSQField(s string) string {
	parts := strings.Split(s, ",")
	sort.Strings(parts)
	return strings.Join(parts, ",")
}

// csqBCSQByPos runs the Go csq engine over the given VCF/FASTA/GFF and
// returns a map keyed by POS\tREF\tALT whose value is the sorted BCSQ INFO
// tag for that record.
func csqBCSQByPos(t *testing.T, vcfName, faName, gffName string) map[string]string {
	t.Helper()
	var buf bytes.Buffer
	_, err := CSQFile(
		filepath.Join(csqTestdataDir, vcfName),
		&buf,
		CSQOptions{
			FastaRef: filepath.Join(csqTestdataDir, faName),
			GFFAnnot: filepath.Join(csqTestdataDir, gffName),
		},
	)
	if err != nil {
		t.Fatalf("CSQFile(%s): %v", vcfName, err)
	}
	return bcsqMapFromVCF(t, &buf)
}

// upstreamCSQBCSQByPos runs the real `bcftools csq` binary over the same
// fixture and returns the BCSQ INFO tag keyed the same way.
func upstreamCSQBCSQByPos(t *testing.T, bin, vcfName, faName, gffName string) map[string]string {
	t.Helper()
	cmd := exec.Command(bin, "csq",
		"-f", filepath.Join(csqTestdataDir, faName),
		"-g", filepath.Join(csqTestdataDir, gffName),
		filepath.Join(csqTestdataDir, vcfName),
		"-Ov",
	)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools csq %s: %v\n%s", vcfName, err, errBuf.String())
	}
	return bcsqMapFromVCF(t, &out)
}

// bcsqMapFromVCF parses a VCF stream and returns the sorted BCSQ INFO tag
// keyed by POS\tREF\tALT.
func bcsqMapFromVCF(t *testing.T, r *bytes.Buffer) map[string]string {
	t.Helper()
	vr := vcf.NewReader(r)
	if _, err := vr.ReadHeader(); err != nil {
		t.Fatalf("read header: %v", err)
	}
	m := make(map[string]string)
	for {
		v, err := vr.Read()
		if err != nil {
			break
		}
		key := itoa(v.Pos) + "\t" + v.Ref + "\t" + strings.Join(v.Alt, ",")
		m[key] = sortCSQField(v.Info["BCSQ"])
	}
	return m
}

// TestCSQ_UpstreamParity runs both the Go csq engine and the live upstream
// `bcftools csq` binary on the same fixtures and asserts the INFO/BCSQ
// haplotype-engine output is identical record-for-record. This replaces the
// former golden-file comparison: the upstream binary is built on demand and
// a build failure is fatal (never skipped).
func TestCSQ_UpstreamParity(t *testing.T) {
	bin := upstreamBcftools(t)
	cases := []struct {
		name string
		vcf  string
		fa   string
		gff  string
	}{
		{"csq.1", "csq.vcf", "csq.fa", "csq.gff3"},
		{"csq.oob-codon", "csq.oob-codon.vcf", "csq.oob-codon.fa", "csq.oob-codon.gff"},
		{"csq.splice.issue-2543", "csq.splice.issue-2543.vcf", "csq.splice.issue-2543.fa", "csq.splice.issue-2543.gff"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := csqBCSQByPos(t, tc.vcf, tc.fa, tc.gff)
			want := upstreamCSQBCSQByPos(t, bin, tc.vcf, tc.fa, tc.gff)
			if len(got) != len(want) {
				t.Errorf("record count: got %d, want %d", len(got), len(want))
			}
			for key, wantBCSQ := range want {
				gotBCSQ, ok := got[key]
				if !ok {
					t.Errorf("missing record %q (upstream BCSQ=%q)", key, wantBCSQ)
					continue
				}
				if gotBCSQ != wantBCSQ {
					t.Errorf("BCSQ mismatch at %q:\n got:  %s\n want: %s", key, gotBCSQ, wantBCSQ)
				}
			}
			for key := range got {
				if _, ok := want[key]; !ok {
					t.Errorf("extra record %q not emitted by upstream", key)
				}
			}
		})
	}
}
