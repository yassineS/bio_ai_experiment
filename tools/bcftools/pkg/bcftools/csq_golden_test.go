package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// csqTestdataDir is the vendored upstream csq fixture directory.
const csqTestdataDir = "../../testdata/csq"

// sortCSQField replicates the upstream test/csq/sort-csq behaviour:
// comma-separated consequence values within one BCSQ/EXP tag are sorted
// lexicographically before comparison.
func sortCSQField(s string) string {
	parts := strings.Split(s, ",")
	sort.Strings(parts)
	return strings.Join(parts, ",")
}

// runCSQToGolden runs the csq engine over the given VCF/FASTA/GFF and
// renders the result as the upstream golden format:
//
//	POS<TAB>REF<TAB>ALT<TAB>EXP
//	POS<TAB>REF<TAB>ALT<TAB>BCSQ
//	<blank line>
//
// matching `bcftools csq ... | sort-csq | bcftools query
// -f'%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n'`.
func runCSQToGolden(t *testing.T, vcfName, faName, gffName string) string {
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

	vr := vcf.NewReader(&buf)
	if _, err := vr.ReadHeader(); err != nil {
		t.Fatalf("read header: %v", err)
	}
	var sb strings.Builder
	for {
		v, err := vr.Read()
		if err != nil {
			break
		}
		alt := strings.Join(v.Alt, ",")
		exp := v.Info["EXP"]
		bcsq := v.Info["BCSQ"]
		sb.WriteString(v.Chrom[:0]) // keep import tidy; POS uses fields below
		sb.WriteString(itoa(v.Pos))
		sb.WriteByte('\t')
		sb.WriteString(v.Ref)
		sb.WriteByte('\t')
		sb.WriteString(alt)
		sb.WriteByte('\t')
		sb.WriteString(sortCSQField(exp))
		sb.WriteByte('\n')
		sb.WriteString(itoa(v.Pos))
		sb.WriteByte('\t')
		sb.WriteString(v.Ref)
		sb.WriteByte('\t')
		sb.WriteString(alt)
		sb.WriteByte('\t')
		sb.WriteString(sortCSQField(bcsq))
		sb.WriteString("\n\n")
	}
	return sb.String()
}

// TestCSQGoldenINFO checks the INFO/BCSQ haplotype-engine output
// byte-for-byte against the vendored upstream goldens.
func TestCSQGoldenINFO(t *testing.T) {
	cases := []struct {
		name   string
		vcf    string
		fa     string
		gff    string
		golden string
	}{
		{"csq.1", "csq.vcf", "csq.fa", "csq.gff3", "csq.1.out"},
		{"csq.oob-codon", "csq.oob-codon.vcf", "csq.oob-codon.fa", "csq.oob-codon.gff", "csq.oob-codon.out"},
		{"csq.splice.issue-2543", "csq.splice.issue-2543.vcf", "csq.splice.issue-2543.fa", "csq.splice.issue-2543.gff", "csq.splice.issue-2543.1.out"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want, err := os.ReadFile(filepath.Join(csqTestdataDir, tc.golden))
			if err != nil {
				t.Fatalf("read golden: %v", err)
			}
			got := runCSQToGolden(t, tc.vcf, tc.fa, tc.gff)
			if got != string(want) {
				t.Errorf("%s: BCSQ output does not match golden %s\n--- got ---\n%s\n--- want ---\n%s",
					tc.name, tc.golden, got, want)
			}
		})
	}
}
