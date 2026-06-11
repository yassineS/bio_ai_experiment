package bcftools

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// goLocalCSQBCSQByPos runs the Go csq engine with -l/--local-csq over the
// given fixture and returns the sorted INFO/BCSQ tag keyed by POS\tREF\tALT.
func goLocalCSQBCSQByPos(t *testing.T, vcfName, faName, gffName string) map[string]string {
	t.Helper()
	var buf bytes.Buffer
	_, err := CSQFile(
		filepath.Join(csqTestdataDir, vcfName),
		&buf,
		CSQOptions{
			FastaRef: filepath.Join(csqTestdataDir, faName),
			GFFAnnot: filepath.Join(csqTestdataDir, gffName),
			LocalCSQ: true,
			Phase:    'a',
		},
	)
	if err != nil {
		t.Fatalf("CSQFile(-l, %s): %v", vcfName, err)
	}
	return bcsqMapFromVCF(t, &buf)
}

// upstreamLocalCSQBCSQByPos runs the upstream `bcftools csq -l` binary over
// the same fixture and returns the INFO/BCSQ tag keyed the same way.
func upstreamLocalCSQBCSQByPos(t *testing.T, bin, vcfName, faName, gffName string) map[string]string {
	t.Helper()
	cmd := exec.Command(bin, "csq", "-l", "-p", "a",
		"-f", filepath.Join(csqTestdataDir, faName),
		"-g", filepath.Join(csqTestdataDir, gffName),
		filepath.Join(csqTestdataDir, vcfName),
		"-Ov",
	)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools csq -l %s: %v\n%s", vcfName, err, errBuf.String())
	}
	return bcsqMapFromVCF(t, &out)
}

// TestCSQ_LocalUpstreamParity validates `bcftools csq -l/--local-csq` (the
// per-record, non-haplotype-aware caller, upstream test_cds_local) against
// the live upstream C binary. Both implementations run on the same csq
// fixture and their INFO/BCSQ output is compared record-for-record. The
// upstream binary is built once on demand; a build failure is fatal
// (t.Fatalf), never skipped.
func TestCSQ_LocalUpstreamParity(t *testing.T) {
	bin := upstreamBcftools(t)
	cases := []struct {
		name string
		vcf  string
		fa   string
		gff  string
	}{
		{"csq.local", "csq.vcf", "csq.fa", "csq.gff3"},
		{"csq.local.oob-codon", "csq.oob-codon.vcf", "csq.oob-codon.fa", "csq.oob-codon.gff"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := upstreamLocalCSQBCSQByPos(t, bin, tc.vcf, tc.fa, tc.gff)
			// Guard against a vacuous pass: the fixture carries many
			// BCSQ-bearing records, so an empty upstream map would mean the
			// -l invocation silently produced nothing.
			nonEmpty := 0
			for _, v := range want {
				if strings.TrimSpace(v) != "" {
					nonEmpty++
				}
			}
			if nonEmpty == 0 {
				t.Fatalf("%s: upstream produced no INFO/BCSQ; cannot parity-check", tc.name)
			}
			got := goLocalCSQBCSQByPos(t, tc.vcf, tc.fa, tc.gff)
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

// TestCSQ_LocalDiffersFromHaplotype confirms the -l path is genuinely the
// per-record caller and not a silent alias for the haplotype-aware engine:
// the two must disagree on at least one compound consequence in the
// fixture (e.g. POS 107, where the haplotype engine joins 107G>A+108T>A
// into a single missense but the local caller reports only 107's own
// missense). This locks the dispatch wiring so a regression that routed
// -l back through test_cds would be caught.
func TestCSQ_LocalDiffersFromHaplotype(t *testing.T) {
	hap := csqBCSQByPos(t, "csq.vcf", "csq.fa", "csq.gff3")
	local := goLocalCSQBCSQByPos(t, "csq.vcf", "csq.fa", "csq.gff3")
	differ := false
	for key, h := range hap {
		if l, ok := local[key]; ok && l != h {
			differ = true
			break
		}
	}
	if !differ {
		t.Fatalf("local and haplotype-aware BCSQ are identical for every record; -l is not selecting test_cds_local")
	}
}
