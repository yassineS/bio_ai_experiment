package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// csqParityDir holds the csq parity fixtures and the expected
// INFO/BCSQ renderings captured from the upstream bcftools binary
// (`bcftools csq <args> | bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'`).
// The goldens were produced with bcftools 1.23.1-73-ge0ec6ab0; see
// docs/PARITY_ROADMAP.md#bcftools.
const csqParityDir = "../../testdata/parity/csq"

// renderCSQBCSQ runs the haplotype-aware engine over the parity fixture
// with the given options and renders one
//
//	POS<TAB>REF<TAB>ALT<TAB>BCSQ
//
// line per record, matching the upstream
// `bcftools query -f'%POS\t%REF\t%ALT\t%INFO/BCSQ\n'` layout used to
// capture the goldens. Records without a BCSQ tag are skipped (upstream
// query prints nothing for them under this format string... in practice
// every parity record carries one).
func renderCSQBCSQ(t *testing.T, opts CSQOptions) string {
	t.Helper()
	opts.FastaRef = filepath.Join(csqParityDir, "csq.fa")
	opts.GFFAnnot = filepath.Join(csqParityDir, "csq.gff3")
	var buf bytes.Buffer
	if _, err := CSQFile(filepath.Join(csqParityDir, "csq.vcf"), &buf, opts); err != nil {
		t.Fatalf("CSQFile: %v", err)
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
		bcsq, ok := v.Info["BCSQ"]
		if !ok {
			continue
		}
		sb.WriteString(itoa(v.Pos))
		sb.WriteByte('\t')
		sb.WriteString(v.Ref)
		sb.WriteByte('\t')
		sb.WriteString(strings.Join(v.Alt, ","))
		sb.WriteByte('\t')
		sb.WriteString(sortCSQField(bcsq))
		sb.WriteByte('\n')
	}
	return sb.String()
}

// sortGoldenBCSQ applies the same per-tag consequence sort to the
// captured upstream golden so a benign reordering of comma-separated
// consequences within a single BCSQ tag does not fail the comparison
// (mirrors the upstream test/csq sort-csq helper).
func sortGoldenBCSQ(golden string) string {
	var sb strings.Builder
	for _, line := range strings.Split(golden, "\n") {
		if line == "" {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) == 4 {
			cols[3] = sortCSQField(cols[3])
		}
		sb.WriteString(strings.Join(cols, "\t"))
		sb.WriteByte('\n')
	}
	return sb.String()
}

// TestCSQParityBriefGencode validates `bcftools csq -b/--brief-predictions`,
// `-B/--trim-protein-seq` and `-C/--genetic-code` against goldens
// captured from the upstream C binary. The test auto-skips when the
// golden file is absent (the goldens are committed, so this only fires
// if the fixture directory is incomplete).
func TestCSQParityBriefGencode(t *testing.T) {
	cases := []struct {
		name   string
		opts   CSQOptions
		golden string
	}{
		// -b is upstream's alias for -B 1: abbreviate each amino-acid
		// prediction to its first residue plus "..<index>".
		{"brief", CSQOptions{TrimProteinSeq: 1}, "csq.brief.bcsq.expected"},
		// -B 2 keeps two residues before the "..".
		{"trim2", CSQOptions{TrimProteinSeq: 2}, "csq.trim2.bcsq.expected"},
		// -C 1 is the full standard table (identical AAs to table 0 here).
		{"gencode1", CSQOptions{GeneticCode: 1}, "csq.gc1.bcsq.expected"},
		// -C 2 is the vertebrate mitochondrial table (e.g. ATA->M).
		{"gencode2", CSQOptions{GeneticCode: 2}, "csq.gc2.bcsq.expected"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			raw, err := os.ReadFile(filepath.Join(csqParityDir, tc.golden))
			if err != nil {
				if os.IsNotExist(err) {
					t.Skipf("golden %s absent; regenerate from upstream bcftools (see docs/PARITY_ROADMAP.md#bcftools)", tc.golden)
				}
				t.Fatalf("read golden: %v", err)
			}
			want := sortGoldenBCSQ(string(raw))
			got := renderCSQBCSQ(t, tc.opts)
			if got != want {
				t.Errorf("%s: BCSQ output does not match upstream golden %s\n--- got ---\n%s\n--- want ---\n%s",
					tc.name, tc.golden, got, want)
			}
		})
	}
}

// TestGencodeTablesWellFormed asserts each transcribed NCBI table has
// the 64-entry code/stop strings the codon index requires.
func TestGencodeTablesWellFormed(t *testing.T) {
	for _, gc := range gencodeTables {
		if len(gc.Code) != 64 {
			t.Errorf("gencode %d (%s): len(Code)=%d, want 64", gc.ID, gc.Name, len(gc.Code))
		}
		if len(gc.Stop) != 64 {
			t.Errorf("gencode %d (%s): len(Stop)=%d, want 64", gc.ID, gc.Name, len(gc.Stop))
		}
	}
	if _, ok := gencodeByID(0); !ok {
		t.Error("gencodeByID(0): standard table must be present")
	}
	if _, ok := gencodeByID(99); ok {
		t.Error("gencodeByID(99): expected absent")
	}
	if !GeneticCodeKnown(2) || GeneticCodeKnown(99) {
		t.Error("GeneticCodeKnown: 2 should be known and 99 unknown")
	}
}

// TestKprintAAPrediction exercises the brief-prediction truncation
// directly, mirroring upstream's kprint_aa_prediction edge cases.
func TestKprintAAPrediction(t *testing.T) {
	cases := []struct {
		name  string
		brief int
		beg   int
		aa    string
		stop  string
		want  string
	}{
		// brief==0 => full prediction.
		{"full", 0, 12, "YVRT", "----", "YVRT"},
		// short prediction (len-brief < 3) => not abbreviated.
		{"too-short", 1, 12, "YVR", "---", "YVR"},
		// -b (brief==1): first residue + ".." + index past the prediction.
		{"brief1", 1, 12, "YVRT", "----", "Y..16"},
		// -B 2: keep two residues (needs len-brief>=3 to abbreviate).
		{"brief2", 2, 11, "VVRTY", "-----", "VV..16"},
		// trailing stop codon is dropped from the length.
		{"drop-stop", 1, 11, "VVRTY*", "-----*", "V..16"},
	}
	e := &hapEngine{}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			e.brief = tc.brief
			var sb strings.Builder
			e.kprintAAPrediction(&sb, tc.beg, []byte(tc.aa), []byte(tc.stop))
			if sb.String() != tc.want {
				t.Errorf("kprintAAPrediction(%q): got %q, want %q", tc.aa, sb.String(), tc.want)
			}
		})
	}
}
