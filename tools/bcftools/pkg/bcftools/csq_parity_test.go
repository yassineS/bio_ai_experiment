package bcftools

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// csqParityDir holds the csq parity fixtures: the input VCF, the GFF3
// annotation and the reference FASTA. Unlike a golden-file setup, the
// expected INFO/BCSQ renderings are NOT committed — the live
// upstream-parity test (TestCSQ_BriefGencodeUpstreamParity) computes
// them in-process by running the upstream bcftools binary.
const csqParityDir = "../../testdata/parity/csq"

// renderCSQBCSQ runs the haplotype-aware Go engine over the parity
// fixture with the given options and renders one
//
//	POS<TAB>REF<TAB>ALT<TAB>BCSQ
//
// line per record that carries an INFO/BCSQ tag, matching the upstream
// `bcftools query -f'%POS\t%REF\t%ALT\t%INFO/BCSQ\n'` layout. Records
// without a BCSQ tag are skipped (upstream query prints nothing for them
// under this format string).
func renderCSQBCSQ(t *testing.T, opts CSQOptions) string {
	t.Helper()
	opts.FastaRef = filepath.Join(csqParityDir, "csq.fa")
	opts.GFFAnnot = filepath.Join(csqParityDir, "csq.gff3")
	var buf bytes.Buffer
	if _, err := CSQFile(filepath.Join(csqParityDir, "csq.vcf"), &buf, opts); err != nil {
		t.Fatalf("CSQFile: %v", err)
	}
	return renderBCSQFromVCF(t, buf.Bytes())
}

// renderBCSQFromVCF parses a VCF byte stream (from either the Go engine
// or the upstream binary) and renders the per-record
// POS<TAB>REF<TAB>ALT<TAB>BCSQ summary, applying sortCSQField so a benign
// reordering of comma-separated consequences within a single BCSQ tag
// does not cause a spurious mismatch (mirrors upstream test/csq sort-csq).
func renderBCSQFromVCF(t *testing.T, vcfBytes []byte) string {
	t.Helper()
	vr := vcf.NewReader(bytes.NewReader(vcfBytes))
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

// =====================================================================
// Unit tests (always run): translation tables, brief truncation,
// -B/--trim-protein-seq validation.
// =====================================================================

// TestGencodeTablesWellFormed asserts each transcribed NCBI table
// (ids 0,1,2,3,5) has the 64-entry code/stop strings the codon index
// requires, and that lookup helpers agree on which ids are known.
func TestGencodeTablesWellFormed(t *testing.T) {
	wantIDs := map[int]bool{0: true, 1: true, 2: true, 3: true, 5: true}
	gotIDs := map[int]bool{}
	for _, gc := range gencodeTables {
		gotIDs[gc.ID] = true
		if len(gc.Code) != 64 {
			t.Errorf("gencode %d (%s): len(Code)=%d, want 64", gc.ID, gc.Name, len(gc.Code))
		}
		if len(gc.Stop) != 64 {
			t.Errorf("gencode %d (%s): len(Stop)=%d, want 64", gc.ID, gc.Name, len(gc.Stop))
		}
	}
	for id := range wantIDs {
		if !gotIDs[id] {
			t.Errorf("gencode table %d missing", id)
		}
		if _, ok := gencodeByID(id); !ok {
			t.Errorf("gencodeByID(%d): expected present", id)
		}
		if !GeneticCodeKnown(id) {
			t.Errorf("GeneticCodeKnown(%d): expected known", id)
		}
	}
	if _, ok := gencodeByID(99); ok {
		t.Error("gencodeByID(99): expected absent")
	}
	if GeneticCodeKnown(99) {
		t.Error("GeneticCodeKnown(99): expected unknown")
	}
}

// TestGencodeTranslation exercises codon->amino-acid translation for the
// supported tables, locking in the residues that distinguish each table
// from the standard one (e.g. table 2 vertebrate-mitochondrial translates
// ATA->M and AGA->* where the standard table gives I and R).
func TestGencodeTranslation(t *testing.T) {
	type aa struct {
		c0, c1, c2 byte
		want       byte
	}
	cases := []struct {
		id     int
		codons []aa
	}{
		// Standard simplified (0) and Standard (1) agree on these.
		{0, []aa{{'A', 'T', 'G', 'M'}, {'A', 'T', 'A', 'I'}, {'A', 'G', 'A', 'R'}, {'T', 'G', 'G', 'W'}, {'T', 'A', 'A', '*'}}},
		{1, []aa{{'A', 'T', 'G', 'M'}, {'A', 'T', 'A', 'I'}, {'A', 'G', 'A', 'R'}, {'T', 'G', 'G', 'W'}, {'T', 'A', 'A', '*'}}},
		// Vertebrate mitochondrial (2): ATA->M, AGA->*, TGA->W.
		{2, []aa{{'A', 'T', 'A', 'M'}, {'A', 'G', 'A', '*'}, {'T', 'G', 'A', 'W'}}},
		// Yeast mitochondrial (3): CTN->T (e.g. CTG->T), TGA->W.
		{3, []aa{{'C', 'T', 'G', 'T'}, {'A', 'T', 'A', 'M'}, {'T', 'G', 'A', 'W'}}},
		// Invertebrate mitochondrial (5): AGA->S, ATA->M, TGA->W.
		{5, []aa{{'A', 'G', 'A', 'S'}, {'A', 'T', 'A', 'M'}, {'T', 'G', 'A', 'W'}}},
	}
	for _, tc := range cases {
		gc, ok := gencodeByID(tc.id)
		if !ok {
			t.Fatalf("gencodeByID(%d): not found", tc.id)
		}
		for _, c := range tc.codons {
			if got := gc.dna2aa(c.c0, c.c1, c.c2); got != c.want {
				t.Errorf("table %d: dna2aa(%c%c%c)=%c, want %c", tc.id, c.c0, c.c1, c.c2, got, c.want)
			}
		}
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

// TestKprintAAPredictionGencodeAware confirms the brief truncation is
// applied independently of the genetic-code table — the truncation
// operates on the already-translated residue string, so any table that
// produces a >=3 residue prediction abbreviates the same way.
func TestKprintAAPredictionGencodeAware(t *testing.T) {
	e := &hapEngine{brief: 1}
	var sb strings.Builder
	// A vertebrate-mitochondrial run starting at residue 8: first residue
	// kept, index points one past the prediction (8+4 == 12).
	e.kprintAAPrediction(&sb, 8, []byte("MWSV"), []byte("M---"))
	if got := sb.String(); got != "M..12" {
		t.Errorf("brief mito prediction: got %q, want %q", got, "M..12")
	}
}

// =====================================================================
// Live upstream-binary parity (always runs; builds upstream on demand).
// =====================================================================

// runUpstreamCSQ invokes the upstream bcftools binary
// (`bcftools csq -p a -f ... -g ... <args> csq.vcf`) and returns the
// emitted VCF bytes. -p a (take GTs as phased) matches the Go engine's
// per-haplotype interpretation of the fixture.
func runUpstreamCSQ(t *testing.T, bin string, extraArgs ...string) []byte {
	t.Helper()
	args := []string{"csq", "-p", "a",
		"-f", filepath.Join(csqParityDir, "csq.fa"),
		"-g", filepath.Join(csqParityDir, "csq.gff3"),
	}
	args = append(args, extraArgs...)
	args = append(args, filepath.Join(csqParityDir, "csq.vcf"))
	cmd := exec.Command(bin, args...)
	// Capture stdout (the VCF) separately from stderr, which carries
	// upstream's progress diagnostics ("Parsing ... gff3", deprecation
	// warnings) that would otherwise corrupt the VCF stream.
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools %v: %v\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// TestCSQ_BriefGencodeUpstreamParity validates `bcftools csq`
// -b/--brief-predictions, -B/--trim-protein-seq and -C/--genetic-code
// against the upstream C binary, live and in-process: both upstream and
// the Go port run on the same fixture, and their INFO/BCSQ output is
// compared directly with no committed snapshot. The upstream binary is
// built once on demand (htslib + bcftools); a genuine build failure is a
// hard error (t.Fatalf), never a skip.
func TestCSQ_BriefGencodeUpstreamParity(t *testing.T) {
	bin, err := ensureUpstreamBcftools(t)
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}

	cases := []struct {
		name string
		args []string   // upstream CLI flags
		opts CSQOptions // equivalent Go-port options
	}{
		// -b is upstream's alias for -B 1.
		{"brief", []string{"-b"}, CSQOptions{TrimProteinSeq: 1}},
		// -B 2 keeps two residues before the "..".
		{"trim2", []string{"-B", "2"}, CSQOptions{TrimProteinSeq: 2}},
		// -C 1: full standard table.
		{"gencode1", []string{"-C", "1"}, CSQOptions{GeneticCode: 1}},
		// -C 2: vertebrate mitochondrial table.
		{"gencode2", []string{"-C", "2"}, CSQOptions{GeneticCode: 2}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			wantVCF := runUpstreamCSQ(t, bin, tc.args...)
			want := renderBCSQFromVCF(t, wantVCF)
			// Guard against a vacuous pass: the fixture is known to carry
			// many BCSQ-bearing records, so empty upstream output means the
			// invocation silently produced nothing (e.g. a flag was rejected
			// or the binary changed its tag), which must not match an
			// equally-empty Go run.
			if strings.TrimSpace(want) == "" {
				t.Fatalf("%s: upstream produced no INFO/BCSQ output; cannot parity-check", tc.name)
			}
			got := renderCSQBCSQ(t, tc.opts)
			if got != want {
				t.Errorf("%s: INFO/BCSQ does not match upstream bcftools %v\n--- got (go) ---\n%s\n--- want (upstream) ---\n%s",
					tc.name, tc.args, got, want)
			}
		})
	}
}
