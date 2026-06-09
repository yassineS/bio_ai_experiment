package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"

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

var (
	upstreamBcftoolsOnce sync.Once
	upstreamBcftoolsPath string
	upstreamBcftoolsErr  error
)

// repoRoot returns the absolute path of the repository root relative to
// this test package (tools/bcftools/pkg/bcftools).
func repoRoot(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", ".."))
	if err != nil {
		t.Fatalf("abs repo root: %v", err)
	}
	return abs
}

// runWithRetry runs cmd, retrying on failure with exponential backoff
// (2/4/8/16s) — used for network-bound submodule fetches. It returns the
// combined output and the final error.
func runWithRetry(dir string, name string, args ...string) ([]byte, error) {
	backoffs := []time.Duration{2 * time.Second, 4 * time.Second, 8 * time.Second, 16 * time.Second}
	var out []byte
	var err error
	for attempt := 0; ; attempt++ {
		cmd := exec.Command(name, args...)
		cmd.Dir = dir
		out, err = cmd.CombinedOutput()
		if err == nil {
			return out, nil
		}
		if attempt >= len(backoffs) {
			return out, err
		}
		time.Sleep(backoffs[attempt])
	}
}

// buildUpstreamBcftools locates or builds the upstream bcftools binary
// and returns its absolute path. The build is memoised across all tests
// in the package via sync.Once.
func buildUpstreamBcftools(t *testing.T) (string, error) {
	t.Helper()
	upstreamBcftoolsOnce.Do(func() {
		root := repoRoot(t)
		htslibDir := filepath.Join(root, "reference_code", "htslib")
		bcftoolsDir := filepath.Join(root, "reference_code", "bcftools")
		binPath := filepath.Join(bcftoolsDir, "bcftools")

		// Fast path: already built.
		if st, err := os.Stat(binPath); err == nil && st.Mode()&0o111 != 0 {
			upstreamBcftoolsPath = binPath
			return
		}

		// Ensure both submodules are checked out. The sentinels are files
		// each submodule genuinely commits at its top level (htslib ships
		// configure.ac, not a generated configure; bcftools ships a static
		// Makefile) — checking a file that the submodule does not commit
		// (e.g. htslib has no Makefile.am) would force an unnecessary
		// network fetch even when the checkout is already complete, which
		// can wrongly fail the build offline. RECURSIVE is required: htslib
		// nests the htscodecs submodule, and a non-recursive init makes
		// htslib's ./configure abort.
		_, htsErr := os.Stat(filepath.Join(htslibDir, "configure.ac"))
		_, bcfErr := os.Stat(filepath.Join(bcftoolsDir, "Makefile"))
		if htsErr != nil || bcfErr != nil {
			if out, err := runWithRetry(root, "git", "submodule", "update", "--init", "--recursive",
				"reference_code/htslib", "reference_code/bcftools"); err != nil {
				upstreamBcftoolsErr = wrapBuildErr("git submodule update", out, err)
				return
			}
		}

		// Build htslib first (autoreconf + configure + make). bcftools
		// links against this in-tree htslib.
		if _, err := os.Stat(filepath.Join(htslibDir, "configure")); err != nil {
			if out, err := runCmd(htslibDir, "autoreconf", "-i"); err != nil {
				upstreamBcftoolsErr = wrapBuildErr("htslib autoreconf", out, err)
				return
			}
		}
		if out, err := runCmd(htslibDir, "./configure"); err != nil {
			upstreamBcftoolsErr = wrapBuildErr("htslib configure", out, err)
			return
		}
		if out, err := runCmd(htslibDir, "make", "-j"); err != nil {
			upstreamBcftoolsErr = wrapBuildErr("htslib make", out, err)
			return
		}

		// Build bcftools. Use a plain `make` — do NOT run bcftools'
		// ./configure, which can clobber htslib's config.mk.
		if out, err := runCmd(bcftoolsDir, "make", "-j"); err != nil {
			upstreamBcftoolsErr = wrapBuildErr("bcftools make", out, err)
			return
		}

		if st, err := os.Stat(binPath); err != nil || st.Mode()&0o111 == 0 {
			upstreamBcftoolsErr = wrapBuildErr("bcftools binary", nil, err)
			return
		}
		upstreamBcftoolsPath = binPath
	})
	return upstreamBcftoolsPath, upstreamBcftoolsErr
}

// runCmd runs name with args in dir and returns combined output and error.
func runCmd(dir string, name string, args ...string) ([]byte, error) {
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	return cmd.CombinedOutput()
}

// wrapBuildErr decorates a build-step failure with the captured output so
// the t.Fatalf message is actionable.
func wrapBuildErr(step string, out []byte, err error) error {
	tail := string(out)
	if len(tail) > 4000 {
		tail = "...[truncated]\n" + tail[len(tail)-4000:]
	}
	if err == nil {
		err = os.ErrInvalid
	}
	return &buildError{step: step, out: tail, err: err}
}

type buildError struct {
	step string
	out  string
	err  error
}

func (e *buildError) Error() string {
	return e.step + ": " + e.err.Error() + "\n--- output ---\n" + e.out
}

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
	bin, err := buildUpstreamBcftools(t)
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
