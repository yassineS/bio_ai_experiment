package samtools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// tviewParityRef is a 120 bp reference used by the tview live-parity tests.
// It is plain repetitive enough that the per-position consensus is
// deterministic yet has enough variety to drive non-trivial dot/mismatch
// rendering.
const tviewParityRef = ">chr1\n" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"

// tviewParitySAM exercises every rendering path the text/HTML backends must
// match: forward/reverse matches and mismatches, a deletion (D), a reference
// skip (N), an insertion (I), an improper-pair read (underline attribute),
// and reads packed onto shared rows after a gap.
const tviewParitySAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:120
@RG	ID:rg1	SM:sampleA
@RG	ID:rg2	SM:sampleB
m1	0	chr1	1	60	12M	*	0	0	ACGTACGTACGT	IIIIIIIIIIII	RG:Z:rg1
m2	16	chr1	3	40	10M	*	0	0	GTACGTTCGT	IIIIIIIIII	RG:Z:rg1
del1	0	chr1	5	60	4M2D6M	*	0	0	ACGTGTACGT	IIIIIIIIII	RG:Z:rg2
skip1	0	chr1	8	50	4M3N4M	*	0	0	GTACACGT	IIIIIIII	RG:Z:rg2
ins1	0	chr1	14	60	5M2I5M	*	0	0	ACGTAGGTACGT	IIIIIIIIIIII	RG:Z:rg1
imp1	99	chr1	18	30	8M	*	0	0	GTACGTAC	IIIIIIII	RG:Z:rg2
gapfill	0	chr1	30	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	RG:Z:rg1
rev1	16	chr1	40	45	12M	*	0	0	ACGTACGTACGT	IIIIIIIIIIII	RG:Z:rg2
`

// buildTviewFixture writes the SAM/reference fixtures and uses the upstream
// samtools binary to produce a coordinate-sorted BAM, its .bai index, and the
// reference .fai. It returns (bamPath, refPath). It never t.Skip: an inability
// to build the inputs is a hard failure.
func buildTviewFixture(t *testing.T, bin, samText, refText string) (string, string) {
	t.Helper()
	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	bamPath := filepath.Join(dir, "in.bam")
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(samPath, []byte(samText), 0o600); err != nil {
		t.Fatalf("write SAM: %v", err)
	}
	if err := os.WriteFile(refPath, []byte(refText), 0o600); err != nil {
		t.Fatalf("write reference: %v", err)
	}
	// SAM -> BAM (the fixture is already coordinate-sorted).
	if out, err := exec.Command(bin, "view", "-b", "-o", bamPath, samPath).CombinedOutput(); err != nil {
		t.Fatalf("upstream samtools view -b failed: %v\n%s", err, out)
	}
	if out, err := exec.Command(bin, "index", bamPath).CombinedOutput(); err != nil {
		t.Fatalf("upstream samtools index failed: %v\n%s", err, out)
	}
	if out, err := exec.Command(bin, "faidx", refPath).CombinedOutput(); err != nil {
		t.Fatalf("upstream samtools faidx failed: %v\n%s", err, out)
	}
	return bamPath, refPath
}

// runUpstreamTview invokes upstream `samtools tview` and returns its stdout.
func runUpstreamTview(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, append([]string{"tview"}, args...)...)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream samtools tview %v failed: %v", args, err)
	}
	return out
}

// TestTviewLiveParity is the gold-standard live check: for a fixture that
// exercises matches/mismatches/deletions/refskips/insertions/strand and
// improper-pair underlining, it compares the Go TviewText / TviewHTML output
// byte-for-byte against the vendored upstream `samtools tview -d T|H` over a
// sweep of start positions and widths. Per the project's testing rules it
// builds (never skips) the upstream binary.
func TestTviewLiveParity(t *testing.T) {
	bin := upstreamSamtoolsBinary(t)
	bamPath, refPath := buildTviewFixture(t, bin, tviewParitySAM, tviewParityRef)

	positions := []string{"chr1:1", "chr1:3", "chr1:5", "chr1:14", "chr1:18", "chr1:30", "chr1:40"}
	widths := []int{20, 40, 80}
	modes := []struct {
		flag string
		mode TviewMode
	}{
		{"T", TviewText},
		{"H", TviewHTML},
	}
	inserts := []struct {
		hide bool
		flag []string
	}{
		{false, nil},           // default: expand insertion sub-columns
		{true, []string{"-i"}}, // -i: hide insertions
	}

	for _, pos := range positions {
		for _, w := range widths {
			for _, m := range modes {
				for _, ins := range inserts {
					name := pos + "/w" + itoa(w) + "/-d" + m.flag
					if ins.hide {
						name += "/-i"
					}
					t.Run(name, func(t *testing.T) {
						upArgs := []string{"-d", m.flag, "-w", itoa(w), "-p", pos}
						upArgs = append(upArgs, ins.flag...)
						upArgs = append(upArgs, bamPath, refPath)
						want := runUpstreamTview(t, bin, upArgs...)

						var got bytes.Buffer
						opts := TviewOptions{
							Input:       bamPath,
							Reference:   refPath,
							Position:    pos,
							Width:       w,
							Mode:        m.mode,
							HideInserts: ins.hide,
						}
						if err := Tview(opts, &got); err != nil {
							t.Fatalf("Tview: %v", err)
						}
						if !bytes.Equal(got.Bytes(), want) {
							t.Fatalf("tview -d %s parity mismatch at %s w%d (hide=%v)\n--- got ---\n%s\n--- want ---\n%s",
								m.flag, pos, w, ins.hide, got.Bytes(), want)
						}
					})
				}
			}
		}
	}
}

// TestTviewSampleFilterParity checks that -s (sample / read-group filtering)
// matches upstream byte-for-byte for both an @RG SM: value and an @RG ID:.
func TestTviewSampleFilterParity(t *testing.T) {
	bin := upstreamSamtoolsBinary(t)
	bamPath, refPath := buildTviewFixture(t, bin, tviewParitySAM, tviewParityRef)

	for _, sample := range []string{"sampleA", "sampleB", "rg1", "rg2"} {
		t.Run(sample, func(t *testing.T) {
			want := runUpstreamTview(t, bin, "-d", "T", "-w", "40", "-p", "chr1:1", "-s", sample, bamPath, refPath)
			var got bytes.Buffer
			opts := TviewOptions{
				Input:     bamPath,
				Reference: refPath,
				Position:  "chr1:1",
				Width:     40,
				Sample:    sample,
				Mode:      TviewText,
			}
			if err := Tview(opts, &got); err != nil {
				t.Fatalf("Tview: %v", err)
			}
			if !bytes.Equal(got.Bytes(), want) {
				t.Fatalf("-s %s parity mismatch\n--- got ---\n%s\n--- want ---\n%s", sample, got.Bytes(), want)
			}
		})
	}
}

// TestTviewNoReference checks that, without a reference FASTA, the reference
// line is all 'N' and matches upstream's reference-free rendering.
func TestTviewNoReferenceParity(t *testing.T) {
	bin := upstreamSamtoolsBinary(t)
	bamPath, _ := buildTviewFixture(t, bin, tviewParitySAM, tviewParityRef)

	want := runUpstreamTview(t, bin, "-d", "T", "-w", "40", "-p", "chr1:1", bamPath)
	var got bytes.Buffer
	opts := TviewOptions{
		Input:    bamPath,
		Position: "chr1:1",
		Width:    40,
		Mode:     TviewText,
	}
	if err := Tview(opts, &got); err != nil {
		t.Fatalf("Tview: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("no-reference parity mismatch\n--- got ---\n%s\n--- want ---\n%s", got.Bytes(), want)
	}
}

// TestTviewCharMapping is a focused unit test of the per-read base-to-character
// mapping: forward/reverse match -> '.'/',', mismatch -> UPPER/lower base,
// deletion -> '*', reference skip -> '>'/'<'. It drives drawTviewScreen
// directly with hand-built records and asserts the read rows.
func TestTviewCharMapping(t *testing.T) {
	// ref: AAAAAAAAAA (positions 1..10), so an 'A' read base is a match.
	ref := []byte("AAAAAAAAAA")
	recs := []*sam.Record{
		mkRec(t, "fwd", 0, 1, "AACAA", "5M"),  // match,match,mismatch C,match,match (forward)
		mkRec(t, "rev", 16, 1, "AACAA", "5M"), // same, reverse strand
		mkRec(t, "del", 0, 1, "AAAA", "2M2D2M"),
		mkRec(t, "skip", 0, 1, "AAAA", "2M3N2M"),
	}
	em := errmod.Init(tvErrModDepCorr)
	s := drawTviewScreen(recs, ref, 0, 10, em, true)

	// Row layout: 0 ruler, 1 reference, 2 consensus, 3.. reads. The four
	// reads all start at pos 1 (overlapping), so each lands on its own level
	// (rows 3,4,5,6) in input order.
	rows := renderRows(s)
	if len(rows) < 7 {
		t.Fatalf("expected at least 7 rows, got %d:\n%v", len(rows), rows)
	}
	// fwd read: match '.', match '.', mismatch C -> upper 'C', match, match.
	checkRow(t, rows[3], "..C..", "forward match/mismatch")
	// rev read: matches ',', mismatch 'c' lower.
	checkRow(t, rows[4], ",,c,,", "reverse match/mismatch")
	// del read: 2 match, 2 deletion '*', 2 match.
	checkRow(t, rows[5], "..**..", "deletion")
	// skip read: 2 match, 3 refskip '>', 2 match.
	checkRow(t, rows[6], "..>>>..", "reference skip")
}

// TestTviewGreedyRowPacking checks the level-packing: two non-overlapping
// reads separated by a 2-column gap share a row, while overlapping reads do
// not. We verify by counting the number of distinct read rows produced.
func TestTviewGreedyRowPacking(t *testing.T) {
	ref := bytes.Repeat([]byte("A"), 30)
	// r1 covers pos 1-5; r2 covers pos 8-12 (gap of 2 cols at pos 6,7) so it
	// reuses r1's row; r3 overlaps r1 (pos 1-10) so it needs its own row.
	recs := []*sam.Record{
		mkRec(t, "r1", 0, 1, "AAAAA", "5M"),
		mkRec(t, "r3", 0, 1, "AAAAAAAAAA", "10M"),
		mkRec(t, "r2", 0, 8, "AAAAA", "5M"),
	}
	em := errmod.Init(tvErrModDepCorr)
	s := drawTviewScreen(recs, ref, 0, 30, em, true)
	rows := renderRows(s)
	// Rows 0-2 are ruler/ref/consensus. Read rows are 3+. r1 and r2 share a
	// row (they don't overlap and the gap is >= TV_GAP); r3 gets its own.
	// So there should be exactly 2 read rows (5 total rows).
	readRows := 0
	for i := 3; i < len(rows); i++ {
		if hasNonSpace(rows[i]) {
			readRows++
		}
	}
	if readRows != 2 {
		t.Fatalf("expected 2 read rows (r1+r2 packed, r3 separate), got %d:\n%s", readRows, rows)
	}
	// Confirm r1 and r2 are on the same row: that row has read content at
	// cols 1-5 and 8-12 with a gap at 6-7.
	var packed string
	for i := 3; i < len(rows); i++ {
		if len(rows[i]) >= 12 && rows[i][0] != ' ' && rows[i][7] != ' ' && rows[i][5] == ' ' {
			packed = rows[i]
			break
		}
	}
	if packed == "" {
		t.Fatalf("did not find the packed row (r1 cols1-5, gap cols6-7, r2 cols8-12):\n%s", rows)
	}
}

// --- small test helpers ---
//
// itoa (int -> decimal string) is shared from depth_test.go.

// mkRec builds a mapped sam.Record on chr1 with the given QNAME, flag,
// 1-based position, query sequence, and CIGAR. Per-base qualities default to
// Phred 40 ('I'), which keeps the consensus caller above its min-baseQ floor.
func mkRec(t *testing.T, qname string, flag uint16, pos int, seq, cigar string) *sam.Record {
	t.Helper()
	cig, err := sam.ParseCigar(cigar)
	if err != nil {
		t.Fatalf("ParseCigar(%q): %v", cigar, err)
	}
	qual := make([]byte, len(seq))
	for i := range qual {
		qual[i] = 40
	}
	return &sam.Record{
		QName: qname,
		Flag:  flag,
		RName: "chr1",
		Pos:   int32(pos),
		MapQ:  60,
		Cigar: cig,
		Seq:   seq,
		Qual:  qual,
	}
}

// renderRows returns the screen's rows as strings (mcol chars each, trailing
// spaces preserved) for assertion.
func renderRows(s *tvScreen) []string {
	out := make([]string, len(s.rows))
	for i, row := range s.rows {
		b := make([]byte, s.mcol)
		for x := 0; x < s.mcol; x++ {
			b[x] = row[x].ch
		}
		out[i] = string(b)
	}
	return out
}

// checkRow asserts that row begins with want (ignoring trailing spaces).
func checkRow(t *testing.T, row, want, what string) {
	t.Helper()
	if len(row) < len(want) || row[:len(want)] != want {
		t.Errorf("%s: row prefix = %q, want %q", what, trimTrailingSpaces(row), want)
	}
}

func trimTrailingSpaces(s string) string {
	i := len(s)
	for i > 0 && s[i-1] == ' ' {
		i--
	}
	return s[:i]
}

func hasNonSpace(s string) bool {
	for i := 0; i < len(s); i++ {
		if s[i] != ' ' {
			return true
		}
	}
	return false
}
