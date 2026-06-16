package cram

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// shuffledHeaderText is a SAM header whose @-lines are deliberately NOT in
// htslib's canonical emission order: a comment first, then @PG, @RG, an @SQ,
// the @HD, and another @SQ. A faithful CRAM writer must re-emit these in
// htslib's grouped order (@HD, @CO, @PG, @RG, @SQ) so the embedded header
// byte-matches what samtools writes.
const shuffledHeaderText = "@CO\ta comment first\n" +
	"@PG\tID:prog1\tPN:prog1\n" +
	"@RG\tID:rg1\tSM:s1\n" +
	"@SQ\tSN:chr2\tLN:50000\n" +
	"@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:100000\n"

// TestCRAMHeaderCanonicalOrder asserts the CRAM writer serialises the SAM
// header in htslib's canonical @-line order, so that upstream
// `samtools view -H` on a CRAM WE wrote yields exactly the same header (and
// thus the same order) as `samtools view -H` on a CRAM samtools wrote from the
// identical input SAM.
func TestCRAMHeaderCanonicalOrder(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	dir := t.TempDir()

	// 1. Our CRAM, written from the shuffled header.
	h, err := sam.ParseHeaderText(shuffledHeaderText)
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	ourCRAM := filepath.Join(dir, "ours.cram")
	f, err := os.Create(ourCRAM)
	if err != nil {
		t.Fatalf("create ours.cram: %v", err)
	}
	rw, err := NewRecordWriter(f, h)
	if err != nil {
		f.Close()
		t.Fatalf("NewRecordWriter: %v", err)
	}
	if err := rw.Close(); err != nil {
		f.Close()
		t.Fatalf("RecordWriter.Close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("close ours.cram: %v", err)
	}

	// 2. samtools' own CRAM from the identical SAM input.
	srcSAM := filepath.Join(dir, "src.sam")
	if err := os.WriteFile(srcSAM, []byte(shuffledHeaderText), 0o644); err != nil {
		t.Fatalf("write src.sam: %v", err)
	}
	upstreamCRAM := filepath.Join(dir, "upstream.cram")
	cmd := exec.Command(samtools, "view", "-C", "--no-PG", "-o", upstreamCRAM, srcSAM)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools view -C: %v\n%s", err, out)
	}

	// 3. Decode both headers with the SAME tool (samtools view -H --no-PG) so
	// only the stored @-line order can differ.
	ourHdr := samtoolsHeader(t, samtools, ourCRAM)
	upHdr := samtoolsHeader(t, samtools, upstreamCRAM)

	if ourHdr != upHdr {
		t.Fatalf("CRAM header order differs from upstream.\n--- ours ---\n%s\n--- upstream ---\n%s", ourHdr, upHdr)
	}

	// Pin the expected canonical order explicitly so a future regression is
	// obvious even if upstream's behaviour drifts.
	wantOrder := []string{"@HD", "@CO", "@PG", "@RG", "@SQ", "@SQ"}
	gotOrder := headerLineTags(ourHdr)
	if len(gotOrder) != len(wantOrder) {
		t.Fatalf("header line count = %d (%v), want %d (%v)", len(gotOrder), gotOrder, len(wantOrder), wantOrder)
	}
	for i := range wantOrder {
		if gotOrder[i] != wantOrder[i] {
			t.Fatalf("header line %d = %s, want %s (full order %v)", i, gotOrder[i], wantOrder[i], gotOrder)
		}
	}
}

// samtoolsHeader runs `samtools view -H --no-PG file` and returns the header
// text. --no-PG keeps samtools from injecting its own @PG line on read.
func samtoolsHeader(t *testing.T, samtools, file string) string {
	t.Helper()
	out, err := exec.Command(samtools, "view", "-H", "--no-PG", file).Output()
	if err != nil {
		t.Fatalf("samtools view -H %s: %v", file, err)
	}
	return string(out)
}

// headerLineTags returns the leading @XY tag of each non-empty line.
func headerLineTags(header string) []string {
	var tags []string
	for _, line := range strings.Split(header, "\n") {
		if line == "" {
			continue
		}
		if len(line) >= 3 {
			tags = append(tags, line[:3])
		}
	}
	return tags
}

// TestCRAMEmbeddedHeaderTextCanonical checks the writer's embedded SAM header
// block round-trips through our own reader in canonical order, without needing
// samtools — a fast guard that complements the upstream parity test.
func TestCRAMEmbeddedHeaderTextCanonical(t *testing.T) {
	h, err := sam.ParseHeaderText(shuffledHeaderText)
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, nil); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	got := rr.Header().Text()
	wantOrder := []string{"@HD", "@CO", "@PG", "@RG", "@SQ", "@SQ"}
	gotOrder := headerLineTags(got)
	if len(gotOrder) != len(wantOrder) {
		t.Fatalf("embedded header line count = %d (%v), want %d", len(gotOrder), gotOrder, len(wantOrder))
	}
	for i := range wantOrder {
		if gotOrder[i] != wantOrder[i] {
			t.Fatalf("embedded header line %d = %s, want %s (full %v)", i, gotOrder[i], wantOrder[i], gotOrder)
		}
	}
}
