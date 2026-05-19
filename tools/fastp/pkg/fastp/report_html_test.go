package fastp

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// makeRichStats builds a ProcessStats populated with per-cycle quality,
// composition, and length data so the SVG renderers have non-empty
// inputs to exercise. Both R1 and R2 are populated when paired is true.
func makeRichStats(paired bool) *ProcessStats {
	s := &ProcessStats{}
	s.TotalReads = 100
	s.TotalBases = 10000
	s.CleanReads = 80
	s.CleanBases = 7800
	s.LowQualityReads = 10
	s.TooShortReads = 5
	s.TooManyNReads = 2
	s.TooLongReads = 0
	s.AdapterTrimmedReads = 30
	s.AdapterTrimmedBases = 600
	s.Q20BasesBefore = 9000
	s.Q30BasesBefore = 7000
	s.GCBasesBefore = 5000
	s.Q20BasesAfter = 7700
	s.Q30BasesAfter = 6000
	s.GCBasesAfter = 4000
	s.TotalReadsR1 = 50
	s.TotalBasesR1 = 5000
	s.CleanReadsR1 = 40
	s.CleanBasesR1 = 3900
	if paired {
		s.TotalReadsR2 = 50
		s.TotalBasesR2 = 5000
		s.CleanReadsR2 = 40
		s.CleanBasesR2 = 3900
	}

	cycles := 100
	for r := 0; r < 2; r++ {
		if r == 1 && !paired {
			continue
		}
		s.QualSumByCycle[r] = make([]int64, cycles)
		s.QualCountByCycle[r] = make([]int64, cycles)
		for b := 0; b < 5; b++ {
			s.BaseCountByCycle[r][b] = make([]int64, cycles)
		}
		for i := 0; i < cycles; i++ {
			s.QualSumByCycle[r][i] = int64(50 * (35 - i/20))
			s.QualCountByCycle[r][i] = 50
			// Composition: balanced 25%.
			for b := 0; b < 4; b++ {
				s.BaseCountByCycle[r][b][i] = 12
			}
			s.BaseCountByCycle[r][4][i] = 2
		}
		s.LengthHistAfter[r] = map[int]int64{100: 40, 99: 5, 95: 3}
		s.LengthHistBefore[r] = map[int]int64{100: 50}
	}
	return s
}

func TestWriteHTMLReportContainsAllSections(t *testing.T) {
	stats := makeRichStats(true)
	dir := t.TempDir()
	path := filepath.Join(dir, "report.html")
	if err := WriteHTMLReport(path, stats); err != nil {
		t.Fatalf("WriteHTMLReport: %v", err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read report: %v", err)
	}
	body := string(data)

	required := []string{
		"<title>fastp report</title>",
		"<h2>Summary</h2>",
		"<h2>Per-base quality</h2>",
		"<h2>Per-base composition</h2>",
		"<h2>Length distribution</h2>",
		"<h2>Filtering reasons</h2>",
		"<h2>Adapter trimming</h2>",
		"<svg",
		"Read 1 per-base mean quality",
		"Read 2 per-base mean quality",
		"fastp (Go) v",
	}
	for _, want := range required {
		if !strings.Contains(body, want) {
			t.Errorf("HTML report missing %q", want)
		}
	}
	// No JS, no CDN.
	if strings.Contains(strings.ToLower(body), "<script") {
		t.Errorf("HTML report contains <script>; should be JS-free")
	}
	if strings.Contains(body, "https://cdn") || strings.Contains(body, "http://cdn") {
		t.Errorf("HTML report references a CDN")
	}
}

func TestWriteHTMLReportSingleEndSkipsR2Sections(t *testing.T) {
	stats := makeRichStats(false)
	dir := t.TempDir()
	path := filepath.Join(dir, "se.html")
	if err := WriteHTMLReport(path, stats); err != nil {
		t.Fatalf("WriteHTMLReport: %v", err)
	}
	data, _ := os.ReadFile(path)
	body := string(data)
	if !strings.Contains(body, "Read 1 per-base mean quality") {
		t.Errorf("missing R1 quality chart")
	}
	if strings.Contains(body, "Read 2 per-base mean quality") {
		t.Errorf("SE report should not contain R2 chart")
	}
}

func TestWriteHTMLReportEmptyStats(t *testing.T) {
	stats := &ProcessStats{}
	dir := t.TempDir()
	path := filepath.Join(dir, "empty.html")
	if err := WriteHTMLReport(path, stats); err != nil {
		t.Fatalf("WriteHTMLReport on empty stats: %v", err)
	}
	data, _ := os.ReadFile(path)
	if !bytes.Contains(data, []byte("<svg")) {
		t.Errorf("empty-stats report should still contain placeholder SVGs")
	}
}

func TestWriteHTMLReportNilStats(t *testing.T) {
	if err := WriteHTMLReport("/tmp/never.html", nil); err == nil {
		t.Errorf("expected error for nil stats")
	}
}

func TestWriteHTMLReportDuplicationSection(t *testing.T) {
	stats := makeRichStats(false)
	stats.DupRate = 0.25
	stats.DupTotal = 100
	stats.DedupDropped = 7
	stats.DupHist = map[int]int64{1: 75, 2: 20, 5: 5}

	dir := t.TempDir()
	path := filepath.Join(dir, "dup.html")
	if err := WriteHTMLReport(path, stats); err != nil {
		t.Fatalf("WriteHTMLReport: %v", err)
	}
	body, _ := os.ReadFile(path)
	for _, want := range []string{
		"<h2>Duplication</h2>",
		"Duplication rate",
		"25.00%",
		"Reads dropped by --dedup",
		"Duplication levels",
	} {
		if !strings.Contains(string(body), want) {
			t.Errorf("duplication section missing %q", want)
		}
	}
}

func TestWriteHTMLReportNoDuplicationSectionWhenDisabled(t *testing.T) {
	stats := makeRichStats(false)
	// stats.DupTotal == 0 -> section should be hidden.
	dir := t.TempDir()
	path := filepath.Join(dir, "nodup.html")
	if err := WriteHTMLReport(path, stats); err != nil {
		t.Fatalf("WriteHTMLReport: %v", err)
	}
	body, _ := os.ReadFile(path)
	if strings.Contains(string(body), "<h2>Duplication</h2>") {
		t.Errorf("Duplication section rendered without DupTotal > 0")
	}
}

func TestWriteHTMLReportRendersFromActualRun(t *testing.T) {
	// Build a tiny SE input, run the processor, then render. This is the
	// end-to-end path the CLI uses; it catches stat-collection regressions.
	in := `@r1
ACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIII
@r2
ACGTACGTACGTACGTACGTAAAA
+
IIIIIIIIIIIIIIIIIIIIIIII
`
	var out bytes.Buffer
	opts := DefaultProcessOptions()
	opts.MinLength = 5
	stats, err := ProcessSingleEnd(strings.NewReader(in), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd: %v", err)
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "end2end.html")
	if err := WriteHTMLReport(path, stats); err != nil {
		t.Fatalf("WriteHTMLReport: %v", err)
	}
	data, _ := os.ReadFile(path)
	if !strings.Contains(string(data), "Read 1 per-base mean quality") {
		t.Errorf("expected quality section in end-to-end HTML report")
	}
}
