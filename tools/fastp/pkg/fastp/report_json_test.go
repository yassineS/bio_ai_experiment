package fastp

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

func TestWriteJSONReportTopLevelKeys(t *testing.T) {
	stats := makeRichStats(true)
	stats.DetectedAdapter = "AGATCGGAAGAGC"
	stats.DetectedAdapterR1 = "AGATCGGAAGAGC"
	stats.DetectedAdapterR2 = "CTGTCTCTTATA"

	dir := t.TempDir()
	path := filepath.Join(dir, "report.json")
	if err := WriteJSONReport(path, stats); err != nil {
		t.Fatalf("WriteJSONReport: %v", err)
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	var m map[string]interface{}
	if err := json.Unmarshal(raw, &m); err != nil {
		t.Fatalf("unmarshal: %v\n%s", err, raw)
	}
	wantKeys := []string{
		"summary",
		"filtering_result",
		"duplication",
		"adapter_cutting",
		"read1_before_filtering",
		"read2_before_filtering",
		"read1_after_filtering",
		"read2_after_filtering",
		"tool",
	}
	for _, k := range wantKeys {
		if _, ok := m[k]; !ok {
			t.Errorf("missing top-level key %q", k)
		}
	}

	// Spot-check the summary section's nested keys.
	summary, ok := m["summary"].(map[string]interface{})
	if !ok {
		t.Fatalf("summary is not an object")
	}
	for _, sub := range []string{"before_filtering", "after_filtering"} {
		section, ok := summary[sub].(map[string]interface{})
		if !ok {
			t.Errorf("summary.%s is not an object", sub)
			continue
		}
		for _, field := range []string{
			"total_reads", "total_bases", "q20_bases", "q30_bases",
			"q20_rate", "q30_rate", "gc_content", "read1_mean_length",
		} {
			if _, ok := section[field]; !ok {
				t.Errorf("summary.%s missing %q", sub, field)
			}
		}
	}

	// adapter_cutting should include adapter_trimmed_reads/bases plus the
	// detected adapter sequence.
	ac, ok := m["adapter_cutting"].(map[string]interface{})
	if !ok {
		t.Fatalf("adapter_cutting is not an object")
	}
	if _, ok := ac["adapter_trimmed_reads"]; !ok {
		t.Errorf("adapter_cutting missing adapter_trimmed_reads")
	}
	if _, ok := ac["adapter_trimmed_bases"]; !ok {
		t.Errorf("adapter_cutting missing adapter_trimmed_bases")
	}
	if got, _ := ac["detected_adapter_sequence"].(string); got != "AGATCGGAAGAGC" {
		t.Errorf("detected_adapter_sequence = %q, want AGATCGGAAGAGC", got)
	}
}

func TestWriteJSONReportSingleEndOmitsR2(t *testing.T) {
	stats := makeRichStats(false)
	dir := t.TempDir()
	path := filepath.Join(dir, "se.json")
	if err := WriteJSONReport(path, stats); err != nil {
		t.Fatalf("WriteJSONReport: %v", err)
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	var m map[string]interface{}
	if err := json.Unmarshal(raw, &m); err != nil {
		t.Fatalf("unmarshal: %v", err)
	}
	if _, ok := m["read2_before_filtering"]; ok {
		t.Errorf("SE report should not include read2_before_filtering")
	}
	if _, ok := m["read2_after_filtering"]; ok {
		t.Errorf("SE report should not include read2_after_filtering")
	}
}

func TestWriteJSONReportNilStats(t *testing.T) {
	if err := WriteJSONReport("/tmp/never.json", nil); err == nil {
		t.Errorf("expected error for nil stats")
	}
}

func TestWriteJSONReportEndToEnd(t *testing.T) {
	in := `@r1
ACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIII
@r2
ACGTACGTACGTACGTACGTAAAA
+
IIIIIIIIIIIIIIIIIIIIIIII
`
	opts := DefaultProcessOptions()
	opts.MinLength = 5
	stats, err := ProcessSingleEnd(strings.NewReader(in), &discardWriter{}, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd: %v", err)
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "end2end.json")
	if err := WriteJSONReport(path, stats); err != nil {
		t.Fatalf("WriteJSONReport: %v", err)
	}
	raw, _ := os.ReadFile(path)
	var m map[string]interface{}
	if err := json.Unmarshal(raw, &m); err != nil {
		t.Fatalf("unmarshal: %v", err)
	}
	summary, _ := m["summary"].(map[string]interface{})
	before, _ := summary["before_filtering"].(map[string]interface{})
	if before["total_reads"].(float64) != 2 {
		t.Errorf("expected 2 reads before filtering, got %v", before["total_reads"])
	}
	tool, ok := m["tool"].(map[string]interface{})
	if !ok || tool["name"] != "fastp" {
		t.Errorf("tool.name missing/wrong: %v", m["tool"])
	}
}

// discardWriter implements io.Writer and throws data away. Used as the
// output sink for tests that only care about stats.
type discardWriter struct{}

func (d *discardWriter) Write(p []byte) (int, error) { return len(p), nil }
