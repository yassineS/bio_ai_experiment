// JSON report generator. The schema is intentionally close to upstream
// fastp's fastp.json layout so existing tooling (MultiQC etc.) can
// consume our output. Notable top-level keys: summary, filtering_result,
// duplication, adapter_cutting, read1_before_filtering,
// read2_before_filtering, read1_after_filtering, read2_after_filtering.

package fastp

import (
	"encoding/json"
	"fmt"
	"os"
	"time"
)

// reportJSON is the top-level JSON structure written by WriteJSONReport.
type reportJSON struct {
	Summary           jsonSummary        `json:"summary"`
	FilteringResult   jsonFiltering      `json:"filtering_result"`
	Duplication       jsonDuplication    `json:"duplication"`
	AdapterCutting    jsonAdapterCutting `json:"adapter_cutting"`
	Read1BeforeFilter jsonReadStats      `json:"read1_before_filtering"`
	Read2BeforeFilter *jsonReadStats     `json:"read2_before_filtering,omitempty"`
	Read1AfterFilter  jsonReadStats      `json:"read1_after_filtering"`
	Read2AfterFilter  *jsonReadStats     `json:"read2_after_filtering,omitempty"`
	Tool              jsonTool           `json:"tool"`
}

type jsonSummary struct {
	Before jsonSummarySection `json:"before_filtering"`
	After  jsonSummarySection `json:"after_filtering"`
}

type jsonSummarySection struct {
	TotalReads      int64   `json:"total_reads"`
	TotalBases      int64   `json:"total_bases"`
	Q20Bases        int64   `json:"q20_bases"`
	Q30Bases        int64   `json:"q30_bases"`
	Q20Rate         float64 `json:"q20_rate"`
	Q30Rate         float64 `json:"q30_rate"`
	GCContent       float64 `json:"gc_content"`
	Read1MeanLength int     `json:"read1_mean_length"`
	Read2MeanLength int     `json:"read2_mean_length,omitempty"`
}

type jsonFiltering struct {
	PassedFilterReads  int64 `json:"passed_filter_reads"`
	CorrectedReads     int64 `json:"corrected_reads,omitempty"`
	CorrectedBases     int64 `json:"corrected_bases,omitempty"`
	LowQualityReads    int64 `json:"low_quality_reads"`
	TooManyNReads      int64 `json:"too_many_N_reads"`
	LowComplexityReads int64 `json:"low_complexity_reads"`
	TooShortReads      int64 `json:"too_short_reads"`
	TooLongReads       int64 `json:"too_long_reads"`
}

type jsonDuplication struct {
	Rate      float64          `json:"rate"`
	Histogram map[string]int64 `json:"histogram,omitempty"`
}

type jsonAdapterCutting struct {
	AdapterTrimmedReads int64  `json:"adapter_trimmed_reads"`
	AdapterTrimmedBases int64  `json:"adapter_trimmed_bases"`
	Read1Adapter        string `json:"read1_adapter_sequence,omitempty"`
	Read2Adapter        string `json:"read2_adapter_sequence,omitempty"`
	DetectedAdapter     string `json:"detected_adapter_sequence,omitempty"`
}

type jsonReadStats struct {
	TotalReads               int64                `json:"total_reads"`
	TotalBases               int64                `json:"total_bases"`
	Q20Bases                 int64                `json:"q20_bases"`
	Q30Bases                 int64                `json:"q30_bases"`
	TotalCycles              int                  `json:"total_cycles"`
	QualityCurves            map[string][]float64 `json:"quality_curves,omitempty"`
	ContentCurves            map[string][]float64 `json:"content_curves,omitempty"`
	LengthDistribution       map[string]int64     `json:"length_distribution,omitempty"`
	OverrepresentedSequences map[string]int64     `json:"overrepresented_sequences,omitempty"`
}

type jsonTool struct {
	Name    string `json:"name"`
	Version string `json:"version"`
	Time    string `json:"time"`
}

// ToolVersion is the fastp Go implementation version string. Bumped here
// (rather than in a build flag) so the value is captured in the JSON and
// HTML reports.
const ToolVersion = "1.1.0"

// WriteJSONReport writes a fastp-compatible JSON report to path.
//
// The report aggregates the per-read counters in stats into the
// upstream fastp schema. It is safe to call after a streaming run from
// stdin/stdout because all required data is collected in stats during
// processing.
func WriteJSONReport(path string, stats *ProcessStats) error {
	if stats == nil {
		return fmt.Errorf("stats is nil")
	}

	report := buildJSONReport(stats)
	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return fmt.Errorf("marshal JSON report: %w", err)
	}
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create JSON report %q: %w", path, err)
	}
	defer f.Close()
	if _, err := f.Write(data); err != nil {
		return fmt.Errorf("write JSON report: %w", err)
	}
	return nil
}

// buildJSONReport composes the JSON structure for WriteJSONReport. Kept
// public-ish for tests in the same package.
func buildJSONReport(stats *ProcessStats) reportJSON {
	pe := stats.TotalReadsR2 > 0

	totalBefore := stats.TotalBasesR1 + stats.TotalBasesR2
	if totalBefore == 0 {
		totalBefore = stats.TotalBases
	}
	totalAfter := stats.CleanBasesR1 + stats.CleanBasesR2
	if totalAfter == 0 {
		totalAfter = stats.CleanBases
	}

	r1MeanBefore := meanLen(stats.LengthHistBefore[0])
	r2MeanBefore := meanLen(stats.LengthHistBefore[1])
	r1MeanAfter := meanLen(stats.LengthHistAfter[0])
	r2MeanAfter := meanLen(stats.LengthHistAfter[1])

	r := reportJSON{
		Summary: jsonSummary{
			Before: jsonSummarySection{
				TotalReads:      int64(stats.TotalReads),
				TotalBases:      totalBefore,
				Q20Bases:        stats.Q20BasesBefore,
				Q30Bases:        stats.Q30BasesBefore,
				Q20Rate:         safeDiv(stats.Q20BasesBefore, totalBefore),
				Q30Rate:         safeDiv(stats.Q30BasesBefore, totalBefore),
				GCContent:       safeDiv(stats.GCBasesBefore, totalBefore),
				Read1MeanLength: r1MeanBefore,
				Read2MeanLength: r2MeanBefore,
			},
			After: jsonSummarySection{
				TotalReads:      int64(stats.CleanReads),
				TotalBases:      totalAfter,
				Q20Bases:        stats.Q20BasesAfter,
				Q30Bases:        stats.Q30BasesAfter,
				Q20Rate:         safeDiv(stats.Q20BasesAfter, totalAfter),
				Q30Rate:         safeDiv(stats.Q30BasesAfter, totalAfter),
				GCContent:       safeDiv(stats.GCBasesAfter, totalAfter),
				Read1MeanLength: r1MeanAfter,
				Read2MeanLength: r2MeanAfter,
			},
		},
		FilteringResult: jsonFiltering{
			PassedFilterReads:  int64(stats.CleanReads),
			CorrectedReads:     stats.CorrectedReads,
			CorrectedBases:     stats.BasesCorrected,
			LowQualityReads:    int64(stats.LowQualityReads),
			TooManyNReads:      int64(stats.TooManyNReads),
			LowComplexityReads: int64(stats.LowComplexityReads),
			TooShortReads:      int64(stats.TooShortReads),
			TooLongReads:       int64(stats.TooLongReads),
		},
		Duplication: jsonDuplication{
			Rate:      stats.DupRate,
			Histogram: dupHistToJSON(stats.DupHist),
		},
		AdapterCutting: jsonAdapterCutting{
			AdapterTrimmedReads: int64(stats.AdapterTrimmedReads),
			AdapterTrimmedBases: stats.AdapterTrimmedBases,
			Read1Adapter:        stats.DetectedAdapterR1,
			Read2Adapter:        stats.DetectedAdapterR2,
			DetectedAdapter:     stats.DetectedAdapter,
		},
		Read1BeforeFilter: buildReadStats(stats, 0, true),
		Read1AfterFilter:  buildReadStats(stats, 0, false),
		Tool: jsonTool{
			Name:    "fastp",
			Version: ToolVersion,
			Time:    time.Now().UTC().Format(time.RFC3339),
		},
	}
	if pe {
		r2b := buildReadStats(stats, 1, true)
		r2a := buildReadStats(stats, 1, false)
		r.Read2BeforeFilter = &r2b
		r.Read2AfterFilter = &r2a
	}
	return r
}

// buildReadStats fills jsonReadStats for the given read index. before
// selects the BEFORE-filtering histograms; otherwise AFTER.
func buildReadStats(s *ProcessStats, readIdx int, before bool) jsonReadStats {
	hist := s.LengthHistBefore[readIdx]
	if !before {
		hist = s.LengthHistAfter[readIdx]
	}
	out := jsonReadStats{
		TotalReads:         readTotalReads(s, readIdx, before),
		TotalBases:         readTotalBases(s, readIdx, before),
		Q20Bases:           0, // approximate per-read Q20/Q30 not tracked
		Q30Bases:           0,
		TotalCycles:        len(s.QualSumByCycle[readIdx]),
		LengthDistribution: lengthHistToJSON(hist),
	}
	if before {
		// Per-cycle quality + composition only captured BEFORE filtering.
		qc, cc := cycleCurves(s, readIdx)
		out.QualityCurves = qc
		out.ContentCurves = cc
		// Overrepresented sequences are sampled on the before-filtering
		// stream (upstream emits them under read{1,2}_before_filtering).
		if s.overrep[readIdx] != nil {
			out.OverrepresentedSequences = s.overrep[readIdx].passedSequences()
		}
	}
	return out
}

func readTotalReads(s *ProcessStats, readIdx int, before bool) int64 {
	if before {
		if readIdx == 0 {
			return int64(s.TotalReadsR1)
		}
		return int64(s.TotalReadsR2)
	}
	if readIdx == 0 {
		return int64(s.CleanReadsR1)
	}
	return int64(s.CleanReadsR2)
}

func readTotalBases(s *ProcessStats, readIdx int, before bool) int64 {
	if before {
		if readIdx == 0 {
			return s.TotalBasesR1
		}
		return s.TotalBasesR2
	}
	if readIdx == 0 {
		return s.CleanBasesR1
	}
	return s.CleanBasesR2
}

// cycleCurves returns the per-cycle mean quality and per-cycle base
// fraction curves, keyed by metric name.
func cycleCurves(s *ProcessStats, readIdx int) (quality map[string][]float64, content map[string][]float64) {
	n := len(s.QualSumByCycle[readIdx])
	if n == 0 {
		return nil, nil
	}
	mean := make([]float64, n)
	for i := 0; i < n; i++ {
		c := s.QualCountByCycle[readIdx][i]
		if c > 0 {
			mean[i] = float64(s.QualSumByCycle[readIdx][i]) / float64(c)
		}
	}
	quality = map[string][]float64{"mean": mean}

	content = make(map[string][]float64, 5)
	for b, name := range []string{"A", "C", "G", "T", "N"} {
		row := make([]float64, n)
		for i := 0; i < n; i++ {
			c := s.QualCountByCycle[readIdx][i]
			if c > 0 {
				row[i] = float64(s.BaseCountByCycle[readIdx][b][i]) / float64(c)
			}
		}
		content[name] = row
	}
	return quality, content
}

// dupHistToJSON converts a duplication count -> reads-at-that-count map
// into a string-keyed map suitable for JSON. Returns nil if the input is
// empty so the field is omitted from the report.
func dupHistToJSON(hist map[int]int64) map[string]int64 {
	if len(hist) == 0 {
		return nil
	}
	out := make(map[string]int64, len(hist))
	for k, v := range hist {
		out[fmt.Sprintf("%d", k)] = v
	}
	return out
}

// lengthHistToJSON converts a length->count map into a string-keyed
// map suitable for JSON. Keys are stringified ints so the output is a
// regular JSON object (Go forbids non-string JSON object keys).
func lengthHistToJSON(hist map[int]int64) map[string]int64 {
	if len(hist) == 0 {
		return nil
	}
	out := make(map[string]int64, len(hist))
	for k, v := range hist {
		out[fmt.Sprintf("%d", k)] = v
	}
	return out
}

// safeDiv divides a/b and returns 0 when b is zero, avoiding NaN in JSON
// output.
func safeDiv(a, b int64) float64 {
	if b == 0 {
		return 0
	}
	return float64(a) / float64(b)
}

// meanLen returns the integer mean of a length histogram.
func meanLen(hist map[int]int64) int {
	var sum int64
	var n int64
	for l, c := range hist {
		sum += int64(l) * c
		n += c
	}
	if n == 0 {
		return 0
	}
	return int(sum / n)
}
