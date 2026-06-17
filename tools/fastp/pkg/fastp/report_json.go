// JSON report generator. The schema is intentionally close to upstream
// fastp's fastp.json layout so existing tooling (MultiQC etc.) can
// consume our output. Notable top-level keys: summary, filtering_result,
// duplication, adapter_cutting, read1_before_filtering,
// read2_before_filtering, read1_after_filtering, read2_after_filtering.

package fastp

import (
	"bytes"
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
	InsertSize        *jsonInsertSize    `json:"insert_size,omitempty"`
	AdapterCutting    jsonAdapterCutting `json:"adapter_cutting"`
	Read1BeforeFilter jsonReadStats      `json:"read1_before_filtering"`
	Read2BeforeFilter *jsonReadStats     `json:"read2_before_filtering,omitempty"`
	Read1AfterFilter  jsonReadStats      `json:"read1_after_filtering"`
	Read2AfterFilter  *jsonReadStats     `json:"read2_after_filtering,omitempty"`
	Tool              jsonTool           `json:"tool"`
}

type jsonSummary struct {
	// Sequencing describes the layout and cycle counts, e.g.
	// "single end (100 cycles)" or "paired end (100 cycles + 100 cycles)".
	// It is deterministic given the input, matching upstream's
	// summary.sequencing (jsonreporter.cpp:28-33).
	Sequencing string             `json:"sequencing"`
	Before     jsonSummarySection `json:"before_filtering"`
	After      jsonSummarySection `json:"after_filtering"`
}

// jsonInsertSize mirrors upstream's insert_size block (jsonreporter.cpp:121-134),
// emitted only for paired-end input.
type jsonInsertSize struct {
	Peak      int     `json:"peak"`
	Unknown   int64   `json:"unknown"`
	Histogram []int64 `json:"histogram"`
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
	TotalReads    int64                `json:"total_reads"`
	TotalBases    int64                `json:"total_bases"`
	Q20Bases      int64                `json:"q20_bases"`
	Q30Bases      int64                `json:"q30_bases"`
	Q40Bases      int64                `json:"q40_bases"`
	TotalCycles   int                  `json:"total_cycles"`
	QualityCurves map[string][]float64 `json:"quality_curves,omitempty"`
	ContentCurves map[string][]float64 `json:"content_curves,omitempty"`
	// KmerCount is the 5-mer histogram (1024 entries) matching upstream's
	// kmer_count block. It is computed for both the before and after streams.
	KmerCount map[string]int64 `json:"kmer_count,omitempty"`
	// OverrepresentedSequences is emitted for both streams when -p/-P is on.
	OverrepresentedSequences map[string]int64 `json:"overrepresented_sequences,omitempty"`
	// LengthDistribution is a Go extension (upstream's JSON does not include
	// it; the data is carried in the HTML report there). Kept for downstream
	// consumers; omitted when empty.
	LengthDistribution map[string]int64 `json:"length_distribution,omitempty"`
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
	if stats.MergeEnabled {
		// Upstream renames the after-filtering block to `merged_and_filtered`
		// and omits `read2_after_filtering` when -m/--merge is set
		// (jsonreporter.cpp:158-167). Drop read2 here; the read1 key is
		// renamed on the marshalled bytes below to preserve key order.
		report.Read2AfterFilter = nil
	}
	data, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return fmt.Errorf("marshal JSON report: %w", err)
	}
	if stats.MergeEnabled {
		data = bytes.Replace(data,
			[]byte(`"read1_after_filtering":`),
			[]byte(`"merged_and_filtered":`), 1)
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
			Sequencing: sequencingInfo(stats, pe),
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
		if is := buildInsertSize(stats); is != nil {
			r.InsertSize = is
		}
	}
	return r
}

// sequencingInfo builds the summary.sequencing descriptor, deterministic given
// the input cycle counts, matching jsonreporter.cpp:28-33. For PE it reports
// both read1 and read2 cycle counts.
func sequencingInfo(s *ProcessStats, pe bool) string {
	c1 := 0
	if s.curvesBefore[0] != nil {
		c1 = s.curvesBefore[0].cycles()
	}
	if pe {
		c2 := 0
		if s.curvesBefore[1] != nil {
			c2 = s.curvesBefore[1].cycles()
		}
		return fmt.Sprintf("paired end (%d cycles + %d cycles)", c1, c2)
	}
	return fmt.Sprintf("single end (%d cycles)", c1)
}

// buildInsertSize emits upstream's insert_size block from the accumulated
// insert-size histogram. Returns nil when no histogram was collected.
func buildInsertSize(s *ProcessStats) *jsonInsertSize {
	if len(s.InsertHist) == 0 {
		return nil
	}
	max := len(s.InsertHist) - 1 // last bucket is "unknown"
	peak := 0
	peakCount := int64(-1)
	for d := 0; d < max; d++ {
		if s.InsertHist[d] > peakCount {
			peakCount = s.InsertHist[d]
			peak = d
		}
	}
	hist := make([]int64, max)
	copy(hist, s.InsertHist[:max])
	return &jsonInsertSize{
		Peak:      peak,
		Unknown:   s.InsertHist[max],
		Histogram: hist,
	}
}

// buildReadStats fills jsonReadStats for the given read index. before
// selects the BEFORE-filtering stream; otherwise AFTER. The per-cycle
// quality_curves / content_curves, the kmer_count histogram, the q40 total
// and (before) the overrepresented sequences are taken from the readCurves
// accumulator, reproducing upstream Stats::reportJson for BOTH streams.
func buildReadStats(s *ProcessStats, readIdx int, before bool) jsonReadStats {
	hist := s.LengthHistBefore[readIdx]
	curves := s.curvesBefore[readIdx]
	streamIdx := 0
	if !before {
		hist = s.LengthHistAfter[readIdx]
		curves = s.curvesAfter[readIdx]
		streamIdx = 1
	}
	out := jsonReadStats{
		TotalReads:         readTotalReads(s, readIdx, before),
		TotalBases:         readTotalBases(s, readIdx, before),
		Q20Bases:           s.Q20ByRead[streamIdx][readIdx],
		Q30Bases:           s.Q30ByRead[streamIdx][readIdx],
		LengthDistribution: lengthHistToJSON(hist),
	}
	if curves != nil {
		out.Q40Bases = curves.q40
		out.TotalCycles = curves.cycles()
		out.QualityCurves = curves.qualityCurves()
		out.ContentCurves = curves.contentCurves()
		out.KmerCount = curves.kmerCounts()
	}
	if before {
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
