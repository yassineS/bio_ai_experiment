package prinseq

import (
	"bytes"
	"fmt"
	"io"
	"time"
)

// BenchmarkResult holds benchmark timing and performance data
type BenchmarkResult struct {
	Operation     string        `json:"operation"`
	Duration      time.Duration `json:"duration_ns"`
	DurationMs    float64       `json:"duration_ms"`
	ThroughputMBs float64       `json:"throughput_mbs,omitempty"`
	ReadsPerSec   float64       `json:"reads_per_sec,omitempty"`
	NumReads      int           `json:"num_reads,omitempty"`
	FileSize      int64         `json:"file_size_bytes,omitempty"`
}

// BenchmarkStats measures the time to calculate statistics
func BenchmarkStats(reader io.Reader, isFastq bool) (*BenchmarkResult, *Stats, error) {
	// Read all data into memory first to isolate computation time
	data, err := io.ReadAll(reader)
	if err != nil {
		return nil, nil, fmt.Errorf("error reading data: %w", err)
	}

	start := time.Now()
	stats, err := CalculateEnhancedStats(bytes.NewReader(data), isFastq)
	duration := time.Since(start)

	if err != nil {
		return nil, nil, err
	}

	result := &BenchmarkResult{
		Operation:  "stats",
		Duration:   duration,
		DurationMs: float64(duration.Milliseconds()),
		NumReads:   stats.NumReads,
		FileSize:   int64(len(data)),
	}

	// Calculate throughput
	if duration.Seconds() > 0 {
		result.ThroughputMBs = float64(len(data)) / (1024 * 1024) / duration.Seconds()
		result.ReadsPerSec = float64(stats.NumReads) / duration.Seconds()
	}

	return result, stats, nil
}

// BenchmarkFilter measures the time to filter sequences
func BenchmarkFilter(reader io.Reader, isFastq bool, opts FilterOptions) (*BenchmarkResult, error) {
	// Read all data into memory first
	data, err := io.ReadAll(reader)
	if err != nil {
		return nil, fmt.Errorf("error reading data: %w", err)
	}

	// Use a buffer to capture output without I/O overhead
	var output bytes.Buffer

	start := time.Now()
	err = Filter(bytes.NewReader(data), &output, isFastq, opts)
	duration := time.Since(start)

	if err != nil {
		return nil, err
	}

	result := &BenchmarkResult{
		Operation:  "filter",
		Duration:   duration,
		DurationMs: float64(duration.Milliseconds()),
		FileSize:   int64(len(data)),
	}

	// Calculate throughput
	if duration.Seconds() > 0 {
		result.ThroughputMBs = float64(len(data)) / (1024 * 1024) / duration.Seconds()
	}

	return result, nil
}

// RunBenchmarkSuite runs a comprehensive benchmark suite
func RunBenchmarkSuite(reader io.Reader, isFastq bool) ([]*BenchmarkResult, error) {
	// Read data once
	data, err := io.ReadAll(reader)
	if err != nil {
		return nil, fmt.Errorf("error reading data: %w", err)
	}

	results := make([]*BenchmarkResult, 0)

	// Benchmark basic stats
	result1, _, err := BenchmarkStats(bytes.NewReader(data), isFastq)
	if err != nil {
		return nil, err
	}
	results = append(results, result1)

	// Benchmark enhanced stats
	start := time.Now()
	_, err = CalculateEnhancedStats(bytes.NewReader(data), isFastq)
	duration := time.Since(start)
	if err == nil {
		results = append(results, &BenchmarkResult{
			Operation:  "enhanced_stats",
			Duration:   duration,
			DurationMs: float64(duration.Milliseconds()),
			FileSize:   int64(len(data)),
		})
	}

	// Benchmark different filter operations
	filterOpts := []struct {
		name string
		opts FilterOptions
	}{
		{"filter_length", FilterOptions{MinLen: 50}},
		{"filter_gc", FilterOptions{MinGC: 40, MaxGC: 60}},
		{"filter_quality", FilterOptions{MinQualMean: 20}},
		{"filter_combined", FilterOptions{MinLen: 50, MinGC: 40, MaxGC: 60, MinQualMean: 20}},
	}

	for _, f := range filterOpts {
		var output bytes.Buffer
		start := time.Now()
		err := Filter(bytes.NewReader(data), &output, isFastq, f.opts)
		duration := time.Since(start)

		if err == nil {
			result := &BenchmarkResult{
				Operation:  f.name,
				Duration:   duration,
				DurationMs: float64(duration.Milliseconds()),
				FileSize:   int64(len(data)),
			}
			if duration.Seconds() > 0 {
				result.ThroughputMBs = float64(len(data)) / (1024 * 1024) / duration.Seconds()
			}
			results = append(results, result)
		}
	}

	return results, nil
}

// FormatBenchmarkResults formats benchmark results for display
func FormatBenchmarkResults(results []*BenchmarkResult) string {
	var buf bytes.Buffer

	fmt.Fprintf(&buf, "Benchmark Results\n")
	fmt.Fprintf(&buf, "=================\n\n")
	fmt.Fprintf(&buf, "%-20s | %12s | %12s | %15s\n", "Operation", "Duration", "Throughput", "Reads/sec")
	fmt.Fprintf(&buf, "%s\n", "--------------------------------------------------------------------------------")

	for _, r := range results {
		throughput := "-"
		if r.ThroughputMBs > 0 {
			throughput = fmt.Sprintf("%.2f MB/s", r.ThroughputMBs)
		}

		readsPerSec := "-"
		if r.ReadsPerSec > 0 {
			readsPerSec = fmt.Sprintf("%.0f reads/s", r.ReadsPerSec)
		}

		fmt.Fprintf(&buf, "%-20s | %10.2f ms | %12s | %15s\n",
			r.Operation, r.DurationMs, throughput, readsPerSec)
	}

	return buf.String()
}
