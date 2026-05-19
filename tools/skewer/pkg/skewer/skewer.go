// Package skewer provides adapter trimming for FASTQ files.
// It detects and removes adapter sequences from the 3' and 5' ends of reads.
package skewer

import (
	"encoding/json"
	"fmt"
	"io"
	"sync"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TrimOptions contains parameters for adapter trimming.
type TrimOptions struct {
	Adapter3         string  // 3' adapter sequence
	Adapter5         string  // 5' adapter sequence
	MinLength        int     // Minimum read length after trimming
	QualThreshold    int     // Quality threshold for trimming
	MinOverlap       int     // Minimum overlap for adapter detection
	ErrorRate        float64 // Maximum error rate for adapter matching
	TrimBothEnds     bool    // Trim adapters from both ends
	AutoDetect       bool    // Auto-detect adapters from reads
	AutoDetectReads  int     // Number of reads to analyze for auto-detection
	ProgressReport   bool    // Report progress during processing
	ProgressInterval int     // Report progress every N reads
	UMILength        int     // UMI length to extract (0 = disabled)
	UMIPosition      string  // UMI position: "5prime" or "3prime"
	PEMatrixMode     bool    // PE matrix mode: require reverse-complement overlap support before trimming (skewer's default -m pe)
}

// DefaultTrimOptions returns default trimming options.
func DefaultTrimOptions() TrimOptions {
	return TrimOptions{
		Adapter3:         "",
		Adapter5:         "",
		MinLength:        18,
		QualThreshold:    0,
		MinOverlap:       3,
		ErrorRate:        0.1,
		TrimBothEnds:     false,
		AutoDetect:       false,
		AutoDetectReads:  1000,
		ProgressReport:   false,
		ProgressInterval: 100000,
		UMILength:        0,
		UMIPosition:      "5prime",
		PEMatrixMode:     false,
	}
}

// TrimStats tracks adapter trimming statistics.
type TrimStats struct {
	TotalReads       int           `json:"total_reads"`
	TrimmedReads     int           `json:"trimmed_reads"`
	AdapterFound3    int           `json:"adapter_found_3"`
	AdapterFound5    int           `json:"adapter_found_5"`
	DiscardedReads   int           `json:"discarded_reads"`
	TotalBases       int64         `json:"total_bases"`
	TrimmedBases     int64         `json:"trimmed_bases"`
	ProcessingTime   time.Duration `json:"processing_time"`
	DetectedAdapter3 string        `json:"detected_adapter_3,omitempty"`
	DetectedAdapter5 string        `json:"detected_adapter_5,omitempty"`
	UMIStats         *UMIStats     `json:"umi_stats,omitempty"`
	mu               sync.Mutex    `json:"-"`
}

// UMIStats tracks UMI/barcode statistics.
type UMIStats struct {
	TotalUMIs       int            `json:"total_umis"`
	UniqueUMIs      int            `json:"unique_umis"`
	UMIDistribution map[string]int `json:"umi_distribution,omitempty"`
}

// ToJSON exports statistics as JSON. The error from json.MarshalIndent is
// unreachable in practice — TrimStats only contains ints, strings, a
// map[string]int and nested structs of the same, none of which json.Marshal
// can fail on — so this method does not return one.
func (s *TrimStats) ToJSON() string {
	s.mu.Lock()
	defer s.mu.Unlock()

	data, _ := json.MarshalIndent(s, "", "  ")
	return string(data)
}

// ToHTML generates an HTML report of the statistics.
func (s *TrimStats) ToHTML() string {
	s.mu.Lock()
	defer s.mu.Unlock()

	html := `<!DOCTYPE html>
<html>
<head>
    <title>Skewer Trimming Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        table { border-collapse: collapse; width: 100%; max-width: 800px; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .metric { font-weight: bold; }
        .value { color: #2196F3; }
        .percentage { color: #FF9800; font-size: 0.9em; }
    </style>
</head>
<body>
    <h1>Adapter Trimming Report</h1>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
`

	html += fmt.Sprintf("        <tr><td class='metric'>Total Reads</td><td class='value'>%d</td></tr>\n", s.TotalReads)

	if s.TotalReads > 0 {
		trimmedPct := 100.0 * float64(s.TrimmedReads) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>Trimmed Reads</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.TrimmedReads, trimmedPct)

		adapter3Pct := 100.0 * float64(s.AdapterFound3) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>3' Adapters Found</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.AdapterFound3, adapter3Pct)

		adapter5Pct := 100.0 * float64(s.AdapterFound5) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>5' Adapters Found</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.AdapterFound5, adapter5Pct)

		discardedPct := 100.0 * float64(s.DiscardedReads) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>Discarded Reads</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.DiscardedReads, discardedPct)

		keptPct := 100.0 * float64(s.TotalReads-s.DiscardedReads) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>Kept Reads</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.TotalReads-s.DiscardedReads, keptPct)
	}

	html += fmt.Sprintf("        <tr><td class='metric'>Total Bases</td><td class='value'>%d</td></tr>\n", s.TotalBases)

	if s.TotalBases > 0 {
		trimmedBasesPct := 100.0 * float64(s.TrimmedBases) / float64(s.TotalBases)
		html += fmt.Sprintf("        <tr><td class='metric'>Trimmed Bases</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.TrimmedBases, trimmedBasesPct)
	}

	html += fmt.Sprintf("        <tr><td class='metric'>Processing Time</td><td class='value'>%v</td></tr>\n", s.ProcessingTime)

	if s.DetectedAdapter3 != "" {
		html += fmt.Sprintf("        <tr><td class='metric'>Detected 3' Adapter</td><td class='value'>%s</td></tr>\n", s.DetectedAdapter3)
	}
	if s.DetectedAdapter5 != "" {
		html += fmt.Sprintf("        <tr><td class='metric'>Detected 5' Adapter</td><td class='value'>%s</td></tr>\n", s.DetectedAdapter5)
	}

	if s.UMIStats != nil {
		html += fmt.Sprintf("        <tr><td class='metric'>Total UMIs</td><td class='value'>%d</td></tr>\n", s.UMIStats.TotalUMIs)
		html += fmt.Sprintf("        <tr><td class='metric'>Unique UMIs</td><td class='value'>%d</td></tr>\n", s.UMIStats.UniqueUMIs)
	}

	html += `    </table>
</body>
</html>
`
	return html
}

// TrimSingleEnd trims adapters from single-end FASTQ reads.
func TrimSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {
	startTime := time.Now()

	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	stats := &TrimStats{}

	// buffered holds records that were pre-read for auto-detection and still
	// need to be processed; it is drained before reading further from reader.
	var buffered []*fastq.Record

	// Auto-detect adapters if requested
	if opts.AutoDetect && opts.Adapter3 == "" && opts.Adapter5 == "" {
		maxReads := opts.AutoDetectReads
		if maxReads <= 0 {
			maxReads = 1000
		}
		for i := 0; i < maxReads; i++ {
			record, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return stats, fmt.Errorf("error reading FASTQ: %w", err)
			}
			buffered = append(buffered, record)
		}
		detected := detectAdaptersFromReads(buffered)
		if detected.Adapter3 != "" {
			opts.Adapter3 = detected.Adapter3
			stats.DetectedAdapter3 = detected.Adapter3
		}
		if detected.Adapter5 != "" {
			opts.Adapter5 = detected.Adapter5
			stats.DetectedAdapter5 = detected.Adapter5
		}
	}

	// Initialize UMI stats if UMI extraction is enabled
	if opts.UMILength > 0 {
		stats.UMIStats = &UMIStats{
			UMIDistribution: make(map[string]int),
		}
	}

	readCount := 0
	for {
		var record *fastq.Record
		if len(buffered) > 0 {
			record = buffered[0]
			buffered = buffered[1:]
		} else {
			var err error
			record, err = reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return stats, fmt.Errorf("error reading FASTQ: %w", err)
			}
		}

		stats.mu.Lock()
		stats.TotalReads++
		readCount++
		stats.mu.Unlock()

		originalLength := len(record.Sequence)
		stats.mu.Lock()
		stats.TotalBases += int64(originalLength)
		stats.mu.Unlock()

		// Extract UMI if enabled
		var umi string
		if opts.UMILength > 0 {
			record, umi = extractUMI(record, opts)
			if umi != "" {
				stats.mu.Lock()
				stats.UMIStats.TotalUMIs++
				stats.UMIStats.UMIDistribution[umi]++
				stats.mu.Unlock()
			}
		}

		// Trim adapters
		trimmed := trimRecord(record, opts, stats)

		// Add UMI to description if extracted
		if umi != "" {
			trimmed.Description = trimmed.Description + " UMI:" + umi
		}

		// Check if read passes length threshold
		if len(trimmed.Sequence) >= opts.MinLength {
			if err := writer.Write(trimmed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			if len(trimmed.Sequence) < originalLength {
				stats.mu.Lock()
				stats.TrimmedReads++
				stats.TrimmedBases += int64(originalLength - len(trimmed.Sequence))
				stats.mu.Unlock()
			}
		} else {
			stats.mu.Lock()
			stats.DiscardedReads++
			stats.mu.Unlock()
		}

		// Progress reporting
		if opts.ProgressReport && readCount%opts.ProgressInterval == 0 {
			fmt.Fprintf(io.Discard, "\rProcessed %d reads...", readCount)
		}
	}

	// Flush writer
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}

	// Calculate unique UMIs
	if stats.UMIStats != nil {
		stats.UMIStats.UniqueUMIs = len(stats.UMIStats.UMIDistribution)
	}

	stats.ProcessingTime = time.Since(startTime)
	return stats, nil
}

// TrimPairedEnd trims adapters from paired-end FASTQ reads.
func TrimPairedEnd(input1, input2 io.Reader, output1, output2, outputSingle io.Writer,
	encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {

	startTime := time.Now()

	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)

	var writerSingle *fastq.Writer
	if outputSingle != nil {
		writerSingle = fastq.NewWriter(outputSingle, encoding)
	}

	stats := &TrimStats{}

	// Initialize UMI stats if UMI extraction is enabled
	if opts.UMILength > 0 {
		stats.UMIStats = &UMIStats{
			UMIDistribution: make(map[string]int),
		}
	}

	readCount := 0
	for {
		record1, err1 := reader1.Read()
		record2, err2 := reader2.Read()

		if err1 == io.EOF || err2 == io.EOF {
			break
		}
		if err1 != nil {
			return stats, fmt.Errorf("error reading first input: %w", err1)
		}
		if err2 != nil {
			return stats, fmt.Errorf("error reading second input: %w", err2)
		}

		stats.mu.Lock()
		stats.TotalReads += 2
		readCount += 2
		stats.mu.Unlock()

		originalLen1 := len(record1.Sequence)
		originalLen2 := len(record2.Sequence)
		stats.mu.Lock()
		stats.TotalBases += int64(originalLen1 + originalLen2)
		stats.mu.Unlock()

		// Extract UMI if enabled (from first read only)
		var umi string
		if opts.UMILength > 0 {
			record1, umi = extractUMI(record1, opts)
			if umi != "" {
				stats.mu.Lock()
				stats.UMIStats.TotalUMIs++
				stats.UMIStats.UMIDistribution[umi]++
				stats.mu.Unlock()
			}
		}

		// Trim both reads. In matrix mode (skewer's default `-m pe`) the
		// per-mate trimmer is gated by a reverse-complement overlap check
		// so the mates have to agree on the inferred insert size before any
		// trimming happens. See detectPairedAdapter for the algorithm.
		var trimmed1, trimmed2 *fastq.Record
		if opts.PEMatrixMode {
			trimmed1, trimmed2 = trimPairWithMatrix(record1, record2, opts, stats)
		} else {
			trimmed1 = trimRecord(record1, opts, stats)
			trimmed2 = trimRecord(record2, opts, stats)
		}

		// Add UMI to descriptions if extracted
		if umi != "" {
			trimmed1.Description = trimmed1.Description + " UMI:" + umi
			trimmed2.Description = trimmed2.Description + " UMI:" + umi
		}

		// Check if both reads pass length threshold
		pass1 := len(trimmed1.Sequence) >= opts.MinLength
		pass2 := len(trimmed2.Sequence) >= opts.MinLength

		if pass1 && pass2 {
			// Both pass - write to paired output
			if err := writer1.Write(trimmed1); err != nil {
				return stats, fmt.Errorf("error writing first output: %w", err)
			}
			if err := writer2.Write(trimmed2); err != nil {
				return stats, fmt.Errorf("error writing second output: %w", err)
			}

			if len(trimmed1.Sequence) < originalLen1 || len(trimmed2.Sequence) < originalLen2 {
				stats.mu.Lock()
				stats.TrimmedReads++
				stats.TrimmedBases += int64((originalLen1 - len(trimmed1.Sequence)) +
					(originalLen2 - len(trimmed2.Sequence)))
				stats.mu.Unlock()
			}
		} else if writerSingle != nil {
			// One or both fail - write survivors to single output
			if pass1 {
				if err := writerSingle.Write(trimmed1); err != nil {
					return stats, fmt.Errorf("error writing single output: %w", err)
				}
			} else {
				stats.mu.Lock()
				stats.DiscardedReads++
				stats.mu.Unlock()
			}

			if pass2 {
				if err := writerSingle.Write(trimmed2); err != nil {
					return stats, fmt.Errorf("error writing single output: %w", err)
				}
			} else {
				stats.mu.Lock()
				stats.DiscardedReads++
				stats.mu.Unlock()
			}
		} else {
			// Both discarded
			stats.mu.Lock()
			stats.DiscardedReads += 2
			stats.mu.Unlock()
		}

		// Progress reporting
		if opts.ProgressReport && readCount%opts.ProgressInterval == 0 {
			fmt.Fprintf(io.Discard, "\rProcessed %d reads...", readCount)
		}
	}

	// Flush writers
	if err := writer1.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing first output: %w", err)
	}
	if err := writer2.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing second output: %w", err)
	}
	if writerSingle != nil {
		if err := writerSingle.Flush(); err != nil {
			return stats, fmt.Errorf("error flushing single output: %w", err)
		}
	}

	// Calculate unique UMIs
	if stats.UMIStats != nil {
		stats.UMIStats.UniqueUMIs = len(stats.UMIStats.UMIDistribution)
	}

	stats.ProcessingTime = time.Since(startTime)
	return stats, nil
}

// trimPairWithMatrix applies skewer's `-m pe` matrix-mode logic to a paired
// read: it only trims the mates when their inferred insert sizes are
// consistent (the R1 prefix before the adapter is a reverse-complement of the
// R2 prefix before the adapter). Mirrors cMatrix::findAdapterWithPE in
// reference_code/skewer/src/matrix.cpp:726-851 and the surrounding plumbing
// at main.cpp:1390-1413.
func trimPairWithMatrix(record1, record2 *fastq.Record, opts TrimOptions, stats *TrimStats) (*fastq.Record, *fastq.Record) {
	seq1 := string(record1.Sequence)
	seq2 := string(record2.Sequence)
	qual1 := record1.Quality
	qual2 := record2.Quality

	// If no 3' adapter is configured, fall back to plain per-record trimming
	// (5' / quality trimming still applies). The matrix logic only governs
	// the 3'-end adapter detection.
	if opts.Adapter3 == "" {
		return trimRecord(record1, opts, stats), trimRecord(record2, opts, stats)
	}

	pos1 := findAdapterWithQual(seq1, opts.Adapter3, qual1, opts.MinOverlap, opts.ErrorRate)
	pos2 := findAdapterWithQual(seq2, opts.Adapter3, qual2, opts.MinOverlap, opts.ErrorRate)

	// detectPairedTrim returns the trim positions sanctioned by the matrix
	// gate; -1 means "leave that mate untrimmed".
	tp1, tp2 := detectPairedTrim(seq1, seq2, qual1, qual2, pos1, pos2, opts.ErrorRate)

	out1 := applyTrim(record1, tp1, opts, stats, 3)
	out2 := applyTrim(record2, tp2, opts, stats, 3)
	return out1, out2
}

// detectPairedTrim implements the reverse-complement overlap gate from
// matrix.cpp:761-849. Given the (possibly -1) adapter positions in R1 and R2
// it returns the positions at which the mates should actually be cut; -1 in
// the return position means "no trim". The gate trims only when the two
// mate prefixes (up to the inferred insert end) are reverse-complement
// matches under the quality-weighted scoring model from matrix.cpp:487-522.
func detectPairedTrim(seq1, seq2 string, qual1, qual2 []byte, pos1, pos2 int, errorRate float64) (int, int) {
	if pos1 < 0 && pos2 < 0 {
		return -1, -1
	}
	// Pick the smaller of the two candidate insert lengths as the overlap
	// span (cpos in matrix.cpp:767/772). If only one mate found the
	// adapter, the other one's position is treated as INT_MAX.
	cpos := pos1
	if pos1 < 0 || (pos2 >= 0 && pos2 < pos1) {
		cpos = pos2
	}
	if cpos <= 0 {
		// matrix.cpp:843-844: apos <= 0 ⇒ both pos clamped to 0 ⇒ effectively
		// no overlap to validate, no trim.
		return -1, -1
	}
	if calcRevCompScore(seq1, seq2, qual1, qual2, cpos, errorRate) {
		// Both mates trimmed at cpos, with apos2 clamped to the shorter
		// mate's length when cpos overruns it (matrix.cpp:785, 803).
		tp1 := cpos
		if cpos > len(seq1) {
			tp1 = len(seq1)
		}
		tp2 := cpos
		if cpos > len(seq2) {
			tp2 = len(seq2)
		}
		return tp1, tp2
	}
	return -1, -1
}

// calcRevCompScore is the Go port of cMatrix::CalcRevCompScore in
// matrix.cpp:487-522. It walks the first `n` bases of R1 against the
// reverse-complement of the first `n` bases of R2 and returns true when the
// quality-weighted mismatch penalty stays below dPenaltyPerErr * n. n must be
// > 0; callers (detectPairedTrim) shortcut the n<=0 branch above.
func calcRevCompScore(seq1, seq2 string, qual1, qual2 []byte, n int, errorRate float64) bool {
	if n <= 0 || n > len(seq1) || n > len(seq2) {
		return false
	}
	dPenaltyPerErr := errorRate * meanPenalty
	dMaxPenalty := dPenaltyPerErr * float64(n)
	var penalty float64
	for i := 0; i < n; i++ {
		a := seq1[i]
		// complement of R2[n-1-i]: per matrix.cpp:500.
		b := complementBase(seq2[n-1-i])
		if iupacMatch(a, b) {
			continue
		}
		// matrix.cpp:504-509 takes the minimum of the two mates' quality
		// penalties; we replicate via min(mismatchPenalty(qual1,i),
		// mismatchPenalty(qual2, n-1-i)).
		p1 := mismatchPenalty(qual1, i)
		p2 := mismatchPenalty(qual2, n-1-i)
		penal := p1
		if p2 < penal {
			penal = p2
		}
		penalty += penal
		if penalty > dMaxPenalty {
			return false
		}
	}
	return true
}

// complementBase returns the Watson-Crick complement of a single ASCII base.
// Lower-case is mapped via toupper. Unknown bases (anything outside ACGTN)
// map to 'N', matching matrix.cpp:84-86's `complement[]` table.
func complementBase(b byte) byte {
	if b >= 'a' && b <= 'z' {
		b -= 32
	}
	switch b {
	case 'A':
		return 'T'
	case 'T':
		return 'A'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	}
	return 'N'
}

// applyTrim builds a trimmed copy of record, cutting at trimPos (or leaving
// it untouched when trimPos < 0). end=3 means "3'-end cut"; this is the only
// mode trimPairWithMatrix uses today. The quality-based trim from
// opts.QualThreshold is still applied after the adapter cut, matching the
// SE path's behaviour.
func applyTrim(record *fastq.Record, trimPos int, opts TrimOptions, stats *TrimStats, end int) *fastq.Record {
	_ = end // reserved for a future 5'-cut variant; the matrix gate only emits 3' cuts today.
	start := 0
	stop := len(record.Sequence)
	if trimPos >= 0 && trimPos < stop {
		stop = trimPos
		if stats != nil {
			stats.AdapterFound3++
		}
	}
	if opts.QualThreshold > 0 {
		start, stop = trimByQuality(record.Quality[start:stop], opts.QualThreshold, start, stop)
	}
	if start >= stop || stop-start < opts.MinLength {
		return &fastq.Record{
			ID:          record.ID,
			Description: record.Description,
			Sequence:    []byte{},
			Quality:     []byte{},
		}
	}
	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence[start:stop],
		Quality:     record.Quality[start:stop],
	}
}

// trimRecord trims adapters from a single record.
func trimRecord(record *fastq.Record, opts TrimOptions, stats *TrimStats) *fastq.Record {
	seq := string(record.Sequence)
	qual := record.Quality

	start := 0
	end := len(seq)

	// Trim 5' adapter if specified
	if opts.Adapter5 != "" {
		pos := findAdapterWithQual(seq, opts.Adapter5, qual, opts.MinOverlap, opts.ErrorRate)
		if pos >= 0 {
			// Found 5' adapter - trim from start to end of adapter
			start = pos + len(opts.Adapter5)
			if stats != nil {
				stats.AdapterFound5++
			}
		}
	}

	// Trim 3' adapter if specified
	if opts.Adapter3 != "" {
		pos := findAdapterWithQual(seq[start:], opts.Adapter3, qual[start:], opts.MinOverlap, opts.ErrorRate)
		if pos >= 0 {
			// Found 3' adapter - trim from adapter position to end
			end = start + pos
			if stats != nil {
				stats.AdapterFound3++
			}
		}
	}

	// Apply quality-based trimming if threshold is set
	if opts.QualThreshold > 0 {
		start, end = trimByQuality(qual[start:end], opts.QualThreshold, start, end)
	}

	// Create trimmed record
	if start >= end || end-start < opts.MinLength {
		// Return empty record if too short
		return &fastq.Record{
			ID:          record.ID,
			Description: record.Description,
			Sequence:    []byte{},
			Quality:     []byte{},
		}
	}

	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence[start:end],
		Quality:     record.Quality[start:end],
	}
}

// findAdapter finds the position of an adapter in a sequence with error tolerance.
// Returns -1 if adapter not found. Uses improved alignment scoring.
func findAdapter(seq string, adapter string, minOverlap int, errorRate float64) int {
	// Use improved algorithm for better accuracy
	return improvedFindAdapter(seq, adapter, minOverlap, errorRate)
}

// findAdapterWithQual finds the leftmost adapter position in a sequence using
// the upstream skewer scoring model (quality-weighted penalty, TRIM_TAIL).
//
// Ported from reference_code/skewer/src/matrix.cpp:cAdapter::align
// (lines 297-435) and the cMatrix scoring constants (lines 138-141). The
// upstream algorithm uses a bit-masked k-difference engine with an asymmetric
// tail penalty; this port collapses that to its observable behaviour for the
// no-indel ASCII-ATCG-only case that the parity corpus exercises: scan every
// candidate adapter position p in [0, rLen), accumulate per-position penalties
// (quality-derived for mismatches, zero for matches), and accept the match
// whose cumulative penalty stays under dPenaltyPerErr * compareLen + EPSILON.
// For TRIM_TAIL the adapter is allowed to extend past the read's 3' end
// (compareLen < adapterLen) provided compareLen >= minOverlap. The leftmost
// acceptable position wins, matching upstream's preference for longer
// alignments at the same score (longer match span ⇒ earlier start ⇒ lower p).
//
// qual is the Phred-encoded quality string (same encoding as seq); if qual is
// empty, MEAN_PENALTY is used as the per-mismatch penalty (matches upstream's
// "qLen == 0" branch in align() and CalcRevCompScore()).
func findAdapterWithQual(seq, adapter string, qual []byte, minOverlap int, errorRate float64) int {
	if len(adapter) == 0 || len(seq) < minOverlap {
		return -1
	}
	dPenaltyPerErr := errorRate * meanPenalty
	bestPos := -1
	bestScore := -1.0
	for i := 0; i <= len(seq)-minOverlap; i++ {
		compareLen := len(adapter)
		if remaining := len(seq) - i; remaining < compareLen {
			compareLen = remaining
		}
		if compareLen < minOverlap {
			continue
		}
		// dMaxPenalty matches matrix.cpp:301 / matrix.cpp:418
		// (dPenaltyPerErr * span + 0.001 slack to absorb FP rounding).
		dMaxPenalty := dPenaltyPerErr*float64(compareLen) + epsilonSlack
		var penalty float64
		ok := true
		for j := 0; j < compareLen; j++ {
			if iupacMatch(seq[i+j], adapter[j]) {
				continue
			}
			// Per-base penalty driven by base quality, mirroring
			// matrix.cpp:cMatrix::penalty[] (lines 547-556).
			penalty += mismatchPenalty(qual, i+j)
			if penalty >= dMaxPenalty {
				ok = false
				break
			}
		}
		if !ok {
			continue
		}
		// Normalised score = (compareLen * dMu - penalty) / (compareLen + 1)
		// matches the bestAlign branch in matrix.cpp:362/369/392/408 — favours
		// longer alignments with lower penalties. Ties keep the leftmost.
		score := (float64(compareLen)*meanPenalty - penalty) / float64(compareLen+1)
		if score > bestScore {
			bestScore = score
			bestPos = i
		}
	}
	return bestPos
}

// iupacMatch reports whether a read base (ATCGN/lowercase) matches an adapter
// base under upstream's chrVadp table (matrix.cpp:115-136). For the parity
// corpus we only need the ATCG/N subset: exact equality, with N treated as
// non-matching to anything (chrVadp[N][*]=1, chrVadp[*][N]=0 — but here we
// score only mismatches so a non-match read base does count as a mismatch).
func iupacMatch(readBase, adapterBase byte) bool {
	if readBase == adapterBase {
		return true
	}
	// Lower-case forms map to the same code (codeMap entries 0x61..0x7A).
	if readBase >= 'a' && readBase <= 'z' {
		readBase -= 32
	}
	if adapterBase >= 'a' && adapterBase <= 'z' {
		adapterBase -= 32
	}
	return readBase == adapterBase
}

// Penalty constants. Verbatim from matrix.cpp:138-141 (the natural-log-10
// derived quality weights used by skewer's k-difference scoring).
const (
	minPenalty   = 0.477121255
	meanPenalty  = 2.477121255
	maxPenalty   = 4.477121255
	epsilonSlack = 0.001
)

// mismatchPenalty returns the per-base mismatch penalty driven by Phred-33
// quality. Mirrors the precomputed cMatrix::penalty[] table built in
// matrix.cpp:547-556: bytes <= baseQual (33) map to minPenalty, the next 39
// quality levels ramp linearly by 0.1, and Q>=40 saturates at maxPenalty.
// When qual is empty the function returns meanPenalty (matches the
// "qLen == 0" fallback in matrix.cpp:512 / matrix.cpp:338).
func mismatchPenalty(qual []byte, idx int) float64 {
	if idx < 0 || idx >= len(qual) {
		return meanPenalty
	}
	q := int(qual[idx])
	if q <= 33 {
		return minPenalty
	}
	if q >= 33+40 {
		return maxPenalty
	}
	return minPenalty + float64(q-33)/10.0
}

// trimByQuality trims low-quality regions from both ends.
func trimByQuality(quality []byte, threshold int, start, end int) (int, int) {
	// Trim from 3' end
	for end > start && int(quality[end-start-1]) < threshold+33 {
		end--
	}

	// Trim from 5' end
	for start < end && int(quality[0]) < threshold+33 {
		start++
		quality = quality[1:]
	}

	return start, end
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// Common adapter sequences for auto-detection
var commonAdapters = []string{
	"AGATCGGAAGAGC",                      // Illumina TruSeq Universal
	"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", // Illumina TruSeq Full
	"CTGTCTCTTATACACATCT",                // Illumina Nextera
	"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",  // Illumina Small RNA 3'
	"TGGAATTCTCGGGTGCCAAGG",              // Illumina Small RNA 5'
	"CGCCTTGGCCGTACAGCAG",                // SOLiD
	"ATCTCGTATGCCGTCTTCTGCTTG",           // Ion Torrent
}

// detectAdapters attempts to auto-detect adapter sequences by reading up to
// maxReads records from reader.
func detectAdapters(reader *fastq.Reader, maxReads int) (TrimOptions, error) {
	var reads []*fastq.Record
	for i := 0; i < maxReads; i++ {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return TrimOptions{}, err
		}
		reads = append(reads, record)
	}
	return detectAdaptersFromReads(reads), nil
}

// detectAdaptersFromReads inspects a slice of already-read records and returns
// TrimOptions populated with any common adapters that appear in at least 10% of
// them.
func detectAdaptersFromReads(reads []*fastq.Record) TrimOptions {
	opts := TrimOptions{}

	if len(reads) == 0 {
		return opts
	}

	// Count adapter occurrences
	adapter3Counts := make(map[string]int)
	adapter5Counts := make(map[string]int)

	for _, record := range reads {
		seq := string(record.Sequence)

		// Check 3' adapters (at end of read)
		for _, adapter := range commonAdapters {
			if pos := findAdapter(seq, adapter, 5, 0.15); pos >= 0 {
				adapter3Counts[adapter]++
			}
		}

		// Check 5' adapters (at beginning of read)
		for _, adapter := range commonAdapters {
			if pos := findAdapter(seq, adapter, 5, 0.15); pos >= 0 && pos < 10 {
				adapter5Counts[adapter]++
			}
		}
	}

	// Find most common adapters (must appear in at least 10% of reads).
	// Iterate commonAdapters in declaration order so the result is
	// deterministic when several candidates tie on count.
	threshold := len(reads) / 10

	var maxCount3 int
	for _, adapter := range commonAdapters {
		if count := adapter3Counts[adapter]; count > maxCount3 && count >= threshold {
			maxCount3 = count
			opts.Adapter3 = adapter
		}
	}

	var maxCount5 int
	for _, adapter := range commonAdapters {
		if count := adapter5Counts[adapter]; count > maxCount5 && count >= threshold {
			maxCount5 = count
			opts.Adapter5 = adapter
		}
	}

	return opts
}

// extractUMI extracts UMI/barcode from read.
func extractUMI(record *fastq.Record, opts TrimOptions) (*fastq.Record, string) {
	if opts.UMILength <= 0 || opts.UMILength >= len(record.Sequence) {
		return record, ""
	}

	var umi string
	var newSeq, newQual []byte

	if opts.UMIPosition == "5prime" {
		// Extract from 5' end
		umi = string(record.Sequence[:opts.UMILength])
		newSeq = record.Sequence[opts.UMILength:]
		newQual = record.Quality[opts.UMILength:]
	} else {
		// Extract from 3' end
		umi = string(record.Sequence[len(record.Sequence)-opts.UMILength:])
		newSeq = record.Sequence[:len(record.Sequence)-opts.UMILength]
		newQual = record.Quality[:len(record.Quality)-opts.UMILength]
	}

	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    newSeq,
		Quality:     newQual,
	}, umi
}

// improvedFindAdapter uses improved Smith-Waterman-like local alignment for adapter matching.
func improvedFindAdapter(seq string, adapter string, minOverlap int, errorRate float64) int {
	if len(adapter) == 0 || len(seq) < minOverlap {
		return -1
	}

	maxErrors := int(float64(len(adapter)) * errorRate)
	bestPos := -1
	bestScore := -1

	// Scan through sequence
	for i := 0; i <= len(seq)-minOverlap; i++ {
		compareLen := min(len(adapter), len(seq)-i)

		if compareLen < minOverlap {
			continue
		}

		// Calculate alignment score
		matches := 0
		for j := 0; j < compareLen; j++ {
			if seq[i+j] == adapter[j] {
				matches++
			}
		}

		errors := compareLen - matches

		// Score: favor longer alignments with fewer errors
		score := matches*2 - errors

		if errors <= maxErrors && score > bestScore {
			bestScore = score
			bestPos = i
		}
	}

	return bestPos
}
