// Package fastp provides all-in-one preprocessing for FASTQ files.
// It combines quality filtering, adapter trimming, and various other preprocessing steps.
package fastp

import (
	"fmt"
	"io"
	"strings"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// Common adapter sequences for automatic detection
var CommonAdapters = map[string]string{
	"TruSeq":     "AGATCGGAAGAGC",
	"Nextera":    "CTGTCTCTTATA",
	"SmallRNA":   "TGGAATTCTCGG",
	"TruSeq_R2":  "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
	"Nextera_R2": "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
}

// ProcessOptions contains all preprocessing parameters.
type ProcessOptions struct {
	// Adapter trimming
	Adapter3      string
	Adapter5      string
	DetectAdapter bool

	// Quality filtering
	QualThreshold int
	MinLength     int
	MaxLength     int
	QualPercent   int // Percentage of bases that must meet quality threshold

	// Complexity filtering
	LowComplexity       bool
	ComplexityThreshold float64

	// Poly-tail trimming
	TrimPolyG   bool
	TrimPolyX   bool
	PolyGMinLen int

	// Sliding-window quality trimming
	CutFront       bool // Slide a window from the 5' end; trim leading bases while window mean quality < threshold
	CutTail        bool // Slide a window from the 3' end; trim trailing bases while window mean quality < threshold
	CutRight       bool // Slide a window 5'->3'; the moment a window's mean quality drops below threshold, cut there and to its right
	CutWindowSize  int  // Window size for cut_front/cut_tail/cut_right (default 4)
	CutMeanQuality int  // Mean-quality (Phred) threshold for the sliding window (default 20)

	// N filtering
	MaxNCount   int
	MaxNPercent float64

	// Length filtering
	LengthRequired int
	LengthLimit    int

	// UMI processing (legacy fields, kept for backward compatibility with
	// the older umi-length/umi-location/umi-skip flag names).
	UMILength   int
	UMILocation string // "read1", "read2", "index"

	// UMI processing — current API matching upstream fastp's --umi family.
	UMI       bool   // Enable UMI processing
	UMILoc    string // One of: read1, read2, per_read, index1, index2, per_index
	UMILen    int    // Number of UMI bases on the read (sequence-located modes only)
	UMIPrefix string // Optional prefix prepended to the UMI in the read name
	UMISkip   int    // Bases to skip immediately after the UMI bases

	// Duplication evaluation
	DupCalcAccuracy int  // 1..6; controls dup hash-table size. 0 means "disabled".
	Dedup           bool // When true, drop duplicate reads from the output stream.

	// Base correction
	BaseCorrection      bool
	CorrectionThreshold int

	// Overlap analysis (paired-end)
	MergeOverlap bool
	MinOverlap   int
	MaxMismatch  int

	// Multi-threading
	Threads int

	// HTML report
	HTMLReport string // Path to HTML report file

	// JSON report
	JSONReport string // Path to JSON report file

	// Automatic adapter detection
	DetectAdapterPE bool // Enable overlap-based adapter detection for paired-end
	DetectAdapterSE bool // Enable kmer-frequency-based adapter detection for single-end
}

// DefaultProcessOptions returns default processing options.
func DefaultProcessOptions() ProcessOptions {
	return ProcessOptions{
		Adapter3:            "",
		Adapter5:            "",
		DetectAdapter:       false,
		QualThreshold:       15,
		MinLength:           15,
		MaxLength:           0, // no limit
		QualPercent:         40,
		LowComplexity:       false,
		ComplexityThreshold: 0.3,
		TrimPolyG:           false,
		TrimPolyX:           false,
		PolyGMinLen:         10,
		CutFront:            false,
		CutTail:             false,
		CutRight:            false,
		CutWindowSize:       4,
		CutMeanQuality:      20,
		MaxNCount:           5,
		MaxNPercent:         20.0,
		LengthRequired:      15,
		LengthLimit:         0,
		UMILength:           0,
		UMILocation:         "",
		UMI:                 false,
		UMILoc:              "",
		UMILen:              0,
		UMIPrefix:           "",
		UMISkip:             0,
		DupCalcAccuracy:     0,
		Dedup:               false,
		BaseCorrection:      false,
		CorrectionThreshold: 20,
		MergeOverlap:        false,
		MinOverlap:          30,
		MaxMismatch:         5,
		Threads:             1,
		HTMLReport:          "",
		JSONReport:          "",
		DetectAdapterPE:     false,
		DetectAdapterSE:     false,
	}
}

// ProcessStats tracks preprocessing statistics.
//
// In addition to top-line counters used by the CLI, ProcessStats records
// per-cycle quality and base-composition histograms (separately for read 1
// and read 2) plus length distributions, which the HTML and JSON reports
// consume. The Mu field is used to guard updates from parallel workers.
type ProcessStats struct {
	// Mu guards all fields below when stats are accumulated from multiple
	// goroutines. Sequential callers may ignore it.
	Mu sync.Mutex

	TotalReads          int
	TotalBases          int64
	CleanReads          int
	CleanBases          int64
	LowQualityReads     int
	TooShortReads       int
	TooLongReads        int
	TooManyNReads       int
	AdapterTrimmedReads int
	AdapterTrimmedBases int64
	PolyGTrimmedReads   int
	PolyGTrimmedBases   int64
	QualityCutReads     int
	QualityCutBases     int64
	DetectedAdapter     string
	DetectedAdapterR1   string // Adapter detected for read 1 in PE mode
	DetectedAdapterR2   string // Adapter detected for read 2 in PE mode
	UMIExtracted        int    // Legacy counter (kept for back-compat with older tests/scripts).
	UMIProcessed        int    // Count of records that had a UMI extracted.
	BasesCorrected      int64
	OverlappingReads    int
	MergedReads         int

	// Duplication evaluation results, populated when DupCalcAccuracy > 0.
	DupRate      float64
	DupHist      map[int]int64
	DedupDropped int   // Reads removed from the output stream when Dedup is true.
	DupTotal     int64 // Total reads observed by the dup tracker.

	// Per-base quality + composition tracking, captured BEFORE filtering.
	// Index 0 is read 1; index 1 (when populated) is read 2.
	// QualSumByCycle[r][i] is the sum of phred qualities seen at cycle i;
	// QualCountByCycle[r][i] is the number of reads that reached cycle i.
	// BaseCountByCycle[r][b][i] is the count of base b at cycle i, where
	// b indexes into "ACGTN" via baseIndex.
	QualSumByCycle   [2][]int64
	QualCountByCycle [2][]int64
	BaseCountByCycle [2][5][]int64

	// Length histograms BEFORE and AFTER filtering, indexed by read (0/1).
	LengthHistBefore [2]map[int]int64
	LengthHistAfter  [2]map[int]int64

	// Aggregate quality buckets BEFORE filtering, totals across both reads.
	Q20BasesBefore int64
	Q30BasesBefore int64
	GCBasesBefore  int64

	// Aggregate quality buckets AFTER filtering.
	Q20BasesAfter int64
	Q30BasesAfter int64
	GCBasesAfter  int64
	TotalReadsR1  int
	TotalReadsR2  int
	TotalBasesR1  int64
	TotalBasesR2  int64
	CleanReadsR1  int
	CleanReadsR2  int
	CleanBasesR1  int64
	CleanBasesR2  int64
}

// baseIndex maps a base byte to its slot in BaseCountByCycle. Anything not
// A/C/G/T is treated as N.
func baseIndex(b byte) int {
	switch b {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	default:
		return 4
	}
}

// growCycles ensures the per-cycle slices in stats can hold at least n
// cycles for the given read index (0 or 1). It allocates lazily.
func (s *ProcessStats) growCycles(readIdx, n int) {
	if cap(s.QualSumByCycle[readIdx]) < n {
		nq := make([]int64, n)
		copy(nq, s.QualSumByCycle[readIdx])
		s.QualSumByCycle[readIdx] = nq
		nc := make([]int64, n)
		copy(nc, s.QualCountByCycle[readIdx])
		s.QualCountByCycle[readIdx] = nc
		for b := 0; b < 5; b++ {
			nb := make([]int64, n)
			copy(nb, s.BaseCountByCycle[readIdx][b])
			s.BaseCountByCycle[readIdx][b] = nb
		}
	} else if len(s.QualSumByCycle[readIdx]) < n {
		s.QualSumByCycle[readIdx] = s.QualSumByCycle[readIdx][:n]
		s.QualCountByCycle[readIdx] = s.QualCountByCycle[readIdx][:n]
		for b := 0; b < 5; b++ {
			s.BaseCountByCycle[readIdx][b] = s.BaseCountByCycle[readIdx][b][:n]
		}
	}
}

// recordBefore updates the BEFORE-filtering histograms for a single record.
// readIdx is 0 for R1 / SE, 1 for R2.
func (s *ProcessStats) recordBefore(record *fastq.Record, readIdx int, encoding fastq.QualityEncoding) {
	if record == nil {
		return
	}
	offset := phredOffset(encoding)
	n := len(record.Sequence)
	s.growCycles(readIdx, n)
	if s.LengthHistBefore[readIdx] == nil {
		s.LengthHistBefore[readIdx] = make(map[int]int64)
	}
	s.LengthHistBefore[readIdx][n]++
	if readIdx == 0 {
		s.TotalReadsR1++
		s.TotalBasesR1 += int64(n)
	} else {
		s.TotalReadsR2++
		s.TotalBasesR2 += int64(n)
	}
	for i := 0; i < n; i++ {
		b := record.Sequence[i]
		s.BaseCountByCycle[readIdx][baseIndex(b)][i]++
		if b == 'G' || b == 'g' || b == 'C' || b == 'c' {
			s.GCBasesBefore++
		}
		q := int(record.Quality[i]) - offset
		s.QualSumByCycle[readIdx][i] += int64(q)
		s.QualCountByCycle[readIdx][i]++
		if q >= 20 {
			s.Q20BasesBefore++
		}
		if q >= 30 {
			s.Q30BasesBefore++
		}
	}
}

// recordAfter updates the AFTER-filtering histograms for a single record
// that passed all filters.
func (s *ProcessStats) recordAfter(record *fastq.Record, readIdx int, encoding fastq.QualityEncoding) {
	if record == nil {
		return
	}
	offset := phredOffset(encoding)
	n := len(record.Sequence)
	if s.LengthHistAfter[readIdx] == nil {
		s.LengthHistAfter[readIdx] = make(map[int]int64)
	}
	s.LengthHistAfter[readIdx][n]++
	if readIdx == 0 {
		s.CleanReadsR1++
		s.CleanBasesR1 += int64(n)
	} else {
		s.CleanReadsR2++
		s.CleanBasesR2 += int64(n)
	}
	for i := 0; i < n; i++ {
		b := record.Sequence[i]
		if b == 'G' || b == 'g' || b == 'C' || b == 'c' {
			s.GCBasesAfter++
		}
		q := int(record.Quality[i]) - offset
		if q >= 20 {
			s.Q20BasesAfter++
		}
		if q >= 30 {
			s.Q30BasesAfter++
		}
	}
}

// OverlapResult represents the result of paired-end overlap analysis
type OverlapResult struct {
	HasOverlap    bool
	OverlapLength int
	Mismatches    int
	MergedSeq     string
	MergedQual    []byte
}

// ProcessPairedEnd processes paired-end FASTQ reads with all filters.
func ProcessPairedEnd(input1, input2 io.Reader, output1, output2 io.Writer, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	opts = normalizeUMIOptions(opts)
	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)

	stats := &ProcessStats{}
	var dupTracker *DupTracker
	if opts.DupCalcAccuracy > 0 || opts.Dedup {
		acc := opts.DupCalcAccuracy
		if acc <= 0 {
			acc = dupAccuracyDefault
		}
		dupTracker = NewDupTracker(acc)
	}

	// Process with multi-threading if enabled. UMI and duplication
	// tracking serialize through the input pipeline because they need a
	// deterministic per-record view, so we fall back to single-threaded
	// mode when either is active.
	if opts.Threads > 1 && !opts.UMI && opts.UMILength == 0 && dupTracker == nil {
		return processPairedEndParallel(reader1, reader2, writer1, writer2, encoding, opts, stats)
	}

	// Buffer the first batch of reads so we can run overlap-based PE
	// adapter detection on them before processing. The detected adapters
	// (if any) feed into Adapter3 / Adapter5 for the remainder of the run.
	type readPair struct {
		r1 *fastq.Record
		r2 *fastq.Record
	}
	var detectBuffer []readPair
	if opts.DetectAdapterPE {
		detectBuffer = make([]readPair, 0, adapterDetectSampleSize)
		for i := 0; i < adapterDetectSampleSize; i++ {
			r1, err1 := reader1.Read()
			r2, err2 := reader2.Read()
			if err1 == io.EOF && err2 == io.EOF {
				break
			}
			if err1 == io.EOF || err2 == io.EOF {
				return stats, fmt.Errorf("paired files have different number of reads")
			}
			if err1 != nil {
				return stats, fmt.Errorf("error reading read1: %w", err1)
			}
			if err2 != nil {
				return stats, fmt.Errorf("error reading read2: %w", err2)
			}
			detectBuffer = append(detectBuffer, readPair{r1: r1, r2: r2})
		}
		pairs := make([][2]*fastq.Record, len(detectBuffer))
		for i, p := range detectBuffer {
			pairs[i] = [2]*fastq.Record{p.r1, p.r2}
		}
		r1Adapter, r2Adapter := DetectAdaptersFromPairs(pairs)
		stats.DetectedAdapterR1 = r1Adapter
		stats.DetectedAdapterR2 = r2Adapter
		if r1Adapter != "" {
			stats.DetectedAdapter = r1Adapter
		}
		if opts.Adapter3 == "" && r1Adapter != "" {
			opts.Adapter3 = r1Adapter
		}
		if opts.Adapter5 == "" && r2Adapter != "" {
			opts.Adapter5 = r2Adapter
		}
	}

	processPair := func(record1, record2 *fastq.Record) error {
		stats.TotalReads += 2
		stats.TotalBases += int64(len(record1.Sequence) + len(record2.Sequence))
		stats.recordBefore(record1, 0, encoding)
		stats.recordBefore(record2, 1, encoding)

		// Duplication evaluation: hash the R1 sequence (matches upstream
		// fastp). When --dedup is on and the read is a duplicate, drop
		// both mates and skip further processing.
		if dupTracker != nil {
			isDup := dupTracker.Observe(record1.Sequence)
			if isDup && opts.Dedup {
				stats.DedupDropped += 2
				return nil
			}
		}

		// Extract UMI if configured
		if opts.UMI || opts.UMILength > 0 {
			record1, record2 = applyUMI(record1, record2, opts, stats)
		}

		// Check for overlap and merge if enabled
		if opts.MergeOverlap {
			overlap := analyzeOverlap(record1, record2, opts, encoding)
			if overlap.HasOverlap {
				stats.OverlappingReads++
				if overlap.OverlapLength >= opts.MinOverlap && overlap.Mismatches <= opts.MaxMismatch {
					merged := &fastq.Record{
						ID:          record1.ID,
						Description: record1.Description,
						Sequence:    []byte(overlap.MergedSeq),
						Quality:     overlap.MergedQual,
					}
					processed, pass := processRecord(merged, opts, stats, encoding)
					if pass {
						if err := writer1.Write(processed); err != nil {
							return fmt.Errorf("error writing merged read: %w", err)
						}
						stats.CleanReads++
						stats.CleanBases += int64(len(processed.Sequence))
						stats.MergedReads++
						stats.recordAfter(processed, 0, encoding)
					}
					return nil
				}
			}
		}

		processed1, pass1 := processRecord(record1, opts, stats, encoding)
		processed2, pass2 := processRecord(record2, opts, stats, encoding)

		if pass1 && pass2 {
			if err := writer1.Write(processed1); err != nil {
				return fmt.Errorf("error writing read1: %w", err)
			}
			if err := writer2.Write(processed2); err != nil {
				return fmt.Errorf("error writing read2: %w", err)
			}
			stats.CleanReads += 2
			stats.CleanBases += int64(len(processed1.Sequence) + len(processed2.Sequence))
			stats.recordAfter(processed1, 0, encoding)
			stats.recordAfter(processed2, 1, encoding)
		}
		return nil
	}

	for _, p := range detectBuffer {
		if err := processPair(p.r1, p.r2); err != nil {
			return stats, err
		}
	}

	for {
		record1, err1 := reader1.Read()
		record2, err2 := reader2.Read()
		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF || err2 == io.EOF {
			return stats, fmt.Errorf("paired files have different number of reads")
		}
		if err1 != nil {
			return stats, fmt.Errorf("error reading read1: %w", err1)
		}
		if err2 != nil {
			return stats, fmt.Errorf("error reading read2: %w", err2)
		}
		if err := processPair(record1, record2); err != nil {
			return stats, err
		}
	}

	// Flush writers and return early; the old loop body is removed below.
	if err := writer1.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output1: %w", err)
	}
	if err := writer2.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output2: %w", err)
	}
	finalizeDupStats(stats, dupTracker)
	return stats, nil
}

// finalizeDupStats copies the final duplication metrics out of the
// tracker into stats. Called once at the end of each Process* function.
// It is safe to call with a nil tracker (no-op).
func finalizeDupStats(stats *ProcessStats, tracker *DupTracker) {
	if stats == nil || tracker == nil {
		return
	}
	stats.DupRate = tracker.Rate()
	stats.DupHist = tracker.Histogram()
	stats.DupTotal = tracker.Total()
}

// ProcessSingleEnd processes single-end FASTQ reads with all filters.
func ProcessSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	opts = normalizeUMIOptions(opts)
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	stats := &ProcessStats{}
	var dupTracker *DupTracker
	if opts.DupCalcAccuracy > 0 || opts.Dedup {
		acc := opts.DupCalcAccuracy
		if acc <= 0 {
			acc = dupAccuracyDefault
		}
		dupTracker = NewDupTracker(acc)
	}

	// Process with multi-threading if enabled. UMI and duplication
	// tracking serialize through the input pipeline because they need a
	// deterministic per-record view, so we fall back to single-threaded
	// mode when either is active.
	if opts.Threads > 1 && !opts.UMI && opts.UMILength == 0 && dupTracker == nil {
		return processSingleEndParallel(reader, writer, encoding, opts, stats)
	}

	// Buffer reads when SE adapter detection is requested so the same
	// reads can be inspected before processing. We only buffer up to
	// adapterDetectSampleSize records; remaining reads stream normally.
	var detectBuffer []*fastq.Record
	if opts.DetectAdapterSE && opts.Adapter3 == "" {
		detectBuffer = make([]*fastq.Record, 0, adapterDetectSampleSize)
		for i := 0; i < adapterDetectSampleSize; i++ {
			record, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return stats, fmt.Errorf("error reading FASTQ: %w", err)
			}
			detectBuffer = append(detectBuffer, record)
		}
		stats.DetectedAdapter = DetectAdapterSE(detectBuffer)
		stats.DetectedAdapterR1 = stats.DetectedAdapter
		if stats.DetectedAdapter != "" && opts.Adapter3 == "" {
			opts.Adapter3 = stats.DetectedAdapter
		}
	}

	processOne := func(record *fastq.Record) error {
		stats.TotalReads++
		originalLength := len(record.Sequence)
		stats.TotalBases += int64(originalLength)
		stats.recordBefore(record, 0, encoding)

		// Duplication tracking: hash the read sequence and optionally
		// drop duplicates when --dedup is on.
		if dupTracker != nil {
			isDup := dupTracker.Observe(record.Sequence)
			if isDup && opts.Dedup {
				stats.DedupDropped++
				return nil
			}
		}

		// Extract UMI if configured
		if opts.UMI || opts.UMILength > 0 {
			record, _ = applyUMI(record, nil, opts, stats)
		}

		processed, pass := processRecord(record, opts, stats, encoding)
		if pass {
			if err := writer.Write(processed); err != nil {
				return fmt.Errorf("error writing FASTQ: %w", err)
			}
			stats.CleanReads++
			stats.CleanBases += int64(len(processed.Sequence))
			stats.recordAfter(processed, 0, encoding)
		}
		return nil
	}

	for _, rec := range detectBuffer {
		if err := processOne(rec); err != nil {
			return stats, err
		}
	}

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return stats, fmt.Errorf("error reading FASTQ: %w", err)
		}
		if err := processOne(record); err != nil {
			return stats, err
		}
	}

	// Flush writer
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}

	finalizeDupStats(stats, dupTracker)
	return stats, nil
}

// processRecord applies all processing steps to a single record.
func processRecord(record *fastq.Record, opts ProcessOptions, stats *ProcessStats, encoding fastq.QualityEncoding) (*fastq.Record, bool) {
	seq := string(record.Sequence)
	qual := record.Quality

	start := 0
	end := len(seq)

	// Step 0: Base correction if enabled
	if opts.BaseCorrection {
		seq, qual = correctBases(seq, qual, opts.CorrectionThreshold, stats, encoding)
		record.Sequence = []byte(seq)
		record.Quality = qual
	}

	// Step 1: Trim adapters if specified
	if opts.Adapter5 != "" {
		pos := findAdapter(seq, opts.Adapter5)
		if pos >= 0 {
			start = pos + len(opts.Adapter5)
			stats.AdapterTrimmedReads++
			stats.AdapterTrimmedBases += int64(pos + len(opts.Adapter5))
		}
	}

	if opts.Adapter3 != "" {
		pos := findAdapter(seq[start:], opts.Adapter3)
		if pos >= 0 {
			end = start + pos
			stats.AdapterTrimmedReads++
			stats.AdapterTrimmedBases += int64(len(seq) - end)
		}
	}

	// Step 2: Trim poly-G tails if enabled
	if opts.TrimPolyG {
		polyLen := countPolyTail(seq[start:end], 'G')
		if polyLen >= opts.PolyGMinLen {
			end -= polyLen
			stats.PolyGTrimmedReads++
			stats.PolyGTrimmedBases += int64(polyLen)
		}
	}

	// Step 3: Trim poly-X tails if enabled
	if opts.TrimPolyX {
		for _, base := range []byte{'A', 'T', 'C'} {
			polyLen := countPolyTail(seq[start:end], rune(base))
			if polyLen >= opts.PolyGMinLen {
				end -= polyLen
			}
		}
	}

	// Step 3.5: Sliding-window quality trimming (cut_front / cut_tail / cut_right)
	if opts.CutFront || opts.CutTail || opts.CutRight {
		cutLo, cutHi := slidingWindowCut(qual[start:end], encoding, opts)
		if cutLo != 0 || cutHi != end-start {
			removed := (cutLo) + (end - start - cutHi)
			start += cutLo
			end = start + (cutHi - cutLo)
			if removed > 0 {
				stats.QualityCutReads++
				stats.QualityCutBases += int64(removed)
			}
		}
	}

	// Step 4: Quality-based trimming
	if opts.QualThreshold > 0 {
		start, end = trimByQuality(qual[start:end], opts.QualThreshold, start, end, encoding)
	}

	// Check if read is too short after trimming
	if end-start < opts.MinLength || end-start < opts.LengthRequired {
		stats.TooShortReads++
		return nil, false
	}

	// Check if read is too long
	if opts.MaxLength > 0 && end-start > opts.MaxLength {
		stats.TooLongReads++
		return nil, false
	}

	if opts.LengthLimit > 0 && end-start > opts.LengthLimit {
		end = start + opts.LengthLimit
	}

	// Step 5: Check N content
	nCount := countNs(seq[start:end])
	nPercent := 100.0 * float64(nCount) / float64(end-start)

	if nCount > opts.MaxNCount || nPercent > opts.MaxNPercent {
		stats.TooManyNReads++
		return nil, false
	}

	// Step 6: Check quality (percentage of bases meeting threshold)
	if opts.QualPercent > 0 {
		qualScores := getQualityScores(qual[start:end], encoding)
		passCount := 0
		for _, q := range qualScores {
			if q >= opts.QualThreshold {
				passCount++
			}
		}
		passPercent := 100.0 * float64(passCount) / float64(len(qualScores))
		if passPercent < float64(opts.QualPercent) {
			stats.LowQualityReads++
			return nil, false
		}
	}

	// Step 7: Check complexity if enabled
	if opts.LowComplexity {
		complexity := calculateComplexity(seq[start:end])
		if complexity < opts.ComplexityThreshold {
			// Low complexity read - discard
			return nil, false
		}
	}

	// Create processed record
	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence[start:end],
		Quality:     record.Quality[start:end],
	}, true
}

// findAdapter finds the position of an adapter in a sequence.
func findAdapter(seq string, adapter string) int {
	return strings.Index(seq, adapter)
}

// countPolyTail counts the length of a poly-X tail at the end of a sequence.
func countPolyTail(seq string, base rune) int {
	count := 0
	for i := len(seq) - 1; i >= 0 && rune(seq[i]) == base; i-- {
		count++
	}
	return count
}

// countNs counts the number of N bases in a sequence.
func countNs(seq string) int {
	count := 0
	for _, b := range seq {
		if b == 'N' || b == 'n' {
			count++
		}
	}
	return count
}

// trimByQuality trims low-quality regions from both ends.
func trimByQuality(quality []byte, threshold int, start, end int, encoding fastq.QualityEncoding) (int, int) {
	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}

	// Trim from 3' end
	for end > start && int(quality[end-start-1])-offset < threshold {
		end--
	}

	// Trim from 5' end
	for start < end && int(quality[0])-offset < threshold {
		start++
		quality = quality[1:]
	}

	return start, end
}

// phredOffset returns the ASCII offset for the given quality encoding.
func phredOffset(encoding fastq.QualityEncoding) int {
	if encoding == fastq.Phred64 {
		return 64
	}
	return 33
}

// slidingWindowCut applies fastp's sliding-window quality trimming to a quality
// slice and returns the half-open index range [lo, hi) of the bases that should
// be kept. The three modes (cut_front, cut_tail, cut_right) are applied in that
// order when more than one is enabled, matching upstream fastp.
//
// cut_front: scanning 5'->3', find the first window whose mean quality is >=
// CutMeanQuality and trim everything before it.
// cut_tail: scanning 3'->5', find the first window whose mean quality is >=
// CutMeanQuality and trim everything after it.
// cut_right: scanning 5'->3', the moment a window's mean quality drops below
// CutMeanQuality, cut the read at that window's start (keep [0, windowStart)).
//
// If the window size is larger than the (current) read, the whole read is
// treated as a single short window.
func slidingWindowCut(quality []byte, encoding fastq.QualityEncoding, opts ProcessOptions) (lo, hi int) {
	offset := phredOffset(encoding)
	window := opts.CutWindowSize
	if window < 1 {
		window = 1
	}
	threshold := float64(opts.CutMeanQuality)

	lo, hi = 0, len(quality)

	// meanQual computes the mean quality of quality[a:b].
	meanQual := func(a, b int) float64 {
		if b <= a {
			return 0
		}
		sum := 0
		for i := a; i < b; i++ {
			sum += int(quality[i]) - offset
		}
		return float64(sum) / float64(b-a)
	}

	if opts.CutFront {
		n := hi - lo
		w := window
		if w > n {
			w = n
		}
		if w > 0 {
			newLo := lo
			for newLo+w <= hi {
				if meanQual(newLo, newLo+w) >= threshold {
					break
				}
				newLo++
			}
			if newLo+w > hi {
				// No qualifying window found: drop everything.
				newLo = hi
			}
			lo = newLo
		}
	}

	if opts.CutTail {
		n := hi - lo
		w := window
		if w > n {
			w = n
		}
		if w > 0 {
			newHi := hi
			for newHi-w >= lo {
				if meanQual(newHi-w, newHi) >= threshold {
					break
				}
				newHi--
			}
			if newHi-w < lo {
				// No qualifying window found: drop everything.
				newHi = lo
			}
			hi = newHi
		}
	}

	if opts.CutRight {
		n := hi - lo
		w := window
		if w > n {
			w = n
		}
		if w > 0 {
			for s := lo; s+w <= hi; s++ {
				if meanQual(s, s+w) < threshold {
					hi = s
					break
				}
			}
		}
	}

	if hi < lo {
		hi = lo
	}
	return lo, hi
}

// getQualityScores converts ASCII-encoded quality scores to numeric values.
func getQualityScores(quality []byte, encoding fastq.QualityEncoding) []int {
	scores := make([]int, len(quality))
	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}
	for i, q := range quality {
		scores[i] = int(q) - offset
	}
	return scores
}

// calculateComplexity calculates sequence complexity (0-1, higher is more complex).
func calculateComplexity(seq string) float64 {
	if len(seq) == 0 {
		return 0
	}

	// Simple complexity measure: unique 2-mers / total 2-mers
	if len(seq) < 2 {
		return 1.0
	}

	kmers := make(map[string]bool)
	for i := 0; i < len(seq)-1; i++ {
		kmer := seq[i : i+2]
		kmers[kmer] = true
	}

	return float64(len(kmers)) / float64(len(seq)-1)
}

// detectAdapter detects the most common adapter in a sample of reads
func detectAdapter(input1, input2 io.Reader, encoding fastq.QualityEncoding) string {
	// Sample first 10000 reads to detect adapter
	reader1 := fastq.NewReader(input1, encoding)

	adapterCounts := make(map[string]int)
	sampleSize := 10000

	for i := 0; i < sampleSize; i++ {
		record, err := reader1.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		seq := string(record.Sequence)
		// Check each common adapter
		for name, adapter := range CommonAdapters {
			if strings.Contains(seq, adapter) || strings.Contains(seq, adapter[:len(adapter)/2]) {
				adapterCounts[name]++
			}
		}
	}

	// Find most common adapter
	maxCount := 0
	detectedName := ""
	for name, count := range adapterCounts {
		if count > maxCount {
			maxCount = count
			detectedName = name
		}
	}

	if detectedName != "" && maxCount > sampleSize/100 { // At least 1% detection rate
		return CommonAdapters[detectedName]
	}

	return ""
}

// DetectAdapterFromReads detects the most common adapter from FASTQ reads
// This is a public function that can be called before processing
func DetectAdapterFromReads(reads []string) string {
	adapterCounts := make(map[string]int)

	for _, seq := range reads {
		// Check each common adapter
		for name, adapter := range CommonAdapters {
			if strings.Contains(seq, adapter) || strings.Contains(seq, adapter[:len(adapter)/2]) {
				adapterCounts[name]++
			}
		}
	}

	// Find most common adapter
	maxCount := 0
	detectedName := ""
	for name, count := range adapterCounts {
		if count > maxCount {
			maxCount = count
			detectedName = name
		}
	}

	if detectedName != "" && maxCount > len(reads)/100 { // At least 1% detection rate
		return CommonAdapters[detectedName]
	}

	return ""
}

// normalizeUMIOptions promotes the legacy --umi-length / --umi-location /
// --umi-skip fields into the newer fastp-aligned UMI/UMILen/UMILoc/UMISkip
// fields when the caller has only filled in the legacy set. This keeps
// older programmatic callers working without changing the public API.
func normalizeUMIOptions(opts ProcessOptions) ProcessOptions {
	if !opts.UMI && opts.UMILength > 0 {
		opts.UMI = true
		if opts.UMILen == 0 {
			opts.UMILen = opts.UMILength
		}
		if opts.UMILoc == "" {
			if opts.UMILocation != "" {
				opts.UMILoc = opts.UMILocation
			} else {
				opts.UMILoc = UMILocRead1
			}
		}
	}
	return opts
}

// extractUMI is a thin shim retained for the existing call sites that
// invoke UMI extraction with two records. It delegates to applyUMI so all
// supported --umi_loc modes are honoured. Legacy callers that only set
// the UMILength/UMILocation/UMISkip fields are accommodated by
// normalizeUMIOptions.
func extractUMI(record1, record2 *fastq.Record, opts ProcessOptions, stats *ProcessStats) (*fastq.Record, *fastq.Record) {
	return applyUMI(record1, record2, normalizeUMIOptions(opts), stats)
}

// correctBases performs base correction based on quality scores
func correctBases(seq string, qual []byte, threshold int, stats *ProcessStats, encoding fastq.QualityEncoding) (string, []byte) {
	if threshold == 0 {
		return seq, qual
	}

	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}

	corrected := []byte(seq)
	correctedCount := int64(0)

	for i := 0; i < len(seq); i++ {
		q := int(qual[i]) - offset

		// Correct low-quality bases to N
		if q < threshold && seq[i] != 'N' {
			corrected[i] = 'N'
			correctedCount++
		}
	}

	stats.BasesCorrected += correctedCount
	return string(corrected), qual
}

// analyzeOverlap analyzes overlap between paired-end reads
func analyzeOverlap(record1, record2 *fastq.Record, opts ProcessOptions, encoding fastq.QualityEncoding) OverlapResult {
	result := OverlapResult{HasOverlap: false}

	seq1 := string(record1.Sequence)
	seq2 := reverseComplement(string(record2.Sequence))
	qual1 := record1.Quality
	qual2 := reverseSlice(record2.Quality)

	minLen := opts.MinOverlap
	if minLen < 10 {
		minLen = 10
	}

	// Try different overlap positions
	bestOverlap := 0
	bestMismatches := len(seq1)
	bestPos := -1

	for pos := len(seq1) - minLen; pos >= 0; pos-- {
		overlapLen := len(seq1) - pos
		if overlapLen > len(seq2) {
			overlapLen = len(seq2)
		}

		mismatches := 0
		for i := 0; i < overlapLen; i++ {
			if seq1[pos+i] != seq2[i] {
				mismatches++
			}
		}

		if mismatches < bestMismatches && overlapLen >= minLen {
			bestOverlap = overlapLen
			bestMismatches = mismatches
			bestPos = pos
		}
	}

	if bestPos >= 0 && bestOverlap >= minLen {
		result.HasOverlap = true
		result.OverlapLength = bestOverlap
		result.Mismatches = bestMismatches

		// Merge sequences using quality scores
		merged := make([]byte, bestPos+len(seq2))
		mergedQual := make([]byte, bestPos+len(seq2))

		// Copy non-overlapping part of seq1
		copy(merged[:bestPos], seq1[:bestPos])
		copy(mergedQual[:bestPos], qual1[:bestPos])

		// Merge overlapping region
		offset := 33
		if encoding == fastq.Phred64 {
			offset = 64
		}

		for i := 0; i < len(seq2); i++ {
			pos := bestPos + i
			if i < bestOverlap {
				// In overlap - use higher quality base
				q1 := int(qual1[pos]) - offset
				q2 := int(qual2[i]) - offset

				if q1 >= q2 {
					merged[pos] = seq1[pos]
					mergedQual[pos] = qual1[pos]
				} else {
					merged[pos] = seq2[i]
					mergedQual[pos] = qual2[i]
				}
			} else {
				// After overlap - use seq2
				merged[pos] = seq2[i]
				mergedQual[pos] = qual2[i]
			}
		}

		result.MergedSeq = string(merged)
		result.MergedQual = mergedQual
	}

	return result
}

// reverseComplement returns the reverse complement of a DNA sequence
func reverseComplement(seq string) string {
	complement := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
		'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
		'N': 'N', 'n': 'n',
	}

	result := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		if comp, ok := complement[seq[len(seq)-1-i]]; ok {
			result[i] = comp
		} else {
			result[i] = 'N'
		}
	}

	return string(result)
}

// reverseSlice reverses a byte slice
func reverseSlice(s []byte) []byte {
	result := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		result[i] = s[len(s)-1-i]
	}
	return result
}

// processPairedEndParallel processes paired-end reads in parallel
func processPairedEndParallel(reader1, reader2 *fastq.Reader, writer1, writer2 *fastq.Writer, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats) (*ProcessStats, error) {
	type readPair struct {
		record1 *fastq.Record
		record2 *fastq.Record
	}

	type resultPair struct {
		processed1 *fastq.Record
		processed2 *fastq.Record
		pass       bool
	}

	inputChan := make(chan readPair, opts.Threads*2)
	outputChan := make(chan resultPair, opts.Threads*2)

	var wg sync.WaitGroup
	var statsMu sync.Mutex

	// Start worker goroutines
	for i := 0; i < opts.Threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			localStats := &ProcessStats{}

			for pair := range inputChan {
				processed1, pass1 := processRecord(pair.record1, opts, localStats, encoding)
				processed2, pass2 := processRecord(pair.record2, opts, localStats, encoding)

				outputChan <- resultPair{
					processed1: processed1,
					processed2: processed2,
					pass:       pass1 && pass2,
				}
			}

			statsMu.Lock()
			stats.LowQualityReads += localStats.LowQualityReads
			stats.TooShortReads += localStats.TooShortReads
			stats.TooLongReads += localStats.TooLongReads
			stats.TooManyNReads += localStats.TooManyNReads
			statsMu.Unlock()
		}()
	}

	// Reader goroutine
	go func() {
		for {
			record1, err1 := reader1.Read()
			record2, err2 := reader2.Read()

			if err1 == io.EOF && err2 == io.EOF {
				break
			}
			if err1 != nil || err2 != nil {
				break
			}

			stats.TotalReads += 2
			stats.TotalBases += int64(len(record1.Sequence) + len(record2.Sequence))
			stats.recordBefore(record1, 0, encoding)
			stats.recordBefore(record2, 1, encoding)

			inputChan <- readPair{record1: record1, record2: record2}
		}
		close(inputChan)
	}()

	// Wait for workers to finish
	go func() {
		wg.Wait()
		close(outputChan)
	}()

	// Writer goroutine
	for result := range outputChan {
		if result.pass {
			if err := writer1.Write(result.processed1); err != nil {
				return stats, fmt.Errorf("error writing read1: %w", err)
			}
			if err := writer2.Write(result.processed2); err != nil {
				return stats, fmt.Errorf("error writing read2: %w", err)
			}
			stats.CleanReads += 2
			stats.CleanBases += int64(len(result.processed1.Sequence) + len(result.processed2.Sequence))
			stats.recordAfter(result.processed1, 0, encoding)
			stats.recordAfter(result.processed2, 1, encoding)
		}
	}

	return stats, nil
}

// processSingleEndParallel processes single-end reads in parallel
func processSingleEndParallel(reader *fastq.Reader, writer *fastq.Writer, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats) (*ProcessStats, error) {
	type result struct {
		processed *fastq.Record
		pass      bool
	}

	inputChan := make(chan *fastq.Record, opts.Threads*2)
	outputChan := make(chan result, opts.Threads*2)

	var wg sync.WaitGroup
	var statsMu sync.Mutex

	// Start worker goroutines
	for i := 0; i < opts.Threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			localStats := &ProcessStats{}

			for record := range inputChan {
				processed, pass := processRecord(record, opts, localStats, encoding)
				outputChan <- result{processed: processed, pass: pass}
			}

			statsMu.Lock()
			stats.LowQualityReads += localStats.LowQualityReads
			stats.TooShortReads += localStats.TooShortReads
			stats.TooLongReads += localStats.TooLongReads
			stats.TooManyNReads += localStats.TooManyNReads
			statsMu.Unlock()
		}()
	}

	// Reader goroutine
	go func() {
		for {
			record, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				break
			}

			stats.TotalReads++
			stats.TotalBases += int64(len(record.Sequence))
			stats.recordBefore(record, 0, encoding)
			inputChan <- record
		}
		close(inputChan)
	}()

	// Wait for workers to finish
	go func() {
		wg.Wait()
		close(outputChan)
	}()

	// Writer goroutine
	for res := range outputChan {
		if res.pass {
			if err := writer.Write(res.processed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			stats.CleanReads++
			stats.CleanBases += int64(len(res.processed.Sequence))
			stats.recordAfter(res.processed, 0, encoding)
		}
	}

	return stats, nil
}
