// Package fastp provides all-in-one preprocessing for FASTQ files.
// It combines quality filtering, adapter trimming, and various other preprocessing steps.
package fastp

import (
	"fmt"
	"io"
	"strings"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
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

	// DisableAdapterTrimming maps to upstream --disable_adapter_trimming
	// (-A). When set, all adapter trimming (explicit Adapter3/Adapter5,
	// auto-detection, and AdapterFasta) is skipped.
	DisableAdapterTrimming bool

	// DisableQualityFiltering maps to upstream --disable_quality_filtering
	// (-Q). Upstream gates the whole quality-filter block on
	// qualfilter.enabled (filter.cpp:43-50): the low-quality-base
	// percentage limit, the average-quality requirement, AND the
	// N-base-count limit. When this is set the Go port skips all three so
	// every read survives the quality filter, matching upstream.
	DisableQualityFiltering bool

	// DisableLengthFiltering maps to upstream --disable_length_filtering
	// (-L). Upstream gates the length-filter block on lengthFilter.enabled
	// (filter.cpp:52-57): the length_required (too-short) check and the
	// length_limit (too-long) check. When this is set the Go port skips
	// both, matching upstream.
	DisableLengthFiltering bool

	// DisableTrimPolyG maps to upstream --disable_trim_poly_g (-G). When
	// set it force-disables poly-G tail trimming even if TrimPolyG was
	// requested (upstream errors if both --trim_poly_g and
	// --disable_trim_poly_g are given; here it simply wins, mirroring the
	// "disabled" outcome).
	DisableTrimPolyG bool

	// AdapterFasta holds adapter sequences loaded from --adapter_fasta.
	// Every read (R1 and, in PE mode, R2) is trimmed against each of these
	// sequences after any explicit/auto-detected adapter pass. Ordered by
	// sorted FASTA contig name to mirror upstream's std::map iteration.
	AdapterFasta []string

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
	PolyXMinLen int // Minimum poly-X length (upstream --poly_x_min_len, default 10).

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

	// Overlap-based base correction (paired-end), upstream --correction.
	// When enabled, mismatched bases inside the detected PE overlap are
	// corrected using the higher-quality mate (see overlap.go).
	Correction bool

	// Overlap analysis knobs shared by --correction and merge. These map
	// directly to upstream's --overlap_len_require / --overlap_diff_limit /
	// --overlap_diff_percent_limit.
	OverlapRequire          int // Minimum overlap length (default 30).
	OverlapDiffLimit        int // Maximum mismatched bases in overlap (default 5).
	OverlapDiffPercentLimit int // Maximum mismatch percentage in overlap (default 20).
	InsertSizeMax           int // Insert-size histogram cap / unknown bucket (upstream --insert_size_max, default 512).

	// Overlap analysis (paired-end)
	MergeOverlap bool
	MinOverlap   int
	MaxMismatch  int

	// Overlap-driven merge writer, upstream -m/--merge family. When Merge
	// is set, overlapping pairs are merged into a single read (via the
	// upstream OverlapAnalysis::merge port) and written to the merge
	// stream; non-overlapping pairs are dropped unless IncludeUnmerged is
	// set, in which case each surviving mate is written to the merge stream.
	Merge           bool
	IncludeUnmerged bool

	// Overrepresentation analysis, upstream -p / -P.
	OverrepAnalysis bool // Enable overrepresented-sequence analysis.
	OverrepSampling int  // 1-in-N sampling rate (default 20).

	// Output splitting, upstream -s / -S / -d.
	SplitNumber       int // --split: split into this many files (2-999).
	SplitByLines      int // --split_by_lines: max lines per output file.
	SplitPrefixDigits int // --split_prefix_digits: zero-pad width (default 4).

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
		Adapter3:                "",
		Adapter5:                "",
		DetectAdapter:           false,
		QualThreshold:           15,
		MinLength:               15,
		MaxLength:               0, // no limit
		QualPercent:             40,
		LowComplexity:           false,
		ComplexityThreshold:     0.3,
		TrimPolyG:               false,
		TrimPolyX:               false,
		PolyGMinLen:             10,
		PolyXMinLen:             10,
		CutFront:                false,
		CutTail:                 false,
		CutRight:                false,
		CutWindowSize:           4,
		CutMeanQuality:          20,
		MaxNCount:               5,
		MaxNPercent:             20.0,
		LengthRequired:          15,
		LengthLimit:             0,
		UMILength:               0,
		UMILocation:             "",
		UMI:                     false,
		UMILoc:                  "",
		UMILen:                  0,
		UMIPrefix:               "",
		UMISkip:                 0,
		DupCalcAccuracy:         0,
		Dedup:                   false,
		BaseCorrection:          false,
		CorrectionThreshold:     20,
		Correction:              false,
		OverlapRequire:          30,
		OverlapDiffLimit:        5,
		OverlapDiffPercentLimit: 20,
		InsertSizeMax:           defaultInsertSizeMax,
		MergeOverlap:            false,
		MinOverlap:              30,
		MaxMismatch:             5,
		Merge:                   false,
		IncludeUnmerged:         false,
		OverrepAnalysis:         false,
		OverrepSampling:         20,
		SplitNumber:             0,
		SplitByLines:            0,
		SplitPrefixDigits:       4,
		Threads:                 1,
		HTMLReport:              "",
		JSONReport:              "",
		DetectAdapterPE:         false,
		DetectAdapterSE:         false,
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
	LowComplexityReads  int
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
	CorrectedReads      int64 // Reads with >=1 base corrected by overlap analysis (--correction).
	OverlappingReads    int
	MergedReads         int
	// MergeEnabled records whether -m/--merge was active. Upstream renames the
	// JSON `read1_after_filtering` block to `merged_and_filtered` (and omits
	// `read2_after_filtering`) when the merge flag is set, regardless of how
	// many pairs actually merged (jsonreporter.cpp:158).
	MergeEnabled bool

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

	// Base-resolved per-cycle curves + 5-mer histograms, indexed by read
	// (0/1), for the BEFORE and AFTER streams. These drive the JSON report's
	// quality_curves / content_curves / kmer_count / q40_bases sub-fields and
	// the per-read q20/q30 totals (reproducing upstream Stats; see
	// report_curves.go). Allocated lazily by recordBefore/recordAfter.
	curvesBefore [2]*readCurves
	curvesAfter  [2]*readCurves

	// Per-read aggregate quality buckets, split by stream and read index, so
	// the per-read JSON blocks report the real q20/q30 totals upstream emits
	// (rather than the previous 0 placeholders). Index [0]=before, [1]=after.
	Q20ByRead [2][2]int64
	Q30ByRead [2][2]int64

	// Insert-size histogram for paired-end overlap analysis, reproducing
	// upstream's mInsertHist (peprocessor.cpp). Index d in [0, insertSizeMax)
	// counts pairs whose inferred insert size is d; index insertSizeMax is the
	// "unknown" bucket (no detectable overlap). Populated only in PE mode.
	InsertHist []int64

	// Overrepresented-sequence analyzers, indexed by read (0/1), populated
	// when OverrepAnalysis is enabled. Tracks the before-filtering stream.
	overrep [2]*overrepAnalyzer

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
			s.Q20ByRead[0][readIdx]++
		}
		if q >= 30 {
			s.Q30BasesBefore++
			s.Q30ByRead[0][readIdx]++
		}
	}
	if s.curvesBefore[readIdx] == nil {
		s.curvesBefore[readIdx] = &readCurves{}
	}
	s.curvesBefore[readIdx].stat(record, offset)
	if s.overrep[readIdx] != nil {
		s.overrep[readIdx].sampleRead(string(record.Sequence))
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
			s.Q20ByRead[1][readIdx]++
		}
		if q >= 30 {
			s.Q30BasesAfter++
			s.Q30ByRead[1][readIdx]++
		}
	}
	if s.curvesAfter[readIdx] == nil {
		s.curvesAfter[readIdx] = &readCurves{}
	}
	s.curvesAfter[readIdx].stat(record, offset)
}

// defaultInsertSizeMax is upstream fastp's default --insert_size_max (512):
// the histogram length and the cap/"unknown" bucket index.
const defaultInsertSizeMax = 512

// statInsertSize reproduces upstream PairEndProcessor::statInsertSize
// (peprocessor.cpp:698-711): it runs the overlap analysis on the pair and
// bins the inferred insert size into InsertHist. frontTrimmed1/2 account for
// bases removed by UMI/front trimming before this point (0 in the common
// case). A pair with no detectable overlap is binned at insertSizeMax (the
// "unknown" bucket); any computed size above insertSizeMax is clamped there.
func (s *ProcessStats) statInsertSize(record1, record2 *fastq.Record, opts ProcessOptions, frontTrimmed1, frontTrimmed2 int) {
	if record1 == nil || record2 == nil {
		return
	}
	max := opts.InsertSizeMax
	if max <= 0 {
		max = defaultInsertSizeMax
	}
	if s.InsertHist == nil {
		s.InsertHist = make([]int64, max+1)
	}
	rcSeq2 := reverseComplement(string(record2.Sequence))
	ov := analyzeOverlapPair(string(record1.Sequence), rcSeq2,
		opts.OverlapDiffLimit, opts.OverlapRequire,
		float64(opts.OverlapDiffPercentLimit)/100.0)
	isize := max
	if ov.Overlapped {
		if ov.Offset > 0 {
			isize = len(record1.Sequence) + len(record2.Sequence) - ov.OverlapLen + frontTrimmed1 + frontTrimmed2
		} else {
			isize = ov.OverlapLen + frontTrimmed1 + frontTrimmed2
		}
	}
	if isize > max {
		isize = max
	}
	s.InsertHist[isize]++
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
	return ProcessPairedEndMerge(input1, input2, output1, output2, nil, encoding, opts)
}

// ProcessPairedEndMerge processes paired-end FASTQ reads with all filters,
// additionally honouring the upstream merge writer (-m/--merge family).
// When opts.Merge is set, overlapping pairs are merged and written to
// mergeOutput; non-overlapping pairs go to output1/output2 only when
// --include_unmerged routes them there (in which case they are written to
// mergeOutput, matching upstream). When mergeOutput is nil and opts.Merge is
// set, merged reads are written to output1 (upstream's
// "compatible with fastp 0.19.8" fallback). When opts.Merge is unset,
// mergeOutput is ignored and behaviour is identical to ProcessPairedEnd.
func ProcessPairedEndMerge(input1, input2 io.Reader, output1, output2, mergeOutput io.Writer, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	opts = normalizeUMIOptions(opts)
	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)

	// Resolve the merge writer. In merge mode without an explicit
	// --merged_out, upstream falls back to writing merged reads to out1.
	var mergeWriter *fastq.Writer
	if opts.Merge {
		if mergeOutput != nil {
			mergeWriter = fastq.NewWriter(mergeOutput, encoding)
		} else {
			mergeWriter = writer1
		}
	}

	stats := &ProcessStats{MergeEnabled: opts.Merge}
	dupTracker := newDupTrackerForOpts(opts)

	// Process with multi-threading if enabled. UMI and duplication
	// tracking serialize through the input pipeline because they need a
	// deterministic per-record view, so we fall back to single-threaded
	// mode when either is active.
	if opts.Threads > 1 && !opts.UMI && opts.UMILength == 0 && dupTracker == nil &&
		!opts.Correction && !opts.OverrepAnalysis && !opts.Merge {
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

	// Overrepresentation analysis (upstream -p): buffer remaining pairs so
	// each stream's candidate hot-sequence map is built before sampling.
	if opts.OverrepAnalysis {
		for {
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
		seqs1 := make([]*fastq.Record, len(detectBuffer))
		seqs2 := make([]*fastq.Record, len(detectBuffer))
		for i, p := range detectBuffer {
			seqs1[i] = p.r1
			seqs2[i] = p.r2
		}
		stats.overrep[0] = newOverrepAnalyzer(recordSeqs(seqs1), opts.OverrepSampling, evaluatedSeqLen(seqs1))
		stats.overrep[1] = newOverrepAnalyzer(recordSeqs(seqs2), opts.OverrepSampling, evaluatedSeqLen(seqs2))
	}

	for _, p := range detectBuffer {
		if err := processPairOnce(p.r1, p.r2, writer1, writer2, mergeWriter, encoding, opts, stats, dupTracker); err != nil {
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
		if err := processPairOnce(record1, record2, writer1, writer2, mergeWriter, encoding, opts, stats, dupTracker); err != nil {
			return stats, err
		}
	}

	if err := writer1.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output1: %w", err)
	}
	if err := writer2.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output2: %w", err)
	}
	// Flush the merge writer when it is a distinct stream (when it aliases
	// writer1 it was already flushed above).
	if mergeWriter != nil && mergeWriter != writer1 {
		if err := mergeWriter.Flush(); err != nil {
			return stats, fmt.Errorf("error flushing merged output: %w", err)
		}
	}
	finalizeDupStats(stats, dupTracker)
	return stats, nil
}

// processPairOnce runs the full paired-end pipeline for one read pair,
// writing surviving reads to writer1/writer2 (or, in merge mode, the merged
// read to writer1). Shared by the streaming and split code paths.
func processPairOnce(record1, record2 *fastq.Record, writer1, writer2, mergeWriter recordWriter, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats, dupTracker *DupTracker) error {
	stats.TotalReads += 2
	stats.TotalBases += int64(len(record1.Sequence) + len(record2.Sequence))
	stats.recordBefore(record1, 0, encoding)
	stats.recordBefore(record2, 1, encoding)

	// Insert-size histogram (upstream's mInsertSizeHist). Upstream runs the
	// overlap analysis once per pair (after poly-G trimming, before adapter
	// trimming / poly-X) and bins the inferred insert size. We reproduce the
	// same overlap analysis and binning here on the pre-trim reads; the common
	// case (no poly-G, no front trim) matches upstream byte-for-byte. See
	// statInsertSize and peprocessor.cpp:698-711.
	stats.statInsertSize(record1, record2, opts, 0, 0)

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

	// Overlap-based base correction (upstream --correction). Mismatched
	// bases inside the detected PE overlap are corrected using the
	// higher-quality mate, in place, before any per-read processing.
	if opts.Correction {
		rcSeq2 := reverseComplement(string(record2.Sequence))
		ov := analyzeOverlapPair(string(record1.Sequence), rcSeq2,
			opts.OverlapDiffLimit, opts.OverlapRequire,
			float64(opts.OverlapDiffPercentLimit)/100.0)
		if ov.Overlapped {
			stats.OverlappingReads++
			if n, reads := correctByOverlapAnalysis(record1, record2, ov, encoding); n > 0 {
				stats.BasesCorrected += int64(n)
				stats.CorrectedReads += int64(reads)
			}
		}
	}

	// Overlap-driven merge writer (upstream -m/--merge). Trim each mate,
	// then run the upstream overlap analysis on the trimmed reads, merge an
	// overlapping pair into a single read, and write it to mergeWriter. When
	// the pair does not overlap, the surviving mates are dropped unless
	// --include_unmerged is set, in which case each passing mate is written
	// to the merge stream. This branch fully owns the pair's output.
	if opts.Merge {
		return processMergePair(record1, record2, mergeWriter, encoding, opts, stats)
	}

	// Legacy overlap merge path (the project's pre-upstream-port merge).
	// Kept for backward compatibility with the --merge-overlap flag.
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

// processMergePair implements upstream fastp's merging mode for one read
// pair (peprocessor.cpp:488-523). Each mate is trimmed, the upstream
// overlap analysis is run on the trimmed reads, and an overlapping pair is
// merged via mergeOverlappedPair and written to mergeWriter after passing
// the filters. A non-overlapping pair is dropped unless IncludeUnmerged is
// set, in which case each surviving mate is written to mergeWriter.
func processMergePair(record1, record2 *fastq.Record, mergeWriter recordWriter, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats) error {
	t1 := trimRecord(record1, opts, stats, encoding)
	t2 := trimRecord(record2, opts, stats, encoding)

	// Upstream auto-enables base correction whenever merging is on
	// (options.cpp:115-117), and the peprocessor applies it to the trimmed
	// reads before the merge overlap analysis. Mirror that here so the
	// merged sequence reflects corrected bases.
	if opts.Correction || opts.Merge {
		cov := analyzeOverlapPair(string(t1.Sequence), reverseComplement(string(t2.Sequence)),
			opts.OverlapDiffLimit, opts.OverlapRequire,
			float64(opts.OverlapDiffPercentLimit)/100.0)
		if cov.Overlapped {
			if n, reads := correctByOverlapAnalysis(t1, t2, cov, encoding); n > 0 {
				stats.BasesCorrected += int64(n)
				stats.CorrectedReads += int64(reads)
			}
		}
	}

	ov := analyzeOverlapPair(string(t1.Sequence), reverseComplement(string(t2.Sequence)),
		opts.OverlapDiffLimit, opts.OverlapRequire,
		float64(opts.OverlapDiffPercentLimit)/100.0)

	if ov.Overlapped {
		stats.OverlappingReads++
		merged := mergeOverlappedPair(t1, t2, ov)
		processed, pass := filterRecord(merged, opts, stats, encoding)
		if pass {
			if err := mergeWriter.Write(processed); err != nil {
				return fmt.Errorf("error writing merged read: %w", err)
			}
			stats.CleanReads++
			stats.CleanBases += int64(len(processed.Sequence))
			stats.MergedReads++
			stats.recordAfter(processed, 0, encoding)
		}
		return nil
	}

	if opts.IncludeUnmerged {
		p1, pass1 := filterRecord(t1, opts, stats, encoding)
		if pass1 {
			if err := mergeWriter.Write(p1); err != nil {
				return fmt.Errorf("error writing unmerged read1: %w", err)
			}
			stats.CleanReads++
			stats.CleanBases += int64(len(p1.Sequence))
			stats.recordAfter(p1, 0, encoding)
		}
		p2, pass2 := filterRecord(t2, opts, stats, encoding)
		if pass2 {
			if err := mergeWriter.Write(p2); err != nil {
				return fmt.Errorf("error writing unmerged read2: %w", err)
			}
			stats.CleanReads++
			stats.CleanBases += int64(len(p2.Sequence))
			stats.recordAfter(p2, 0, encoding)
		}
	}
	return nil
}

// ProcessPairedEndSplit processes paired-end FASTQ reads with all filters,
// routing surviving read pairs across numbered split files derived from
// outputBase1/outputBase2 (e.g. 0001.out1.fq / 0001.out2.fq). Split
// parameters come from opts. It buffers the input to size the splits,
// matching upstream's read-number evaluation. Merge mode is not supported
// alongside splitting (upstream rejects that combination).
func ProcessPairedEndSplit(input1, input2 io.Reader, outputBase1, outputBase2 string, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	opts = normalizeUMIOptions(opts)
	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	stats := &ProcessStats{}
	dupTracker := newDupTrackerForOpts(opts)

	records1, err := drainRecords(reader1)
	if err != nil {
		return stats, err
	}
	records2, err := drainRecords(reader2)
	if err != nil {
		return stats, err
	}
	if len(records1) != len(records2) {
		return stats, fmt.Errorf("paired files have different number of reads")
	}

	if opts.DetectAdapterPE {
		sample := len(records1)
		if sample > adapterDetectSampleSize {
			sample = adapterDetectSampleSize
		}
		pairs := make([][2]*fastq.Record, sample)
		for i := 0; i < sample; i++ {
			pairs[i] = [2]*fastq.Record{records1[i], records2[i]}
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
	if opts.OverrepAnalysis {
		stats.overrep[0] = newOverrepAnalyzer(recordSeqs(records1), opts.OverrepSampling, evaluatedSeqLen(records1))
		stats.overrep[1] = newOverrepAnalyzer(recordSeqs(records2), opts.OverrepSampling, evaluatedSeqLen(records2))
	}

	cfg := resolveSplitConfig(opts, len(records1))
	sw1 := newSplitWriter(outputBase1, cfg, encoding)
	sw2 := newSplitWriter(outputBase2, cfg, encoding)

	for i := range records1 {
		// Announce the input read index so the split writers can attribute the
		// surviving mates to the upstream worker thread (pack i/256 % threads)
		// that owns them, matching upstream's multi-threaded file boundaries.
		sw1.SetInputPos(i)
		sw2.SetInputPos(i)
		if err := processPairOnce(records1[i], records2[i], sw1, sw2, nil, encoding, opts, stats, dupTracker); err != nil {
			return stats, err
		}
	}
	if err := sw1.Close(); err != nil {
		return stats, fmt.Errorf("error closing split output1: %w", err)
	}
	if err := sw2.Close(); err != nil {
		return stats, fmt.Errorf("error closing split output2: %w", err)
	}
	finalizeDupStats(stats, dupTracker)
	return stats, nil
}

// drainRecords reads all remaining records from r into a slice. Used when
// overrepresentation analysis needs to seed a hot-sequence map from the
// full input before the main processing pass.
func drainRecords(r *fastq.Reader) ([]*fastq.Record, error) {
	var out []*fastq.Record
	for {
		rec, err := r.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, fmt.Errorf("error reading FASTQ: %w", err)
		}
		out = append(out, rec)
	}
}

// evaluatedSeqLen returns the maximum sequence length over records, the
// per-stream "evaluated read length" upstream uses to size the
// overrepresentation step windows. Returns 0 for an empty slice.
func evaluatedSeqLen(records []*fastq.Record) int {
	maxLen := 0
	for _, r := range records {
		if r != nil && len(r.Sequence) > maxLen {
			maxLen = len(r.Sequence)
		}
	}
	return maxLen
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

// ProcessSingleEnd processes single-end FASTQ reads with all filters,
// writing the surviving reads to a single output writer.
func ProcessSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	opts = normalizeUMIOptions(opts)
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	stats := &ProcessStats{MergeEnabled: opts.Merge}
	dupTracker := newDupTrackerForOpts(opts)

	// Process with multi-threading if enabled. UMI and duplication
	// tracking serialize through the input pipeline because they need a
	// deterministic per-record view, so we fall back to single-threaded
	// mode when either is active.
	if opts.Threads > 1 && !opts.UMI && opts.UMILength == 0 && dupTracker == nil && !opts.OverrepAnalysis {
		return processSingleEndParallel(reader, writer, encoding, opts, stats)
	}

	if err := runSingleEnd(reader, writer, encoding, opts, stats, dupTracker); err != nil {
		return stats, err
	}
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}
	finalizeDupStats(stats, dupTracker)
	return stats, nil
}

// ProcessSingleEndSplit processes single-end FASTQ reads with all filters,
// routing the surviving reads across numbered split files derived from
// outputBase (e.g. 0001.out.fq, 0002.out.fq). The split parameters come
// from opts (SplitNumber / SplitByLines / SplitPrefixDigits). It buffers
// the input to size the splits, matching upstream's read-number evaluation.
func ProcessSingleEndSplit(input io.Reader, outputBase string, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	opts = normalizeUMIOptions(opts)
	reader := fastq.NewReader(input, encoding)
	stats := &ProcessStats{}
	dupTracker := newDupTrackerForOpts(opts)

	records, err := drainRecords(reader)
	if err != nil {
		return stats, err
	}

	// SE adapter auto-detect over the buffered records, mirroring the
	// streaming path.
	if opts.DetectAdapterSE && opts.Adapter3 == "" {
		sample := records
		if len(sample) > adapterDetectSampleSize {
			sample = sample[:adapterDetectSampleSize]
		}
		stats.DetectedAdapter = DetectAdapterSE(sample)
		stats.DetectedAdapterR1 = stats.DetectedAdapter
		if stats.DetectedAdapter != "" {
			opts.Adapter3 = stats.DetectedAdapter
		}
	}
	if opts.OverrepAnalysis {
		stats.overrep[0] = newOverrepAnalyzer(recordSeqs(records), opts.OverrepSampling, evaluatedSeqLen(records))
	}

	cfg := resolveSplitConfig(opts, len(records))
	sw := newSplitWriter(outputBase, cfg, encoding)

	for i, rec := range records {
		// Announce the input read index so the split writer can attribute the
		// surviving read to the upstream worker thread (pack i/256 % threads)
		// that owns it, matching upstream's multi-threaded file boundaries.
		sw.SetInputPos(i)
		if err := processOneSE(rec, sw, encoding, opts, stats, dupTracker); err != nil {
			return stats, err
		}
	}
	if err := sw.Close(); err != nil {
		return stats, fmt.Errorf("error closing split output: %w", err)
	}
	finalizeDupStats(stats, dupTracker)
	return stats, nil
}

// newDupTrackerForOpts builds a duplication tracker when dup evaluation or
// dedup is requested, or returns nil.
func newDupTrackerForOpts(opts ProcessOptions) *DupTracker {
	if opts.DupCalcAccuracy <= 0 && !opts.Dedup {
		return nil
	}
	acc := opts.DupCalcAccuracy
	if acc <= 0 {
		acc = dupAccuracyDefault
	}
	return NewDupTracker(acc)
}

// runSingleEnd runs the single-end processing loop against reader, writing
// surviving records to writer. It does not flush or finalize stats; the
// caller owns those steps so split and non-split paths can share the body.
func runSingleEnd(reader *fastq.Reader, writer recordWriter, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats, dupTracker *DupTracker) error {
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
				return fmt.Errorf("error reading FASTQ: %w", err)
			}
			detectBuffer = append(detectBuffer, record)
		}
		stats.DetectedAdapter = DetectAdapterSE(detectBuffer)
		stats.DetectedAdapterR1 = stats.DetectedAdapter
		if stats.DetectedAdapter != "" && opts.Adapter3 == "" {
			opts.Adapter3 = stats.DetectedAdapter
		}
	}

	// Overrepresentation analysis (upstream -p): buffer the remaining
	// records so the candidate hot-sequence map can be built from a sample
	// of the full input before the main pass samples against it.
	if opts.OverrepAnalysis {
		rest, err := drainRecords(reader)
		if err != nil {
			return err
		}
		all := append(detectBuffer, rest...)
		detectBuffer = all
		stats.overrep[0] = newOverrepAnalyzer(recordSeqs(all), opts.OverrepSampling, evaluatedSeqLen(all))
	}

	for _, rec := range detectBuffer {
		if err := processOneSE(rec, writer, encoding, opts, stats, dupTracker); err != nil {
			return err
		}
	}

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("error reading FASTQ: %w", err)
		}
		if err := processOneSE(record, writer, encoding, opts, stats, dupTracker); err != nil {
			return err
		}
	}
	return nil
}

// processOneSE runs the full single-end pipeline for one record and writes
// it to writer if it survives filtering. Shared by the streaming and split
// code paths.
func processOneSE(record *fastq.Record, writer recordWriter, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats, dupTracker *DupTracker) error {
	stats.TotalReads++
	stats.TotalBases += int64(len(record.Sequence))
	stats.recordBefore(record, 0, encoding)

	// Duplication tracking: hash the read sequence and optionally drop
	// duplicates when --dedup is on.
	if dupTracker != nil {
		isDup := dupTracker.Observe(record.Sequence)
		if isDup && opts.Dedup {
			stats.DedupDropped++
			return nil
		}
	}

	// Extract UMI if configured.
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

// processRecord applies all processing steps (trimming + filtering) to a
// single record, returning the surviving record and whether it passed.
func processRecord(record *fastq.Record, opts ProcessOptions, stats *ProcessStats, encoding fastq.QualityEncoding) (*fastq.Record, bool) {
	trimmed := trimRecord(record, opts, stats, encoding)
	return filterRecord(trimmed, opts, stats, encoding)
}

// trimRecord applies the per-read trimming steps (base correction, adapter
// trimming including --adapter_fasta, poly-G, poly-X, sliding-window
// quality cut, and quality-based end trimming) and returns the trimmed
// record. It never rejects a read; length/N/quality/complexity filtering is
// the job of filterRecord. Splitting the two mirrors upstream fastp, which
// trims each mate before the merge/overlap analysis and only then applies
// passFilter.
func trimRecord(record *fastq.Record, opts ProcessOptions, stats *ProcessStats, encoding fastq.QualityEncoding) *fastq.Record {
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

	// Step 1: Sliding-window quality trimming (cut_front / cut_tail /
	// cut_right) — upstream's Filter::trimAndCut.
	//
	// ORDER IS LOAD-BEARING: upstream runs trimAndCut FIRST, before poly-G,
	// adapter trimming, and poly-X (seprocessor.cpp:235 / peprocessor.cpp).
	// Running it after adapter trimming (as a previous version of this port
	// did) makes the sliding-window see an already-shortened read, which
	// shifts the window-boundary math (the `s+w<l-tail` bound and the
	// last-partial-window handling) and diverges from upstream by ~1bp on
	// any read where adapter/poly trimming fired. Keep this step first.
	if opts.CutFront || opts.CutTail || opts.CutRight {
		cutLo, cutHi := slidingWindowCut([]byte(seq[start:end]), qual[start:end], encoding, opts)
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

	// Step 2: Trim poly-G tails if enabled.
	//
	// ORDER: upstream applies poly-G AFTER the sliding-window cut but
	// BEFORE adapter trimming (seprocessor.cpp:239 runs before the adapter
	// block at :242). Uses trimPolyG, a verbatim port of upstream's
	// PolyX::trimPolyG (reference_code/fastp/src/polyx.cpp:16-42). The
	// upstream algorithm tolerates 1 mismatch per 8 bases scanned (capped
	// at 5 total) and anchors the trim at the last-G position seen.
	if opts.TrimPolyG && !opts.DisableTrimPolyG {
		newEnd := trimPolyG(seq[start:end], opts.PolyGMinLen)
		if newEnd < end-start {
			polyLen := (end - start) - newEnd
			end = start + newEnd
			stats.PolyGTrimmedReads++
			stats.PolyGTrimmedBases += int64(polyLen)
		}
	}

	// Step 3: Trim adapters if specified. Suppressed entirely when
	// --disable_adapter_trimming is set (upstream gates the whole adapter
	// block on adapter.enabled == !disable_adapter_trimming).
	if !opts.DisableAdapterTrimming {
		explicitTrimmed := false
		if opts.Adapter5 != "" {
			pos := findAdapter(seq, opts.Adapter5)
			if pos >= 0 {
				start = pos + len(opts.Adapter5)
				stats.AdapterTrimmedReads++
				stats.AdapterTrimmedBases += int64(pos + len(opts.Adapter5))
				explicitTrimmed = true
			}
		}

		if opts.Adapter3 != "" {
			// Route the single configured/auto-detected 3' adapter through
			// the verbatim upstream AdapterTrimmer::trimBySequence algorithm
			// (adaptertrimmer.cpp:71-170) rather than a plain substring
			// search. This matters for parity: trimBySequence tolerates one
			// mismatch per 8 overlapping bases, matches a PARTIAL adapter
			// prefix at the 3' end (cmplen = min(rlen-pos, alen)), and uses
			// the negative "A-tailing" start offset. A plain strings.Index
			// missed read-through cases where only the adapter prefix fits
			// inside the read. Upstream calls this with the default
			// matchReq = 4 (seprocessor.cpp:245).
			window := &fastq.Record{
				Sequence: []byte(seq[start:end]),
				Quality:  qual[start:end],
			}
			before := len(window.Sequence)
			if trimBySequenceUpstream(window, opts.Adapter3, false, 4) {
				removed := before - len(window.Sequence)
				end -= removed
				stats.AdapterTrimmedReads++
				stats.AdapterTrimmedBases += int64(removed)
				explicitTrimmed = true
			}
		}

		// --adapter_fasta: trim the current window against every loaded
		// adapter, mirroring upstream's AdapterTrimmer::trimByMultiSequences.
		// The read counter is incremented only when the explicit-adapter
		// pass did NOT already trim (and thus already count) this read —
		// upstream passes incTrimmedCounter = !trimmed (seprocessor.cpp:246,
		// filterresult.cpp:124-128). Bases are always accumulated.
		if len(opts.AdapterFasta) > 0 && end > start {
			window := &fastq.Record{
				Sequence: []byte(seq[start:end]),
				Quality:  qual[start:end],
			}
			if n, trimmed := trimByMultiSequences(window, opts.AdapterFasta, false); trimmed {
				end -= n
				if !explicitTrimmed {
					stats.AdapterTrimmedReads++
				}
				stats.AdapterTrimmedBases += int64(n)
			}
		}
	}

	// Step 4: Trim poly-X tails if enabled.
	//
	// Uses trimPolyXUpstream, a verbatim port of upstream's PolyX::trimPolyX
	// (reference_code/fastp/src/polyx.cpp:49-116): it tolerates 1 mismatch
	// per 8 bases scanned (capped at 5), treats N as every base, and picks
	// the most frequent of A/T/C/G as the poly base. The length threshold is
	// --poly_x_min_len, defaulting to 10 but independent of poly-G's knob.
	if opts.TrimPolyX {
		polyXMin := opts.PolyXMinLen
		if polyXMin <= 0 {
			polyXMin = 10
		}
		newEnd, trimmed := trimPolyXUpstream(seq[start:end], polyXMin)
		if trimmed > 0 {
			end = start + newEnd
		}
	}

	// NOTE: there is deliberately NO standalone end quality-trim here. Upstream
	// fastp's only quality-based trimming is the cut_front/cut_tail/cut_right
	// sliding-window cut, handled in Step 1 above (gated on those flags, off by
	// default). qualified_quality_phred (-q, default 15) is a FILTER threshold
	// used by the quality-percentage filter in filterRecord — NOT a trim
	// threshold. A previous standalone trimByQuality step here ran by default
	// (QualThreshold defaults to 15) and trimmed low-quality 3' tails upstream
	// leaves intact, which both shifted lengths and let too-many-N reads slip
	// through the N filter (their N-laden tails were trimmed away). See upstream
	// seprocessor.cpp:234-262, which has no such step.

	// Return the trimmed window. Length/N/quality/complexity filtering is
	// deferred to filterRecord so the merge path can interpose its overlap
	// analysis between trimming and filtering (as upstream does).
	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence[start:end],
		Quality:     record.Quality[start:end],
	}
}

// filterRecord applies the length, N-content, quality-percentage, and
// complexity filters to an already-trimmed record. It returns the record
// (possibly resized by --length_limit) and whether it passed all filters.
func filterRecord(record *fastq.Record, opts ProcessOptions, stats *ProcessStats, encoding fastq.QualityEncoding) (*fastq.Record, bool) {
	seq := string(record.Sequence)
	qual := record.Quality
	start := 0
	end := len(seq)

	// Length filtering (upstream filter.cpp:52-57, gated on
	// lengthFilter.enabled). --disable_length_filtering (-L) skips the
	// too-short and too-long checks. Note --length_limit truncation is a
	// distinct trimming step (max_len in upstream's Read::resize path) and
	// is intentionally NOT gated by the length filter toggle.
	if !opts.DisableLengthFiltering {
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
	}

	if opts.LengthLimit > 0 && end-start > opts.LengthLimit {
		end = start + opts.LengthLimit
	}

	// Quality filtering (upstream filter.cpp:43-50, gated on
	// qualfilter.enabled). --disable_quality_filtering (-Q) skips the
	// N-base-count limit AND the low-quality-base percentage limit (and
	// the average-quality requirement, once ported).
	if !opts.DisableQualityFiltering {
		// Step 5: Check N content (nBaseLimit / N-percent).
		nCount := countNs(seq[start:end])
		nPercent := 100.0 * float64(nCount) / float64(end-start)

		if nCount > opts.MaxNCount || nPercent > opts.MaxNPercent {
			stats.TooManyNReads++
			return nil, false
		}

		// Step 6: Check quality (percentage of bases meeting threshold).
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
	}

	// Step 7: Check complexity if enabled.
	//
	// Upstream fastp defines sequence complexity as the fraction of
	// adjacent base pairs that differ: `diff/(length-1)` where diff is
	// the count of indices i where seq[i] != seq[i+1]. Threshold defaults
	// to 0.3 (== upstream's --complexity_threshold 30, which is interpreted
	// as a percentage). We accept the value as a fraction in [0,1]; the
	// CLI also accepts the percentage form via cliflag's float parser, and
	// we normalise values > 1 down by dividing by 100 so `--complexity-threshold 30`
	// behaves like `--complexity-threshold 0.3` (matching upstream's `-Y 30`).
	if opts.LowComplexity {
		complexity := calculateComplexity(seq[start:end])
		threshold := opts.ComplexityThreshold
		if threshold > 1 {
			threshold = threshold / 100
		}
		if complexity < threshold {
			// Low complexity read — discard. Mirrors upstream fastp which
			// also counts these under a dedicated `low_complexity_reads`
			// bucket in the filtering result.
			stats.LowComplexityReads++
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

// phredOffset returns the ASCII offset for the given quality encoding.
func phredOffset(encoding fastq.QualityEncoding) int {
	if encoding == fastq.Phred64 {
		return 64
	}
	return 33
}

// slidingWindowCut applies fastp's sliding-window quality trimming and
// returns the half-open index range [lo, hi) of bases that should be kept
// in seq / quality.
//
// This is a verbatim port of upstream fastp's Filter::trimAndCut
// (reference_code/fastp/src/filter.cpp:83-222), with PE pre-trim front/tail
// hard-wired to 0 (the Go caller passes the already-clipped slice). The
// three modes are applied in upstream's order:
//
//   - cut_front (filter.cpp:111-142): walk a window 5'->3' until the
//     window mean quality is >= CutMeanQuality. Then advance past the
//     window minus one base, then skip any leading 'N's. Keep [s, end).
//
//   - cut_right (filter.cpp:144-178): walk a window 5'->3'; on the first
//     window with mean below threshold, walk the bad window keeping the
//     high-Q prefix (qualstr[s] >= threshold). Trim from there onward.
//
//   - cut_tail (filter.cpp:180-209): only runs when cut_right is OFF.
//     Walks a window 3'->5'; on the first window with mean >= threshold,
//     trim everything past the window's START (then skip trailing N's).
//
// Critically, the upstream loop bound is strictly `s + w < l` (not
// `s + w <= l`), so the very last window of the read is never scanned.
// Replicated here for byte-for-byte parity with the C++ implementation.
//
// The cut_front asymmetric trim - front jumps to s+w-1 once we find the
// qualifying window, dropping (w-1) bases of that window - is also part
// of the upstream behavior and is preserved verbatim.
func slidingWindowCut(seq, quality []byte, encoding fastq.QualityEncoding, opts ProcessOptions) (lo, hi int) {
	offset := phredOffset(encoding)
	window := opts.CutWindowSize
	if window < 1 {
		window = 1
	}
	l := len(quality)
	front := 0
	tail := 0
	rlen := l

	// q returns the integer Phred score at position i in the (full)
	// quality slice. Mirrors upstream's `qualstr[i] - 33` math; we use
	// the configured phred offset (33 or 64).
	q := func(i int) int { return int(quality[i]) - offset }

	// CUT FRONT - filter.cpp:111-142.
	if opts.CutFront {
		w := window
		// l - front - tail - w <= 0 -> nothing to do; upstream returns NULL.
		if l-front-tail-w > 0 {
			totalQual := 0
			s := front
			// preparing rolling: sum w-1 leading qualities.
			for i := 0; i < w-1; i++ {
				totalQual += q(s + i)
			}
			thresh := float64(opts.CutMeanQuality)
			for s = front; s+w < l-tail; s++ {
				totalQual += q(s + w - 1)
				if s > front {
					totalQual -= q(s - 1)
				}
				if float64(totalQual)/float64(w) >= thresh {
					break
				}
			}
			// the trimming in front is forwarded and rlen is recalculated
			if s > 0 {
				s = s + w - 1
			}
			for s < l && seq[s] == 'N' {
				s++
			}
			front = s
			rlen = l - front - tail
		}
	}

	// CUT RIGHT - filter.cpp:144-178.
	if opts.CutRight {
		w := window
		if l-front-tail-w > 0 {
			totalQual := 0
			s := front
			for i := 0; i < w-1; i++ {
				totalQual += q(s + i)
			}
			thresh := float64(opts.CutMeanQuality)
			foundLowQualWindow := false
			for s = front; s+w < l-tail; s++ {
				totalQual += q(s + w - 1)
				if s > front {
					totalQual -= q(s - 1)
				}
				if float64(totalQual)/float64(w) < thresh {
					foundLowQualWindow = true
					break
				}
			}
			if foundLowQualWindow {
				// keep the good bases in the (bad) window
				for s < l-1 && q(s) >= opts.CutMeanQuality {
					s++
				}
				rlen = s - front
			}
		}
	}

	// CUT TAIL - filter.cpp:180-209. Suppressed when cut_right is on.
	if !opts.CutRight && opts.CutTail {
		w := window
		if l-front-tail-w > 0 {
			totalQual := 0
			t := l - tail - 1
			// preparing rolling: sum w-1 trailing qualities.
			for i := 0; i < w-1; i++ {
				totalQual += q(t - i)
			}
			thresh := float64(opts.CutMeanQuality)
			for t = l - tail - 1; t-w >= front; t-- {
				totalQual += q(t - w + 1)
				if t < l-tail-1 {
					totalQual -= q(t + 1)
				}
				if float64(totalQual)/float64(w) >= thresh {
					break
				}
			}
			if t < l-1 {
				t = t - w + 1
			}
			for t >= 0 && seq[t] == 'N' {
				t--
			}
			rlen = t - front + 1
		}
	}

	// Upstream returns NULL (i.e. the read is dropped) when rlen <= 0 or
	// front >= l-1. We model that by emitting a degenerate range, which
	// the caller's MinLength check then rejects.
	if rlen <= 0 || front >= l-1 {
		return l, l
	}
	return front, front + rlen
}

// trimPolyG trims a 3' poly-G run from seq and returns the new length
// (the index at which to truncate). It is a verbatim port of upstream
// fastp's PolyX::trimPolyG (reference_code/fastp/src/polyx.cpp:16-42).
//
// The algorithm scans bases right-to-left, allowing one mismatch per 8
// bases scanned (capped at 5 total). It tracks `firstGPos`, the
// leftmost-G seen so far in the run; once we accumulate enough
// mismatches that the run can't reasonably continue (and we've already
// scanned at least compareReq bases) we stop and truncate at firstGPos.
//
// Returns len(seq) if no trim should be applied.
func trimPolyG(seq string, compareReq int) int {
	const allowOneMismatchForEach = 8
	const maxMismatch = 5

	rlen := len(seq)
	mismatch := 0
	i := 0
	firstGPos := rlen - 1
	for i = 0; i < rlen; i++ {
		b := seq[rlen-i-1]
		if b != 'G' && b != 'g' {
			mismatch++
		} else {
			firstGPos = rlen - i - 1
		}
		allowedMismatch := (i + 1) / allowOneMismatchForEach
		if mismatch > maxMismatch || (mismatch > allowedMismatch && i >= compareReq-1) {
			break
		}
	}
	if i >= compareReq {
		return firstGPos
	}
	return rlen
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

// calculateComplexity returns the fraction of adjacent base pairs that
// differ — i.e. count(i : seq[i] != seq[i+1]) / (len(seq)-1). This matches
// upstream fastp's passLowComplexityFilter definition (see filter.cpp).
// A run of identical bases returns 0; a perfectly alternating ATAT...
// sequence returns 1.0. Sequences shorter than 2 bases return 0 (upstream
// also rejects them).
func calculateComplexity(seq string) float64 {
	if len(seq) <= 1 {
		return 0
	}
	diff := 0
	for i := 0; i < len(seq)-1; i++ {
		if seq[i] != seq[i+1] {
			diff++
		}
	}
	return float64(diff) / float64(len(seq)-1)
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
