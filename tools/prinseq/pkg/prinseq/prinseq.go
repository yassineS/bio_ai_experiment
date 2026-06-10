package prinseq

import (
	"bufio"
	"fmt"
	"io"
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// Stats holds sequence statistics
type Stats struct {
	NumReads            int            `json:"num_reads"`
	TotalBases          int            `json:"total_bases"`
	MinLength           int            `json:"min_length"`
	MaxLength           int            `json:"max_length"`
	AvgLength           float64        `json:"avg_length"`
	GCContent           float64        `json:"gc_content"`
	AvgQuality          float64        `json:"avg_quality,omitempty"` // Only for FASTQ
	NumNs               int            `json:"num_ns"`
	LengthDistribution  map[int]int    `json:"length_distribution,omitempty"`
	QualityDistribution map[int]int    `json:"quality_distribution,omitempty"`
	BaseComposition     map[string]int `json:"base_composition,omitempty"`
	Dinucleotides       map[string]int `json:"dinucleotides,omitempty"`
	PositionalQuality   []float64      `json:"positional_quality,omitempty"`
}

// CalculateStats computes statistics for FASTA or FASTQ files
func CalculateStats(reader io.Reader, isFastq bool) (*Stats, error) {
	return CalculateStatsWithEncoding(reader, isFastq, "sanger")
}

// CalculateStatsWithEncoding computes statistics with a specific quality encoding
func CalculateStatsWithEncoding(reader io.Reader, isFastq bool, qualType string) (*Stats, error) {
	stats := &Stats{
		MinLength: -1,
	}

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 1024*1024), 1024*1024) // 1MB buffer

	if isFastq {
		offset := 33
		if qualType == "illumina" {
			offset = 64
		}
		return calculateFastqStatsWithOffset(scanner, stats, offset)
	}
	return calculateFastaStats(scanner, stats)
}

func calculateFastaStats(scanner *bufio.Scanner, stats *Stats) (*Stats, error) {
	var currentSeq string
	gcCount := 0

	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			// Process previous sequence
			if currentSeq != "" {
				processSequence(currentSeq, &gcCount, stats)
				currentSeq = ""
			}
		} else {
			currentSeq += line
		}
	}

	// Process last sequence
	if currentSeq != "" {
		processSequence(currentSeq, &gcCount, stats)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
	}

	return stats, nil
}

func calculateFastqStats(scanner *bufio.Scanner, stats *Stats) (*Stats, error) {
	return calculateFastqStatsWithOffset(scanner, stats, 33)
}

func calculateFastqStatsWithOffset(scanner *bufio.Scanner, stats *Stats, offset int) (*Stats, error) {
	gcCount := 0
	totalQuality := 0.0
	lineNum := 0

	var seq string
	var qual string

	for scanner.Scan() {
		line := scanner.Text()
		mod := lineNum % 4

		switch mod {
		case 0: // Header line
			if len(line) == 0 || line[0] != '@' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 1: // Sequence line
			seq = line
		case 2: // Plus line
			if len(line) == 0 || line[0] != '+' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 3: // Quality line
			qual = line
			if len(qual) != len(seq) {
				return nil, fmt.Errorf("quality length (%d) doesn't match sequence length (%d) at line %d",
					len(qual), len(seq), lineNum+1)
			}
			// Process the complete FASTQ record
			processSequence(seq, &gcCount, stats)
			totalQuality += calculateAvgQualityScoreWithOffset(qual, offset)
			seq = ""
			qual = ""
		}

		lineNum++
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Check if we have incomplete records
	if lineNum%4 != 0 {
		return nil, fmt.Errorf("incomplete FASTQ record")
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
		stats.AvgQuality = totalQuality / float64(stats.NumReads)
	}

	return stats, nil
}

func processSequence(seq string, gcCount *int, stats *Stats) {
	seqLen := len(seq)
	stats.NumReads++
	stats.TotalBases += seqLen

	if stats.MinLength == -1 || seqLen < stats.MinLength {
		stats.MinLength = seqLen
	}
	if seqLen > stats.MaxLength {
		stats.MaxLength = seqLen
	}

	// Count GC and Ns
	for _, base := range seq {
		switch base {
		case 'G', 'C', 'g', 'c':
			*gcCount++
		case 'N', 'n':
			stats.NumNs++
		}
	}
}

func calculateAvgQualityScore(qual string) float64 {
	return calculateAvgQualityScoreWithOffset(qual, 33)
}

func calculateAvgQualityScoreWithOffset(qual string, offset int) float64 {
	if len(qual) == 0 {
		return 0.0
	}

	total := 0
	for _, q := range qual {
		// Use specified offset (33 for Phred+33, 64 for Phred+64)
		total += int(q) - offset
	}
	return float64(total) / float64(len(qual))
}

// trimSequence applies various trimming operations to a sequence
func trimSequence(seq, qual string, opts FilterOptions) (string, string) {
	seqBytes := []byte(seq)
	var qualBytes []byte
	if qual != "" {
		qualBytes = []byte(qual)
	}

	// Apply trim_left (fixed position from left)
	if opts.TrimLeft > 0 && len(seqBytes) > opts.TrimLeft {
		seqBytes = seqBytes[opts.TrimLeft:]
		if len(qualBytes) > 0 {
			qualBytes = qualBytes[opts.TrimLeft:]
		}
	}

	// Apply trim_right (fixed position from right)
	if opts.TrimRight > 0 && len(seqBytes) > opts.TrimRight {
		seqBytes = seqBytes[:len(seqBytes)-opts.TrimRight]
		if len(qualBytes) > 0 {
			qualBytes = qualBytes[:len(qualBytes)-opts.TrimRight]
		}
	}

	// Apply trim_left_p (percentage from left)
	if opts.TrimLeftP > 0 && len(seqBytes) > 0 {
		trimPos := (len(seqBytes) * opts.TrimLeftP) / 100
		if trimPos > 0 && trimPos < len(seqBytes) {
			seqBytes = seqBytes[trimPos:]
			if len(qualBytes) > 0 {
				qualBytes = qualBytes[trimPos:]
			}
		}
	}

	// Apply trim_right_p (percentage from right)
	if opts.TrimRightP > 0 && len(seqBytes) > 0 {
		trimPos := (len(seqBytes) * opts.TrimRightP) / 100
		if trimPos > 0 && trimPos < len(seqBytes) {
			seqBytes = seqBytes[:len(seqBytes)-trimPos]
			if len(qualBytes) > 0 {
				qualBytes = qualBytes[:len(qualBytes)-trimPos]
			}
		}
	}

	// Apply trim_ns_left (trim poly-N from left)
	if opts.TrimNsLeft > 0 {
		seqBytes, qualBytes = trimPolyNLeft(seqBytes, qualBytes, opts.TrimNsLeft)
	}

	// Apply trim_ns_right (trim poly-N from right)
	if opts.TrimNsRight > 0 {
		seqBytes, qualBytes = trimPolyNRight(seqBytes, qualBytes, opts.TrimNsRight)
	}

	// Apply trim_tail_left (trim poly-A/T from left)
	if opts.TrimTailLeft > 0 {
		seqBytes, qualBytes = trimPolyATLeft(seqBytes, qualBytes, opts.TrimTailLeft)
	}

	// Apply trim_tail_right (trim poly-A/T from right)
	if opts.TrimTailRight > 0 {
		seqBytes, qualBytes = trimPolyATRight(seqBytes, qualBytes, opts.TrimTailRight)
	}

	// Apply sliding-window quality-based trimming. Upstream evaluates
	// the left side first, then the right side over the already
	// left-trimmed read (prinseq-lite.pl:3215-3287).
	if (opts.TrimQualL > 0 || opts.TrimQualR > 0) && len(qualBytes) > 0 {
		seqBytes, qualBytes = trimQualityWindow(seqBytes, qualBytes, opts)
	}

	// Apply trim_to_len: hard-trim to a maximum length from the 5' end.
	// Upstream runs this AFTER quality trimming and BEFORE the length/GC
	// filters (prinseq-lite.pl:3382-3385). The `length > trim_to_len`
	// guard means equal-or-shorter reads are left untouched.
	if opts.TrimToLen > 0 && len(seqBytes) > opts.TrimToLen {
		seqBytes = seqBytes[:opts.TrimToLen]
		if len(qualBytes) > 0 {
			qualBytes = qualBytes[:opts.TrimToLen]
		}
	}

	return string(seqBytes), string(qualBytes)
}

// trimQualityWindow performs prinseq's sliding-window quality trimming on
// both ends of a read. It mirrors prinseq-lite.pl:3215-3287 exactly:
//
//   - The window starts at the read end and advances by TrimQualStep each
//     time the rule is satisfied. The window shrinks to the remaining bases
//     near the far end (so the last partial window can be size 1).
//   - The window's score is aggregated by TrimQualType: "min" (default),
//     "max", "mean", or "sum". A window of size 1 always uses the raw value.
//     For "sum", a partial final window (smaller than TrimQualWindow) ends
//     the loop without trimming further.
//   - The rule compares the score against the per-end threshold:
//     "lt" trims while score < threshold, "gt" while score > threshold,
//     "et" while score == threshold.
//
// The left side is processed first; the right side then operates on the
// already left-trimmed read using a reversed view of the quality scores,
// matching upstream's `reverse(@$qualsnums)` and `$length -= $begintmp`.
func trimQualityWindow(seq, qual []byte, opts FilterOptions) ([]byte, []byte) {
	nums := decodeQual(qual, phredOffset(opts.QualType))

	window := opts.TrimQualWindow
	if window <= 0 {
		window = 1
	}
	step := opts.TrimQualStep
	if step <= 0 {
		step = 1
	}
	qtype := opts.TrimQualType
	if qtype == "" {
		qtype = "min"
	}
	rule := opts.TrimQualRule
	if rule == "" {
		rule = "lt"
	}

	length := len(nums)

	// Left.
	beginTmp := 0
	if opts.TrimQualL > 0 {
		beginTmp = scanQualTrim(nums, length, window, step, qtype, rule, opts.TrimQualL)
		if beginTmp >= length {
			// Whole read trimmed away; the length filter will drop it.
			beginTmp = length
		}
		if beginTmp > 0 {
			seq = seq[beginTmp:]
			qual = qual[beginTmp:]
		}
	}

	// Right. Operate on the read remaining after the left trim, using a
	// reversed copy of the quality scores (upstream reverses the FULL
	// quality array, then limits the scan to length-beginTmp).
	if opts.TrimQualR > 0 {
		rightLen := length - beginTmp
		rev := make([]int, len(nums))
		for i := range nums {
			rev[i] = nums[len(nums)-1-i]
		}
		endTmp := scanQualTrim(rev, rightLen, window, step, qtype, rule, opts.TrimQualR)
		if endTmp >= rightLen {
			endTmp = rightLen
		}
		if endTmp > 0 && endTmp <= len(seq) {
			seq = seq[:len(seq)-endTmp]
			qual = qual[:len(qual)-endTmp]
		}
	}

	return seq, qual
}

// scanQualTrim returns the number of bases to trim from the start of nums
// (limited to the first `length` entries) under prinseq's window/step/type/
// rule scheme. It is shared by the left and right passes; the right pass
// supplies a reversed score array. See trimQualityWindow for the semantics.
func scanQualTrim(nums []int, length, window, step int, qtype, rule string, threshold int) int {
	trimmed := 0
	i := 0
	for i < length {
		// The effective window shrinks near the far end.
		w := window
		if i+window > length {
			w = length - i
		}
		var val float64
		switch {
		case w == 1:
			val = float64(nums[i])
		case qtype == "min":
			val = float64(minInt(nums[i : i+w]))
		case qtype == "max":
			val = float64(maxInt(nums[i : i+w]))
		case qtype == "mean":
			val = meanInt(nums[i : i+w])
		case qtype == "sum":
			if w < window {
				// Partial final window: upstream `last`s without trimming.
				return trimmed
			}
			val = float64(sumInt(nums[i : i+w]))
		default:
			return trimmed
		}
		if qualRuleMatches(rule, val, threshold) {
			trimmed += step
			i += step
		} else {
			break
		}
	}
	return trimmed
}

// qualRuleMatches implements prinseq's trim_qual_rule comparison: "lt" is
// true when val < threshold, "gt" when val > threshold, "et" when
// val == threshold (prinseq-lite.pl:3236).
func qualRuleMatches(rule string, val float64, threshold int) bool {
	t := float64(threshold)
	switch rule {
	case "gt":
		return val > t
	case "et":
		return val == t
	default: // "lt"
		return val < t
	}
}

// minInt returns the smallest element of a non-empty slice.
func minInt(xs []int) int {
	m := xs[0]
	for _, x := range xs[1:] {
		if x < m {
			m = x
		}
	}
	return m
}

// maxInt returns the largest element of a non-empty slice.
func maxInt(xs []int) int {
	m := xs[0]
	for _, x := range xs[1:] {
		if x > m {
			m = x
		}
	}
	return m
}

// sumInt returns the sum of the slice's elements.
func sumInt(xs []int) int {
	s := 0
	for _, x := range xs {
		s += x
	}
	return s
}

// meanInt returns the arithmetic mean of a non-empty slice as a float, the
// same value upstream's getArrayMean produces.
func meanInt(xs []int) float64 {
	return float64(sumInt(xs)) / float64(len(xs))
}

// trimQualityLeft trims bases from the 5' end while each base's decoded
// quality score is below threshold. It is the window=1/step=1/type=min/
// rule=lt special case of trimQualityWindow and is retained as a small,
// directly-testable entry point. The whole-read case (everything below
// threshold) leaves the read untouched, matching the original behaviour
// the package's filter loop relies on.
func trimQualityLeft(seq, qual []byte, threshold, offset int) ([]byte, []byte) {
	nums := decodeQual(qual, offset)
	trimmed := scanQualTrim(nums, len(nums), 1, 1, "min", "lt", threshold)
	if trimmed > 0 && trimmed < len(seq) {
		seq = seq[trimmed:]
		qual = qual[trimmed:]
	}
	return seq, qual
}

// trimQualityRight trims bases from the 3' end while each base's decoded
// quality score is below threshold. It mirrors trimQualityLeft on the
// reversed read.
func trimQualityRight(seq, qual []byte, threshold, offset int) ([]byte, []byte) {
	nums := decodeQual(qual, offset)
	rev := make([]int, len(nums))
	for i := range nums {
		rev[i] = nums[len(nums)-1-i]
	}
	trimmed := scanQualTrim(rev, len(rev), 1, 1, "min", "lt", threshold)
	if trimmed > 0 && trimmed < len(seq) {
		seq = seq[:len(seq)-trimmed]
		qual = qual[:len(qual)-trimmed]
	}
	return seq, qual
}

// decodeQual converts ASCII quality bytes to integer Phred scores under the
// given offset.
func decodeQual(qual []byte, offset int) []int {
	nums := make([]int, len(qual))
	for i := range qual {
		nums[i] = int(qual[i]) - offset
	}
	return nums
}

func trimPolyNLeft(seq, qual []byte, minLen int) ([]byte, []byte) {
	nCount := 0
	for i := 0; i < len(seq); i++ {
		if seq[i] == 'N' || seq[i] == 'n' {
			nCount++
		} else {
			break
		}
	}

	if nCount >= minLen {
		seq = seq[nCount:]
		if len(qual) > 0 {
			qual = qual[nCount:]
		}
	}
	return seq, qual
}

func trimPolyNRight(seq, qual []byte, minLen int) ([]byte, []byte) {
	nCount := 0
	for i := len(seq) - 1; i >= 0; i-- {
		if seq[i] == 'N' || seq[i] == 'n' {
			nCount++
		} else {
			break
		}
	}

	if nCount >= minLen {
		seq = seq[:len(seq)-nCount]
		if len(qual) > 0 {
			qual = qual[:len(qual)-nCount]
		}
	}
	return seq, qual
}

// trimPolyATLeft trims a poly-A or poly-T tail from the 5'-end of seq when
// the run is at least minLen bases long. The logic mirrors upstream prinseq:
// (1) match either `^A{minLen}` or `^T{minLen}`, picking whichever exists;
// (2) extend the run with more of the SAME base or with Ns; (3) if the run
// covers the entire read, trim everything (the caller will length-filter it
// out separately, matching upstream's `good = 0` branch).
//
// Before this fix we treated A and T as interchangeable, so a run that
// alternated A and T (or that started with A and continued with T) would
// over-trim by one or more bases relative to upstream.
func trimPolyATLeft(seq, qual []byte, minLen int) ([]byte, []byte) {
	if len(seq) < minLen {
		return seq, qual
	}
	// Detect a homogeneous A-run or T-run of length minLen at the head.
	var anchor byte
	switch {
	case allEqualCase(seq[:minLen], 'A'):
		anchor = 'A'
	case allEqualCase(seq[:minLen], 'T'):
		anchor = 'T'
	default:
		return seq, qual
	}
	// Extend with the same base or N.
	end := minLen
	for end < len(seq) && (matchesCase(seq[end], anchor) || seq[end] == 'N' || seq[end] == 'n') {
		end++
	}
	seq = seq[end:]
	if len(qual) > 0 {
		qual = qual[end:]
	}
	return seq, qual
}

// trimPolyATRight is the 3'-end counterpart of trimPolyATLeft.
func trimPolyATRight(seq, qual []byte, minLen int) ([]byte, []byte) {
	if len(seq) < minLen {
		return seq, qual
	}
	tail := seq[len(seq)-minLen:]
	var anchor byte
	switch {
	case allEqualCase(tail, 'A'):
		anchor = 'A'
	case allEqualCase(tail, 'T'):
		anchor = 'T'
	default:
		return seq, qual
	}
	start := len(seq) - minLen
	for start > 0 && (matchesCase(seq[start-1], anchor) || seq[start-1] == 'N' || seq[start-1] == 'n') {
		start--
	}
	seq = seq[:start]
	if len(qual) > 0 {
		qual = qual[:start]
	}
	return seq, qual
}

// matchesCase reports whether b equals base or its lowercase form.
func matchesCase(b, base byte) bool {
	return b == base || b == (base|0x20)
}

// allEqualCase reports whether every byte in s equals base or its
// lowercase form.
func allEqualCase(s []byte, base byte) bool {
	for _, b := range s {
		if !matchesCase(b, base) {
			return false
		}
	}
	return true
}

// phredOffset returns the ASCII offset used to decode quality characters for
// the given prinseq quality-type string. It returns 64 for Illumina 1.3-1.7
// style encodings ("illumina") and 33 for Sanger/Phred+33 (the default).
func phredOffset(qualType string) int {
	if qualType == "illumina" {
		return 64
	}
	return 33
}

// FilterOptions holds filtering parameters
type FilterOptions struct {
	MinLen        int
	MaxLen        int
	MinGC         float64
	MaxGC         float64
	MinQual       float64
	MaxNsP        float64 // Max percentage of Ns (upstream --ns_max_p)
	MaxNsN        int     // Max number of Ns (upstream --ns_max_n)
	TrimLeft      int
	TrimRight     int
	TrimLeftP     int // Trim percentage from left
	TrimRightP    int // Trim percentage from right
	TrimQualL     int // Quality threshold for left trimming
	TrimQualR     int // Quality threshold for right trimming
	TrimNsLeft    int // Trim poly-N from left
	TrimNsRight   int // Trim poly-N from right
	TrimTailLeft  int // Trim poly-A/T from left
	TrimTailRight int // Trim poly-A/T from right
	MinQualMean   float64
	MaxQualMean   float64
	Derep         int       // Duplicate removal mode (1=exact, 4=reverse complement)
	DerepMin      int       // Minimum occurrences to keep
	QualType      string    // Quality encoding type: "sanger" (Phred+33) or "illumina" (Phred+64)
	OutBad        io.Writer // Writer for rejected sequences (optional)
	LcMethod      string    // Low complexity method: "dust" or "entropy"
	LcThreshold   float64   // Low complexity threshold (7 for dust, 70 for entropy)

	// Upstream PRINSEQ-lite parity flags (added in PR #prinseq-missing-flags):
	//
	//   NonIUPAC    — implements `-noniupac` (reject any base outside ACGTN,
	//                 case-insensitive; upstream upper-cases before the
	//                 [^ACGTN] check at prinseq-lite.pl line 3478).
	//   SeqID       — implements `-seq_id <prefix>` (rename headers to
	//                 "<prefix><counter>"; counter starts at 1 and only
	//                 increments for records that pass all filters, matching
	//                 the upstream behaviour at line 3648).
	//   SeqIDMap    — implements `-seq_id_mappings <file>`; emits
	//                 "<orig_id>\t<new_id>\n" for each renamed record
	//                 (upstream line 3646). Requires SeqID to be non-empty;
	//                 callers should validate that at the CLI layer.
	//   OutFormat   — implements `-out_format` (1=FASTA, 2=FASTA+QUAL,
	//                 3=FASTQ, 4=FASTQ+FASTA, 5=FASTQ+FASTA+QUAL). When
	//                 zero, the output format mirrors the input format
	//                 (matches the inferred default at lines 785-789).
	//   FastaOut    — secondary FASTA writer for `--out_format 4 / 5`.
	//   QualOut     — QUAL writer for `--out_format 2 / 5`. QUAL records
	//                 follow upstream `convertQualArrayToString`
	//                 (lines 2531-2546): two-character space-padded decimal
	//                 phred scores separated by spaces, wrapped every
	//                 QualLineWidth values; default 60 (LINE_WIDTH at
	//                 line 45).
	//   QualLineWidth — overrides the default 60-column QUAL wrap (mirrors
	//                 upstream `-line_width`); 0 keeps the default.
	NonIUPAC      bool
	SeqID         string
	SeqIDMap      io.Writer
	OutFormat     int
	FastaOut      io.Writer
	QualOut       io.Writer
	QualLineWidth int

	// Sliding-window quality trimming and range filters (mirrors
	// prinseq-lite.pl 0.20.4):
	//
	//   TrimQualWindow — `-trim_qual_window`: window size in bases used
	//                    when evaluating the quality trimming rule. The
	//                    upstream default is 1 (per-base trimming);
	//                    callers should default this to 1 whenever
	//                    TrimQualL/TrimQualR are in use to match upstream.
	//   TrimQualStep   — `-trim_qual_step`: how many bases the window
	//                    advances after a window satisfies the rule.
	//                    Upstream default 1.
	//   TrimQualType   — `-trim_qual_type`: how the window's score is
	//                    aggregated before comparison: "min" (default),
	//                    "max", "mean" or "sum" (prinseq-lite.pl:3219-3232).
	//   TrimQualRule   — `-trim_qual_rule`: comparison between the window
	//                    score and the TrimQualL/TrimQualR threshold:
	//                    "lt" (default; trim while score < threshold),
	//                    "gt" (trim while score > threshold) or "et"
	//                    (trim while score == threshold)
	//                    (prinseq-lite.pl:3236, 3276).
	//   TrimToLen      — `-trim_to_len`: hard-trim each read to at most
	//                    this many bases from the 5' end, applied AFTER all
	//                    other trimming and BEFORE the length/GC filters
	//                    (prinseq-lite.pl:3382-3385). Zero disables it.
	//   RangeLen       — `-range_len`: comma-separated "min-max" ranges; a
	//                    read passes only when its trimmed length lies
	//                    within every listed range (prinseq-lite.pl:3403,
	//                    checkRange at 2548). Empty disables it.
	//   RangeGC        — `-range_gc`: comma-separated "min-max" ranges
	//                    applied to the integer GC percentage
	//                    (prinseq-lite.pl:3458). Empty disables it.
	TrimQualWindow int
	TrimQualStep   int
	TrimQualType   string
	TrimQualRule   string
	TrimToLen      int
	RangeLen       string
	RangeGC        string

	// Sequence/header transforms and misc knobs (this PR; semantics from
	// prinseq-lite.pl 0.20.4):
	//
	//   SeqCase     — implements `-seq_case <upper|lower>` (force the
	//                 emitted sequence to upper- or lower-case;
	//                 prinseq-lite.pl:3664-3671). Empty leaves the case
	//                 untouched.
	//   DNARNA      — implements `-dna_rna <dna|rna>`: "dna" maps U->T /
	//                 u->t, "rna" maps T->U / t->u, preserving case
	//                 (prinseq-lite.pl:3672-3679). Empty disables.
	//   RmHeader    — implements `-rm_header`: drop the trailing FASTA/FASTQ
	//                 comment, leaving only the identifier
	//                 (prinseq-lite.pl:3651-3653, `$header = undef`).
	//   NoQualHeader — implements `-no_qual_header`: emit a bare "+" line
	//                 on FASTQ output instead of repeating the id
	//                 (prinseq-lite.pl:3686).
	//   SeqNum      — implements `-seq_num <int>`: keep only the first N
	//                 records that pass all other filters
	//                 (prinseq-lite.pl:3147-3150). Zero disables.
	//   ExactOnly   — accepts `-exact_only` for CLI compatibility. Upstream
	//                 uses it to restrict duplicate detection to exact (and,
	//                 with derep&4, reverse-complement) duplicates, skipping
	//                 the heuristic 5'/3' prefix/suffix dedup
	//                 (prinseq-lite.pl:3604-3628). This Go port never
	//                 implemented that prefix/suffix heuristic — its derep is
	//                 already exact-only — so the flag is a documented no-op
	//                 here: it is validated against the derep mode at the CLI
	//                 layer but does not alter Filter's behaviour.
	//   CustomParams — implements `-custom_params <string>` /
	//                 `-params <file>`: a list of complexity rules
	//                 evaluated against the upper-cased sequence
	//                 (prinseq-lite.pl:1069-1085, 3484-3503).
	SeqCase      string
	DNARNA       string
	RmHeader     bool
	NoQualHeader bool
	SeqNum       int
	ExactOnly    bool
	CustomParams []CustomParam
}

// CustomParam is one parsed `-custom_params` rule (prinseq-lite.pl:1069-1085).
// Bases is an upper-cased A/C/G/T/N pattern. When Percent is false the rule
// rejects a read whose upper-cased sequence contains Bases repeated Count
// times in a row. When Percent is true the rule rejects a read in which the
// fraction of Bases occurrences exceeds Count percent of the read length.
type CustomParam struct {
	Bases   string
	Count   int
	Percent bool
}

// ParseCustomParams parses a `-custom_params` string into a slice of
// CustomParam rules following upstream semantics (prinseq-lite.pl:1071-1085):
// rules are separated by ';'; each rule is "<bases> <count>" separated by
// whitespace. The bases token must consist solely of UPPER-CASE A/C/G/T/N:
// upstream validates with `$bases = ($tmp[0] =~ tr/ACGTN//); next if($bases <
// length($tmp[0]))`, and `tr/ACGTN//` only counts upper-case characters, so a
// token containing any lower-case letter fails the length check and is
// skipped. A '%' anywhere in the count marks a percentage rule. Tokens that
// fail validation are silently skipped, exactly as upstream's `next`
// statements do.
func ParseCustomParams(s string) []CustomParam {
	var rules []CustomParam
	for _, rule := range splitTrim(s, ';') {
		fields := splitWhitespace(rule)
		if len(fields) != 2 {
			continue
		}
		bases := fields[0]
		if !allUpperACGTN(bases) {
			continue
		}
		countStr := fields[1]
		percent := false
		if containsByte(countStr, '%') {
			countStr = removeByte(countStr, '%')
			percent = true
		}
		n, ok := parseUint(countStr)
		if !ok {
			continue
		}
		rules = append(rules, CustomParam{Bases: upperASCII(bases), Count: n, Percent: percent})
	}
	return rules
}

// Filter filters a FASTA/FASTQ file based on the given options
func Filter(reader io.Reader, writer io.Writer, isFastq bool, opts FilterOptions) error {
	if isFastq {
		return filterFastq(reader, writer, opts)
	}
	return filterFasta(reader, writer, opts)
}

// getQualityEncoding returns the appropriate quality encoding based on options
func getQualityEncoding(qualType string) fastq.QualityEncoding {
	if qualType == "illumina" {
		return fastq.Phred64
	}
	return fastq.Phred33 // Default to sanger
}

func filterFasta(reader io.Reader, writer io.Writer, opts FilterOptions) error {
	fastaReader := fasta.NewReader(reader)
	// We write good records via writePrinseqFastaWidth so that --line_width,
	// --seq_id comment preservation, and the case/DNA-RNA transforms all
	// apply; the bad stream keeps the plain fasta.Writer layout.
	fastaBW := bufio.NewWriter(writer)

	var badWriter *fasta.Writer
	if opts.OutBad != nil {
		badWriter = fasta.NewWriter(opts.OutBad, 80)
	}

	seenSeqs := make(map[string]int) // For duplicate tracking
	seqCount := 0                    // For --seq_id renumbering

	for {
		record, err := fastaReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		origDesc := record.Description
		origID, origComment := splitHeader(origDesc)
		if origID == "" {
			origID = firstToken(origDesc)
		}
		seq := string(record.Sequence)

		// Apply trimming
		seq, _ = trimSequence(seq, "", opts)

		// Enforce --seq_num: keep only the first N passing records.
		// Upstream checks `seq_num <= seqcount` before emitting
		// (prinseq-lite.pl:3147-3150); seqCount counts records already
		// written, so once we have emitted SeqNum the rest are rejected.
		if opts.SeqNum > 0 && seqCount >= opts.SeqNum {
			if badWriter != nil {
				record.Sequence = []byte(seq)
				if err := badWriter.Write(record); err != nil {
					return err
				}
			}
			continue
		}

		// Check for duplicates if derep is enabled
		if opts.Derep > 0 {
			if shouldFilterDuplicate(seq, seenSeqs, opts) {
				if badWriter != nil {
					record.Sequence = []byte(seq)
					if err := badWriter.Write(record); err != nil {
						return err
					}
				}
				continue
			}
		}

		// Apply filters
		if shouldFilterSequence(seq, "", opts) {
			// Write to bad output if enabled
			if badWriter != nil {
				record.Sequence = []byte(seq)
				if err := badWriter.Write(record); err != nil {
					return err
				}
			}
			continue
		}

		// This record passes; it counts towards --seq_num and gets a
		// fresh --seq_id counter value.
		seqCount++

		// Resolve the output identifier and trailing comment.
		id := origID
		comment := origComment
		if opts.SeqID != "" {
			id = renameDescription(opts.SeqID, seqCount)
			if err := writeSeqIDMapping(opts.SeqIDMap, origID, id); err != nil {
				return err
			}
		}
		// --rm_header drops the trailing comment (prinseq-lite.pl:3651).
		if opts.RmHeader {
			comment = ""
		}

		// Apply sequence transforms (case, DNA<->RNA).
		seq = applySeqCase(seq, opts.SeqCase)
		seq = applyDNARNA(seq, opts.DNARNA)

		record.Description = joinHeader(id, comment)
		record.Sequence = []byte(seq)

		// Write filtered record honouring --line_width.
		if err := writePrinseqFastaWidth(fastaBW, record.Description, record.Sequence, opts.QualLineWidth); err != nil {
			return err
		}
	}

	if err := fastaBW.Flush(); err != nil {
		return err
	}
	if badWriter != nil {
		return badWriter.Flush()
	}
	return nil
}

// firstToken returns the substring of s up to (but not including) the
// first ASCII whitespace byte, or s itself when none is present. It is
// the Go equivalent of upstream prinseq's `/^(\S+)/` extraction of the
// sequence identifier from a FASTA/FASTQ header.
func firstToken(s string) string {
	for i := 0; i < len(s); i++ {
		if s[i] == ' ' || s[i] == '\t' {
			return s[:i]
		}
	}
	return s
}

// splitHeader splits a full FASTA/FASTQ description line into the upstream
// `$sid` (the first whitespace-delimited token) and `$header` (everything
// after the first run of whitespace, with that leading whitespace removed).
// This mirrors upstream's `/^[\@\>](\S+)\s*(.*)$/` capture
// (prinseq-lite.pl:1660, 1678): the comment is the trailing remainder.
func splitHeader(desc string) (id, comment string) {
	i := 0
	for i < len(desc) && desc[i] != ' ' && desc[i] != '\t' {
		i++
	}
	id = desc[:i]
	for i < len(desc) && (desc[i] == ' ' || desc[i] == '\t') {
		i++
	}
	comment = desc[i:]
	return id, comment
}

// joinHeader recombines an identifier and trailing comment into the
// upstream output form `$sid.($header ? ' '.$header : ”)`
// (prinseq-lite.pl:3685, 3704). An empty comment yields just the id.
func joinHeader(id, comment string) string {
	if comment == "" {
		return id
	}
	return id + " " + comment
}

// applySeqCase forces the sequence case for `-seq_case`
// (prinseq-lite.pl:3664-3671). An unrecognised mode leaves seq untouched.
func applySeqCase(seq, mode string) string {
	switch mode {
	case "upper":
		return upperASCII(seq)
	case "lower":
		return lowerASCII(seq)
	default:
		return seq
	}
}

// applyDNARNA converts between DNA and RNA for `-dna_rna`
// (prinseq-lite.pl:3672-3679). "dna" maps U/u to T/t; "rna" maps T/t to
// U/u; case is preserved. Other modes leave seq untouched.
func applyDNARNA(seq, mode string) string {
	var from, to byte
	switch mode {
	case "dna":
		from, to = 'U', 'T'
	case "rna":
		from, to = 'T', 'U'
	default:
		return seq
	}
	fromLower := from | 0x20
	toLower := to | 0x20
	b := []byte(seq)
	for i := range b {
		switch b[i] {
		case from:
			b[i] = to
		case fromLower:
			b[i] = toLower
		}
	}
	return string(b)
}

// shouldFilterCustomParams reports whether the upper-cased sequence seqUpper
// matches any custom-params rejection rule (prinseq-lite.pl:3484-3503).
func shouldFilterCustomParams(seqUpper string, rules []CustomParam) bool {
	length := len(seqUpper)
	for _, p := range rules {
		if p.Percent {
			if length == 0 {
				continue
			}
			occ := countOccurrences(seqUpper, p.Bases)
			if float64(100*occ)/float64(length) > float64(p.Count) {
				return true
			}
		} else {
			if p.Count <= 0 {
				continue
			}
			if containsSubstr(seqUpper, repeatString(p.Bases, p.Count)) {
				return true
			}
		}
	}
	return false
}

// upperASCII upper-cases ASCII letters in s without touching other bytes.
func upperASCII(s string) string {
	b := []byte(s)
	for i := range b {
		if b[i] >= 'a' && b[i] <= 'z' {
			b[i] -= 0x20
		}
	}
	return string(b)
}

// lowerASCII lower-cases ASCII letters in s without touching other bytes.
func lowerASCII(s string) string {
	b := []byte(s)
	for i := range b {
		if b[i] >= 'A' && b[i] <= 'Z' {
			b[i] += 0x20
		}
	}
	return string(b)
}

// allUpperACGTN reports whether every byte of s is one of the UPPER-CASE
// bases A/C/G/T/N. It models upstream's custom_params validation, where
// `tr/ACGTN//` counts only upper-case characters, so any lower-case letter
// causes the rule to be skipped (prinseq-lite.pl:1075-1076).
func allUpperACGTN(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case 'A', 'C', 'G', 'T', 'N':
		default:
			return false
		}
	}
	return true
}

// splitTrim splits s on sep and trims surrounding ASCII whitespace from each
// field, mirroring upstream's `split(/\s*;\s*/, ...)`-style behaviour.
func splitTrim(s string, sep byte) []string {
	var out []string
	start := 0
	for i := 0; i <= len(s); i++ {
		if i == len(s) || s[i] == sep {
			out = append(out, trimSpace(s[start:i]))
			start = i + 1
		}
	}
	return out
}

// splitWhitespace splits s into non-empty runs separated by ASCII whitespace.
func splitWhitespace(s string) []string {
	var out []string
	start := -1
	for i := 0; i < len(s); i++ {
		if s[i] == ' ' || s[i] == '\t' {
			if start >= 0 {
				out = append(out, s[start:i])
				start = -1
			}
		} else if start < 0 {
			start = i
		}
	}
	if start >= 0 {
		out = append(out, s[start:])
	}
	return out
}

// trimSpace removes leading and trailing ASCII spaces and tabs from s.
func trimSpace(s string) string {
	i := 0
	for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
		i++
	}
	j := len(s)
	for j > i && (s[j-1] == ' ' || s[j-1] == '\t') {
		j--
	}
	return s[i:j]
}

// containsByte reports whether s contains byte c.
func containsByte(s string, c byte) bool {
	for i := 0; i < len(s); i++ {
		if s[i] == c {
			return true
		}
	}
	return false
}

// removeByte returns s with every occurrence of byte c removed.
func removeByte(s string, c byte) string {
	b := make([]byte, 0, len(s))
	for i := 0; i < len(s); i++ {
		if s[i] != c {
			b = append(b, s[i])
		}
	}
	return string(b)
}

// parseUint parses a non-negative base-10 integer. ok is false when s is
// empty or contains a non-digit, matching upstream's `^\d+$` guard.
func parseUint(s string) (int, bool) {
	if s == "" {
		return 0, false
	}
	n := 0
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			return 0, false
		}
		n = n*10 + int(s[i]-'0')
	}
	return n, true
}

// repeatString returns s concatenated n times (n <= 0 yields "").
func repeatString(s string, n int) string {
	if n <= 0 {
		return ""
	}
	b := make([]byte, 0, len(s)*n)
	for i := 0; i < n; i++ {
		b = append(b, s...)
	}
	return string(b)
}

// containsSubstr reports whether s contains sub. An empty sub never matches
// (upstream skips zero-repeat rules).
func containsSubstr(s, sub string) bool {
	if sub == "" {
		return false
	}
	if len(sub) > len(s) {
		return false
	}
	for i := 0; i+len(sub) <= len(s); i++ {
		if s[i:i+len(sub)] == sub {
			return true
		}
	}
	return false
}

// countOccurrences counts non-overlapping occurrences of sub in s, matching
// the behaviour of Perl's `$n++ while($seq =~ /$v/g)` for a literal pattern.
func countOccurrences(s, sub string) int {
	if sub == "" {
		return 0
	}
	count := 0
	for i := 0; i+len(sub) <= len(s); {
		if s[i:i+len(sub)] == sub {
			count++
			i += len(sub)
		} else {
			i++
		}
	}
	return count
}

// writePrinseqFastq writes one FASTQ record in the upstream PRINSEQ-lite
// layout, which repeats the sequence header on the "+" separator line
// (e.g. "+read1\n" rather than the bare "+\n" emitted by sickle / fastp).
// Upstream switches to a bare "+" only when -no_qual_header is supplied.
// Centralising the write here lets us swap encodings later without
// touching the filter loop.
func writePrinseqFastq(w *bufio.Writer, rec *fastq.Record) error {
	return writePrinseqFastqHeader(w, rec.Description, rec.Sequence, rec.Quality, false)
}

// writePrinseqFastqHeader writes one FASTQ record with explicit control over
// the "+" separator line. When noQualHeader is true the separator is the bare
// "+" upstream emits under -no_qual_header (prinseq-lite.pl:3686); otherwise
// it repeats the description, matching the default layout.
func writePrinseqFastqHeader(w *bufio.Writer, desc string, seq, qual []byte, noQualHeader bool) error {
	if _, err := fmt.Fprintf(w, "@%s\n", desc); err != nil {
		return err
	}
	if _, err := w.Write(seq); err != nil {
		return err
	}
	if err := w.WriteByte('\n'); err != nil {
		return err
	}
	if noQualHeader {
		if _, err := w.WriteString("+\n"); err != nil {
			return err
		}
	} else {
		if _, err := fmt.Fprintf(w, "+%s\n", desc); err != nil {
			return err
		}
	}
	if _, err := w.Write(qual); err != nil {
		return err
	}
	return w.WriteByte('\n')
}

// writePrinseqFasta writes one FASTA record under the upstream prinseq
// convention: ">desc\n<seq>\n" with no soft-wrap. The seq is already
// supplied as a single line; we deliberately do not split it here to
// match the byte-layout used by upstream when no `-line_width` is set
// (prinseq-lite.pl:3704-3708, where the substitution `s/(.{$linelen})/$1\n/g`
// is gated on `$linelen` being non-zero).
func writePrinseqFasta(w *bufio.Writer, desc string, seq []byte) error {
	return writePrinseqFastaWidth(w, desc, seq, 0)
}

// writePrinseqFastaWidth writes one FASTA record, wrapping the sequence every
// lineWidth bytes when lineWidth > 0. Upstream's wrap is the substitution
// `s/(.{$linelen})/$1\n/g` gated on a non-zero `$linelen`
// (prinseq-lite.pl:3699-3702); lineWidth <= 0 emits the sequence unwrapped.
func writePrinseqFastaWidth(w *bufio.Writer, desc string, seq []byte, lineWidth int) error {
	if _, err := fmt.Fprintf(w, ">%s\n", desc); err != nil {
		return err
	}
	if lineWidth > 0 {
		for off := 0; off < len(seq); off += lineWidth {
			end := off + lineWidth
			if end > len(seq) {
				end = len(seq)
			}
			if _, err := w.Write(seq[off:end]); err != nil {
				return err
			}
			if err := w.WriteByte('\n'); err != nil {
				return err
			}
		}
		if len(seq) == 0 {
			return w.WriteByte('\n')
		}
		return nil
	}
	if _, err := w.Write(seq); err != nil {
		return err
	}
	return w.WriteByte('\n')
}

// defaultQualLineWidth is the upstream LINE_WIDTH constant
// (prinseq-lite.pl:45). It controls how many phred values are written
// per line of a QUAL record.
const defaultQualLineWidth = 60

// writePrinseqQual writes one record to a QUAL file using upstream's
// convertQualArrayToString layout (prinseq-lite.pl:2531-2546):
//
//   - decode each ASCII quality byte to a decimal phred score under the
//     given offset (33 for Sanger / Phred+33, 64 for Illumina /
//     Phred+64);
//   - emit each score as a two-character field (leading space when
//     <10), followed by a single space, wrapping every lineWidth
//     values to a new line; and
//   - strip any trailing space/newline from the last token, then end
//     the record with a single "\n".
//
// linelen <= 0 selects the default of 60.
func writePrinseqQual(w *bufio.Writer, desc string, qual []byte, offset, lineWidth int) error {
	if lineWidth <= 0 {
		lineWidth = defaultQualLineWidth
	}
	if _, err := fmt.Fprintf(w, ">%s\n", desc); err != nil {
		return err
	}
	count := 0
	for i, q := range qual {
		score := int(q) - offset
		// Upstream caps to 93 on write (line 2487); on read it
		// accepts arbitrary Phred values. We don't cap here
		// because the values are read back as ASCII via the
		// FASTQ reader which has already validated the encoding.
		if score < 10 {
			if err := w.WriteByte(' '); err != nil {
				return err
			}
			if err := w.WriteByte(byte('0' + score)); err != nil {
				return err
			}
		} else {
			if _, err := fmt.Fprintf(w, "%d", score); err != nil {
				return err
			}
		}
		count++
		if count >= lineWidth && i != len(qual)-1 {
			if err := w.WriteByte('\n'); err != nil {
				return err
			}
			count = 0
		} else if i != len(qual)-1 {
			if err := w.WriteByte(' '); err != nil {
				return err
			}
		}
	}
	return w.WriteByte('\n')
}

// writeSeqIDMapping writes one row of the seq_id_mappings TSV in upstream's
// "<orig_id>\t<new_id>\n" format (prinseq-lite.pl:3646). Returns nil if
// the writer is nil so callers can unconditionally call it.
func writeSeqIDMapping(w io.Writer, origID, newID string) error {
	if w == nil {
		return nil
	}
	_, err := fmt.Fprintf(w, "%s\t%s\n", origID, newID)
	return err
}

// renameDescription returns the new IDENTIFIER used by `--seq_id`:
// "<SeqID><counter>". It replaces only the upstream `$sid` token; the
// caller re-attaches any trailing comment via joinHeader so that a record
// like "@read1 sample=A" becomes "@<prefix>N sample=A", matching upstream's
// `$sid.($header ? ' '.$header : ”)` (prinseq-lite.pl:3685-3691). The
// earlier divergence that dropped the comment has been resolved.
func renameDescription(prefix string, counter int) string {
	return fmt.Sprintf("%s%d", prefix, counter)
}

func filterFastq(reader io.Reader, writer io.Writer, opts FilterOptions) error {
	encoding := getQualityEncoding(opts.QualType)
	fastqReader := fastq.NewReader(reader, encoding)
	bw := bufio.NewWriter(writer)

	var bbw *bufio.Writer
	if opts.OutBad != nil {
		bbw = bufio.NewWriter(opts.OutBad)
	}

	// out_format == 2 (FASTA+QUAL) is the only mode where the FASTQ
	// input must be silently demoted to a FASTA stream on the primary
	// writer, with QUAL going to opts.QualOut. For modes 3/4/5 the
	// primary writer carries FASTQ as usual; for mode 1 it carries
	// FASTA. Mode 0 (unset) preserves the input format (FASTQ here).
	primaryFasta := opts.OutFormat == 1 || opts.OutFormat == 2
	emitFasta := opts.OutFormat == 4 || opts.OutFormat == 5
	emitQual := opts.OutFormat == 2 || opts.OutFormat == 5

	var fastaBW, qualBW *bufio.Writer
	if emitFasta && opts.FastaOut != nil {
		fastaBW = bufio.NewWriter(opts.FastaOut)
	}
	if emitQual && opts.QualOut != nil {
		qualBW = bufio.NewWriter(opts.QualOut)
	}

	qualOffset := phredOffset(opts.QualType)
	seenSeqs := make(map[string]int) // For duplicate tracking
	seqCount := 0                    // For --seq_id renumbering

	for {
		record, err := fastqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		origID, origComment := splitHeader(record.Description)
		if origID == "" {
			origID = record.ID
		}
		seq := string(record.Sequence)
		qual := string(record.Quality)

		// Apply trimming
		seq, qual = trimSequence(seq, qual, opts)

		// Enforce --seq_num (prinseq-lite.pl:3147-3150): once SeqNum good
		// records are emitted, reject the rest.
		if opts.SeqNum > 0 && seqCount >= opts.SeqNum {
			if bbw != nil {
				record.Sequence = []byte(seq)
				record.Quality = []byte(qual)
				if err := writePrinseqFastq(bbw, record); err != nil {
					return err
				}
			}
			continue
		}

		// Check for duplicates if derep is enabled
		if opts.Derep > 0 {
			if shouldFilterDuplicate(seq, seenSeqs, opts) {
				if bbw != nil {
					record.Sequence = []byte(seq)
					record.Quality = []byte(qual)
					if err := writePrinseqFastq(bbw, record); err != nil {
						return err
					}
				}
				continue
			}
		}

		// Apply filters
		if shouldFilterSequence(seq, qual, opts) {
			// Write to bad output if enabled
			if bbw != nil {
				record.Sequence = []byte(seq)
				record.Quality = []byte(qual)
				if err := writePrinseqFastq(bbw, record); err != nil {
					return err
				}
			}
			continue
		}

		// This record passes; advance the --seq_num / --seq_id counter.
		seqCount++

		// Resolve output identifier and trailing comment.
		id := origID
		comment := origComment
		if opts.SeqID != "" {
			id = renameDescription(opts.SeqID, seqCount)
			if err := writeSeqIDMapping(opts.SeqIDMap, origID, id); err != nil {
				return err
			}
		}
		if opts.RmHeader {
			comment = ""
		}

		// Apply sequence transforms (case, DNA<->RNA).
		seq = applySeqCase(seq, opts.SeqCase)
		seq = applyDNARNA(seq, opts.DNARNA)

		desc := joinHeader(id, comment)
		record.Description = desc
		record.Sequence = []byte(seq)
		record.Quality = []byte(qual)

		// Write primary output. Format depends on --out_format.
		if primaryFasta {
			if err := writePrinseqFastaWidth(bw, desc, record.Sequence, opts.QualLineWidth); err != nil {
				return err
			}
		} else {
			if err := writePrinseqFastqHeader(bw, desc, record.Sequence, record.Quality, opts.NoQualHeader); err != nil {
				return err
			}
		}

		// Optional secondary FASTA output (out_format 4/5).
		if fastaBW != nil {
			if err := writePrinseqFastaWidth(fastaBW, desc, record.Sequence, opts.QualLineWidth); err != nil {
				return err
			}
		}

		// Optional QUAL output (out_format 2/5).
		if qualBW != nil {
			if err := writePrinseqQual(qualBW, desc, record.Quality, qualOffset, opts.QualLineWidth); err != nil {
				return err
			}
		}
	}

	if err := bw.Flush(); err != nil {
		return err
	}
	if fastaBW != nil {
		if err := fastaBW.Flush(); err != nil {
			return err
		}
	}
	if qualBW != nil {
		if err := qualBW.Flush(); err != nil {
			return err
		}
	}
	if bbw != nil {
		return bbw.Flush()
	}
	return nil
}

func shouldFilterSequence(seq, qual string, opts FilterOptions) bool {
	seqLen := len(seq)

	// Zero-length reads are always dropped. Upstream marks any read that
	// trimming reduced to length 0 as bad via its `zero_length` filter
	// (prinseq-lite.pl:3389-3392), independent of any min_len setting.
	if seqLen == 0 {
		return true
	}

	// Length filters
	if opts.MinLen > 0 && seqLen < opts.MinLen {
		return true
	}
	if opts.MaxLen > 0 && seqLen > opts.MaxLen {
		return true
	}
	// range_len: reject when the length is outside the requested
	// range(s) (prinseq-lite.pl:3403). Upstream applies this after the
	// min_len/max_len checks.
	if opts.RangeLen != "" && !checkRange(opts.RangeLen, seqLen) {
		return true
	}

	// GC content filter
	gcCount := 0
	nCount := 0
	for _, base := range seq {
		switch base {
		case 'G', 'C', 'g', 'c':
			gcCount++
		case 'N', 'n':
			nCount++
		}
	}

	if seqLen > 0 {
		gcContent := float64(gcCount) / float64(seqLen) * 100.0
		if opts.MinGC > 0 && gcContent < opts.MinGC {
			return true
		}
		if opts.MaxGC > 0 && gcContent > opts.MaxGC {
			return true
		}
		// range_gc: reject when the integer GC percentage is outside
		// the requested range(s). Upstream truncates the percentage to
		// an integer via sprintf("%d", gc*100/len) before the
		// checkRange comparison (prinseq-lite.pl:3450-3460), so we
		// truncate toward zero here too rather than reusing the float
		// gcContent above.
		if opts.RangeGC != "" {
			gcInt := (gcCount * 100) / seqLen
			if !checkRange(opts.RangeGC, gcInt) {
				return true
			}
		}
	}

	// N content filters
	if opts.MaxNsN > 0 && nCount > opts.MaxNsN {
		return true
	}
	if opts.MaxNsP > 0 && seqLen > 0 {
		nPercent := float64(nCount) / float64(seqLen) * 100.0
		if nPercent > opts.MaxNsP {
			return true
		}
	}

	// Non-IUPAC strict filter. Upstream's check is
	// `uc($seq) =~ /[^ACGTN]/o` (prinseq-lite.pl:3478), i.e. any base
	// other than A/C/G/T/N (case-insensitive) marks the read for
	// rejection. We match that semantics here.
	if opts.NonIUPAC {
		for _, base := range seq {
			switch base {
			case 'A', 'C', 'G', 'T', 'N',
				'a', 'c', 'g', 't', 'n':
				// ok
			default:
				return true
			}
		}
	}

	// Quality filters (only for FASTQ). Honour the per-options Phred
	// offset so that `-phred64` inputs are decoded against ASCII 64.
	// Before this fix the mean was always computed under Phred+33,
	// silently mis-classifying Phred+64 reads as much higher quality
	// than they actually were.
	if qual != "" {
		avgQual := calculateAvgQualityScoreWithOffset(qual, phredOffset(opts.QualType))
		if opts.MinQualMean > 0 && avgQual < opts.MinQualMean {
			return true
		}
		if opts.MaxQualMean > 0 && avgQual > opts.MaxQualMean {
			return true
		}
	}

	// Custom-params complexity rules (upstream `-custom_params` /
	// `-params`; prinseq-lite.pl:3484-3503). Rules are matched against the
	// upper-cased sequence.
	if len(opts.CustomParams) > 0 {
		if shouldFilterCustomParams(upperASCII(seq), opts.CustomParams) {
			return true
		}
	}

	// Complexity filter
	if opts.LcMethod != "" && opts.LcThreshold > 0 {
		if opts.LcMethod == "dust" {
			score := calculateDustScore(seq)
			if score > opts.LcThreshold {
				return true
			}
		} else if opts.LcMethod == "entropy" {
			score := calculateEntropy(seq)
			if score < opts.LcThreshold {
				return true
			}
		}
	}

	return false
}

// calculateDustScore computes DUST score for low-complexity filtering
// Lower scores indicate higher complexity; higher scores indicate low complexity
func calculateDustScore(seq string) float64 {
	if len(seq) < 3 {
		return 0.0
	}

	// Count triplet frequencies
	triplets := make(map[string]int)
	for i := 0; i <= len(seq)-3; i++ {
		triplet := seq[i : i+3]
		triplets[triplet]++
	}

	// Calculate DUST score as sum of (count * (count-1) / 2) for each triplet
	score := 0.0
	for _, count := range triplets {
		if count > 1 {
			score += float64(count * (count - 1) / 2)
		}
	}

	// Normalize by sequence length
	if len(seq) > 0 {
		score = score / float64(len(seq)-2) * 10.0
	}

	return score
}

// calculateEntropy computes Shannon entropy for complexity filtering
// Higher scores indicate higher complexity; lower scores indicate low complexity
func calculateEntropy(seq string) float64 {
	if len(seq) == 0 {
		return 0.0
	}

	// Count base frequencies
	counts := make(map[rune]int)
	for _, base := range seq {
		counts[base]++
	}

	// Calculate Shannon entropy
	entropy := 0.0
	length := float64(len(seq))
	for _, count := range counts {
		if count > 0 {
			p := float64(count) / length
			entropy -= p * math.Log2(p)
		}
	}

	// Normalize to percentage (0-100)
	// Maximum entropy for DNA is log2(4) = 2.0, so normalize to 100
	return (entropy / 2.0) * 100.0
}

func shouldFilterDuplicate(seq string, seenSeqs map[string]int, opts FilterOptions) bool {
	if opts.Derep == 0 {
		return false
	}

	// Mode 1: exact duplicate
	if opts.Derep&1 != 0 {
		seenSeqs[seq]++
		if seenSeqs[seq] >= opts.DerepMin {
			return true
		}
	}

	// Mode 4: reverse complement exact duplicate
	if opts.Derep&4 != 0 {
		revComp := reverseComplement(seq)
		seenSeqs[revComp]++
		if seenSeqs[revComp] >= opts.DerepMin {
			return true
		}
	}

	return false
}

// checkRange reports whether val lies within the comma-separated set of
// "min-max" ranges, faithfully reproducing prinseq-lite.pl's checkRange
// (lines 2548-2557). Note the upstream semantics: the value must lie
// within EVERY listed range — the helper returns false as soon as val is
// below the lower bound or above the upper bound of any range. With
// disjoint ranges (e.g. "10-20,40-50") this means nothing passes; we
// mirror that behaviour exactly rather than "fixing" it to an OR.
//
// Malformed range tokens are treated leniently, matching Perl's loose
// numeric coercion: a token without a "-" parses its single number as the
// lower bound with an implicit upper bound of 0 (so any positive value
// fails), and non-numeric fields coerce to 0.
func checkRange(ranges string, val int) bool {
	for _, r := range splitOn(ranges, ',') {
		parts := splitOn(r, '-')
		lo := atoiOrZero(parts[0])
		hi := 0
		if len(parts) > 1 {
			hi = atoiOrZero(parts[1])
		}
		if val < lo || val > hi {
			return false
		}
	}
	return true
}

// splitOn splits s on the single byte sep, returning at least one element
// (the empty string when s is empty). It avoids pulling in strings.Split's
// behaviour differences for the simple cases checkRange needs and keeps the
// Perl `split(/\,/)` / `split(/\-/)` semantics of dropping nothing.
func splitOn(s string, sep byte) []string {
	var out []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == sep {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	out = append(out, s[start:])
	return out
}

// atoiOrZero parses a leading signed integer from s (after trimming ASCII
// spaces), returning 0 when no digits are present. This matches Perl's
// numeric coercion of a string in a comparison context, which is how
// upstream's split fields are used.
func atoiOrZero(s string) int {
	i := 0
	for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
		i++
	}
	neg := false
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		neg = s[i] == '-'
		i++
	}
	n := 0
	got := false
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		n = n*10 + int(s[i]-'0')
		i++
		got = true
	}
	if !got {
		return 0
	}
	if neg {
		return -n
	}
	return n
}

func reverseComplement(seq string) string {
	complement := map[rune]rune{
		'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
		'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
		'N': 'N', 'n': 'n',
	}

	result := make([]rune, len(seq))
	for i, base := range seq {
		if comp, ok := complement[base]; ok {
			result[len(seq)-1-i] = comp
		} else {
			result[len(seq)-1-i] = base
		}
	}
	return string(result)
}

// FilterPaired filters paired-end FASTA/FASTQ files
func FilterPaired(reader1, reader2 io.Reader, writer1, writer2 io.Writer, isFastq bool, opts FilterOptions) error {
	if isFastq {
		return filterPairedFastq(reader1, reader2, writer1, writer2, opts)
	}
	return filterPairedFasta(reader1, reader2, writer1, writer2, opts)
}

func filterPairedFasta(reader1, reader2 io.Reader, writer1, writer2 io.Writer, opts FilterOptions) error {
	fastaReader1 := fasta.NewReader(reader1)
	fastaReader2 := fasta.NewReader(reader2)
	fastaWriter1 := fasta.NewWriter(writer1, 80)
	fastaWriter2 := fasta.NewWriter(writer2, 80)

	seenSeqs := make(map[string]int)

	for {
		record1, err1 := fastaReader1.Read()
		record2, err2 := fastaReader2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF || err2 == io.EOF {
			return fmt.Errorf("paired files have different number of sequences")
		}
		if err1 != nil {
			return err1
		}
		if err2 != nil {
			return err2
		}

		seq1 := string(record1.Sequence)
		seq2 := string(record2.Sequence)

		// Apply trimming
		seq1, _ = trimSequence(seq1, "", opts)
		seq2, _ = trimSequence(seq2, "", opts)

		// Check duplicates (consider both reads together)
		if opts.Derep > 0 {
			combinedSeq := seq1 + "|" + seq2
			if shouldFilterDuplicate(combinedSeq, seenSeqs, opts) {
				continue
			}
		}

		// Filter: if either read fails, both are filtered
		if shouldFilterSequence(seq1, "", opts) || shouldFilterSequence(seq2, "", opts) {
			continue
		}

		// Update records
		record1.Sequence = []byte(seq1)
		record2.Sequence = []byte(seq2)

		// Write both records
		if err := fastaWriter1.Write(record1); err != nil {
			return err
		}
		if err := fastaWriter2.Write(record2); err != nil {
			return err
		}
	}

	if err := fastaWriter1.Flush(); err != nil {
		return err
	}
	return fastaWriter2.Flush()
}

func filterPairedFastq(reader1, reader2 io.Reader, writer1, writer2 io.Writer, opts FilterOptions) error {
	encoding := getQualityEncoding(opts.QualType)
	fastqReader1 := fastq.NewReader(reader1, encoding)
	fastqReader2 := fastq.NewReader(reader2, encoding)
	bw1 := bufio.NewWriter(writer1)
	bw2 := bufio.NewWriter(writer2)

	seenSeqs := make(map[string]int)

	for {
		record1, err1 := fastqReader1.Read()
		record2, err2 := fastqReader2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF || err2 == io.EOF {
			return fmt.Errorf("paired files have different number of sequences")
		}
		if err1 != nil {
			return err1
		}
		if err2 != nil {
			return err2
		}

		seq1 := string(record1.Sequence)
		qual1 := string(record1.Quality)
		seq2 := string(record2.Sequence)
		qual2 := string(record2.Quality)

		// Apply trimming
		seq1, qual1 = trimSequence(seq1, qual1, opts)
		seq2, qual2 = trimSequence(seq2, qual2, opts)

		// Check duplicates (consider both reads together)
		if opts.Derep > 0 {
			combinedSeq := seq1 + "|" + seq2
			if shouldFilterDuplicate(combinedSeq, seenSeqs, opts) {
				continue
			}
		}

		// Filter: if either read fails, both are filtered
		if shouldFilterSequence(seq1, qual1, opts) || shouldFilterSequence(seq2, qual2, opts) {
			continue
		}

		// Update records
		record1.Sequence = []byte(seq1)
		record1.Quality = []byte(qual1)
		record2.Sequence = []byte(seq2)
		record2.Quality = []byte(qual2)

		// Write both records in upstream-compatible FASTQ format.
		if err := writePrinseqFastq(bw1, record1); err != nil {
			return err
		}
		if err := writePrinseqFastq(bw2, record2); err != nil {
			return err
		}
	}

	if err := bw1.Flush(); err != nil {
		return err
	}
	return bw2.Flush()
}

// CalculateEnhancedStats computes enhanced statistics including distributions and dinucleotides
func CalculateEnhancedStats(reader io.Reader, isFastq bool) (*Stats, error) {
	return CalculateEnhancedStatsWithEncoding(reader, isFastq, "sanger")
}

// CalculateEnhancedStatsWithEncoding computes enhanced statistics with specific quality encoding
func CalculateEnhancedStatsWithEncoding(reader io.Reader, isFastq bool, qualType string) (*Stats, error) {
	stats := &Stats{
		MinLength:           -1,
		LengthDistribution:  make(map[int]int),
		QualityDistribution: make(map[int]int),
		BaseComposition:     make(map[string]int),
		Dinucleotides:       make(map[string]int),
	}

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 1024*1024), 1024*1024) // 1MB buffer

	if isFastq {
		offset := 33
		if qualType == "illumina" {
			offset = 64
		}
		return calculateEnhancedFastqStats(scanner, stats, offset)
	}
	return calculateEnhancedFastaStats(scanner, stats)
}

func calculateEnhancedFastaStats(scanner *bufio.Scanner, stats *Stats) (*Stats, error) {
	var currentSeq string
	gcCount := 0

	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			// Process previous sequence
			if currentSeq != "" {
				processEnhancedSequence(currentSeq, "", &gcCount, stats)
				currentSeq = ""
			}
		} else {
			currentSeq += line
		}
	}

	// Process last sequence
	if currentSeq != "" {
		processEnhancedSequence(currentSeq, "", &gcCount, stats)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
	}

	return stats, nil
}

func calculateEnhancedFastqStats(scanner *bufio.Scanner, stats *Stats, offset int) (*Stats, error) {
	gcCount := 0
	totalQuality := 0.0
	lineNum := 0
	maxLen := 0

	var seq string
	var qual string

	for scanner.Scan() {
		line := scanner.Text()
		mod := lineNum % 4

		switch mod {
		case 0: // Header line
			if len(line) == 0 || line[0] != '@' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 1: // Sequence line
			seq = line
		case 2: // Plus line
			if len(line) == 0 || line[0] != '+' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 3: // Quality line
			qual = line
			if len(qual) != len(seq) {
				return nil, fmt.Errorf("quality length (%d) doesn't match sequence length (%d) at line %d",
					len(qual), len(seq), lineNum+1)
			}
			// Process the complete FASTQ record
			processEnhancedSequence(seq, qual, &gcCount, stats)
			avgQual := calculateAvgQualityScoreWithOffset(qual, offset)
			totalQuality += avgQual

			// Update quality distribution
			qualInt := int(avgQual)
			stats.QualityDistribution[qualInt]++

			// Update positional quality
			if len(seq) > maxLen {
				maxLen = len(seq)
				// Extend positional quality array
				for len(stats.PositionalQuality) < maxLen {
					stats.PositionalQuality = append(stats.PositionalQuality, 0)
				}
			}
			for i, q := range qual {
				qualVal := float64(int(q) - offset)
				if i < len(stats.PositionalQuality) {
					stats.PositionalQuality[i] += qualVal
				}
			}

			seq = ""
			qual = ""
		}

		lineNum++
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Check if we have incomplete records
	if lineNum%4 != 0 {
		return nil, fmt.Errorf("incomplete FASTQ record")
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
		stats.AvgQuality = totalQuality / float64(stats.NumReads)

		// Average positional quality
		for i := range stats.PositionalQuality {
			stats.PositionalQuality[i] /= float64(stats.NumReads)
		}
	}

	return stats, nil
}

func processEnhancedSequence(seq, qual string, gcCount *int, stats *Stats) {
	seqLen := len(seq)
	stats.NumReads++
	stats.TotalBases += seqLen

	if stats.MinLength == -1 || seqLen < stats.MinLength {
		stats.MinLength = seqLen
	}
	if seqLen > stats.MaxLength {
		stats.MaxLength = seqLen
	}

	// Update length distribution
	stats.LengthDistribution[seqLen]++

	// Count bases, GC, Ns, and dinucleotides
	var prevBase byte
	for i, base := range seq {
		// Base composition
		baseStr := string(base)
		stats.BaseComposition[baseStr]++

		switch base {
		case 'G', 'C', 'g', 'c':
			*gcCount++
		case 'N', 'n':
			stats.NumNs++
		}

		// Dinucleotides
		if i > 0 {
			dinuc := string([]byte{prevBase, byte(base)})
			stats.Dinucleotides[dinuc]++
		}
		prevBase = byte(base)
	}
}
