// Package sickle provides quality-based trimming for FASTQ files using a sliding window approach.
package sickle

import (
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TrimOptions contains parameters for quality-based trimming.
type TrimOptions struct {
	QualThreshold   int  // Minimum quality score threshold
	LengthThreshold int  // Minimum length threshold after trimming
	NoFivePrime     bool // Don't trim 5' end
	TruncateN       bool // Truncate at first N
	// WindowSize is the size of the sliding window used for quality
	// assessment. A value <= 0 selects upstream sickle's dynamic per-read
	// window of int(0.1*read_length) (falling back to the full read length
	// when that rounds to 0). A positive value is a Go-port extension that
	// pins a fixed window for every read; upstream sickle has no -w flag.
	WindowSize  int
	Progress    bool // Show progress reporting
	Recalibrate bool // Recalibrate quality scores
}

// DefaultTrimOptions returns default trimming options. WindowSize defaults to 0
// (dynamic), so the sliding window is computed per read as int(0.1*read_length),
// matching upstream sickle exactly. Upstream has no -w flag and always uses this
// dynamic window; the fixed-window override (WindowSize > 0) is a Go extension.
func DefaultTrimOptions() TrimOptions {
	return TrimOptions{
		QualThreshold:   20,
		LengthThreshold: 20,
		NoFivePrime:     false,
		TruncateN:       false,
		WindowSize:      0,
	}
}

// TrimStats tracks trimming statistics.
type TrimStats struct {
	TotalReads     int
	TrimmedReads   int
	DiscardedReads int
	TotalBases     int64
	TrimmedBases   int64
}

// TrimSingleEnd trims a single-end FASTQ file based on quality scores.
func TrimSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	stats := &TrimStats{}

	// Progress reporting setup
	var progressCounter int
	const progressInterval = 10000

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return stats, fmt.Errorf("error reading FASTQ: %w", err)
		}

		stats.TotalReads++
		stats.TotalBases += int64(len(record.Sequence))

		// Apply recalibration if requested
		if opts.Recalibrate {
			record = recalibrateRecord(record, encoding)
		}

		// Apply trimming (upstream-faithful sliding window; needs the
		// encoding so the Phred offset is right for sanger vs illumina).
		trimmed := trimRecord(record, opts, encoding)

		// Check if read passes length threshold
		if len(trimmed.Sequence) >= opts.LengthThreshold {
			if err := writer.Write(trimmed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			if len(trimmed.Sequence) < len(record.Sequence) {
				stats.TrimmedReads++
				stats.TrimmedBases += int64(len(record.Sequence) - len(trimmed.Sequence))
			}
		} else {
			stats.DiscardedReads++
		}

		// Progress reporting
		if opts.Progress && stats.TotalReads%progressInterval == 0 {
			progressCounter++
			fmt.Fprintf(os.Stderr, "\rProcessed %d reads...", stats.TotalReads)
		}
	}

	// Clear progress line
	if opts.Progress {
		fmt.Fprintf(os.Stderr, "\r")
	}

	// Flush writer
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}

	return stats, nil
}

// TrimPairedEnd trims paired-end FASTQ files, maintaining synchronization.
func TrimPairedEnd(input1, input2 io.Reader, output1, output2, outputSingle io.Writer,
	encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {

	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)

	var writerSingle *fastq.Writer
	if outputSingle != nil {
		writerSingle = fastq.NewWriter(outputSingle, encoding)
	}

	stats := &TrimStats{}

	// Progress reporting setup
	const progressInterval = 10000

	for {
		record1, err1 := reader1.Read()
		record2, err2 := reader2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 != nil && err1 != io.EOF {
			return stats, fmt.Errorf("error reading first FASTQ: %w", err1)
		}
		if err2 != nil && err2 != io.EOF {
			return stats, fmt.Errorf("error reading second FASTQ: %w", err2)
		}
		if (err1 == io.EOF) != (err2 == io.EOF) {
			return stats, fmt.Errorf("paired-end files have different number of reads")
		}

		stats.TotalReads += 2
		stats.TotalBases += int64(len(record1.Sequence) + len(record2.Sequence))

		// Apply recalibration if requested
		if opts.Recalibrate {
			record1 = recalibrateRecord(record1, encoding)
			record2 = recalibrateRecord(record2, encoding)
		}

		// Trim both reads (upstream-faithful sliding window).
		trimmed1 := trimRecord(record1, opts, encoding)
		trimmed2 := trimRecord(record2, opts, encoding)

		pass1 := len(trimmed1.Sequence) >= opts.LengthThreshold
		pass2 := len(trimmed2.Sequence) >= opts.LengthThreshold

		// Both reads pass - write to paired output
		if pass1 && pass2 {
			if err := writer1.Write(trimmed1); err != nil {
				return stats, fmt.Errorf("error writing first FASTQ: %w", err)
			}
			if err := writer2.Write(trimmed2); err != nil {
				return stats, fmt.Errorf("error writing second FASTQ: %w", err)
			}

			if len(trimmed1.Sequence) < len(record1.Sequence) {
				stats.TrimmedReads++
				stats.TrimmedBases += int64(len(record1.Sequence) - len(trimmed1.Sequence))
			}
			if len(trimmed2.Sequence) < len(record2.Sequence) {
				stats.TrimmedReads++
				stats.TrimmedBases += int64(len(record2.Sequence) - len(trimmed2.Sequence))
			}
		} else if writerSingle != nil {
			// One read passes - write to single output if available
			if pass1 {
				if err := writerSingle.Write(trimmed1); err != nil {
					return stats, fmt.Errorf("error writing single FASTQ: %w", err)
				}
				if len(trimmed1.Sequence) < len(record1.Sequence) {
					stats.TrimmedReads++
					stats.TrimmedBases += int64(len(record1.Sequence) - len(trimmed1.Sequence))
				}
			}
			if pass2 {
				if err := writerSingle.Write(trimmed2); err != nil {
					return stats, fmt.Errorf("error writing single FASTQ: %w", err)
				}
				if len(trimmed2.Sequence) < len(record2.Sequence) {
					stats.TrimmedReads++
					stats.TrimmedBases += int64(len(record2.Sequence) - len(trimmed2.Sequence))
				}
			}
			if !pass1 {
				stats.DiscardedReads++
			}
			if !pass2 {
				stats.DiscardedReads++
			}
		} else {
			// Discard both if either fails and no single output
			stats.DiscardedReads += 2
		}

		// Progress reporting
		if opts.Progress && stats.TotalReads%progressInterval == 0 {
			fmt.Fprintf(os.Stderr, "\rProcessed %d reads...", stats.TotalReads)
		}
	}

	// Clear progress line
	if opts.Progress {
		fmt.Fprintf(os.Stderr, "\r")
	}

	// Flush all writers
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

	return stats, nil
}

// qualityOffset returns the ASCII offset (Q = ASCII byte - offset) for the
// given fastq.QualityEncoding. Upstream sickle interprets quality bytes as
// decoded integer Phred values, so the offset must match its per-encoding
// table exactly: 33 for sanger, 64 for illumina/solexa.
func qualityOffset(enc fastq.QualityEncoding) int {
	if enc == fastq.Phred64 {
		return 64
	}
	return 33
}

// resolveWindowSize returns the sliding-window length to use for a read of
// seqLen bases, given the user-requested optWindow.
//
// When optWindow <= 0 (the default), it reproduces upstream sickle's dynamic
// per-read rule from reference_code/sickle/src/sliding.c exactly:
//
//	window_size = (int)(0.1 * seq.l);
//	if (window_size == 0) window_size = seq.l;
//
// i.e. the window is one tenth of the read length truncated toward zero, and
// falls back to the full read length when that truncates to 0 (reads < 10 bp).
// Upstream never clamps this dynamic value, because int(0.1*L) can never exceed
// L. Upstream has no -w flag, so this dynamic window is upstream's only mode.
//
// When optWindow > 0 (a Go-port extension), that fixed window is used for every
// read, clamped down to seqLen so a window longer than the read still works.
func resolveWindowSize(seqLen, optWindow int) int {
	if optWindow > 0 {
		if optWindow > seqLen {
			return seqLen
		}
		return optWindow
	}
	w := int(0.1 * float64(seqLen))
	if w == 0 {
		w = seqLen
	}
	return w
}

// trimRecord applies upstream-sickle-faithful sliding-window quality trimming
// to a FASTQ record. The algorithm mirrors reference_code/sickle/src/sliding.c:
//
//  1. Window size is resolved by resolveWindowSize: by default (opts.WindowSize
//     <= 0) it is the per-read dynamic int(0.1*len(seq)), falling back to the
//     full read length when that rounds to 0 (read < 10 bp). Callers may pin a
//     fixed window by setting opts.WindowSize > 0 (a Go-port extension —
//     upstream has no -w flag and always uses the dynamic window).
//  2. Walk the window left-to-right. Per position:
//     a. If no 5' cut yet and the window's *average* decoded quality is
//     >= QualThreshold, set five_prime_cut to the first index within the
//     window whose decoded quality is >= QualThreshold.
//     b. If a 5' cut has been found (or NoFivePrime is set) AND either the
//     window's average drops below threshold OR the window now extends
//     past the read end, set three_prime_cut to the first index within
//     the window whose decoded quality is < threshold, then break.
//  3. If TruncateN is set and seq contains an N (case-insensitive), override
//     three_prime_cut with the index of that first N (upstream's
//     `strstr(seq, "N") || strstr(seq, "n")` does this unconditionally
//     *after* the sliding pass).
//  4. If no 5' cut was found (and NoFivePrime is off) OR
//     three_prime_cut - five_prime_cut < LengthThreshold, the read is
//     discarded — represented here by returning an empty sequence so the
//     calling loop can drop it (matches upstream's three_prime_cut = -1
//     sentinel).
//
// qualityEnc selects the ASCII offset for decoding qual bytes (33 vs 64).
func trimRecord(record *fastq.Record, opts TrimOptions, qualityEnc fastq.QualityEncoding) *fastq.Record {
	seq := record.Sequence
	qual := record.Quality

	emptyOut := func() *fastq.Record {
		return &fastq.Record{
			ID:          record.ID,
			Description: record.Description,
			Sequence:    []byte{},
			Quality:     []byte{},
		}
	}

	if len(seq) == 0 {
		return emptyOut()
	}

	// Upstream discards reads shorter than the length threshold up front.
	if len(seq) < opts.LengthThreshold {
		return emptyOut()
	}

	offset := qualityOffset(qualityEnc)

	// Resolve the sliding-window size for this read.
	windowSize := resolveWindowSize(len(seq), opts.WindowSize)

	threshold := opts.QualThreshold
	threePrimeCut := len(seq)
	fivePrimeCut := 0
	foundFivePrime := false

	// Seed the window.
	windowTotal := 0
	for i := 0; i < windowSize; i++ {
		windowTotal += int(qual[i]) - offset
	}

	windowStart := 0
	endLoop := len(qual) - windowSize
	for i := 0; i <= endLoop; i++ {
		windowAvg := float64(windowTotal) / float64(windowSize)

		// 5' cut: first window whose average >= threshold.
		if !opts.NoFivePrime && !foundFivePrime && windowAvg >= float64(threshold) {
			for j := windowStart; j < windowStart+windowSize && j < len(qual); j++ {
				if int(qual[j])-offset >= threshold {
					fivePrimeCut = j
					break
				}
			}
			foundFivePrime = true
		}

		// 3' cut: window average dipped below threshold, OR we've reached
		// the last window — but only after the 5' cut is settled (or
		// 5'-trimming is disabled).
		if (windowAvg < float64(threshold) || windowStart+windowSize > len(qual)) && (foundFivePrime || opts.NoFivePrime) {
			for j := windowStart; j < windowStart+windowSize && j < len(qual); j++ {
				if int(qual[j])-offset < threshold {
					threePrimeCut = j
					break
				}
			}
			break
		}

		// Slide one position: subtract the byte leaving the window, add
		// the byte entering it (if one exists).
		windowTotal -= int(qual[windowStart]) - offset
		if windowStart+windowSize < len(qual) {
			windowTotal += int(qual[windowStart+windowSize]) - offset
		}
		windowStart++
	}

	// trunc-N override (case-insensitive): if any N appears, force 3' cut
	// to its position, mirroring upstream's behavior.
	if opts.TruncateN {
		for i := 0; i < len(seq); i++ {
			if seq[i] == 'N' || seq[i] == 'n' {
				threePrimeCut = i
				break
			}
		}
	}

	// Discard rules: no 5' cut found (when 5'-trimming is on), or the
	// kept span is shorter than the length threshold.
	if (!foundFivePrime && !opts.NoFivePrime) || threePrimeCut-fivePrimeCut < opts.LengthThreshold {
		return emptyOut()
	}

	if fivePrimeCut < 0 {
		fivePrimeCut = 0
	}
	if threePrimeCut > len(seq) {
		threePrimeCut = len(seq)
	}
	if fivePrimeCut >= threePrimeCut {
		return emptyOut()
	}

	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    seq[fivePrimeCut:threePrimeCut],
		Quality:     qual[fivePrimeCut:threePrimeCut],
	}
}

// recalibrateRecord recalibrates quality scores using empirical base quality score recalibration.
// This is a simplified version that adjusts quality scores based on sequence context.
func recalibrateRecord(record *fastq.Record, encoding fastq.QualityEncoding) *fastq.Record {
	if len(record.Quality) == 0 {
		return record
	}

	// Get quality offset
	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}

	// Create a copy of the quality scores
	newQuality := make([]byte, len(record.Quality))
	copy(newQuality, record.Quality)

	// Simple recalibration: adjust quality scores based on position and context
	// This is a simplified algorithm - real recalibration would use machine learning
	for i := 0; i < len(newQuality); i++ {
		currentQual := int(newQuality[i]) - offset

		// Position-based adjustment: quality tends to degrade toward read ends
		positionFactor := 1.0
		readLength := len(newQuality)
		if i < readLength/10 || i > 9*readLength/10 {
			// First and last 10% of read: slight quality penalty
			positionFactor = 0.95
		}

		// Context-based adjustment: homopolymers are more error-prone
		if i > 0 && i < len(record.Sequence)-1 {
			prevBase := record.Sequence[i-1]
			currBase := record.Sequence[i]
			nextBase := record.Sequence[i+1]

			// Penalize homopolymer runs
			if prevBase == currBase || currBase == nextBase {
				positionFactor *= 0.95
			}
		}

		// Apply adjustments
		adjustedQual := int(float64(currentQual) * positionFactor)

		// Clamp to valid range
		if adjustedQual < 0 {
			adjustedQual = 0
		}
		maxQual := 93
		if encoding == fastq.Phred64 {
			maxQual = 62
		}
		if adjustedQual > maxQual {
			adjustedQual = maxQual
		}

		newQuality[i] = byte(adjustedQual + offset)
	}

	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence,
		Quality:     newQuality,
	}
}
