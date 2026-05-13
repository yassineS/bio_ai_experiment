// Encoding detection for FASTQ quality scores.
//
// Sickle accepts three named encodings (sanger/Phred+33, illumina/Phred+64,
// solexa/Solexa+64) but the underlying fastq library only models the two Phred
// variants. DetectEncoding reads up to the first ~10000 quality bytes from a
// FASTQ stream WITHOUT consuming the underlying io.Reader (it uses
// bufio.Reader.Peek), so the same stream can be passed to fastq.NewReader
// afterwards and processed in full.
package sickle

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// DetectionResult describes what DetectEncoding inferred from a FASTQ stream.
// Encoding is the matching fastq.QualityEncoding (always one of fastq.Phred33
// or fastq.Phred64 — the fastq library does not model Solexa separately, so
// Solexa+64 is reported as fastq.Phred64 with Name set to "solexa"). Name is
// the human-readable label ("sanger", "illumina", "solexa"). Ambiguous is true
// when the byte range did not unambiguously match one encoding and a fallback
// was chosen.
type DetectionResult struct {
	Encoding  fastq.QualityEncoding
	Name      string
	MinByte   int
	MaxByte   int
	Sampled   int
	Ambiguous bool
}

// detectionPeekSize is the maximum number of raw bytes DetectEncoding will
// peek from the input. 256 KiB comfortably contains the first 10000 quality
// characters of typical short-read FASTQ files (~150 bp reads × 4 lines ≈ 600
// bytes/record → ~430 records in 256 KiB → tens of thousands of qual bytes).
const detectionPeekSize = 256 * 1024

// detectionMaxQualBytes caps how many quality characters we consider when
// computing min/max. Matches the spec; more is wasteful, fewer is less robust.
const detectionMaxQualBytes = 10000

// DetectEncoding inspects the head of a FASTQ stream and returns the inferred
// quality encoding. It uses bufio.Reader.Peek so the bytes it inspects remain
// available for downstream readers — the typical pattern is:
//
//	br := bufio.NewReaderSize(r, 256*1024)
//	res, err := sickle.DetectEncoding(br)
//	// ... then pass br to fastq.NewReader; the quality data peeked at is
//	// still in the buffer and will be re-read normally.
//
// Detection algorithm (matches what most QC tools do):
//   - Scan up to the first 10000 quality bytes.
//   - If min < 33: invalid FASTQ — return an error.
//   - If min < 64 && max <= 73: Phred+33 (sanger).
//   - If min >= 64 && max <= 104: Phred+64 (illumina).
//   - If min >= 59 && min < 64: Solexa+64 (reported as fastq.Phred64,
//     Name="solexa").
//   - min < 64 && max > 73: Phred+33 (Illumina 1.8+ where ASCII can rise above
//     73). This case is also flagged Ambiguous=false.
//   - Anything else falls back to Phred+33 with Ambiguous=true.
//
// If the stream contains no quality bytes at all (empty / not a FASTQ), an
// error is returned.
func DetectEncoding(br *bufio.Reader) (DetectionResult, error) {
	res := DetectionResult{Encoding: fastq.Phred33, Name: "sanger"}

	peeked, err := br.Peek(detectionPeekSize)
	// io.EOF / ErrBufferFull are fine — we just work with whatever we got.
	if err != nil && err != io.EOF && err != bufio.ErrBufferFull {
		return res, fmt.Errorf("peeking FASTQ for encoding detection: %w", err)
	}
	if len(peeked) == 0 {
		return res, fmt.Errorf("empty input: no FASTQ data to detect encoding from")
	}

	minQ, maxQ, sampled, scanErr := scanQualityBytes(peeked, detectionMaxQualBytes)
	if scanErr != nil {
		return res, scanErr
	}
	if sampled == 0 {
		return res, fmt.Errorf("no quality bytes found in first %d bytes of input", len(peeked))
	}
	res.MinByte = minQ
	res.MaxByte = maxQ
	res.Sampled = sampled

	if minQ < 33 {
		return res, fmt.Errorf("invalid FASTQ: quality byte %d is below ASCII 33", minQ)
	}

	switch {
	case minQ < 64 && maxQ <= 73:
		res.Encoding = fastq.Phred33
		res.Name = "sanger"
	case minQ >= 64 && maxQ <= 104:
		res.Encoding = fastq.Phred64
		res.Name = "illumina"
	case minQ >= 59 && minQ < 64:
		// Solexa+64. The fastq library doesn't have a dedicated constant, so
		// we treat it as Phred+64 for read offsetting but label it "solexa".
		res.Encoding = fastq.Phred64
		res.Name = "solexa"
	case minQ < 64 && maxQ > 73:
		// Illumina 1.8+ allows quality up to 41, which encodes to ASCII 74 ('J').
		// Still Phred+33.
		res.Encoding = fastq.Phred33
		res.Name = "sanger"
	default:
		res.Encoding = fastq.Phred33
		res.Name = "sanger"
		res.Ambiguous = true
	}

	return res, nil
}

// scanQualityBytes walks the peeked buffer line-by-line treating it as FASTQ
// and accumulates min/max byte values across up to maxBytes quality
// characters. It tolerates a partial trailing record (it just stops at EOF).
// Returns minByte (255 if none seen), maxByte (0 if none seen), and the count
// of quality bytes inspected.
func scanQualityBytes(buf []byte, maxBytes int) (int, int, int, error) {
	minQ := 255
	maxQ := 0
	sampled := 0

	// Walk buffer line-by-line.
	lineNo := 0 // 0-indexed within the current 4-line FASTQ record
	i := 0
	for i < len(buf) {
		// Find the end of this line.
		lineStart := i
		for i < len(buf) && buf[i] != '\n' {
			i++
		}
		// We've either found '\n' or hit end of buffer (partial line).
		complete := i < len(buf)
		line := buf[lineStart:i]
		// Trim trailing '\r' if present (Windows line endings).
		if len(line) > 0 && line[len(line)-1] == '\r' {
			line = line[:len(line)-1]
		}
		if complete {
			i++ // skip the '\n'
		}

		switch lineNo {
		case 0: // header — must start with '@'
			if !complete {
				// Partial header at end of buffer; stop.
				return minQ, maxQ, sampled, nil
			}
			if len(line) == 0 || line[0] != '@' {
				return minQ, maxQ, sampled, fmt.Errorf("invalid FASTQ: expected '@' header, got %q", string(line))
			}
		case 1: // sequence
			if !complete {
				return minQ, maxQ, sampled, nil
			}
		case 2: // separator — must start with '+'
			if !complete {
				return minQ, maxQ, sampled, nil
			}
			if len(line) == 0 || line[0] != '+' {
				return minQ, maxQ, sampled, fmt.Errorf("invalid FASTQ: expected '+' separator, got %q", string(line))
			}
		case 3: // quality
			// Even if this line is partial (no newline), the bytes we have
			// are still valid quality bytes — count them.
			for _, b := range line {
				bv := int(b)
				if bv < minQ {
					minQ = bv
				}
				if bv > maxQ {
					maxQ = bv
				}
				sampled++
				if sampled >= maxBytes {
					return minQ, maxQ, sampled, nil
				}
			}
			if !complete {
				return minQ, maxQ, sampled, nil
			}
		}

		lineNo = (lineNo + 1) % 4
	}

	return minQ, maxQ, sampled, nil
}

// EncodingFromName maps a CLI quality-type string ("sanger", "illumina",
// "solexa", "phred33", "phred64") to the underlying fastq.QualityEncoding. It
// returns an error for unknown names so callers can surface a clear message
// instead of silently defaulting.
func EncodingFromName(name string) (fastq.QualityEncoding, error) {
	switch name {
	case "sanger", "phred33":
		return fastq.Phred33, nil
	case "illumina", "phred64":
		return fastq.Phred64, nil
	case "solexa":
		// Solexa uses a different formula but shares the +64 ASCII offset for
		// the practical range we encounter; treat as Phred+64 for I/O.
		return fastq.Phred64, nil
	default:
		return fastq.Phred33, fmt.Errorf("unknown quality type %q (expected sanger, illumina, solexa, or auto)", name)
	}
}
