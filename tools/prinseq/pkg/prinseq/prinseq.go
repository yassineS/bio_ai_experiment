package prinseq

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// Stats holds sequence statistics
type Stats struct {
	NumReads   int
	TotalBases int
	MinLength  int
	MaxLength  int
	AvgLength  float64
	GCContent  float64
	AvgQuality float64 // Only for FASTQ
	NumNs      int
}

// CalculateStats computes statistics for FASTA or FASTQ files
func CalculateStats(reader io.Reader, isFastq bool) (*Stats, error) {
	stats := &Stats{
		MinLength: -1,
	}

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 1024*1024), 1024*1024) // 1MB buffer

	if isFastq {
		return calculateFastqStats(scanner, stats)
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
			totalQuality += calculateAvgQualityScore(qual)
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
	if len(qual) == 0 {
		return 0.0
	}

	total := 0
	for _, q := range qual {
		// Assume Phred+33 encoding
		total += int(q) - 33
	}
	return float64(total) / float64(len(qual))
}

// FilterOptions holds filtering parameters
type FilterOptions struct {
	MinLen      int
	MaxLen      int
	MinGC       float64
	MaxGC       float64
	MinQual     float64
	MaxNsP      float64 // Max percentage of Ns
	MaxNsN      int     // Max number of Ns
	TrimLeft    int
	TrimRight   int
	TrimQualL   int
	TrimQualR   int
	MinQualMean float64
	MaxQualMean float64
}

// Filter filters a FASTA/FASTQ file based on the given options
func Filter(reader io.Reader, writer io.Writer, isFastq bool, opts FilterOptions) error {
	if isFastq {
		return filterFastq(reader, writer, opts)
	}
	return filterFasta(reader, writer, opts)
}

func filterFasta(reader io.Reader, writer io.Writer, opts FilterOptions) error {
	fastaReader := fasta.NewReader(reader)
	fastaWriter := fasta.NewWriter(writer, 80)

	for {
		record, err := fastaReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		// Apply filters
		if shouldFilterSequence(string(record.Sequence), "", opts) {
			continue
		}

		// Write filtered record
		if err := fastaWriter.Write(record); err != nil {
			return err
		}
	}

	return fastaWriter.Flush()
}

func filterFastq(reader io.Reader, writer io.Writer, opts FilterOptions) error {
	fastqReader := fastq.NewReader(reader, fastq.Phred33)
	fastqWriter := fastq.NewWriter(writer, fastq.Phred33)

	for {
		record, err := fastqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		// Apply filters
		if shouldFilterSequence(string(record.Sequence), string(record.Quality), opts) {
			continue
		}

		// Write filtered record
		if err := fastqWriter.Write(record); err != nil {
			return err
		}
	}

	return fastqWriter.Flush()
}

func shouldFilterSequence(seq, qual string, opts FilterOptions) bool {
	seqLen := len(seq)

	// Length filters
	if opts.MinLen > 0 && seqLen < opts.MinLen {
		return true
	}
	if opts.MaxLen > 0 && seqLen > opts.MaxLen {
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

	// Quality filters (only for FASTQ)
	if qual != "" {
		avgQual := calculateAvgQualityScore(qual)
		if opts.MinQualMean > 0 && avgQual < opts.MinQualMean {
			return true
		}
		if opts.MaxQualMean > 0 && avgQual > opts.MaxQualMean {
			return true
		}
	}

	return false
}
