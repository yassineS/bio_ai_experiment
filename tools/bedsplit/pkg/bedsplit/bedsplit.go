// Package bedsplit implements `bedtools split`: it splits a single BED file
// into N approximately-equal-sized output files.
//
// Two algorithms are supported (matching upstream's `-a {simple|size}`):
//
//   - simple: split by record count. Records are partitioned in input order
//     so each file gets ceil(total/N) records (the last file may have less).
//   - size:   split by total bp. Records are sorted by length descending and
//     greedily assigned to the currently-smallest bin (LPT scheduling), so
//     bp totals across files are balanced.
//
// The CLI writes each shard to `<prefix>.NNNNN.bed` (5-digit zero-padded,
// 1-based) and emits a manifest TSV to stdout: `filename<TAB>total_bp<TAB>num_records`.
package bedsplit

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cppsort"
)

// Algorithm selects the partitioning strategy.
type Algorithm int

const (
	// AlgSimple partitions records in input order, equal counts per file.
	AlgSimple Algorithm = iota
	// AlgSize partitions records by total bp, balancing bp across files.
	AlgSize
)

// ParseAlgorithm maps a -a flag value to an Algorithm.
func ParseAlgorithm(s string) (Algorithm, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "simple":
		return AlgSimple, nil
	case "size":
		return AlgSize, nil
	default:
		return 0, fmt.Errorf("unknown -a algorithm %q (want \"simple\" or \"size\")", s)
	}
}

// Options configures Split.
type Options struct {
	N         int       // number of output files (must be >= 1)
	Prefix    string    // output filename prefix
	Algorithm Algorithm // partitioning strategy
}

// ManifestRow describes one output shard.
type ManifestRow struct {
	Filename   string
	TotalBP    int
	NumRecords int
}

// record holds one parsed input line.
type record struct {
	length int
	line   string
}

// Split reads BED records from r, partitions them per opts, writes one shard
// file per partition to disk, and emits a manifest TSV to manifestW.
//
// Output filenames are `<prefix>.NNNNN.bed` (1-based, zero-padded to 5 digits),
// matching upstream behaviour.
func Split(r io.Reader, manifestW io.Writer, opts Options) ([]ManifestRow, error) {
	if opts.N < 1 {
		return nil, fmt.Errorf("split count -n must be >= 1, got %d", opts.N)
	}
	if opts.Prefix == "" {
		return nil, fmt.Errorf("output prefix -p is required")
	}

	records, err := readAll(r)
	if err != nil {
		return nil, err
	}
	if len(records) == 0 {
		return nil, nil
	}

	// Cap N at the number of records — empty files would otherwise be
	// produced; upstream skips them as well.
	n := opts.N
	if n > len(records) {
		n = len(records)
	}

	var bins [][]int // bins[i] = indices into `records` going into shard i
	switch opts.Algorithm {
	case AlgSimple:
		bins = simpleBins(len(records), n)
	case AlgSize:
		bins = sizeBins(records, n)
	default:
		return nil, fmt.Errorf("unknown algorithm")
	}

	// Write each shard.
	manifest := make([]ManifestRow, 0, len(bins))
	for i, idxs := range bins {
		fname := fmt.Sprintf("%s.%05d.bed", opts.Prefix, i+1)
		f, err := os.Create(fname)
		if err != nil {
			return manifest, fmt.Errorf("creating %s: %w", fname, err)
		}
		bw := bufio.NewWriter(f)
		totalBP := 0
		for _, ix := range idxs {
			rec := records[ix]
			totalBP += rec.length
			if _, err := bw.WriteString(rec.line); err != nil {
				_ = f.Close()
				return manifest, fmt.Errorf("writing %s: %w", fname, err)
			}
			if err := bw.WriteByte('\n'); err != nil {
				_ = f.Close()
				return manifest, fmt.Errorf("writing %s: %w", fname, err)
			}
		}
		if err := bw.Flush(); err != nil {
			_ = f.Close()
			return manifest, fmt.Errorf("flushing %s: %w", fname, err)
		}
		if err := f.Close(); err != nil {
			return manifest, fmt.Errorf("closing %s: %w", fname, err)
		}
		manifest = append(manifest, ManifestRow{
			Filename:   fname,
			TotalBP:    totalBP,
			NumRecords: len(idxs),
		})
	}

	// Emit manifest TSV.
	if manifestW != nil {
		bw := bufio.NewWriter(manifestW)
		for _, row := range manifest {
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\n", row.Filename, row.TotalBP, row.NumRecords); err != nil {
				return manifest, err
			}
		}
		if err := bw.Flush(); err != nil {
			return manifest, err
		}
	}
	return manifest, nil
}

// simpleBins assigns the i-th input record to bin (i % n), matching upstream
// `bedtools split -a simple` (a strict round-robin distribution).
func simpleBins(total, n int) [][]int {
	bins := make([][]int, n)
	for i := 0; i < total; i++ {
		bins[i%n] = append(bins[i%n], i)
	}
	return bins
}

// sizeBins ports upstream `bedtools split -a size` ("doEuristicSplitOnTotalSize"):
//
//  1. Sort all records by length descending (stable).
//  2. Open n empty bins.
//  3. For each record, pick the bin where placing it would yield the
//     lowest sum-of-absolute-deviations from the current expected mean.
//     The first bin to achieve that minimum wins (no further tie-breaker).
//  4. Within a bin, records are emitted in the order they were added —
//     i.e. size-descending order, since records are processed largest-first.
//     Upstream stores each BED* in the split's items vector as it is added
//     and writes them in that same vector order (splitBed.cpp saveBedItems),
//     so we must NOT re-sort the bin back into input order.
//
// This is O(records * n^2) but n is small in practice (number of shards).
func sizeBins(records []record, n int) [][]int {
	type lenIdx struct {
		length int
		idx    int
	}
	byLen := make([]lenIdx, len(records))
	for i, r := range records {
		byLen[i] = lenIdx{length: r.length, idx: i}
	}
	// Upstream sorts with std::sort(items, sortBySizeDesc) — a length-ONLY
	// comparator (aLen > bLen), so equal-length records come out in introsort's
	// artifact order, which determines the greedy bin assignment below. Use the
	// libstdc++ introsort port so the per-file record assignment matches
	// byte-for-byte (a stable sort would tie equal-length records in input order
	// and diverge).
	cppsort.Sort(byLen, func(a, b lenIdx) bool { return a.length > b.length })

	bins := make([][]int, 0, n)
	binTotals := make([]int, 0, n)
	totalBases := 0

	for _, li := range byLen {
		// Phase 1: each new bin starts with one record (no choice to make).
		if len(bins) < n {
			bins = append(bins, []int{li.idx})
			binTotals = append(binTotals, li.length)
			totalBases += li.length
			continue
		}
		if n == 1 {
			bins[0] = append(bins[0], li.idx)
			binTotals[0] += li.length
			totalBases += li.length
			continue
		}
		// Phase 2: pick the bin whose insertion minimises sum |bin_size - mean|.
		mean := float64(totalBases+li.length) / float64(n)
		bestIdx := 0
		bestDev := -1.0
		for try := 0; try < len(bins); try++ {
			dev := 0.0
			for i, t := range binTotals {
				size := float64(t)
				if i == try {
					size += float64(li.length)
				}
				diff := size - mean
				if diff < 0 {
					diff = -diff
				}
				dev += diff
			}
			if try == 0 || dev < bestDev {
				bestDev = dev
				bestIdx = try
			}
		}
		bins[bestIdx] = append(bins[bestIdx], li.idx)
		binTotals[bestIdx] += li.length
		totalBases += li.length
	}
	return bins
}

func readAll(r io.Reader) ([]record, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []record
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("BED record needs >=3 columns, got %d in line: %q", len(fields), line)
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		out = append(out, record{length: end - start, line: line})
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("reading BED: %w", err)
	}
	return out, nil
}
