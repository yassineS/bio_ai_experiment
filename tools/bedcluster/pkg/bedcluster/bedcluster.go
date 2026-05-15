// Package bedcluster ports `bedtools cluster`: it groups overlapping
// intervals into clusters and tags each input record with a cluster ID.
//
// Unlike `bedmerge`, no records are collapsed — every input record appears
// in the output, with an extra trailing column holding its cluster ID.
//
// Two intervals end up in the same cluster when they are within MaxDistance
// bp of each other (default 0 = require overlap, matching upstream's
// `-d 0` default). With StrandSpec, only same-strand intervals cluster.
package bedcluster

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// Options configures Cluster.
type Options struct {
	// MaxDistance is the largest gap (in bp) between two intervals that still
	// puts them in the same cluster. Default 0 = require overlap (book-ends
	// included: an interval ending at 100 and one starting at 100 cluster
	// when MaxDistance=0, matching upstream).
	MaxDistance int
	// StrandSpec, when true, only clusters intervals with the same strand.
	StrandSpec bool
}

// row is one parsed input line, kept with all its raw fields so we can
// re-emit them verbatim in the original column count.
type row struct {
	chrom  string
	start  int
	end    int
	strand string
	fields []string
}

// Cluster reads BED records from r, assigns each one a cluster ID, and writes
// each record back to w with the cluster ID appended as the final column.
//
// Records are emitted in genomic-sort order (chrom, then — for StrandSpec —
// strand, then start, then end), which matches upstream's behaviour: the
// upstream tool requires sorted input and groups by strand internally when
// `-s` is set, producing per-strand contiguous output.
//
// Returns the number of records written.
func Cluster(r io.Reader, w io.Writer, opts Options) (int, error) {
	rows, err := readAll(r)
	if err != nil {
		return 0, err
	}
	if len(rows) == 0 {
		return 0, nil
	}

	// Stable sort by (chrom [, strand], start) only. Equal-start records
	// keep their input order — matches upstream's behaviour, which only
	// requires that input be sorted by chrom+start (not end).
	sort.SliceStable(rows, func(i, j int) bool {
		if rows[i].chrom != rows[j].chrom {
			return rows[i].chrom < rows[j].chrom
		}
		if opts.StrandSpec && rows[i].strand != rows[j].strand {
			return rows[i].strand < rows[j].strand
		}
		return rows[i].start < rows[j].start
	})

	clusterIDs := assignClusters(rows, opts)

	bw := bufio.NewWriter(w)
	defer bw.Flush()
	for i, r := range rows {
		fields := append(append([]string(nil), r.fields...), strconv.Itoa(clusterIDs[i]))
		if _, err := fmt.Fprintln(bw, strings.Join(fields, "\t")); err != nil {
			return 0, err
		}
	}
	return len(rows), nil
}

// assignClusters returns a slice (parallel to the sorted `rows`) of cluster
// IDs. Cluster IDs are 1-based and assigned in walk order.
func assignClusters(rows []row, opts Options) []int {
	ids := make([]int, len(rows))
	if len(rows) == 0 {
		return ids
	}
	curID := 1
	curChrom := rows[0].chrom
	curStrand := rows[0].strand
	curEnd := rows[0].end
	ids[0] = curID

	for i := 1; i < len(rows); i++ {
		r := rows[i]
		newCluster := false
		switch {
		case r.chrom != curChrom:
			newCluster = true
		case opts.StrandSpec && r.strand != curStrand:
			newCluster = true
		case r.start > curEnd+opts.MaxDistance:
			// Strict gap: a gap > MaxDistance opens a new cluster.
			newCluster = true
		}
		if newCluster {
			curID++
			curChrom = r.chrom
			curStrand = r.strand
			curEnd = r.end
		} else if r.end > curEnd {
			curEnd = r.end
		}
		ids[i] = curID
	}
	return ids
}

func readAll(r io.Reader) ([]row, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var rows []row
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
		strand := ""
		if len(fields) >= 6 {
			strand = fields[5]
		}
		rows = append(rows, row{
			chrom:  fields[0],
			start:  start,
			end:    end,
			strand: strand,
			fields: fields,
		})
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("reading BED: %w", err)
	}
	return rows, nil
}
