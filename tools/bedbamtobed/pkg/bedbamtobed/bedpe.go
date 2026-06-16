package bedbamtobed

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// convertBedPE implements the -bedpe conversion. It assumes the input is
// grouped or sorted by query name (as upstream requires) and joins consecutive
// records with the same QNAME into a single BEDPE line.
//
// The pairing logic mirrors upstream ConvertBamToBedpe: it reads two records at
// a time; when their names differ it advances bam1 until it finds a matching
// mate (warning on stderr about orphaned paired reads), then emits the pair;
// when bam1 and bam2 share a name and both are paired, it emits them directly.
//
// Upstream's BamReader::GetNextAlignment leaves the destination record
// unchanged when it returns false at EOF. We reproduce that "stale on EOF"
// behaviour precisely so the terminal pairing matches byte-for-byte.
func convertBedPE(sr sam.Reader, w io.Writer, opts Options) (int, error) {
	n := 0
	// next returns (record, gotOne, err). gotOne is false at EOF.
	next := func() (*sam.Record, bool, error) {
		rec, err := sr.Read()
		if err == io.EOF {
			return nil, false, nil
		}
		if err != nil {
			return nil, false, err
		}
		n++
		return rec, true, nil
	}

	var bam1, bam2 *sam.Record
	got1, err := false, error(nil)
	bam1, got1, err = next()
	if err != nil {
		return n, err
	}
	for got1 {
		// GetNextAlignment(bam2); on EOF bam2 keeps its previous value.
		if rec, got2, err := next(); err != nil {
			return n, err
		} else if got2 {
			bam2 = rec
		}

		name1, name2 := nameOf(bam1), nameOf(bam2)
		if name1 != name2 {
			for nameOf(bam1) != nameOf(bam2) {
				// Upstream emits a "marked as paired but mate not adjacent"
				// warning to stderr here; it does not affect the BED output, so
				// we omit it to keep the converter pure (the CLI layer reports
				// it if desired).
				bam1 = bam2
				rec, got2, err := next()
				if err != nil {
					return n, err
				}
				if got2 {
					bam2 = rec
				} else {
					// EOF: bam2 stays stale; names now equal -> loop exits.
					break
				}
			}
			if err := printBedPE(w, bam1, bam2, opts); err != nil {
				return n, err
			}
		} else if bam1 != nil && bam2 != nil && bam1.IsPaired() && bam2.IsPaired() {
			if err := printBedPE(w, bam1, bam2, opts); err != nil {
				return n, err
			}
		}

		bam1, got1, err = next()
		if err != nil {
			return n, err
		}
	}
	return n, nil
}

// nameOf returns a record's QNAME, or "" for a nil record.
func nameOf(rec *sam.Record) string {
	if rec == nil {
		return ""
	}
	return rec.QName
}

// peMate holds the per-mate fields extracted for a BEDPE line. An unmapped mate
// keeps the upstream sentinels: chrom/strand "." and start/end -1.
type peMate struct {
	chrom  string
	start  int
	end    int
	strand string
	editND int64
	mapped bool
	mapQ   uint8
}

// extractMate pulls the BEDPE fields from a single alignment, returning the
// upstream sentinel values when the read is unmapped.
func extractMate(rec *sam.Record, useED bool) (peMate, error) {
	m := peMate{chrom: ".", start: -1, end: -1, strand: ".", mapped: false}
	if rec == nil || rec.IsUnmapped() {
		return m, nil
	}
	m.mapped = true
	m.chrom = rec.RName
	m.start = int(rec.Pos) - 1
	m.end = m.start + rec.Cigar.ReferenceLength()
	m.strand = strandOf(rec)
	m.mapQ = rec.MapQ
	if useED {
		a, ok := rec.GetAux("NM")
		if !ok {
			return m, fmt.Errorf("The edit distance tag (NM) was not found in the BAM file.  Please disable -ed.  Exiting")
		}
		v, _ := a.Int()
		m.editND = v
	}
	return m, nil
}

// printBedPE writes a single BEDPE line from two mates, applying the upstream
// end-swap rules and score selection (min MAPQ or summed edit distance).
func printBedPE(w io.Writer, bam1, bam2 *sam.Record, opts Options) error {
	m1, err := extractMate(bam1, opts.UseEditDistance)
	if err != nil {
		return err
	}
	m2, err := extractMate(bam2, opts.UseEditDistance)
	if err != nil {
		return err
	}

	if !opts.Mate1First {
		// Order by (chrom, start): swap so the lexicographically-smaller end
		// (and, on ties, the smaller start) is reported first.
		if m1.chrom > m2.chrom || (m1.chrom == m2.chrom && m1.start > m2.start) {
			m1, m2 = m2, m1
		}
	} else {
		// Always report mate 1 (FlagRead1) first.
		if bam1 != nil && !bam1.IsRead1() {
			m1, m2 = m2, m1
		}
	}

	name := nameOf(bam1)

	var score int64
	if !opts.UseEditDistance {
		// Minimum MAPQ across the two mapped mates (0 if either is unmapped).
		if m1.mapped && m2.mapped {
			if m1.mapQ < m2.mapQ {
				score = int64(m1.mapQ)
			} else {
				score = int64(m2.mapQ)
			}
		}
	} else {
		switch {
		case m1.mapped && m2.mapped:
			score = m1.editND + m2.editND
		case m1.mapped:
			score = m1.editND
		case m2.mapped:
			score = m2.editND
		}
	}

	_, err = fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n",
		m1.chrom, m1.start, m1.end,
		m2.chrom, m2.start, m2.end,
		name, score, m1.strand, m2.strand)
	return err
}
