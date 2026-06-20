package samtools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// FlagstatCounts is the tally of flag combinations needed to render the
// 13-line samtools-style flagstat report. Each counter is a [QC-passed,
// QC-failed] pair where QC-failed sums records with the 0x200 flag bit set.
type FlagstatCounts struct {
	Total             [2]int
	Primary           [2]int
	Secondary         [2]int
	Supplementary     [2]int
	Duplicates        [2]int
	PrimaryDuplicates [2]int
	Mapped            [2]int
	PrimaryMapped     [2]int
	Paired            [2]int
	Read1             [2]int
	Read2             [2]int
	ProperlyPaired    [2]int
	WithItselfAndMate [2]int
	Singletons        [2]int
	MateDiffChr       [2]int
	MateDiffChrMapq5  [2]int
}

// CountFlagstat consumes a SAM/BAM stream from r and returns the flagstat
// tallies. The header is parsed first via sam.NewReader.
func CountFlagstat(r io.Reader) (*FlagstatCounts, error) {
	return CountFlagstatThreaded(r, 0)
}

// CountFlagstatThreaded is CountFlagstat with block-parallel BGZF input decode
// wired to a thread count (`-@/--threads`). When threads >= 2 and the input is
// a BGZF-wrapped BAM, its blocks are inflated concurrently; the tallies are
// identical for any thread count because the decoded record stream is identical.
func CountFlagstatThreaded(r io.Reader, threads int) (*FlagstatCounts, error) {
	rd, err := alnio.NewReaderThreaded(r, "", threads)
	if err != nil {
		return nil, err
	}
	if rc, ok := rd.(io.Closer); ok {
		defer rc.Close()
	}
	return countFromReader(rd)
}

// add updates the counters for one record.
func (c *FlagstatCounts) add(r *sam.Record) {
	idx := 0
	if r.IsQCFail() {
		idx = 1
	}
	c.Total[idx]++

	primary := r.IsPrimary()
	switch {
	case r.IsSecondary():
		c.Secondary[idx]++
	case r.IsSupplementary():
		c.Supplementary[idx]++
	default:
		c.Primary[idx]++
	}

	if r.IsDuplicate() {
		c.Duplicates[idx]++
		if primary {
			c.PrimaryDuplicates[idx]++
		}
	}

	if r.IsMapped() {
		c.Mapped[idx]++
		if primary {
			c.PrimaryMapped[idx]++
		}
	}

	if !r.IsPaired() {
		return
	}
	c.Paired[idx]++
	if r.IsRead1() {
		c.Read1[idx]++
	}
	if r.IsRead2() {
		c.Read2[idx]++
	}
	if r.IsProperPair() {
		c.ProperlyPaired[idx]++
	}
	// Mate-related accounting only applies to mapped records.
	if r.IsUnmapped() {
		return
	}
	if !r.IsMateUnmapped() {
		c.WithItselfAndMate[idx]++
		// Mate on a different chromosome. Upstream (bam_stat.c
		// flagstat_loop) compares decoded reference ids, c->mtid != c->tid,
		// inside the both-mates-mapped branch. Here the read is mapped
		// (RNAME != "*") and the mate is mapped (FMUNMAP clear), so the mate
		// id differs from the read id unless RNEXT names the read's own
		// reference — either via the "=" reflexive marker or by spelling the
		// same RNAME. Every other RNEXT value, including "*"/"" (which htslib
		// decodes to mtid=-1, distinct from a mapped read's tid>=0), counts as
		// a different chr. Comparing ref ids rather than the raw RNEXT string
		// fixes both the RNEXT==RNAME over-count and the RNEXT=="*"-with-mate-
		// mapped under-count.
		if r.RNext != "=" && r.RNext != r.RName {
			c.MateDiffChr[idx]++
			if r.MapQ >= 5 {
				c.MateDiffChrMapq5[idx]++
			}
		}
	} else {
		c.Singletons[idx]++
	}
}

// Format writes the 13-line samtools-style flagstat report to w.
func (c *FlagstatCounts) Format(w io.Writer) error {
	pct := func(num, denom int) string {
		if denom == 0 {
			return "N/A"
		}
		return fmt.Sprintf("%.2f%%", 100.0*float64(num)/float64(denom))
	}

	lines := []string{
		fmt.Sprintf("%d + %d in total (QC-passed reads + QC-failed reads)", c.Total[0], c.Total[1]),
		fmt.Sprintf("%d + %d primary", c.Primary[0], c.Primary[1]),
		fmt.Sprintf("%d + %d secondary", c.Secondary[0], c.Secondary[1]),
		fmt.Sprintf("%d + %d supplementary", c.Supplementary[0], c.Supplementary[1]),
		fmt.Sprintf("%d + %d duplicates", c.Duplicates[0], c.Duplicates[1]),
		fmt.Sprintf("%d + %d primary duplicates", c.PrimaryDuplicates[0], c.PrimaryDuplicates[1]),
		fmt.Sprintf("%d + %d mapped (%s : %s)", c.Mapped[0], c.Mapped[1], pct(c.Mapped[0], c.Total[0]), pct(c.Mapped[1], c.Total[1])),
		fmt.Sprintf("%d + %d primary mapped (%s : %s)", c.PrimaryMapped[0], c.PrimaryMapped[1], pct(c.PrimaryMapped[0], c.Primary[0]), pct(c.PrimaryMapped[1], c.Primary[1])),
		fmt.Sprintf("%d + %d paired in sequencing", c.Paired[0], c.Paired[1]),
		fmt.Sprintf("%d + %d read1", c.Read1[0], c.Read1[1]),
		fmt.Sprintf("%d + %d read2", c.Read2[0], c.Read2[1]),
		fmt.Sprintf("%d + %d properly paired (%s : %s)", c.ProperlyPaired[0], c.ProperlyPaired[1], pct(c.ProperlyPaired[0], c.Paired[0]), pct(c.ProperlyPaired[1], c.Paired[1])),
		fmt.Sprintf("%d + %d with itself and mate mapped", c.WithItselfAndMate[0], c.WithItselfAndMate[1]),
		fmt.Sprintf("%d + %d singletons (%s : %s)", c.Singletons[0], c.Singletons[1], pct(c.Singletons[0], c.Paired[0]), pct(c.Singletons[1], c.Paired[1])),
		fmt.Sprintf("%d + %d with mate mapped to a different chr", c.MateDiffChr[0], c.MateDiffChr[1]),
		fmt.Sprintf("%d + %d with mate mapped to a different chr (mapQ>=5)", c.MateDiffChrMapq5[0], c.MateDiffChrMapq5[1]),
	}
	for _, ln := range lines {
		if _, err := fmt.Fprintln(w, ln); err != nil {
			return err
		}
	}
	return nil
}

// Flagstat is the high-level entry point used by the CLI: reads from r,
// writes the report to w.
func Flagstat(r io.Reader, w io.Writer) error {
	return FlagstatThreaded(r, w, 0)
}

// FlagstatThreaded is Flagstat with block-parallel BGZF input decode wired to a
// thread count (`-@/--threads`). The report is byte-identical for any thread
// count; threading only changes decode throughput.
func FlagstatThreaded(r io.Reader, w io.Writer, threads int) error {
	c, err := CountFlagstatThreaded(r, threads)
	if err != nil {
		return err
	}
	return c.Format(w)
}

// FlagstatFile opens the alignment file at path (SAM, BAM, or CRAM; "-" is
// stdin) and writes its flagstat report to w. When threads >= 2 and the input
// is a BGZF-wrapped BAM, the BGZF blocks are inflated concurrently across up to
// threads worker goroutines; the report is identical to the single-threaded
// path. Opening the raw file here (rather than a pre-decompressed stream) is
// what lets the parallel BGZF reader engage.
func FlagstatFile(path string, w io.Writer, threads int) error {
	rd, err := alnio.OpenReaderThreaded(path, "", threads)
	if err != nil {
		return err
	}
	defer rd.Close()
	c, err := countFromReader(rd)
	if err != nil {
		return err
	}
	return c.Format(w)
}

// countFromReader tallies flagstat counts from rd, preferring the reader's
// shallow decode (fixed-prefix fields only) when available, then the
// allocation-free ReadInto, then a plain Read. flagstat never retains a record
// past add, so reusing one record across calls is safe. The shallow path is
// the common BAM case and skips the variable-length region (name/CIGAR/SEQ/
// QUAL/aux) that flagstat never reads.
func countFromReader(rd interface {
	Read() (*sam.Record, error)
}) (*FlagstatCounts, error) {
	c := &FlagstatCounts{}
	if rs, ok := rd.(interface {
		ReadShallowInto(*sam.Record) error
	}); ok {
		var rec sam.Record
		for {
			err := rs.ReadShallowInto(&rec)
			if err == io.EOF {
				break
			}
			if err != nil {
				return nil, err
			}
			c.add(&rec)
		}
		return c, nil
	}
	if ri, ok := rd.(interface {
		ReadInto(*sam.Record) error
	}); ok {
		var rec sam.Record
		for {
			err := ri.ReadInto(&rec)
			if err == io.EOF {
				break
			}
			if err != nil {
				return nil, err
			}
			c.add(&rec)
		}
		return c, nil
	}
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		c.add(rec)
	}
	return c, nil
}
