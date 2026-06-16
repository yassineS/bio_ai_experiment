package samtools

import (
	"bytes"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// IdxstatsRow is one line of `samtools idxstats` output: reference name,
// reference length in bases, number of mapped reads aligned to it, and
// number of unmapped reads anchored on it (plus, for the trailing "*"
// row, the count of fully unplaced unmapped records).
type IdxstatsRow struct {
	Name     string
	Length   int32
	Mapped   uint64
	Unmapped uint64
}

// Idxstats reads the BAI index sitting next to bamPath (i.e. <bamPath>.bai)
// and returns one row per @SQ entry from the BAM header, plus a final
// "*\t0\t0\t<n_no_coor>" row covering the truly unplaced records. The BAI
// meta pseudo-bin (37450) carries the per-reference (mapped, unmapped)
// counts that we already populate during indexing, so the operation is
// O(refs) rather than O(records).
//
// When the .bai file is missing Idxstats falls back to a full streaming
// scan of the BAM body (matching upstream behaviour when index is absent
// is to error out — but we keep parity with htslib's behaviour by
// returning a clear error pointing at how to build the index).
func Idxstats(bamPath string) ([]IdxstatsRow, error) {
	return IdxstatsThreaded(bamPath, 0)
}

// IdxstatsThreaded is Idxstats with block-parallel BGZF decode wired to a
// thread count (`-@/--threads`). The index path (the common case) does not
// decode the BAM body and is unaffected by threads; only the no-index fallback
// scan inflates the whole BAM body, and there threads >= 2 spreads the BGZF
// decode across worker goroutines. The rows are identical for any thread count.
func IdxstatsThreaded(bamPath string, threads int) ([]IdxstatsRow, error) {
	in, err := openSeekable(bamPath)
	if err != nil {
		return nil, err
	}
	defer in.Close()
	br, err := sam.NewBAMReader(in)
	if err != nil {
		return nil, err
	}
	defer br.Close()
	hdr := br.Header()

	baiPath := bamPath + ".bai"
	if baiBytes, berr := hfile.ReadFile(baiPath); berr == nil {
		idx, ierr := bam.ReadBAI(bytes.NewReader(baiBytes))
		if ierr == nil {
			return idxstatsFromIndex(hdr, idx), nil
		}
		// fall through to scan if BAI is corrupt
	}

	// No usable index — count records directly. This matches upstream's
	// fallback when the BAI is missing or corrupt: samtools 1.x errors
	// out, but htslib 1.10+ scans linearly. We pick the more useful
	// behaviour and surface a stderr note via the caller. Re-open through the
	// thread-aware reader so the linear scan benefits from -@ parallel BGZF
	// decode; the already-open br is only used above for the header.
	if threads >= 2 {
		if rows, ok, err := idxstatsScanThreaded(bamPath, threads); ok {
			return rows, err
		}
	}
	return idxstatsFromScan(br)
}

// idxstatsScanThreaded performs the no-index linear scan through the thread-aware
// alnio reader so BGZF blocks are inflated in parallel. It returns ok=false when
// the file cannot be re-opened (e.g. a non-seekable stream), so the caller can
// fall back to the single-threaded scan over the already-open reader.
func idxstatsScanThreaded(bamPath string, threads int) ([]IdxstatsRow, bool, error) {
	rd, err := alnio.OpenReaderThreaded(bamPath, "", threads)
	if err != nil {
		return nil, false, nil
	}
	defer rd.Close()
	hdr := rd.Header()
	rows := make([]IdxstatsRow, len(hdr.Refs)+1)
	for i, r := range hdr.Refs {
		rows[i].Name = r.Name
		rows[i].Length = r.Length
	}
	rows[len(hdr.Refs)].Name = "*"
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, true, err
		}
		if rec.RName == "" || rec.RName == "*" {
			rows[len(hdr.Refs)].Unmapped++
			continue
		}
		idx := hdr.RefIndex(rec.RName)
		if idx < 0 {
			rows[len(hdr.Refs)].Unmapped++
			continue
		}
		if rec.IsUnmapped() {
			rows[idx].Unmapped++
		} else {
			rows[idx].Mapped++
		}
	}
	return rows, true, nil
}

// idxstatsFromIndex builds the per-reference rows by reading each
// reference's BAI meta pseudo-bin.
func idxstatsFromIndex(hdr *sam.Header, idx *bam.BAIIndex) []IdxstatsRow {
	out := make([]IdxstatsRow, 0, len(hdr.Refs)+1)
	for i, ref := range hdr.Refs {
		var mapped, unmapped uint64
		if mp, um, ok := idx.MetaCounts(i); ok {
			mapped = mp
			unmapped = um
		}
		out = append(out, IdxstatsRow{
			Name:     ref.Name,
			Length:   ref.Length,
			Mapped:   mapped,
			Unmapped: unmapped,
		})
	}
	out = append(out, IdxstatsRow{
		Name:     "*",
		Length:   0,
		Mapped:   0,
		Unmapped: idx.NoCoor,
	})
	return out
}

// idxstatsFromScan walks every record in the BAM and tallies (mapped,
// unmapped) per reference. Used when no .bai is available.
func idxstatsFromScan(br *sam.BAMReader) ([]IdxstatsRow, error) {
	hdr := br.Header()
	rows := make([]IdxstatsRow, len(hdr.Refs)+1)
	for i, r := range hdr.Refs {
		rows[i].Name = r.Name
		rows[i].Length = r.Length
	}
	rows[len(hdr.Refs)].Name = "*"
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if rec.RName == "" || rec.RName == "*" {
			rows[len(hdr.Refs)].Unmapped++
			continue
		}
		idx := hdr.RefIndex(rec.RName)
		if idx < 0 {
			rows[len(hdr.Refs)].Unmapped++
			continue
		}
		if rec.IsUnmapped() {
			rows[idx].Unmapped++
		} else {
			rows[idx].Mapped++
		}
	}
	return rows, nil
}

// WriteIdxstats serialises rows as the canonical text table:
// "name\tlen\tmapped\tunmapped\n" per row.
func WriteIdxstats(w io.Writer, rows []IdxstatsRow) error {
	for _, r := range rows {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\n", r.Name, r.Length, r.Mapped, r.Unmapped); err != nil {
			return err
		}
	}
	return nil
}

// IdxstatsFile is the high-level CLI entry point: reads the BAM at path,
// computes the rows, and writes them to w. threads (`-@/--threads`) parallelises
// the BGZF decode of the no-index fallback scan; the output is identical for any
// thread count.
func IdxstatsFile(path string, w io.Writer, threads int) error {
	rows, err := IdxstatsThreaded(path, threads)
	if err != nil {
		return err
	}
	return WriteIdxstats(w, rows)
}
