// bcftools mpileup — indexed region fast path.
//
// `bcftools mpileup -r/-R` with an indexed, coordinate-sorted BAM must NOT
// drain the whole file into memory: a whole-genome BAM is tens of GB of
// records, and bucketing every record by chromosome (mpileupReadBAM) blew
// RSS to ~11 GB where upstream stays at ~106 MB. This file provides a
// region-restricted reader that seeks only the BAI/CSI index chunks
// covering the requested regions, so a `-r` query reads only the region's
// reads — bounded memory and wall.
//
// It is built directly on the shared pkg/htsgo primitives
// (bam.UnionChunks / bam.UnionChunksCSI, bgzf, sam, hfile, region) rather
// than the samtools-package helpers, so it adds zero churn to the
// just-merged samtools `mpileup -r` / `consensus -r` / `view` paths. The
// read+filter logic mirrors samtools' bamRegionReader exactly: chunks come
// out coordinate-ordered (UnionChunks returns them ascending and merged,
// and each chunk is itself coordinate-sorted), and the per-record overlap
// filter keeps left-overlapping reads (reads that start before the region
// but span into it) — the indel-context invariant that makes an indexed
// `-r` query over a whole-genome BAM byte-identical to upstream's, not to a
// `samtools view`-extracted sub-BAM.
package bcftools

import (
	"bytes"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
	hregion "github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// mpileupRegionChunks loads the BAM's .csi (preferred) or .bai index and
// returns the merged, coordinate-sorted file chunks covering the resolved
// regions, or ok=false when neither sibling index exists (the caller then
// linear-scans).
func mpileupRegionChunks(path string, regions []hregion.ResolvedRegion) (chunks []bam.BAIChunk, ok bool, err error) {
	if csiBytes, cerr := hfile.ReadFile(path + ".csi"); cerr == nil {
		idx, ierr := bam.ReadCSI(bytes.NewReader(csiBytes))
		if ierr != nil {
			return nil, false, ierr
		}
		return bam.UnionChunksCSI(idx, regions), true, nil
	}
	if baiBytes, berr := hfile.ReadFile(path + ".bai"); berr == nil {
		idx, ierr := bam.ReadBAI(bytes.NewReader(baiBytes))
		if ierr != nil {
			return nil, false, ierr
		}
		return bam.UnionChunks(idx, regions), true, nil
	}
	return nil, false, nil
}

// mpileupHeaderIsCoordinateSorted reports whether the @HD line declares
// SO:coordinate, the precondition for an index seek (records arrive in
// reference order). A BAM that is not coordinate-sorted routes to the
// linear fallback.
func mpileupHeaderIsCoordinateSorted(hdr *sam.Header) bool {
	for _, f := range hdr.HDFields {
		if f.Tag == "SO" {
			return f.Value == "coordinate"
		}
	}
	return false
}

// mpileupBuildRegionFilter returns nil when no regions are configured;
// otherwise returns a predicate that keeps records overlapping any region's
// range on the matching reference. It keeps a read whose alignment starts
// before the region but ends inside it (left-overlap), preserving the
// indel-context reads upstream's region query also keeps.
func mpileupBuildRegionFilter(regions []hregion.ResolvedRegion, hdr *sam.Header) func(*sam.Record) bool {
	if len(regions) == 0 {
		return nil
	}
	return func(rec *sam.Record) bool {
		if rec.RName == "" || rec.RName == "*" {
			return false
		}
		rid := hdr.RefIndex(rec.RName)
		if rid < 0 {
			return false
		}
		pos0 := int(rec.Pos) - 1
		if pos0 < 0 {
			pos0 = 0
		}
		refLen := rec.Cigar.ReferenceLength()
		if refLen <= 0 {
			refLen = 1
		}
		for _, r := range regions {
			if r.RefID != rid {
				continue
			}
			recEnd := pos0 + refLen
			if pos0 < r.End0 && recEnd > r.Beg0 {
				return true
			}
		}
		return false
	}
}

// mpileupChunkBoundedReader stops returning bytes once its underlying BGZF
// reader has advanced past the chunk-end virtual offset. We watch the
// virtual offset of the bgzip.Reader after every read and report io.EOF as
// soon as we cross the end boundary.
type mpileupChunkBoundedReader struct {
	r   *bgzip.Reader
	end uint64
}

// Read reads from the bounded BGZF stream, returning io.EOF once the
// chunk-end virtual offset is reached.
func (c *mpileupChunkBoundedReader) Read(p []byte) (int, error) {
	if c.r.VirtualOffset() >= c.end {
		return 0, io.EOF
	}
	return c.r.Read(p)
}

// bamRegionReader is a sam.Reader that yields only the records overlapping a
// set of regions by seeking the BAM's index chunks rather than scanning the
// whole file. It is the indexed fast path for `bcftools mpileup -r`: a
// region query against a whole-genome BAM then reads only the region's
// reads — bounded memory and wall — instead of bucketing every record.
// Records come out in coordinate order (UnionChunks returns the chunks
// ascending and merged, and each chunk is itself coordinate-sorted),
// exactly what the per-chrom bucketing downstream expects.
type bamRegionReader struct {
	f        hfile.SeekHandle
	hdr      *sam.Header
	chunks   []bam.BAIChunk
	ci       int
	bgz      *bgzip.Reader
	br       *sam.BAMReader
	regionOK func(*sam.Record) bool
}

// openBAMRegionReader opens path, reads its header, loads the .csi/.bai
// index and returns a region-restricted reader. It returns (nil, nil) — so
// the caller falls back to a linear scan — when the input is not a
// coordinate-sorted BAM or no index file is present; a non-nil error is a
// genuine read/parse failure.
func openBAMRegionReader(path string, regionStrs []string) (*bamRegionReader, error) {
	f, err := hfile.OpenSeekable(path)
	if err != nil {
		return nil, err
	}
	// sam.NewBAMReader parses the header and fails on a non-BAM (CRAM/SAM)
	// stream, so this also routes those to the linear fallback.
	hr, herr := sam.NewBAMReader(f)
	if herr != nil {
		_ = f.Close()
		return nil, nil
	}
	hdr := hr.Header()
	if !mpileupHeaderIsCoordinateSorted(hdr) {
		_ = f.Close()
		return nil, nil
	}
	resolved, _, rerr := hregion.ResolveRegions(regionStrs, func(name string) int { return hdr.RefIndex(name) })
	if rerr != nil {
		_ = f.Close()
		return nil, rerr
	}
	chunks, ok, cerr := mpileupRegionChunks(path, resolved)
	if cerr != nil {
		_ = f.Close()
		return nil, cerr
	}
	if !ok {
		_ = f.Close()
		return nil, nil // no index — linear fallback
	}
	return &bamRegionReader{
		f:        f,
		hdr:      hdr,
		chunks:   chunks,
		regionOK: mpileupBuildRegionFilter(resolved, hdr),
	}, nil
}

// Header returns the input's SAM header.
func (r *bamRegionReader) Header() *sam.Header { return r.hdr }

// Read returns the next record overlapping the regions, in coordinate
// order, or io.EOF when every chunk is exhausted.
func (r *bamRegionReader) Read() (*sam.Record, error) {
	for {
		if r.br == nil {
			if r.ci >= len(r.chunks) {
				return nil, io.EOF
			}
			c := r.chunks[r.ci]
			r.ci++
			if c.Beg >= c.End {
				continue
			}
			startBlock := int64(c.Beg >> 16)
			if _, err := r.f.Seek(startBlock, io.SeekStart); err != nil {
				return nil, err
			}
			// NewReaderAt(startBlock) makes the reader's virtual offset
			// absolute so the chunk-bounded wrapper stops exactly at c.End.
			bgz, err := bgzip.NewReaderAt(r.f, startBlock)
			if err != nil {
				return nil, err
			}
			// Discard the in-block bytes before the chunk's first record.
			if uoff := int(c.Beg & 0xFFFF); uoff > 0 {
				if _, err := io.CopyN(io.Discard, bgz, int64(uoff)); err != nil {
					_ = bgz.Close()
					return nil, err
				}
			}
			r.bgz = bgz
			r.br = sam.NewBAMBodyReader(&mpileupChunkBoundedReader{r: bgz, end: c.End}, r.hdr)
		}
		rec, err := r.br.Read()
		if err == io.EOF {
			_ = r.bgz.Close()
			r.bgz, r.br = nil, nil
			continue
		}
		if err != nil {
			return nil, err
		}
		// Chunks are coarse; keep only records that actually overlap a region.
		if r.regionOK != nil && !r.regionOK(rec) {
			continue
		}
		return rec, nil
	}
}

// Close releases the open BGZF reader and the file handle.
func (r *bamRegionReader) Close() error {
	if r.bgz != nil {
		_ = r.bgz.Close()
		r.bgz = nil
	}
	if r.f != nil {
		err := r.f.Close()
		r.f = nil
		return err
	}
	return nil
}
