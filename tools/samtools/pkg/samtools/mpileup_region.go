package samtools

import (
	"bytes"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// mpileupRegionChunks loads the BAM's .csi (preferred) or .bai index and returns
// the merged, coordinate-sorted file chunks covering the resolved regions, or
// ok=false when neither sibling index exists (the caller then linear-scans).
func mpileupRegionChunks(path string, regions []region.ResolvedRegion) (chunks []bam.BAIChunk, ok bool, err error) {
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

// bamRegionReader is a sam.Reader that yields only the records overlapping a set
// of regions by seeking the BAM's index chunks rather than scanning the whole
// file. It is the indexed fast path for `samtools mpileup -r`: a region query
// against a whole-genome BAM then reads only the region's reads — bounded memory
// and wall — instead of linearly churning through every record. Records come out
// in coordinate order (UnionChunks returns the chunks ascending and merged, and
// each chunk is itself coordinate-sorted), exactly what the streaming pileup
// expects.
type bamRegionReader struct {
	f        hfile.SeekHandle
	hdr      *sam.Header
	chunks   []bam.BAIChunk
	ci       int
	bgz      *bgzip.Reader
	br       *sam.BAMReader
	regionOK func(*sam.Record) bool
}

// openBAMRegionReader opens path, reads its header, loads the .csi/.bai index
// and returns a region-restricted reader. It returns (nil, nil) — so the caller
// falls back to a linear scan — when the input is not a coordinate-sorted BAM or
// no index file is present; a non-nil error is a genuine read/parse failure.
func openBAMRegionReader(path string, regionStrs []string) (*bamRegionReader, error) {
	f, err := openSeekable(path)
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
	// A present .csi/.bai index implies the file is physically coordinate-sorted
	// (samtools index refuses to index anything else), and the index yields the
	// region's records in coordinate order regardless of the @HD SO tag. So we
	// accept a mis-tagged SO:unsorted (or missing SO) BAM here too — the common
	// real-world case — rather than falling back to a whole-file linear drain.
	// SO:queryname (and other non-coordinate values) are still rejected.
	if !headerConsensusCanStream(hdr) {
		_ = f.Close()
		return nil, nil
	}
	resolved, _, rerr := region.ResolveRegions(regionStrs, func(name string) int { return hdr.RefIndex(name) })
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
		regionOK: buildRegionFilter(resolved, hdr),
	}, nil
}

// Header returns the input's SAM header.
func (r *bamRegionReader) Header() *sam.Header { return r.hdr }

// Read returns the next record overlapping the regions, in coordinate order, or
// io.EOF when every chunk is exhausted.
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
			// NewReaderAt(startBlock) makes the reader's virtual offset absolute so
			// the chunk-bounded wrapper stops exactly at c.End.
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
			r.br = sam.NewBAMBodyReader(&chunkBoundedReader{r: bgz, end: c.End}, r.hdr)
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
