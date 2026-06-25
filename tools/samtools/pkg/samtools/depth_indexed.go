package samtools

import (
	"bytes"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// DepthFile runs depth over one or more on-disk inputs. When the caller
// requested one or more -r regions and an input is a seekable BAM with a
// sibling .csi / .bai index, it seeks straight to the BGZF chunks the regions
// overlap and feeds only those records into the streaming engine — inflating a
// small fraction of the file instead of the whole genome. This mirrors the
// indexed seek-and-scan that samtools view (ViewFile) and upstream bam2depth
// already use, and is what closes the large gap against upstream on a
// single-chromosome query.
//
// When no regions are requested, an input is not a local indexed BAM (CRAM,
// SAM, stdin, a BAM with no sibling index), or the index cannot be read, it
// falls back to the linear streaming Depth over the supplied readers, so the
// emitted depth is identical to the non-indexed path in every case. The
// per-position output is byte-for-byte the same as a whole-file scan because
// the streaming engine still clamps every position to the requested intervals
// (inInterval / refBeg / refEnd); the index only changes which BGZF blocks are
// inflated, never which positions are emitted.
//
// paths and fallback are positional-parallel: fallback[i] is the already-open
// reader for paths[i], used when paths[i] cannot take the indexed path.
func DepthFile(paths []string, fallback []io.Reader, out io.Writer, opts DepthOptions) error {
	if len(paths) == 0 {
		return fmt.Errorf("samtools depth: no inputs")
	}
	if len(paths) != len(fallback) {
		return fmt.Errorf("samtools depth: internal: %d paths but %d readers", len(paths), len(fallback))
	}
	// The indexed path only applies to an -r region query (not -b BED, not -A
	// all-trans, not the plain whole-file scan). Anything else streams.
	if len(opts.Regions) == 0 || opts.BedPath != "" {
		return Depth(fallback, out, opts)
	}

	readers := make([]sam.Reader, len(paths))
	handles := make([]io.Closer, 0, len(paths))
	cleanup := func() {
		for _, c := range handles {
			_ = c.Close()
		}
	}

	for i, path := range paths {
		rd, hs, ok, err := openIndexedDepthReader(path, opts)
		if err != nil {
			cleanup()
			return err
		}
		if !ok {
			// This input cannot take the indexed path; fall back to a full
			// streaming run over every reader so the merge stays simple and the
			// output is identical. (Mixing indexed and streaming inputs in one
			// merge is never needed in practice — depth inputs share a layout.)
			cleanup()
			return Depth(fallback, out, opts)
		}
		readers[i] = rd
		handles = append(handles, hs...)
	}
	defer cleanup()
	return depthFromReaders(readers, out, opts)
}

// openIndexedDepthReader tries to build a region-bounded sam.Reader over the
// indexed BAM at path. ok is false (with no error) when path is not a local
// indexed BAM and the caller should fall back to streaming. The returned
// closers must be closed when the scan is done.
func openIndexedDepthReader(path string, opts DepthOptions) (sam.Reader, []io.Closer, bool, error) {
	if path == "" || path == "-" {
		return nil, nil, false, nil
	}
	// CRAM uses a different (container) seek model; leave it to the streaming
	// path, which already decodes CRAM correctly.
	if isCRAM, _ := pathIsCRAM(path); isCRAM {
		return nil, nil, false, nil
	}

	// Prefer .csi (wider coordinate range) over .bai, matching ViewFile.
	csiPath := path + ".csi"
	if csiBytes, err := hfile.ReadFile(csiPath); err == nil {
		idx, ierr := bam.ReadCSI(bytes.NewReader(csiBytes))
		if ierr != nil {
			return nil, nil, false, fmt.Errorf("samtools depth: read %s: %w", csiPath, ierr)
		}
		return newChunkScanReaderCSI(path, idx, opts)
	}
	baiPath := path + ".bai"
	baiBytes, err := hfile.ReadFile(baiPath)
	if err != nil {
		// No sibling index — stream.
		return nil, nil, false, nil
	}
	idx, ierr := bam.ReadBAI(strings.NewReader(string(baiBytes)))
	if ierr != nil {
		return nil, nil, false, fmt.Errorf("samtools depth: read %s: %w", baiPath, ierr)
	}
	return newChunkScanReaderBAI(path, idx, opts)
}

func newChunkScanReaderBAI(path string, idx *bam.BAIIndex, opts DepthOptions) (sam.Reader, []io.Closer, bool, error) {
	return newChunkScanReader(path, opts, func(hdr *sam.Header, resolved []region.ResolvedRegion) []bam.BAIChunk {
		return bam.UnionChunks(idx, resolved)
	})
}

func newChunkScanReaderCSI(path string, idx *bam.CSIIndex, opts DepthOptions) (sam.Reader, []io.Closer, bool, error) {
	return newChunkScanReader(path, opts, func(hdr *sam.Header, resolved []region.ResolvedRegion) []bam.BAIChunk {
		return bam.UnionChunksCSI(idx, resolved)
	})
}

// newChunkScanReader opens path, reads its header, resolves the requested
// regions to a merged BGZF chunk list, and returns a sam.Reader that yields the
// records inside those chunks in coordinate order. It returns ok=false when no
// region resolves to a known reference (the caller then streams).
func newChunkScanReader(path string, opts DepthOptions, unionFn func(*sam.Header, []region.ResolvedRegion) []bam.BAIChunk) (sam.Reader, []io.Closer, bool, error) {
	f, err := openSeekable(path)
	if err != nil {
		return nil, nil, false, nil
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		_ = f.Close()
		return nil, nil, false, err
	}
	hdrReader, err := sam.NewBAMReader(f)
	if err != nil {
		// Not a BAM (or unreadable header) — let the streaming path handle it.
		_ = f.Close()
		return nil, nil, false, nil
	}
	hdr := hdrReader.Header()
	resolved, _, perr := region.ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if perr != nil {
		_ = f.Close()
		return nil, nil, false, perr
	}
	chunks := unionFn(hdr, resolved)
	r := &chunkScanReader{
		f:      f,
		hdr:    hdr,
		chunks: chunks,
		idx:    -1,
	}
	return r, []io.Closer{f}, true, nil
}

// chunkScanReader is a sam.Reader that walks a sorted, merged list of BGZF
// virtual-offset chunks, decoding the records inside them in coordinate order.
// It is the seek-and-scan counterpart of the linear depthSource and exposes the
// same depth-tailored ReadDepthInto fast decode, so the engine reads identical
// records either way — only fewer BGZF blocks are inflated.
type chunkScanReader struct {
	f      hfile.SeekHandle
	hdr    *sam.Header
	chunks []bam.BAIChunk

	idx  int            // index of the chunk currently being scanned, -1 before first
	bgz  *bgzip.Reader  // BGZF reader positioned inside the current chunk
	body *sam.BAMReader // record decoder over the bounded current chunk
}

// Header returns the BAM header parsed from the file.
func (c *chunkScanReader) Header() *sam.Header { return c.hdr }

// Read returns the next record across all chunks, or io.EOF when exhausted.
func (c *chunkScanReader) Read() (*sam.Record, error) {
	for {
		if c.body == nil {
			if !c.nextChunk() {
				return nil, io.EOF
			}
		}
		rec, err := c.body.Read()
		if err == io.EOF {
			c.closeChunk()
			continue
		}
		if err != nil {
			return nil, err
		}
		return rec, nil
	}
}

// ReadInto decodes the next record into dst without allocating, advancing
// across chunk boundaries as needed.
func (c *chunkScanReader) ReadInto(dst *sam.Record) error {
	for {
		if c.body == nil {
			if !c.nextChunk() {
				return io.EOF
			}
		}
		err := c.body.ReadInto(dst)
		if err == io.EOF {
			c.closeChunk()
			continue
		}
		return err
	}
}

// ReadDepthInto decodes only the depth-relevant fields of the next record into
// dst (the same fast path depthSource uses on the linear scan), advancing
// across chunk boundaries as needed.
func (c *chunkScanReader) ReadDepthInto(dst *sam.Record, needQual bool) error {
	for {
		if c.body == nil {
			if !c.nextChunk() {
				return io.EOF
			}
		}
		err := c.body.ReadDepthInto(dst, needQual)
		if err == io.EOF {
			c.closeChunk()
			continue
		}
		return err
	}
}

// nextChunk seeks to the next non-empty chunk and primes a record decoder over
// it. It returns false when there are no more chunks.
func (c *chunkScanReader) nextChunk() bool {
	for {
		c.idx++
		if c.idx >= len(c.chunks) {
			return false
		}
		ch := c.chunks[c.idx]
		if ch.Beg >= ch.End {
			continue
		}
		startBlock := int64(ch.Beg >> 16)
		if _, err := c.f.Seek(startBlock, io.SeekStart); err != nil {
			// A seek failure leaves body nil so Read surfaces EOF; the caller's
			// stream still produces correct (if region-empty) output. Treat it
			// as no-more-chunks.
			return false
		}
		// NewReaderAt(startBlock) so VirtualOffset is absolute and the
		// chunk-bounded wrapper stops at the right place (see view.go / bgzf).
		bgz, err := bgzip.NewReaderAt(c.f, startBlock)
		if err != nil {
			return false
		}
		// Skip the in-block bytes before the first wanted record.
		if uoff := int(ch.Beg & 0xFFFF); uoff > 0 {
			if _, err := io.CopyN(io.Discard, bgz, int64(uoff)); err != nil {
				_ = bgz.Close()
				return false
			}
		}
		c.bgz = bgz
		c.body = sam.NewBAMBodyReader(&chunkBoundedReader{r: bgz, end: ch.End}, c.hdr)
		return true
	}
}

// closeChunk releases the current chunk's BGZF reader and clears the decoder so
// nextChunk advances on the following Read.
func (c *chunkScanReader) closeChunk() {
	if c.bgz != nil {
		_ = c.bgz.Close()
		c.bgz = nil
	}
	c.body = nil
}
