package cram

import (
	"bufio"
	"io"
	"os"
)

// countingReader wraps an io.Reader and tracks the total number of bytes
// it has yielded. It is used by the .crai index builder to record the
// absolute byte offset of every container as the CRAM stream is walked.
type countingReader struct {
	r io.Reader
	n int64
}

// Read implements io.Reader, advancing the byte counter.
func (c *countingReader) Read(p []byte) (int, error) {
	n, err := c.r.Read(p)
	c.n += int64(n)
	return n, err
}

// offset returns the number of bytes consumed from the stream so far.
func (c *countingReader) offset() int64 { return c.n }

// BuildCRAI walks the CRAM stream in r and returns one CRAIEntry per
// slice, in file order. Each entry records the slice's reference span
// (taken from the slice header) and the byte offsets a seeking reader
// needs: the absolute offset of the container that holds the slice, the
// slice-header block's offset within the container body, and the byte
// size of the slice's blocks (the slice-header block plus its data
// blocks).
//
// BuildCRAI is the index-building counterpart of ReadCRAI/WriteCRAI:
// WriteCRAI(BuildCRAI(cram)) produces the .crai a samtools-style index
// command writes. The file-header container and the EOF marker carry no
// alignment slices and contribute no entries.
func BuildCRAI(r io.Reader) ([]CRAIEntry, error) {
	cr := &countingReader{r: r}
	br := bufio.NewReader(cr)
	def, err := readFileDefinition(br)
	if err != nil {
		return nil, err
	}
	// readFileDefinition consumed the 26-byte file definition; the first
	// container begins where the buffered reader's cursor now sits.
	var entries []CRAIEntry
	for {
		// The absolute offset of the next container is the bytes the
		// counting reader has handed out, minus whatever the bufio.Reader
		// has buffered but not yet served to a parser.
		containerOffset := cr.offset() - int64(br.Buffered())
		if _, err := br.Peek(1); err == io.EOF {
			// A well-formed CRAM file ends with the EOF marker, which the
			// loop below consumes; reaching a clean stream end here means
			// the marker was already seen and the file is exhausted.
			return entries, nil
		}
		hdr, err := readContainerHeader(br, def)
		if err == io.EOF {
			return entries, nil
		}
		if err != nil {
			return nil, err
		}
		if hdr.IsEOF {
			return entries, nil
		}
		// The container body holds NumBlocks blocks. Walk them, recording
		// the per-slice index entries. A slice is a slice-header block
		// followed by that slice's data blocks.
		blockEnds, err := craiContainerBlocks(br, def, hdr)
		if err != nil {
			return nil, err
		}
		entries = append(entries, craiSliceEntries(blockEnds, containerOffset)...)
	}
}

// craiBlock is one block walked while building a .crai: its parsed form
// plus the byte offset, relative to the start of the container body, at
// which the block begins.
type craiBlock struct {
	block       Block
	bodyOffset  int64 // offset of the block within the container body.
	encodedSize int64 // the block's total on-disk byte size.
}

// craiContainerBlocks reads every block of one container body, recording
// each block's body-relative offset and on-disk size. It is the index
// builder's structural walk: it mirrors Reader.Next's block loop but
// keeps the per-block offsets the .crai needs.
func craiContainerBlocks(br *bufio.Reader, def FileDefinition, hdr ContainerHeader) ([]craiBlock, error) {
	bodyReader := &countingReader{r: io.LimitReader(br, int64(hdr.Length))}
	limited := &io.LimitedReader{R: bodyReader, N: int64(hdr.Length)}
	var blocks []craiBlock
	for i := int32(0); i < hdr.NumBlocks; i++ {
		start := bodyReader.offset()
		blk, err := readBlock(limited, def, limited.N)
		if err != nil {
			return nil, wrapf(err, "container block %d", i)
		}
		blocks = append(blocks, craiBlock{
			block:       blk,
			bodyOffset:  start,
			encodedSize: bodyReader.offset() - start,
		})
	}
	return blocks, nil
}

// craiSliceEntries turns the walked blocks of one container into the
// .crai entries for its slices. The container body is a compression-
// header block followed, per slice, by a slice-header block and that
// slice's NumBlocks data blocks. Each slice yields one entry: its
// reference span from the slice header, the slice-header block's
// body-relative offset as SliceOffset, and the summed byte size of the
// slice-header block and its data blocks as SliceSize.
func craiSliceEntries(blocks []craiBlock, containerOffset int64) []CRAIEntry {
	var entries []CRAIEntry
	for i := 0; i < len(blocks); i++ {
		b := &blocks[i]
		if b.block.ContentType != ContentMappedSlice {
			continue // the compression-header block and stray blocks.
		}
		payload, err := b.block.Decompress()
		if err != nil {
			continue // an unreadable slice header contributes no entry.
		}
		sh, err := parseSliceHeader(payload)
		if err != nil {
			continue
		}
		sliceSize := b.encodedSize
		// The slice owns the next NumBlocks blocks; sum their sizes.
		for j := 0; j < int(sh.NumBlocks) && i+1+j < len(blocks); j++ {
			sliceSize += blocks[i+1+j].encodedSize
		}
		entries = append(entries, CRAIEntry{
			RefID:           sh.RefSeqID,
			AlignmentStart:  sh.AlignmentStart,
			AlignmentSpan:   sh.AlignmentSpan,
			ContainerOffset: containerOffset,
			SliceOffset:     b.bodyOffset,
			SliceSize:       sliceSize,
		})
	}
	return entries
}

// BuildCRAIFile walks the named CRAM file and returns its .crai entries.
func BuildCRAIFile(path string) ([]CRAIEntry, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return BuildCRAI(f)
}

// CreateCRAI walks the CRAM file at cramPath and writes its .crai index
// to craiPath. It is the file-level convenience wrapper a samtools-style
// `index` command uses for CRAM input.
func CreateCRAI(cramPath, craiPath string) error {
	entries, err := BuildCRAIFile(cramPath)
	if err != nil {
		return err
	}
	f, err := os.Create(craiPath)
	if err != nil {
		return err
	}
	if err := WriteCRAI(f, entries); err != nil {
		f.Close()
		os.Remove(craiPath)
		return err
	}
	return f.Close()
}
