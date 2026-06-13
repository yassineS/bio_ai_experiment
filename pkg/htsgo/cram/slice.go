package cram

import (
	"fmt"
)

// Slice is one parsed CRAM slice: its slice header plus the data blocks
// the slice owns (the CORE block and the external blocks). It is the
// unit a data-series decode operates over — every series of every
// record in the slice is decoded from this block set.
type Slice struct {
	// Header is the parsed slice-header block.
	Header *SliceHeader
	// core is the slice's CORE (bit-packed) data block, or nil when the
	// slice has no CORE block.
	core *Block
	// external maps each external data block's content id to the block.
	external map[int32]*Block
	// major is the CRAM major version of the container the slice belongs
	// to. It is threaded into the seriesSource so EXTERNAL/VARINT integer
	// values are read in the matching format (ITF-8 for v2/v3, uint7 for
	// v4).
	major uint8
}

// DataContainer is a CRAM data container parsed one level deeper than
// the structural Container: its compression header and every slice are
// decoded into Go structures, ready for data-series extraction.
//
// It is produced by ParseDataContainer from a structural Container. The
// first container of a CRAM file holds the SAM header and is not a data
// container; ParseDataContainer rejects it.
type DataContainer struct {
	// Compression is the container's parsed compression header.
	Compression *CompressionHeader
	// Slices holds the container's slices in file order.
	Slices []*Slice
}

// ParseDataContainer interprets a structural Container as a CRAM data
// container: it parses the compression-header block and, for each
// slice, the slice-header block, grouping the slice's data blocks by
// content id. It returns an error for a container that is not a data
// container (for example the file-header container or the EOF marker).
func ParseDataContainer(c *Container) (*DataContainer, error) {
	if c == nil {
		return nil, fmt.Errorf("cram: cannot parse a nil container")
	}
	if len(c.Blocks) == 0 {
		return nil, fmt.Errorf("cram: container %d has no blocks (not a data container)", c.Index)
	}
	first := &c.Blocks[0]
	if first.ContentType != ContentCompressionHeader {
		return nil, fmt.Errorf("cram: container %d's first block is %s, not a compression header",
			c.Index, first.ContentType)
	}
	// The integer encoding throughout the container — compression-header
	// maps, per-encoding parameters, slice headers and data-series values
	// — depends on the CRAM major version: ITF-8/LTF-8 for v2/v3, uint7
	// varints for v4. Containers from Reader.Next carry the version; a
	// hand-built Container leaves Major zero, which we treat as v3+ to
	// preserve the historical default.
	major := c.Major
	if major == 0 {
		major = 3
	}

	payload, err := first.Decompress()
	if err != nil {
		return nil, fmt.Errorf("cram: container %d compression header: %w", c.Index, err)
	}
	ch, err := parseCompressionHeader(payload, major)
	if err != nil {
		return nil, fmt.Errorf("cram: container %d compression header: %w", c.Index, err)
	}
	dc := &DataContainer{Compression: ch}

	// Walk the blocks after the compression header. They alternate: a
	// slice-header block (content type MAPPED_SLICE) followed by that
	// slice's NumBlocks data blocks.
	i := 1
	for i < len(c.Blocks) {
		hdrBlk := &c.Blocks[i]
		if hdrBlk.ContentType != ContentMappedSlice {
			return nil, fmt.Errorf("cram: container %d block %d is %s, expected a slice header",
				c.Index, i, hdrBlk.ContentType)
		}
		hdrPayload, err := hdrBlk.Decompress()
		if err != nil {
			return nil, fmt.Errorf("cram: container %d slice %d header: %w", c.Index, len(dc.Slices), err)
		}
		sh, err := parseSliceHeader(hdrPayload, major)
		if err != nil {
			return nil, fmt.Errorf("cram: container %d slice %d header: %w", c.Index, len(dc.Slices), err)
		}
		i++
		sl := &Slice{Header: sh, external: make(map[int32]*Block), major: major}
		if int(sh.NumBlocks) > len(c.Blocks)-i {
			return nil, fmt.Errorf("cram: container %d slice %d declares %d data blocks but only %d remain",
				c.Index, len(dc.Slices), sh.NumBlocks, len(c.Blocks)-i)
		}
		for j := int32(0); j < sh.NumBlocks; j++ {
			b := &c.Blocks[i]
			i++
			switch b.ContentType {
			case ContentCoreData:
				sl.core = b
			case ContentExternal:
				sl.external[b.ContentID] = b
			default:
				return nil, fmt.Errorf("cram: container %d slice %d data block %d has unexpected content type %s",
					c.Index, len(dc.Slices), j, b.ContentType)
			}
		}
		dc.Slices = append(dc.Slices, sl)
	}
	return dc, nil
}

// newSeriesSource decompresses the slice's CORE block and every external
// block once, returning a seriesSource the slice's data series can be
// decoded through. Decoding more than one series from a slice should
// reuse a single seriesSource so series sharing an external block
// advance the same cursor.
func (sl *Slice) newSeriesSource() (*seriesSource, error) {
	var core *bitReader
	if sl.core != nil {
		payload, err := sl.core.Decompress()
		if err != nil {
			return nil, fmt.Errorf("cram: slice CORE block: %w", err)
		}
		core = newBitReader(payload)
	} else {
		// A slice with no CORE block still needs a (empty) bit reader so
		// an encoding that unexpectedly reaches for the CORE stream
		// errors cleanly rather than dereferencing nil.
		core = newBitReader(nil)
	}
	blocks := make(map[int32][]byte, len(sl.external))
	for id, b := range sl.external {
		payload, err := b.Decompress()
		if err != nil {
			return nil, fmt.Errorf("cram: slice external block (content id %d): %w", id, err)
		}
		blocks[id] = payload
	}
	return &seriesSource{
		core:     core,
		external: make(map[int32]*byteCursor),
		blocks:   blocks,
		reader:   newIntReader(sl.major),
	}, nil
}

// SeriesSource is an opaque handle to a slice's decompressed CORE
// bitstream and external blocks. A single SeriesSource is created per
// slice with Slice.NewSource and threaded through the slice's data
// series decodes so that series sharing an external block advance the
// same read cursor.
type SeriesSource struct {
	s *seriesSource
}

// NewSource decompresses the slice's CORE and external blocks and
// returns a SeriesSource for decoding the slice's data series. Each call
// produces a fresh, independent source positioned at the start of every
// block.
func (sl *Slice) NewSource() (*SeriesSource, error) {
	s, err := sl.newSeriesSource()
	if err != nil {
		return nil, err
	}
	return &SeriesSource{s: s}, nil
}

// HasEmbeddedReference reports whether the slice carries an embedded
// reference block — a per-slice copy of the reference span the slice's
// reads align to, written by samtools' embed_ref mode. A slice with an
// embedded reference needs no external FASTA, REF_CACHE or REF_PATH to
// reconstruct its sequences.
func (sl *Slice) HasEmbeddedReference() bool {
	return sl.Header != nil && sl.Header.EmbeddedRefID >= 0 && sl.external[sl.Header.EmbeddedRefID] != nil
}

// EmbeddedReference decompresses and returns the slice's embedded
// reference bases. The returned span begins at the slice header's
// AlignmentStart and must be at least AlignmentSpan bases long
// (matching htslib's cram_decode_slice embed_ref handling). It returns
// an error if the slice declares an embedded reference id whose block is
// absent or decompresses to a span shorter than the slice covers.
//
// Unlike an external reference, an embedded reference is not MD5-checked
// against the slice header: the bases come from the file itself, so
// there is no separate source to disagree with — htslib likewise trusts
// the embedded block verbatim.
func (sl *Slice) EmbeddedReference() ([]byte, error) {
	if sl.Header == nil || sl.Header.EmbeddedRefID < 0 {
		return nil, fmt.Errorf("cram: slice carries no embedded reference")
	}
	b := sl.external[sl.Header.EmbeddedRefID]
	if b == nil {
		return nil, fmt.Errorf("cram: slice declares embedded reference block id %d but it is absent",
			sl.Header.EmbeddedRefID)
	}
	bases, err := b.Decompress()
	if err != nil {
		return nil, fmt.Errorf("cram: decompressing embedded reference block (content id %d): %w",
			sl.Header.EmbeddedRefID, err)
	}
	if int64(len(bases)) < int64(sl.Header.AlignmentSpan) {
		return nil, fmt.Errorf("cram: embedded reference is %d bases, too small for slice span %d-%d",
			len(bases), sl.Header.AlignmentStart,
			int64(sl.Header.AlignmentStart)+int64(sl.Header.AlignmentSpan)-1)
	}
	return bases, nil
}

// HasSeriesData reports whether the named data series has on-disk data
// in this slice. A series whose encoding draws from an external block
// (EXTERNAL, BYTE_ARRAY_STOP) carries no values when that block is
// absent: CRAM omits a series' block when the series contributes
// nothing. A series read from the CORE bitstream, or one absent from
// the encoding map, always reports true (its presence cannot be
// determined from block membership alone).
func (sl *Slice) HasSeriesData(h *CompressionHeader, key string) bool {
	enc := h.Encoding(key)
	if enc == nil {
		return false
	}
	switch enc.ID {
	case EncodingExternal, EncodingByteArrayStop,
		EncodingVarintUnsigned, EncodingVarintSigned:
		// CRAM v4 VARINT codecs draw from an external block exactly as
		// EXTERNAL does, so their presence is the block's presence.
		return sl.external[enc.ExternalID] != nil
	case EncodingConstByte, EncodingConstInt:
		// A constant series carries no block but always has a (constant)
		// value, so it is always present.
		return true
	case EncodingXPack, EncodingXRLE, EncodingXDelta:
		// A transform codec has data when the block(s) its sub-codec(s)
		// draw from are present.
		return sl.encodingBlocksPresent(enc)
	default:
		return true
	}
}

// encodingBlocksPresent reports whether every external block an encoding
// ultimately reads from is present in the slice. It descends through the
// CRAM v4 transform codecs (XPACK / XRLE / XDELTA) to their wrapped
// sub-encoding(s); a CORE-bitstream or constant leaf is treated as present
// since it needs no external block.
func (sl *Slice) encodingBlocksPresent(enc *Encoding) bool {
	if enc == nil || enc.ID == EncodingNull {
		return true
	}
	switch enc.ID {
	case EncodingExternal, EncodingByteArrayStop,
		EncodingVarintUnsigned, EncodingVarintSigned:
		return sl.external[enc.ExternalID] != nil
	case EncodingConstByte, EncodingConstInt:
		return true
	case EncodingByteArrayLen:
		return sl.encodingBlocksPresent(enc.LenEnc) && sl.encodingBlocksPresent(enc.ValEnc)
	case EncodingXPack, EncodingXDelta:
		return sl.encodingBlocksPresent(enc.SubEnc)
	case EncodingXRLE:
		return sl.encodingBlocksPresent(enc.LenSubEnc) && sl.encodingBlocksPresent(enc.LitSubEnc)
	default:
		// CORE-bitstream encodings (HUFFMAN, BETA, …) carry no external block.
		return true
	}
}

// DecodeIntSeries decodes exactly n integer values of the named
// two-character data series. Pass the slice's record count for the
// series that store one value per record; feature and nested series
// store a data-dependent count and need C4b's record traversal to
// supply n. The source argument is threaded so series sharing an
// external block stay positionally consistent.
func (sl *Slice) DecodeIntSeries(h *CompressionHeader, source *SeriesSource, key string, n int) ([]int32, error) {
	enc := h.Encoding(key)
	if enc == nil {
		return nil, fmt.Errorf("cram: data series %q is not in the compression header's encoding map", key)
	}
	return enc.decodeInts(source.s, n)
}

// DecodeByteArraySeries decodes exactly n byte-array values of the named
// data series.
func (sl *Slice) DecodeByteArraySeries(h *CompressionHeader, source *SeriesSource, key string, n int) ([][]byte, error) {
	enc := h.Encoding(key)
	if enc == nil {
		return nil, fmt.Errorf("cram: data series %q is not in the compression header's encoding map", key)
	}
	return enc.decodeByteArrays(source.s, n)
}

// DrainIntSeries decodes every integer value of the named EXTERNAL data
// series, reading its external block to exhaustion. It is the way to
// extract a whole series without first knowing the value count, and it
// doubles as a correctness check: a clean drain means every byte of the
// block decoded as a well-formed ITF-8 integer. It returns an error for
// a series whose encoding is not EXTERNAL, since only an external block
// is self-delimiting.
func (sl *Slice) DrainIntSeries(h *CompressionHeader, source *SeriesSource, key string) ([]int32, error) {
	enc := h.Encoding(key)
	if enc == nil {
		return nil, fmt.Errorf("cram: data series %q is not in the compression header's encoding map", key)
	}
	return enc.drainInts(source.s)
}

// DrainByteArraySeries decodes every byte-array value of the named
// BYTE_ARRAY_STOP data series, reading its external block to exhaustion.
func (sl *Slice) DrainByteArraySeries(h *CompressionHeader, source *SeriesSource, key string) ([][]byte, error) {
	enc := h.Encoding(key)
	if enc == nil {
		return nil, fmt.Errorf("cram: data series %q is not in the compression header's encoding map", key)
	}
	return enc.drainByteArrays(source.s)
}

// DrainResult is the outcome of draining one data series with
// Slice.DrainSeries. Exactly one of Ints, Bytes or ByteArrays is
// populated according to the series' encoding and value kind; Count is
// the number of values decoded.
type DrainResult struct {
	// Ints holds the values of an integer EXTERNAL series.
	Ints []int32
	// Bytes holds the bytes of a byte-valued EXTERNAL series, one byte
	// per value.
	Bytes []byte
	// ByteArrays holds the values of a BYTE_ARRAY_STOP series.
	ByteArrays [][]byte
	// Count is the number of values decoded (len of whichever field is
	// populated).
	Count int
}

// DrainSeries decodes the named data series of the slice in full,
// reading its external block(s) to exhaustion. It is defined for series
// whose encoding is self-delimiting — EXTERNAL, BYTE_ARRAY_STOP, and
// BYTE_ARRAY_LEN with an EXTERNAL length sub-encoding — and picks the
// integer, raw-byte or byte-array interpretation from the series'
// encoding and the data-series value catalogue. It returns an error,
// never a wrong answer, for a CORE-bitstream encoding (HUFFMAN, BETA,
// GAMMA, SUBEXP, GOLOMB, GOLOMB_RICE), whose values have no per-series
// boundary in the shared block and so need C4b's record traversal to
// bound them.
//
// A clean return is itself the correctness check: every byte of the
// series' block decoded as a well-formed value with nothing left over.
func (sl *Slice) DrainSeries(h *CompressionHeader, source *SeriesSource, key string) (DrainResult, error) {
	enc := h.Encoding(key)
	if enc == nil {
		return DrainResult{}, fmt.Errorf("cram: data series %q is not in the compression header's encoding map", key)
	}
	switch enc.ID {
	case EncodingByteArrayStop:
		ba, err := enc.drainByteArrays(source.s)
		return DrainResult{ByteArrays: ba, Count: len(ba)}, err
	case EncodingByteArrayLen:
		ba, err := enc.drainByteArrayLen(source.s)
		return DrainResult{ByteArrays: ba, Count: len(ba)}, err
	case EncodingExternal, EncodingVarintUnsigned, EncodingVarintSigned:
		if enc.ID == EncodingExternal && SeriesValueKind(key) == SeriesByte {
			b, err := enc.drainRawBytes(source.s)
			return DrainResult{Bytes: b, Count: len(b)}, err
		}
		iv, err := enc.drainInts(source.s)
		return DrainResult{Ints: iv, Count: len(iv)}, err
	case EncodingXPack, EncodingXRLE:
		// CRAM v4 transform codecs expand a self-delimiting external block,
		// so the whole expanded byte series can be drained.
		b, err := enc.drainRawBytes(source.s)
		return DrainResult{Bytes: b, Count: len(b)}, err
	default:
		return DrainResult{}, fmt.Errorf("cram: %s series %q is not drainable; it is read from the shared CORE bitstream and needs a record count", enc.ID, key)
	}
}

// DrainTag decodes every value of the named three-byte tag series of
// the slice in full. Tag series carry per-record auxiliary SAM tag
// values and, like data series, are drainable only when self-delimiting
// (EXTERNAL, BYTE_ARRAY_STOP, or BYTE_ARRAY_LEN with an EXTERNAL length
// sub-encoding). It returns an error for a CORE-bitstream tag encoding
// or an absent block.
func (sl *Slice) DrainTag(h *CompressionHeader, source *SeriesSource, tag string) (DrainResult, error) {
	if len(tag) != 3 {
		return DrainResult{}, fmt.Errorf("cram: tag key %q must be three bytes", tag)
	}
	enc := h.Tags[tagKey{tag[0], tag[1], tag[2]}]
	if enc == nil {
		return DrainResult{}, fmt.Errorf("cram: tag %q is not in the compression header's tag-encoding map", tag)
	}
	switch enc.ID {
	case EncodingByteArrayStop:
		ba, err := enc.drainByteArrays(source.s)
		return DrainResult{ByteArrays: ba, Count: len(ba)}, err
	case EncodingByteArrayLen:
		ba, err := enc.drainByteArrayLen(source.s)
		return DrainResult{ByteArrays: ba, Count: len(ba)}, err
	case EncodingExternal, EncodingVarintUnsigned, EncodingVarintSigned:
		iv, err := enc.drainInts(source.s)
		return DrainResult{Ints: iv, Count: len(iv)}, err
	case EncodingXPack, EncodingXRLE:
		b, err := enc.drainRawBytes(source.s)
		return DrainResult{Bytes: b, Count: len(b)}, err
	default:
		return DrainResult{}, fmt.Errorf("cram: %s tag %q is not drainable from a self-delimiting block", enc.ID, tag)
	}
}

// TagDrainable reports whether the named three-byte tag series can be
// decoded in full with DrainTag.
func (sl *Slice) TagDrainable(h *CompressionHeader, tag string) bool {
	if len(tag) != 3 {
		return false
	}
	enc := h.Tags[tagKey{tag[0], tag[1], tag[2]}]
	if enc == nil {
		return false
	}
	switch enc.ID {
	case EncodingExternal, EncodingByteArrayStop,
		EncodingVarintUnsigned, EncodingVarintSigned:
		return sl.external[enc.ExternalID] != nil
	case EncodingXPack, EncodingXRLE:
		return sl.encodingBlocksPresent(enc)
	case EncodingByteArrayLen:
		if enc.LenEnc == nil || enc.LenEnc.ID != EncodingExternal {
			return false
		}
		if sl.external[enc.LenEnc.ExternalID] == nil {
			return false
		}
		if enc.ValEnc != nil && enc.ValEnc.ID == EncodingExternal {
			return sl.external[enc.ValEnc.ExternalID] != nil
		}
		return true
	default:
		return false
	}
}

// Drainable reports whether the named data series can be decoded in full
// with DrainSeries: true for an EXTERNAL or BYTE_ARRAY_STOP series whose
// external block is present, or a BYTE_ARRAY_LEN series whose EXTERNAL
// length and values blocks are both present; false for a CORE-bitstream
// series or a series whose block is absent.
func (sl *Slice) Drainable(h *CompressionHeader, key string) bool {
	enc := h.Encoding(key)
	if enc == nil {
		return false
	}
	switch enc.ID {
	case EncodingExternal, EncodingByteArrayStop,
		EncodingVarintUnsigned, EncodingVarintSigned:
		return sl.external[enc.ExternalID] != nil
	case EncodingXPack, EncodingXRLE:
		return sl.encodingBlocksPresent(enc)
	case EncodingByteArrayLen:
		if enc.LenEnc == nil || enc.LenEnc.ID != EncodingExternal {
			return false
		}
		if sl.external[enc.LenEnc.ExternalID] == nil {
			return false
		}
		// The values block must be present (or be a CORE encoding,
		// which decodeRawBytes can still pull from).
		if enc.ValEnc != nil && enc.ValEnc.ID == EncodingExternal {
			return sl.external[enc.ValEnc.ExternalID] != nil
		}
		return true
	default:
		return false
	}
}
