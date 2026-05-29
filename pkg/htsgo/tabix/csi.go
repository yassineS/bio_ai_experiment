package tabix

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"sort"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// CSIMagic is the four-byte signature at the start of every CSI file.
var CSIMagic = [4]byte{'C', 'S', 'I', 1}

// Sentinel errors for the CSI parser.
var (
	ErrCSIBadMagic  = errors.New("csi: bad magic")
	ErrCSITruncated = errors.New("csi: truncated input")
)

// CSI is the in-memory representation of a `.csi` index. It mirrors the BAI/
// TBI shape but is parameterised by (MinShift, Depth) so it can address
// chromosomes longer than the BAI 2^29 bp ceiling.
type CSI struct {
	MinShift int32
	Depth    int32
	Aux      []byte
	Names    []string // optional: filled from a tabix-style aux block
	Refs     []CSIRef
	NoCoor   uint64
}

// CSIRef is the per-reference index portion.
type CSIRef struct {
	Bins []CSIBin

	// offBeg/offEnd track the virtual-offset span covered by this
	// reference's records, and nMapped counts records placed here. These
	// feed the meta/pseudo-bin (META_BIN) htslib writes after a reference's
	// data bins. metaDone guards against appending the pseudo-bin twice.
	offBeg   VOffset
	offEnd   VOffset
	nMapped  uint64
	haveSpan bool
	metaDone bool
}

// metaBin is the CSI meta/pseudo-bin number. htslib defines it as
// META_BIN = n_bins + 1 where n_bins = ((1<<(3*depth+3))-1)/7. BinLimit()
// returns exactly that n_bins value, so META_BIN is BinLimit()+1. For the
// default depth of 5 this is 37449+1 = 37450 — identical to the TBI/BAI scheme.
func (c *CSI) metaBin() uint32 { return c.BinLimit() + 1 }

// CSIBin is one bin in a CSI index. LOffset is the linear-index "earliest
// virtual offset of any record reaching this bin" used to short-circuit
// region queries.
type CSIBin struct {
	ID      uint32
	LOffset VOffset
	Chunks  []CSIChunk
}

// CSIChunk is one (virtual offset, virtual offset) record interval.
type CSIChunk struct{ Beg, End VOffset }

// NewCSI returns a freshly initialised CSI. minShift==0 and depth==0 select
// the htslib defaults (14 and 5 — matching BAI bins).
func NewCSI(minShift, depth int32) *CSI {
	if minShift <= 0 {
		minShift = 14
	}
	if depth <= 0 {
		depth = 5
	}
	return &CSI{MinShift: minShift, Depth: depth}
}

// TBXMaxShift mirrors htslib's TBX_MAX_SHIFT: the maximum interval shift the
// tabix CSI scheme addresses before adjusting min_shift.
const TBXMaxShift = 31

// csiDepthForGeneric reproduces htslib tbx_index's choice of n_lvls (CSI
// depth) for a tabix index built with the given positive min_shift when no
// per-reference maximum length is known up front — the situation for the BED
// and other generic presets, whose header lines carry no contig lengths.
// htslib uses n_lvls = max_n_lvls - (min_shift-10)/3 for 10 <= min_shift < 25
// (with max_n_lvls = 9), n_lvls = 9 for min_shift < 10, and n_lvls = 4 for
// min_shift >= 25. With the conventional min_shift of 14 this yields 8.
func csiDepthForGeneric(minShift int32) int32 {
	const maxNLvls = 9
	switch {
	case minShift < 10:
		return maxNLvls
	case minShift < 25:
		return maxNLvls - (minShift-10)/3
	default:
		return 4
	}
}

// MaxPos returns the largest 0-based position addressable at the current
// (MinShift, Depth) combination.
func (c *CSI) MaxPos() int64 { return int64(1) << uint32(c.MinShift+c.Depth*3) }

// BinLimit returns the count of bins implied by (MinShift, Depth): sum of
// 8^k for k = 0..depth.
func (c *CSI) BinLimit() uint32 {
	var n uint32
	for k := int32(0); k <= c.Depth; k++ {
		n += 1 << uint32(3*k)
	}
	return n
}

// Reg2bin returns the smallest CSI bin that fully contains the half-open
// interval [beg, end). The math is the BAI scheme generalised to a
// configurable (MinShift, Depth).
func (c *CSI) Reg2bin(beg, end int64) uint32 {
	end--
	l := c.Depth
	s := uint32(c.MinShift)
	t := uint32(0)
	for k := int32(0); k < c.Depth; k++ {
		t += 1 << uint32(3*k)
	}
	for ; l > 0; l-- {
		if beg>>s == end>>s {
			return t + uint32(beg>>s)
		}
		s += 3
		t -= 1 << uint32(3*l-3)
	}
	return 0
}

// Reg2bins returns every bin number whose span overlaps [beg, end).
func (c *CSI) Reg2bins(beg, end int64) []uint32 {
	if end <= beg {
		return []uint32{0}
	}
	end--
	out := []uint32{0}
	s := uint32(c.MinShift + c.Depth*3)
	t := uint32(0)
	for l := int32(1); l <= c.Depth; l++ {
		s -= 3
		t += 1 << uint32(3*(l-1))
		b := t + uint32(beg>>s)
		e := t + uint32(end>>s)
		for k := b; k <= e; k++ {
			out = append(out, k)
		}
	}
	return out
}

// AddRecord registers a record under (refID, beg, end). It is the imperative
// counterpart of a CSI Build for streaming callers.
func (c *CSI) AddRecord(refID int, beg, end int64, vbeg, vend VOffset) {
	for len(c.Refs) <= refID {
		c.Refs = append(c.Refs, CSIRef{})
	}
	ref := &c.Refs[refID]
	// Track the per-reference virtual-offset span and mapped-record count
	// for the meta/pseudo-bin emitted by Finalize.
	if !ref.haveSpan {
		ref.offBeg = vbeg
		ref.haveSpan = true
	}
	ref.offEnd = vend
	ref.nMapped++
	binID := c.Reg2bin(beg, end)
	var bin *CSIBin
	for i := range ref.Bins {
		if ref.Bins[i].ID == binID {
			bin = &ref.Bins[i]
			break
		}
	}
	if bin == nil {
		ref.Bins = append(ref.Bins, CSIBin{ID: binID, LOffset: vbeg})
		bin = &ref.Bins[len(ref.Bins)-1]
	}
	// Track the minimum begin offset reaching this bin. Do NOT treat a
	// stored LOffset of 0 as "unset": 0 is a legitimate virtual offset (the
	// first byte of the first block), and htslib keeps it as the loff for a
	// bin whose first record starts at the file head.
	if vbeg < bin.LOffset {
		bin.LOffset = vbeg
	}
	if n := len(bin.Chunks); n > 0 && bin.Chunks[n-1].End >= vbeg {
		if vend > bin.Chunks[n-1].End {
			bin.Chunks[n-1].End = vend
		}
	} else {
		bin.Chunks = append(bin.Chunks, CSIChunk{Beg: vbeg, End: vend})
	}
}

// Finalize appends each non-empty reference's meta/pseudo-bin, mirroring
// htslib's hts_idx_finish. The pseudo-bin (id metaBin()) carries two chunks:
// {off_beg, off_end} (the reference's virtual-offset span) and
// {n_mapped, n_unmapped} (record counts, stored verbatim in the chunk's
// begin/end slots). Per htslib update_loff, a CSI bin with id >= n_bins (the
// meta-bin) has loff == 0. Finalize is idempotent per reference.
func (c *CSI) Finalize() {
	for r := range c.Refs {
		ref := &c.Refs[r]
		if !ref.haveSpan || ref.metaDone {
			continue
		}
		ref.metaDone = true
		ref.Bins = append(ref.Bins, CSIBin{
			ID:      c.metaBin(),
			LOffset: 0,
			Chunks: []CSIChunk{
				{Beg: ref.offBeg, End: ref.offEnd},
				{Beg: VOffset(ref.nMapped), End: VOffset(0)},
			},
		})
	}
}

// Write serialises c to w using the CSI binary layout. The output is little-
// endian and matches htslib byte ordering.
func (c *CSI) Write(w io.Writer) error {
	bw := bufio.NewWriter(w)
	if _, err := bw.Write(CSIMagic[:]); err != nil {
		return err
	}
	if err := binary.Write(bw, binary.LittleEndian, c.MinShift); err != nil {
		return err
	}
	if err := binary.Write(bw, binary.LittleEndian, c.Depth); err != nil {
		return err
	}
	if err := binary.Write(bw, binary.LittleEndian, int32(len(c.Aux))); err != nil {
		return err
	}
	if len(c.Aux) > 0 {
		if _, err := bw.Write(c.Aux); err != nil {
			return err
		}
	}
	if err := binary.Write(bw, binary.LittleEndian, int32(len(c.Refs))); err != nil {
		return err
	}
	for i := range c.Refs {
		ref := c.Refs[i]
		sort.Slice(ref.Bins, func(a, b int) bool { return ref.Bins[a].ID < ref.Bins[b].ID })
		if err := binary.Write(bw, binary.LittleEndian, int32(len(ref.Bins))); err != nil {
			return err
		}
		for _, bin := range ref.Bins {
			if err := binary.Write(bw, binary.LittleEndian, bin.ID); err != nil {
				return err
			}
			if err := binary.Write(bw, binary.LittleEndian, uint64(bin.LOffset)); err != nil {
				return err
			}
			if err := binary.Write(bw, binary.LittleEndian, int32(len(bin.Chunks))); err != nil {
				return err
			}
			for _, ch := range bin.Chunks {
				if err := binary.Write(bw, binary.LittleEndian, uint64(ch.Beg)); err != nil {
					return err
				}
				if err := binary.Write(bw, binary.LittleEndian, uint64(ch.End)); err != nil {
					return err
				}
			}
		}
	}
	// htslib's idx_save_core always writes the trailing n_no_coor (uint64)
	// for CSI indexes, even when zero.
	if err := binary.Write(bw, binary.LittleEndian, c.NoCoor); err != nil {
		return err
	}
	return bw.Flush()
}

// WriteFile serialises c to path inside a BGZF stream (the canonical CSI
// on-disk form).
func (c *CSI) WriteFile(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bgzip.NewWriter(f)
	if err := c.Write(bw); err != nil {
		return err
	}
	return bw.Close()
}

// ReadCSI parses a CSI binary stream (already BGZF-decoded) into a fresh CSI.
func ReadCSI(r io.Reader) (*CSI, error) {
	br := bufio.NewReader(r)
	var magic [4]byte
	if _, err := io.ReadFull(br, magic[:]); err != nil {
		return nil, csiWrapEOF(err)
	}
	if magic != CSIMagic {
		return nil, fmt.Errorf("%w: got %v", ErrCSIBadMagic, magic)
	}
	c := &CSI{}
	if err := binary.Read(br, binary.LittleEndian, &c.MinShift); err != nil {
		return nil, csiWrapEOF(err)
	}
	if err := binary.Read(br, binary.LittleEndian, &c.Depth); err != nil {
		return nil, csiWrapEOF(err)
	}
	var lAux int32
	if err := binary.Read(br, binary.LittleEndian, &lAux); err != nil {
		return nil, csiWrapEOF(err)
	}
	if lAux > 0 {
		c.Aux = make([]byte, lAux)
		if _, err := io.ReadFull(br, c.Aux); err != nil {
			return nil, csiWrapEOF(err)
		}
		c.Names = parseCSIAuxNames(c.Aux)
	}
	var nRef int32
	if err := binary.Read(br, binary.LittleEndian, &nRef); err != nil {
		return nil, csiWrapEOF(err)
	}
	if nRef < 0 {
		return nil, ErrCSIBadMagic
	}
	c.Refs = make([]CSIRef, nRef)
	for r := int32(0); r < nRef; r++ {
		var nBin int32
		if err := binary.Read(br, binary.LittleEndian, &nBin); err != nil {
			return nil, csiWrapEOF(err)
		}
		bins := make([]CSIBin, nBin)
		for i := range bins {
			var lo uint64
			if err := binary.Read(br, binary.LittleEndian, &bins[i].ID); err != nil {
				return nil, csiWrapEOF(err)
			}
			if err := binary.Read(br, binary.LittleEndian, &lo); err != nil {
				return nil, csiWrapEOF(err)
			}
			bins[i].LOffset = VOffset(lo)
			var nChunk int32
			if err := binary.Read(br, binary.LittleEndian, &nChunk); err != nil {
				return nil, csiWrapEOF(err)
			}
			chunks := make([]CSIChunk, nChunk)
			for j := range chunks {
				var beg, end uint64
				if err := binary.Read(br, binary.LittleEndian, &beg); err != nil {
					return nil, csiWrapEOF(err)
				}
				if err := binary.Read(br, binary.LittleEndian, &end); err != nil {
					return nil, csiWrapEOF(err)
				}
				chunks[j] = CSIChunk{Beg: VOffset(beg), End: VOffset(end)}
			}
			bins[i].Chunks = chunks
		}
		c.Refs[r] = CSIRef{Bins: bins}
	}
	var trailer uint64
	if err := binary.Read(br, binary.LittleEndian, &trailer); err == nil {
		c.NoCoor = trailer
	}
	return c, nil
}

// ReadCSIFile is the on-disk counterpart of ReadCSI.
func ReadCSIFile(path string) (*CSI, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	br, err := bgzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer br.Close()
	return ReadCSI(br)
}

// SetAuxFromTabix encodes a tabix-style auxiliary block: 6 int32 config
// values followed by an int32 l_nm and a NUL-terminated reference-name list.
// We use this so CSI indexes built from .vcf.gz or .bcf can carry the same
// column-config payload `.tbi` files do.
func (c *CSI) SetAuxFromTabix(cfg Config, names []string) {
	var buf bytes.Buffer
	for _, v := range []int32{cfg.Format, cfg.ColSeq, cfg.ColBeg, cfg.ColEnd, cfg.Meta, cfg.Skip} {
		binary.Write(&buf, binary.LittleEndian, v)
	}
	var nmBuf bytes.Buffer
	for _, n := range names {
		nmBuf.WriteString(n)
		nmBuf.WriteByte(0)
	}
	binary.Write(&buf, binary.LittleEndian, int32(nmBuf.Len()))
	buf.Write(nmBuf.Bytes())
	c.Aux = buf.Bytes()
	c.Names = append(c.Names[:0], names...)
}

// parseCSIAuxNames extracts the trailing reference-name list from the tabix-
// style aux block. If the block doesn't look like the tabix layout we return
// an empty slice and let the caller fall back to a higher-level dictionary.
func parseCSIAuxNames(aux []byte) []string {
	// 6 int32 config + l_nm int32 + names → at least 28 bytes.
	if len(aux) < 28 {
		return nil
	}
	off := 24 // skip 6 int32 config values
	lNm := int32(binary.LittleEndian.Uint32(aux[off : off+4]))
	off += 4
	if lNm < 0 || off+int(lNm) > len(aux) {
		return nil
	}
	nameBuf := aux[off : off+int(lNm)]
	var out []string
	start := 0
	for i := 0; i < len(nameBuf); i++ {
		if nameBuf[i] == 0 {
			out = append(out, string(nameBuf[start:i]))
			start = i + 1
		}
	}
	return out
}

// RegionChunks returns the merged chunks that may contain records overlapping
// [beg, end) on refID. The linear-index trim that BAI does is replaced by the
// per-bin LOffset minimum: skip every chunk whose End ≤ min(LOffset across
// bins overlapping the query).
func (c *CSI) RegionChunks(refID int, beg, end int64) []CSIChunk {
	if refID < 0 || refID >= len(c.Refs) || end <= beg {
		return nil
	}
	ref := c.Refs[refID]
	bins := c.Reg2bins(beg, end)
	binSet := make(map[uint32]struct{}, len(bins))
	for _, b := range bins {
		binSet[b] = struct{}{}
	}
	var minOff VOffset
	have := false
	for _, b := range ref.Bins {
		if _, ok := binSet[b.ID]; !ok {
			continue
		}
		if !have || b.LOffset < minOff {
			minOff = b.LOffset
			have = true
		}
	}
	var chunks []CSIChunk
	for _, bin := range ref.Bins {
		if _, ok := binSet[bin.ID]; !ok {
			continue
		}
		for _, ch := range bin.Chunks {
			if ch.End <= minOff {
				continue
			}
			if ch.Beg < minOff {
				ch.Beg = minOff
			}
			chunks = append(chunks, ch)
		}
	}
	if len(chunks) == 0 {
		return nil
	}
	sort.Slice(chunks, func(a, b int) bool { return chunks[a].Beg < chunks[b].Beg })
	merged := chunks[:0]
	cur := chunks[0]
	for i := 1; i < len(chunks); i++ {
		if chunks[i].Beg <= cur.End {
			if chunks[i].End > cur.End {
				cur.End = chunks[i].End
			}
		} else {
			merged = append(merged, cur)
			cur = chunks[i]
		}
	}
	merged = append(merged, cur)
	return merged
}

// ChromIDInCSI returns the dictionary index of name in c.Names, or -1 if
// the name is unknown.
func ChromIDInCSI(c *CSI, name string) int {
	for i, n := range c.Names {
		if n == name {
			return i
		}
	}
	return -1
}

// BuildCSIFromDataFile streams a bgzipped tab-delimited file (preset-config
// driven), records per-record virtual offsets, and returns a populated CSI
// index ready to be written or queried. It mirrors the existing TBI Build
// path with the bin scheme parameterised by minShift.
func BuildCSIFromDataFile(path string, cfg Config, minShift int32) (*CSI, error) {
	if err := validateConfig(cfg); err != nil {
		return nil, err
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	offsets, err := bgzip.Scan(f)
	if err != nil && !errors.Is(err, bgzip.ErrTruncated) {
		return nil, err
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		return nil, err
	}
	br, err := bgzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer br.Close()
	data, err := io.ReadAll(br)
	if err != nil {
		return nil, err
	}

	uoffToV := func(pos int64) VOffset { return VOffsetAt(offsets, pos) }

	// Match htslib tbx_index's CSI parameter choice: a positive min_shift
	// selects CSI, and for presets without header-supplied contig lengths
	// (BED and other generic formats) the depth is derived from min_shift.
	ms := minShift
	if ms <= 0 {
		ms = 14
	}
	csi := NewCSI(ms, csiDepthForGeneric(ms))
	idx := NewIndex(cfg) // we reuse the existing parser bits
	skip := int(cfg.Skip)
	scanner := bufio.NewScanner(bytes.NewReader(data))
	scanner.Buffer(make([]byte, 0, 1<<16), 1<<24)
	var pos int64
	var lineNo int
	var fieldBuf [64][]byte
	for scanner.Scan() {
		raw := scanner.Bytes()
		lineStart := pos
		pos += int64(len(raw)) + 1
		lineNo++
		if len(raw) == 0 {
			continue
		}
		if lineNo <= skip {
			continue
		}
		if raw[0] == byte(cfg.Meta) {
			continue
		}
		fields := splitFields(raw, fieldBuf[:0])
		if int(cfg.ColSeq) < 1 || int(cfg.ColSeq) > len(fields) {
			return nil, ErrBadRecord
		}
		chrom := string(fields[cfg.ColSeq-1])
		if chrom == "" || chrom == "*" {
			csi.NoCoor++
			continue
		}
		beg, end, err := idx.recordEnd(fields)
		if err != nil {
			return nil, fmt.Errorf("tabix: line %d: %w", lineNo, err)
		}
		refID := idx.ensureRef(chrom)
		v := uoffToV(lineStart)
		vEnd := uoffToV(pos)
		csi.AddRecord(refID, int64(beg), int64(end), v, vEnd)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	// Carry the contig name list in the aux block (tabix-style).
	csi.SetAuxFromTabix(cfg, idx.Names)
	csi.Finalize()
	return csi, nil
}

func csiWrapEOF(err error) error {
	if err == io.EOF || err == io.ErrUnexpectedEOF {
		return ErrCSITruncated
	}
	return err
}
