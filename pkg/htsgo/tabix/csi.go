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
type CSIRef struct{ Bins []CSIBin }

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
	if vbeg < bin.LOffset || bin.LOffset == 0 {
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
	if c.NoCoor > 0 {
		if err := binary.Write(bw, binary.LittleEndian, c.NoCoor); err != nil {
			return err
		}
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
	return ReadCSIGz(f)
}

// ReadCSIGz parses a BGZF-compressed .csi index read from r. It is the
// reader-based counterpart of ReadCSIFile, used when the index bytes come from
// a source other than the local filesystem (for example a remote sibling index
// downloaded through hfile.ReadFile and wrapped in a bytes.Reader).
func ReadCSIGz(r io.Reader) (*CSI, error) {
	br, err := bgzip.NewReader(r)
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

// ConfigFromAux reconstructs the column Config stored in a tabix-style CSI
// aux block (the 6 leading int32 values written by SetAuxFromTabix). It
// returns the zero Config and false when the aux block is too short to hold
// a config header.
func (c *CSI) ConfigFromAux() (Config, bool) {
	if len(c.Aux) < 24 {
		return Config{}, false
	}
	g := func(off int) int32 { return int32(binary.LittleEndian.Uint32(c.Aux[off : off+4])) }
	return Config{
		Format: g(0),
		ColSeq: g(4),
		ColBeg: g(8),
		ColEnd: g(12),
		Meta:   g(16),
		Skip:   g(20),
	}, true
}

// QueryBytes returns every raw record line in dataPath whose interval
// overlaps [beg, end) on chrom, using c as the index. beg and end are
// 0-based half-open. It is the CSI counterpart of Index.QueryBytes and is
// driven entirely by the column Config carried in the CSI aux block. If
// chrom is not in the index the result is an empty slice and no error.
func (c *CSI) QueryBytes(dataPath, chrom string, beg, end int) ([][]byte, error) {
	cfg, ok := c.ConfigFromAux()
	if !ok {
		return nil, fmt.Errorf("csi: aux block lacks a column config")
	}
	refID := ChromIDInCSI(c, chrom)
	if refID < 0 {
		return nil, nil
	}
	chunks := c.RegionChunks(refID, int64(beg), int64(end))
	if len(chunks) == 0 {
		return nil, nil
	}
	f, err := os.Open(dataPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return c.queryChunks(f, cfg, chrom, beg, end, chunks)
}

// QueryBytesReader is the reader-based counterpart of QueryBytes: it returns
// every raw record line read from src (an already-open seekable bgzipped stream)
// whose interval overlaps [beg, end) on chrom. It is the entry point used for
// transparent remote (http(s)/s3/gs) access; the caller owns src and must close
// it. beg and end are 0-based half-open.
func (c *CSI) QueryBytesReader(src io.ReadSeeker, chrom string, beg, end int) ([][]byte, error) {
	cfg, ok := c.ConfigFromAux()
	if !ok {
		return nil, fmt.Errorf("csi: aux block lacks a column config")
	}
	refID := ChromIDInCSI(c, chrom)
	if refID < 0 {
		return nil, nil
	}
	chunks := c.RegionChunks(refID, int64(beg), int64(end))
	if len(chunks) == 0 {
		return nil, nil
	}
	return c.queryChunks(src, cfg, chrom, beg, end, chunks)
}

// queryChunks decodes the given CSI chunks from src and returns the matching
// record lines. It is shared by QueryBytes (local file) and QueryBytesReader
// (already-open stream).
func (c *CSI) queryChunks(src io.ReadSeeker, cfg Config, chrom string, beg, end int, chunks []CSIChunk) ([][]byte, error) {
	// Reuse the TBI record decoder via an ephemeral Index that shares the
	// column Config and name dictionary.
	idx := &Index{Config: cfg, Names: c.Names}
	tbiChunks := make([]Chunk, len(chunks))
	for i, ch := range chunks {
		tbiChunks[i] = Chunk{Beg: ch.Beg, End: ch.End}
	}
	recs, err := idx.readRegionRecords(src, chrom, beg, end, tbiChunks)
	if err != nil || len(recs) == 0 {
		return nil, err
	}
	out := make([][]byte, len(recs))
	for i := range recs {
		out[i] = recs[i].Line
	}
	return out, nil
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

	uoffToV := func(pos int64) VOffset {
		lo, hi := 0, len(offsets)
		for lo < hi {
			mid := (lo + hi) / 2
			if int64(offsets[mid].UncompressedOffset) <= pos {
				lo = mid + 1
			} else {
				hi = mid
			}
		}
		i := lo - 1
		if i < 0 {
			i = 0
		}
		blk := offsets[i]
		uoff := int(pos - blk.UncompressedOffset)
		return MakeVOffset(blk.CompressedOffset, uoff)
	}

	csi := NewCSI(minShift, 5)
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
	return csi, nil
}

func csiWrapEOF(err error) error {
	if err == io.EOF || err == io.ErrUnexpectedEOF {
		return ErrCSITruncated
	}
	return err
}
