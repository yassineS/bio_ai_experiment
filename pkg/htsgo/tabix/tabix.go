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
	"strconv"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// Magic is the 4-byte file signature at the start of every .tbi file.
var Magic = [4]byte{'T', 'B', 'I', 1}

// Format constants. The low 16 bits of Index.Format select a preset; the
// high 16 bits are coordinate-system flags.
const (
	// FormatGeneric covers GFF/BED/anything table-shaped with explicit
	// columns for chrom/begin/end.
	FormatGeneric = 0
	// FormatSAM marks a SAM-style file (the index treats unmapped reads
	// specially when this bit is set).
	FormatSAM = 1
	// FormatVCF marks a VCF file: a record begins at POS and the length is
	// derived from the REF allele.
	FormatVCF = 2

	// FlagZeroBased indicates BED-style 0-based half-open coordinates.
	FlagZeroBased = 0x10000
	// FlagEndInclusive is reserved by the spec; htslib defines it but does
	// not exercise it widely. Kept here for completeness.
	FlagEndInclusive = 0x20000
)

// Preset names accepted by --preset.
const (
	PresetGFF = "gff"
	PresetBED = "bed"
	PresetSAM = "sam"
	PresetVCF = "vcf"
)

// Config describes how to parse a data file's records. It corresponds to
// the header fields stored in a `.tbi` file.
type Config struct {
	// Format is the format/flag word stored verbatim at byte offset 8 of
	// the .tbi header. The low 16 bits are one of FormatGeneric/SAM/VCF
	// and the high 16 bits hold coordinate-system flags such as
	// FlagZeroBased.
	Format int32
	// ColSeq is the 1-based column index of the sequence (chrom) name.
	ColSeq int32
	// ColBeg is the 1-based column index of the begin / single position.
	ColBeg int32
	// ColEnd is the 1-based column index of the end position. Zero means
	// "use ColBeg" (single-position records).
	ColEnd int32
	// Meta is the first byte of comment / header lines. Defaults to '#'.
	Meta int32
	// Skip is the number of header lines to skip past before parsing
	// records.
	Skip int32
}

// PresetConfig returns the canonical Config for one of the standard preset
// names ("gff", "bed", "sam", "vcf").
func PresetConfig(name string) (Config, error) {
	switch name {
	case PresetGFF:
		return Config{Format: FormatGeneric, ColSeq: 1, ColBeg: 4, ColEnd: 5, Meta: '#', Skip: 0}, nil
	case PresetBED:
		return Config{Format: FormatGeneric | FlagZeroBased, ColSeq: 1, ColBeg: 2, ColEnd: 3, Meta: '#', Skip: 0}, nil
	case PresetSAM:
		return Config{Format: FormatSAM, ColSeq: 3, ColBeg: 4, ColEnd: 0, Meta: '@', Skip: 0}, nil
	case PresetVCF:
		return Config{Format: FormatVCF, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}, nil
	}
	return Config{}, fmt.Errorf("%w: %s", ErrBadPreset, name)
}

// ZeroBased reports whether the index uses 0-based half-open coordinates
// (bit 16 of Format).
func (c Config) ZeroBased() bool { return c.Format&FlagZeroBased != 0 }

// FormatBase returns the low 16 bits of Format (the preset selector).
func (c Config) FormatBase() int32 { return c.Format & 0xFFFF }

// Chunk is a half-open [Beg, End) range of virtual offsets.
type Chunk struct {
	Beg, End VOffset
}

// Bin is one entry in the bin-index portion of an Index.
type Bin struct {
	ID     uint32
	Chunks []Chunk
}

// RefIndex holds every bin and the linear index for one reference
// sequence.
type RefIndex struct {
	Bins   []Bin
	Linear []VOffset
}

// Index is the in-memory representation of a `.tbi` file plus the
// associated parsing Config.
type Index struct {
	Config Config
	Names  []string
	Refs   []RefIndex
	// NoCoor counts records without coordinates (optional trailing field).
	NoCoor uint64
}

// Errors returned by the tabix parser.
var (
	ErrBadMagic  = errors.New("tabix: bad .tbi magic")
	ErrTruncated = errors.New("tabix: .tbi truncated")
	ErrBadConfig = errors.New("tabix: invalid configuration")
	ErrBadRecord = errors.New("tabix: malformed record (too few columns)")
	ErrBadPreset = errors.New("tabix: unknown preset")
)

// NewIndex returns a freshly-initialised Index with the given Config.
func NewIndex(cfg Config) *Index {
	return &Index{Config: cfg}
}

// ChromID returns the index of chrom in Names, or -1 if not present.
func (idx *Index) ChromID(chrom string) int {
	for i, n := range idx.Names {
		if n == chrom {
			return i
		}
	}
	return -1
}

// Chroms returns the chromosome names recorded in the index, in their
// canonical order.
func (idx *Index) Chroms() []string {
	out := make([]string, len(idx.Names))
	copy(out, idx.Names)
	return out
}

// ensureRef grows the Refs / Names slices so that index id is addressable
// and returns it.
func (idx *Index) ensureRef(name string) int {
	if id := idx.ChromID(name); id >= 0 {
		return id
	}
	idx.Names = append(idx.Names, name)
	idx.Refs = append(idx.Refs, RefIndex{})
	return len(idx.Names) - 1
}

// recordEnd returns the half-open end coordinate (0-based, exclusive) for a
// record's parsed fields. For VCF and single-pos generic records, end is
// beg+1. For SAM, end is beg+cigarRefLen. For generic with an explicit end
// column, end is taken from that column directly (already 0-based half-open
// for BED preset; converted from 1-based-inclusive otherwise).
func (idx *Index) recordEnd(fields [][]byte) (beg, end int, err error) {
	cfg := idx.Config
	if int(cfg.ColBeg) < 1 || int(cfg.ColBeg) > len(fields) {
		return 0, 0, ErrBadRecord
	}
	bv, err := strconv.Atoi(string(fields[cfg.ColBeg-1]))
	if err != nil {
		return 0, 0, fmt.Errorf("tabix: invalid begin value: %w", err)
	}
	if !cfg.ZeroBased() {
		// 1-based, inclusive -> 0-based, half-open begin.
		bv--
	}
	beg = bv
	end = beg + 1

	if cfg.ColEnd > 0 {
		if int(cfg.ColEnd) > len(fields) {
			return 0, 0, ErrBadRecord
		}
		ev, err := strconv.Atoi(string(fields[cfg.ColEnd-1]))
		if err != nil {
			return 0, 0, fmt.Errorf("tabix: invalid end value: %w", err)
		}
		// For BED-style (0-based half-open) end is already exclusive; for
		// GFF (1-based inclusive) end is inclusive so the half-open form
		// is unchanged numerically (end_inclusive == end_exclusive_0based
		// because we subtracted from beg only).
		end = ev
		if cfg.ZeroBased() {
			// BED preset: end is exclusive already.
		}
	} else if cfg.FormatBase() == FormatVCF {
		// VCF: REF allele is column 4 (1-based). End is beg + len(REF).
		if len(fields) < 4 {
			return 0, 0, ErrBadRecord
		}
		end = beg + len(fields[3])
		if end <= beg {
			end = beg + 1
		}
	}
	if end <= beg {
		end = beg + 1
	}
	if beg < 0 {
		beg = 0
	}
	return beg, end, nil
}

// splitFields slices line by tab characters into pre-allocated `out` if
// possible. It is allocation-light: returns sub-slices of the input line.
func splitFields(line []byte, out [][]byte) [][]byte {
	out = out[:0]
	start := 0
	for i := 0; i < len(line); i++ {
		if line[i] == '\t' {
			out = append(out, line[start:i])
			start = i + 1
		}
	}
	out = append(out, line[start:])
	return out
}

// Build streams the bgzipped data file at `path` and constructs an in-memory
// Index using cfg. The data file must be sorted by (chrom, beg). Build
// returns the populated Index ready to be Written or Queried.
func Build(path string, cfg Config) (*Index, error) {
	if err := validateConfig(cfg); err != nil {
		return nil, err
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	idx := NewIndex(cfg)

	// We need both the raw decoded text *and* the virtual offset at which
	// each line begins. The bgzip Reader does not expose offsets, so we
	// scan the BGZF blocks once and decode them sequentially, accumulating
	// each block's offset.
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

	// uncompressedToVOffset maps an absolute uncompressed byte position to
	// the virtual offset of that byte. We need this every time a line
	// starts.
	uoffToV := func(pos int64) VOffset {
		// Binary search offsets for the block whose
		// [UncompressedOffset, UncompressedOffset+UncompressedSize)
		// contains pos.
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

	skip := int(cfg.Skip)
	scanner := bufio.NewScanner(bytes.NewReader(data))
	scanner.Buffer(make([]byte, 0, 1<<16), 1<<24)
	var pos int64
	var lineNo int
	var fieldBuf [64][]byte
	for scanner.Scan() {
		raw := scanner.Bytes()
		lineStart := pos
		pos += int64(len(raw)) + 1 // +1 for the consumed '\n'
		lineNo++
		if len(raw) == 0 {
			continue
		}
		// Skip header / comment lines.
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
			idx.NoCoor++
			continue
		}
		beg, end, err := idx.recordEnd(fields)
		if err != nil {
			return nil, fmt.Errorf("tabix: line %d: %w", lineNo, err)
		}
		refID := idx.ensureRef(chrom)
		v := uoffToV(lineStart)
		vEnd := uoffToV(pos) // virtual offset just past this line
		idx.addRecord(refID, beg, end, v, vEnd)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	idx.finalize()
	return idx, nil
}

// addRecord registers one record with bin and linear-index slots.
func (idx *Index) addRecord(refID, beg, end int, v, vEnd VOffset) {
	ref := &idx.Refs[refID]
	binID := uint32(Reg2bin(beg, end))

	// Find or create the bin entry.
	var bin *Bin
	for i := range ref.Bins {
		if ref.Bins[i].ID == binID {
			bin = &ref.Bins[i]
			break
		}
	}
	if bin == nil {
		ref.Bins = append(ref.Bins, Bin{ID: binID})
		bin = &ref.Bins[len(ref.Bins)-1]
	}
	// Coalesce with previous chunk if they touch / overlap (sorted input).
	if n := len(bin.Chunks); n > 0 && bin.Chunks[n-1].End >= v {
		if vEnd > bin.Chunks[n-1].End {
			bin.Chunks[n-1].End = vEnd
		}
	} else {
		bin.Chunks = append(bin.Chunks, Chunk{Beg: v, End: vEnd})
	}

	// Linear index: every 16-kbp tile spanned by the record gets the
	// minimum virtual offset of any record that touches it. We use
	// ^VOffset(0) as the "no record yet" sentinel because VOffset(0) is
	// a perfectly legal offset (the first byte of the first block).
	first := LinearTile(beg)
	last := LinearTile(end - 1)
	if last < first {
		last = first
	}
	for last >= len(ref.Linear) {
		ref.Linear = append(ref.Linear, ^VOffset(0))
	}
	for t := first; t <= last; t++ {
		if ref.Linear[t] == ^VOffset(0) || v < ref.Linear[t] {
			ref.Linear[t] = v
		}
	}
}

// finalize fills any "no record" sentinel slots in the linear index using
// the htslib convention "carry forward the previous value". This makes
// the index lookup work for queries that start in a tile with no records.
// Sentinel-prefixed tiles (those before the first record) become 0.
func (idx *Index) finalize() {
	for r := range idx.Refs {
		lin := idx.Refs[r].Linear
		var last VOffset
		seen := false
		for i := range lin {
			if lin[i] == ^VOffset(0) {
				if seen {
					lin[i] = last
				} else {
					lin[i] = 0
				}
			} else {
				last = lin[i]
				seen = true
			}
		}
	}
}

// validateConfig checks the Config for the minimum invariants required to
// build a meaningful index.
func validateConfig(cfg Config) error {
	if cfg.ColSeq < 1 || cfg.ColBeg < 1 {
		return ErrBadConfig
	}
	if cfg.ColEnd < 0 {
		return ErrBadConfig
	}
	return nil
}

// Write serialises idx to w in the canonical `.tbi` byte layout.
// The output is little-endian and matches htslib byte-for-byte.
func (idx *Index) Write(w io.Writer) error {
	bw := bufio.NewWriter(w)
	if _, err := bw.Write(Magic[:]); err != nil {
		return err
	}
	hdr := []int32{
		int32(len(idx.Names)),
		idx.Config.Format,
		idx.Config.ColSeq,
		idx.Config.ColBeg,
		idx.Config.ColEnd,
		idx.Config.Meta,
		idx.Config.Skip,
	}
	for _, v := range hdr {
		if err := binary.Write(bw, binary.LittleEndian, v); err != nil {
			return err
		}
	}
	// Concatenated NUL-terminated names.
	var nameBuf bytes.Buffer
	for _, n := range idx.Names {
		nameBuf.WriteString(n)
		nameBuf.WriteByte(0)
	}
	if err := binary.Write(bw, binary.LittleEndian, int32(nameBuf.Len())); err != nil {
		return err
	}
	if _, err := bw.Write(nameBuf.Bytes()); err != nil {
		return err
	}
	for i := range idx.Refs {
		ref := idx.Refs[i]
		// Sort bins by ID for stable output.
		sort.Slice(ref.Bins, func(a, b int) bool { return ref.Bins[a].ID < ref.Bins[b].ID })
		if err := binary.Write(bw, binary.LittleEndian, int32(len(ref.Bins))); err != nil {
			return err
		}
		for _, bin := range ref.Bins {
			if err := binary.Write(bw, binary.LittleEndian, bin.ID); err != nil {
				return err
			}
			if err := binary.Write(bw, binary.LittleEndian, int32(len(bin.Chunks))); err != nil {
				return err
			}
			for _, c := range bin.Chunks {
				if err := binary.Write(bw, binary.LittleEndian, uint64(c.Beg)); err != nil {
					return err
				}
				if err := binary.Write(bw, binary.LittleEndian, uint64(c.End)); err != nil {
					return err
				}
			}
		}
		if err := binary.Write(bw, binary.LittleEndian, int32(len(ref.Linear))); err != nil {
			return err
		}
		for _, v := range ref.Linear {
			if err := binary.Write(bw, binary.LittleEndian, uint64(v)); err != nil {
				return err
			}
		}
	}
	if idx.NoCoor > 0 {
		if err := binary.Write(bw, binary.LittleEndian, idx.NoCoor); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// WriteFile encodes idx to a .tbi-formatted bgzipped file at path. It is
// what tabix uses to emit `<data>.tbi` on disk.
func (idx *Index) WriteFile(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bgzip.NewWriter(f)
	if err := idx.Write(bw); err != nil {
		return err
	}
	return bw.Close()
}

// Read parses a .tbi byte stream (already BGZF-decoded) into a fresh Index.
func Read(r io.Reader) (*Index, error) {
	br := bufio.NewReader(r)
	var magic [4]byte
	if _, err := io.ReadFull(br, magic[:]); err != nil {
		if err == io.EOF || err == io.ErrUnexpectedEOF {
			return nil, ErrTruncated
		}
		return nil, err
	}
	if magic != Magic {
		return nil, ErrBadMagic
	}
	var hdr [7]int32
	for i := range hdr {
		if err := binary.Read(br, binary.LittleEndian, &hdr[i]); err != nil {
			return nil, wrapEOF(err)
		}
	}
	nRef := hdr[0]
	if nRef < 0 {
		return nil, ErrBadMagic
	}
	cfg := Config{
		Format: hdr[1],
		ColSeq: hdr[2],
		ColBeg: hdr[3],
		ColEnd: hdr[4],
		Meta:   hdr[5],
		Skip:   hdr[6],
	}
	idx := NewIndex(cfg)

	var lNm int32
	if err := binary.Read(br, binary.LittleEndian, &lNm); err != nil {
		return nil, wrapEOF(err)
	}
	if lNm < 0 {
		return nil, ErrTruncated
	}
	nameBuf := make([]byte, lNm)
	if _, err := io.ReadFull(br, nameBuf); err != nil {
		return nil, wrapEOF(err)
	}
	names := splitNullStrings(nameBuf)
	if int32(len(names)) != nRef {
		// Some encoders emit a trailing empty entry after the final NUL;
		// be lenient when that happens.
		if int32(len(names)) == nRef+1 && names[len(names)-1] == "" {
			names = names[:nRef]
		} else {
			return nil, fmt.Errorf("tabix: name count %d does not match n_ref %d", len(names), nRef)
		}
	}
	idx.Names = names
	idx.Refs = make([]RefIndex, nRef)

	for r := int32(0); r < nRef; r++ {
		var nBin int32
		if err := binary.Read(br, binary.LittleEndian, &nBin); err != nil {
			return nil, wrapEOF(err)
		}
		bins := make([]Bin, nBin)
		for i := range bins {
			if err := binary.Read(br, binary.LittleEndian, &bins[i].ID); err != nil {
				return nil, wrapEOF(err)
			}
			var nChunk int32
			if err := binary.Read(br, binary.LittleEndian, &nChunk); err != nil {
				return nil, wrapEOF(err)
			}
			chunks := make([]Chunk, nChunk)
			for j := range chunks {
				var beg, end uint64
				if err := binary.Read(br, binary.LittleEndian, &beg); err != nil {
					return nil, wrapEOF(err)
				}
				if err := binary.Read(br, binary.LittleEndian, &end); err != nil {
					return nil, wrapEOF(err)
				}
				chunks[j] = Chunk{Beg: VOffset(beg), End: VOffset(end)}
			}
			bins[i].Chunks = chunks
		}
		var nIntv int32
		if err := binary.Read(br, binary.LittleEndian, &nIntv); err != nil {
			return nil, wrapEOF(err)
		}
		linear := make([]VOffset, nIntv)
		for i := range linear {
			var v uint64
			if err := binary.Read(br, binary.LittleEndian, &v); err != nil {
				return nil, wrapEOF(err)
			}
			linear[i] = VOffset(v)
		}
		idx.Refs[r] = RefIndex{Bins: bins, Linear: linear}
	}
	// Optional trailing n_no_coor.
	var trailer uint64
	if err := binary.Read(br, binary.LittleEndian, &trailer); err == nil {
		idx.NoCoor = trailer
	}
	return idx, nil
}

func wrapEOF(err error) error {
	if err == io.EOF || err == io.ErrUnexpectedEOF {
		return ErrTruncated
	}
	return err
}

// ReadFile opens path (which must be a BGZF-compressed .tbi) and parses it.
func ReadFile(path string) (*Index, error) {
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
	return Read(br)
}

// splitNullStrings splits a NUL-delimited byte run into strings, dropping a
// trailing empty entry that arises from a terminator at end of buffer.
func splitNullStrings(buf []byte) []string {
	var out []string
	start := 0
	for i := 0; i < len(buf); i++ {
		if buf[i] == 0 {
			out = append(out, string(buf[start:i]))
			start = i + 1
		}
	}
	if start < len(buf) {
		out = append(out, string(buf[start:]))
	}
	return out
}

// QueryBytes returns every raw line whose interval overlaps [beg, end) on
// chrom. beg and end are 0-based half-open. The returned byte slices each
// represent one full record (no trailing newline). If chrom is not in the
// index the result is an empty slice and no error.
func (idx *Index) QueryBytes(dataPath, chrom string, beg, end int) ([][]byte, error) {
	chunks, refID, err := idx.RegionChunks(chrom, beg, end)
	if err != nil {
		return nil, err
	}
	if len(chunks) == 0 || refID < 0 {
		return nil, nil
	}
	return idx.readRegionRecords(dataPath, chrom, beg, end, chunks)
}

// RegionChunks returns the set of merged chunks that may contain records
// overlapping [beg, end) on chrom. It also returns the chrom's reference
// ID for callers that want to short-circuit further work when the chrom is
// not indexed (refID will be -1 in that case).
func (idx *Index) RegionChunks(chrom string, beg, end int) ([]Chunk, int, error) {
	if beg < 0 {
		beg = 0
	}
	if end <= beg {
		return nil, -1, nil
	}
	refID := idx.ChromID(chrom)
	if refID < 0 {
		return nil, -1, nil
	}
	ref := idx.Refs[refID]

	// Linear-index trim: the first chunk that can contain a record
	// overlapping [beg, end) is the one whose Beg ≥ Linear[beg>>14].
	minOff := VOffset(0)
	tile := LinearTile(beg)
	if tile < len(ref.Linear) {
		minOff = ref.Linear[tile]
	}

	bins := Reg2bins(beg, end)
	binSet := make(map[uint32]struct{}, len(bins))
	for _, b := range bins {
		binSet[uint32(b)] = struct{}{}
	}

	var chunks []Chunk
	for _, bin := range ref.Bins {
		if _, ok := binSet[bin.ID]; !ok {
			continue
		}
		for _, c := range bin.Chunks {
			if c.End <= minOff {
				continue
			}
			if c.Beg < minOff {
				c.Beg = minOff
			}
			chunks = append(chunks, c)
		}
	}
	if len(chunks) == 0 {
		return nil, refID, nil
	}
	sort.Slice(chunks, func(a, b int) bool { return chunks[a].Beg < chunks[b].Beg })
	// Merge overlapping / adjacent chunks.
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
	return merged, refID, nil
}

// readRegionRecords decodes the requested chunks from dataPath and returns
// every record line whose interval overlaps [beg, end) on chrom.
func (idx *Index) readRegionRecords(dataPath, chrom string, beg, end int, chunks []Chunk) ([][]byte, error) {
	f, err := os.Open(dataPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var out [][]byte
	var fieldBuf [64][]byte
	for _, c := range chunks {
		if c.Beg >= c.End {
			continue
		}
		data, err := readVirtualRange(f, c.Beg, c.End)
		if err != nil {
			return nil, err
		}
		// Each block boundary may split a line; the BGZF reader already
		// handles that because we read a contiguous decompressed slice.
		// We then split on newlines and filter records.
		start := 0
		for i := 0; i <= len(data); i++ {
			if i == len(data) || data[i] == '\n' {
				line := data[start:i]
				start = i + 1
				if len(line) == 0 {
					continue
				}
				if line[0] == byte(idx.Config.Meta) {
					continue
				}
				fields := splitFields(line, fieldBuf[:0])
				if int(idx.Config.ColSeq) > len(fields) {
					continue
				}
				if string(fields[idx.Config.ColSeq-1]) != chrom {
					continue
				}
				rbeg, rend, err := idx.recordEnd(fields)
				if err != nil {
					continue
				}
				if rend <= beg || rbeg >= end {
					continue
				}
				// Copy because the underlying buffer is reused.
				cp := make([]byte, len(line))
				copy(cp, line)
				out = append(out, cp)
			}
		}
	}
	return out, nil
}

// readVirtualRange seeks to the compressed block that contains begV and
// reads sequentially until the byte at endV (exclusive) has been decoded,
// returning the decompressed bytes from begV.Uoff() through endV.
func readVirtualRange(f *os.File, begV, endV VOffset) ([]byte, error) {
	startBlock := begV.Coff()
	startInBlock := begV.Uoff()
	endBlock := endV.Coff()
	endInBlock := endV.Uoff()

	if _, err := f.Seek(startBlock, io.SeekStart); err != nil {
		return nil, err
	}
	br, err := bgzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer br.Close()

	// Read enough to span past endBlock. We read until we have decoded
	// (endBlock - startBlock) compressed bytes' worth plus the in-block
	// offset endInBlock. The simplest correct approach: read until the
	// decompressed length is at least (endBlock - startBlock)*expansion
	// + endInBlock. But we don't know the expansion factor up front.
	//
	// Instead: read in 64 KiB decoded chunks until the total compressed
	// position is past endBlock+endInBlock, or until EOF. Since we cannot
	// inspect the underlying file position from inside bgzip.Reader (it
	// owns the consumption), we use the simpler strategy of reading the
	// raw byte span needed: we know endV.Coff() is the start of the
	// block AFTER the last record, plus an in-block offset within it.
	// Read until we have read (endBlock - startBlock) blocks worth of
	// data plus endInBlock more bytes.
	//
	// We approximate by reading until io.EOF or until our cumulative
	// decoded-byte count crosses the boundary. The exact boundary is
	// tracked by walking the bgzip Scan output we already keep.
	want := int64(0)
	if endBlock > startBlock {
		// We need every byte from offset startInBlock of the first
		// block through endInBlock of the end block. The decompressed
		// length is at most the sum of all block ISIZEs in between
		// plus endInBlock minus startInBlock. We do not have ISIZEs
		// handy here, so read everything until we have read enough
		// data to span the compressed range — bgzip blocks are at
		// most MaxBlockSize=65280 bytes uncompressed each. As a safe
		// upper bound, read all available data and trim afterwards.
		_ = want
	}

	// Read up to (endBlock - startBlock + 1) * MaxBlockSize bytes.
	// That covers every block from start through end inclusive.
	upper := int(endBlock-startBlock+1) * bgzip.MaxBlockSize
	if upper <= 0 || upper > 1<<30 {
		upper = 1 << 30
	}
	buf := make([]byte, 0, 1<<16)
	tmp := make([]byte, 1<<16)
	total := 0
	for total < upper {
		n, err := br.Read(tmp)
		if n > 0 {
			buf = append(buf, tmp[:n]...)
			total += n
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			// ErrTruncated is benign here because we may be reading
			// past the indexed region into the trailing portion of
			// the file.
			if errors.Is(err, bgzip.ErrTruncated) {
				break
			}
			return nil, err
		}
	}

	// Trim to [startInBlock, decompressed-end-position]. The end position
	// in the buffer is the sum of (uncompressed sizes of all blocks
	// strictly between startBlock and endBlock) + endInBlock. To compute
	// that we re-scan the BGZF blocks within the slice... but we don't
	// have offsets handy here. Instead, fall back to using the entire
	// decoded buffer past startInBlock; the caller already filters by
	// chrom and coordinates, so trailing extra bytes are harmless.
	_ = endInBlock
	if startInBlock > len(buf) {
		return nil, nil
	}
	return buf[startInBlock:], nil
}
