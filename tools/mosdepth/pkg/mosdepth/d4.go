package mosdepth

import (
	"encoding/binary"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// D4 (Dense Depth Data Dump) is the binary depth-track container format used
// by mosdepth's `-d/--d4` mode. This file implements a pure-Go writer that
// reproduces, byte-for-byte, the on-disk format the upstream mosdepth binary
// emits when built with D4 support (which links the Rust `d4` crate via the
// `d4binding` C API). It also provides a matching reader used by the
// round-trip parity tests.
//
// Format overview
//
// A D4 file is a "d4-framefile": a small hierarchical container holding named
// streams, blobs and sub-directories. The container layout is:
//
//	[4-byte magic "d4\xdd\xdd"][4 zero bytes]
//	[root directory stream]            (the table of contents)
//	[".metadata" stream]               (JSON header: chroms, dictionary, denominator)
//	[".ptab" blob]                     (bit-packed dense primary depth table)
//	[".stab" sub-directory]            (secondary table: per-chrom overflow streams)
//	[".index" sub-directory]           (secondary frame index blob)
//
// Streams are stored as singly-linked lists of fixed-size frames. Each frame
// begins with a 16-byte FrameHeader: an int64 little-endian relative offset to
// the next frame (0 = last frame) followed by a uint64 little-endian size of
// that next frame. The remaining bytes of the frame are payload.
//
// The primary table is the dense, bit-width-packed depth track. mosdepth uses
// the fixed dictionary SimpleRange{low:0, high:128}, i.e. a 7-bit code per
// base (128 distinct codes). Depth d in [0,128) is stored as code d. Depth
// >= 128 is clamped to the all-ones code 127 in the primary table, matching
// upstream mosdepth's per-base D4 output, which caps over-dictionary depths and
// does NOT populate the secondary table (the chromosome summary still reports
// the true maximum depth). Codes are packed least-significant-bit first into a
// contiguous byte array, one chromosome after another in chrom-list order;
// chromosome c occupies (size*7 + 7)/8 bytes.
//
// The secondary table (".stab") is a SimpleKV/range sparse array: a
// ".metadata" stream describing one partition per chromosome, plus one data
// stream per chromosome (named by its index "0", "1", ...) that would hold
// out-of-dictionary overflow records. For mosdepth's per-base output every
// overflow stream is empty (a single zero-filled frame). The ".index"
// sub-directory holds a "secondary_frame_index" blob. The reader below still
// understands populated overflow streams so it can decode general D4 files.
//
// This writer reproduces the exact frame sizes, padding and directory layout
// the upstream encoder produces, so that for the same BAM the resulting
// `.per-base.d4` is byte-identical to upstream's.

// d4Magic is the 4-byte D4 file signature, matching upstream's FILE_MAGIC_NUM.
var d4Magic = [4]byte{'d', '4', 0xdd, 0xdd}

// d4DictHigh is the exclusive upper bound of the SimpleRange dictionary
// mosdepth uses (low=0, high=128). It yields a 7-bit primary-table code width.
const d4DictHigh = 128

// d4Bits is the per-value bit width of the dense primary encoding, derived
// from the dictionary size (floor(log2(128)) = 7).
const d4Bits = 7

// d4SentinelCode is the all-ones primary code (2^d4Bits - 1 = 127) written
// when a depth is too large for the dictionary; the real value then lives in
// the secondary table.
const d4SentinelCode = (1 << d4Bits) - 1

// d4InitBlockSize is the directory/stream init frame size used by the upstream
// framefile (Directory::INIT_BLOCK_SIZE and the metadata/stab stream frame
// size).
const d4InitBlockSize = 512

// d4FrameHeaderSize is the size of a stream FrameHeader (int64 link + uint64
// size).
const d4FrameHeaderSize = 16

// d4Chrom is one chromosome entry: a name and a length in bases.
type d4Chrom struct {
	Name   string
	Length int64
}

// d4Writer streams a dense per-base depth track to a D4 file in the upstream
// on-disk format. Chromosomes are written one at a time via writeChrom in the
// declared order; the bit-packed primary bytes are appended to the ".ptab"
// blob region as they are produced, and the directory, metadata, secondary
// table and index are finalised in Close.
type d4Writer struct {
	f      *os.File
	path   string
	chroms []d4Chrom

	// order/next track the expected chromosome write sequence.
	next int

	// ptabStart is the absolute file offset of the ".ptab" blob payload.
	ptabStart int64
	// pos is the current write offset within the ".ptab" blob.
	pos int64
}

// d4Overflow is a single secondary-table record: a half-open base range
// [Left,Right) on one chromosome that all share an out-of-dictionary Value.
type d4Overflow struct {
	Left  uint32
	Right uint32
	Value int32
}

// newD4Writer creates a D4 file at path whose metadata declares the given
// chromosomes in order. The magic, root directory, metadata stream and the
// (zero-filled) primary-table blob region are written up front; callers then
// invoke writeChrom once per chromosome in order before calling Close.
func newD4Writer(path string, chroms []d4Chrom) (*d4Writer, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	w := &d4Writer{
		f:      f,
		path:   path,
		chroms: chroms,
	}
	if err := w.writeHeaderAndReservePtab(); err != nil {
		_ = f.Close()
		return nil, err
	}
	return w, nil
}

// ptabBytesForLen returns the number of primary-table bytes a chromosome of
// the given base length occupies: ceil(length * bits / 8).
func ptabBytesForLen(length int64) int64 {
	return (length*int64(d4Bits) + 7) / 8
}

// primaryTableSize returns the total size of the ".ptab" blob: the sum of the
// per-chromosome packed byte counts.
func (w *d4Writer) primaryTableSize() int64 {
	var total int64
	for _, c := range w.chroms {
		total += ptabBytesForLen(c.Length)
	}
	return total
}

// metadataJSON builds the exact ".metadata" JSON header string the upstream
// encoder writes: chrom_list, then the SimpleRange dictionary, then the
// denominator ("One").
func (w *d4Writer) metadataJSON() []byte {
	var b strings.Builder
	b.WriteString(`{"chrom_list":[`)
	for i, c := range w.chroms {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(`{"name":`)
		b.WriteString(strconv.Quote(c.Name))
		b.WriteString(`,"size":`)
		b.WriteString(strconv.FormatInt(c.Length, 10))
		b.WriteByte('}')
	}
	b.WriteString(`],"dictionary":{"SimpleRange":{"low":0,"high":`)
	b.WriteString(strconv.Itoa(d4DictHigh))
	b.WriteString(`}},"denominator":"One"}`)
	return []byte(b.String())
}

// writeStream writes a framefile stream whose payload is data, using the given
// frame size, at the current end of file. It returns the absolute offset of
// the stream's primary (first) frame and the size of that primary frame.
//
// The stream is a chain of contiguous frames, each frameSize bytes: a 16-byte
// header (link to next frame, or 0 for the last) followed by up to
// frameSize-16 payload bytes. The final frame is zero-padded to frameSize.
func writeStream(f *os.File, data []byte, frameSize int) (primaryOffset, primarySize int64, err error) {
	start, err := f.Seek(0, 2)
	if err != nil {
		return 0, 0, err
	}
	payloadCap := frameSize - d4FrameHeaderSize
	// At least one frame is always emitted (even for empty payload), matching
	// the upstream Stream::create behaviour.
	nFrames := 1
	if len(data) > 0 {
		nFrames = (len(data) + payloadCap - 1) / payloadCap
	}
	buf := make([]byte, frameSize)
	off := 0
	for i := 0; i < nFrames; i++ {
		for j := range buf {
			buf[j] = 0
		}
		last := i == nFrames-1
		if !last {
			// link = +frameSize (next frame immediately follows), size = frameSize
			binary.LittleEndian.PutUint64(buf[0:8], uint64(frameSize))
			binary.LittleEndian.PutUint64(buf[8:16], uint64(frameSize))
		}
		n := copy(buf[d4FrameHeaderSize:], data[off:])
		off += n
		if _, werr := f.Write(buf); werr != nil {
			return 0, 0, werr
		}
	}
	return start, int64(frameSize), nil
}

// reserveBlob appends a zero-filled region of size bytes at the current end of
// file and returns its starting offset. Used for the ".ptab" blob and the
// secondary frame index blob.
func reserveBlob(f *os.File, size int64) (int64, error) {
	start, err := f.Seek(0, 2)
	if err != nil {
		return 0, err
	}
	if size == 0 {
		return start, nil
	}
	if err := f.Truncate(start + size); err != nil {
		return 0, err
	}
	if _, err := f.Seek(0, 2); err != nil {
		return 0, err
	}
	return start, nil
}

// dirEntry is one directory table entry.
type dirEntry struct {
	kind          byte // 0=stream, 1=subdir, 2=blob
	primaryOffset int64
	primarySize   int64
	name          string
}

// encodeDirEntries serialises directory entries into the byte payload of a
// directory stream, exactly as upstream's append_directory does. For each
// entry it writes a 0x01 "has next" flag, the kind byte, the int64
// little-endian offset *relative to the directory base*, the uint64
// little-endian size, and the NUL-terminated name. Upstream additionally
// writes a trailing 0x00 "no more entries" sentinel after each entry, but the
// next entry's flag rewinds one byte and overwrites that sentinel — so on disk
// only a single 0x00 separates the name from the next entry's flag, and one
// final 0x00 sentinel follows the last entry.
func encodeDirEntries(base int64, entries []dirEntry) []byte {
	var out []byte
	var tmp [8]byte
	for _, e := range entries {
		out = append(out, 1, e.kind)
		binary.LittleEndian.PutUint64(tmp[:], uint64(e.primaryOffset-base))
		out = append(out, tmp[:]...)
		binary.LittleEndian.PutUint64(tmp[:], uint64(e.primarySize))
		out = append(out, tmp[:]...)
		out = append(out, []byte(e.name)...)
		out = append(out, 0)
	}
	// Final "no more entries" sentinel after the last entry.
	out = append(out, 0)
	return out
}

// writeHeaderAndReservePtab writes the magic, root directory frame, metadata
// stream and reserves the primary-table blob region. The root directory frame
// is written after the metadata stream and ptab blob are laid out, because its
// entries reference their absolute offsets; to keep a single 512-byte root
// frame (as upstream does) we write the magic, leave a 512-byte gap for the
// root frame, then write metadata + ptab, and finally patch the root frame in
// Close once stab/index offsets are known. We compute all offsets up front
// here since they are fully determined by the chrom list.
func (w *d4Writer) writeHeaderAndReservePtab() error {
	// Magic (4) + 4 zero bytes.
	if _, err := w.f.Write(d4Magic[:]); err != nil {
		return err
	}
	if _, err := w.f.Write([]byte{0, 0, 0, 0}); err != nil {
		return err
	}
	// Reserve the root directory frame (single 512-byte frame). It is patched
	// in Close once all child offsets are known. Upstream's root directory
	// fits all four entries in one 512-byte frame for any realistic chrom
	// count.
	if _, err := reserveBlob(w.f, d4InitBlockSize); err != nil {
		return err
	}
	// ".metadata" stream.
	if _, _, err := writeStream(w.f, w.metadataJSON(), d4InitBlockSize); err != nil {
		return err
	}
	// ".ptab" blob, zero-filled; we will fill it in chromosome by chromosome.
	start, err := reserveBlob(w.f, w.primaryTableSize())
	if err != nil {
		return err
	}
	w.ptabStart = start
	w.pos = start
	return nil
}

// writeChrom packs the dense depth codes for one chromosome and appends them
// to the primary-table blob. depths must contain exactly the chromosome's
// declared length of values, and writeChrom must be called once per chromosome
// in chrom-list order.
//
// Depths >= the dictionary high bound (128) are clamped to the all-ones
// sentinel code (127) in the primary table, and NO secondary overflow record
// is written — this matches the upstream mosdepth per-base D4 output exactly,
// where the d4 C-binding writer caps per-base values at the dictionary maximum
// and leaves the secondary table empty. (The chromosome-wide summary still
// reports the true maximum depth; only the D4 per-base track is capped.)
func (w *d4Writer) writeChrom(name string, depths []int32) error {
	if w.next >= len(w.chroms) {
		return fmt.Errorf("mosdepth: D4 writeChrom for %q after all chromosomes written", name)
	}
	want := w.chroms[w.next]
	if want.Name != name {
		return fmt.Errorf("mosdepth: D4 writeChrom out of order: got %q, expected %q", name, want.Name)
	}
	if int64(len(depths)) != want.Length {
		return fmt.Errorf("mosdepth: D4 writeChrom %q: got %d depths, expected %d", name, len(depths), want.Length)
	}

	nbytes := ptabBytesForLen(want.Length)
	packed := make([]byte, nbytes)
	for i, d := range depths {
		var code uint32
		switch {
		case d < 0:
			code = 0
		case int(d) >= d4DictHigh:
			code = d4SentinelCode
		default:
			code = uint32(d)
		}
		bit := int64(i) * int64(d4Bits)
		byteOff := bit / 8
		shift := uint(bit % 8)
		// LSB-first packing into the little-endian byte array. A 7-bit code
		// spans at most two bytes.
		packed[byteOff] |= byte(code << shift)
		if shift+d4Bits > 8 {
			packed[byteOff+1] |= byte(code >> (8 - shift))
		}
	}

	if _, err := w.f.WriteAt(packed, w.pos); err != nil {
		return err
	}
	w.pos += nbytes
	w.next++
	return nil
}

// Close finalises the secondary table (".stab"), the index (".index") and the
// root directory, then closes the file. It is an error to Close before every
// declared chromosome has been written.
func (w *d4Writer) Close() error {
	if w.next != len(w.chroms) {
		err := fmt.Errorf("mosdepth: D4 file closed with %d of %d chromosomes written", w.next, len(w.chroms))
		_ = w.f.Close()
		return err
	}

	// ".stab" sub-directory. Build its directory frame, a ".metadata" stream
	// and one data stream per chromosome, all laid out contiguously after the
	// ptab blob.
	stabOffset, stabSize, err := w.writeStabDirectory()
	if err != nil {
		_ = w.f.Close()
		return err
	}

	// ".index" sub-directory: a directory frame plus a "secondary_frame_index"
	// blob.
	indexOffset, indexSize, err := w.writeIndexDirectory()
	if err != nil {
		_ = w.f.Close()
		return err
	}

	// Patch the root directory frame (single 512-byte frame at offset 8) with
	// the four entries now that all offsets are known.
	metaOffset := int64(8 + d4InitBlockSize)
	metaSize := int64(d4InitBlockSize)
	ptabSize := w.primaryTableSize()
	entries := []dirEntry{
		{kind: 0, primaryOffset: metaOffset, primarySize: metaSize, name: ".metadata"},
		{kind: 2, primaryOffset: w.ptabStart, primarySize: ptabSize, name: ".ptab"},
		{kind: 1, primaryOffset: stabOffset, primarySize: stabSize, name: ".stab"},
		{kind: 1, primaryOffset: indexOffset, primarySize: indexSize, name: ".index"},
	}
	rootBase := int64(8)
	payload := encodeDirEntries(rootBase, entries)
	if len(payload) > d4InitBlockSize-d4FrameHeaderSize {
		_ = w.f.Close()
		return fmt.Errorf("mosdepth: D4 root directory too large for one frame (%d bytes)", len(payload))
	}
	frame := make([]byte, d4InitBlockSize)
	// link=0, size=0 (single frame); header already zero.
	copy(frame[d4FrameHeaderSize:], payload)
	if _, err := w.f.WriteAt(frame, rootBase); err != nil {
		_ = w.f.Close()
		return err
	}
	return w.f.Close()
}

// writeStabDirectory writes the ".stab" sub-directory and returns its offset
// and total size. The sub-directory holds a ".metadata" stream (the SimpleKV
// partition description) and one data stream per chromosome.
func (w *d4Writer) writeStabDirectory() (offset, size int64, err error) {
	dirBase, err := w.f.Seek(0, 2)
	if err != nil {
		return 0, 0, err
	}
	// Reserve the sub-directory's own frame (single 512-byte frame). Its
	// entries reference offsets relative to dirBase; patch after children.
	if _, err := reserveBlob(w.f, d4InitBlockSize); err != nil {
		return 0, 0, err
	}

	var entries []dirEntry

	// ".metadata" stream for the sparse array.
	stabMeta := w.stabMetadataJSON()
	mOff, mSize, err := writeStream(w.f, stabMeta, d4InitBlockSize)
	if err != nil {
		return 0, 0, err
	}
	entries = append(entries, dirEntry{kind: 0, primaryOffset: mOff, primarySize: mSize, name: ".metadata"})

	// One data stream per chromosome (named by index). mosdepth's per-base D4
	// output clamps over-dictionary depths into the primary table and never
	// populates the secondary table, so every overflow stream is empty (a
	// single zero-filled frame), matching upstream byte-for-byte.
	for i := range w.chroms {
		sOff, sSize, err := writeStream(w.f, nil, d4InitBlockSize)
		if err != nil {
			return 0, 0, err
		}
		entries = append(entries, dirEntry{kind: 0, primaryOffset: sOff, primarySize: sSize, name: strconv.Itoa(i)})
	}

	// Patch the sub-directory frame.
	payload := encodeDirEntries(dirBase, entries)
	if len(payload) > d4InitBlockSize-d4FrameHeaderSize {
		return 0, 0, fmt.Errorf("mosdepth: D4 stab directory too large for one frame (%d bytes)", len(payload))
	}
	frame := make([]byte, d4InitBlockSize)
	copy(frame[d4FrameHeaderSize:], payload)
	if _, err := w.f.WriteAt(frame, dirBase); err != nil {
		return 0, 0, err
	}

	end, err := w.f.Seek(0, 2)
	if err != nil {
		return 0, 0, err
	}
	return dirBase, end - dirBase, nil
}

// stabMetadataJSON builds the ".stab/.metadata" JSON: SimpleKV / range format
// with one partition per chromosome spanning its full length, no compression.
func (w *d4Writer) stabMetadataJSON() []byte {
	var b strings.Builder
	b.WriteString(`{"format":"SimpleKV","record_format":"range","partitions":[`)
	for i, c := range w.chroms {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteByte('[')
		b.WriteString(strconv.Quote(c.Name))
		b.WriteString(`,0,`)
		b.WriteString(strconv.FormatInt(c.Length, 10))
		b.WriteByte(']')
	}
	b.WriteString(`],"compression":"NoCompression"}`)
	return []byte(b.String())
}

// writeIndexDirectory writes the ".index" sub-directory (a directory frame
// plus the "secondary_frame_index" blob) and returns its offset and size.
func (w *d4Writer) writeIndexDirectory() (offset, size int64, err error) {
	dirBase, err := w.f.Seek(0, 2)
	if err != nil {
		return 0, 0, err
	}
	if _, err := reserveBlob(w.f, d4InitBlockSize); err != nil {
		return 0, 0, err
	}
	// "secondary_frame_index" blob: 8 bytes (a single uint64 count, 0 when no
	// secondary frames need indexing).
	blobOff, err := reserveBlob(w.f, 8)
	if err != nil {
		return 0, 0, err
	}
	entries := []dirEntry{
		{kind: 2, primaryOffset: blobOff, primarySize: 8, name: "secondary_frame_index"},
	}
	payload := encodeDirEntries(dirBase, entries)
	frame := make([]byte, d4InitBlockSize)
	copy(frame[d4FrameHeaderSize:], payload)
	if _, err := w.f.WriteAt(frame, dirBase); err != nil {
		return 0, 0, err
	}
	end, err := w.f.Seek(0, 2)
	if err != nil {
		return 0, 0, err
	}
	return dirBase, end - dirBase, nil
}

// d4DenseDepths materialises the full per-base depth array for one
// accumulator, length == refLen. Bases past the highest event keep depth 0.
func d4DenseDepths(a *covAccum) []int32 {
	n := a.refLen
	if n < 0 {
		n = 0
	}
	out := make([]int32, n)
	a.emit(func(pos int, depth int32) {
		if pos >= 0 && pos < n {
			out[pos] = depth
		}
	})
	return out
}

// d4Reader reads back a D4 file produced by d4Writer (or upstream mosdepth),
// decoding the bit-packed primary table and secondary overflow records. It is
// used by the round-trip parity tests to confirm the depth track reproduces
// the per-base depths exactly.
type d4Reader struct {
	chroms    []d4Chrom
	data      []byte
	ptabStart int64
	overflow  map[string][]d4Overflow
}

// openD4Reader reads the whole file at path into memory and parses its header,
// primary-table location and secondary overflow records.
func openD4Reader(path string) (*d4Reader, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	if len(raw) < 8 || raw[0] != d4Magic[0] || raw[1] != d4Magic[1] || raw[2] != d4Magic[2] || raw[3] != d4Magic[3] {
		return nil, fmt.Errorf("mosdepth: bad D4 magic")
	}
	root, err := parseDirectory(raw, 8, d4InitBlockSize)
	if err != nil {
		return nil, err
	}
	var metaEnt, ptabEnt, stabEnt *dirEntry
	for i := range root {
		switch root[i].name {
		case ".metadata":
			metaEnt = &root[i]
		case ".ptab":
			ptabEnt = &root[i]
		case ".stab":
			stabEnt = &root[i]
		}
	}
	if metaEnt == nil || ptabEnt == nil {
		return nil, fmt.Errorf("mosdepth: D4 file missing .metadata or .ptab")
	}
	metaBytes := readStream(raw, metaEnt.primaryOffset, metaEnt.primarySize)
	chroms, err := parseChromList(metaBytes)
	if err != nil {
		return nil, err
	}
	r := &d4Reader{
		chroms:    chroms,
		data:      raw,
		ptabStart: ptabEnt.primaryOffset,
		overflow:  map[string][]d4Overflow{},
	}
	if stabEnt != nil {
		if err := r.parseStab(stabEnt); err != nil {
			return nil, err
		}
	}
	return r, nil
}

// parseDirectory decodes a framefile directory stream starting at offset with
// the given frame size, returning its entries with absolute offsets.
func parseDirectory(data []byte, offset, frameSize int64) ([]dirEntry, error) {
	payload := readStream(data, offset, frameSize)
	var entries []dirEntry
	i := 0
	for i < len(payload) {
		if payload[i] == 0 {
			break
		}
		i++ // has-next flag
		if i+1+8+8 > len(payload) {
			return nil, fmt.Errorf("mosdepth: truncated D4 directory entry")
		}
		kind := payload[i]
		i++
		off := int64(binary.LittleEndian.Uint64(payload[i : i+8]))
		i += 8
		size := int64(binary.LittleEndian.Uint64(payload[i : i+8]))
		i += 8
		start := i
		for i < len(payload) && payload[i] != 0 {
			i++
		}
		name := string(payload[start:i])
		i++ // name terminator
		entries = append(entries, dirEntry{kind: kind, primaryOffset: off + offset, primarySize: size, name: name})
	}
	return entries, nil
}

// readStream reassembles a framefile stream's payload by following the frame
// link chain, starting at offset with the given primary frame size.
func readStream(data []byte, offset, frameSize int64) []byte {
	var out []byte
	cur := offset
	size := frameSize
	for {
		if cur < 0 || cur+16 > int64(len(data)) {
			break
		}
		link := int64(binary.LittleEndian.Uint64(data[cur : cur+8]))
		end := cur + size
		if end > int64(len(data)) {
			end = int64(len(data))
		}
		out = append(out, data[cur+16:end]...)
		if link == 0 {
			break
		}
		nextSize := int64(binary.LittleEndian.Uint64(data[cur+8 : cur+16]))
		cur += link
		size = nextSize
	}
	return out
}

// parseChromList extracts the chromosome list from the metadata JSON. It does a
// lightweight scan for "name"/"size" pairs rather than pulling in a JSON
// dependency, matching the exact shape this writer (and upstream) emits.
func parseChromList(meta []byte) ([]d4Chrom, error) {
	s := strings.TrimRight(string(meta), "\x00")
	idx := strings.Index(s, `"chrom_list":[`)
	if idx < 0 {
		return nil, fmt.Errorf("mosdepth: D4 metadata missing chrom_list")
	}
	rest := s[idx+len(`"chrom_list":[`):]
	var chroms []d4Chrom
	for {
		ni := strings.Index(rest, `"name":`)
		if ni < 0 {
			break
		}
		rest = rest[ni+len(`"name":`):]
		name, err := strconv.Unquote(scanJSONString(rest))
		if err != nil {
			return nil, fmt.Errorf("mosdepth: D4 metadata bad name: %w", err)
		}
		si := strings.Index(rest, `"size":`)
		if si < 0 {
			break
		}
		rest = rest[si+len(`"size":`):]
		j := 0
		for j < len(rest) && (rest[j] == '-' || (rest[j] >= '0' && rest[j] <= '9')) {
			j++
		}
		sz, err := strconv.ParseInt(rest[:j], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("mosdepth: D4 metadata bad size: %w", err)
		}
		chroms = append(chroms, d4Chrom{Name: name, Length: sz})
		rest = rest[j:]
		if strings.HasPrefix(rest, "]") {
			break
		}
	}
	return chroms, nil
}

// scanJSONString returns the leading JSON string literal (including quotes)
// from s, which must start at a '"'.
func scanJSONString(s string) string {
	if len(s) == 0 || s[0] != '"' {
		return `""`
	}
	for i := 1; i < len(s); i++ {
		if s[i] == '\\' {
			i++
			continue
		}
		if s[i] == '"' {
			return s[:i+1]
		}
	}
	return s
}

// parseStab decodes the secondary table sub-directory into per-chromosome
// overflow records keyed by chromosome name.
func (r *d4Reader) parseStab(stab *dirEntry) error {
	entries, err := parseDirectory(r.data, stab.primaryOffset, d4InitBlockSize)
	if err != nil {
		return err
	}
	for _, e := range entries {
		if e.name == ".metadata" {
			continue
		}
		idx, err := strconv.Atoi(e.name)
		if err != nil || idx < 0 || idx >= len(r.chroms) {
			continue
		}
		payload := readStream(r.data, e.primaryOffset, e.primarySize)
		recs := decodeStabStream(payload)
		if len(recs) > 0 {
			r.overflow[r.chroms[idx].Name] = recs
		}
	}
	return nil
}

// decodeStabStream decodes RangeRecords from a secondary-table stream payload.
// Each record is the 10-byte packed (left:u32, size_enc:u16, value:i32) tuple.
// A record with stored left==0 is invalid and terminates the stream (this is
// how the zero-padded tail of the final frame signals end-of-records). The
// decoded half-open base range is [left-1, left+size_enc).
func decodeStabStream(payload []byte) []d4Overflow {
	var out []d4Overflow
	for off := 0; off+10 <= len(payload); off += 10 {
		left := binary.LittleEndian.Uint32(payload[off : off+4])
		if left == 0 {
			break
		}
		sizeEnc := binary.LittleEndian.Uint16(payload[off+4 : off+6])
		val := int32(binary.LittleEndian.Uint32(payload[off+6 : off+10]))
		actualLeft := left - 1
		out = append(out, d4Overflow{Left: actualLeft, Right: actualLeft + uint32(sizeEnc) + 1, Value: val})
	}
	return out
}

// chromDepths returns the dense per-base depth array for the named chromosome,
// decoded from the bit-packed primary table and the secondary overflow table.
func (r *d4Reader) chromDepths(name string) ([]int32, error) {
	base := r.ptabStart
	for _, c := range r.chroms {
		nbytes := ptabBytesForLen(c.Length)
		if c.Name != name {
			base += nbytes
			continue
		}
		if base+nbytes > int64(len(r.data)) {
			return nil, fmt.Errorf("mosdepth: D4 primary table truncated for %q", name)
		}
		seg := r.data[base : base+nbytes]
		out := make([]int32, c.Length)
		mask := uint32((1 << d4Bits) - 1)
		for i := range out {
			bit := int64(i) * int64(d4Bits)
			byteOff := bit / 8
			shift := uint(bit % 8)
			var v uint32
			v = uint32(seg[byteOff]) >> shift
			if shift+d4Bits > 8 {
				v |= uint32(seg[byteOff+1]) << (8 - shift)
			}
			out[i] = int32(v & mask)
		}
		// Apply overflow records.
		for _, ov := range r.overflow[name] {
			for p := ov.Left; p < ov.Right && int64(p) < c.Length; p++ {
				out[p] = ov.Value
			}
		}
		return out, nil
	}
	return nil, fmt.Errorf("mosdepth: chromosome %q not in D4 file", name)
}
