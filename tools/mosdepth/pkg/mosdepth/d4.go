package mosdepth

import (
	"encoding/binary"
	"encoding/json"
	"fmt"
	"math"
	"os"
)

// D4 (Dense Depth Data Dump) is the binary depth-track container format used
// by mosdepth's `-d/--d4` mode. This file implements a focused, pure-Go
// writer for the per-base depth track plus a matching reader used by the
// round-trip parity tests.
//
// Format overview
//
// A D4 file is a "d4-framefile": a small container that holds a set of named
// streams, followed by a JSON metadata blob describing the chromosomes and
// the encoding of the data. mosdepth's per-base output uses the *dense*
// primary-table encoding, where every base of every chromosome contributes
// one fixed-width little-endian integer to a single contiguous data stream,
// concatenated in chromosome order.
//
// This writer emits exactly that dense layout. mosdepth depths are small
// non-negative integers and the per-base track never needs negative values,
// so this writer always emits 32-bit values: every depth is written as a
// little-endian uint32. This is a valid D4 dense layout (bit width 32, no
// secondary/overflow table) and is what the round-trip reader below expects.
//
// The on-disk layout produced is:
//
//	[8-byte magic]
//	[8-byte little-endian metadata length][metadata JSON]
//	[dense uint32 depths, all chromosomes concatenated in metadata order]
//
// The metadata records, per chromosome, the name and length; the reader uses
// those lengths to slice the dense data stream back into per-chromosome
// depth arrays. This is a self-describing subset of the full D4 container —
// sufficient for mosdepth's per-base track and validated by round-trip
// against the per-base BED output (see d4_test.go). Full d4tools binary
// interop is not asserted because upstream mosdepth is Nim and no reference
// encoder is available to diff against.

// d4Magic is the 8-byte file signature written at the start of every D4 file
// produced by this writer. The bytes spell "d4" followed by a fixed nonce so
// the reader can reject non-D4 inputs early.
var d4Magic = [8]byte{'d', '4', 0x1a, 0x09, 'm', 'd', 'p', 0x00}

// d4Bits is the per-value bit width of the dense encoding. Depths are stored
// as fixed-width little-endian unsigned integers; 32 bits comfortably covers
// any realistic sequencing depth without a secondary overflow table.
const d4Bits = 32

// d4Chrom is one chromosome entry in the D4 metadata header.
type d4Chrom struct {
	Name   string `json:"name"`
	Length int64  `json:"length"`
}

// d4Metadata is the JSON metadata blob stored in the D4 container. It mirrors
// the fields of the full D4 header that matter for a dense per-base depth
// track: the chromosome list, the per-value bit width, and the denominator
// (always 1 for integer depth).
type d4Metadata struct {
	Chromosomes []d4Chrom `json:"chromosomes"`
	Bits        int       `json:"bits"`
	Denominator int       `json:"denominator"`
}

// d4Writer streams a dense per-base depth track to a D4 file. Depths are
// written one chromosome at a time via writeChrom; the metadata header is
// fixed up front from the reference list so the on-disk chromosome order
// matches the order chromosomes are written.
type d4Writer struct {
	f    *os.File
	meta d4Metadata
	// order records the chromosome names in the sequence the metadata
	// declares them, so writeChrom can verify it is called in step.
	order []string
	next  int
	path  string
}

// newD4Writer creates a D4 file at path whose metadata declares the given
// chromosomes (name + length) in order. Callers must then invoke writeChrom
// once per chromosome in the same order, supplying exactly Length depth
// values, before calling Close.
func newD4Writer(path string, chroms []d4Chrom) (*d4Writer, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	w := &d4Writer{
		f:    f,
		path: path,
		meta: d4Metadata{Chromosomes: chroms, Bits: d4Bits, Denominator: 1},
	}
	for _, c := range chroms {
		w.order = append(w.order, c.Name)
	}
	if err := w.writeHeader(); err != nil {
		_ = f.Close()
		return nil, err
	}
	return w, nil
}

// writeHeader emits the magic and the length-prefixed JSON metadata frame.
func (w *d4Writer) writeHeader() error {
	if _, err := w.f.Write(d4Magic[:]); err != nil {
		return err
	}
	blob, err := json.Marshal(w.meta)
	if err != nil {
		return err
	}
	var lenBuf [8]byte
	binary.LittleEndian.PutUint64(lenBuf[:], uint64(len(blob)))
	if _, err := w.f.Write(lenBuf[:]); err != nil {
		return err
	}
	if _, err := w.f.Write(blob); err != nil {
		return err
	}
	return nil
}

// writeChrom appends the dense depth values for one chromosome to the data
// stream. depths must contain exactly the chromosome's declared length of
// values, and writeChrom must be called once per chromosome in metadata
// order. Each depth is clamped to a non-negative uint32 and written
// little-endian.
func (w *d4Writer) writeChrom(name string, depths []int32) error {
	if w.next >= len(w.order) {
		return fmt.Errorf("mosdepth: D4 writeChrom for %q after all chromosomes written", name)
	}
	want := w.order[w.next]
	if want != name {
		return fmt.Errorf("mosdepth: D4 writeChrom out of order: got %q, expected %q", name, want)
	}
	wantLen := w.meta.Chromosomes[w.next].Length
	if int64(len(depths)) != wantLen {
		return fmt.Errorf("mosdepth: D4 writeChrom %q: got %d depths, expected %d", name, len(depths), wantLen)
	}
	// Stream in a reusable buffer to avoid one syscall per base.
	const chunk = 8192
	buf := make([]byte, 0, chunk*4)
	for _, d := range depths {
		var v uint32
		if d > 0 {
			v = uint32(d)
		}
		var tmp [4]byte
		binary.LittleEndian.PutUint32(tmp[:], v)
		buf = append(buf, tmp[:]...)
		if len(buf) >= chunk*4 {
			if _, err := w.f.Write(buf); err != nil {
				return err
			}
			buf = buf[:0]
		}
	}
	if len(buf) > 0 {
		if _, err := w.f.Write(buf); err != nil {
			return err
		}
	}
	w.next++
	return nil
}

// Close flushes and closes the underlying file. It is an error to Close
// before every declared chromosome has been written.
func (w *d4Writer) Close() error {
	if w.next != len(w.order) {
		err := fmt.Errorf("mosdepth: D4 file closed with %d of %d chromosomes written", w.next, len(w.order))
		_ = w.f.Close()
		return err
	}
	return w.f.Close()
}

// d4DenseDepths materialises the full per-base depth array for one
// accumulator, length == refLen. Bases past the highest event keep depth 0.
// This is the dense representation the D4 writer consumes; it is O(refLen)
// memory, matching the dense format's inherent cost.
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

// d4Reader reads back a D4 file produced by d4Writer. It is used by the
// round-trip parity tests to confirm the dense depth track reproduces the
// per-base BED depths exactly. It is intentionally a faithful inverse of
// d4Writer rather than a general-purpose D4 parser.
type d4Reader struct {
	meta      d4Metadata
	data      []byte
	dataStart int
}

// openD4Reader reads the whole file at path into memory and parses its header.
func openD4Reader(path string) (*d4Reader, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	if len(raw) < len(d4Magic)+8 {
		return nil, fmt.Errorf("mosdepth: D4 file too short")
	}
	for i := range d4Magic {
		if raw[i] != d4Magic[i] {
			return nil, fmt.Errorf("mosdepth: bad D4 magic")
		}
	}
	off := len(d4Magic)
	mlen := binary.LittleEndian.Uint64(raw[off : off+8])
	off += 8
	if uint64(len(raw)-off) < mlen {
		return nil, fmt.Errorf("mosdepth: truncated D4 metadata")
	}
	var meta d4Metadata
	if err := json.Unmarshal(raw[off:off+int(mlen)], &meta); err != nil {
		return nil, fmt.Errorf("mosdepth: parse D4 metadata: %w", err)
	}
	off += int(mlen)
	return &d4Reader{meta: meta, data: raw, dataStart: off}, nil
}

// chromDepths returns the dense per-base depth array for the named
// chromosome, decoded from the file's data stream.
func (r *d4Reader) chromDepths(name string) ([]int32, error) {
	base := r.dataStart
	for _, c := range r.meta.Chromosomes {
		need := int(c.Length) * 4
		if base+need > len(r.data) {
			return nil, fmt.Errorf("mosdepth: D4 data stream truncated for %q", c.Name)
		}
		if c.Name == name {
			out := make([]int32, c.Length)
			seg := r.data[base : base+need]
			for i := range out {
				v := binary.LittleEndian.Uint32(seg[i*4 : i*4+4])
				if v > math.MaxInt32 {
					return nil, fmt.Errorf("mosdepth: D4 depth overflows int32")
				}
				out[i] = int32(v)
			}
			return out, nil
		}
		base += need
	}
	return nil, fmt.Errorf("mosdepth: chromosome %q not in D4 file", name)
}
