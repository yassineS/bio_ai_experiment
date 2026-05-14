package samtools

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/tools/tabix/pkg/tabix"
)

// BAIMagic is the 4-byte signature at the head of every well-formed .bai
// index file ("BAI\1").
var BAIMagic = [4]byte{'B', 'A', 'I', 0x01}

// MetaBin is the htslib "pseudo-bin" used to record the (first, last) virtual
// offsets covering a reference and the (mapped, unmapped) record counts. The
// number ((1<<18)-1)/7 + 1 = 37450 is one past the last legal regular bin
// (tabix.MaxBin == 37449).
const MetaBin uint32 = 37450

// ErrBadBAIMagic indicates the input did not start with the BAI magic bytes.
var ErrBadBAIMagic = errors.New("samtools: not a BAI index (bad magic)")

// BAIChunk is one virtual-offset range — exactly two uint64s on disk.
type BAIChunk struct {
	Beg, End uint64
}

// BAIBin is one (binNumber, []chunks) entry within a reference. The pseudo-
// bin (BinID == MetaBin) holds two chunks: the first is
// (firstVOff, lastVOff) and the second is (mapped, unmapped) counts.
type BAIBin struct {
	BinID  uint32
	Chunks []BAIChunk
}

// BAIRef is the per-reference index payload: a list of bins plus the linear
// index of tile-min virtual offsets.
type BAIRef struct {
	Bins   []BAIBin
	Linear []uint64
}

// BAIIndex is the in-memory form of a BAI file.
type BAIIndex struct {
	// Refs is one entry per reference sequence in @SQ order.
	Refs []BAIRef
	// NoCoor is the number of unmapped (refID == -1) records, recorded
	// after the last reference per the BAI spec's optional trailer.
	NoCoor uint64
}

// baiRefBuilder is the in-progress per-reference state used by the BAI
// builder while it streams records.
type baiRefBuilder struct {
	bins   map[uint32]*BAIBin
	linear []uint64
	// firstVOff and lastVOff record the smallest and largest virtual offsets
	// of any mapped record on this reference. The "pseudo-bin" stores them
	// for the runtime to bracket per-ref reads.
	firstVOff uint64
	lastVOff  uint64
	mapped    uint64
	unmapped  uint64
	// haveOff tracks whether we have seen any record on this ref yet, so we
	// can distinguish "first record" from "no records yet".
	haveOff bool
}

// BAIBuilder collects the bin/linear/meta state for an in-progress index.
// Callers add records sequentially in BAM order via AddRecord, then call
// Finish to obtain the populated BAIIndex.
type BAIBuilder struct {
	refs      []*baiRefBuilder
	noCoor    uint64
	totalRefs int
}

// NewBAIBuilder creates an empty builder for an index over numRefs reference
// sequences. The reference IDs in AddRecord refer to positions in [0, numRefs).
func NewBAIBuilder(numRefs int) *BAIBuilder {
	b := &BAIBuilder{totalRefs: numRefs, refs: make([]*baiRefBuilder, numRefs)}
	for i := range b.refs {
		b.refs[i] = newBaiRefBuilder()
	}
	return b
}

func newBaiRefBuilder() *baiRefBuilder {
	return &baiRefBuilder{bins: map[uint32]*BAIBin{}, firstVOff: ^uint64(0)}
}

// AddRecord registers one BAM record with the index. refID is the 0-based
// reference index (or -1 for unmapped); beg/end are the 0-based half-open
// reference coordinates that the record spans (end - beg == CIGAR
// reference length, with end clamped to beg+1 for empty CIGARs). vBeg is the
// virtual offset of the record's first byte and vEnd is the virtual offset
// just past the record (i.e. the start of the next record).
//
// Records that are unmapped (refID == -1) bump the trailing n_no_coor
// counter. Records with refID >= 0 but the unmapped flag (mapped == false)
// still contribute to the per-ref unmapped meta count.
func (b *BAIBuilder) AddRecord(refID, beg, end int, vBeg, vEnd uint64, mapped bool) error {
	if refID < 0 {
		b.noCoor++
		return nil
	}
	if refID >= b.totalRefs {
		return fmt.Errorf("samtools: BAI refID %d >= %d", refID, b.totalRefs)
	}
	if end <= beg {
		end = beg + 1
	}
	ref := b.refs[refID]

	// Update meta bookkeeping.
	if !ref.haveOff || vBeg < ref.firstVOff {
		ref.firstVOff = vBeg
		ref.haveOff = true
	}
	if vEnd > ref.lastVOff {
		ref.lastVOff = vEnd
	}
	if mapped {
		ref.mapped++
	} else {
		ref.unmapped++
	}

	// Bin lookup + chunk merge.
	binID := uint32(tabix.Reg2bin(beg, end))
	bin, ok := ref.bins[binID]
	if !ok {
		bin = &BAIBin{BinID: binID}
		ref.bins[binID] = bin
	}
	if n := len(bin.Chunks); n > 0 && bin.Chunks[n-1].End >= vBeg {
		// Merge contiguous / overlapping chunks. This both keeps the chunk
		// list small and matches the htslib convention.
		if vEnd > bin.Chunks[n-1].End {
			bin.Chunks[n-1].End = vEnd
		}
	} else {
		bin.Chunks = append(bin.Chunks, BAIChunk{Beg: vBeg, End: vEnd})
	}

	// Linear index update — every 16-Kbp tile spanned by the record gets
	// the smallest virtual offset of any record touching it.
	first := tabix.LinearTile(beg)
	last := tabix.LinearTile(end - 1)
	if last < first {
		last = first
	}
	for last >= len(ref.linear) {
		ref.linear = append(ref.linear, ^uint64(0))
	}
	for t := first; t <= last; t++ {
		if ref.linear[t] == ^uint64(0) || vBeg < ref.linear[t] {
			ref.linear[t] = vBeg
		}
	}
	return nil
}

// Finish finalises the linear indices (carry-forward the last seen offset
// over sentinel slots, per htslib) and returns the assembled BAIIndex.
func (b *BAIBuilder) Finish() *BAIIndex {
	idx := &BAIIndex{NoCoor: b.noCoor, Refs: make([]BAIRef, b.totalRefs)}
	for i, ref := range b.refs {
		// Linear-index finalize.
		lin := ref.linear
		var last uint64
		seen := false
		for j := range lin {
			if lin[j] == ^uint64(0) {
				if seen {
					lin[j] = last
				} else {
					lin[j] = 0
				}
			} else {
				last = lin[j]
				seen = true
			}
		}
		// Bin slice — sort by ID for deterministic output.
		bins := make([]BAIBin, 0, len(ref.bins)+1)
		for _, bb := range ref.bins {
			bins = append(bins, *bb)
		}
		sort.Slice(bins, func(a, b int) bool { return bins[a].BinID < bins[b].BinID })
		// Inject the meta pseudo-bin if we saw at least one record on this
		// ref. htslib emits it unconditionally for refs touched by any
		// record (mapped or unmapped on this ref).
		if ref.haveOff || ref.mapped > 0 || ref.unmapped > 0 {
			meta := BAIBin{BinID: MetaBin, Chunks: []BAIChunk{
				{Beg: ref.firstVOff, End: ref.lastVOff},
				{Beg: ref.mapped, End: ref.unmapped},
			}}
			// Insert in sorted position (MetaBin > all regular bins, so
			// just append).
			bins = append(bins, meta)
		}
		idx.Refs[i] = BAIRef{Bins: bins, Linear: lin}
	}
	return idx
}

// WriteBAI serialises idx to w in the canonical .bai byte layout (little-
// endian throughout).
func WriteBAI(w io.Writer, idx *BAIIndex) error {
	bw := bufio.NewWriter(w)
	if _, err := bw.Write(BAIMagic[:]); err != nil {
		return err
	}
	if err := binary.Write(bw, binary.LittleEndian, int32(len(idx.Refs))); err != nil {
		return err
	}
	for _, ref := range idx.Refs {
		if err := binary.Write(bw, binary.LittleEndian, int32(len(ref.Bins))); err != nil {
			return err
		}
		for _, bin := range ref.Bins {
			if err := binary.Write(bw, binary.LittleEndian, bin.BinID); err != nil {
				return err
			}
			if err := binary.Write(bw, binary.LittleEndian, int32(len(bin.Chunks))); err != nil {
				return err
			}
			for _, c := range bin.Chunks {
				if err := binary.Write(bw, binary.LittleEndian, c.Beg); err != nil {
					return err
				}
				if err := binary.Write(bw, binary.LittleEndian, c.End); err != nil {
					return err
				}
			}
		}
		if err := binary.Write(bw, binary.LittleEndian, int32(len(ref.Linear))); err != nil {
			return err
		}
		for _, v := range ref.Linear {
			if err := binary.Write(bw, binary.LittleEndian, v); err != nil {
				return err
			}
		}
	}
	// Optional trailer — htslib emits n_no_coor whenever it is non-zero.
	if idx.NoCoor > 0 {
		if err := binary.Write(bw, binary.LittleEndian, idx.NoCoor); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// ReadBAI parses a .bai byte stream into a fresh BAIIndex. The reader is
// expected to be raw (not BGZF-compressed) — BAI files are stored as plain
// little-endian bytes.
func ReadBAI(r io.Reader) (*BAIIndex, error) {
	br := bufio.NewReader(r)
	var magic [4]byte
	if _, err := io.ReadFull(br, magic[:]); err != nil {
		return nil, err
	}
	if magic != BAIMagic {
		return nil, ErrBadBAIMagic
	}
	var nRef int32
	if err := binary.Read(br, binary.LittleEndian, &nRef); err != nil {
		return nil, err
	}
	if nRef < 0 {
		return nil, fmt.Errorf("samtools: BAI negative n_ref %d", nRef)
	}
	idx := &BAIIndex{Refs: make([]BAIRef, nRef)}
	for i := int32(0); i < nRef; i++ {
		var nBin int32
		if err := binary.Read(br, binary.LittleEndian, &nBin); err != nil {
			return nil, err
		}
		if nBin < 0 {
			return nil, fmt.Errorf("samtools: BAI negative n_bin %d", nBin)
		}
		bins := make([]BAIBin, nBin)
		for j := range bins {
			if err := binary.Read(br, binary.LittleEndian, &bins[j].BinID); err != nil {
				return nil, err
			}
			var nChunk int32
			if err := binary.Read(br, binary.LittleEndian, &nChunk); err != nil {
				return nil, err
			}
			if nChunk < 0 {
				return nil, fmt.Errorf("samtools: BAI negative n_chunk %d", nChunk)
			}
			chunks := make([]BAIChunk, nChunk)
			for k := range chunks {
				if err := binary.Read(br, binary.LittleEndian, &chunks[k].Beg); err != nil {
					return nil, err
				}
				if err := binary.Read(br, binary.LittleEndian, &chunks[k].End); err != nil {
					return nil, err
				}
			}
			bins[j].Chunks = chunks
		}
		var nIntv int32
		if err := binary.Read(br, binary.LittleEndian, &nIntv); err != nil {
			return nil, err
		}
		if nIntv < 0 {
			return nil, fmt.Errorf("samtools: BAI negative n_intv %d", nIntv)
		}
		linear := make([]uint64, nIntv)
		for k := range linear {
			if err := binary.Read(br, binary.LittleEndian, &linear[k]); err != nil {
				return nil, err
			}
		}
		idx.Refs[i] = BAIRef{Bins: bins, Linear: linear}
	}
	// Optional trailer.
	var trailer uint64
	if err := binary.Read(br, binary.LittleEndian, &trailer); err == nil {
		idx.NoCoor = trailer
	}
	return idx, nil
}

// MetaCounts returns the (mapped, unmapped) counts stored in the BAI meta
// pseudo-bin for the given reference, or (0, 0, false) if the pseudo-bin
// is absent or the refID is out of range.
func (idx *BAIIndex) MetaCounts(refID int) (mapped, unmapped uint64, ok bool) {
	if refID < 0 || refID >= len(idx.Refs) {
		return 0, 0, false
	}
	for _, bin := range idx.Refs[refID].Bins {
		if bin.BinID == MetaBin && len(bin.Chunks) >= 2 {
			return bin.Chunks[1].Beg, bin.Chunks[1].End, true
		}
	}
	return 0, 0, false
}

// MetaBounds returns the (firstVOff, lastVOff) virtual offsets recorded in
// the BAI meta pseudo-bin for the given reference, or (0, 0, false) when
// absent.
func (idx *BAIIndex) MetaBounds(refID int) (first, last uint64, ok bool) {
	if refID < 0 || refID >= len(idx.Refs) {
		return 0, 0, false
	}
	for _, bin := range idx.Refs[refID].Bins {
		if bin.BinID == MetaBin && len(bin.Chunks) >= 1 {
			return bin.Chunks[0].Beg, bin.Chunks[0].End, true
		}
	}
	return 0, 0, false
}

// RegionChunks returns the set of merged chunks that may contain records
// overlapping [beg, end) on the given reference. beg and end are 0-based
// half-open. Returns nil chunks when the refID is out of range.
func (idx *BAIIndex) RegionChunks(refID, beg, end int) []BAIChunk {
	if refID < 0 || refID >= len(idx.Refs) || end <= beg {
		return nil
	}
	if beg < 0 {
		beg = 0
	}
	ref := idx.Refs[refID]

	// Linear-index lower bound.
	var minOff uint64
	tile := tabix.LinearTile(beg)
	if tile < len(ref.Linear) {
		minOff = ref.Linear[tile]
	}

	bins := tabix.Reg2bins(beg, end)
	want := make(map[uint32]struct{}, len(bins))
	for _, b := range bins {
		want[uint32(b)] = struct{}{}
	}

	var chunks []BAIChunk
	for _, bin := range ref.Bins {
		if bin.BinID == MetaBin {
			continue
		}
		if _, ok := want[bin.BinID]; !ok {
			continue
		}
		for _, c := range bin.Chunks {
			if c.End <= minOff {
				continue
			}
			cc := c
			if cc.Beg < minOff {
				cc.Beg = minOff
			}
			chunks = append(chunks, cc)
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
