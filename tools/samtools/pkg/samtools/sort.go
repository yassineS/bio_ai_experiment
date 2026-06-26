package samtools

import (
	"container/heap"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// SortOrder names the sort key for the output records.
type SortOrder int

const (
	// SortCoordinate sorts by (refID, 0-based pos, qname-as-tiebreaker).
	SortCoordinate SortOrder = iota
	// SortByName sorts by QName using a plain (case-sensitive) byte compare.
	SortByName
	// SortByNameNatural sorts by QName with numeric runs compared as
	// integers ("r2" < "r10").
	SortByNameNatural
	// SortByTag sorts by the value of an aux tag. Integer tags compare as
	// numbers; string tags compare lexicographically; floats compare
	// numerically; absent tags sort last.
	SortByTag
)

// SortOptions configures the external-merge sort pipeline.
type SortOptions struct {
	// Order picks the sort key.
	Order SortOrder
	// Tag is the two-letter aux tag used when Order == SortByTag.
	Tag string
	// MaxMemBytes is the per-shard memory budget. When the running estimate
	// of in-memory record bytes crosses this threshold the shard is sorted
	// and spilled to a temp file. 0 means "use SortDefaultMem".
	MaxMemBytes int64
	// TmpPrefix is the path prefix used for spill files. Empty selects the
	// OS temp directory.
	TmpPrefix string
	// CompressLevel selects the BGZF deflate level for the output file.
	// 0..9 are honoured; -1 (or any value out of range) selects the bgzip
	// default.
	CompressLevel int
	// OutputBAM forces BAM output; otherwise output format is chosen by the
	// output file extension (.sam / .bam) or defaults to BAM.
	OutputBAM bool
	// OutputSAM forces text SAM output (overrides OutputBAM when both set —
	// SAM wins).
	OutputSAM bool
	// Threads sets the BGZF compression worker count from `-@/--threads`. A
	// value above 1 compresses both the temporary sort shards and the final
	// BAM output in parallel; the decoded records are byte-identical to the
	// single-threaded path. It has no effect on SAM output.
	Threads int
	// NoPG is accepted but is a no-op (we don't inject @PG lines anyway).
	NoPG bool
}

// SortDefaultMem is the default per-shard memory budget when SortOptions.
// MaxMemBytes is zero. 768 MiB matches upstream samtools sort's default.
const SortDefaultMem int64 = 768 << 20

// ErrEmptyTag is returned when SortByTag is selected but Tag is empty or
// not exactly 2 characters long.
var ErrEmptyTag = errors.New("samtools sort: --by-tag requires a 2-character tag")

// Sort reads BAM/SAM records from in, sorts them per opts, and writes them
// out as BAM (or SAM if opts.OutputSAM is set). It performs an external-
// merge sort so the working set never exceeds opts.MaxMemBytes.
func Sort(in io.Reader, out io.Writer, opts SortOptions) error {
	if opts.Order == SortByTag && len(opts.Tag) != 2 {
		return ErrEmptyTag
	}
	if opts.MaxMemBytes <= 0 {
		opts.MaxMemBytes = SortDefaultMem
	}

	// The per-shard spill budget is the flat MaxMemBytes (default 768 MiB),
	// deliberately NOT scaled by -@. Scaling the BYTE budget by thread count
	// (mirroring upstream's `max_mem = _max_mem * n_threads`) was tried and
	// regressed peak RSS to ~2x upstream: a buffered Go *sam.Record's real heap
	// footprint is several times its packed byte size, so a thread-scaled byte
	// budget lets the live buffer's resident set blow past upstream's, breaking
	// task #32's memory bound — and it did not improve the full-sort wall (the
	// per-spill decode/re-encode round-trip, not the shard count, dominates).
	// Keeping the flat budget holds peak RSS within bound. Reducing the
	// full-sort shard count without inflating RSS needs a packed/arena spill
	// format that avoids re-decoding records — tracked as a follow-up in
	// docs/PARITY_ROADMAP.md.
	budget := opts.MaxMemBytes

	// Decode the input with block-parallel BGZF inflation when -@/--threads
	// asks for it. sort retains every record (it buffers and spills them), so
	// it reads via Read() — which returns a fresh record each call — and never
	// uses the allocation-free ReadInto fast path; the threaded reader is a
	// drop-in here. The decoded record stream is byte-for-byte identical for any
	// thread count, so the sorted output is unchanged; -@ only overlaps input
	// inflation with the sort/spill work. opts.Threads is honoured verbatim
	// (not widened to NumCPU): with no -@ it is 0/1 and NewReaderThreaded falls
	// back to the existing single-threaded sam.NewReader path, so the default is
	// behaviour-, output-, and performance-identical to before.
	r, err := alnio.NewReaderThreaded(in, "", opts.Threads)
	if err != nil {
		return err
	}
	if rc, ok := r.(io.Closer); ok {
		defer rc.Close()
	}
	hdr := r.Header()
	// Build a reference-name → refID map for coordinate sort.
	refIndex := make(map[string]int, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refIndex[ref.Name] = i
	}

	tmpDir := opts.TmpPrefix
	if tmpDir == "" {
		tmpDir = os.TempDir()
	}
	// If TmpPrefix is a file-path prefix we treat its directory portion as
	// the tmpdir and reserve the basename for the file prefix.
	tmpBaseDir := tmpDir
	tmpName := "sort"
	if fi, ferr := os.Stat(tmpDir); ferr != nil || !fi.IsDir() {
		tmpBaseDir = filepath.Dir(tmpDir)
		tmpName = filepath.Base(tmpDir)
		if tmpBaseDir == "" || tmpBaseDir == "." {
			tmpBaseDir = "."
		}
	}

	// PACKED FAST PATH. When the input is BGZF-wrapped (or raw) BAM the reader
	// exposes ReadRaw, which hands back the verbatim on-disk record body. The
	// sort path never mutates a record (RG/PG translation is merge-only; Sort's
	// writeOutput applies none), so we buffer/spill/merge those raw bodies and
	// emit them unchanged — eliminating the per-record decode→re-encode round
	// trip that dominated the profile. The first raw read also doubles as the
	// probe: ErrNoRawRead (returned without consuming input) means SAM/CRAM
	// input, so we fall through to the decode path below with the stream intact.
	if rr, ok := r.(rawReader); ok {
		first, ferr := rr.ReadRaw()
		if ferr == nil || ferr == io.EOF {
			return sortPacked(rr, first, ferr, out, hdr, opts, budget, tmpBaseDir, tmpName)
		}
		if ferr != alnio.ErrNoRawRead {
			return ferr
		}
		// ErrNoRawRead: no bytes consumed, fall through to the decode path.
	}

	cmp := makeRecordLess(opts, refIndex)

	var shards []string
	defer func() {
		for _, p := range shards {
			_ = os.Remove(p)
		}
	}()

	var buffer []*sam.Record
	var bufBytes int64
	flush := func() error {
		if len(buffer) == 0 {
			return nil
		}
		sort.SliceStable(buffer, func(a, b int) bool { return cmp(buffer[a], buffer[b]) })
		path, ferr := writeShard(tmpBaseDir, tmpName, len(shards), hdr, buffer)
		if ferr != nil {
			return ferr
		}
		shards = append(shards, path)
		buffer = buffer[:0]
		bufBytes = 0
		return nil
	}

	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		size := recordSize(rec)
		// If we are about to overflow and we already have content, spill
		// first so the new record always lands in a fresh buffer.
		if len(buffer) > 0 && bufBytes+size > budget {
			if err := flush(); err != nil {
				return err
			}
		}
		buffer = append(buffer, rec)
		bufBytes += size
	}
	// Fast path: zero shards → sort the buffer in place and write to out.
	if len(shards) == 0 {
		sort.SliceStable(buffer, func(a, b int) bool { return cmp(buffer[a], buffer[b]) })
		return writeOutput(out, hdr, sliceIter(buffer), opts)
	}
	// Otherwise: flush the tail and k-way merge.
	if err := flush(); err != nil {
		return err
	}
	readers := make([]*sam.BAMReader, 0, len(shards))
	defer func() {
		for _, br := range readers {
			_ = br.Close()
		}
	}()
	for _, p := range shards {
		f, oerr := os.Open(p)
		if oerr != nil {
			return oerr
		}
		// Hold the file open via the reader's underlying io.Reader; closing
		// the reader does not close the file, so we close the file when we
		// hit EOF inside the iterator.
		br, oerr := sam.NewBAMReader(f)
		if oerr != nil {
			_ = f.Close()
			return oerr
		}
		readers = append(readers, br)
	}
	merged := newMergeIterator(readers, cmp)
	return writeOutput(out, hdr, merged, opts)
}

// sortPacked is the packed-spill external-merge sort: it buffers verbatim BAM
// record bodies (packedRec), spills sorted shards of raw bodies, and merges them
// straight to the output, decoding only each record's sort key — never the full
// record. firstRaw/firstErr carry the result of the probe ReadRaw the caller
// already performed (firstErr is nil for a record or io.EOF for empty input).
//
// The sorted output is byte-identical to the decode path: the comparators are
// field-for-field mirrors of coordCmp/nameLess/tagCmp, the in-memory sort is
// sort.SliceStable and the merge heap tie-breaks on source index, so equal-key
// records keep input order exactly as the *sam.Record path does. For SAM text
// output (not perf-critical) the raw bodies are decoded on the way out.
func sortPacked(rr rawReader, firstRaw []byte, firstErr error, out io.Writer, hdr *sam.Header, opts SortOptions, budget int64, tmpBaseDir, tmpName string) error {
	cmp := makePackedCmp(opts)

	var shards []string
	defer func() {
		for _, p := range shards {
			_ = os.Remove(p)
		}
	}()

	var buffer []packedRec
	var bufBytes int64
	flush := func() error {
		if len(buffer) == 0 {
			return nil
		}
		sort.SliceStable(buffer, func(a, b int) bool { return cmp(&buffer[a], &buffer[b]) < 0 })
		path, ferr := writePackedShard(tmpBaseDir, tmpName, len(shards), hdr, buffer)
		if ferr != nil {
			return ferr
		}
		shards = append(shards, path)
		buffer = buffer[:0]
		bufBytes = 0
		return nil
	}

	// add ingests one raw record body into the buffer, spilling first when it
	// would overflow the budget so the new record lands in a fresh buffer.
	add := func(raw []byte) error {
		p, perr := packRecord(raw)
		if perr != nil {
			return perr
		}
		size := packedSize(raw)
		if len(buffer) > 0 && bufBytes+size > budget {
			if err := flush(); err != nil {
				return err
			}
		}
		buffer = append(buffer, p)
		bufBytes += size
		return nil
	}

	// Consume the probe record, then the rest of the stream.
	if firstErr == nil {
		if err := add(firstRaw); err != nil {
			return err
		}
	}
	if firstErr == nil {
		for {
			raw, err := rr.ReadRaw()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err := add(raw); err != nil {
				return err
			}
		}
	}

	// Fast path: zero shards → sort the buffer in place and write to out.
	if len(shards) == 0 {
		sort.SliceStable(buffer, func(a, b int) bool { return cmp(&buffer[a], &buffer[b]) < 0 })
		return writePackedOutput(out, hdr, &packedSliceIter{recs: buffer}, opts)
	}
	// Otherwise: flush the tail and k-way merge.
	if err := flush(); err != nil {
		return err
	}
	readers := make([]*sam.BAMReader, 0, len(shards))
	defer func() {
		for _, br := range readers {
			_ = br.Close()
		}
	}()
	for _, p := range shards {
		f, oerr := os.Open(p)
		if oerr != nil {
			return oerr
		}
		br, oerr := sam.NewBAMReader(f)
		if oerr != nil {
			_ = f.Close()
			return oerr
		}
		readers = append(readers, br)
	}
	merged, merr := newPackedMergeIterator(readers, cmp)
	if merr != nil {
		return merr
	}
	return writePackedOutput(out, hdr, merged, opts)
}

// rawReader is the optional interface a sam.Reader exposes when it can hand back
// undecoded BAM record bodies. The BAM reader (and the alnio wrappers around it)
// implement it; SAM and CRAM readers do not, so sort detects the packed fast
// path by type-asserting to this interface and probing its first read.
type rawReader interface {
	ReadRaw() ([]byte, error)
}

// packedRec is a buffered/spilled record on the packed sort fast path: the
// verbatim on-disk BAM record body plus the few sort-key fields decoded once
// at read time. It carries NO *sam.Record — the comparator works straight off
// the cached key fields and the raw bytes, and the raw bytes are emitted to the
// output unchanged, eliminating the per-record decode→re-encode round-trip.
//
// refID/pos are the on-disk BAM values (refID == -1 and pos == -1 are the
// unmapped/unset sentinels). name is a sub-slice of raw (the read name WITHOUT
// its trailing NUL) used by the -n/-N name comparators. For -t the matching tag
// is decoded lazily from raw's aux region by the comparator (see packedTagCmp).
type packedRec struct {
	raw   []byte
	refID int32
	pos   int32
	flag  uint16
	name  []byte
}

// packRecord builds a packedRec from an owned BAM record body, decoding only the
// fixed-offset sort-key fields (refID@0, pos@4, flag@14) and slicing out the
// read name (read_name@32, l_read_name@8 bytes minus the trailing NUL). It does
// NOT decode CIGAR/SEQ/QUAL/aux — those bytes ride along in raw untouched.
func packRecord(raw []byte) (packedRec, error) {
	if len(raw) < 32 {
		return packedRec{}, fmt.Errorf("samtools sort: BAM record body too small (%d)", len(raw))
	}
	refID := int32(binary.LittleEndian.Uint32(raw[0:4]))
	pos := int32(binary.LittleEndian.Uint32(raw[4:8]))
	lReadName := int(raw[8])
	flag := binary.LittleEndian.Uint16(raw[14:16])
	nameEnd := 32 + lReadName
	if nameEnd > len(raw) {
		return packedRec{}, fmt.Errorf("samtools sort: truncated read name")
	}
	name := raw[32:nameEnd]
	// Drop the trailing NUL the BAM read name carries, mirroring decodeInto.
	if len(name) > 0 && name[len(name)-1] == 0 {
		name = name[:len(name)-1]
	}
	return packedRec{raw: raw, refID: refID, pos: pos, flag: flag, name: name}, nil
}

// packedRefRank maps an on-disk BAM refID to the rank the coordinate comparator
// orders on: real refs keep their id, the unmapped sentinel (-1) becomes a value
// that sorts after every real reference — identical to the refIndex map's
// unmapped fallback (1<<30) used by the *sam.Record comparator.
func packedRefRank(refID int32) int {
	if refID < 0 {
		return 1 << 30
	}
	return int(refID)
}

// packedCoordCmp reproduces coordCmp field-for-field on packed records: order by
// reference rank, then 0-based position (monotonic with the 1-based pos coordCmp
// compares), then the reverse-strand FLAG bit (forward before reverse). Returns
// -1, 0, or 1.
func packedCoordCmp(a, b *packedRec) int {
	ra := packedRefRank(a.refID)
	rb := packedRefRank(b.refID)
	if ra != rb {
		if ra < rb {
			return -1
		}
		return 1
	}
	if a.pos != b.pos {
		if a.pos < b.pos {
			return -1
		}
		return 1
	}
	arev := a.flag&sam.FlagReverse != 0
	brev := b.flag&sam.FlagReverse != 0
	if arev != brev {
		if !arev {
			return -1
		}
		return 1
	}
	return 0
}

// packedNameCmp reproduces nameLess on packed records: QNAME (natural-numeric or
// plain byte compare) then the FLAG-derived secondary key. Returns -1, 0, or 1.
func packedNameCmp(a, b *packedRec, natural bool) int {
	c := strnumCmpBytes(a.name, b.name, natural)
	if c != 0 {
		return c
	}
	ka := flagSortKey(a.flag)
	kb := flagSortKey(b.flag)
	if ka != kb {
		if ka < kb {
			return -1
		}
		return 1
	}
	return 0
}

// packedTagCmp reproduces tagCmp on packed records for `samtools sort -t TAG`.
// It decodes only the matching tag from each record's raw aux region (via
// sam.FindRawAuxTag), then mirrors tagCmp's normalised-type/value ordering and
// its coordinate fallback exactly. On a malformed aux stream it errors via the
// returned bool=false sentinel — callers using this on the packed path have
// already validated each record body's fixed region, so a decode error here is
// surfaced by treating the tag as absent (which matches "tags equal → coord
// fallback" only when both error; to stay safe we fall back to coord). Returns
// -1, 0, or 1.
func packedTagCmp(a, b *packedRec, tag string) int {
	av, aok := packedGetAux(a, tag)
	bv, bok := packedGetAux(b, tag)
	switch {
	case !aok && bok:
		return -1
	case aok && !bok:
		return 1
	case !aok && !bok:
		return packedCoordCmp(a, b)
	}

	atype := normalizeAuxType(av.Type)
	btype := normalizeAuxType(bv.Type)
	if atype != btype {
		if atype == 'c' && btype == 'f' {
			atype = 'f'
		} else if atype == 'f' && btype == 'c' {
			btype = 'f'
		} else {
			if atype < btype {
				return -1
			}
			return 1
		}
	}

	switch atype {
	case 'c':
		ai, _ := av.Int()
		bi, _ := bv.Int()
		if ai != bi {
			if ai < bi {
				return -1
			}
			return 1
		}
	case 'f':
		af := auxFloat(av)
		bf := auxFloat(bv)
		if af != bf {
			if af < bf {
				return -1
			}
			return 1
		}
	case 'A':
		as, _ := av.Value.(string)
		bs, _ := bv.Value.(string)
		ac, bc := byte(0), byte(0)
		if len(as) > 0 {
			ac = as[0]
		}
		if len(bs) > 0 {
			bc = bs[0]
		}
		if ac != bc {
			if ac < bc {
				return -1
			}
			return 1
		}
	case 'H':
		as, _ := av.Value.(string)
		bs, _ := bv.Value.(string)
		if as != bs {
			if as < bs {
				return -1
			}
			return 1
		}
	}
	return packedCoordCmp(a, b)
}

// packedGetAux decodes the single aux field with the given tag from a packed
// record's raw aux region. A malformed aux stream or a missing tag both report
// "absent" (ok=false); the comparator then orders the record exactly as tagCmp
// orders a tag-less record (tags-absent sort first, with a coordinate fallback
// when both are absent).
func packedGetAux(p *packedRec, tag string) (sam.Aux, bool) {
	off, err := sam.RawAuxOffset(p.raw)
	if err != nil {
		return sam.Aux{}, false
	}
	a, ok, ferr := sam.FindRawAuxTag(p.raw[off:], tag)
	if ferr != nil || !ok {
		return sam.Aux{}, false
	}
	return a, true
}

// strnumCmpBytes is strnumCmp over byte slices (the packed read name is a raw
// sub-slice of the BAM body, avoiding a per-record string allocation). It is a
// byte-for-byte mirror of strnumCmp; the two share their semantics exactly.
func strnumCmpBytes(a, b []byte, natural bool) int {
	if !natural {
		return bytesCompare(a, b)
	}
	pa, pb := 0, 0
	for pa < len(a) && pb < len(b) {
		if !isDigit(a[pa]) || !isDigit(b[pb]) {
			if a[pa] != b[pb] {
				return int(a[pa]) - int(b[pb])
			}
			pa++
			pb++
			continue
		}
		for pa < len(a) && a[pa] == '0' {
			pa++
		}
		for pb < len(b) && b[pb] == '0' {
			pb++
		}
		for pa < len(a) && pb < len(b) && isDigit(a[pa]) && a[pa] == b[pb] {
			pa++
			pb++
		}
		var ca, cb int
		if pa < len(a) {
			ca = int(a[pa])
		}
		if pb < len(b) {
			cb = int(b[pb])
		}
		diff := ca - cb
		for pa < len(a) && pb < len(b) && isDigit(a[pa]) && isDigit(b[pb]) {
			pa++
			pb++
		}
		if pa < len(a) && isDigit(a[pa]) {
			return 1
		}
		if pb < len(b) && isDigit(b[pb]) {
			return -1
		}
		if diff != 0 {
			return diff
		}
	}
	if pa < len(a) {
		return 1
	}
	if pb < len(b) {
		return -1
	}
	return 0
}

// bytesCompare returns the sign of a lexicographic byte comparison, matching the
// plain (non-natural) branch of strnumCmp where Go's string `<`/`>` reduces to a
// byte-wise compare.
func bytesCompare(a, b []byte) int {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			if a[i] < b[i] {
				return -1
			}
			return 1
		}
	}
	if len(a) < len(b) {
		return -1
	}
	if len(a) > len(b) {
		return 1
	}
	return 0
}

// makePackedLess builds the three-way packed comparator for opts, mirroring
// makeRecordLess. It returns a func reporting a<0 / 0 / a>0 so both the stable
// in-memory sort (as a bool via <0) and the merge heap's equal-key tie-break can
// share one comparator.
func makePackedCmp(opts SortOptions) func(a, b *packedRec) int {
	switch opts.Order {
	case SortByName:
		return func(a, b *packedRec) int { return packedNameCmp(a, b, false) }
	case SortByNameNatural:
		return func(a, b *packedRec) int { return packedNameCmp(a, b, true) }
	case SortByTag:
		tag := opts.Tag
		return func(a, b *packedRec) int { return packedTagCmp(a, b, tag) }
	default:
		return func(a, b *packedRec) int { return packedCoordCmp(a, b) }
	}
}

// packedSize estimates a buffered raw record's real Go heap footprint, used to
// decide when to spill on the packed path. A packedRec holds one heap allocation
// — the raw []byte body — plus the small struct/slot bookkeeping; it carries no
// expanded *sam.Record, so the footprint is close to the on-disk record size.
// This is why the packed path packs far more records per shard than the fat
// *sam.Record budget: the same 768 MiB holds many more raw bodies, cutting the
// shard count without raising peak RSS above upstream's packed bam1_t budget.
func packedSize(raw []byte) int64 {
	const structBytes = 64 // packedRec struct + its slot in the []packedRec buffer
	n := (len(raw) + 15) &^ 15
	return int64(structBytes + n)
}

// packedIter is the packed-path analogue of recordIter: it yields the next raw
// BAM record body (and the input index it came from) or (nil, _, EOF) when
// exhausted. The src index lets the merge tie-break on input order so equal-key
// records keep their relative position, exactly as the *sam.Record merge does.
type packedIter interface {
	Next() (packedRec, error)
}

// packedSliceIter iterates an in-memory []packedRec (the zero-shard fast path).
type packedSliceIter struct {
	recs []packedRec
	i    int
}

func (s *packedSliceIter) Next() (packedRec, error) {
	if s.i >= len(s.recs) {
		return packedRec{}, io.EOF
	}
	r := s.recs[s.i]
	s.i++
	return r, nil
}

// packedMergeIterator does a k-way heap-merge across n BAM shard readers, pulling
// raw record bodies (ReadRaw) and decoding only each record's sort key for the
// heap comparison. The merged bytes are emitted verbatim by writePackedOutput.
type packedMergeIterator struct {
	h *packedHeap
}

func newPackedMergeIterator(readers []*sam.BAMReader, cmp func(a, b *packedRec) int) (packedIter, error) {
	h := &packedHeap{cmp: cmp}
	heap.Init(h)
	for i, rd := range readers {
		raw, err := rd.ReadRaw()
		if err == io.EOF {
			continue
		}
		if err != nil {
			return nil, err
		}
		p, perr := packRecord(raw)
		if perr != nil {
			return nil, perr
		}
		heap.Push(h, packedEntry{rec: p, src: i, reader: rd})
	}
	return &packedMergeIterator{h: h}, nil
}

func (m *packedMergeIterator) Next() (packedRec, error) {
	if m.h.Len() == 0 {
		return packedRec{}, io.EOF
	}
	top := heap.Pop(m.h).(packedEntry)
	raw, err := top.reader.ReadRaw()
	if err == nil {
		p, perr := packRecord(raw)
		if perr != nil {
			return packedRec{}, perr
		}
		heap.Push(m.h, packedEntry{rec: p, src: top.src, reader: top.reader})
	} else if err != io.EOF {
		return top.rec, err
	}
	return top.rec, nil
}

type packedEntry struct {
	rec    packedRec
	src    int
	reader *sam.BAMReader
}

type packedHeap struct {
	items []packedEntry
	cmp   func(a, b *packedRec) int
}

func (h *packedHeap) Len() int { return len(h.items) }
func (h *packedHeap) Less(i, j int) bool {
	c := h.cmp(&h.items[i].rec, &h.items[j].rec)
	if c != 0 {
		return c < 0
	}
	// Tie-break on source index so the merge is deterministic and equal-key
	// records keep input order, mirroring recordHeap.Less.
	return h.items[i].src < h.items[j].src
}
func (h *packedHeap) Swap(i, j int)      { h.items[i], h.items[j] = h.items[j], h.items[i] }
func (h *packedHeap) Push(x interface{}) { h.items = append(h.items, x.(packedEntry)) }
func (h *packedHeap) Pop() interface{} {
	n := len(h.items)
	x := h.items[n-1]
	h.items = h.items[:n-1]
	return x
}

// recordIter is a lazy iterator that returns the next record or (nil, EOF)
// when exhausted.
type recordIter interface {
	Next() (*sam.Record, error)
}

type sliceIterT struct {
	recs []*sam.Record
	i    int
}

func sliceIter(recs []*sam.Record) recordIter { return &sliceIterT{recs: recs} }

func (s *sliceIterT) Next() (*sam.Record, error) {
	if s.i >= len(s.recs) {
		return nil, io.EOF
	}
	r := s.recs[s.i]
	s.i++
	return r, nil
}

// mergeIterator does a k-way heap-merge across n BAM readers.
type mergeIterator struct {
	h    *recordHeap
	less func(a, b *sam.Record) bool
}

func newMergeIterator(readers []*sam.BAMReader, less func(a, b *sam.Record) bool) recordIter {
	h := &recordHeap{less: less}
	heap.Init(h)
	for i, rd := range readers {
		rec, err := rd.Read()
		if err == io.EOF {
			continue
		}
		if err != nil {
			// Surface the error on the next Next() call by stashing it.
			return &errIter{err: err}
		}
		heap.Push(h, mergeEntry{rec: rec, src: i, reader: rd})
	}
	return &mergeIterator{h: h, less: less}
}

func (m *mergeIterator) Next() (*sam.Record, error) {
	rec, _, err := m.nextWithSrc()
	return rec, err
}

// nextWithSrc behaves like Next but also reports the input index the record
// came from, which `merge` needs to apply its per-input @RG ID translation.
func (m *mergeIterator) nextWithSrc() (*sam.Record, int, error) {
	if m.h.Len() == 0 {
		return nil, -1, io.EOF
	}
	top := heap.Pop(m.h).(mergeEntry)
	// Pull the next record from that source.
	next, err := top.reader.Read()
	if err == nil {
		heap.Push(m.h, mergeEntry{rec: next, src: top.src, reader: top.reader})
	} else if err != io.EOF {
		return top.rec, top.src, err
	}
	return top.rec, top.src, nil
}

type errIter struct{ err error }

func (e *errIter) Next() (*sam.Record, error) { return nil, e.err }

type mergeEntry struct {
	rec    *sam.Record
	src    int
	reader *sam.BAMReader
}

type recordHeap struct {
	items []mergeEntry
	less  func(a, b *sam.Record) bool
}

func (h *recordHeap) Len() int { return len(h.items) }
func (h *recordHeap) Less(i, j int) bool {
	if h.less(h.items[i].rec, h.items[j].rec) {
		return true
	}
	if h.less(h.items[j].rec, h.items[i].rec) {
		return false
	}
	// Tie-break on source index so the merge is deterministic.
	return h.items[i].src < h.items[j].src
}
func (h *recordHeap) Swap(i, j int)      { h.items[i], h.items[j] = h.items[j], h.items[i] }
func (h *recordHeap) Push(x interface{}) { h.items = append(h.items, x.(mergeEntry)) }
func (h *recordHeap) Pop() interface{} {
	n := len(h.items)
	x := h.items[n-1]
	h.items = h.items[:n-1]
	return x
}

// recordSize estimates a record's ACTUAL Go heap footprint (not just its data
// bytes), used to decide when to spill an in-memory chunk. A *sam.Record is many
// separate heap allocations — the struct itself plus a backing array per
// string/slice field, each rounded up to a size class — so the real footprint
// is roughly twice the raw field bytes. Charging only the field bytes (the old
// estimate) let the in-memory chunk grow to ~2x MaxMemBytes before spilling,
// doubling peak RSS versus upstream's packed bam1_t budget.
func recordSize(r *sam.Record) int64 {
	// alloc approximates a single Go heap allocation of n payload bytes: the
	// allocator rounds small objects up to a size class (modelled as the next
	// multiple of 16) and there is per-object bookkeeping; an empty field still
	// costs nothing.
	alloc := func(n int) int64 {
		if n == 0 {
			return 0
		}
		return int64((n + 15) &^ 15)
	}
	const structBytes = 208 // the Record struct + its slot in the []*sam.Record buffer
	n := int64(structBytes) +
		alloc(len(r.QName)) + alloc(len(r.Seq)) + alloc(len(r.Qual)) +
		alloc(len(r.Cigar)*4) + alloc(len(r.RName)) + alloc(len(r.RNext))
	for _, a := range r.Aux {
		n += 48 // Aux struct + interface value box
		if s, ok := a.Value.(string); ok {
			n += alloc(len(s))
		}
		n += alloc(len(a.ArrayValues) * 8)
	}
	return n
}

// makeRecordLess constructs the comparison function that drives both the
// in-memory sort and the merge heap, per opts.
func makeRecordLess(opts SortOptions, refIndex map[string]int) func(a, b *sam.Record) bool {
	switch opts.Order {
	case SortByName:
		// Plain lexicographic QNAME (upstream `-N`, natural_sort=0), tie-broken
		// by the FLAG-derived key just like upstream's heap_lt/bam1_cmp_core.
		return func(a, b *sam.Record) bool { return nameLess(a, b, false) }
	case SortByNameNatural:
		// Natural-numeric QNAME (upstream `-n`, natural_sort=1), same FLAG
		// secondary key.
		return func(a, b *sam.Record) bool { return nameLess(a, b, true) }
	case SortByTag:
		tag := opts.Tag
		return func(a, b *sam.Record) bool { return tagLess(a, b, tag, refIndex) }
	default:
		return func(a, b *sam.Record) bool { return coordLess(a, b, refIndex) }
	}
}

// flagSortKey maps a record's FLAG to the integer upstream samtools uses as
// the secondary key for name sorts. It reorders the READ1/READ2, SECONDARY
// and SUPPLEMENTARY bits so that a plain numeric subtraction yields the
// order READ1, READ2, (PRIMARY), SUPPLEMENTARY, SECONDARY. This is a direct
// port of the bit shuffle in bam_sort.c's heap_lt and bam1_cmp_core:
//
//	f = ((flag&0xc0)<<8) | ((flag&0x100)<<3) | ((flag&0x800)>>3)
func flagSortKey(flag uint16) int {
	f := int(flag)
	return ((f & 0xc0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3)
}

// nameLess orders two records the way upstream samtools sort does for a name
// sort: first by QNAME (natural-numeric when natural is true, else a plain
// byte compare) and then, on a QNAME tie, by the FLAG-derived secondary key.
// This matches bam1_cmp_core's QueryName branch byte-for-byte.
func nameLess(a, b *sam.Record, natural bool) bool {
	c := strnumCmp(a.QName, b.QName, natural)
	if c != 0 {
		return c < 0
	}
	return flagSortKey(a.Flag) < flagSortKey(b.Flag)
}

// coordLess sorts records by the exact key upstream samtools uses for a
// coordinate sort (bam_sort.c bam1_cmp_core, non-QueryName branch):
// reference ID, then 1-based position, then the reverse-strand flag (forward
// before reverse). Unmapped records (tid == -1, mapped to UINT32_MAX as an
// unsigned compare) sort after every mapped record. Records identical in all
// three keys compare equal — their relative order is then fixed by the
// stable sort, reproducing htslib's "equal records keep input order" result
// (samtools sort spills/merges deterministically and our stable sort + heap
// tie-break on source index give the same record order on decode).
//
// NOTE: upstream tie-breaks on the reverse-strand flag, NOT on QNAME. An
// earlier port used a lexicographic QNAME tie-break here, which reordered
// equal-(rname,pos) records relative to htslib.
func coordLess(a, b *sam.Record, refIndex map[string]int) bool {
	return coordCmp(a, b, refIndex) < 0
}

// coordCmp is the three-way coordinate comparator (refID, 1-based pos,
// reverse-strand flag) shared with the stable-sort path. It returns -1, 0, or
// 1, mirroring bam1_cmp_core.
func coordCmp(a, b *sam.Record, refIndex map[string]int) int {
	ra := refID(a.RName, refIndex)
	rb := refID(b.RName, refIndex)
	if ra != rb {
		if ra < rb {
			return -1
		}
		return 1
	}
	if a.Pos != b.Pos {
		if a.Pos < b.Pos {
			return -1
		}
		return 1
	}
	arev := a.Flag&sam.FlagReverse != 0
	brev := b.Flag&sam.FlagReverse != 0
	if arev != brev {
		if !arev { // forward (0) sorts before reverse (1)
			return -1
		}
		return 1
	}
	return 0
}

// naturalLess implements numeric-aware string comparison: runs of digits in
// both strings are compared as integers, everything else byte-by-byte. So
// "r10" sorts after "r2" but before "r20". It is a thin wrapper over
// strnumCmp's natural mode, retained for callers that only need a bool.
func naturalLess(a, b string) bool { return strnumCmp(a, b, true) < 0 }

// strnumCmp is a byte-for-byte port of upstream samtools' strnum_cmp
// (bam_sort.c). It returns a negative, zero, or positive value if a sorts
// before, equal to, or after b respectively. When natural is false it
// degrades to a plain byte compare (the `samtools sort -N` path, where
// upstream sets natural_sort=0 and calls strcmp). When natural is true it
// compares maximal digit runs numerically: leading zeros are skipped, then
// matching digits are skipped, and whichever number keeps producing digits
// the longest is the larger; equal-length runs fall back to the first
// differing digit.
func strnumCmp(a, b string, natural bool) int {
	if !natural {
		if a < b {
			return -1
		}
		if a > b {
			return 1
		}
		return 0
	}
	pa, pb := 0, 0
	for pa < len(a) && pb < len(b) {
		if !isDigit(a[pa]) || !isDigit(b[pb]) {
			if a[pa] != b[pb] {
				return int(a[pa]) - int(b[pb])
			}
			pa++
			pb++
			continue
		}
		// Both positions are digits: compare the two numeric runs.
		// Skip leading zeros.
		for pa < len(a) && a[pa] == '0' {
			pa++
		}
		for pb < len(b) && b[pb] == '0' {
			pb++
		}
		// Skip matching digits.
		for pa < len(a) && pb < len(b) && isDigit(a[pa]) && a[pa] == b[pb] {
			pa++
			pb++
		}
		// Capture the byte difference at the current cursor, exactly as
		// upstream does ((int)*pa - (int)*pb). Past the end of either string
		// the C code reads the NUL terminator, i.e. 0.
		var ca, cb int
		if pa < len(a) {
			ca = int(a[pa])
		}
		if pb < len(b) {
			cb = int(b[pb])
		}
		diff := ca - cb
		// Whichever number is still emitting digits is the larger one.
		for pa < len(a) && pb < len(b) && isDigit(a[pa]) && isDigit(b[pb]) {
			pa++
			pb++
		}
		if pa < len(a) && isDigit(a[pa]) {
			return 1 // a still going, so larger
		}
		if pb < len(b) && isDigit(b[pb]) {
			return -1 // b still going, so larger
		}
		if diff != 0 {
			return diff // same length, so earlier diff decides
		}
	}
	if pa < len(a) {
		return 1
	}
	if pb < len(b) {
		return -1
	}
	return 0
}

func isDigit(b byte) bool { return b >= '0' && b <= '9' }

// tagLess reports whether record a sorts before record b under `samtools
// sort -t TAG` (the TagCoordinate order). It is a direct port of upstream's
// bam1_cmp_by_tag (bam_sort.c): records lacking the tag sort *first*; present
// tags are compared first by their normalised aux type, then by value
// (integers numerically, floats as floats — including the int-vs-float
// coercion upstream performs — single chars by code point, strings/hex
// byte-wise); on a full tie it falls back to the coordinate secondary key
// (refID, 1-based pos, reverse-strand flag) via tagCoordCmp.
func tagLess(a, b *sam.Record, tag string, refIndex map[string]int) bool {
	return tagCmp(a, b, tag, refIndex) < 0
}

// tagCmp returns a negative, zero, or positive value mirroring upstream's
// bam1_cmp_by_tag three-way result.
func tagCmp(a, b *sam.Record, tag string, refIndex map[string]int) int {
	av, aok := a.GetAux(tag)
	bv, bok := b.GetAux(tag)
	// Reads not carrying the tag sort first (return -1 means a < b).
	switch {
	case !aok && bok:
		return -1
	case aok && !bok:
		return 1
	case !aok && !bok:
		return tagCoordCmp(a, b, refIndex)
	}

	atype := normalizeAuxType(av.Type)
	btype := normalizeAuxType(bv.Type)
	if atype != btype {
		// Upstream coerces int↔float so the numeric ordering is total;
		// any other type mismatch falls back to comparing the type bytes.
		if atype == 'c' && btype == 'f' {
			atype = 'f'
		} else if atype == 'f' && btype == 'c' {
			btype = 'f'
		} else {
			if atype < btype {
				return -1
			}
			return 1
		}
	}

	switch atype {
	case 'c': // any integer width, normalised
		ai, _ := av.Int()
		bi, _ := bv.Int()
		if ai != bi {
			if ai < bi {
				return -1
			}
			return 1
		}
	case 'f':
		af := auxFloat(av)
		bf := auxFloat(bv)
		if af != bf {
			if af < bf {
				return -1
			}
			return 1
		}
	case 'A':
		as, _ := av.Value.(string)
		bs, _ := bv.Value.(string)
		ac, bc := byte(0), byte(0)
		if len(as) > 0 {
			ac = as[0]
		}
		if len(bs) > 0 {
			bc = bs[0]
		}
		if ac != bc {
			if ac < bc {
				return -1
			}
			return 1
		}
	case 'H': // Z and H both normalise to 'H'; compare byte-wise like strcmp
		as, _ := av.Value.(string)
		bs, _ := bv.Value.(string)
		if as != bs {
			if as < bs {
				return -1
			}
			return 1
		}
	}
	// Tags are equal — fall back to the coordinate secondary key.
	return tagCoordCmp(a, b, refIndex)
}

// normalizeAuxType folds a SAM aux type byte to the canonical letter upstream
// uses for total-ordering tag comparisons (bam_sort.c normalize_type): all
// integer widths become 'c', float/double become 'f', printable string and
// hex become 'H', and any other type keeps its own byte.
func normalizeAuxType(t byte) byte {
	switch t {
	case 'c', 'C', 's', 'S', 'i', 'I':
		return 'c'
	case 'f', 'd':
		return 'f'
	case 'H', 'Z':
		return 'H'
	default:
		return t
	}
}

// auxFloat returns an aux value as a float64, accepting either a stored
// float or an integer (matching upstream's bam_aux2f int-to-float read).
func auxFloat(a sam.Aux) float64 {
	if f, ok := a.Value.(float64); ok {
		return f
	}
	if i, ok := a.Int(); ok {
		return float64(i)
	}
	return 0
}

// tagCoordCmp is the secondary key used by `samtools sort -t TAG` when two
// records' tag values are equal: it orders by reference ID, then 1-based
// position, then the reverse-strand flag, exactly as upstream's
// bam1_cmp_core does for the (non-QueryName) coordinate path. Returns -1, 0,
// or 1.
func tagCoordCmp(a, b *sam.Record, refIndex map[string]int) int {
	return coordCmp(a, b, refIndex)
}

// refID resolves a reference name to its 0-based header index, mapping the
// unmapped sentinel ("*"/empty/unknown) to a value that sorts after every
// real reference — matching the unsigned tid ordering upstream relies on
// (tid == -1 becomes UINT32_MAX).
func refID(name string, refIndex map[string]int) int {
	if name != "" && name != "*" {
		if id, ok := refIndex[name]; ok {
			return id
		}
	}
	return 1 << 30
}

// writeShard creates a tmp BAM file and writes the given (already sorted)
// records to it. Returns the file path.
//
// FIX 3a: spill shards always use the SINGLE-THREADED BGZF writer, never the
// parallel MultiWriter back end. A -@N run spills many shards over its lifetime;
// previously each one spun up (and tore down) a fresh parallel-BGZF pool — N
// worker goroutines, a collector, and two 256-deep channels — for a small,
// serially written file. That create/destroy churn (plus the fd traffic) made
// the -@4 sys-time explode, so -@4 ran slower than -@1. The parallel BGZF writer
// is now reserved for the FINAL output only (see writeOutput), where the stream
// is large enough to amortise the pool. OUTPUT-COMPAT: only the spill shards'
// BGZF deflate threading changes; BGZF is block-parallel and decodes to
// identical plaintext regardless of how many goroutines compressed it, and the
// merge reads those identical records back, so the merged output is byte-for-
// byte unchanged.
func writeShard(dir, prefix string, idx int, hdr *sam.Header, recs []*sam.Record) (string, error) {
	path := filepath.Join(dir, fmt.Sprintf("%s.%d.%d.bam", prefix, os.Getpid(), idx))
	f, err := os.Create(path)
	if err != nil {
		return "", err
	}
	bw := sam.NewBAMWriter(f)
	if err := bw.WriteHeader(hdr); err != nil {
		_ = f.Close()
		_ = os.Remove(path)
		return "", err
	}
	for _, rec := range recs {
		if err := bw.Write(rec); err != nil {
			_ = bw.Close()
			_ = f.Close()
			_ = os.Remove(path)
			return "", err
		}
	}
	if err := bw.Close(); err != nil {
		_ = f.Close()
		_ = os.Remove(path)
		return "", err
	}
	if err := f.Close(); err != nil {
		_ = os.Remove(path)
		return "", err
	}
	return path, nil
}

// writePackedShard is writeShard for the packed path: it writes the spill file's
// header then each record's size prefix + raw body verbatim via WriteRaw, with
// no decode/re-encode. Like writeShard it always uses the single-threaded BGZF
// writer (the parallel pool is reserved for the final output). The spilled bytes
// decode to the same records the merge reads back, so the merged output is
// byte-for-byte unchanged from a re-encoded shard.
func writePackedShard(dir, prefix string, idx int, hdr *sam.Header, recs []packedRec) (string, error) {
	path := filepath.Join(dir, fmt.Sprintf("%s.%d.%d.bam", prefix, os.Getpid(), idx))
	f, err := os.Create(path)
	if err != nil {
		return "", err
	}
	bw := sam.NewBAMWriter(f)
	if err := bw.WriteHeader(hdr); err != nil {
		_ = f.Close()
		_ = os.Remove(path)
		return "", err
	}
	for i := range recs {
		if err := bw.WriteRaw(recs[i].raw); err != nil {
			_ = bw.Close()
			_ = f.Close()
			_ = os.Remove(path)
			return "", err
		}
	}
	if err := bw.Close(); err != nil {
		_ = f.Close()
		_ = os.Remove(path)
		return "", err
	}
	if err := f.Close(); err != nil {
		_ = os.Remove(path)
		return "", err
	}
	return path, nil
}

// writeOutput writes the iterator's records to out as BAM or SAM per opts,
// stamping the header's @HD SO field on the way.
func writeOutput(out io.Writer, hdr *sam.Header, it recordIter, opts SortOptions) error {
	// Stamp SO in @HD according to the chosen sort order.
	stampSortOrder(hdr, opts.Order)

	var w sam.Writer
	if opts.OutputSAM {
		w = sam.NewSAMWriter(out)
	} else {
		// Default: BAM. The parallel BGZF back end (-@/--threads) compresses
		// blocks concurrently; the decoded records are identical to the
		// single-threaded path.
		bw := sam.NewBAMWriterThreads(out, opts.Threads)
		w = bw
	}
	// closeOnErr ensures a parallel BGZF back end's worker goroutines drain on
	// any early return; Close is idempotent so the success path's explicit
	// Close is still the one whose error is reported.
	closed := false
	closeOnErr := func() {
		if !closed {
			_ = w.Close()
		}
	}
	if err := w.WriteHeader(hdr); err != nil {
		closeOnErr()
		return err
	}
	for {
		rec, err := it.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			closeOnErr()
			return err
		}
		if err := w.Write(rec); err != nil {
			closeOnErr()
			return err
		}
	}
	closed = true
	return w.Close()
}

// writePackedOutput writes a packed iterator's records to out, stamping @HD SO
// like writeOutput. For BAM output (the perf-critical, common case) it emits
// each raw record body verbatim via WriteRaw — no decode/re-encode — keeping the
// parallel BGZF back end (-@/--threads). For SAM text output (not perf-critical)
// it decodes each raw body to a *sam.Record and writes it through the SAM
// writer, exactly as the decode path's writeOutput does.
func writePackedOutput(out io.Writer, hdr *sam.Header, it packedIter, opts SortOptions) error {
	stampSortOrder(hdr, opts.Order)

	if opts.OutputSAM {
		w := sam.NewSAMWriter(out)
		if err := w.WriteHeader(hdr); err != nil {
			_ = w.Close()
			return err
		}
		// A header-only BAM reader used purely to decode buffered raw bodies.
		dec := sam.NewBAMBodyReader(nil, hdr)
		for {
			p, err := it.Next()
			if err == io.EOF {
				break
			}
			if err != nil {
				_ = w.Close()
				return err
			}
			rec, derr := dec.DecodeRecordBody(p.raw)
			if derr != nil {
				_ = w.Close()
				return derr
			}
			if err := w.Write(rec); err != nil {
				_ = w.Close()
				return err
			}
		}
		return w.Close()
	}

	// BAM output: raw passthrough through the (optionally parallel) BGZF writer.
	bw := sam.NewBAMWriterThreads(out, opts.Threads)
	closed := false
	closeOnErr := func() {
		if !closed {
			_ = bw.Close()
		}
	}
	if err := bw.WriteHeader(hdr); err != nil {
		closeOnErr()
		return err
	}
	for {
		p, err := it.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			closeOnErr()
			return err
		}
		if err := bw.WriteRaw(p.raw); err != nil {
			closeOnErr()
			return err
		}
	}
	closed = true
	return bw.Close()
}

// stampSortOrder rewrites the @HD line's SO: and SS: fields (creating an
// @HD line if none exists) to reflect the sort order in use. The SS
// (sub-sort) tag is what upstream samtools writes for queryname order —
// "queryname:natural" for the default name-sort, or
// "queryname:lexicographical" for the -N variant.
func stampSortOrder(hdr *sam.Header, order SortOrder) {
	var so, ss string
	switch order {
	case SortByName:
		so = "queryname"
		ss = "queryname:lexicographical"
	case SortByNameNatural:
		so = "queryname"
		ss = "queryname:natural"
	case SortByTag:
		so = "unknown" // SAM spec does not define a tag-sort SO value.
	default:
		so = "coordinate"
	}
	setField := func(line *sam.HeaderLine, tag, value string) {
		for i := range line.Fields {
			if line.Fields[i].Tag == tag {
				line.Fields[i].Value = value
				return
			}
		}
		line.Fields = append(line.Fields, sam.HeaderField{Tag: tag, Value: value})
	}
	removeField := func(line *sam.HeaderLine, tag string) {
		out := line.Fields[:0]
		for _, f := range line.Fields {
			if f.Tag != tag {
				out = append(out, f)
			}
		}
		line.Fields = out
	}
	// Find an existing @HD line.
	for i := range hdr.Lines {
		if hdr.Lines[i].Tag != "HD" {
			continue
		}
		setField(&hdr.Lines[i], "SO", so)
		if ss != "" {
			setField(&hdr.Lines[i], "SS", ss)
		} else {
			removeField(&hdr.Lines[i], "SS")
		}
		hdr.HDFields = hdr.Lines[i].Fields
		return
	}
	// No @HD line — create one and place it at the top.
	fields := []sam.HeaderField{
		{Tag: "VN", Value: "1.6"},
		{Tag: "SO", Value: so},
	}
	if ss != "" {
		fields = append(fields, sam.HeaderField{Tag: "SS", Value: ss})
	}
	newLine := sam.HeaderLine{Tag: "HD", Fields: fields}
	hdr.Lines = append([]sam.HeaderLine{newLine}, hdr.Lines...)
	hdr.HDFields = newLine.Fields
}

// ParseSortOrder recognises the upstream samtools `-O` order strings.
func ParseSortOrder(s string) (SortOrder, error) {
	switch strings.ToLower(s) {
	case "", "coordinate":
		return SortCoordinate, nil
	case "queryname", "name":
		return SortByName, nil
	case "natural":
		return SortByNameNatural, nil
	case "tag":
		return SortByTag, nil
	}
	return 0, fmt.Errorf("samtools sort: unknown sort order %q", s)
}

// ParseMemBudget parses a memory budget string in the form N[K|M|G]. Plain
// integers are treated as bytes.
func ParseMemBudget(s string) (int64, error) {
	if s == "" {
		return 0, nil
	}
	mult := int64(1)
	body := s
	switch s[len(s)-1] {
	case 'k', 'K':
		mult = 1 << 10
		body = s[:len(s)-1]
	case 'm', 'M':
		mult = 1 << 20
		body = s[:len(s)-1]
	case 'g', 'G':
		mult = 1 << 30
		body = s[:len(s)-1]
	}
	n, err := strconv.ParseInt(body, 10, 64)
	if err != nil {
		return 0, fmt.Errorf("samtools sort: bad memory budget %q: %w", s, err)
	}
	if n < 0 {
		return 0, fmt.Errorf("samtools sort: negative memory budget %q", s)
	}
	return n * mult, nil
}
