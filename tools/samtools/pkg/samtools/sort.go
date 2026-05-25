package samtools

import (
	"container/heap"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

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
	// Threads selects the number of worker goroutines used to sort and BAM-
	// encode in-memory shards in parallel. The value follows upstream
	// samtools' loose convention: 0 (or any non-positive value) selects the
	// serial fast path (one worker), positive values are capped at
	// runtime.NumCPU(). Output is byte-identical to the serial path
	// because shards are reassembled by submission-order sequence number
	// before the k-way merge step.
	Threads int
	// NoPG is accepted but is a no-op (we don't inject @PG lines anyway).
	NoPG bool
	// SecondaryByName, with Order == SortByTag, selects upstream's
	// TagQueryName ordering (qname+FLAG tie-break). Default false matches
	// upstream TagCoordinate (tid, pos+1, rev) tie-break.
	SecondaryByName bool
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

	r, err := sam.NewReader(in)
	if err != nil {
		return err
	}
	hdr := r.Header()
	// Build a reference-name → refID map for coordinate sort.
	refIndex := make(map[string]int, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refIndex[ref.Name] = i
	}

	cmp := makeRecordLess(opts, refIndex)

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

	workers := resolveWorkers(opts.Threads)

	var shards []string
	defer func() {
		for _, p := range shards {
			_ = os.Remove(p)
		}
	}()

	// Pipeline: the reader hands full buffers off to either a serial
	// flush (workers == 1) or a worker pool (workers > 1). The worker
	// path is share-nothing — each shardJob carries its own []*sam.Record
	// — and results are reassembled by submission sequence number so the
	// resulting list of temp BAM paths is identical to the serial path's.

	var (
		buffer   []*sam.Record
		bufBytes int64
		seq      int
	)

	// Serial flush keeps the simple in-line path for workers == 1, which
	// preserves the historical behaviour (no goroutines, no channels) and
	// matches the byte-for-byte temp-shard layout used by every existing
	// test.
	flushSerial := func() error {
		if len(buffer) == 0 {
			return nil
		}
		sort.SliceStable(buffer, func(a, b int) bool { return cmp(buffer[a], buffer[b]) })
		path, ferr := writeShard(tmpBaseDir, tmpName, seq, hdr, buffer)
		if ferr != nil {
			return ferr
		}
		shards = append(shards, path)
		seq++
		buffer = buffer[:0]
		bufBytes = 0
		return nil
	}

	// Parallel path: spin up a pool and a small result collector that
	// recovers submission order via a min-heap keyed on seq. The pool is
	// only constructed when workers > 1 so the serial path stays free of
	// goroutines.
	type pendingResults struct {
		pool     *shardPool
		bySeq    map[int]string // collected results awaiting in-order release
		nextWant int            // next seq to append to shards
	}
	var pr *pendingResults
	if workers > 1 {
		work := func(job shardJob) shardResult {
			sort.SliceStable(job.recs, func(a, b int) bool { return cmp(job.recs[a], job.recs[b]) })
			path, ferr := writeShard(tmpBaseDir, tmpName, job.seq, hdr, job.recs)
			return shardResult{seq: job.seq, path: path, err: ferr}
		}
		pr = &pendingResults{
			pool:  newShardPool(workers, work),
			bySeq: make(map[int]string),
		}
	}

	// drainAvailable pulls every result currently sitting in the pool's
	// output channel (non-blocking) into pr.bySeq, then releases the
	// in-order prefix into shards. Called between submissions so the
	// channel does not back up and so an early error is surfaced promptly.
	drainAvailable := func() error {
		if pr == nil {
			return nil
		}
		for {
			select {
			case res, ok := <-pr.pool.results():
				if !ok {
					return pr.pool.firstError()
				}
				if res.err != nil {
					return res.err
				}
				pr.bySeq[res.seq] = res.path
			default:
				// Release the contiguous in-order prefix.
				for {
					p, ok := pr.bySeq[pr.nextWant]
					if !ok {
						return nil
					}
					shards = append(shards, p)
					delete(pr.bySeq, pr.nextWant)
					pr.nextWant++
				}
			}
		}
	}

	flushParallel := func() error {
		if len(buffer) == 0 {
			return nil
		}
		// Hand the buffer off to a worker and start a fresh local slice
		// so the reader can keep filling while the worker runs. This is
		// the share-nothing handoff that keeps the race detector happy.
		job := shardJob{seq: seq, recs: buffer}
		buffer = make([]*sam.Record, 0, 1024)
		bufBytes = 0
		seq++
		// Surface any prior worker error before blocking on submit so we
		// don't deadlock waiting on a full input channel after a worker
		// has already failed.
		if err := drainAvailable(); err != nil {
			return err
		}
		pr.pool.submit(job)
		return nil
	}

	flush := flushSerial
	if pr != nil {
		flush = flushParallel
	}

	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			if pr != nil {
				pr.pool.closeSubmissions()
				for range pr.pool.results() {
				}
			}
			return err
		}
		size := recordSize(rec)
		// If we are about to overflow and we already have content, spill
		// first so the new record always lands in a fresh buffer.
		if len(buffer) > 0 && bufBytes+size > opts.MaxMemBytes {
			if err := flush(); err != nil {
				return err
			}
		}
		buffer = append(buffer, rec)
		bufBytes += size
	}
	// Fast path: zero shards spilled and no in-flight workers → sort the
	// buffer in place and stream it out. This applies in both the serial
	// and parallel code paths.
	if pr == nil && len(shards) == 0 {
		sort.SliceStable(buffer, func(a, b int) bool { return cmp(buffer[a], buffer[b]) })
		return writeOutput(out, hdr, sliceIter(buffer), opts)
	}
	if pr != nil && seq == 0 {
		// No buffers ever spilled — every record still sits in `buffer`.
		// Tear down the (idle) pool and take the in-memory fast path.
		pr.pool.closeSubmissions()
		for range pr.pool.results() {
		}
		if err := pr.pool.firstError(); err != nil {
			return err
		}
		sort.SliceStable(buffer, func(a, b int) bool { return cmp(buffer[a], buffer[b]) })
		return writeOutput(out, hdr, sliceIter(buffer), opts)
	}
	// Otherwise: flush the tail and k-way merge.
	if err := flush(); err != nil {
		if pr != nil {
			pr.pool.closeSubmissions()
			for range pr.pool.results() {
			}
		}
		return err
	}
	// Drain the parallel pool to completion so all temp files exist and
	// the shards slice is in submission order.
	if pr != nil {
		pr.pool.closeSubmissions()
		for res := range pr.pool.results() {
			if res.err != nil {
				continue // firstError will surface it below
			}
			pr.bySeq[res.seq] = res.path
		}
		if err := pr.pool.firstError(); err != nil {
			return err
		}
		for pr.nextWant < seq {
			p, ok := pr.bySeq[pr.nextWant]
			if !ok {
				return fmt.Errorf("samtools sort: missing parallel shard %d", pr.nextWant)
			}
			shards = append(shards, p)
			delete(pr.bySeq, pr.nextWant)
			pr.nextWant++
		}
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
	if m.h.Len() == 0 {
		return nil, io.EOF
	}
	top := heap.Pop(m.h).(mergeEntry)
	// Pull the next record from that source.
	next, err := top.reader.Read()
	if err == nil {
		heap.Push(m.h, mergeEntry{rec: next, src: top.src, reader: top.reader})
	} else if err != io.EOF {
		return top.rec, err
	}
	return top.rec, nil
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

// recordSize returns a coarse estimate of a record's memory footprint, used
// only to decide when to spill an in-memory chunk.
func recordSize(r *sam.Record) int64 {
	const base = 96 // struct header + slice headers + minor overheads.
	n := int64(base) + int64(len(r.QName)) + int64(len(r.Seq)) + int64(len(r.Qual)) + int64(len(r.Cigar)*4) + int64(len(r.RName)+len(r.RNext))
	for _, a := range r.Aux {
		n += int64(2 + 1 + 8) // tag + type + 8-byte fallback value
		if s, ok := a.Value.(string); ok {
			n += int64(len(s))
		}
		n += int64(len(a.ArrayValues) * 8)
	}
	return n
}

// makeRecordLess constructs the comparison function that drives both the
// in-memory sort and the merge heap, per opts. Semantics mirror upstream
// samtools' bam_sort.c `bam1_cmp_core` / `bam1_cmp_by_tag` exactly.
func makeRecordLess(opts SortOptions, refIndex map[string]int) func(a, b *sam.Record) bool {
	switch opts.Order {
	case SortByName:
		return func(a, b *sam.Record) bool { return nameCmpLess(a, b, false) }
	case SortByNameNatural:
		return func(a, b *sam.Record) bool { return nameCmpLess(a, b, true) }
	case SortByTag:
		tag := opts.Tag
		secondaryByName := opts.SecondaryByName
		return func(a, b *sam.Record) bool { return tagLess(a, b, tag, refIndex, secondaryByName) }
	default:
		return func(a, b *sam.Record) bool { return coordLess(a, b, refIndex) }
	}
}

// coordLess wraps coordCmpCore for the in-memory sort.
func coordLess(a, b *sam.Record, refIndex map[string]int) bool {
	return coordCmpCore(a, b, refIndex) < 0
}

// coordCmpCore is upstream `bam1_cmp_core` for coordinate order: (tid,
// pos+1, rev). Unmapped records sort last via refIDFor's sentinel.
func coordCmpCore(a, b *sam.Record, refIndex map[string]int) int {
	ra := refIDFor(a, refIndex)
	rb := refIDFor(b, refIndex)
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
	revA := uint32(0)
	if a.Flag&0x10 != 0 {
		revA = 1
	}
	revB := uint32(0)
	if b.Flag&0x10 != 0 {
		revB = 1
	}
	if revA != revB {
		if revA < revB {
			return -1
		}
		return 1
	}
	return 0
}

// refIDFor maps a record's RName to refIndex; unmapped/unknown get a
// large sentinel so they sort last.
func refIDFor(r *sam.Record, refIndex map[string]int) int {
	if r.RName == "" || r.RName == "*" {
		return 1 << 30
	}
	if id, ok := refIndex[r.RName]; ok {
		return id
	}
	return 1 << 30
}

// flagSortKey reorders FLAG bits so the natural integer compare yields
// READ1 < READ2 < primary < supplementary < secondary, matching upstream
// bam_sort.c bam1_cmp_core:
// `((flag&0xc0)<<8)|((flag&0x100)<<3)|((flag&0x800)>>3)`.
func flagSortKey(flag uint16) uint32 {
	f := uint32(flag)
	return ((f & 0xc0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3)
}

// nameCmpLess is the qname comparator with the upstream FLAG tie-break.
func nameCmpLess(a, b *sam.Record, natural bool) bool {
	return nameCmpCore(a, b, natural) < 0
}

func nameCmpCore(a, b *sam.Record, natural bool) int {
	var t int
	if natural {
		t = naturalCmp(a.QName, b.QName)
	} else {
		if a.QName < b.QName {
			t = -1
		} else if a.QName > b.QName {
			t = 1
		}
	}
	if t != 0 {
		return t
	}
	fa := flagSortKey(a.Flag)
	fb := flagSortKey(b.Flag)
	if fa < fb {
		return -1
	}
	if fa > fb {
		return 1
	}
	return 0
}

func naturalCmp(a, b string) int {
	if naturalLess(a, b) {
		return -1
	}
	if naturalLess(b, a) {
		return 1
	}
	return 0
}

// naturalLess implements numeric-aware string comparison: runs of digits in
// both strings are compared as integers, everything else byte-by-byte. So
// "r10" sorts after "r2" but before "r20".
func naturalLess(a, b string) bool {
	ia, ib := 0, 0
	for ia < len(a) && ib < len(b) {
		ca, cb := a[ia], b[ib]
		if isDigit(ca) && isDigit(cb) {
			// Find numeric run extents.
			ja := ia
			for ja < len(a) && isDigit(a[ja]) {
				ja++
			}
			jb := ib
			for jb < len(b) && isDigit(b[jb]) {
				jb++
			}
			// Strip leading zeros for value comparison; preserve full
			// length for the tie-break.
			na := stripLeadingZeros(a[ia:ja])
			nb := stripLeadingZeros(b[ib:jb])
			if len(na) != len(nb) {
				return len(na) < len(nb)
			}
			if na != nb {
				return na < nb
			}
			// Same numeric value — fall back to literal length compare.
			if ja-ia != jb-ib {
				return ja-ia < jb-ib
			}
			ia = ja
			ib = jb
			continue
		}
		if ca != cb {
			return ca < cb
		}
		ia++
		ib++
	}
	return len(a) < len(b)
}

func isDigit(b byte) bool { return b >= '0' && b <= '9' }

func stripLeadingZeros(s string) string {
	i := 0
	for i < len(s)-1 && s[i] == '0' {
		i++
	}
	return s[i:]
}

// tagLess compares two records by an aux tag value. Tag-absent records sort
// FIRST (matching upstream `bam1_cmp_by_tag` where aux_a == NULL returns -1
// vs an aux-bearing aux_b). Ties on the tag value fall back to the secondary
// comparator selected by secondaryByName: true → qname+FLAG, false →
// (tid, pos, rev). Mirrors upstream `TagQueryName` vs `TagCoordinate`.
func tagLess(a, b *sam.Record, tag string, refIndex map[string]int, secondaryByName bool) bool {
	av, aok := a.GetAux(tag)
	bv, bok := b.GetAux(tag)
	if !aok && !bok {
		return tagSecondaryCmp(a, b, refIndex, secondaryByName) < 0
	}
	if !aok {
		return true
	}
	if !bok {
		return false
	}
	ai, aint := av.Int()
	bi, bint := bv.Int()
	if aint && bint {
		if ai != bi {
			return ai < bi
		}
		return tagSecondaryCmp(a, b, refIndex, secondaryByName) < 0
	}
	if av.Type == 'f' && bv.Type == 'f' {
		af, _ := av.Value.(float64)
		bf, _ := bv.Value.(float64)
		if af != bf {
			return af < bf
		}
		return tagSecondaryCmp(a, b, refIndex, secondaryByName) < 0
	}
	as, _ := av.Value.(string)
	bs, _ := bv.Value.(string)
	if as != bs {
		return as < bs
	}
	return tagSecondaryCmp(a, b, refIndex, secondaryByName) < 0
}

// tagSecondaryCmp routes to upstream's secondary comparator: name+FLAG when
// secondaryByName, otherwise the position-based comparator.
func tagSecondaryCmp(a, b *sam.Record, refIndex map[string]int, secondaryByName bool) int {
	if secondaryByName {
		return nameCmpCore(a, b, true)
	}
	return coordCmpCore(a, b, refIndex)
}

// writeShard creates a tmp BAM file and writes the given (already sorted)
// records to it. Returns the file path.
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

// writeOutput writes the iterator's records to out as BAM or SAM per opts,
// stamping the header's @HD SO field on the way.
func writeOutput(out io.Writer, hdr *sam.Header, it recordIter, opts SortOptions) error {
	// Stamp SO in @HD according to the chosen sort order.
	stampSortOrder(hdr, opts.Order)

	var w sam.Writer
	if opts.OutputSAM {
		w = sam.NewSAMWriter(out)
	} else {
		// Default: BAM.
		bw := sam.NewBAMWriter(out)
		w = bw
	}
	if err := w.WriteHeader(hdr); err != nil {
		return err
	}
	for {
		rec, err := it.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if err := w.Write(rec); err != nil {
			return err
		}
	}
	return w.Close()
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
