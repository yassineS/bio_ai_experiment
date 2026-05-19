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
	// Threads is accepted but ignored — the v1 pipeline is single-threaded.
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
		if len(buffer) > 0 && bufBytes+size > opts.MaxMemBytes {
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
// in-memory sort and the merge heap, per opts.
func makeRecordLess(opts SortOptions, refIndex map[string]int) func(a, b *sam.Record) bool {
	switch opts.Order {
	case SortByName:
		return func(a, b *sam.Record) bool { return a.QName < b.QName }
	case SortByNameNatural:
		return func(a, b *sam.Record) bool { return naturalLess(a.QName, b.QName) }
	case SortByTag:
		tag := opts.Tag
		return func(a, b *sam.Record) bool { return tagLess(a, b, tag) }
	default:
		return func(a, b *sam.Record) bool { return coordLess(a, b, refIndex) }
	}
}

// coordLess sorts records by (refID, 0-based pos). Unmapped records (refID
// == -1) sort after every mapped record. Within a reference, ties on pos
// fall back to lexicographic QName so the order is fully deterministic.
func coordLess(a, b *sam.Record, refIndex map[string]int) bool {
	ra := -1
	if a.RName != "" && a.RName != "*" {
		if id, ok := refIndex[a.RName]; ok {
			ra = id
		}
	}
	rb := -1
	if b.RName != "" && b.RName != "*" {
		if id, ok := refIndex[b.RName]; ok {
			rb = id
		}
	}
	// Convert -1 to a sentinel that sorts last.
	const last = 1 << 30
	if ra < 0 {
		ra = last
	}
	if rb < 0 {
		rb = last
	}
	if ra != rb {
		return ra < rb
	}
	if a.Pos != b.Pos {
		return a.Pos < b.Pos
	}
	return a.QName < b.QName
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

// tagLess compares two records by an aux tag value. Missing tags sort last;
// integer tags compare numerically; float tags compare as floats; string
// tags compare lexicographically.
func tagLess(a, b *sam.Record, tag string) bool {
	av, aok := a.GetAux(tag)
	bv, bok := b.GetAux(tag)
	if !aok && !bok {
		return a.QName < b.QName
	}
	if !aok {
		return false
	}
	if !bok {
		return true
	}
	// Integer comparison (works for c/C/s/S/i/I).
	ai, aint := av.Int()
	bi, bint := bv.Int()
	if aint && bint {
		if ai != bi {
			return ai < bi
		}
		return a.QName < b.QName
	}
	if av.Type == 'f' && bv.Type == 'f' {
		af, _ := av.Value.(float64)
		bf, _ := bv.Value.(float64)
		if af != bf {
			return af < bf
		}
		return a.QName < b.QName
	}
	as, _ := av.Value.(string)
	bs, _ := bv.Value.(string)
	if as != bs {
		return as < bs
	}
	return a.QName < b.QName
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
