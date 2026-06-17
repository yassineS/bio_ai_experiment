package bedmerge

import "container/heap"

// strandQueue ports upstream's StrandQueue (src/utils/FileRecordTools/Records/
// StrandQueue.{h,cpp}): a set of per-strand priority queues used by the merge
// state machine to defer records that cannot join the current merged group yet
// (wrong strand under -s, a different chromosome, or out of range). Each queue
// is a min-heap keyed on (chrom, start, end) — upstream's Record::operator< —
// so a record pulled back out is always the smallest-positioned one available.
//
// Upstream keeps three queues (FORWARD, REVERSE, UNKNOWN). top()/pop() return
// the global minimum across all queues; top(strand)/pop(strand) restrict to a
// single strand. We index queues by strand string ("+", "-", and "" for the
// unknown/'.' bucket).
type strandQueue struct {
	queues map[string]*recHeap
}

// strandKey maps a record's strand value to its queue key. '+' and '-' map to
// themselves; everything else ('.' or empty) is the unknown bucket "".
func strandKey(strand string) string {
	if strand == "+" || strand == "-" {
		return strand
	}
	return ""
}

func (q *strandQueue) ensure(key string) *recHeap {
	if q.queues == nil {
		q.queues = make(map[string]*recHeap, 3)
	}
	h := q.queues[key]
	if h == nil {
		h = &recHeap{}
		q.queues[key] = h
	}
	return h
}

// push stores a record in the queue for its strand.
func (q *strandQueue) push(r record) {
	heap.Push(q.ensure(strandKey(r.strand)), r)
}

// top returns the globally smallest-positioned stored record across all
// strands, without removing it.
func (q *strandQueue) top() (record, bool) {
	key, ok := q.minKey()
	if !ok {
		return record{}, false
	}
	return (*q.queues[key])[0], true
}

// pop removes the globally smallest-positioned stored record.
func (q *strandQueue) pop() {
	key, ok := q.minKey()
	if !ok {
		return
	}
	heap.Pop(q.queues[key])
}

// topStrand returns the smallest-positioned stored record for the given strand.
func (q *strandQueue) topStrand(strand string) (record, bool) {
	h := q.queues[strandKey(strand)]
	if h == nil || h.Len() == 0 {
		return record{}, false
	}
	return (*h)[0], true
}

// popStrand removes the smallest-positioned stored record for the given strand.
func (q *strandQueue) popStrand(strand string) {
	h := q.queues[strandKey(strand)]
	if h == nil || h.Len() == 0 {
		return
	}
	heap.Pop(h)
}

// minKey returns the queue key holding the global minimum record. Ties between
// queues are broken by record order via recLess, matching upstream getMinIdx
// (which keeps the first queue whose top is strictly less than the running min).
func (q *strandQueue) minKey() (string, bool) {
	var bestKey string
	var best record
	found := false
	// Iterate in upstream's queue order: FORWARD, REVERSE, UNKNOWN.
	for _, key := range []string{"+", "-", ""} {
		h := q.queues[key]
		if h == nil || h.Len() == 0 {
			continue
		}
		top := (*h)[0]
		if !found || recLess(top, best) {
			best = top
			bestKey = key
			found = true
		}
	}
	return bestKey, found
}

// recLess is upstream Record::operator<: order by (chrom, start, end).
func recLess(a, b record) bool {
	if a.chrom != b.chrom {
		return a.chrom < b.chrom
	}
	if a.start != b.start {
		return a.start < b.start
	}
	return a.end < b.end
}

// recHeap is a min-heap of records ordered by recLess.
type recHeap []record

func (h recHeap) Len() int            { return len(h) }
func (h recHeap) Less(i, j int) bool  { return recLess(h[i], h[j]) }
func (h recHeap) Swap(i, j int)       { h[i], h[j] = h[j], h[i] }
func (h *recHeap) Push(x interface{}) { *h = append(*h, x.(record)) }
func (h *recHeap) Pop() interface{} {
	old := *h
	n := len(old)
	r := old[n-1]
	*h = old[:n-1]
	return r
}
