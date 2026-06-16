package cram

import (
	"io"
	"runtime"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// SetThreads enables block/slice-parallel CRAM decode for this RecordReader,
// decoding up to threads data containers concurrently while still yielding
// records in strict file order. A threads value below 2 leaves the reader on
// its single-threaded path (the default); any higher value spreads the
// expensive per-container codec decode and record reconstruction across a
// worker pool.
//
// The non-negotiable invariant: the sequence of records Read returns — and
// therefore every downstream byte of text or binary output — is identical for
// any thread count. Threading changes only decode throughput, never the data.
// This mirrors the BGZF MultiReader used for BAM input: a feeder reads the
// cheap structural container framing sequentially, workers run the costly
// decode out of order, and an ordering collector reassembles the records in
// the original container order.
//
// SetThreads must be called before the first Read. Calling it after decoding
// has begun is a no-op.
func (rr *RecordReader) SetThreads(threads int) {
	if rr.next > 0 || rr.done || rr.par != nil {
		return
	}
	if threads < 2 {
		rr.threads = 0
		return
	}
	rr.threads = threads
}

// parallelDriver decodes CRAM data containers concurrently and delivers their
// reconstructed records to the RecordReader in strict file order. It owns the
// feeder, the worker pool and the ordering collector; the RecordReader pulls
// ordered per-container record batches from out.
//
// The structural Reader (rd) is advanced only by the single feeder goroutine,
// so its non-concurrent Next contract is upheld. The expensive work —
// ParseDataContainer, block decompression, data-series decode and record
// reconstruction, all of which operate on a self-contained Container — runs in
// the workers. Reference resolution is funnelled through the RecordReader's
// mutex-guarded helper so the shared resolver memo and reference handles are
// touched by one goroutine at a time.
type parallelDriver struct {
	rr      *RecordReader
	threads int

	jobs    chan *containerJob
	results chan *containerJob
	out     chan *containerJob

	wg     sync.WaitGroup
	closed chan struct{}
	once   sync.Once
}

// containerJob carries one structural data container through the worker
// pipeline. seq fixes its file position so the collector can emit batches in
// order even though workers finish out of order. recs holds the container's
// reconstructed records once decoded; eof marks the end of the stream.
type containerJob struct {
	seq            int64
	container      *Container
	recs           []*sam.Record
	needsReference bool
	eof            bool
	err            error
}

// containerJobBuffer bounds the number of in-flight containers (queued plus
// awaiting collection), keeping memory modest while leaving slack to keep
// workers busy. A CRAM container is far larger than a BGZF block, so a smaller
// bound than the BGZF reader's suffices.
const containerJobBuffer = 32

// startParallel spins up the feeder, worker pool and ordering collector. It is
// called lazily on the first Read once SetThreads has requested parallel
// decode.
func (rr *RecordReader) startParallel() {
	threads := rr.threads
	if max := runtime.NumCPU(); threads > max+1 {
		// Cap workers at the available parallelism plus a small slack; more
		// goroutines than cores only adds scheduling overhead for a
		// CPU-bound decode.
		threads = max + 1
	}
	pd := &parallelDriver{
		rr:      rr,
		threads: threads,
		jobs:    make(chan *containerJob, containerJobBuffer),
		results: make(chan *containerJob, containerJobBuffer),
		out:     make(chan *containerJob, containerJobBuffer),
		closed:  make(chan struct{}),
	}
	rr.par = pd
	for i := 0; i < threads; i++ {
		pd.wg.Add(1)
		go pd.worker()
	}
	go func() {
		pd.wg.Wait()
		close(pd.results)
	}()
	go pd.feed()
	go pd.collect()
}

// feed reads structural data containers sequentially and dispatches each to the
// worker pool. Reading a container is cheap framing work (no codec decode); the
// expensive decode runs in the workers. Non-data containers (none normally
// follow the file-header container, but the format permits them) are skipped
// here so the workers only ever see decodable containers.
func (pd *parallelDriver) feed() {
	defer close(pd.jobs)
	var seq int64
	for {
		c, err := pd.rr.rd.Next()
		if err == io.EOF {
			pd.send(&containerJob{seq: seq, eof: true})
			return
		}
		if err != nil {
			pd.send(&containerJob{seq: seq, err: err})
			return
		}
		if len(c.Blocks) == 0 || c.Blocks[0].ContentType != ContentCompressionHeader {
			continue // a non-data container; do not advance the sequence.
		}
		if !pd.send(&containerJob{seq: seq, container: c}) {
			return
		}
		seq++
	}
}

// send queues a job for the workers, returning false if the driver has been
// closed (the consumer abandoned the stream) so the feeder can stop early.
func (pd *parallelDriver) send(job *containerJob) bool {
	select {
	case pd.jobs <- job:
		return true
	case <-pd.closed:
		return false
	}
}

// worker decodes one container's records and forwards the batch to the
// collector tagged with its sequence number. The decode is identical to the
// sequential decodeContainerInto path; only the reference resolution is routed
// through the mutex-guarded helper so the shared resolver is safe.
func (pd *parallelDriver) worker() {
	defer pd.wg.Done()
	for job := range pd.jobs {
		if job.err != nil || job.eof {
			pd.emit(job)
			continue
		}
		recs, needsRef, err := pd.rr.decodeContainerParallel(job.container)
		job.recs = recs
		job.needsReference = needsRef
		job.err = err
		pd.emit(job)
	}
}

// emit forwards a decoded job to the collector, unblocking on close so a worker
// never stalls on a full results channel after the consumer abandoned the
// stream.
func (pd *parallelDriver) emit(job *containerJob) {
	select {
	case pd.results <- job:
	case <-pd.closed:
	}
}

// collect receives decoded containers out of order, buffers them by sequence
// number, and forwards their record batches to the consumer strictly in
// sequence. The first error (or the EOF marker) terminates the ordered stream.
func (pd *parallelDriver) collect() {
	defer close(pd.out)
	var want int64
	pending := make(map[int64]*containerJob)
	for job := range pd.results {
		pending[job.seq] = job
		for {
			next, ok := pending[want]
			if !ok {
				break
			}
			delete(pending, want)
			select {
			case pd.out <- next:
			case <-pd.closed:
				return
			}
			if next.err != nil || next.eof {
				return
			}
			want++
		}
	}
}

// fillNextSliceParallel pulls the next ordered container batch from the driver
// into the RecordReader's pending buffer. It is the parallel analogue of
// fillNextSlice: each call yields one container's worth of records (in file
// order), or sets done at end of stream. Containers that decoded to zero
// records are skipped so Read always makes progress.
func (rr *RecordReader) fillNextSliceParallel() error {
	if rr.par == nil {
		rr.startParallel()
	}
	rr.pending = rr.pending[:0]
	rr.next = 0
	for {
		job, ok := <-rr.par.out
		if !ok {
			rr.done = true
			return nil
		}
		if job.err != nil {
			rr.done = true
			rr.par.stop()
			return job.err
		}
		if job.eof {
			rr.done = true
			return nil
		}
		if job.needsReference {
			rr.needsReference = true
		}
		if len(job.recs) == 0 {
			continue // an empty container; keep pulling.
		}
		rr.pending = job.recs
		return nil
	}
}

// stop tears down the driver, unblocking the feeder, workers and collector if
// the consumer abandons the stream before EOF, and waits for the worker pool to
// drain. Waiting matters because Close may release the reference FASTA handle
// right after; a worker must not still be fetching reference bases through it.
// It is idempotent and safe to call from the collector's own goroutine path
// (the wait only blocks on the workers, never on the collector).
func (pd *parallelDriver) stop() {
	pd.once.Do(func() {
		close(pd.closed)
		pd.wg.Wait()
	})
}

// decodeContainerParallel decodes one structural data container into its
// reconstructed records in file order, returning the records, whether any
// record needed an external reference, and the first error. It mirrors
// decodeContainerInto + decodeSlice but takes no shared RecordReader mutable
// state: needsReference is returned rather than written to the shared field,
// and reference resolution is funnelled through the mutex-guarded helper. This
// keeps the per-container decode safe to run on a worker goroutine.
func (rr *RecordReader) decodeContainerParallel(c *Container) ([]*sam.Record, bool, error) {
	dc, err := ParseDataContainer(c)
	if err != nil {
		return nil, false, wrapf(err, "container %d", c.Index)
	}
	var out []*sam.Record
	needsRef := false
	for si, sl := range dc.Slices {
		recs, sliceNeedsRef, err := rr.decodeSliceParallel(dc.Compression, sl, c.Index, si)
		if err != nil {
			return nil, false, err
		}
		if sliceNeedsRef {
			needsRef = true
		}
		out = append(out, recs...)
	}
	return out, needsRef, nil
}

// decodeSliceParallel decodes one slice's records on a worker goroutine. It is
// decodeSlice with two changes for concurrency safety: it returns whether the
// slice needed an external reference (rather than writing the shared
// needsReference field), and it resolves the slice's reference span through the
// mutex-guarded resolveSliceReferenceLocked. The per-slice decode itself
// (newRecordDecoder, decodeSliceRecords, regenerateMDNM) operates entirely on
// self-contained inputs and shares no mutable state with other slices.
func (rr *RecordReader) decodeSliceParallel(h *CompressionHeader, sl *Slice, containerIdx, sliceIdx int) ([]*sam.Record, bool, error) {
	if sl.Header == nil {
		return nil, false, errFormat("container %d slice %d has no header", containerIdx, sliceIdx)
	}
	if sl.Header.NumRecords < 0 {
		return nil, false, errFormat("container %d slice %d declares a negative record count %d",
			containerIdx, sliceIdx, sl.Header.NumRecords)
	}
	src, err := sl.NewSource()
	if err != nil {
		return nil, false, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	refBases, refStart, err := rr.resolveSliceReferenceLocked(sl)
	if err != nil {
		return nil, false, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	dec, err := newRecordDecoder(h, sl.Header, src, rr.refNames, rr.readGroups, refBases, refStart)
	if err != nil {
		return nil, false, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	dec.namePrefix = rr.namePrefix
	recs, err := dec.decodeSliceRecords(sl.Header.NumRecords)
	if err != nil {
		return nil, false, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	if !sl.HasEmbeddedReference() {
		regenerateMDNM(recs, refBases, refStart)
	}
	return recs, dec.needsReference, nil
}

// resolveSliceReferenceLocked is resolveSliceReference guarded by the
// RecordReader's reference mutex. The mutex serialises access to the shared
// referenceResolver's memo (lastContig/lastBases) and any stateful reference
// source so the parallel workers resolve spans safely. The faidx ReadAt fetch,
// the REF_CACHE file read and the REF_PATH source are themselves stateless or
// already mutex-guarded; the lock here protects the memo and keeps a single
// consistent resolution path shared with the sequential code.
func (rr *RecordReader) resolveSliceReferenceLocked(sl *Slice) ([]byte, int32, error) {
	// An embedded reference is self-contained and touches no shared state, so
	// the common embed_ref case resolves without contending on the mutex.
	if sl.Header != nil && sl.Header.RefSeqID >= 0 && sl.HasEmbeddedReference() {
		return rr.resolveSliceReference(sl)
	}
	rr.refMu.Lock()
	defer rr.refMu.Unlock()
	return rr.resolveSliceReference(sl)
}
