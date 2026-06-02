package bgzf

import (
	"errors"
	"io"
	"runtime"
	"sync"
)

// ParallelWriter is a drop-in replacement for Writer that compresses BGZF
// blocks in parallel using a fixed worker pool. The block boundaries and
// emission order are identical to a serial Writer at the same compression
// level: bytes written are partitioned into MaxBlockSize-byte uncompressed
// chunks (with explicit Flush boundaries honored), and each chunk is
// encoded into one gzip member. Workers compress chunks concurrently, but
// finished blocks are drained in submission order through a sequenced
// reassembly stage so the on-disk byte sequence is byte-identical to the
// serial path. Close emits the canonical EOFBlock just like Writer.
//
// The parallel path passes the race detector by construction: each shard
// owns its own input slice and its own libdeflate scratch buffer, and the
// only shared state is the worker channels and the sequenced output
// drainer.
type ParallelWriter struct {
	w     io.Writer
	level int

	buf [MaxBlockSize]byte
	n   int

	// nextSeq is the sequence number to assign the next block submitted.
	nextSeq int

	// jobs feeds workers. Each parWorkJob is a freshly-allocated payload
	// slice (workers may keep references after the producer has moved on)
	// plus the sequence number.
	jobs chan parWorkJob
	// done is closed by the drainer when it has finished writing all
	// results.
	done chan struct{}

	// drainerErr is the first error the drainer hit while writing a
	// result block. workerErr is the first error any worker hit while
	// compressing. err is the sticky surface error returned to the
	// caller; consulted by Write/Flush/Close.
	mu         sync.Mutex
	err        error
	workerWG   sync.WaitGroup
	drainerWG  sync.WaitGroup
	closed     bool
	submitting bool
}

type parWorkJob struct {
	seq     int
	payload []byte
}

type parWorkResult struct {
	seq   int
	block []byte // full gzip member, ready to write
	err   error
}

// NewParallelWriter returns a ParallelWriter using DefaultCompression and
// the requested number of worker goroutines. workers <= 1 still spawns one
// background worker and one drainer (so behaviour is uniform), but the
// caller should generally just use NewWriter for the serial path.
func NewParallelWriter(w io.Writer, workers int) *ParallelWriter {
	pw, _ := NewParallelWriterLevel(w, DefaultCompression, workers)
	return pw
}

// NewParallelWriterLevel returns a ParallelWriter compressing at the given
// level with the requested worker count. The returned writer must be
// Close()d to flush buffered bytes, drain workers, and emit the BGZF EOF
// block.
func NewParallelWriterLevel(w io.Writer, level, workers int) (*ParallelWriter, error) {
	// Validate level by trying to construct a serial Writer with it.
	if _, err := NewWriterLevel(io.Discard, level); err != nil {
		return nil, err
	}
	if workers < 1 {
		workers = 1
	}
	if cpu := runtime.NumCPU(); workers > cpu && cpu > 0 {
		workers = cpu
	}
	p := &ParallelWriter{
		w:     w,
		level: level,
		// Buffer jobs slightly so producers don't always block on
		// every block boundary; sized at 2x worker count.
		jobs: make(chan parWorkJob, workers*2),
		done: make(chan struct{}),
	}
	// Each worker has its own results channel; the drainer multiplexes
	// them into a min-heap keyed on seq so output goes out in order.
	results := make(chan parWorkResult, workers*4)
	p.workerWG.Add(workers)
	for i := 0; i < workers; i++ {
		go p.worker(results)
	}
	// When all workers exit, close the results channel so the drainer
	// can finish.
	go func() {
		p.workerWG.Wait()
		close(results)
	}()
	p.drainerWG.Add(1)
	go p.drainer(results)
	return p, nil
}

// worker pulls jobs and compresses them, posting results to the shared
// results channel. Each worker owns its own serial Writer used solely as
// a per-block encoder via writeBlockTo.
func (p *ParallelWriter) worker(results chan<- parWorkResult) {
	defer p.workerWG.Done()
	enc := &blockEncoder{level: p.level}
	for job := range p.jobs {
		block, err := enc.encode(job.payload)
		results <- parWorkResult{seq: job.seq, block: block, err: err}
	}
}

// drainer pulls results from workers, reorders by seq using a small heap,
// and writes the contiguous prefix to the underlying io.Writer.
func (p *ParallelWriter) drainer(results <-chan parWorkResult) {
	defer p.drainerWG.Done()
	defer close(p.done)
	next := 0
	pending := map[int]parWorkResult{}
	for r := range results {
		if r.err != nil {
			p.setErr(r.err)
			// Keep draining to let workers exit cleanly.
			continue
		}
		pending[r.seq] = r
		for {
			cur, ok := pending[next]
			if !ok {
				break
			}
			delete(pending, next)
			if _, err := p.w.Write(cur.block); err != nil {
				p.setErr(err)
				next++
				continue
			}
			next++
		}
	}
}

func (p *ParallelWriter) setErr(err error) {
	p.mu.Lock()
	if p.err == nil {
		p.err = err
	}
	p.mu.Unlock()
}

func (p *ParallelWriter) getErr() error {
	p.mu.Lock()
	defer p.mu.Unlock()
	return p.err
}

// Write appends p to the writer's current block, flushing blocks to the
// worker pool whenever a block fills. Mirrors Writer.Write semantically.
func (p *ParallelWriter) Write(buf []byte) (int, error) {
	if err := p.getErr(); err != nil {
		return 0, err
	}
	if p.closed {
		return 0, errors.New("bgzf: write on closed ParallelWriter")
	}
	total := 0
	for len(buf) > 0 {
		space := MaxBlockSize - p.n
		if space == 0 {
			if err := p.submitBlock(); err != nil {
				return total, err
			}
			space = MaxBlockSize
		}
		take := len(buf)
		if take > space {
			take = space
		}
		copy(p.buf[p.n:], buf[:take])
		p.n += take
		total += take
		buf = buf[take:]
	}
	return total, nil
}

// Flush emits the buffered bytes as a BGZF block (even if not full), like
// Writer.Flush. Safe on an empty buffer.
func (p *ParallelWriter) Flush() error {
	if err := p.getErr(); err != nil {
		return err
	}
	if p.n == 0 {
		return nil
	}
	return p.submitBlock()
}

// FlushTry mirrors Writer.FlushTry: flushes the current partial block iff
// appending size more bytes would overflow MaxBlockSize.
func (p *ParallelWriter) FlushTry(size int) error {
	if err := p.getErr(); err != nil {
		return err
	}
	if p.n+size > MaxBlockSize {
		return p.Flush()
	}
	return nil
}

// WriteRaw flushes any buffered bytes (so block ordering is preserved)
// and then writes p directly to the underlying stream without
// compressing it. p is expected to be one or more complete BGZF members
// (e.g. blocks obtained from ReadRawBlock). Mirrors Writer.WriteRaw.
//
// To preserve sequence order across mixed compressed/verbatim members,
// the buffered bytes are flushed, the worker pool is drained, and only
// then are the raw bytes written. After WriteRaw returns the pool is
// quiescent and subsequent Write/Flush calls resume normally.
func (p *ParallelWriter) WriteRaw(buf []byte) error {
	if err := p.getErr(); err != nil {
		return err
	}
	if p.closed {
		return errors.New("bgzf: WriteRaw on closed ParallelWriter")
	}
	if p.n > 0 {
		if err := p.submitBlock(); err != nil {
			return err
		}
	}
	// Drain the worker pool so all preceding submitted blocks have been
	// written before the raw bytes appear in the output.
	if err := p.quiesce(); err != nil {
		return err
	}
	if _, err := p.w.Write(buf); err != nil {
		p.setErr(err)
		return err
	}
	return nil
}

// quiesce waits for the worker pool to drain all currently-submitted
// blocks to the underlying writer. It does NOT close the jobs channel;
// after quiesce returns the producer may submit further blocks and the
// workers/drainer remain running.
func (p *ParallelWriter) quiesce() error {
	// Synchronise via a marker job: submit a sentinel and wait for the
	// drainer to observe a contiguous prefix up through it. The simplest
	// approach is to send a zero-payload job that the worker turns into
	// nothing on the wire; but the drainer's order-recovery state must
	// see the seq slot. We instead repurpose the existing pipeline by
	// sending a tiny dummy block IF the producer has buffered work to
	// flush.
	//
	// In practice WriteRaw is rare (used only by reheader/cat) and
	// performance here is not critical, so we use the heavy-handed
	// option: close and recreate the worker pool. That guarantees all
	// in-flight work has been drained before we write the raw bytes.
	close(p.jobs)
	<-p.done
	p.workerWG.Wait()
	p.drainerWG.Wait()
	if err := p.getErr(); err != nil {
		return err
	}
	// Recreate the pool so the writer is reusable.
	workers := cap(p.jobs) / 2
	if workers < 1 {
		workers = 1
	}
	p.jobs = make(chan parWorkJob, workers*2)
	p.done = make(chan struct{})
	results := make(chan parWorkResult, workers*4)
	p.workerWG.Add(workers)
	for i := 0; i < workers; i++ {
		go p.worker(results)
	}
	go func() {
		p.workerWG.Wait()
		close(results)
	}()
	p.drainerWG.Add(1)
	go p.drainer(results)
	// Reset nextSeq for the new pool — the drainer's expected-next
	// counter starts at 0 again because it's a fresh drainer.
	p.nextSeq = 0
	return nil
}

// submitBlock hands the current buffer off to the worker pool. The buffer
// payload is copied so the producer can immediately reuse p.buf.
func (p *ParallelWriter) submitBlock() error {
	if p.n == 0 {
		return nil
	}
	payload := make([]byte, p.n)
	copy(payload, p.buf[:p.n])
	p.n = 0
	seq := p.nextSeq
	p.nextSeq++
	// The send may block if workers are busy — that's the desired
	// back-pressure. Also surface any drainer-side error early.
	if err := p.getErr(); err != nil {
		return err
	}
	p.jobs <- parWorkJob{seq: seq, payload: payload}
	return nil
}

// Close flushes any buffered bytes, drains the worker pool, emits the
// BGZF EOF block, and releases resources. Does not close the underlying
// writer.
func (p *ParallelWriter) Close() error {
	if p.closed {
		return nil
	}
	p.closed = true
	if err := p.getErr(); err != nil {
		// Still need to drain pool to let goroutines exit.
		close(p.jobs)
		<-p.done
		return err
	}
	// Flush any trailing partial block before closing the jobs channel.
	if p.n > 0 {
		if err := p.submitBlock(); err != nil {
			close(p.jobs)
			<-p.done
			return err
		}
	}
	close(p.jobs)
	// Wait for drainer to finish writing everything.
	<-p.done
	if err := p.getErr(); err != nil {
		return err
	}
	// Emit the canonical EOF block.
	if _, err := p.w.Write(EOFBlock); err != nil {
		p.setErr(err)
		return err
	}
	return nil
}

// blockEncoder is a per-worker encoder that compresses a single BGZF block.
// It is a thin wrapper around the same code path as Writer.encodeBlock, but
// emits to a private buffer instead of the shared io.Writer.
type blockEncoder struct {
	level int
	// Reuse a single Writer instance per worker for its deflatePayload
	// path. Its w field is unused (we never call encodeBlock through it).
	tmp *Writer
}

// encode builds a single complete gzip member for payload and returns it
// as a fresh byte slice.
func (e *blockEncoder) encode(payload []byte) ([]byte, error) {
	if e.tmp == nil {
		wr, err := NewWriterLevel(io.Discard, e.level)
		if err != nil {
			return nil, err
		}
		e.tmp = wr
	}
	deflated, err := e.tmp.deflatePayload(payload)
	if err != nil {
		return nil, err
	}
	return buildBGZFMember(payload, deflated)
}
