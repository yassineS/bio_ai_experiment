package bgzf

import (
	"bytes"
	"compress/flate"
	"errors"
	"io"
	"sync"
)

// MultiWriter is a BGZF writer that compresses blocks concurrently across a
// pool of worker goroutines while preserving the exact on-disk block order of a
// single-threaded Writer. BGZF is block-parallel by construction: every block
// is an independent gzip member, so distributing the deflate work across N
// goroutines and reassembling the framed blocks in their original order yields
// a fully valid BGZF stream (and the same decoded plaintext) regardless of the
// thread count.
//
// MultiWriter satisfies io.WriteCloser. The caller MUST call Close to flush the
// trailing partial block, append the BGZF EOF marker, and surface any deferred
// I/O error. A MultiWriter with a thread count of 1 is equivalent to the
// single-threaded Writer and is handled without spawning goroutines.
type MultiWriter struct {
	w     io.Writer
	level int

	// blockSize is the uncompressed payload size of every full block. It
	// matches the single-threaded Writer's MaxBlockSize so the two paths
	// produce identical block boundaries for the same input.
	blockSize int

	// buf accumulates input until it reaches blockSize, at which point it is
	// dispatched as one job.
	buf []byte

	jobs    chan *blockJob
	results chan *blockJob
	pool    sync.Pool // *flate.Writer

	// nextSeq is the sequence number assigned to the next dispatched block.
	nextSeq int64

	// writerErr is set (under mu) by the collector goroutine when a downstream
	// write fails; subsequent Writes short-circuit with it.
	mu        sync.Mutex
	writerErr error

	wg     sync.WaitGroup
	collWg sync.WaitGroup

	closed  bool
	started bool
}

// blockJob carries one block through the worker pipeline. seq fixes its
// position in the output stream so the collector can emit blocks in order even
// though workers finish out of order.
type blockJob struct {
	seq     int64
	payload []byte // uncompressed input for this block (owned by the job)
	frame   []byte // framed BGZF block, filled in by a worker
	err     error
}

// DefaultBlockJobBuffer bounds the number of in-flight blocks (queued plus
// awaiting collection). It keeps memory use modest while leaving enough slack to
// keep every worker busy.
const DefaultBlockJobBuffer = 256

// NewMultiWriter returns a MultiWriter that compresses to w at the given level
// using up to threads concurrent goroutines. A threads value below 1 is treated
// as 1. Valid levels match flate's: HuffmanOnly, NoCompression (0),
// BestSpeed (1) through BestCompression (9), and DefaultCompression (-1).
func NewMultiWriter(w io.Writer, level, threads int) (*MultiWriter, error) {
	if threads < 1 {
		threads = 1
	}
	// Validate the level eagerly so callers fail fast rather than on first
	// block.
	if _, err := flate.NewWriter(io.Discard, level); err != nil {
		return nil, err
	}
	mw := &MultiWriter{
		w:         w,
		level:     level,
		blockSize: MaxBlockSize,
		buf:       make([]byte, 0, MaxBlockSize),
	}
	mw.pool.New = func() any {
		fw, _ := flate.NewWriter(io.Discard, level)
		return fw
	}
	mw.startWorkers(threads)
	return mw, nil
}

// startWorkers spins up the worker pool and the ordering collector. It is a
// no-op after the first call.
func (mw *MultiWriter) startWorkers(threads int) {
	if mw.started {
		return
	}
	mw.started = true
	mw.jobs = make(chan *blockJob, DefaultBlockJobBuffer)
	mw.results = make(chan *blockJob, DefaultBlockJobBuffer)

	for i := 0; i < threads; i++ {
		mw.wg.Add(1)
		go mw.worker()
	}
	mw.collWg.Add(1)
	go mw.collector()

	// Close results once every worker has exited.
	go func() {
		mw.wg.Wait()
		close(mw.results)
	}()
}

// worker pulls jobs, compresses each payload, and frames the block. The framed
// bytes are attached to the job and forwarded to the collector.
func (mw *MultiWriter) worker() {
	defer mw.wg.Done()
	for job := range mw.jobs {
		fw := mw.pool.Get().(*flate.Writer)
		var deflated bytes.Buffer
		// A block payload is at most MaxBlockSize; its deflate output is
		// bounded by MaxCompressedBlockSize, so preallocate to avoid growth.
		deflated.Grow(MaxCompressedBlockSize)
		fw.Reset(&deflated)
		if _, err := fw.Write(job.payload); err != nil {
			job.err = err
		} else if err := fw.Close(); err != nil {
			job.err = err
		} else {
			frame, err := frameBlock(job.payload, deflated.Bytes(), nil)
			job.frame = frame
			job.err = err
		}
		mw.pool.Put(fw)
		mw.results <- job
	}
}

// collector receives framed blocks out of order, buffers them by sequence
// number, and writes them to the underlying stream strictly in order. It
// records the first downstream write error in writerErr.
func (mw *MultiWriter) collector() {
	defer mw.collWg.Done()
	var want int64
	pending := make(map[int64]*blockJob)
	for job := range mw.results {
		pending[job.seq] = job
		for {
			next, ok := pending[want]
			if !ok {
				break
			}
			delete(pending, want)
			want++
			if next.err != nil {
				mw.setErr(next.err)
				continue
			}
			mw.mu.Lock()
			err := mw.writerErr
			mw.mu.Unlock()
			if err != nil {
				continue
			}
			if _, err := mw.w.Write(next.frame); err != nil {
				mw.setErr(err)
			}
		}
	}
}

// setErr records the first error seen by the collector.
func (mw *MultiWriter) setErr(err error) {
	mw.mu.Lock()
	if mw.writerErr == nil {
		mw.writerErr = err
	}
	mw.mu.Unlock()
}

// err returns the recorded downstream error, if any.
func (mw *MultiWriter) err() error {
	mw.mu.Lock()
	defer mw.mu.Unlock()
	return mw.writerErr
}

// Write buffers p and dispatches full blocks to the worker pool.
func (mw *MultiWriter) Write(p []byte) (int, error) {
	if mw.closed {
		return 0, errors.New("bgzf: write on closed MultiWriter")
	}
	if err := mw.err(); err != nil {
		return 0, err
	}
	total := 0
	for len(p) > 0 {
		space := mw.blockSize - len(mw.buf)
		n := space
		if n > len(p) {
			n = len(p)
		}
		mw.buf = append(mw.buf, p[:n]...)
		p = p[n:]
		total += n
		if len(mw.buf) == mw.blockSize {
			mw.dispatch()
		}
	}
	return total, nil
}

// dispatch hands the currently buffered block to the worker pool and resets the
// buffer. The payload is copied into a fresh slice owned by the job so the
// shared buffer can be reused for the next block.
func (mw *MultiWriter) dispatch() {
	payload := make([]byte, len(mw.buf))
	copy(payload, mw.buf)
	mw.buf = mw.buf[:0]
	job := &blockJob{seq: mw.nextSeq, payload: payload}
	mw.nextSeq++
	mw.jobs <- job
}

// Flush emits any buffered bytes as a BGZF block. It does not emit the EOF
// block. Flush blocks until earlier blocks have drained only insofar as the job
// queue is bounded; it does not wait for in-flight compression to finish.
func (mw *MultiWriter) Flush() error {
	if mw.closed {
		return errors.New("bgzf: flush on closed MultiWriter")
	}
	if err := mw.err(); err != nil {
		return err
	}
	if len(mw.buf) > 0 {
		mw.dispatch()
	}
	return nil
}

// Close flushes the trailing partial block, drains the worker pool, writes the
// BGZF EOF marker, and returns the first deferred I/O error. Close does not
// close the underlying writer. It is safe to call Close more than once.
func (mw *MultiWriter) Close() error {
	if mw.closed {
		return mw.err()
	}
	mw.closed = true

	if len(mw.buf) > 0 {
		mw.dispatch()
	}
	// Signal workers there is no more input and wait for the collector to
	// finish writing every block in order.
	close(mw.jobs)
	mw.collWg.Wait()

	if err := mw.err(); err != nil {
		return err
	}
	// Append the canonical EOF block.
	if _, err := mw.w.Write(EOFBlock); err != nil {
		mw.setErr(err)
		return err
	}
	return nil
}
