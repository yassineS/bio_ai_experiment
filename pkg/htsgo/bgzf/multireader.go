package bgzf

import (
	"bytes"
	"encoding/binary"
	kflate "github.com/klauspost/compress/flate"
	"hash/crc32"
	"io"
	"sync"
)

// MultiReader is a BGZF reader that decompresses blocks concurrently across a
// pool of worker goroutines while delivering the decoded bytes in their exact
// original order. BGZF is block-parallel by construction: every block is an
// independent gzip member, so the deflate work can be spread across N
// goroutines and the decoded payloads reassembled in sequence to reproduce the
// single-threaded byte stream exactly.
//
// MultiReader satisfies io.Reader (and io.Closer). Its output is byte-for-byte
// identical to the sequential Reader's regardless of the worker count, so it is
// a drop-in accelerator for the BAM record parser, which still consumes the
// decoded stream sequentially. A worker count below 2 falls back to the
// sequential Reader so callers never pay the goroutine overhead for `-t 1`.
//
// MultiReader does not support seeking or virtual offsets; it is meant for the
// streaming full-file decode path (e.g. mosdepth's whole-BAM sweep). Callers
// that need region queries should use the seekable Reader.
//
// It assumes a source that always makes progress to EOF (a file or any reader
// whose Read eventually returns io.EOF). Close unblocks a collector stalled on
// an abandoned pipe, but it does NOT interrupt the feeder goroutine if the
// underlying source blocks forever; on such a never-terminating streaming
// source the worker goroutines would not unwind. This matches the intended
// file-backed BAM use and is not a concern for it.
type MultiReader struct {
	// seq is the sequential reader used when threads < 2.
	seq *Reader

	r       io.Reader
	jobs    chan *decodeJob
	results chan *decodeJob
	pool    sync.Pool // *flate.Reader (io.ReadCloser)

	// pending buffers out-of-order decoded blocks by sequence number; the
	// collector goroutine writes them into the ordered pipe in sequence.
	pr *io.PipeReader
	pw *io.PipeWriter

	wg      sync.WaitGroup
	started bool
	closed  bool
}

// decodeJob carries one compressed block through the worker pipeline. seq fixes
// its position in the decoded stream so the collector can emit payloads in
// order even though workers finish out of order. An eof job (eof=true) marks
// the BGZF EOF block and carries no payload.
type decodeJob struct {
	seq      int64
	deflated []byte
	wantCRC  uint32
	wantISZ  uint32
	decoded  []byte
	eof      bool
	err      error
}

// NewMultiReader returns a MultiReader that decodes BGZF bytes from r using up
// to threads concurrent goroutines. A threads value below 2 yields a reader
// backed by the sequential Reader (no goroutines spawned). The decoded output
// is identical for any threads value.
func NewMultiReader(r io.Reader, threads int) (*MultiReader, error) {
	if threads < 2 {
		seq, err := NewReader(r)
		if err != nil {
			return nil, err
		}
		return &MultiReader{seq: seq}, nil
	}
	mr := &MultiReader{r: r}
	mr.pool.New = func() any { return kflate.NewReader(bytes.NewReader(nil)) }
	mr.pr, mr.pw = io.Pipe()
	mr.start(threads)
	return mr, nil
}

// decodeJobBuffer bounds the number of in-flight blocks (queued plus awaiting
// collection). It keeps memory modest while leaving slack to keep workers busy.
const decodeJobBuffer = 256

// start spins up the reader, the worker pool, and the ordering collector.
func (mr *MultiReader) start(threads int) {
	mr.started = true
	mr.jobs = make(chan *decodeJob, decodeJobBuffer)
	mr.results = make(chan *decodeJob, decodeJobBuffer)

	for i := 0; i < threads; i++ {
		mr.wg.Add(1)
		go mr.worker()
	}
	go func() {
		mr.wg.Wait()
		close(mr.results)
	}()
	go mr.feed()
	go mr.collect()
}

// feed reads framed blocks from the compressed stream sequentially and
// dispatches each one's deflate payload to the worker pool. Reading the framed
// blocks is cheap (no decompression); the expensive inflate runs in workers.
func (mr *MultiReader) feed() {
	defer close(mr.jobs)
	var seq int64
	for {
		hdr, err := readBlockHeader(mr.r)
		if err != nil {
			if err == io.EOF {
				return
			}
			mr.jobs <- &decodeJob{seq: seq, err: err}
			return
		}
		deflatedLen := hdr.compressedSize - hdr.headerLen - 8
		if deflatedLen < 0 {
			mr.jobs <- &decodeJob{seq: seq, err: ErrBadBSIZE}
			return
		}
		deflated := make([]byte, deflatedLen)
		if _, err := io.ReadFull(mr.r, deflated); err != nil {
			mr.jobs <- &decodeJob{seq: seq, err: ioErrUnexpected(err)}
			return
		}
		var footer [8]byte
		if _, err := io.ReadFull(mr.r, footer[:]); err != nil {
			mr.jobs <- &decodeJob{seq: seq, err: ioErrUnexpected(err)}
			return
		}
		job := &decodeJob{
			seq:      seq,
			deflated: deflated,
			wantCRC:  binary.LittleEndian.Uint32(footer[0:4]),
			wantISZ:  binary.LittleEndian.Uint32(footer[4:8]),
		}
		if job.wantISZ == 0 && deflatedLen == 2 {
			// The canonical empty BGZF EOF block; the stream ends here.
			job.eof = true
			mr.jobs <- job
			return
		}
		mr.jobs <- job
		seq++
	}
}

// worker inflates one block's deflate payload, verifies its CRC32/ISIZE, and
// forwards the decoded bytes to the collector tagged with its sequence number.
func (mr *MultiReader) worker() {
	defer mr.wg.Done()
	for job := range mr.jobs {
		if job.err != nil || job.eof {
			mr.results <- job
			continue
		}
		fr := mr.pool.Get().(io.ReadCloser)
		if rs, ok := fr.(kflate.Resetter); ok {
			_ = rs.Reset(bytes.NewReader(job.deflated), nil)
		}
		var buf bytes.Buffer
		buf.Grow(int(job.wantISZ))
		if _, err := io.Copy(&buf, fr); err != nil {
			job.err = err
			mr.pool.Put(fr)
			mr.results <- job
			continue
		}
		mr.pool.Put(fr)
		decoded := buf.Bytes()
		if uint32(len(decoded)) != job.wantISZ {
			job.err = ErrISIZE
		} else if crc32.ChecksumIEEE(decoded) != job.wantCRC {
			job.err = ErrChecksum
		} else {
			job.decoded = decoded
		}
		mr.results <- job
	}
}

// collect receives decoded blocks out of order, buffers them by sequence
// number, and writes their payloads into the ordered pipe strictly in
// sequence. It closes the pipe writer with the first error encountered (or nil
// on a clean EOF).
func (mr *MultiReader) collect() {
	var want int64
	pending := make(map[int64]*decodeJob)
	var finalErr error
	sawEOF := false
	for job := range mr.results {
		pending[job.seq] = job
		for {
			next, ok := pending[want]
			if !ok {
				break
			}
			delete(pending, want)
			if next.err != nil {
				finalErr = next.err
				goto done
			}
			if next.eof {
				sawEOF = true
				goto done
			}
			if _, err := mr.pw.Write(next.decoded); err != nil {
				finalErr = err
				goto done
			}
			want++
		}
	}
done:
	if finalErr == nil && !sawEOF {
		// The stream ended without the canonical EOF block.
		finalErr = ErrTruncated
	}
	mr.pw.CloseWithError(finalErr)
	// Drain any straggler results so workers/feed never block on a full
	// channel after we have stopped reading from the pipe.
	for range mr.results {
	}
}

// Read delivers the next decoded bytes in order.
func (mr *MultiReader) Read(p []byte) (int, error) {
	if mr.seq != nil {
		return mr.seq.Read(p)
	}
	return mr.pr.Read(p)
}

// Close releases resources. For the sequential fallback it closes the inner
// Reader; for the parallel path it tears down the pipe so the collector and
// workers unwind. Close does not close the underlying source reader.
func (mr *MultiReader) Close() error {
	if mr.seq != nil {
		return mr.seq.Close()
	}
	if mr.closed {
		return nil
	}
	mr.closed = true
	if mr.pr != nil {
		// Unblock the collector if a reader abandoned the stream early.
		_ = mr.pr.CloseWithError(io.ErrClosedPipe)
	}
	return nil
}
