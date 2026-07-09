package skewer

import (
	"errors"
	"fmt"
	"io"
	"runtime"
	"sync"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// isEOF reports whether err is io.EOF, matching the sequential loop's
// termination condition.
func isEOF(err error) bool {
	return errors.Is(err, io.EOF)
}

// peJob carries one input pair (with its original lengths) plus a monotonic
// sequence id so the collector can reassemble the results in strict input
// order regardless of which worker finishes first.
type peJob struct {
	id       uint64
	rec1     *fastq.Record
	rec2     *fastq.Record
	origLen1 int
	origLen2 int
}

// peResult is the trimmed outcome for one pair, tagged with the same id as its
// job so the collector can order the writes deterministically.
type peResult struct {
	id       uint64
	rec1     *fastq.Record
	rec2     *fastq.Record
	pos1     int
	pos2     int
	origLen1 int
	origLen2 int
	keep     bool // passes the minLen filter (both mates)
	adapter3 int  // AdapterFound3 delta from trimPairWithPE for this pair (0 or 1)
}

// trimPairedEndParallel is the multithreaded PE matrix-mode path. It runs one
// reader goroutine (tagging each pair with a monotonic id), fans the pure
// per-pair CPU work (trimPairWithPE) out over opts.Threads workers, and a
// single collector reassembles results in id order and writes writer1/writer2
// in strict input order. The output is byte-identical to the sequential path
// for any thread count.
//
// Precondition (enforced by the caller): opts.PEMatrixMode is set and
// opts.UMILength == 0. trimPairWithPE operates on local copies of each pair's
// buffers and never mutates shared state, so the only shared mutation is the
// stats merge, which is done per-worker and combined deterministically here.
func trimPairedEndParallel(reader1, reader2 *fastq.Reader, writer1, writer2 *fastq.Writer,
	stats *TrimStats, opts TrimOptions, startTime time.Time) (*TrimStats, error) {

	nWorkers := opts.Threads
	if nWorkers < 1 {
		nWorkers = 1
	}
	if max := runtime.NumCPU(); nWorkers > max*4 {
		// Guard against absurd thread counts; work is CPU-bound.
		nWorkers = max * 4
	}

	// Buffer the channels generously so the reader and collector rarely block
	// the workers. The ordering guarantee comes from the id-keyed reassembly in
	// the collector, not from channel order.
	jobs := make(chan peJob, nWorkers*8)
	results := make(chan peResult, nWorkers*8)

	// readErr is set by the reader goroutine; guarded by readErrMu because the
	// collector reports it after draining.
	var readErr error
	var readErrMu sync.Mutex

	// Reader goroutine: read R1/R2 pairs, tag with a monotonic id, and enqueue.
	go func() {
		defer close(jobs)
		var id uint64
		for {
			record1, err1 := reader1.Read()
			record2, err2 := reader2.Read()
			if err1 != nil || err2 != nil {
				// Either input exhausted => normal termination, mirroring the
				// sequential loop which stops on the first EOF from either mate.
				if isEOF(err1) || isEOF(err2) {
					return
				}
				readErrMu.Lock()
				if err1 != nil {
					readErr = fmt.Errorf("error reading first input: %w", err1)
				} else {
					readErr = fmt.Errorf("error reading second input: %w", err2)
				}
				readErrMu.Unlock()
				return
			}
			jobs <- peJob{
				id:       id,
				rec1:     record1,
				rec2:     record2,
				origLen1: len(record1.Sequence),
				origLen2: len(record2.Sequence),
			}
			id++
		}
	}()

	// Workers: pure per-pair CPU. trimPairWithPE works on copies of the buffers
	// and, when passed a nil stats, never mutates shared state — so no worker
	// touches the TrimStats. All counting happens in the single collector below,
	// making the totals order-independent (thus thread-count-invariant).
	var wg sync.WaitGroup
	wg.Add(nWorkers)
	for w := 0; w < nWorkers; w++ {
		go func() {
			defer wg.Done()
			// A per-worker scratch TrimStats captures trimPairWithPE's only
			// side effect (the AdapterFound3 increment) without touching the
			// shared stats; its delta is carried on the result and merged in the
			// collector so the counter matches the sequential path exactly.
			var local TrimStats
			for job := range jobs {
				before := local.AdapterFound3
				out1, out2, pos1, pos2 := trimPairWithPE(job.rec1, job.rec2, opts, &local)
				results <- peResult{
					id:       job.id,
					rec1:     out1,
					rec2:     out2,
					pos1:     pos1,
					pos2:     pos2,
					origLen1: job.origLen1,
					origLen2: job.origLen2,
					keep:     pos1 >= opts.MinLength && pos2 >= opts.MinLength,
					adapter3: local.AdapterFound3 - before,
				}
			}
		}()
	}

	// Close results once all workers are done, so the collector's range ends.
	go func() {
		wg.Wait()
		close(results)
	}()

	// Collector: reassemble results in id order using a ring/pending map keyed
	// on nextID, writing writer1/writer2 in strict input order. Stats are merged
	// here as each pair commits, so the totals are order-independent.
	pending := make(map[uint64]peResult)
	var nextID uint64
	var writeErr error

	commit := func(r peResult) error {
		// Every pair counts toward the totals, exactly as the sequential loop
		// does before the keep filter.
		stats.TotalReads += 2
		stats.TotalBases += int64(r.origLen1 + r.origLen2)
		stats.AdapterFound3 += r.adapter3
		if !r.keep {
			stats.DiscardedReads += 2
			return nil
		}
		if err := writer1.Write(r.rec1); err != nil {
			return fmt.Errorf("error writing first output: %w", err)
		}
		if err := writer2.Write(r.rec2); err != nil {
			return fmt.Errorf("error writing second output: %w", err)
		}
		if r.pos1 < r.origLen1 || r.pos2 < r.origLen2 {
			stats.TrimmedReads++
			stats.TrimmedBases += int64((r.origLen1 - r.pos1) + (r.origLen2 - r.pos2))
		}
		return nil
	}

	for r := range results {
		pending[r.id] = r
		for {
			nr, ok := pending[nextID]
			if !ok {
				break
			}
			delete(pending, nextID)
			nextID++
			if writeErr == nil {
				if err := commit(nr); err != nil {
					writeErr = err
				}
			}
		}
	}

	// Any straggler still buffered (should not happen once results is drained,
	// but flush defensively in id order).
	for len(pending) > 0 {
		nr, ok := pending[nextID]
		if !ok {
			break
		}
		delete(pending, nextID)
		nextID++
		if writeErr == nil {
			if err := commit(nr); err != nil {
				writeErr = err
			}
		}
	}

	readErrMu.Lock()
	rErr := readErr
	readErrMu.Unlock()
	if rErr != nil {
		return stats, rErr
	}
	if writeErr != nil {
		return stats, writeErr
	}

	if err := writer1.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing first output: %w", err)
	}
	if err := writer2.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing second output: %w", err)
	}

	stats.ProcessingTime = time.Since(startTime)
	return stats, nil
}
