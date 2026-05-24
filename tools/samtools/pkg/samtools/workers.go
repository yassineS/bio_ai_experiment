package samtools

import (
	"runtime"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// resolveWorkers maps a CLI -@/--threads value to a concrete worker count.
//
// Upstream samtools' convention for -@ N is "use N additional worker
// threads" (so N=0 means the main thread only, N=4 means 5 threads total).
// We follow that convention loosely: N <= 0 selects 1 worker (the serial
// fast path), N > 0 selects N workers, capped at runtime.NumCPU() so a
// user passing -@ 1024 on a 4-core box doesn't spawn an absurd number of
// goroutines.
func resolveWorkers(n int) int {
	if n <= 0 {
		return 1
	}
	cap := runtime.NumCPU()
	if cap < 1 {
		cap = 1
	}
	if n > cap {
		return cap
	}
	return n
}

// shardJob is one in-flight sort-and-encode unit for the parallel sort
// path: the caller hands the worker a buffer of records and the shard's
// sequence number, the worker sorts the buffer in place and writes the
// records to a fresh temp BAM file, then returns the temp file's path.
type shardJob struct {
	seq  int
	recs []*sam.Record
}

// shardResult is the per-job product of the parallel sort path. Path is the
// temp BAM file written by the worker; err is non-nil when the sort or
// encode step failed.
type shardResult struct {
	seq  int
	path string
	err  error
}

// shardPool is a fixed-size worker pool that sorts and encodes shard
// buffers in parallel while preserving deterministic output ordering: the
// caller assigns each shard a sequence number and the consumer drains
// shardResult values into a small priority buffer to recover submission
// order.
//
// The pool owns no shared mutable state across workers — each shardJob
// carries its own []*sam.Record slice — so the parallel path is
// share-nothing and passes the race detector by construction.
type shardPool struct {
	in     chan shardJob
	out    chan shardResult
	wg     sync.WaitGroup
	doneCh chan struct{}

	errOnce sync.Once
	err     error
}

// newShardPool starts workers worker goroutines, each of which reads
// shardJob values from the input channel, runs work, and posts the result
// to the output channel.
func newShardPool(workers int, work func(shardJob) shardResult) *shardPool {
	if workers < 1 {
		workers = 1
	}
	p := &shardPool{
		in:     make(chan shardJob, workers),
		out:    make(chan shardResult, workers*2),
		doneCh: make(chan struct{}),
	}
	p.wg.Add(workers)
	for i := 0; i < workers; i++ {
		go func() {
			defer p.wg.Done()
			for job := range p.in {
				res := work(job)
				if res.err != nil {
					p.errOnce.Do(func() { p.err = res.err })
				}
				p.out <- res
			}
		}()
	}
	go func() {
		p.wg.Wait()
		close(p.out)
		close(p.doneCh)
	}()
	return p
}

// submit posts a job into the pool. seq is the caller-assigned position
// used to re-order results.
func (p *shardPool) submit(job shardJob) { p.in <- job }

// closeSubmissions signals to the workers that no more jobs will arrive.
// Workers drain the channel and then exit.
func (p *shardPool) closeSubmissions() { close(p.in) }

// results returns the result channel. Callers must drain it (typically
// into a small priority buffer keyed on seq).
func (p *shardPool) results() <-chan shardResult { return p.out }

// firstError returns the first error any worker reported, if any. Safe to
// call after the result channel is fully drained.
func (p *shardPool) firstError() error { return p.err }
