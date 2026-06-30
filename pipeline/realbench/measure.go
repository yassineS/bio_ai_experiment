package realbench

import (
	"io"
	"os"
	"os/exec"
	"syscall"
	"time"
)

// Measurement is the resource usage of one process execution. CPU and memory
// come from wait4(2)'s struct rusage for the child, captured by the Go runtime
// in (*os.ProcessState).SysUsage; wall-clock is measured around the run.
//
// It mirrors pipeline/cmd/realparity.Measurement and pipeline/bench.Measurement
// deliberately: those types' capture cores are unexported, and the CLAUDE.md
// guidance explicitly permits replicating the ~20-line rusage capture rather
// than exporting new surface. The reduction rule (wall/CPU = min, RSS = max) is
// the same one those harnesses document.
type Measurement struct {
	Wall     time.Duration // wall-clock elapsed
	CPUUser  time.Duration // user-mode CPU (ru_utime)
	CPUSys   time.Duration // kernel-mode CPU (ru_stime)
	MaxRSSKB int64         // peak resident set size in KiB (platform-normalised)
}

// CPUTotal is user+sys CPU time.
func (m Measurement) CPUTotal() time.Duration { return m.CPUUser + m.CPUSys }

// timevalDuration converts a syscall.Timeval (seconds + microseconds) to a
// time.Duration. The field types differ across platforms (int32/int64), so the
// conversion goes through int64.
func timevalDuration(tv syscall.Timeval) time.Duration {
	return time.Duration(int64(tv.Sec))*time.Second + time.Duration(int64(tv.Usec))*time.Microsecond
}

// runOnce executes bin with args, STREAMING the child's stdout into sink rather
// than buffering it, and returns the resource usage. This keeps memory bounded
// to whatever sink retains (e.g. an md5 hash + a 64 KiB head window, or
// io.Discard for timed-only reps) so a command whose stdout is many GB — like
// `samtools view` decoding a multi-GB BAM to SAM — never materialises in RAM.
//
// stderr is discarded (symmetric for both sides). When stdinPath is non-empty it
// is opened and fed as the child's stdin; when env is non-nil it replaces the
// child environment (used to pass REF_PATH=/dev/null so CRAM never reaches the
// network). When sink is nil, stdout is discarded. workDir, when non-empty, sets
// the child's working directory (for tools that write output files relative to
// cwd, e.g. fastp / mosdepth).
//
// A non-nil error means the process failed (non-zero exit or spawn error); the
// caller treats that as a hard cell failure, never as a perf signal.
func runOnce(bin string, args []string, stdinPath, workDir string, env []string, sink io.Writer) (Measurement, error) {
	cmd := exec.Command(bin, args...)
	if env != nil {
		cmd.Env = env
	}
	if workDir != "" {
		cmd.Dir = workDir
	}
	if stdinPath != "" {
		f, err := os.Open(stdinPath)
		if err != nil {
			return Measurement{}, err
		}
		defer f.Close()
		cmd.Stdin = f
	}
	if sink == nil {
		sink = io.Discard
	}
	cmd.Stdout = sink
	cmd.Stderr = nil

	start := time.Now()
	err := cmd.Run()
	wall := time.Since(start)

	m := Measurement{Wall: wall}
	if ps := cmd.ProcessState; ps != nil {
		if ru, ok := ps.SysUsage().(*syscall.Rusage); ok {
			m.CPUUser = timevalDuration(ru.Utime)
			m.CPUSys = timevalDuration(ru.Stime)
			m.MaxRSSKB = maxRSSToKiB(int64(ru.Maxrss))
		}
	}
	return m, err
}

// repeatRun runs the command reps times and reduces the samples: wall and CPU
// take the MINIMUM (the run least perturbed by scheduler/IO noise — the standard
// "how fast can it go" choice), RSS takes the MAXIMUM (the true peak). reps is
// clamped to >= 1. The first error aborts.
//
// Output is streamed, never buffered: exactly ONE rep (the last) writes its
// stdout into cmpSink — the provenance-filtering digester whose digest parity
// compares — while every other rep discards its stdout. Because a deterministic
// command emits identical bytes on every rep, comparing the last rep's digest is
// representative; running the digester on only one rep keeps memory bounded and
// avoids hashing the same multi-GB stream repeatedly. cmpSink may be nil (e.g.
// quickcheck cells that compare exit status, or file-producing cells whose
// stdout is not the comparison stream), in which case every rep discards.
func repeatRun(reps int, bin string, args []string, stdinPath, workDir string, env []string, cmpSink io.Writer) (Measurement, error) {
	if reps < 1 {
		reps = 1
	}
	var best Measurement
	for i := 0; i < reps; i++ {
		// Only the final rep feeds the comparison digester; earlier reps discard.
		var sink io.Writer = io.Discard
		if i == reps-1 {
			sink = cmpSink
		}
		m, err := runOnce(bin, args, stdinPath, workDir, env, sink)
		if err != nil {
			return Measurement{}, err
		}
		if i == 0 {
			best = m
			continue
		}
		if m.Wall < best.Wall {
			best.Wall = m.Wall
		}
		if m.CPUTotal() < best.CPUTotal() {
			best.CPUUser, best.CPUSys = m.CPUUser, m.CPUSys
		}
		if m.MaxRSSKB > best.MaxRSSKB {
			best.MaxRSSKB = m.MaxRSSKB
		}
	}
	return best, nil
}
