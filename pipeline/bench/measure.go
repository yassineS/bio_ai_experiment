// Package bench holds the performance side of the parity pipeline. Alongside the
// testing.B micro-benchmarks (bench_test.go), it exposes a standalone
// measurement core and a scale-sweeping runner (cmd/parity-bench) that times OUR
// binary against the vendored UPSTREAM binary and records, for each, the three
// resource axes a manuscript needs: wall-clock, CPU time (user+sys), and peak
// resident memory (max RSS). Running the same matrix across the small / medium /
// large fixture tiers turns those point measurements into scalability curves
// versus input size (read count, coverage, variant count, interval count).
package bench

import (
	"os"
	"os/exec"
	"syscall"
	"time"
)

// Measurement is the resource usage of one process execution. CPU and memory
// come from wait4(2)'s struct rusage for the child, captured by the Go runtime
// in (*os.ProcessState).SysUsage; wall-clock is measured around Run.
type Measurement struct {
	Wall     time.Duration // wall-clock elapsed
	CPUUser  time.Duration // user-mode CPU (ru_utime)
	CPUSys   time.Duration // kernel-mode CPU (ru_stime)
	MaxRSSKB int64         // peak resident set size in KiB (ru_maxrss on Linux)
}

// CPUTotal is user+sys CPU time.
func (m Measurement) CPUTotal() time.Duration { return m.CPUUser + m.CPUSys }

// timevalDuration converts a syscall.Timeval (seconds + microseconds) to a
// time.Duration. The field types differ across platforms (int32/int64), so the
// conversion goes through int64.
func timevalDuration(tv syscall.Timeval) time.Duration {
	return time.Duration(int64(tv.Sec))*time.Second + time.Duration(int64(tv.Usec))*time.Microsecond
}

// runMeasured executes bin with args and returns its resource usage. stdout and
// stderr are sent to the null device (symmetric for both sides) unless
// stdoutPath is non-empty, in which case stdout is written there — used for the
// streaming operations whose result is the stdout payload, so the write cost is
// included on both sides. stdinPath, when non-empty, is opened and fed as stdin.
//
// A non-nil error means the process failed (non-zero exit or spawn error); the
// caller treats that as a hard failure for the cell, never as a perf signal.
func runMeasured(bin string, args []string, stdinPath, stdoutPath string) (Measurement, error) {
	cmd := exec.Command(bin, args...)

	if stdinPath != "" {
		f, err := os.Open(stdinPath)
		if err != nil {
			return Measurement{}, err
		}
		defer f.Close()
		cmd.Stdin = f
	}

	if stdoutPath != "" {
		f, err := os.Create(stdoutPath)
		if err != nil {
			return Measurement{}, err
		}
		defer f.Close()
		cmd.Stdout = f
	} else {
		cmd.Stdout = nil // null device
	}
	cmd.Stderr = nil

	start := time.Now()
	err := cmd.Run()
	wall := time.Since(start)

	m := Measurement{Wall: wall}
	if ps := cmd.ProcessState; ps != nil {
		if ru, ok := ps.SysUsage().(*syscall.Rusage); ok {
			m.CPUUser = timevalDuration(ru.Utime)
			m.CPUSys = timevalDuration(ru.Stime)
			m.MaxRSSKB = int64(ru.Maxrss)
		}
	}
	return m, err
}

// repeatMeasured runs runMeasured reps times and reduces the samples: wall and
// CPU take the MINIMUM (the run least perturbed by scheduler/IO noise — the
// standard choice for "how fast can it go"), while max RSS takes the MAXIMUM
// (the true peak across runs). reps must be >= 1. The first error aborts.
//
// It is a thin convenience wrapper over repeatMeasuredSamples that discards the
// raw per-rep samples; callers that need the full distribution (for median/IQR
// and the ratio CI) call repeatMeasuredSamples directly.
func repeatMeasured(reps int, bin string, args []string, stdinPath, stdoutPath string) (Measurement, error) {
	best, _, err := repeatMeasuredSamples(reps, bin, args, stdinPath, stdoutPath)
	return best, err
}

// repeatMeasuredSamples runs runMeasured reps times and returns BOTH the reduced
// Measurement (min wall, min CPU, max RSS — see repeatMeasured) AND the raw
// per-rep samples in execution order. The raw samples are what the manuscript's
// distribution-level statistics (median, IQR, bootstrap ratio CI) are computed
// from; the reduced Measurement is kept for backward compatibility with the
// existing min/max report fields. reps must be >= 1. The first error aborts.
func repeatMeasuredSamples(reps int, bin string, args []string, stdinPath, stdoutPath string) (Measurement, []Measurement, error) {
	if reps < 1 {
		reps = 1
	}
	var best Measurement
	samples := make([]Measurement, 0, reps)
	for i := 0; i < reps; i++ {
		m, err := runMeasured(bin, args, stdinPath, stdoutPath)
		if err != nil {
			return Measurement{}, nil, err
		}
		samples = append(samples, m)
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
	return best, samples, nil
}
