package main

import (
	"bytes"
	"os"
	"os/exec"
	"syscall"
	"time"
)

// Measurement is the resource usage of one process execution. CPU and memory
// come from wait4(2)'s struct rusage for the child, captured by the Go runtime
// in (*os.ProcessState).SysUsage; wall-clock is measured around the run.
//
// It mirrors pipeline/bench.Measurement deliberately: that type's capture core
// (runMeasured / repeatMeasured) is unexported, and the CLAUDE.md guidance
// explicitly permits replicating its ~20-line rusage capture rather than
// exporting new surface. The reduction rule (wall/CPU = min, RSS = max) is the
// same one bench documents.
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

// runResult bundles a process's resource usage with its captured stdout. The
// stdout payload is what parity compares; the Measurement is the perf signal.
type runResult struct {
	Meas   Measurement
	Stdout []byte
}

// runOnce executes bin with args, captures stdout into memory, and returns the
// resource usage. stderr is discarded (symmetric for both sides). When
// stdinPath is non-empty it is opened and fed as the child's stdin; when env is
// non-nil it replaces the child environment (used to pass REF_PATH=/dev/null so
// CRAM never reaches the network).
//
// A non-nil error means the process failed (non-zero exit or spawn error); the
// caller treats that as a hard cell failure, never as a perf signal.
func runOnce(bin string, args []string, stdinPath string, env []string) (runResult, error) {
	cmd := exec.Command(bin, args...)
	if env != nil {
		cmd.Env = env
	}
	if stdinPath != "" {
		f, err := os.Open(stdinPath)
		if err != nil {
			return runResult{}, err
		}
		defer f.Close()
		cmd.Stdin = f
	}
	var out bytes.Buffer
	cmd.Stdout = &out
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
	return runResult{Meas: m, Stdout: out.Bytes()}, err
}

// repeatRun runs runOnce reps times and reduces the samples: wall and CPU take
// the MINIMUM (the run least perturbed by scheduler/IO noise — the standard
// "how fast can it go" choice), RSS takes the MAXIMUM (the true peak), while the
// captured stdout is taken from the LAST successful run (every run produces the
// same bytes for a deterministic command, so any one is representative; the last
// keeps the code simple). reps is clamped to >= 1. The first error aborts.
func repeatRun(reps int, bin string, args []string, stdinPath string, env []string) (runResult, error) {
	if reps < 1 {
		reps = 1
	}
	var best runResult
	for i := 0; i < reps; i++ {
		r, err := runOnce(bin, args, stdinPath, env)
		if err != nil {
			return runResult{}, err
		}
		if i == 0 {
			best = r
			continue
		}
		if r.Meas.Wall < best.Meas.Wall {
			best.Meas.Wall = r.Meas.Wall
		}
		if r.Meas.CPUTotal() < best.Meas.CPUTotal() {
			best.Meas.CPUUser, best.Meas.CPUSys = r.Meas.CPUUser, r.Meas.CPUSys
		}
		if r.Meas.MaxRSSKB > best.Meas.MaxRSSKB {
			best.Meas.MaxRSSKB = r.Meas.MaxRSSKB
		}
		best.Stdout = r.Stdout
	}
	return best, nil
}
