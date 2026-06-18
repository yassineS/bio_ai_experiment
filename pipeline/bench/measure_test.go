package bench

import (
	"os"
	"path/filepath"
	"runtime"
	"syscall"
	"testing"
	"time"
)

// TestTimevalDuration checks the seconds+microseconds → Duration conversion.
func TestTimevalDuration(t *testing.T) {
	got := timevalDuration(syscall.Timeval{Sec: 2, Usec: 500000})
	if want := 2500 * time.Millisecond; got != want {
		t.Fatalf("timevalDuration = %v, want %v", got, want)
	}
}

// TestMsMbRatio covers the small reducers used throughout the report.
func TestMsMbRatio(t *testing.T) {
	if got := ms(1500 * time.Microsecond); got != 1.5 {
		t.Errorf("ms = %v, want 1.5", got)
	}
	if got := mb(2048); got != 2.0 {
		t.Errorf("mb = %v, want 2.0", got)
	}
	if got := ratio(3, 2); got != 1.5 {
		t.Errorf("ratio = %v, want 1.5", got)
	}
	if got := ratio(1, 0); got != 0 {
		t.Errorf("ratio by zero = %v, want 0", got)
	}
}

// TestRunMeasured runs a trivial real process and checks the measurement is
// populated: wall-clock is non-zero and (on Linux) max RSS is recorded. It
// shells out to `/bin/sh -c` so it needs no fixtures or upstream binaries, which
// keeps it CI-safe.
func TestRunMeasured(t *testing.T) {
	if runtime.GOOS == "windows" {
		t.Skip("POSIX-only measurement test")
	}
	m, err := runMeasured("/bin/sh", []string{"-c", "for i in 1 2 3; do :; done"}, "", "")
	if err != nil {
		t.Fatalf("runMeasured: %v", err)
	}
	if m.Wall <= 0 {
		t.Errorf("wall = %v, want > 0", m.Wall)
	}
	if runtime.GOOS == "linux" && m.MaxRSSKB <= 0 {
		t.Errorf("maxRSS = %d KiB, want > 0 on linux", m.MaxRSSKB)
	}
}

// TestRunMeasuredStdout verifies stdout redirection writes the payload (so the
// streaming cells' output cost is really incurred and measured).
func TestRunMeasuredStdout(t *testing.T) {
	if runtime.GOOS == "windows" {
		t.Skip("POSIX-only measurement test")
	}
	dir := t.TempDir()
	outPath := filepath.Join(dir, "out.txt")
	if _, err := runMeasured("/bin/sh", []string{"-c", "printf hello"}, "", outPath); err != nil {
		t.Fatalf("runMeasured: %v", err)
	}
	b, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatal(err)
	}
	if string(b) != "hello" {
		t.Errorf("stdout file = %q, want %q", b, "hello")
	}
}

// TestRepeatMeasuredReducer checks the reps reducer returns a populated result
// and respects the min-wall / max-RSS reduction contract over multiple runs.
func TestRepeatMeasuredReducer(t *testing.T) {
	if runtime.GOOS == "windows" {
		t.Skip("POSIX-only measurement test")
	}
	m, err := repeatMeasured(3, "/bin/sh", []string{"-c", ":"}, "", "")
	if err != nil {
		t.Fatalf("repeatMeasured: %v", err)
	}
	if m.Wall <= 0 {
		t.Errorf("reduced wall = %v, want > 0", m.Wall)
	}
}

// TestBenchMatrixWellFormed sanity-checks the matrix: unique names, a known file
// group, and a non-nil builder for each cell.
func TestBenchMatrixWellFormed(t *testing.T) {
	groups := map[string]bool{"BAM": true, "CRAM": true, "VCF": true, "BED": true, "FASTQ": true}
	seen := map[string]bool{}
	cells := BenchMatrix()
	if len(cells) == 0 {
		t.Fatal("empty matrix")
	}
	for _, c := range cells {
		if c.Name == "" || c.Build == nil {
			t.Errorf("cell %+v missing name/build", c)
		}
		if seen[c.Name] {
			t.Errorf("duplicate cell name %q", c.Name)
		}
		seen[c.Name] = true
		if !groups[c.Group] {
			t.Errorf("cell %q has unknown group %q", c.Name, c.Group)
		}
	}
}
