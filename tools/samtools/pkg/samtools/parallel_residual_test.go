package samtools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestParallelIndexBAIDeterminism pins samtools index -@ as a perf-only
// knob: the .bai bytes produced with Threads = N must equal the serial
// Threads = 0 path exactly. BAI construction is per-record virtual-offset
// ordered, so -@ today is single-threaded on the index build itself; this
// test is the regression gate if a future patch adds a worker-pool
// feeder on the decompression side.
func TestParallelIndexBAIDeterminism(t *testing.T) {
	src := filepath.Join("..", "..", "testdata", "parity", "test_input_1_a.sam")
	if _, err := os.Stat(src); err != nil {
		t.Skipf("fixture missing: %v", err)
	}

	// SAM → sorted BAM bytes once.
	var bamBuf bytes.Buffer
	in, err := os.Open(src)
	if err != nil {
		t.Fatal(err)
	}
	defer in.Close()
	if _, err := View(in, &bamBuf, ViewOptions{OutputBAM: true}); err != nil {
		t.Fatalf("View -> BAM: %v", err)
	}
	var sortedBAM bytes.Buffer
	if err := Sort(bytes.NewReader(bamBuf.Bytes()), &sortedBAM, SortOptions{}); err != nil {
		t.Fatalf("Sort: %v", err)
	}

	// Serial index.
	var serial bytes.Buffer
	if err := Index(bytes.NewReader(sortedBAM.Bytes()), &serial, IndexOptions{Threads: 0}); err != nil {
		t.Fatalf("Index -@ 0: %v", err)
	}

	for _, threads := range []int{2, 4, 8} {
		var parallel bytes.Buffer
		if err := Index(bytes.NewReader(sortedBAM.Bytes()), &parallel, IndexOptions{Threads: threads}); err != nil {
			t.Fatalf("Index -@ %d: %v", threads, err)
		}
		if !bytes.Equal(serial.Bytes(), parallel.Bytes()) {
			t.Errorf(".bai bytes differ at -@ %d (serial=%d, parallel=%d)",
				threads, serial.Len(), parallel.Len())
		}
	}
}

// TestParallelMpileupDeterminism pins samtools mpileup -@ as a perf-only
// knob: text output with Threads = N must equal Threads = 0 byte-for-byte.
func TestParallelMpileupDeterminism(t *testing.T) {
	bam := filepath.Join("..", "..", "..", "bcftools", "testdata", "mpileup", "mpileup.1.bam")
	fa := filepath.Join("..", "..", "..", "bcftools", "testdata", "mpileup", "mpileup.ref.fa")
	for _, p := range []string{bam, fa} {
		if _, err := os.Stat(p); err != nil {
			t.Skipf("fixture missing: %v", err)
		}
	}

	var serial bytes.Buffer
	if err := MpileupFile(MpileupOptions{Inputs: []string{bam}, FastaRef: fa, Threads: 0}, &serial); err != nil {
		t.Fatalf("MpileupFile -@ 0: %v", err)
	}

	for _, threads := range []int{2, 4} {
		var parallel bytes.Buffer
		if err := MpileupFile(MpileupOptions{Inputs: []string{bam}, FastaRef: fa, Threads: threads}, &parallel); err != nil {
			t.Fatalf("MpileupFile -@ %d: %v", threads, err)
		}
		if !bytes.Equal(serial.Bytes(), parallel.Bytes()) {
			t.Errorf("mpileup output differs at -@ %d (serial=%d, parallel=%d)",
				threads, serial.Len(), parallel.Len())
		}
	}
}

// TestMpileupRejectsBCFOutput pins our port's rejection-parity for
// `samtools mpileup -g/-u`: upstream samtools 1.x REMOVED these flags
// and emits a specific message redirecting to `bcftools mpileup`. Our
// port carries the exact same wording, so this is parity, not a gap.
func TestMpileupRejectsBCFOutput(t *testing.T) {
	wantSub := `using "samtools mpileup" to generate BCF or VCF files has been removed`
	if !strings.Contains(ErrMpileupBCFNotImplemented.Error(), wantSub) {
		t.Errorf("rejection message missing upstream wording: %q",
			ErrMpileupBCFNotImplemented.Error())
	}
}
