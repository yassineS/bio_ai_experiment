package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

func TestCompressHomopolymers_Table(t *testing.T) {
	tests := []struct {
		name string
		in   string
		want string
	}{
		{"empty", "", ""},
		{"single base", "A", "A"},
		{"no runs", "ACGTACGT", "ACGTACGT"},
		{"single run", "AAAACGT", "ACGT"},
		{"multiple runs", "AAACCGGGT", "ACGT"},
		{"all same", "AAAAAA", "A"},
		{"alternating", "ACACAC", "ACACAC"},
		{"mixed case treated as distinct", "AAaaCC", "AaC"},
		{"trailing run", "ACGTTTTT", "ACGT"},
		{"interior and edges", "TTTACGGGAAAAT", "TACGAT"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := CompressHomopolymers([]byte(tt.in))
			if string(got) != tt.want {
				t.Errorf("CompressHomopolymers(%q) = %q, want %q", tt.in, string(got), tt.want)
			}
		})
	}
}

func TestHPC_BasicFasta(t *testing.T) {
	in := ">s\nAAACCGT\n"
	want := ">s\nACGT\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}

func TestHPC_MultipleRecords(t *testing.T) {
	in := ">a\nAAA\n>b desc\nCCGGTT\n>c\nA\n"
	want := ">a\nA\n>b desc\nCGT\n>c\nA\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}

func TestHPC_PreservesNameAndDescription(t *testing.T) {
	in := ">seq1 length=10\nAAAACCCC\n"
	want := ">seq1 length=10\nAC\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}

func TestHPC_EmptySequenceSkipped(t *testing.T) {
	// Empty sequence -> no record emitted, matching upstream seqtk hpc.
	in := ">empty\n\n>nonempty\nAAA\n"
	want := ">nonempty\nA\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}

func TestHPC_MultiLineSequenceIsConcatenated(t *testing.T) {
	// The FASTA reader concatenates wrapped sequence lines before HPC sees
	// them; runs that span line boundaries must still collapse correctly.
	in := ">s\nAAAA\nAACC\nCC\n"
	want := ">s\nAC\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}

func TestHPC_FastqInput(t *testing.T) {
	in := "@r1\nAAACCGT\n+\nIIIIIII\n"
	want := ">r1\nACGT\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}

func TestHPC_FastqMultipleRecords(t *testing.T) {
	in := "@r1\nAAACCG\n+\nIIIIII\n@r2\nTTTTTAAA\n+\nIIIIIIII\n"
	want := ">r1\nACG\n>r2\nTA\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}

func TestHPC_PreservesCaseAtRunStart(t *testing.T) {
	// First byte of each run is kept; case at run start is preserved.
	in := ">s\naaAACC\n"
	want := ">s\naAC\n"
	var buf bytes.Buffer
	if err := HPC(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("HPC: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("HPC output = %q, want %q", got, want)
	}
}
