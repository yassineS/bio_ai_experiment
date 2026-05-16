package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// TestSize_TableDriven covers Size on a few hand-built FASTA/FASTQ
// inputs. Output format is "<n>\t<total_bases>\n" matching upstream.
func TestSize_TableDriven(t *testing.T) {
	cases := []struct {
		name  string
		input string
		want  string
	}{
		{
			name:  "fasta three records",
			input: ">a\nACGT\n>b\nGGG\n>c\nT\n",
			want:  "3\t8\n",
		},
		{
			name:  "fasta one record",
			input: ">one\nACGTACGTAC\n",
			want:  "1\t10\n",
		},
		{
			name:  "fastq three records",
			input: "@r1\nACGT\n+\nIIII\n@r2\nGG\n+\n##\n@r3\nTTTT\n+\n!!!!\n",
			want:  "3\t10\n",
		},
		{
			name:  "empty input fasta",
			input: "",
			want:  "0\t0\n",
		},
		{
			name:  "fasta wrapped multi-line sequence collapses to total length",
			input: ">long\nACGT\nACGT\nACGT\n",
			want:  "1\t12\n",
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var out bytes.Buffer
			if err := Size(strings.NewReader(tc.input), &out); err != nil {
				t.Fatalf("Size: %v", err)
			}
			if got := out.String(); got != tc.want {
				t.Fatalf("Size output: want %q, got %q", tc.want, got)
			}
		})
	}
}
