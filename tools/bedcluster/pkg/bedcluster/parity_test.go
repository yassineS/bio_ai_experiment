package bedcluster

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Parity tests run the canonical upstream bedtools cluster test cases
// (reference_code/bedtools/test/cluster/test-cluster.sh) through this Go
// port and compare byte-for-byte against the recorded expected output.

const expT1 = `chr1	72017	884436	a	1	+	1
chr1	72017	844113	b	2	+	1
chr1	939517	1011278	c	3	+	2
chr1	1142976	1203168	d	4	+	3
chr1	1153667	1298845	e	5	-	3
chr1	1153667	1219633	f	6	+	3
chr1	1155173	1200334	g	7	-	3
chr1	1229798	1500664	h	8	-	3
chr1	1297735	1357056	i	9	+	3
chr1	1844181	1931789	j	10	-	4
`

const expT2 = `chr1	72017	884436	a	1	+	1
chr1	72017	844113	b	2	+	1
chr1	939517	1011278	c	3	+	2
chr1	1142976	1203168	d	4	+	3
chr1	1153667	1219633	f	6	+	3
chr1	1297735	1357056	i	9	+	4
chr1	1153667	1298845	e	5	-	5
chr1	1155173	1200334	g	7	-	5
chr1	1229798	1500664	h	8	-	5
chr1	1844181	1931789	j	10	-	6
`

func loadParityInput(t *testing.T) []byte {
	t.Helper()
	candidates := []string{
		filepath.Join("..", "..", "testdata", "parity", "in.bed"),
		filepath.Join("..", "..", "..", "..", "reference_code", "bedtools", "test", "cluster", "in.bed"),
	}
	for _, p := range candidates {
		data, err := os.ReadFile(p)
		if err == nil {
			return data
		}
	}
	t.Skip("upstream cluster test data not available")
	return nil
}

func TestParity_Basic(t *testing.T) {
	in := loadParityInput(t)
	var out bytes.Buffer
	if _, err := Cluster(bytes.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	if got := out.String(); got != expT1 {
		t.Errorf("basic cluster mismatch.\ngot:\n%s\nwant:\n%s", got, expT1)
	}
}

func TestParity_Stranded(t *testing.T) {
	in := loadParityInput(t)
	var out bytes.Buffer
	if _, err := Cluster(bytes.NewReader(in), &out, Options{StrandSpec: true}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	if got := out.String(); got != expT2 {
		t.Errorf("stranded cluster mismatch.\ngot:\n%s\nwant:\n%s", got, expT2)
	}
}

// TestParity_Idempotent re-runs cluster on its own output and confirms the
// cluster IDs from the second run match the first (after stripping the
// trailing column). This isn't an upstream test but it's a useful smoke-check.
func TestParity_Idempotent(t *testing.T) {
	in := loadParityInput(t)
	var first bytes.Buffer
	if _, err := Cluster(bytes.NewReader(in), &first, Options{}); err != nil {
		t.Fatalf("Cluster #1: %v", err)
	}
	// Strip the trailing column to feed back through.
	stripped := strings.Builder{}
	for _, line := range strings.Split(strings.TrimRight(first.String(), "\n"), "\n") {
		fields := strings.Split(line, "\t")
		stripped.WriteString(strings.Join(fields[:len(fields)-1], "\t"))
		stripped.WriteByte('\n')
	}
	var second bytes.Buffer
	if _, err := Cluster(strings.NewReader(stripped.String()), &second, Options{}); err != nil {
		t.Fatalf("Cluster #2: %v", err)
	}
	if first.String() != second.String() {
		t.Errorf("clustering not idempotent")
	}
}
