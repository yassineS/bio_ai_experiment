package bcftools

import (
	"bufio"
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// liveMatrixDir holds the genuine-bcftools oracle: basic.vcf plus
// expected_matrix.tsv (EXPR<TAB>comma-separated POS kept by `view -i EXPR`).
const liveMatrixDir = "testdata/filter_live"

// keptPositions runs `view -i expr` over the live fixture and returns the
// comma-joined POS column of the surviving records, or an error.
func keptPositions(t *testing.T, vcfText, expr string, exclude bool) (string, error) {
	t.Helper()
	opts := ViewOptions{}
	if exclude {
		opts.ExcludeExpr = expr
	} else {
		opts.IncludeExpr = expr
	}
	var out bytes.Buffer
	if _, err := View(strings.NewReader(vcfText), &out, opts); err != nil {
		return "", err
	}
	var pos []string
	for _, line := range strings.Split(out.String(), "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) >= 2 {
			pos = append(pos, fields[1])
		}
	}
	return strings.Join(pos, ","), nil
}

// TestFilterLiveMatrix checks every expression in the committed oracle matrix
// against the kept-record set produced by genuine bcftools 1.23.1. Rows with
// an empty expected set encode the ambiguous bare-tag cases (DP defined in
// both INFO and FORMAT), where upstream aborts with an error instead of
// emitting records; we assert that we likewise reject the expression.
func TestFilterLiveMatrix(t *testing.T) {
	vcfBytes, err := os.ReadFile(filepath.Join(liveMatrixDir, "basic.vcf"))
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	vcfText := string(vcfBytes)

	f, err := os.Open(filepath.Join(liveMatrixDir, "expected_matrix.tsv"))
	if err != nil {
		t.Fatalf("open matrix: %v", err)
	}
	defer f.Close()

	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := sc.Text()
		if strings.TrimSpace(line) == "" {
			continue
		}
		parts := strings.SplitN(line, "\t", 2)
		expr := parts[0]
		var want string
		if len(parts) == 2 {
			want = strings.TrimRight(strings.TrimSpace(parts[1]), ",")
		}

		t.Run(expr, func(t *testing.T) {
			got, err := keptPositions(t, vcfText, expr, false)
			if want == "" {
				// Ambiguous bare-tag expression: upstream errors out, so a
				// non-nil error (or an empty kept set) matches the oracle.
				if err == nil && got != "" {
					t.Fatalf("include %q: expected upstream error / no records, got %q", expr, got)
				}
				return
			}
			if err != nil {
				t.Fatalf("include %q: unexpected error: %v", expr, err)
			}
			if got != want {
				t.Fatalf("include %q: kept %q, want %q", expr, got, want)
			}
		})
	}
	if err := sc.Err(); err != nil {
		t.Fatalf("scan matrix: %v", err)
	}
}

// TestFilterLiveExcludeComplement verifies that, for the well-defined
// (non-ambiguous) expressions, `-e EXPR` keeps exactly the records that
// `-i EXPR` drops — i.e. exclude is the precise complement of include over
// the same input.
func TestFilterLiveExcludeComplement(t *testing.T) {
	vcfBytes, err := os.ReadFile(filepath.Join(liveMatrixDir, "basic.vcf"))
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	vcfText := string(vcfBytes)

	// All records in the fixture, in file order.
	allPos := []string{"100", "200", "300", "400", "50", "150"}

	exprs := []string{
		`QUAL<30`,
		`AC>1`,
		`AF>0.2`,
		`TYPE="snp"`,
		`AN=6`,
		`AC[0]>1`,
		`N_ALT=1`,
		`INFO/DP>50`,
		`POS>200`,
		`CHROM="chr1"`,
		`REF="A"`,
		`GT="het"`,
		`GT="hom"`,
		`FORMAT/DP>10`,
		`N_PASS(GT="alt")>0`,
		`STRLEN(REF)=1`,
	}
	for _, expr := range exprs {
		t.Run(expr, func(t *testing.T) {
			inc, err := keptPositions(t, vcfText, expr, false)
			if err != nil {
				t.Fatalf("include %q: %v", expr, err)
			}
			exc, err := keptPositions(t, vcfText, expr, true)
			if err != nil {
				t.Fatalf("exclude %q: %v", expr, err)
			}
			incSet := map[string]bool{}
			for _, p := range strings.Split(inc, ",") {
				if p != "" {
					incSet[p] = true
				}
			}
			// The complement, in file order, must equal the exclude result.
			var wantExc []string
			for _, p := range allPos {
				if !incSet[p] {
					wantExc = append(wantExc, p)
				}
			}
			if exc != strings.Join(wantExc, ",") {
				t.Fatalf("exclude %q: got %q, want complement %q (include kept %q)",
					expr, exc, strings.Join(wantExc, ","), inc)
			}
		})
	}
}
