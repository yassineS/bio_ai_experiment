package bedclosest

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// runMulti runs ClosestMulti over the given A and B database strings and
// returns the output, row count, and any error.
func runMulti(t *testing.T, a string, bs []string, opts Options) (string, int, error) {
	t.Helper()
	readers := make([]io.Reader, len(bs))
	for i, b := range bs {
		readers[i] = strings.NewReader(b)
	}
	var buf bytes.Buffer
	n, err := ClosestMulti(strings.NewReader(a), readers, &buf, opts)
	return buf.String(), n, err
}

// TestClosestMultiSelection exercises the per-database (`-mdb each`) and
// combined (`-mdb all`) selection logic together with the database-label
// column.
func TestClosestMultiSelection(t *testing.T) {
	// A = chr1:80-100. DB1's closest is [20,60) at signed -21; DB2's is
	// [120,170) at +21; DB3's is [70,90) which overlaps (distance 0).
	a := "chr1\t80\t100\tq1\t1\t+\n"
	db1 := "chr1\t5\t15\td1.1\t1\t+\nchr1\t20\t60\td1.2\t2\t-\nchr1\t200\t220\td1.3\t3\t-\n"
	db2 := "chr1\t15\t35\tdb2.1\t1\t-\nchr1\t120\t170\tdb2.2\t2\t-\nchr1\t210\t230\tdb3\t3\t+\n"
	db3 := "chr1\t70\t90\td3.1\t3\t-\n"
	dbs := []string{db1, db2, db3}

	tests := []struct {
		name string
		opts Options
		want string
	}{
		{
			name: "each default index labels",
			opts: Options{},
			want: "chr1\t80\t100\tq1\t1\t+\t1\tchr1\t20\t60\td1.2\t2\t-\n" +
				"chr1\t80\t100\tq1\t1\t+\t2\tchr1\t120\t170\tdb2.2\t2\t-\n" +
				"chr1\t80\t100\tq1\t1\t+\t3\tchr1\t70\t90\td3.1\t3\t-\n",
		},
		{
			name: "each with names labels",
			opts: Options{DBLabels: []string{"a", "b", "c"}},
			want: "chr1\t80\t100\tq1\t1\t+\ta\tchr1\t20\t60\td1.2\t2\t-\n" +
				"chr1\t80\t100\tq1\t1\t+\tb\tchr1\t120\t170\tdb2.2\t2\t-\n" +
				"chr1\t80\t100\tq1\t1\t+\tc\tchr1\t70\t90\td3.1\t3\t-\n",
		},
		{
			name: "each with signed -D ref distance column",
			opts: Options{ReportDistance: true, DistanceMode: DistanceSignedRef},
			want: "chr1\t80\t100\tq1\t1\t+\t1\tchr1\t20\t60\td1.2\t2\t-\t-21\n" +
				"chr1\t80\t100\tq1\t1\t+\t2\tchr1\t120\t170\tdb2.2\t2\t-\t21\n" +
				"chr1\t80\t100\tq1\t1\t+\t3\tchr1\t70\t90\td3.1\t3\t-\t0\n",
		},
		{
			name: "all picks single overall closest (DB3 overlap) with its index",
			opts: Options{MultiDBMode: MultiDBAll},
			want: "chr1\t80\t100\tq1\t1\t+\t3\tchr1\t70\t90\td3.1\t3\t-\n",
		},
		{
			name: "all with signed -D ref distance column",
			opts: Options{MultiDBMode: MultiDBAll, ReportDistance: true, DistanceMode: DistanceSignedRef},
			want: "chr1\t80\t100\tq1\t1\t+\t3\tchr1\t70\t90\td3.1\t3\t-\t0\n",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, _, err := runMulti(t, a, dbs, tc.opts)
			if err != nil {
				t.Fatalf("ClosestMulti: %v", err)
			}
			if got != tc.want {
				t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", tc.want, got)
			}
		})
	}
}

// TestClosestMultiAllTie verifies that, in -mdb all mode, equidistant hits from
// different databases are reported as a tie and resolved per TieBreak, keeping
// each surviving hit's source-database label.
func TestClosestMultiAllTie(t *testing.T) {
	// A = chr1:100-110. DB1 has [120,130) (right, +11); DB2 has [79,90)
	// (left, -11). Both are 11 away, so -mdb all reports a cross-database tie.
	a := "chr1\t100\t110\n"
	db1 := "chr1\t120\t130\tr\n"
	db2 := "chr1\t79\t90\tl\n"
	dbs := []string{db1, db2}

	tests := []struct {
		name string
		tb   TieBreak
		want string
	}{
		{
			// Upstream reports tied hits in genomic order of B (start), so the
			// left DB2 hit [79,90) comes before the right DB1 hit [120,130),
			// each keeping its own source-database index label.
			name: "all ties report both, in B genomic order",
			tb:   TieAll,
			want: "chr1\t100\t110\t2\tchr1\t79\t90\tl\n" +
				"chr1\t100\t110\t1\tchr1\t120\t130\tr\n",
		},
		{
			name: "first keeps the genomically-first hit (DB2 left)",
			tb:   TieFirst,
			want: "chr1\t100\t110\t2\tchr1\t79\t90\tl\n",
		},
		{
			name: "last keeps the genomically-last hit (DB1 right)",
			tb:   TieLast,
			want: "chr1\t100\t110\t1\tchr1\t120\t130\tr\n",
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, _, err := runMulti(t, a, dbs, Options{MultiDBMode: MultiDBAll, TieBreak: tc.tb})
			if err != nil {
				t.Fatalf("ClosestMulti: %v", err)
			}
			if got != tc.want {
				t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", tc.want, got)
			}
		})
	}
}

// TestClosestMultiNoHit verifies the no-hit handling: a database that lacks the
// query's chromosome contributes nothing (it does NOT emit its own null row),
// and only when NO database yields any hit is a single null row emitted with a
// "." database column, matching upstream's RecordOutputMgr behaviour.
func TestClosestMultiNoHit(t *testing.T) {
	a := "chr2\t10\t20\n"
	db1 := "chr1\t0\t10\n"
	db2 := "chr2\t100\t110\n"

	// each: DB1 has no chr2 (contributes nothing); DB2 has a hit.
	gotEach, _, err := runMulti(t, a, []string{db1, db2}, Options{})
	if err != nil {
		t.Fatalf("ClosestMulti each: %v", err)
	}
	wantEach := "chr2\t10\t20\t2\tchr2\t100\t110\n"
	if gotEach != wantEach {
		t.Fatalf("each mismatch.\nwant:\n%s\ngot:\n%s", wantEach, gotEach)
	}

	// all: only DB2 has a chr2 hit, so the overall closest is from DB2.
	gotAll, _, err := runMulti(t, a, []string{db1, db2}, Options{MultiDBMode: MultiDBAll})
	if err != nil {
		t.Fatalf("ClosestMulti all: %v", err)
	}
	wantAll := "chr2\t10\t20\t2\tchr2\t100\t110\n"
	if gotAll != wantAll {
		t.Fatalf("all mismatch.\nwant:\n%s\ngot:\n%s", wantAll, gotAll)
	}

	// each with no hit at all -> single null row, labelled ".".
	gotNoneEach, _, err := runMulti(t, "chr9\t1\t2\n", []string{db1, db2}, Options{})
	if err != nil {
		t.Fatalf("ClosestMulti each none: %v", err)
	}
	wantNoneEach := "chr9\t1\t2\t.\t.\t-1\t-1\n"
	if gotNoneEach != wantNoneEach {
		t.Fatalf("each none mismatch.\nwant:\n%s\ngot:\n%s", wantNoneEach, gotNoneEach)
	}

	// all with no hit at all -> single null row, labelled ".".
	gotNone, _, err := runMulti(t, "chr9\t1\t2\n", []string{db1, db2}, Options{MultiDBMode: MultiDBAll})
	if err != nil {
		t.Fatalf("ClosestMulti all none: %v", err)
	}
	wantNone := "chr9\t1\t2\t.\t.\t-1\t-1\n"
	if gotNone != wantNone {
		t.Fatalf("all none mismatch.\nwant:\n%s\ngot:\n%s", wantNone, gotNone)
	}
}

// TestClosestMultiValidation covers the label-count and required-input checks.
func TestClosestMultiValidation(t *testing.T) {
	a := "chr1\t0\t10\n"
	db := "chr1\t0\t10\n"

	t.Run("label count mismatch", func(t *testing.T) {
		_, _, err := runMulti(t, a, []string{db, db}, Options{DBLabels: []string{"only-one"}})
		if err == nil {
			t.Fatal("expected error for mismatched label count, got nil")
		}
		if !strings.Contains(err.Error(), "must match the number of -b files") {
			t.Fatalf("unexpected error: %v", err)
		}
	})

	t.Run("no databases", func(t *testing.T) {
		var buf bytes.Buffer
		_, err := ClosestMulti(strings.NewReader(a), nil, &buf, Options{})
		if err == nil {
			t.Fatal("expected error for zero databases, got nil")
		}
	})

	t.Run("matching label count is accepted", func(t *testing.T) {
		_, _, err := runMulti(t, a, []string{db, db}, Options{DBLabels: []string{"x", "y"}})
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
	})
}
