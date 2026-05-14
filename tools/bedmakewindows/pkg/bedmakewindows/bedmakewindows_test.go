package bedmakewindows

import (
	"bytes"
	"strings"
	"testing"
)

func TestParseNaming(t *testing.T) {
	cases := []struct {
		in   string
		want Naming
		ok   bool
	}{
		{"", NameWinNum, true},
		{"winnum", NameWinNum, true},
		{"WINNUM", NameWinNum, true},
		{"srcwinnum", NameSrcWinNum, true},
		{"src", NameSrc, true},
		{"bogus", NoName, false},
	}
	for _, tc := range cases {
		got, err := ParseNaming(tc.in)
		if tc.ok && err != nil {
			t.Errorf("ParseNaming(%q): unexpected err %v", tc.in, err)
		}
		if !tc.ok && err == nil {
			t.Errorf("ParseNaming(%q): expected error", tc.in)
		}
		if tc.ok && got != tc.want {
			t.Errorf("ParseNaming(%q): want %v got %v", tc.in, tc.want, got)
		}
	}
}

func TestOptionsValidate(t *testing.T) {
	t.Run("missing width and count", func(t *testing.T) {
		o := Options{}
		if err := o.validate(); err == nil {
			t.Fatal("expected error")
		}
	})
	t.Run("width and count both set", func(t *testing.T) {
		o := Options{Width: 10, Count: 3}
		if err := o.validate(); err == nil {
			t.Fatal("expected error")
		}
	})
	t.Run("step without width", func(t *testing.T) {
		o := Options{Count: 3, Step: 5}
		if err := o.validate(); err == nil {
			t.Fatal("expected error")
		}
	})
	t.Run("step defaults to width", func(t *testing.T) {
		o := Options{Width: 100}
		if err := o.validate(); err != nil {
			t.Fatalf("err: %v", err)
		}
		if o.Step != 100 {
			t.Fatalf("step want 100 got %d", o.Step)
		}
	})
}

func TestFromGenome(t *testing.T) {
	in := "chr1\t100\nchr2\t200\n# comment\n\n"
	ivs, err := FromGenome(strings.NewReader(in))
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if len(ivs) != 2 {
		t.Fatalf("want 2 intervals, got %d", len(ivs))
	}
	if ivs[0].Chrom != "chr1" || ivs[0].End != 100 {
		t.Errorf("first interval: %+v", ivs[0])
	}
}

func TestFromGenome_Errors(t *testing.T) {
	if _, err := FromGenome(strings.NewReader("chr1\n")); err == nil {
		t.Errorf("expected error for malformed line")
	}
	if _, err := FromGenome(strings.NewReader("chr1\tabc\n")); err == nil {
		t.Errorf("expected error for non-numeric size")
	}
}

func TestFromBED(t *testing.T) {
	in := "chr5\t60000\t70000\tAAA\nchr5\t73000\t90000\tBBB\nchr5\t100000\t101000\tCCC\n"
	ivs, err := FromBED(strings.NewReader(in))
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if len(ivs) != 3 {
		t.Fatalf("want 3 intervals, got %d", len(ivs))
	}
	if ivs[1].Name != "BBB" {
		t.Fatalf("expected name BBB, got %q", ivs[1].Name)
	}
}

func TestFromBED_Errors(t *testing.T) {
	if _, err := FromBED(strings.NewReader("chr1\tA\t10\n")); err == nil {
		t.Errorf("expected error for bad start")
	}
	if _, err := FromBED(strings.NewReader("chr1\t0\tB\n")); err == nil {
		t.Errorf("expected error for bad end")
	}
	if _, err := FromBED(strings.NewReader("chr1\t0\n")); err == nil {
		t.Errorf("expected error for too few fields")
	}
}

func TestMakeWindows_FixedWidth_Forward(t *testing.T) {
	ivs := mustBED(t, "chr5\t60000\t70000\tAAA\nchr5\t73000\t90000\tBBB\nchr5\t100000\t101000\tCCC\n")
	var out, warn bytes.Buffer
	n, err := MakeWindows(ivs, &out, &warn, Options{Width: 5000, Naming: NameWinNum})
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "chr5\t60000\t65000\t1\n" +
		"chr5\t65000\t70000\t2\n" +
		"chr5\t73000\t78000\t1\n" +
		"chr5\t78000\t83000\t2\n" +
		"chr5\t83000\t88000\t3\n" +
		"chr5\t88000\t90000\t4\n" +
		"chr5\t100000\t101000\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\nwant=%q\ngot =%q", want, got)
	}
	if n != 7 {
		t.Fatalf("count: want 7 got %d", n)
	}
}

func TestMakeWindows_FixedWidth_Reverse(t *testing.T) {
	ivs := mustBED(t, "chr5\t60000\t70000\tAAA\nchr5\t100000\t101000\tCCC\n")
	var out bytes.Buffer
	if _, err := MakeWindows(ivs, &out, nil, Options{Width: 5000, Reverse: true, Naming: NameWinNum}); err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "chr5\t60000\t65000\t2\nchr5\t65000\t70000\t1\nchr5\t100000\t101000\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch: want=%q got=%q", want, got)
	}
}

func TestMakeWindows_StepBelowWidth(t *testing.T) {
	// Matches upstream makewindows.t03 for the BBB interval slice.
	ivs := mustBED(t, "chr5\t73000\t90000\tBBB\n")
	var out bytes.Buffer
	if _, err := MakeWindows(ivs, &out, nil, Options{Width: 5000, Step: 2000, Naming: NameWinNum}); err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "chr5\t73000\t78000\t1\n" +
		"chr5\t75000\t80000\t2\n" +
		"chr5\t77000\t82000\t3\n" +
		"chr5\t79000\t84000\t4\n" +
		"chr5\t81000\t86000\t5\n" +
		"chr5\t83000\t88000\t6\n" +
		"chr5\t85000\t90000\t7\n" +
		"chr5\t87000\t90000\t8\n" +
		"chr5\t89000\t90000\t9\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestMakeWindows_FixedCount(t *testing.T) {
	ivs := mustBED(t, "chr5\t60000\t70000\tAAA\nchr5\t73000\t90000\tBBB\nchr5\t100000\t101000\tCCC\n")
	var out bytes.Buffer
	if _, err := MakeWindows(ivs, &out, nil, Options{Count: 3, Naming: NameWinNum}); err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "chr5\t60000\t63333\t1\n" +
		"chr5\t63333\t66666\t2\n" +
		"chr5\t66666\t70000\t3\n" +
		"chr5\t73000\t78666\t1\n" +
		"chr5\t78666\t84332\t2\n" +
		"chr5\t84332\t90000\t3\n" +
		"chr5\t100000\t100333\t1\n" +
		"chr5\t100333\t100666\t2\n" +
		"chr5\t100666\t101000\t3\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestMakeWindows_FixedCount_SrcWinNum(t *testing.T) {
	ivs := mustBED(t, "1\t11\t44\tA\n")
	var out bytes.Buffer
	if _, err := MakeWindows(ivs, &out, nil, Options{Count: 10, Naming: NameSrcWinNum}); err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "1\t11\t14\tA_1\n" +
		"1\t14\t17\tA_2\n" +
		"1\t17\t20\tA_3\n" +
		"1\t20\t23\tA_4\n" +
		"1\t23\t26\tA_5\n" +
		"1\t26\t29\tA_6\n" +
		"1\t29\t32\tA_7\n" +
		"1\t32\t35\tA_8\n" +
		"1\t35\t38\tA_9\n" +
		"1\t38\t44\tA_10\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestMakeWindows_FixedCount_TooSmallWarns(t *testing.T) {
	ivs := mustBED(t, "1\t10\t19\tC\n")
	var out, warn bytes.Buffer
	n, err := MakeWindows(ivs, &out, &warn, Options{Count: 10, Naming: NameSrcWinNum})
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if n != 0 {
		t.Fatalf("want 0 windows, got %d", n)
	}
	wantWarn := "WARNING: Interval 1:10-19 is smaller than the number of windows requested. Skipping.\n"
	if got := warn.String(); got != wantWarn {
		t.Fatalf("warn mismatch:\nwant=%q\ngot =%q", wantWarn, got)
	}
}

func TestMakeWindows_NoName(t *testing.T) {
	ivs := mustBED(t, "chr1\t0\t10\tFOO\n")
	var out bytes.Buffer
	if _, err := MakeWindows(ivs, &out, nil, Options{Width: 5, Naming: NoName}); err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "chr1\t0\t5\nchr1\t5\t10\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch: want=%q got=%q", want, got)
	}
}

func TestMakeWindows_NameSrc(t *testing.T) {
	ivs := mustBED(t, "chr1\t0\t10\tFOO\n")
	var out bytes.Buffer
	if _, err := MakeWindows(ivs, &out, nil, Options{Width: 5, Naming: NameSrc}); err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "chr1\t0\t5\tFOO\nchr1\t5\t10\tFOO\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch: want=%q got=%q", want, got)
	}
}

func TestMakeWindows_EmptyInput(t *testing.T) {
	var out bytes.Buffer
	n, err := MakeWindows(nil, &out, nil, Options{Width: 100})
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if n != 0 {
		t.Fatalf("want 0, got %d", n)
	}
}

func TestMakeWindows_ValidationErr(t *testing.T) {
	if _, err := MakeWindows(nil, nil, nil, Options{}); err == nil {
		t.Fatal("expected validation error")
	}
}

func TestMakeWindows_GenomeInput(t *testing.T) {
	ivs, err := FromGenome(strings.NewReader("chr1\t10\n"))
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	var out bytes.Buffer
	if _, err := MakeWindows(ivs, &out, nil, Options{Width: 4, Naming: NameWinNum}); err != nil {
		t.Fatalf("err: %v", err)
	}
	want := "chr1\t0\t4\t1\nchr1\t4\t8\t2\nchr1\t8\t10\t3\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch: want=%q got=%q", want, got)
	}
}

func mustBED(t *testing.T, s string) []Interval {
	t.Helper()
	ivs, err := FromBED(strings.NewReader(s))
	if err != nil {
		t.Fatalf("FromBED: %v", err)
	}
	return ivs
}
