package cliflag

import (
	"flag"
	"io"
	"reflect"
	"strings"
	"testing"
)

// newPosixFS builds a FlagSet mirroring the kind of short flags samtools view
// registers: a mix of boolean switches (-b -S -h -H) and value-taking flags
// (-q int, -o string, -@ int — note the non-alphanumeric short name).
func newPosixFS() (*flag.FlagSet, *struct {
	b, s, h, hOnly bool
	q, at          int
	o              string
}) {
	fs := flag.NewFlagSet("posixtest", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	dst := &struct {
		b, s, h, hOnly bool
		q, at          int
		o              string
	}{}
	BoolVar(fs, &dst.b, "b", "bam", false, "bool b")
	BoolVar(fs, &dst.s, "S", "sam", false, "bool S")
	BoolVar(fs, &dst.h, "h", "with-header", false, "bool h")
	BoolVar(fs, &dst.hOnly, "H", "header-only", false, "bool H")
	IntVar(fs, &dst.q, "q", "min-mapq", 0, "int q")
	IntVar(fs, &dst.at, "@", "threads", 0, "int @")
	StringVar(fs, &dst.o, "o", "output", "", "string o")
	return fs, dst
}

func TestNormalize(t *testing.T) {
	tests := []struct {
		name    string
		args    []string
		want    []string
		wantErr bool
	}{
		{name: "empty", args: []string{}, want: []string{}},
		{name: "single bool canonical", args: []string{"-b"}, want: []string{"-b"}},
		{name: "bundling two bools", args: []string{"-bS"}, want: []string{"-b", "-S"}},
		{name: "bundling three bools", args: []string{"-bSH"}, want: []string{"-b", "-S", "-H"}},
		{name: "value concat", args: []string{"-q20"}, want: []string{"-q", "20"}},
		{name: "mixed bundle then value concat", args: []string{"-bSq20"}, want: []string{"-b", "-S", "-q", "20"}},
		{name: "value flag at end takes next arg", args: []string{"-bSq", "20"}, want: []string{"-b", "-S", "-q", "20"}},
		{name: "value flag alone takes next arg", args: []string{"-q", "20"}, want: []string{"-q", "20"}},
		{name: "string value concat", args: []string{"-oout.bam"}, want: []string{"-o", "out.bam"}},
		{name: "hb cluster", args: []string{"-hb"}, want: []string{"-h", "-b"}},
		{name: "non-alnum value short @ concat", args: []string{"-@4"}, want: []string{"-@", "4"}},
		{name: "non-alnum value short @ next arg", args: []string{"-@", "4"}, want: []string{"-@", "4"}},
		{name: "double dash terminator", args: []string{"-bS", "--", "-q20", "in.bam"}, want: []string{"-b", "-S", "--", "-q20", "in.bam"}},
		{name: "bare dash stdin passthrough", args: []string{"-bS", "-"}, want: []string{"-b", "-S", "-"}},
		{name: "long flag passthrough", args: []string{"--bam"}, want: []string{"--bam"}},
		{name: "long flag equals passthrough", args: []string{"--min-mapq=20"}, want: []string{"--min-mapq=20"}},
		{name: "positional passthrough", args: []string{"-b", "in.bam"}, want: []string{"-b", "in.bam"}},
		{name: "value concat then positional", args: []string{"-q20", "in.bam"}, want: []string{"-q", "20", "in.bam"}},
		{name: "negative-looking value as next arg", args: []string{"-q", "-5"}, want: []string{"-q", "-5"}},
		{name: "unknown short char", args: []string{"-bZ"}, wantErr: true},
		{name: "unknown lone short char", args: []string{"-Z"}, wantErr: true},
		{name: "unknown char after value is consumed not flag", args: []string{"-qZ"}, want: []string{"-q", "Z"}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs, _ := newPosixFS()
			got, err := Normalize(fs, tt.args)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("Normalize(%v): want error, got %v", tt.args, got)
				}
				return
			}
			if err != nil {
				t.Fatalf("Normalize(%v): unexpected error: %v", tt.args, err)
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Fatalf("Normalize(%v) = %v, want %v", tt.args, got, tt.want)
			}
		})
	}
}

// TestNormalizeIdempotent verifies that normalizing already-canonical output a
// second time leaves it unchanged.
func TestNormalizeIdempotent(t *testing.T) {
	inputs := [][]string{
		{"-b", "-S", "-q", "20", "in.bam"},
		{"-b", "-S", "--", "-q20"},
		{"--bam", "-", "in.bam"},
		{},
	}
	for _, in := range inputs {
		fs, _ := newPosixFS()
		once, err := Normalize(fs, in)
		if err != nil {
			t.Fatalf("Normalize(%v): %v", in, err)
		}
		fs2, _ := newPosixFS()
		twice, err := Normalize(fs2, once)
		if err != nil {
			t.Fatalf("Normalize(%v) second pass: %v", once, err)
		}
		if !reflect.DeepEqual(once, twice) {
			t.Fatalf("not idempotent for %v: %v != %v", in, once, twice)
		}
	}
}

// TestParse checks that the end-to-end Parse wires normalization into
// fs.Parse and that flag destinations end up with the expected values.
func TestParse(t *testing.T) {
	t.Run("mixed bundle and value concat", func(t *testing.T) {
		fs, dst := newPosixFS()
		if err := Parse(fs, []string{"-bSq20", "in.bam"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if !dst.b || !dst.s {
			t.Fatalf("want b and S set, got b=%v S=%v", dst.b, dst.s)
		}
		if dst.q != 20 {
			t.Fatalf("want q=20, got %d", dst.q)
		}
		if got := fs.Args(); !reflect.DeepEqual(got, []string{"in.bam"}) {
			t.Fatalf("want positional [in.bam], got %v", got)
		}
	})

	t.Run("value flag takes following arg", func(t *testing.T) {
		fs, dst := newPosixFS()
		if err := Parse(fs, []string{"-bSq", "30"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if !dst.b || !dst.s || dst.q != 30 {
			t.Fatalf("got b=%v S=%v q=%d", dst.b, dst.s, dst.q)
		}
	})

	t.Run("double dash stops flag parsing", func(t *testing.T) {
		fs, dst := newPosixFS()
		if err := Parse(fs, []string{"-b", "--", "-q20"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if !dst.b {
			t.Fatalf("want b set")
		}
		if got := fs.Args(); !reflect.DeepEqual(got, []string{"-q20"}) {
			t.Fatalf("want positional [-q20], got %v", got)
		}
	})

	t.Run("bare dash is positional", func(t *testing.T) {
		fs, dst := newPosixFS()
		if err := Parse(fs, []string{"-b", "-"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if !dst.b {
			t.Fatalf("want b set")
		}
		if got := fs.Args(); !reflect.DeepEqual(got, []string{"-"}) {
			t.Fatalf("want positional [-], got %v", got)
		}
	})

	t.Run("long flag still works", func(t *testing.T) {
		fs, dst := newPosixFS()
		if err := Parse(fs, []string{"--bam", "--min-mapq=15"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if !dst.b || dst.q != 15 {
			t.Fatalf("got b=%v q=%d", dst.b, dst.q)
		}
	})

	t.Run("unknown short char errors", func(t *testing.T) {
		fs, _ := newPosixFS()
		err := Parse(fs, []string{"-bZ"})
		if err == nil {
			t.Fatalf("want error for unknown short flag")
		}
		if !strings.Contains(err.Error(), "not defined") {
			t.Fatalf("want 'not defined' error, got %v", err)
		}
	})

	t.Run("string value concat", func(t *testing.T) {
		fs, dst := newPosixFS()
		if err := Parse(fs, []string{"-oout.bam"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if dst.o != "out.bam" {
			t.Fatalf("want o=out.bam, got %q", dst.o)
		}
	})
}

// TestIsBoolFlag exercises the helper directly, including the unknown-flag
// and value-flag branches.
func TestIsBoolFlag(t *testing.T) {
	fs, _ := newPosixFS()
	if !isBoolFlag(fs, "b") {
		t.Fatalf("b should be bool")
	}
	if isBoolFlag(fs, "q") {
		t.Fatalf("q should not be bool")
	}
	if isBoolFlag(fs, "nope") {
		t.Fatalf("unknown flag should not report bool")
	}
}
