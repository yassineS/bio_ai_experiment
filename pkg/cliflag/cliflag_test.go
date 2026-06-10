package cliflag

import (
	"flag"
	"io"
	"strings"
	"testing"
	"time"
)

// newFS returns a FlagSet that discards output and returns errors (so error
// cases can be asserted without exiting the process).
func newFS() *flag.FlagSet {
	fs := flag.NewFlagSet("test", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	return fs
}

func TestStringVar(t *testing.T) {
	tests := []struct {
		name     string
		short    string
		long     string
		def      string
		args     []string
		want     string
		wantErr  bool
		errCheck func(*flag.FlagSet) // optional: assert a name is unregistered
	}{
		{name: "short form", short: "i", long: "input", def: "", args: []string{"-i", "a.txt"}, want: "a.txt"},
		{name: "long form", short: "i", long: "input", def: "", args: []string{"--input", "b.txt"}, want: "b.txt"},
		{name: "default used", short: "i", long: "input", def: "def.txt", args: nil, want: "def.txt"},
		{name: "long overrides short (last wins)", short: "i", long: "input", def: "", args: []string{"-i", "x", "--input", "y"}, want: "y"},
		{name: "short overrides long (last wins)", short: "i", long: "input", def: "", args: []string{"--input", "y", "-i", "x"}, want: "x"},
		{name: "short only registered", short: "i", long: "", def: "", args: []string{"-i", "z"}, want: "z"},
		{name: "long only registered", short: "", long: "input", def: "", args: []string{"--input", "z"}, want: "z"},
		{name: "short-only rejects long", short: "i", long: "", def: "d", args: []string{"--input", "z"}, want: "d", wantErr: true},
		{name: "long-only rejects short", short: "", long: "input", def: "d", args: []string{"-i", "z"}, want: "d", wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := newFS()
			var p string
			StringVar(fs, &p, tt.short, tt.long, tt.def, "usage")
			err := fs.Parse(tt.args)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Parse err=%v, wantErr=%v", err, tt.wantErr)
			}
			if !tt.wantErr && p != tt.want {
				t.Errorf("got %q, want %q", p, tt.want)
			}
		})
	}
}

func TestIntVar(t *testing.T) {
	tests := []struct {
		name    string
		short   string
		long    string
		args    []string
		want    int
		wantErr bool
	}{
		{name: "short", short: "n", long: "num", args: []string{"-n", "5"}, want: 5},
		{name: "long", short: "n", long: "num", args: []string{"--num", "7"}, want: 7},
		{name: "default", short: "n", long: "num", args: nil, want: 3},
		{name: "last wins", short: "n", long: "num", args: []string{"-n", "1", "--num", "9"}, want: 9},
		{name: "invalid value", short: "n", long: "num", args: []string{"-n", "notint"}, wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := newFS()
			var p int
			IntVar(fs, &p, tt.short, tt.long, 3, "usage")
			err := fs.Parse(tt.args)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Parse err=%v, wantErr=%v", err, tt.wantErr)
			}
			if !tt.wantErr && p != tt.want {
				t.Errorf("got %d, want %d", p, tt.want)
			}
		})
	}
}

func TestInt64Var(t *testing.T) {
	tests := []struct {
		name    string
		args    []string
		want    int64
		wantErr bool
	}{
		{name: "short", args: []string{"-s", "42"}, want: 42},
		{name: "long", args: []string{"--seed", "100"}, want: 100},
		{name: "default", args: nil, want: -1},
		{name: "large value", args: []string{"--seed", "9223372036854775807"}, want: 9223372036854775807},
		{name: "last wins", args: []string{"-s", "1", "--seed", "2"}, want: 2},
		{name: "invalid", args: []string{"-s", "x"}, wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := newFS()
			var p int64
			Int64Var(fs, &p, "s", "seed", -1, "usage")
			err := fs.Parse(tt.args)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Parse err=%v, wantErr=%v", err, tt.wantErr)
			}
			if !tt.wantErr && p != tt.want {
				t.Errorf("got %d, want %d", p, tt.want)
			}
		})
	}
}

func TestUint64Var(t *testing.T) {
	tests := []struct {
		name    string
		args    []string
		want    uint64
		wantErr bool
	}{
		{name: "short", args: []string{"-c", "5"}, want: 5},
		{name: "long", args: []string{"--count", "8"}, want: 8},
		{name: "default", args: nil, want: 1},
		{name: "negative rejected", args: []string{"-c", "-1"}, wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := newFS()
			var p uint64
			Uint64Var(fs, &p, "c", "count", 1, "usage")
			err := fs.Parse(tt.args)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Parse err=%v, wantErr=%v", err, tt.wantErr)
			}
			if !tt.wantErr && p != tt.want {
				t.Errorf("got %d, want %d", p, tt.want)
			}
		})
	}
}

func TestFloat64Var(t *testing.T) {
	tests := []struct {
		name    string
		args    []string
		want    float64
		wantErr bool
	}{
		{name: "short", args: []string{"-q", "2.5"}, want: 2.5},
		{name: "long", args: []string{"--qual", "3.5"}, want: 3.5},
		{name: "default", args: nil, want: 1.0},
		{name: "last wins", args: []string{"-q", "1.1", "--qual", "9.9"}, want: 9.9},
		{name: "invalid", args: []string{"-q", "abc"}, wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := newFS()
			var p float64
			Float64Var(fs, &p, "q", "qual", 1.0, "usage")
			err := fs.Parse(tt.args)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Parse err=%v, wantErr=%v", err, tt.wantErr)
			}
			if !tt.wantErr && p != tt.want {
				t.Errorf("got %v, want %v", p, tt.want)
			}
		})
	}
}

func TestBoolVar(t *testing.T) {
	tests := []struct {
		name    string
		short   string
		long    string
		def     bool
		args    []string
		want    bool
		wantErr bool
	}{
		{name: "short", short: "f", long: "force", args: []string{"-f"}, want: true},
		{name: "long", short: "f", long: "force", args: []string{"--force"}, want: true},
		{name: "explicit false", short: "f", long: "force", args: []string{"--force=false"}, want: false},
		{name: "default true", short: "f", long: "force", def: true, args: nil, want: true},
		{name: "default false", short: "f", long: "force", args: nil, want: false},
		{name: "short only", short: "f", long: "", args: []string{"-f"}, want: true},
		{name: "long only", short: "", long: "force", args: []string{"--force"}, want: true},
		{name: "invalid value", short: "f", long: "force", args: []string{"--force=notbool"}, wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := newFS()
			var p bool
			BoolVar(fs, &p, tt.short, tt.long, tt.def, "usage")
			err := fs.Parse(tt.args)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Parse err=%v, wantErr=%v", err, tt.wantErr)
			}
			if !tt.wantErr && p != tt.want {
				t.Errorf("got %v, want %v", p, tt.want)
			}
		})
	}
}

func TestDurationVar(t *testing.T) {
	tests := []struct {
		name    string
		args    []string
		want    time.Duration
		wantErr bool
	}{
		{name: "short", args: []string{"-t", "5s"}, want: 5 * time.Second},
		{name: "long", args: []string{"--timeout", "2m"}, want: 2 * time.Minute},
		{name: "default", args: nil, want: time.Second},
		{name: "invalid", args: []string{"-t", "nope"}, wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := newFS()
			var p time.Duration
			DurationVar(fs, &p, "t", "timeout", time.Second, "usage")
			err := fs.Parse(tt.args)
			if (err != nil) != tt.wantErr {
				t.Fatalf("Parse err=%v, wantErr=%v", err, tt.wantErr)
			}
			if !tt.wantErr && p != tt.want {
				t.Errorf("got %v, want %v", p, tt.want)
			}
		})
	}
}

// stringSlice is a minimal repeatable flag.Value used to exercise Var.
type stringSlice []string

func (s *stringSlice) String() string { return strings.Join(*s, ",") }
func (s *stringSlice) Set(v string) error {
	*s = append(*s, v)
	return nil
}

func TestVar(t *testing.T) {
	t.Run("accumulates across both forms", func(t *testing.T) {
		fs := newFS()
		var got stringSlice
		Var(fs, &got, "r", "region", "usage")
		if err := fs.Parse([]string{"-r", "chr1", "--region", "chr2", "-r", "chr3"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if want := "chr1,chr2,chr3"; got.String() != want {
			t.Errorf("got %q, want %q", got.String(), want)
		}
	})
	t.Run("short only rejects long", func(t *testing.T) {
		fs := newFS()
		var got stringSlice
		Var(fs, &got, "r", "", "usage")
		if err := fs.Parse([]string{"--region", "chr1"}); err == nil {
			t.Errorf("expected error for unregistered long form")
		}
	})
	t.Run("long only", func(t *testing.T) {
		fs := newFS()
		var got stringSlice
		Var(fs, &got, "", "region", "usage")
		if err := fs.Parse([]string{"--region", "chr1"}); err != nil {
			t.Fatalf("Parse: %v", err)
		}
		if got.String() != "chr1" {
			t.Errorf("got %q, want chr1", got.String())
		}
	})
}

func TestFormatUsage(t *testing.T) {
	tests := []struct {
		name      string
		short     string
		long      string
		valueType string
		desc      string
		wantSubs  []string
	}{
		{name: "both with type", short: "i", long: "input", valueType: "FILE", desc: "Input file",
			wantSubs: []string{"-i FILE", "--input FILE", "Input file"}},
		{name: "both no type (bool)", short: "h", long: "help", valueType: "", desc: "Show help",
			wantSubs: []string{"-h", "--help", "Show help"}},
		{name: "short only", short: "v", long: "", valueType: "", desc: "Version",
			wantSubs: []string{"-v", "Version"}},
		{name: "long only", short: "", long: "version", valueType: "", desc: "Version",
			wantSubs: []string{"--version", "Version"}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			out := FormatUsage(tt.short, tt.long, tt.valueType, tt.desc)
			for _, sub := range tt.wantSubs {
				if !strings.Contains(out, sub) {
					t.Errorf("FormatUsage(%q,%q,%q,%q) = %q, missing %q", tt.short, tt.long, tt.valueType, tt.desc, out, sub)
				}
			}
		})
	}
	// long-only with type and the no-name case both exercised for full coverage.
	if out := FormatUsage("", "input", "FILE", "d"); !strings.Contains(out, "--input FILE") {
		t.Errorf("long-only with type: %q", out)
	}
	if out := FormatUsage("", "", "", "d"); !strings.Contains(out, "d") {
		t.Errorf("no-name case: %q", out)
	}
}
