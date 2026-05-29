package main

import (
	"flag"
	"io"
	"os"
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
)

// newNormFlagSet builds a FlagSet mirroring the value-taking and boolean
// short flags that `bcftools norm`/`view` register, so the normalizer
// tests exercise the same flag surface as the real subcommands.
func newNormFlagSet() *flag.FlagSet {
	fs := flag.NewFlagSet("test", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		multiallelics string
		rmDup         string
		outputType    string
		outputPath    string
		fastaRef      string
		atomize       bool
		showHelp      bool
	)
	cliflag.StringVar(fs, &multiallelics, "m", "multiallelics", "", "")
	cliflag.StringVar(fs, &rmDup, "d", "rm-dup", "", "")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "")
	cliflag.StringVar(fs, &fastaRef, "f", "fasta-ref", "", "")
	cliflag.BoolVar(fs, &atomize, "a", "atomize", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	registerNoVersionIfAbsent(fs)
	return fs
}

// TestNormalizeShortFlags verifies that getopt-style attached short-flag
// values (`-Xvalue`) are rewritten into the two-token form Go's flag
// package accepts, while every form that must be preserved is left
// untouched.
func TestNormalizeShortFlags(t *testing.T) {
	tests := []struct {
		name string
		in   []string
		want []string
	}{
		{"attached value-flag bundle", []string{"-mboth"}, []string{"-m", "both"}},
		{"attached minus", []string{"-m-"}, []string{"-m", "-"}},
		{"attached plus", []string{"-m+"}, []string{"-m", "+"}},
		{"output type attached", []string{"-Ob"}, []string{"-O", "b"}},
		{"already split", []string{"-m", "-"}, []string{"-m", "-"}},
		{"equals form preserved", []string{"-m=both"}, []string{"-m=both"}},
		{"long flag untouched", []string{"--multiallelics", "-"}, []string{"--multiallelics", "-"}},
		{"long equals untouched", []string{"--output-type=b"}, []string{"--output-type=b"}},
		{"bool short untouched", []string{"-a"}, []string{"-a"}},
		{"no-version long untouched", []string{"--no-version"}, []string{"--no-version"}},
		{"bare dash is stdin", []string{"-"}, []string{"-"}},
		{"end-of-options marker", []string{"--", "-mboth"}, []string{"--", "-mboth"}},
		{"value after dashdash untouched", []string{"-m+", "--", "-Ob"}, []string{"-m", "+", "--", "-Ob"}},
		{
			"realistic mix",
			[]string{"-m-", "--no-version", "-Ob", "-o", "out.bcf", "file.vcf"},
			[]string{"-m", "-", "--no-version", "-O", "b", "-o", "out.bcf", "file.vcf"},
		},
	}
	fs := newNormFlagSet()
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := normalizeShortFlags(fs, tt.in)
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("normalizeShortFlags(%v) = %v, want %v", tt.in, got, tt.want)
			}
		})
	}
}

// TestParseFlagsAttachedShortValues confirms that parsing the various
// upstream-accepted argument forms through the shared parseFlags helper
// populates the destination variables identically, regardless of whether
// the value was attached (`-Ob`), split (`-O b`), or long-equals
// (`--output-type=b`), and that `--no-version` is accepted in any
// position.
func TestParseFlagsAttachedShortValues(t *testing.T) {
	tests := []struct {
		name           string
		args           []string
		wantOutputType string
		wantMulti      string
		wantNoVersion  bool
		wantRest       []string
	}{
		{"attached -Ob", []string{"-Ob", "in.vcf"}, "b", "", false, []string{"in.vcf"}},
		{"split -O b", []string{"-O", "b", "in.vcf"}, "b", "", false, []string{"in.vcf"}},
		{"long equals", []string{"--output-type=b", "in.vcf"}, "b", "", false, []string{"in.vcf"}},
		{"attached -m-", []string{"-m-", "in.vcf"}, "v", "-", false, []string{"in.vcf"}},
		{"attached -m+", []string{"-m+", "in.vcf"}, "v", "+", false, []string{"in.vcf"}},
		{"attached -mboth", []string{"-mboth", "in.vcf"}, "v", "both", false, []string{"in.vcf"}},
		{"no-version then -Ob", []string{"--no-version", "-Ob", "in.vcf"}, "b", "", true, []string{"in.vcf"}},
		{"-Ob then no-version", []string{"-Ob", "--no-version", "in.vcf"}, "b", "", true, []string{"in.vcf"}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fs := flag.NewFlagSet("test", flag.ContinueOnError)
			fs.SetOutput(io.Discard)
			var outputType, multi string
			cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "")
			cliflag.StringVar(fs, &multi, "m", "multiallelics", "", "")
			registerNoVersionIfAbsent(fs)
			if err := parseFlags(fs, tt.args); err != nil {
				t.Fatalf("parseFlags(%v) returned error: %v", tt.args, err)
			}
			if outputType != tt.wantOutputType {
				t.Errorf("output-type = %q, want %q", outputType, tt.wantOutputType)
			}
			if multi != tt.wantMulti {
				t.Errorf("multiallelics = %q, want %q", multi, tt.wantMulti)
			}
			nv := fs.Lookup("no-version")
			if nv == nil {
				t.Fatal("no-version flag not registered")
			}
			if got := nv.Value.String() == "true"; got != tt.wantNoVersion {
				t.Errorf("no-version = %v, want %v", got, tt.wantNoVersion)
			}
			if !reflect.DeepEqual(fs.Args(), tt.wantRest) {
				t.Errorf("positional args = %v, want %v", fs.Args(), tt.wantRest)
			}
		})
	}
}

// TestRunNormAttachedMultiallelic is an end-to-end check that the norm
// subcommand accepts the attached forms `-m-`/`-m+` (split / join) and
// runs to completion (rc 0). `-mboth` must be rejected after parsing,
// matching upstream bcftools which errors "Expected '+' or '-' with -m".
func TestRunNormAttachedMultiallelic(t *testing.T) {
	in := "../../testdata/parity/basic.vcf"
	if rc := runNorm([]string{"-m-", "-o", t.TempDir() + "/split.vcf", in}); rc != 0 {
		t.Errorf("runNorm -m-: rc=%d, want 0", rc)
	}
	if rc := runNorm([]string{"-m+", "-o", t.TempDir() + "/join.vcf", in}); rc != 0 {
		t.Errorf("runNorm -m+: rc=%d, want 0", rc)
	}
	// -mboth parses fine (split into -m both) but the norm validator
	// rejects a body without a leading +/-, exactly like upstream.
	if rc := runNorm([]string{"-mboth", in}); rc == 0 {
		t.Errorf("runNorm -mboth: rc=0, want non-zero (upstream rejects)")
	}
}

// TestRunViewNoVersionAttachedOutput is an end-to-end check that
// `view --no-version -Ob` (and the reversed order) parse without error
// and emit a non-empty BCF stream.
func TestRunViewNoVersionAttachedOutput(t *testing.T) {
	in := "../../testdata/parity/basic.vcf"
	for _, args := range [][]string{
		{"--no-version", "-Ob", "-o", t.TempDir() + "/a.bcf", in},
		{"-Ob", "--no-version", "-o", t.TempDir() + "/b.bcf", in},
	} {
		out := args[3]
		if rc := runView(args); rc != 0 {
			t.Fatalf("runView %v: rc=%d, want 0", args, rc)
		}
		fi, err := os.Stat(out)
		if err != nil {
			t.Fatalf("stat %s: %v", out, err)
		}
		if fi.Size() == 0 {
			t.Errorf("runView %v produced empty output", args)
		}
	}
}
