// Built-in plugin implementations. Upstream bcftools ships these as
// dlopen-loaded shared objects; this port re-implements the common ones
// natively in Go so their output byte-matches upstream without a
// compiled .so or a subprocess. RunPlugin dispatches to builtinPlugins
// before falling back to the BCFTOOLS_PLUGINS executable search.
package bcftools

import (
	"fmt"
	"io"
	"strings"
)

// builtinPluginFunc runs a native plugin against opts (its Name is the
// plugin name, Args are the verbatim plugin argv, InputFile is the host
// input) and writes the result to out. Diagnostics go to stderr.
type builtinPluginFunc func(opts PluginOptions, out io.Writer, stderr io.Writer) error

// builtinPlugins maps a plugin name to its native Go implementation.
var builtinPlugins = map[string]builtinPluginFunc{
	"mendelian2": runBuiltinMendelian2,
	"fill-tags":  runBuiltinFillTags,
}

// pluginArgs is a tiny getopt-style splitter for the verbatim plugin
// argv. It collects the input positional plus a flat list of
// (flag, value) options. It understands the conventions used by the
// bcftools plugins we re-implement: single-dash short flags and
// double-dash long flags, each optionally taking the next token as its
// value. The first bare (non-dash) token is taken as the input file;
// any further bare tokens are returned in extras (regions).
type pluginArgs struct {
	input  string
	extras []string
	opts   map[string]string
	flags  map[string]bool
}

// parsePluginArgs splits argv given the set of flags that consume a
// value (valued) and the set that are pure booleans (bools). Both
// short and long forms must be listed (e.g. "-P" and "--ped").
func parsePluginArgs(argv []string, valued, bools map[string]bool) (pluginArgs, error) {
	pa := pluginArgs{opts: map[string]string{}, flags: map[string]bool{}}
	for i := 0; i < len(argv); i++ {
		tok := argv[i]
		if tok == "-" || len(tok) == 0 || tok[0] != '-' {
			if pa.input == "" {
				pa.input = tok
			} else {
				pa.extras = append(pa.extras, tok)
			}
			continue
		}
		// Support --flag=value.
		if eq := strings.IndexByte(tok, '='); eq > 0 {
			name := tok[:eq]
			pa.opts[name] = tok[eq+1:]
			continue
		}
		if bools[tok] {
			pa.flags[tok] = true
			continue
		}
		if valued[tok] {
			if i+1 >= len(argv) {
				return pa, fmt.Errorf("flag %s requires a value", tok)
			}
			pa.opts[tok] = argv[i+1]
			i++
			continue
		}
		return pa, fmt.Errorf("unrecognised plugin option %q", tok)
	}
	return pa, nil
}

// optVal returns the value for any of the given flag aliases, or "".
func (pa pluginArgs) optVal(aliases ...string) string {
	for _, a := range aliases {
		if v, ok := pa.opts[a]; ok {
			return v
		}
	}
	return ""
}

// runBuiltinMendelian2 is the native `+mendelian2` implementation. It
// parses the upstream plugin flags (-p/--pfm, -P/--ped, -m/--mode,
// -i/--include, -e/--exclude, -o/--output, -O/--output-type) and the
// input positional, then delegates to Mendelian2File.
func runBuiltinMendelian2(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	valued := map[string]bool{
		"-p": true, "--pfm": true,
		"-P": true, "--ped": true,
		"-m": true, "--mode": true,
		"-i": true, "--include": true,
		"-e": true, "--exclude": true,
		"-o": true, "--output": true,
		"-O": true, "--output-type": true,
		"-r": true, "--regions": true,
		"-R": true, "--regions-file": true,
		"-t": true, "--targets": true,
		"-T": true, "--targets-file": true,
	}
	pa, err := parsePluginArgs(opts.Args, valued, nil)
	if err != nil {
		return fmt.Errorf("bcftools +mendelian2: %w", err)
	}
	input := pa.input
	if input == "" {
		input = opts.InputFile
	}
	if input == "" {
		input = "-"
	}

	mopts := Mendelian2Options{}
	if pfm := pa.optVal("-p", "--pfm"); pfm != "" {
		parsed, perr := ParseMendelian2PFM(pfm)
		if perr != nil {
			return perr
		}
		mopts.PFM = &parsed
	}
	mopts.PEDFile = pa.optVal("-P", "--ped")
	if mode := pa.optVal("-m", "--mode"); mode != "" {
		m, merr := ParseMendelian2Mode(mode)
		if merr != nil {
			return merr
		}
		mopts.Mode = m
	}
	mopts.IncludeExpr = pa.optVal("-i", "--include")
	mopts.ExcludeExpr = pa.optVal("-e", "--exclude")
	if ot := pa.optVal("-O", "--output-type"); ot != "" {
		f, ferr := ParseOutputFormat(ot)
		if ferr != nil {
			return ferr
		}
		mopts.OutputFormat = f
	}

	_, err = Mendelian2File(input, out, mopts)
	return err
}
