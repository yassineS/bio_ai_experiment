// Native (in-process, pure-Go) plugin framework for bcftools.
//
// Upstream bcftools ships its default plugins as compiled C shared objects
// loaded via dlopen through the htslib C API. Each plugin exposes the
// lifecycle init(argc,argv,header) -> process(record) -> destroy. This file
// re-implements that contract natively in Go so the most common plugins run
// without spawning a child process, dispatched directly by `+name`. The
// existing exec-based path (RunPlugin in plugin.go) remains as a fallback for
// user-supplied executables that are not part of the native registry.
//
// The native pipeline mirrors the exec path's three stages: stage 1 reads and
// optionally region-slices the host input to VCF records (ViewFile), stage 2
// runs the plugin lifecycle over those records, and stage 3 re-emits them in
// the requested -O container (openOutput). Per-record plugins are processed
// through a worker pool sized by PluginOptions.Threads with strictly ordered
// output, so the result is byte-identical regardless of the thread count.
package bcftools

import (
	"bytes"
	"fmt"
	"io"
	"runtime"
	"sort"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// NativePlugin is the in-process contract modelled on upstream bcftools'
// plugin lifecycle. A plugin instance is created fresh per invocation by its
// registry constructor, then driven through Init (once, before streaming),
// Process (once per record), and Destroy (once, after the last record).
type NativePlugin interface {
	// Name returns the plugin name as dispatched by `+name`.
	Name() string
	// About returns the one-line description shown by `bcftools plugin -l -v`,
	// matching the upstream plugin's about() string.
	About() string
	// Init parses the plugin's own arguments (everything after `--`) and is
	// given the input header. It returns the output header, which may be the
	// same header with added/removed ##INFO or ##FORMAT lines. Init runs once
	// before any record is processed.
	Init(args []string, hdr *vcf.Header) (*vcf.Header, error)
	// Process transforms a single record. It returns zero records (to drop the
	// record), one record (the common case), or multiple records. The input
	// variant may be mutated in place and returned. Process must be safe to
	// call concurrently from multiple goroutines for parallel plugins; plugins
	// that are not concurrency-safe must report Parallel()==false.
	Process(v *vcf.Variant) ([]*vcf.Variant, error)
	// Destroy releases any resources and performs end-of-stream side effects
	// (for example, counts prints its totals here). It runs once after the
	// last record. Destroy may write to the host stderr passed via SetStderr.
	Destroy() error
}

// parallelPlugin is implemented by native plugins that can have Process called
// concurrently across records. Plugins that need cross-record state or
// multiple passes must NOT implement it (or return false), so they run
// serially. All batch-1 plugins are per-record and parallel.
type parallelPlugin interface {
	// Parallel reports whether Process may be invoked concurrently for
	// different records. When false the pipeline processes records serially.
	Parallel() bool
}

// bufferedPlugin is implemented by stateful plugins that need to see the whole
// record stream at once — for example to annotate each site with the distance
// to its neighbours (look-ahead and look-back). When a plugin implements it,
// the pipeline hands the complete, in-order variant slice to ProcessAll
// instead of calling Process per record, and emits the returned slice. Such a
// plugin is inherently serial; it must not also report Parallel()==true.
type bufferedPlugin interface {
	// ProcessAll transforms the entire ordered record stream in one pass and
	// returns the records to emit, in order. Init has already run and Destroy
	// runs afterwards, exactly as for the per-record path.
	ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error)
}

// runStylePlugin is implemented by native plugins whose upstream counterpart
// exports a `run` symbol rather than the init/process/destroy lifecycle.
// Upstream dispatches such plugins by handing the plugin its full argv and
// letting the plugin's own getopt parse it; the plugin options therefore
// precede the trailing input-file positional and there is no `--` separator
// (e.g. `bcftools +variant-distance -d nearest FILE`). The generic
// init/process plugins, by contrast, require the `+name [host-opts] [FILE] --
// [plugin-opts]` form, because options before the file are parsed as host
// options. The host CLI uses RunStyle to pick the matching argument-splitting
// rule, and FlagTakesValue to know which of the plugin's own flags consume the
// following token while splitting the input file out of the plugin options.
type runStylePlugin interface {
	// RunStyle reports whether the plugin uses the upstream run() dispatch,
	// accepting its options before the input file with no `--` separator.
	RunStyle() bool
	// FlagTakesValue reports whether the given plugin flag (e.g. "-d" or
	// "--direction") consumes the following CLI token as its value. It is used
	// by the host only to separate the lone input-file positional from the
	// plugin options in the run-style form.
	FlagTakesValue(flag string) bool
}

// IsRunStyleNativePlugin reports whether name resolves to a native plugin that
// uses the upstream run() dispatch (plugin options before the input file, no
// `--` separator). The leading "+" of the `+name` shorthand must be stripped
// by the caller. It returns false for non-native names and for generic
// init/process plugins.
func IsRunStyleNativePlugin(name string) bool {
	ctor, ok := nativeRegistry[name]
	if !ok {
		return false
	}
	rs, ok := ctor().(runStylePlugin)
	return ok && rs.RunStyle()
}

// NativePluginFlagTakesValue reports whether the named run-style native
// plugin's flag consumes the following CLI token as its value. It returns
// false for non-run-style or unknown plugins. The host uses it to split the
// input-file positional out of a run-style plugin's options.
func NativePluginFlagTakesValue(name, flag string) bool {
	ctor, ok := nativeRegistry[name]
	if !ok {
		return false
	}
	rs, ok := ctor().(runStylePlugin)
	if !ok || !rs.RunStyle() {
		return false
	}
	return rs.FlagTakesValue(flag)
}

// fullPlugin is implemented by native plugins that take over the entire run
// rather than fitting the read-input / per-record / re-emit pipeline. Such a
// plugin reuses an existing engine that owns both input parsing and output
// writing (for example +mendelian2, which calls the shared Mendelian2 engine
// that can emit a text count summary and/or the filtered VCF/BCF itself). When
// a registered plugin implements fullPlugin, runNativePlugin delegates the
// whole invocation to RunFull and does nothing else.
type fullPlugin interface {
	// RunFull executes the plugin end to end, reading opts.InputFile and
	// writing the result to out (stderr receives any diagnostics). It is
	// responsible for option parsing, input reading and output formatting.
	RunFull(opts PluginOptions, out io.Writer, stderr io.Writer) error
}

// multiOutputPlugin is implemented by native plugins that write multiple output
// files themselves (e.g. +split writes one VCF/BCF per sample, +scatter writes
// one per chunk/region, +variantkey-hex writes its index files). Unlike the
// single-output pipeline, such a plugin owns ALL of its writers: it creates the
// output directory and every per-file handle, choosing filenames to match
// upstream exactly. The host's single `out` writer is used only for any textual
// end-of-run report (variantkey-hex prints its counts to stdout); plugins that
// produce no stdout text leave it untouched. RunMulti receives the full
// PluginOptions (it reads opts.InputFile and parses its own -o/-O directory and
// container options from opts.Args) so the plugin can reproduce upstream's
// file-naming and per-file contents byte-for-byte.
type multiOutputPlugin interface {
	// RunMulti executes the plugin end to end, reading opts.InputFile and writing
	// its own set of output files. Any textual report goes to out; diagnostics go
	// to stderr.
	RunMulti(opts PluginOptions, out io.Writer, stderr io.Writer) error
}

// stderrSink is implemented by plugins that emit end-of-run diagnostics
// (e.g. "Filled N alleles") to stderr. The pipeline wires the host stderr in
// before Destroy is called.
type stderrSink interface {
	// SetStderr provides the writer the plugin should use for diagnostics.
	SetStderr(w io.Writer)
}

// outputSuppressor is implemented by plugins (such as counts) whose upstream
// init() returns 1 to suppress the VCF/BCF output entirely; they produce only
// their own textual report on stdout via SetStdout/Destroy. When SuppressVCF
// returns true the pipeline does not open the -O writer or emit any records.
type outputSuppressor interface {
	// SuppressVCF reports whether the VCF/BCF output should be suppressed.
	SuppressVCF() bool
	// SetStdout provides the host stdout writer the plugin reports to.
	SetStdout(w io.Writer)
}

// cmdLineSink is implemented by plugins (such as smpl-stats, indel-stats and
// guess-ploidy) whose upstream report embeds a verbatim command-line line
// ("CMD\t<name> <opts...> <file>", or the "# The command line was:" banner).
// The pipeline hands the plugin the upstream-equivalent argv so it can
// reproduce that line byte-for-byte. The argv mirrors upstream's run()-style
// dispatch: argv[0] is the plugin name, followed by the plugin options, with
// the input file as the trailing positional.
type cmdLineSink interface {
	// SetArgv provides the upstream-equivalent argv (name, plugin options,
	// then the input file).
	SetArgv(argv []string)
}

// nativeRegistry maps a plugin name to a constructor that returns a fresh
// instance. It is populated by init() functions in each native plugin file.
var nativeRegistry = map[string]func() NativePlugin{}

// registerNativePlugin adds a native plugin constructor to the registry. It is
// called from the init() of each plugin source file. It panics on a duplicate
// name so a registration clash is caught at startup rather than silently
// shadowing a plugin.
func registerNativePlugin(name string, ctor func() NativePlugin) {
	if _, dup := nativeRegistry[name]; dup {
		panic(fmt.Sprintf("bcftools: native plugin %q registered twice", name))
	}
	nativeRegistry[name] = ctor
}

// nativePluginNames returns the registered native plugin names, sorted.
func nativePluginNames() []string {
	names := make([]string, 0, len(nativeRegistry))
	for n := range nativeRegistry {
		names = append(names, n)
	}
	sort.Strings(names)
	return names
}

// isNativePlugin reports whether name resolves to a native plugin. The leading
// "+" of the `+name` shorthand is not present here (callers strip it).
func isNativePlugin(name string) bool {
	_, ok := nativeRegistry[name]
	return ok
}

// nativePluginInfos returns PluginInfo entries for the native plugins, with
// the About line filled from each plugin's About() method. These are merged
// with exec-discoverable plugins by ListPlugins.
func nativePluginInfos() []PluginInfo {
	names := nativePluginNames()
	infos := make([]PluginInfo, 0, len(names))
	for _, n := range names {
		p := nativeRegistry[n]()
		infos = append(infos, PluginInfo{Name: n, About: p.About(), Native: true})
	}
	return infos
}

// runNativePlugin executes a registered native plugin against the host input
// and writes the result in the requested -O format. It mirrors RunPlugin's
// stage structure: read+slice input (ViewFile), run the lifecycle, re-emit
// (openOutput). The plugin name in opts.Name has any leading "+" already
// stripped by the caller, but TrimPrefix is applied again defensively.
func runNativePlugin(ctor func() NativePlugin, opts PluginOptions, out io.Writer, stderr io.Writer) error {
	plugin := ctor()
	// A fullPlugin owns its entire invocation (input reading and output
	// writing), bypassing the read/process/re-emit stages below.
	if fp, ok := plugin.(fullPlugin); ok {
		return fp.RunFull(opts, out, stderr)
	}
	// A multiOutputPlugin owns its own per-file writers (one output file per
	// sample / chunk / region), so it too bypasses the single-output pipeline.
	if mp, ok := plugin.(multiOutputPlugin); ok {
		return mp.RunMulti(opts, out, stderr)
	}
	if s, ok := plugin.(stderrSink); ok {
		s.SetStderr(stderr)
	}
	if s, ok := plugin.(outputSuppressor); ok {
		s.SetStdout(out)
	}
	if s, ok := plugin.(cmdLineSink); ok {
		// Reconstruct the upstream run()-style argv: name, plugin options, then
		// the input file as the trailing positional (matching pluginCLIArgs and
		// how upstream's run() observes its own argv).
		name := opts.Name
		if name == "" {
			name = plugin.Name()
		}
		argv := append([]string{name}, opts.Args...)
		if opts.InputFile != "" && opts.InputFile != "-" {
			argv = append(argv, opts.InputFile)
		}
		s.SetArgv(argv)
	}

	regions := append([]string{}, opts.Regions...)
	if opts.RegionsFile != "" {
		regs, rerr := LoadRegionsFile(opts.RegionsFile)
		if rerr != nil {
			return rerr
		}
		regions = append(regions, regs...)
	}

	// Stage 1: normalise the host input to uncompressed VCF text, then parse
	// header + records. ViewFile transparently handles VCF/VCF.gz/BCF and
	// applies region selection, exactly as the exec path does.
	var vcfText bytes.Buffer
	input := opts.InputFile
	if input == "" {
		input = "-"
	}
	if _, err := ViewFile(input, &vcfText, ViewOptions{OutputFormat: OutputVCF, Regions: regions}, stderr); err != nil {
		return fmt.Errorf("reading plugin input: %w", err)
	}

	r := vcf.NewReader(bytes.NewReader(vcfText.Bytes()))
	hdr, err := r.ReadHeader()
	if err != nil {
		return fmt.Errorf("plugin %q: %w", opts.Name, err)
	}
	variants, err := r.ReadAll()
	if err != nil {
		return fmt.Errorf("plugin %q: malformed VCF input: %w", opts.Name, err)
	}

	// Lifecycle: Init runs once on the header before any record.
	outHdr, err := plugin.Init(opts.Args, hdr)
	if err != nil {
		return &PluginExecError{Name: opts.Name, Err: err}
	}
	if outHdr == nil {
		outHdr = hdr
	}

	// Stage 2: transform the records. A bufferedPlugin sees the whole ordered
	// stream at once (for cross-record look-ahead/look-back); every other
	// plugin runs through the per-record worker pool, preserving input order.
	var results []*vcf.Variant
	if bp, ok := plugin.(bufferedPlugin); ok {
		results, err = bp.ProcessAll(variants)
	} else {
		results, err = processRecords(plugin, variants, opts.Threads)
	}
	if err != nil {
		_ = plugin.Destroy()
		return &PluginExecError{Name: opts.Name, Err: err}
	}

	// A plugin may suppress the VCF/BCF output entirely (upstream init()==1),
	// producing only its own report. In that case skip stage 3.
	if s, ok := plugin.(outputSuppressor); ok && s.SuppressVCF() {
		return plugin.Destroy()
	}

	// Stage 3: re-emit in the requested -O container.
	w, cleanup, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, outHdr)
	if err != nil {
		_ = plugin.Destroy()
		return err
	}
	if err := w.WriteHeader(); err != nil {
		cleanup()
		_ = plugin.Destroy()
		return err
	}
	for _, v := range results {
		if v == nil {
			continue
		}
		if err := w.Write(v); err != nil {
			cleanup()
			_ = plugin.Destroy()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		cleanup()
		_ = plugin.Destroy()
		return err
	}
	cleanup()

	// Destroy runs after the output is fully flushed so its stderr summary
	// (e.g. counts totals) appears after all records were emitted.
	return plugin.Destroy()
}

// processRecords runs plugin.Process over every record. When the plugin is
// parallel and threads > 1, records are dispatched to a worker pool but the
// results are reassembled into strict input order, so the emitted stream is
// byte-identical to the single-threaded run. When the plugin is serial (or
// threads <= 1), records are processed in order on the caller's goroutine.
func processRecords(plugin NativePlugin, variants []*vcf.Variant, threads int) ([]*vcf.Variant, error) {
	parallel := true
	if p, ok := plugin.(parallelPlugin); ok {
		parallel = p.Parallel()
	}

	if !parallel || threads <= 1 || len(variants) < 2 {
		out := make([]*vcf.Variant, 0, len(variants))
		for _, v := range variants {
			res, err := plugin.Process(v)
			if err != nil {
				return nil, err
			}
			out = append(out, res...)
		}
		return out, nil
	}

	nworkers := threads
	if max := runtime.NumCPU(); nworkers > max {
		nworkers = max
	}
	if nworkers > len(variants) {
		nworkers = len(variants)
	}

	// Per-record output slots keep the global order; each record may expand to
	// zero or more output records, flattened in order at the end.
	slots := make([][]*vcf.Variant, len(variants))
	var (
		mu       sync.Mutex
		firstErr error
	)
	idx := make(chan int)
	var wg sync.WaitGroup
	for w := 0; w < nworkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range idx {
				res, err := plugin.Process(variants[i])
				if err != nil {
					mu.Lock()
					if firstErr == nil {
						firstErr = err
					}
					mu.Unlock()
					continue
				}
				slots[i] = res
			}
		}()
	}
	for i := range variants {
		idx <- i
	}
	close(idx)
	wg.Wait()

	if firstErr != nil {
		return nil, firstErr
	}
	out := make([]*vcf.Variant, 0, len(variants))
	for _, s := range slots {
		out = append(out, s...)
	}
	return out, nil
}
