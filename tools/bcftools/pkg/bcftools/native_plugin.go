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
	if s, ok := plugin.(stderrSink); ok {
		s.SetStderr(stderr)
	}
	if s, ok := plugin.(outputSuppressor); ok {
		s.SetStdout(out)
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

	// Stage 2: run Process over every record, preserving input order.
	results, err := processRecords(plugin, variants, opts.Threads)
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
