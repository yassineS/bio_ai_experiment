// Plugin support for bcftools. Unlike upstream bcftools — which loads
// plugins as native shared objects via dlopen — this Go port runs a plugin
// as an ordinary child process and streams variant data over its standard
// streams. A plugin is therefore "a filter from VCF on stdin to VCF on
// stdout", which lets users write plugins in any language. The full
// contract is documented in docs/PLUGIN_PROTOCOL.md.
package bcftools

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// pluginEnvVar is the environment variable that holds a colon-separated
// list of directories searched for plugin executables, mirroring upstream
// bcftools' BCFTOOLS_PLUGINS.
const pluginEnvVar = "BCFTOOLS_PLUGINS"

// PluginOptions configures a single `bcftools plugin` / `bcftools +<name>`
// invocation.
type PluginOptions struct {
	// Name is the plugin name as given on the command line (the executable
	// is looked up by this exact name in the BCFTOOLS_PLUGINS directories).
	Name string
	// Args are the plugin-specific arguments passed verbatim as the child
	// process argv (everything after the `--` separator, or after the
	// recognised host options).
	Args []string
	// InputFile is the host input VCF/BCF path ("-" or "" means stdin).
	InputFile string
	// Regions optionally restricts the records fed to the plugin.
	Regions []string
	// RegionsFile is a BED-like regions file (merged into Regions).
	RegionsFile string
	// OutputFormat selects the container the host writes around the
	// plugin's VCF output (-O v|z|b|u).
	OutputFormat OutputFormat
	// CompressLevel forwards the -l gzip level for OutputVCFGz.
	CompressLevel int
	// Threads is upstream's -@/--threads. When greater than 1 it enables
	// parallel BGZF compression of -O z and -O b output via bgzf.MultiWriter;
	// the framed result decodes byte-identically regardless of thread count.
	Threads int
}

// PluginNotFoundError is returned when a plugin name cannot be resolved to
// an executable in any BCFTOOLS_PLUGINS directory.
type PluginNotFoundError struct {
	Name string
	Dirs []string
}

// Error implements the error interface.
func (e *PluginNotFoundError) Error() string {
	if len(e.Dirs) == 0 {
		return fmt.Sprintf("plugin %q not found: %s is not set", e.Name, pluginEnvVar)
	}
	return fmt.Sprintf("plugin %q not found in %s (%s)", e.Name, pluginEnvVar, strings.Join(e.Dirs, ":"))
}

// PluginExecError wraps a non-zero plugin exit. The plugin's own stderr is
// streamed straight to the host's stderr while it runs, so this error value
// only needs to carry the plugin name and the underlying exec error.
type PluginExecError struct {
	Name string
	Err  error
}

// Error implements the error interface.
func (e *PluginExecError) Error() string {
	return fmt.Sprintf("plugin %q failed: %v", e.Name, e.Err)
}

// Unwrap exposes the underlying exec error.
func (e *PluginExecError) Unwrap() error { return e.Err }

// PluginDirs returns the plugin search directories from BCFTOOLS_PLUGINS,
// in order, with empty entries dropped.
func PluginDirs() []string {
	raw := os.Getenv(pluginEnvVar)
	if raw == "" {
		return nil
	}
	var dirs []string
	for _, d := range strings.Split(raw, string(os.PathListSeparator)) {
		if d != "" {
			dirs = append(dirs, d)
		}
	}
	return dirs
}

// isExecutableFile reports whether path is a regular file with at least one
// executable bit set.
func isExecutableFile(path string) bool {
	info, err := os.Stat(path)
	if err != nil || info.IsDir() {
		return false
	}
	return info.Mode().Perm()&0o111 != 0
}

// ResolvePlugin looks up a plugin name in the BCFTOOLS_PLUGINS directories
// and returns the absolute path of the first executable file matching the
// name. The plugin executable must be named exactly `<name>` (an optional
// leading "+" on the requested name is stripped first).
func ResolvePlugin(name string) (string, error) {
	name = strings.TrimPrefix(name, "+")
	dirs := PluginDirs()
	if name == "" {
		return "", &PluginNotFoundError{Name: name, Dirs: dirs}
	}
	// An explicit path (containing a separator) bypasses the search path.
	if strings.ContainsRune(name, os.PathSeparator) {
		if isExecutableFile(name) {
			abs, err := filepath.Abs(name)
			if err != nil {
				return "", err
			}
			return abs, nil
		}
		return "", &PluginNotFoundError{Name: name, Dirs: dirs}
	}
	for _, dir := range dirs {
		cand := filepath.Join(dir, name)
		if isExecutableFile(cand) {
			abs, err := filepath.Abs(cand)
			if err != nil {
				return "", err
			}
			return abs, nil
		}
	}
	return "", &PluginNotFoundError{Name: name, Dirs: dirs}
}

// PluginInfo describes one discoverable plugin.
type PluginInfo struct {
	// Name is the plugin name (the executable file name).
	Name string
	// Path is the absolute path of the plugin executable.
	Path string
	// About is the first line of the plugin's `--about` output, populated
	// only by ListPlugins when verbose probing is requested. It is empty
	// when the plugin does not support the optional `--about` flag. For
	// native plugins it is always populated from the plugin's About() method.
	About string
	// Native reports whether this is an in-process native plugin (dispatched
	// by the native registry) rather than an external executable.
	Native bool
}

// ListPlugins scans the BCFTOOLS_PLUGINS directories and returns the
// discoverable executable plugins, sorted by name and de-duplicated (the
// first directory wins, matching the resolution order of ResolvePlugin).
// When verbose is true each plugin is additionally probed by running it
// with `--about`; probing is best-effort and never fails the listing.
func ListPlugins(verbose bool) []PluginInfo {
	seen := make(map[string]bool)
	var plugins []PluginInfo
	// Native plugins are always listed first; an exec-discoverable plugin of
	// the same name does not shadow the native one (the native dispatch wins
	// in RunPlugin).
	for _, info := range nativePluginInfos() {
		seen[info.Name] = true
		plugins = append(plugins, info)
	}
	for _, dir := range PluginDirs() {
		entries, err := os.ReadDir(dir)
		if err != nil {
			continue
		}
		for _, e := range entries {
			name := e.Name()
			if seen[name] {
				continue
			}
			path := filepath.Join(dir, name)
			if !isExecutableFile(path) {
				continue
			}
			abs, err := filepath.Abs(path)
			if err != nil {
				abs = path
			}
			seen[name] = true
			info := PluginInfo{Name: name, Path: abs}
			if verbose {
				info.About = probeAbout(abs)
			}
			plugins = append(plugins, info)
		}
	}
	sort.Slice(plugins, func(i, j int) bool { return plugins[i].Name < plugins[j].Name })
	return plugins
}

// probeAbout runs `<plugin> --about` and returns its first non-empty output
// line. Probing is optional in the protocol: a plugin that does not handle
// `--about` simply yields an empty string and is still listed.
func probeAbout(path string) string {
	cmd := exec.Command(path, "--about")
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = io.Discard
	if err := cmd.Run(); err != nil {
		return ""
	}
	for _, line := range strings.Split(out.String(), "\n") {
		if s := strings.TrimSpace(line); s != "" {
			return s
		}
	}
	return ""
}

// FormatPluginList renders the output of ListPlugins for `bcftools plugin -l`.
func FormatPluginList(plugins []PluginInfo, verbose bool) string {
	if len(plugins) == 0 {
		return ""
	}
	var b strings.Builder
	for _, p := range plugins {
		if verbose {
			// Native plugins have no on-disk path; show the name in its place
			// so the verbose listing still has a stable first column.
			if p.Native {
				b.WriteString(p.Name)
			} else {
				b.WriteString(p.Path)
			}
			if p.About != "" {
				b.WriteString("\t")
				b.WriteString(p.About)
			}
			b.WriteString("\n")
		} else {
			b.WriteString(p.Name)
			b.WriteString("\n")
		}
	}
	return b.String()
}

// RunPlugin executes a plugin against the host input and writes the result
// to out. The host normalises any input (VCF/BCF, optionally region-sliced)
// to uncompressed VCF text, feeds it to the plugin's stdin, reads the
// plugin's stdout as VCF text, and re-emits it in the requested -O format.
// The plugin's stderr is forwarded to stderr verbatim. A non-zero plugin
// exit is reported as a *PluginExecError.
func RunPlugin(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	// Native dispatch takes precedence: if the requested name resolves to a
	// registered in-process plugin, run the native pipeline. Otherwise fall
	// back to the unchanged exec path below for user-supplied executables.
	nativeName := strings.TrimPrefix(opts.Name, "+")
	if ctor, ok := nativeRegistry[nativeName]; ok {
		return runNativePlugin(ctor, opts, out, stderr)
	}

	path, err := ResolvePlugin(opts.Name)
	if err != nil {
		return err
	}

	regions := append([]string{}, opts.Regions...)
	if opts.RegionsFile != "" {
		regs, rerr := LoadRegionsFile(opts.RegionsFile)
		if rerr != nil {
			return rerr
		}
		regions = append(regions, regs...)
	}

	// Stage 1: normalise the host input to uncompressed VCF text. ViewFile
	// transparently handles VCF/VCF.gz/BCF and applies region selection.
	var vcfText bytes.Buffer
	viewOpts := ViewOptions{OutputFormat: OutputVCF, Regions: regions}
	input := opts.InputFile
	if input == "" {
		input = "-"
	}
	if _, err := ViewFile(input, &vcfText, viewOpts, stderr); err != nil {
		return fmt.Errorf("reading plugin input: %w", err)
	}

	// Stage 2: run the plugin as a child process, piping VCF text through.
	cmd := exec.Command(path, opts.Args...)
	cmd.Stdin = bytes.NewReader(vcfText.Bytes())
	cmd.Stderr = stderr
	var pluginOut bytes.Buffer
	cmd.Stdout = &pluginOut
	if err := cmd.Run(); err != nil {
		// Both a non-zero exit (*exec.ExitError) and a failure to start the
		// child surface the same way to the user.
		return &PluginExecError{Name: opts.Name, Err: err}
	}

	// Stage 3: re-emit the plugin's VCF stdout in the requested container.
	return writePluginOutput(pluginOut.Bytes(), out, opts)
}

// writePluginOutput parses the plugin's VCF stdout and writes it to out
// using the host's -O output formatting.
func writePluginOutput(vcfBytes []byte, out io.Writer, opts PluginOptions) error {
	r := vcf.NewReader(bytes.NewReader(vcfBytes))
	hdr, err := r.ReadHeader()
	if err != nil {
		return fmt.Errorf("plugin %q produced no valid VCF header: %w", opts.Name, err)
	}
	variants, err := r.ReadAll()
	if err != nil {
		return fmt.Errorf("plugin %q produced malformed VCF: %w", opts.Name, err)
	}
	w, cleanup, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, hdr)
	if err != nil {
		return err
	}
	defer cleanup()
	if err := w.WriteHeader(); err != nil {
		return err
	}
	for _, v := range variants {
		if err := w.Write(v); err != nil {
			return err
		}
	}
	return w.Flush()
}
