// Native port of the upstream `split` plugin (plugins/split.c). It splits a VCF
// by sample, writing one output file per sample (default) or per -S sample-set /
// -G group, into the directory given by -o. Output filenames match upstream
// exactly: the sample (or set) name with the characters from the set
// [ \t:/\] replaced by "_", a numeric "-N" suffix appended on clashes, and the
// container suffix (.vcf, .vcf.gz, .bcf) chosen from -O.
//
// This is a multiOutputPlugin: it owns all of its per-file writers. The default
// per-sample split, the -S samples-file and -G groups-file modes, the -k keep
// tags selection and the -O container/level options are supported byte-for-byte.
// The -i/-e per-output filter expressions need the bcftools filter engine and
// the index-jump -r/-R/-t/-T/-W options are not part of the streaming pipeline;
// these are reported as a clean unsupported Init error rather than producing
// divergent output.
package bcftools

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("split", func() NativePlugin { return &splitPlugin{} })
}

// splitSubset is one output file: the kept sample indices, optional renamed
// sample names, and the chosen base filename (without container suffix).
type splitSubset struct {
	smpl   []int
	rename []string
	fname  string
}

// splitPlugin implements split.
type splitPlugin struct {
	outputDir   string
	samplesFile string
	groupsFile  string
	keepTags    string
	format      OutputFormat
	clevel      int
	args        []string
}

// Name returns the plugin name.
func (p *splitPlugin) Name() string { return "split" }

// About returns the one-line description, matching split.c about().
func (p *splitPlugin) About() string {
	return "Split VCF by sample, creating single- or multi-sample VCFs\n"
}

// RunStyle reports that split is a run()-style plugin: upstream's split.c
// exports a `run` symbol, so it owns its entire argv (including -o/-O) before
// the trailing input-file positional, with no `--` separator
// (e.g. `bcftools +split -o DIR FILE`).
func (p *splitPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of split's flags consumes the following
// CLI token as its value, so the host can separate the lone input-file
// positional from the plugin options.
func (p *splitPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-o", "--output", "-O", "--output-type",
		"-S", "--samples-file", "-G", "--groups-file",
		"-k", "--keep-tags", "-v", "--verbosity", "--hts-opts":
		return true
	}
	return false
}

// Init parses and validates the plugin arguments. It rejects the modes that need
// the filter engine or index-jump machinery.
func (p *splitPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.format = OutputVCF
	p.clevel = -1
	p.args = args
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("split: option %q requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-o", "--output":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.outputDir = v
		case "-O", "--output-type":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if err := p.parseOutputType(v); err != nil {
				return nil, err
			}
		case "-S", "--samples-file":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.samplesFile = v
		case "-G", "--groups-file":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.groupsFile = v
		case "-k", "--keep-tags":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.keepTags = v
		case "-i", "--include", "-e", "--exclude":
			return nil, fmt.Errorf("split: the -i/-e per-output filter expressions require the bcftools filter engine and are not supported in the native plugin; run upstream bcftools for that")
		case "-r", "--regions", "-R", "--regions-file", "-t", "--targets", "-T", "--targets-file":
			return nil, fmt.Errorf("split: region/target selection is not supported in the native split plugin; pre-filter with bcftools view")
		case "-W", "--write-index":
			return nil, fmt.Errorf("split: -W/--write-index is not supported in the native plugin")
		case "-v", "--verbosity":
			if _, err := next(); err != nil {
				return nil, err
			}
		case "--hts-opts":
			if _, err := next(); err != nil {
				return nil, err
			}
		default:
			if strings.HasPrefix(a, "-O") && len(a) > 2 {
				if err := p.parseOutputType(a[2:]); err != nil {
					return nil, err
				}
				continue
			}
			if strings.HasPrefix(a, "-o") && len(a) > 2 {
				p.outputDir = a[2:]
				continue
			}
			if strings.HasPrefix(a, "-W") {
				return nil, fmt.Errorf("split: -W/--write-index is not supported in the native plugin")
			}
			return nil, fmt.Errorf("split: unsupported option %q", a)
		}
	}
	if p.samplesFile != "" && p.groupsFile != "" {
		return nil, fmt.Errorf("split: only one of -S/--samples-file or -G/--groups-file can be given")
	}
	if p.outputDir == "" {
		return nil, fmt.Errorf("split: missing the -o option")
	}
	return hdr, nil
}

// parseOutputType parses upstream's -O u|b|v|z[0-9] spelling.
func (p *splitPlugin) parseOutputType(s string) error {
	if s == "" {
		return fmt.Errorf("split: empty -O argument")
	}
	switch s[0] {
	case 'b':
		p.format = OutputBCF
	case 'u':
		p.format = OutputBCFUncompressed
	case 'z':
		p.format = OutputVCFGz
	case 'v':
		p.format = OutputVCF
	default:
		n, err := strconv.Atoi(s)
		if err != nil || n < 0 || n > 9 {
			return fmt.Errorf("split: the output type %q not recognised", s)
		}
		p.clevel = n
		return nil
	}
	if len(s) > 1 {
		n, err := strconv.Atoi(s[1:])
		if err != nil || n < 0 || n > 9 {
			return fmt.Errorf("split: could not parse compression level %q", s[1:])
		}
		p.clevel = n
	}
	return nil
}

// Process is unused: split is a multiOutputPlugin.
func (p *splitPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *splitPlugin) Destroy() error { return nil }

// RunMulti reads the input, builds the subset list, and writes one file per
// subset.
func (p *splitPlugin) RunMulti(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	if _, err := p.Init(opts.Args, nil); err != nil {
		return err
	}
	hdr, variants, err := readPluginInput(opts, stderr)
	if err != nil {
		return err
	}
	if len(hdr.Samples) == 0 {
		return fmt.Errorf("split: no samples to split")
	}
	uniq := &uniqueNames{}
	subsets, err := p.buildSubsets(hdr, uniq, stderr)
	if err != nil {
		return err
	}
	if err := os.MkdirAll(p.outputDir, 0o777); err != nil {
		return err
	}

	keep := p.parseKeepTags()
	for _, set := range subsets {
		if err := p.writeSubset(hdr, variants, set, keep); err != nil {
			return err
		}
	}
	return nil
}

// splitKeep records which INFO and FORMAT tags survive the -k keep-tags spec.
// keepAllInfo / keepAllFmt mean no pruning for that category; otherwise only the
// IDs in infoTags / fmtTags are kept.
type splitKeep struct {
	keepAllInfo bool
	keepAllFmt  bool
	infoTags    map[string]bool
	fmtTags     map[string]bool
}

// pruneInfo reports whether an INFO ID is kept.
func (k splitKeep) pruneInfo(id string) bool {
	if k.keepAllInfo {
		return true
	}
	return k.infoTags[id]
}

// pruneFmt reports whether a FORMAT ID is kept.
func (k splitKeep) pruneFmt(id string) bool {
	if k.keepAllFmt {
		return true
	}
	return k.fmtTags[id]
}

// buildSubsets constructs the list of output files per the default / -S / -G
// mode, mirroring init_subsets.
func (p *splitPlugin) buildSubsets(hdr *vcf.Header, uniq *uniqueNames, stderr io.Writer) ([]splitSubset, error) {
	sampleIdx := map[string]int{}
	for i, s := range hdr.Samples {
		sampleIdx[s] = i
	}
	switch {
	case p.samplesFile == "" && p.groupsFile == "":
		subsets := make([]splitSubset, len(hdr.Samples))
		for i, s := range hdr.Samples {
			subsets[i] = splitSubset{smpl: []int{i}, fname: uniq.make(s)}
		}
		return subsets, nil
	case p.samplesFile != "":
		return p.buildSamplesFileSubsets(hdr, sampleIdx, uniq, stderr)
	default:
		return p.buildGroupsFileSubsets(hdr, sampleIdx, uniq, stderr)
	}
}

// buildSamplesFileSubsets implements the -S samples-file mode: one output per
// line, with up to three whitespace-separated columns (sample list, optional
// rename list, optional output base name).
func (p *splitPlugin) buildSamplesFileSubsets(hdr *vcf.Header, sampleIdx map[string]int, uniq *uniqueNames, stderr io.Writer) ([]splitSubset, error) {
	lines, err := readNonEmptyLines(p.samplesFile)
	if err != nil {
		return nil, err
	}
	var subsets []splitSubset
	for _, line := range lines {
		cols := strings.Fields(line)
		if len(cols) == 0 {
			continue
		}
		set := splitSubset{}
		for _, name := range strings.Split(cols[0], ",") {
			if idx, ok := sampleIdx[name]; ok {
				set.smpl = append(set.smpl, idx)
			} else {
				fmt.Fprintf(stderr, "Warning: The sample \"%s\" is not present in %s\n", name, p.samplesFile)
			}
		}
		if len(set.smpl) == 0 {
			continue
		}
		if len(cols) > 1 && cols[1] != "-" {
			set.rename = strings.Split(cols[1], ",")
			set.fname = uniq.make(set.rename[0])
		}
		if len(cols) > 2 {
			set.fname = uniq.make(cols[2])
		}
		if set.fname == "" {
			set.fname = uniq.make(hdr.Samples[set.smpl[0]])
		}
		subsets = append(subsets, set)
	}
	return subsets, nil
}

// buildGroupsFileSubsets implements the -G groups-file mode: each line assigns a
// sample (with optional rename) to one or more comma-separated output files.
func (p *splitPlugin) buildGroupsFileSubsets(hdr *vcf.Header, sampleIdx map[string]int, uniq *uniqueNames, stderr io.Writer) ([]splitSubset, error) {
	lines, err := readNonEmptyLines(p.groupsFile)
	if err != nil {
		return nil, err
	}
	var subsets []splitSubset
	byFile := map[string]int{}
	for _, line := range lines {
		cols := strings.Fields(line)
		if len(cols) == 0 {
			continue
		}
		sample := cols[0]
		idx, ok := sampleIdx[sample]
		if !ok {
			fmt.Fprintf(stderr, "Warning: The sample \"%s\" is not present in %s\n", sample, p.groupsFile)
			continue
		}
		rename := sample
		var files string
		switch {
		case len(cols) >= 3:
			if cols[1] != "-" {
				rename = cols[1]
			}
			files = cols[2]
		case len(cols) == 2:
			files = cols[1]
		default:
			files = sample
		}
		for _, f := range strings.Split(files, ",") {
			si, seen := byFile[f]
			if !seen {
				subsets = append(subsets, splitSubset{fname: uniq.make(f)})
				si = len(subsets) - 1
				byFile[f] = si
			}
			subsets[si].smpl = append(subsets[si].smpl, idx)
			subsets[si].rename = append(subsets[si].rename, rename)
		}
	}
	return subsets, nil
}

// writeSubset writes one output file for the given subset.
func (p *splitPlugin) writeSubset(hdr *vcf.Header, variants []*vcf.Variant, set splitSubset, keep splitKeep) error {
	path := filepath.Join(p.outputDir, set.fname+p.suffix(set.fname))
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("split: cannot write to %q: %w", path, err)
	}
	defer f.Close()

	outHdr := p.subsetHeader(hdr, set, keep)
	w, cleanup, err := openOutput(f, ViewOptions{OutputFormat: p.format, CompressLevel: p.clevel}, outHdr)
	if err != nil {
		return err
	}
	if err := w.WriteHeader(); err != nil {
		cleanup()
		return err
	}
	for _, v := range variants {
		rec := p.subsetRecord(v, set, keep)
		if err := w.Write(rec); err != nil {
			cleanup()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		cleanup()
		return err
	}
	cleanup()
	return nil
}

// suffix returns the container suffix for fname, suppressing it when fname
// already ends in a recognised extension (matching upstream).
func (p *splitPlugin) suffix(fname string) string {
	lower := strings.ToLower(fname)
	for _, ext := range []string{".bcf", ".vcf", ".vcf.gz", ".vcf.bgz"} {
		if strings.HasSuffix(lower, ext) {
			return ""
		}
	}
	switch p.format {
	case OutputBCF, OutputBCFUncompressed:
		return ".bcf"
	case OutputVCFGz:
		return ".vcf.gz"
	default:
		return ".vcf"
	}
}

// subsetHeader builds the output header: the renamed sample list, with INFO and
// FORMAT lines optionally pruned by -k.
func (p *splitPlugin) subsetHeader(hdr *vcf.Header, set splitSubset, keep splitKeep) *vcf.Header {
	samples := make([]string, len(set.smpl))
	for j := range set.smpl {
		if set.rename != nil && j < len(set.rename) {
			samples[j] = set.rename[j]
		} else {
			samples[j] = hdr.Samples[set.smpl[j]]
		}
	}
	out := &vcf.Header{Samples: samples}
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##INFO=") && !keep.pruneInfo(headerID(m)) {
			continue
		}
		if strings.HasPrefix(m, "##FORMAT=") && !keep.pruneFmt(headerID(m)) {
			continue
		}
		out.MetaInfo = append(out.MetaInfo, m)
	}
	return out
}

// subsetRecord produces a copy of v restricted to the subset's samples, with the
// INFO and FORMAT fields pruned per -k.
func (p *splitPlugin) subsetRecord(v *vcf.Variant, set splitSubset, keep splitKeep) *vcf.Variant {
	rec := *v
	rec.Samples = make([]vcf.Sample, len(set.smpl))
	for j, si := range set.smpl {
		rec.Samples[j] = v.Samples[si]
	}
	if !keep.keepAllInfo {
		rec.Info = map[string]string{}
		rec.InfoOrder = nil
		for _, k := range v.InfoOrder {
			if keep.pruneInfo(k) {
				rec.Info[k] = v.Info[k]
				rec.InfoOrder = append(rec.InfoOrder, k)
			}
		}
	}
	if !keep.keepAllFmt {
		var newFmt []string
		drop := map[int]bool{}
		for i, tag := range v.Format {
			if keep.pruneFmt(tag) {
				newFmt = append(newFmt, tag)
			} else {
				drop[i] = true
			}
		}
		rec.Format = newFmt
		newSamples := make([]vcf.Sample, len(rec.Samples))
		for si := range rec.Samples {
			data := map[string]string{}
			for _, tag := range newFmt {
				if val, ok := rec.Samples[si].Data[tag]; ok {
					data[tag] = val
				}
			}
			newSamples[si] = vcf.Sample{Name: rec.Samples[si].Name, Data: data}
		}
		rec.Samples = newSamples
	}
	return &rec
}

// parseKeepTags resolves the -k keep-tags spec into a splitKeep, porting the
// keep-tag logic of init_data. The spec is a comma list whose items select INFO
// or FORMAT tags; a bare INFO/FMT keeps that whole category; INFO/x, FMT/x and
// FORMAT/x switch the active category for the following bare names. When no -k is
// given (or it selects nothing concrete) both categories are kept in full.
func (p *splitPlugin) parseKeepTags() splitKeep {
	k := splitKeep{infoTags: map[string]bool{}, fmtTags: map[string]bool{}}
	if p.keepTags == "" {
		k.keepAllInfo, k.keepAllFmt = true, true
		return k
	}
	keepInfo, keepFmt := false, false
	ninfo, nfmt := 0, 0
	isInfo := false
	beg := p.keepTags
	for beg != "" {
		switch {
		case strings.HasPrefix(strings.ToUpper(beg), "INFO/"):
			isInfo = true
			beg = beg[5:]
			continue
		case strings.EqualFold(beg, "INFO"):
			keepInfo = true
			beg = ""
			continue
		case strings.HasPrefix(strings.ToUpper(beg), "INFO,"):
			keepInfo = true
			beg = beg[5:]
			continue
		case strings.HasPrefix(strings.ToUpper(beg), "FMT/"):
			isInfo = false
			beg = beg[4:]
			continue
		case strings.HasPrefix(strings.ToUpper(beg), "FORMAT/"):
			isInfo = false
			beg = beg[7:]
			continue
		case strings.EqualFold(beg, "FMT"), strings.EqualFold(beg, "FORMAT"):
			keepFmt = true
			beg = ""
			continue
		case strings.HasPrefix(strings.ToUpper(beg), "FMT,"):
			keepFmt = true
			beg = beg[4:]
			continue
		case strings.HasPrefix(strings.ToUpper(beg), "FORMAT,"):
			keepFmt = true
			beg = beg[7:]
			continue
		}
		// A bare tag name in the current category.
		name := beg
		rest := ""
		if comma := strings.IndexByte(beg, ','); comma >= 0 {
			name = beg[:comma]
			rest = beg[comma+1:]
		}
		if isInfo {
			k.infoTags[name] = true
			ninfo++
		} else {
			k.fmtTags[name] = true
			nfmt++
		}
		beg = rest
	}
	if !keepInfo && !keepFmt && ninfo == 0 && nfmt == 0 {
		keepInfo, keepFmt = true, true
	}
	if !keepFmt && nfmt == 0 {
		keepFmt = true
	}
	k.keepAllInfo = keepInfo
	k.keepAllFmt = keepFmt
	return k
}

// uniqueNames reproduces create_unique_file_name: it sanitises a template
// (replacing [ \t:/\] with "_") and appends "-N" on clashes.
type uniqueNames struct {
	seen map[string]bool
}

// make returns a unique sanitised filename for template.
func (u *uniqueNames) make(template string) string {
	if u.seen == nil {
		u.seen = map[string]bool{}
	}
	var b strings.Builder
	for i := 0; i < len(template); i++ {
		c := template[i]
		if c == ':' || c == '\\' || c == '/' || c == ' ' || c == '\t' {
			b.WriteByte('_')
		} else {
			b.WriteByte(c)
		}
	}
	base := b.String()
	name := base
	id := 0
	for u.seen[name] {
		id++
		name = base + "-" + strconv.Itoa(id)
	}
	u.seen[name] = true
	return name
}

// readNonEmptyLines reads a plain-text file into trimmed non-empty lines.
func readNonEmptyLines(path string) ([]string, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var out []string
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimRight(line, "\r")
		if strings.TrimSpace(line) == "" {
			continue
		}
		out = append(out, line)
	}
	return out, nil
}
