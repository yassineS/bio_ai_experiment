// Native port of the upstream `scatter` plugin (plugins/scatter.c). It scatters
// a VCF into multiple output files, either by fixed-size chunks (-n N, files
// named "0", "1", ... by chunk index) or by genomic regions (-s/-S, files named
// by the region or its optional second-column label, with an optional -x extra
// file for records overlapping none). Output filenames are prefix + label +
// container suffix, with whitespace in the label replaced by "_", matching
// open_set in scatter.c.
//
// This is a multiOutputPlugin owning all writers. The -n chunk mode, the -s
// region list and -S region file modes (with the optional second naming column),
// the -x extra file, the -p prefix and the -O container/level options are
// supported. The -i/-e options are accepted and validated but, exactly as in
// upstream scatter.c, applied to NOTHING: scatter.c parses filter_str/filter_logic
// (and errors if both -i and -e are given) yet never calls filter_init or
// filter_test, so the expression has no effect on the scattered output. We
// reproduce that no-op faithfully. -W/--write-index indexes every scattered
// output file (a CSI by default, a TBI for -W=tbi on VCF.gz) exactly as upstream
// does; the -r/-R/-t/-T region selection is applied before the scatter.
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
	registerNativePlugin("scatter", func() NativePlugin { return &scatterPlugin{} })
}

// scatterRegion is one parsed scatter region with its 0-based inclusive span and
// the output-file label it routes to.
type scatterRegion struct {
	chrom string
	beg   int // 0-based inclusive
	end   int // 0-based inclusive
	label string
}

// scatterPlugin implements scatter.
type scatterPlugin struct {
	outputDir string
	prefix    string
	extra     string
	scatter   string // -s value
	scatterFn string // -S file
	nsites    int
	format    OutputFormat
	clevel    int

	filterStr string // -i/-e expression; accepted and validated for "only one",
	filterSet bool   // but applied to nothing — upstream scatter.c never filters.

	writeIndex writeIndexFmt // -W/--write-index[=FMT]; writeIndexOff when unset

	rt regionTargetFilter // shared -r/-R/-t/-T selection applied before scattering.
}

// SetRegionTarget records the shared -r/-R/-t/-T selection the framework parsed
// out of scatter's argv; it is applied to the input records before they are
// routed into chunk/region output files.
func (p *scatterPlugin) SetRegionTarget(f regionTargetFilter) { p.rt = f }

// Name returns the plugin name.
func (p *scatterPlugin) Name() string { return "scatter" }

// RegionTargetCaps opts scatter into the shared -r/-R/-t/-T region/target filter,
// applied (via SetRegionTarget) to the records before they are scattered.
func (p *scatterPlugin) RegionTargetCaps() regionTargetCaps { return overlapRegionTargetCaps }

// About returns the one-line description, matching scatter.c about().
func (p *scatterPlugin) About() string {
	return "Scatter VCF by chunks or regions, creating multiple VCFs.\n"
}

// RunStyle reports that scatter is a run()-style plugin: upstream's scatter.c
// exports a `run` symbol, so it owns its entire argv (including -o/-O/-n/-s)
// before the trailing input-file positional, with no `--` separator
// (e.g. `bcftools +scatter -o DIR -n 2 FILE`).
func (p *scatterPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of scatter's flags consumes the following
// CLI token as its value, so the host can separate the lone input-file
// positional from the plugin options.
func (p *scatterPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-o", "--output", "-O", "--output-type",
		"-s", "--scatter", "-S", "--scatter-file",
		"-n", "--nsites-per-chunk", "-x", "--extra",
		"-i", "--include", "-e", "--exclude",
		"-r", "--regions", "-R", "--regions-file",
		"-t", "--targets", "-T", "--targets-file",
		"--regions-overlap", "--targets-overlap",
		"-p", "--prefix", "--threads", "-v", "--verbosity", "--hts-opts":
		return true
	}
	return false
}

// Init parses and validates the plugin arguments.
func (p *scatterPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.format = OutputVCF
	p.clevel = -1
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("scatter: option %q requires an argument", a)
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
		case "-s", "--scatter":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.scatter = v
		case "-S", "--scatter-file":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.scatterFn = v
		case "-n", "--nsites-per-chunk":
			v, err := next()
			if err != nil {
				return nil, err
			}
			n, perr := strconv.Atoi(v)
			if perr != nil || n <= 0 {
				return nil, fmt.Errorf("scatter: positive integer required for --nsites-per-chunk: %q", v)
			}
			p.nsites = n
		case "-x", "--extra":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.extra = v
		case "-p", "--prefix":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.prefix = v
		case "--no-version":
			// Provenance lines are stripped during parity comparison anyway.
		case "--threads", "-v", "--verbosity", "--hts-opts":
			if _, err := next(); err != nil {
				return nil, err
			}
		case "-i", "--include", "-e", "--exclude":
			// Upstream scatter.c parses -i/-e but never applies them (no
			// filter_init/filter_test), so the only behaviour to reproduce is the
			// "only one of -i or -e" guard. The expression itself is a no-op.
			v, err := next()
			if err != nil {
				return nil, err
			}
			if p.filterSet {
				return nil, fmt.Errorf("scatter: only one -i or -e expression can be given, and they cannot be combined")
			}
			p.filterStr = v
			p.filterSet = true
		default:
			if sel, handled, werr := parseWriteIndexArg(a); handled {
				if werr != nil {
					return nil, fmt.Errorf("scatter: %w", werr)
				}
				p.writeIndex = sel
				continue
			}
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
			if strings.HasPrefix(a, "-p") && len(a) > 2 {
				p.prefix = a[2:]
				continue
			}
			if strings.HasPrefix(a, "-n") && len(a) > 2 {
				n, perr := strconv.Atoi(a[2:])
				if perr != nil || n <= 0 {
					return nil, fmt.Errorf("scatter: positive integer required for --nsites-per-chunk: %q", a[2:])
				}
				p.nsites = n
				continue
			}
			return nil, fmt.Errorf("scatter: unsupported option %q", a)
		}
	}
	if p.outputDir == "" {
		return nil, fmt.Errorf("scatter: missing the -o option")
	}
	if p.nsites == 0 && p.scatter == "" && p.scatterFn == "" {
		return nil, fmt.Errorf("scatter: missing either the -n or one of the -s or -S options")
	}
	if p.nsites != 0 && (p.scatter != "" || p.scatterFn != "") {
		return nil, fmt.Errorf("scatter: only one of -n or either -s or -S can be given")
	}
	if p.nsites != 0 && p.extra != "" {
		return nil, fmt.Errorf("scatter: cannot use -x together with -n")
	}
	return hdr, nil
}

// parseOutputType parses upstream's -O u|b|v|z[0-9] spelling.
func (p *scatterPlugin) parseOutputType(s string) error {
	if s == "" {
		return fmt.Errorf("scatter: empty -O argument")
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
			return fmt.Errorf("scatter: the output type %q not recognised", s)
		}
		p.clevel = n
		return nil
	}
	if len(s) > 1 {
		n, err := strconv.Atoi(s[1:])
		if err != nil || n < 0 || n > 9 {
			return fmt.Errorf("scatter: could not parse compression level %q", s[1:])
		}
		p.clevel = n
	}
	return nil
}

// Process is unused: scatter is a multiOutputPlugin.
func (p *scatterPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *scatterPlugin) Destroy() error { return nil }

// RunMulti reads the input and dispatches to the chunk or region scatter path.
func (p *scatterPlugin) RunMulti(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	if _, err := p.Init(opts.Args, nil); err != nil {
		return err
	}
	hdr, variants, err := readPluginInput(opts, stderr)
	if err != nil {
		return err
	}
	variants = p.rt.apply(variants)
	if err := os.MkdirAll(p.outputDir, 0o777); err != nil {
		return err
	}
	if p.nsites != 0 {
		return p.scatterByChunks(hdr, variants)
	}
	return p.scatterByRegions(hdr, variants)
}

// scatterByChunks writes consecutive runs of nsites records into files named by
// the 0-based chunk index, matching the -n path of process().
func (p *scatterPlugin) scatterByChunks(hdr *vcf.Header, variants []*vcf.Variant) error {
	chunk := 0
	for i := 0; i < len(variants); i += p.nsites {
		end := i + p.nsites
		if end > len(variants) {
			end = len(variants)
		}
		if err := p.writeFile(strconv.Itoa(chunk), hdr, variants[i:end]); err != nil {
			return err
		}
		chunk++
	}
	return nil
}

// scatterByRegions routes each record to every region file it overlaps; records
// overlapping no region go to the -x extra file when given. Files are created in
// first-seen order (the order the labels appear in -s/-S), with -x appended
// last, matching init_data + process().
func (p *scatterPlugin) scatterByRegions(hdr *vcf.Header, variants []*vcf.Variant) error {
	regions, labels, err := p.parseRegions()
	if err != nil {
		return err
	}
	// Group records per label in input order.
	recs := map[string][]*vcf.Variant{}
	for _, v := range variants {
		matched := false
		for _, r := range regions {
			if r.chrom == v.Chrom && v.Pos-1 >= r.beg && v.Pos-1 <= r.end {
				recs[r.label] = append(recs[r.label], v)
				matched = true
			}
		}
		if !matched && p.extra != "" {
			recs[p.extra] = append(recs[p.extra], v)
		}
	}
	for _, label := range labels {
		if err := p.writeFile(label, hdr, recs[label]); err != nil {
			return err
		}
	}
	if p.extra != "" {
		if err := p.writeFile(p.extra, hdr, recs[p.extra]); err != nil {
			return err
		}
	}
	return nil
}

// parseRegions parses the -s comma list or -S file into regions and the ordered
// unique list of output labels (first occurrence order).
func (p *scatterPlugin) parseRegions() ([]scatterRegion, []string, error) {
	var lines []string
	if p.scatterFn != "" {
		fl, err := readNonEmptyLines(p.scatterFn)
		if err != nil {
			return nil, nil, err
		}
		lines = fl
	} else {
		lines = strings.Split(p.scatter, ",")
	}
	var regions []scatterRegion
	var labels []string
	seen := map[string]bool{}
	for _, line := range lines {
		r, ok, err := parseScatterRegion(line)
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			continue
		}
		regions = append(regions, r)
		if !seen[r.label] {
			seen[r.label] = true
			labels = append(labels, r.label)
		}
	}
	return regions, labels, nil
}

// parseScatterRegion parses one region spec "chr[:beg[-end]][\tlabel]", porting
// regidx_parse_reg_name. The coordinates are 1-based on input and stored 0-based;
// the label defaults to the whole region token when no second column is present.
func parseScatterRegion(line string) (scatterRegion, bool, error) {
	s := strings.TrimLeft(line, " \t")
	if s == "" || s[0] == '#' {
		return scatterRegion{}, false, nil
	}
	// Split the leading region token from an optional trailing label.
	end := 0
	for end < len(s) && s[end] != ' ' && s[end] != '\t' {
		end++
	}
	regTok := s[:end]
	rest := strings.TrimLeft(s[end:], " \t")

	r := scatterRegion{}
	colon := strings.IndexByte(regTok, ':')
	if colon < 0 {
		r.chrom = regTok
		r.beg = 0
		r.end = scatterMaxCoord
	} else {
		r.chrom = regTok[:colon]
		coords := regTok[colon+1:]
		if dash := strings.IndexByte(coords, '-'); dash >= 0 {
			beg, err1 := strconv.Atoi(coords[:dash])
			if err1 != nil || beg == 0 {
				return scatterRegion{}, false, fmt.Errorf("scatter: could not parse region: %s", line)
			}
			r.beg = beg - 1
			if coords[dash+1:] == "" {
				r.end = scatterMaxCoord
			} else {
				e, err2 := strconv.Atoi(coords[dash+1:])
				if err2 != nil {
					return scatterRegion{}, false, fmt.Errorf("scatter: could not parse region: %s", line)
				}
				r.end = e - 1
			}
		} else {
			beg, err := strconv.Atoi(coords)
			if err != nil || beg == 0 {
				return scatterRegion{}, false, fmt.Errorf("scatter: could not parse region: %s", line)
			}
			r.beg = beg - 1
			r.end = beg - 1
		}
	}
	if rest != "" {
		r.label = rest
	} else {
		r.label = regTok
	}
	return r, true, nil
}

// scatterMaxCoord is the open-ended region upper bound (REGIDX_MAX), large
// enough to include any record position.
const scatterMaxCoord = int(^uint32(0) >> 1)

// writeFile writes a single scatter output file: prefix + sanitised(label) +
// suffix, with whitespace in the label replaced by "_".
func (p *scatterPlugin) writeFile(label string, hdr *vcf.Header, variants []*vcf.Variant) error {
	name := sanitizeScatterLabel(p.prefix + label)
	path := filepath.Join(p.outputDir, name+p.suffix())
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("scatter: cannot write to %q: %w", path, err)
	}
	defer f.Close()
	outHdr := &vcf.Header{Samples: hdr.Samples}
	outHdr.MetaInfo = append(outHdr.MetaInfo, hdr.MetaInfo...)
	w, cleanup, err := openOutput(f, ViewOptions{OutputFormat: p.format, CompressLevel: p.clevel}, outHdr)
	if err != nil {
		return err
	}
	if err := w.WriteHeader(); err != nil {
		cleanup()
		return err
	}
	for _, v := range variants {
		if err := w.Write(v); err != nil {
			cleanup()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		cleanup()
		return err
	}
	cleanup()
	if err := f.Close(); err != nil {
		return err
	}
	// -W/--write-index: index each scattered output file as upstream does.
	if err := writeIndexFor(path, p.format, p.writeIndex); err != nil {
		return fmt.Errorf("scatter: %w", err)
	}
	return nil
}

// suffix returns the container suffix from -O.
func (p *scatterPlugin) suffix() string {
	switch p.format {
	case OutputBCF, OutputBCFUncompressed:
		return ".bcf"
	case OutputVCFGz:
		return ".vcf.gz"
	default:
		return ".vcf"
	}
}

// sanitizeScatterLabel replaces any whitespace in the prefixed label with "_",
// matching open_set's isspace replacement (only the prefix+label portion).
func sanitizeScatterLabel(s string) string {
	var b strings.Builder
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f' {
			b.WriteByte('_')
		} else {
			b.WriteByte(c)
		}
	}
	return b.String()
}
