// Native port of the upstream `fixploidy` plugin (plugins/fixploidy.c). It
// adjusts the ploidy of FORMAT/GT genotypes to a per-region, per-sex ploidy
// (from -p/--ploidy and -s/--sex), to a forced uniform ploidy (-f), or to the
// built-in human X/Y/MT default table.
//
// Only -t GT is supported (as upstream). The region/sex ploidy machinery is
// ported in-tree (no htslib ploidy.h dependency).
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("fixploidy", func() NativePlugin { return &fixploidyPlugin{} }) }

// ploidyRegion is one CHROM,FROM,TO,SEX,PLOIDY row of a ploidy definition.
type ploidyRegion struct {
	chrom    string
	from, to int
	sex      string
	ploidy   int
}

// fixploidyPlugin implements fixploidy. It is per-record and parallel: the
// region/sex tables are read-only after Init.
type fixploidyPlugin struct {
	forcePloidy   int // -1 when not forcing
	defaultPloidy int
	regions       []ploidyRegion
	sexes         []string // ordered distinct sex names
	sample2sex    []int    // index into sexes per sample
	sampleNames   []string
	stderr        io.Writer
}

// Name returns the plugin name.
func (p *fixploidyPlugin) Name() string { return "fixploidy" }

// About returns the one-line description, matching fixploidy.c about().
func (p *fixploidyPlugin) About() string { return "Fix ploidy." }

// Parallel reports true.
func (p *fixploidyPlugin) Parallel() bool { return true }

// SetStderr wires the host stderr.
func (p *fixploidyPlugin) SetStderr(w io.Writer) { p.stderr = w }

// defaultPloidyTable is the built-in human X/Y/MT ploidy table used when
// neither -p nor -f is given (matching the string in fixploidy.c init()).
var defaultPloidyTable = []ploidyRegion{
	{"X", 1, 60000, "M", 1},
	{"X", 2699521, 154931043, "M", 1},
	{"Y", 1, 59373566, "M", 1},
	{"Y", 1, 59373566, "F", 0},
	{"MT", 1, 16569, "M", 1},
	{"MT", 1, 16569, "F", 1},
}

// Init parses the plugin options and builds the ploidy/sex tables.
func (p *fixploidyPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.forcePloidy = -1
	p.defaultPloidy = 2
	tags := "GT"
	var ploidyFile, sexFile string
	for i := 0; i < len(args); i++ {
		a := args[i]
		val := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("fixploidy: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-d", "--default-ploidy":
			v, err := val()
			if err != nil {
				return nil, err
			}
			n, err := strconv.Atoi(v)
			if err != nil {
				return nil, fmt.Errorf("fixploidy: could not parse -d %s", v)
			}
			p.defaultPloidy = n
		case "-f", "--force-ploidy":
			v, err := val()
			if err != nil {
				return nil, err
			}
			n, err := strconv.Atoi(v)
			if err != nil {
				return nil, fmt.Errorf("fixploidy: could not parse -f %s", v)
			}
			p.forcePloidy = n
		case "-p", "--ploidy":
			v, err := val()
			if err != nil {
				return nil, err
			}
			ploidyFile = v
		case "-s", "--sex":
			v, err := val()
			if err != nil {
				return nil, err
			}
			sexFile = v
		case "-t", "--tags":
			v, err := val()
			if err != nil {
				return nil, err
			}
			tags = v
		default:
			return nil, fmt.Errorf("fixploidy: unsupported option %q", a)
		}
	}
	if !strings.EqualFold(tags, "GT") {
		return nil, fmt.Errorf("fixploidy: only -t GT is currently supported, sorry")
	}

	p.sampleNames = hdr.Samples
	p.sample2sex = make([]int, len(hdr.Samples))

	if p.forcePloidy == -1 {
		if ploidyFile != "" {
			regs, err := loadPloidyFile(ploidyFile)
			if err != nil {
				return nil, err
			}
			p.regions = regs
		} else {
			// The built-in X/Y/MT table is initialised with a hardcoded default
			// ploidy of 2 (ploidy_init_string(...,2)), so -d is ignored in this
			// mode — only -p threads -d through to unlisted regions.
			p.regions = defaultPloidyTable
			p.defaultPloidy = 2
		}
		// Collect the distinct sexes from the region table plus the default "F".
		p.addSex("F")
		for _, r := range p.regions {
			p.addSex(r.sex)
		}
		// All samples default to "F"; -s overrides.
		dfltSex := p.sexIndex("F")
		for i := range p.sample2sex {
			p.sample2sex[i] = dfltSex
		}
		if sexFile != "" {
			if err := p.loadSexFile(sexFile); err != nil {
				return nil, err
			}
		}
	}
	return hdr, nil
}

// addSex registers a sex name if not already present.
func (p *fixploidyPlugin) addSex(sex string) {
	for _, s := range p.sexes {
		if s == sex {
			return
		}
	}
	p.sexes = append(p.sexes, sex)
}

// sexIndex returns the index of a sex name (-1 if absent).
func (p *fixploidyPlugin) sexIndex(sex string) int {
	for i, s := range p.sexes {
		if s == sex {
			return i
		}
	}
	return -1
}

// loadSexFile reads "SAMPLE SEX" lines and assigns each named sample its sex.
func (p *fixploidyPlugin) loadSexFile(fname string) error {
	f, err := os.Open(fname)
	if err != nil {
		return fmt.Errorf("fixploidy: could not read %s: %w", fname, err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return fmt.Errorf("fixploidy: could not parse: %s", line)
		}
		smpl, sex := fields[0], fields[1]
		si := -1
		for i, n := range p.sampleNames {
			if n == smpl {
				si = i
				break
			}
		}
		if si < 0 {
			fmt.Fprintf(p.stderrOrDiscard(), "Warning: No such sample in the VCF: %s\n", smpl)
			continue
		}
		p.addSex(sex)
		p.sample2sex[si] = p.sexIndex(sex)
	}
	return sc.Err()
}

// stderrOrDiscard returns the host stderr or io.Discard.
func (p *fixploidyPlugin) stderrOrDiscard() io.Writer {
	if p.stderr != nil {
		return p.stderr
	}
	return io.Discard
}

// Process adjusts every sample genotype to the region/sex/forced ploidy.
func (p *fixploidyPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	// Per-sex ploidy for this record's position, and the maximum ploidy.
	sexPloidy, maxPloidy := p.queryPloidy(v.Chrom, v.Pos)

	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		var target int
		if p.forcePloidy != -1 {
			target = p.forcePloidy
		} else {
			target = sexPloidy[p.sample2sex[i]]
		}
		newGT := fixPloidyGT(gt, target, maxPloidy)
		v.Samples[i].Data["GT"] = newGT.String()
	}
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *fixploidyPlugin) Destroy() error { return nil }

// queryPloidy returns the per-sex ploidy slice and the maximum ploidy at the
// given position. For force mode it returns (nil, forcePloidy).
func (p *fixploidyPlugin) queryPloidy(chrom string, pos int) ([]int, int) {
	if p.forcePloidy != -1 {
		return nil, p.forcePloidy
	}
	sexPloidy := make([]int, len(p.sexes))
	for i := range sexPloidy {
		sexPloidy[i] = p.defaultPloidy
	}
	for _, r := range p.regions {
		if r.chrom == chrom && pos+1 >= r.from && pos+1 <= r.to {
			if si := p.sexIndex(r.sex); si >= 0 {
				sexPloidy[si] = r.ploidy
			}
		}
	}
	maxPloidy := 0
	for _, pl := range sexPloidy {
		if pl > maxPloidy {
			maxPloidy = pl
		}
	}
	return sexPloidy, maxPloidy
}

// fixPloidyGT reproduces fixploidy.c's per-sample ploidy adjustment. The
// genotype is reshaped to `ploidy` alleles (the sample's own ploidy), padded by
// repeating the last allele, then to `maxPloidy` overall (the extra positions
// beyond `ploidy` becoming vector-ends, which render as fewer values in text).
func fixPloidyGT(gt genotype, ploidy, maxPloidy int) genotype {
	src := gt.alleles
	srcPhase := gt.phased
	out := genotype{}
	if ploidy == 0 {
		// Haploid-zero: a single missing allele (upstream sets bcf_gt_missing).
		return genotype{alleles: []int{missingAllele}, phased: []bool{false}}
	}
	j := 0
	for j < len(src) && j < ploidy {
		out.alleles = append(out.alleles, src[j])
		out.phased = append(out.phased, srcPhase[j])
		j++
	}
	// Expand "." to "./." and "0" to "0/0" by repeating the last allele.
	for j < ploidy {
		out.alleles = append(out.alleles, out.alleles[j-1])
		// Repeated alleles keep the phase of the position they extend; upstream
		// copies the encoded value, which carries the previous phase bit.
		out.phased = append(out.phased, out.phased[j-1])
		j++
	}
	// Positions beyond `ploidy` up to `maxPloidy` are vector-ends and are not
	// emitted in the textual representation, so we simply stop here.
	return out
}

// loadPloidyFile parses a space/tab-delimited CHROM,FROM,TO,SEX,PLOIDY file.
func loadPloidyFile(fname string) ([]ploidyRegion, error) {
	f, err := os.Open(fname)
	if err != nil {
		return nil, fmt.Errorf("fixploidy: could not read %s: %w", fname, err)
	}
	defer f.Close()
	var regs []ploidyRegion
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) == 1 && fields[0] == "*" {
			continue // catch-all marker; default handles it
		}
		if len(fields) < 5 {
			return nil, fmt.Errorf("fixploidy: could not parse ploidy line: %s", line)
		}
		from, err1 := strconv.Atoi(fields[1])
		to, err2 := strconv.Atoi(fields[2])
		ploidy, err3 := strconv.Atoi(fields[4])
		if err1 != nil || err2 != nil || err3 != nil {
			return nil, fmt.Errorf("fixploidy: could not parse ploidy line: %s", line)
		}
		regs = append(regs, ploidyRegion{fields[0], from, to, fields[3], ploidy})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	// Stable ordering for deterministic max-ploidy queries.
	sort.SliceStable(regs, func(i, j int) bool {
		if regs[i].chrom != regs[j].chrom {
			return regs[i].chrom < regs[j].chrom
		}
		return regs[i].from < regs[j].from
	})
	return regs, nil
}
