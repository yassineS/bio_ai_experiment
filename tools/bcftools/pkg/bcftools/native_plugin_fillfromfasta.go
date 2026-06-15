// Native port of the upstream `fill-from-fasta` plugin
// (plugins/fill-from-fasta.c). It fills the REF column or an INFO tag from the
// bases found at the record position in a FASTA reference. The FASTA is loaded
// in Init via pkg/htsgo/fasta (a sibling .fai is used when present, otherwise an
// index is built on the fly), and each record fetches the bases spanning its
// REF allele. Bases are uppercased to match the C implementation, and with
// -N/--replace-non-ACGTN any base outside {A,C,G,T,N} becomes N.
package bcftools

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("fill-from-fasta", func() NativePlugin { return &fillFromFastaPlugin{} })
}

// fill-from-fasta annotation targets, mirroring the ANNO_* constants in
// fill-from-fasta.c.
const (
	fffAnnoRef    = iota + 1 // fill the REF column
	fffAnnoString            // fill an INFO tag of Type=String
	fffAnnoInt               // fill an INFO tag of Type=Integer (single base only)
)

// fillFromFastaPlugin implements the `fill-from-fasta` plugin. Each Process
// only reads from the loaded FASTA, so it is per-record and parallel.
type fillFromFastaPlugin struct {
	column         string // the INFO tag name (REF stripped to anno==fffAnnoRef)
	anno           int    // one of fffAnnoRef / fffAnnoString / fffAnnoInt
	replaceNonACGT bool   // -N/--replace-non-ACGTN
	fai            *fasta.RandomAccess
	filter         *Filter
	filterExclude  bool // true for -e/--exclude, false for -i/--include
}

// Name returns the plugin name.
func (p *fillFromFastaPlugin) Name() string { return "fill-from-fasta" }

// About returns the one-line description, matching fill-from-fasta.c about().
func (p *fillFromFastaPlugin) About() string {
	return "Fill INFO or REF field based on values in a fasta file"
}

// Parallel reports true: each record reads independently from the FASTA.
func (p *fillFromFastaPlugin) Parallel() bool { return true }

// Init parses the plugin options, loads the FASTA, appends any extra header
// lines (-h), and resolves the target column's annotation type from the header.
func (p *fillFromFastaPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	var refFname, headerFname, filterStr string
	var haveInclude, haveExclude bool
	for i := 0; i < len(args); i++ {
		a := args[i]
		needVal := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("fill-from-fasta: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-c", "--column":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			p.column = v
		case "-f", "--fasta":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			refFname = v
		case "-h", "--header-lines":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			headerFname = v
		case "-i", "--include":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			filterStr = v
			haveInclude = true
		case "-e", "--exclude":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			filterStr = v
			haveExclude = true
		case "-N", "--replace-non-ACGTN":
			p.replaceNonACGT = true
		default:
			return nil, fmt.Errorf("fill-from-fasta: unsupported option %q", a)
		}
	}

	if haveInclude && haveExclude {
		return nil, fmt.Errorf("fill-from-fasta: only one of -i or -e can be given")
	}
	if p.column == "" {
		return nil, fmt.Errorf("fill-from-fasta: --column option is required")
	}
	if refFname == "" {
		return nil, fmt.Errorf("fill-from-fasta: no fasta given")
	}

	// Build the output header, appending the optional header-lines file first
	// so that a freshly declared INFO tag can be resolved below.
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	if headerFname != "" {
		lines, err := readHeaderLinesFile(headerFname)
		if err != nil {
			return nil, fmt.Errorf("fill-from-fasta: %w", err)
		}
		for _, line := range lines {
			out.MetaInfo = appendInfoHeader(out.MetaInfo, line)
		}
	}

	// Resolve the annotation target: REF column, or an INFO tag whose declared
	// Type decides between string and integer fills.
	if strings.EqualFold(p.column, "REF") {
		p.anno = fffAnnoRef
	} else {
		col := p.column
		if strings.HasPrefix(strings.ToUpper(col), "INFO/") {
			col = col[len("INFO/"):]
			p.column = col
		}
		typ, ok := infoTagType(out.MetaInfo, col)
		if !ok {
			return nil, fmt.Errorf("fill-from-fasta: no header ID found for %s. Header lines can be added with the --header-lines option", col)
		}
		switch strings.ToLower(typ) {
		case "integer":
			p.anno = fffAnnoInt
		case "string", "character":
			p.anno = fffAnnoString
		default:
			return nil, fmt.Errorf("fill-from-fasta: the type of %s not recognised (%s)", col, typ)
		}
	}

	fai, err := fasta.OpenRandomAccess(refFname)
	if err != nil {
		return nil, fmt.Errorf("fill-from-fasta: %w", err)
	}
	p.fai = fai

	if filterStr != "" {
		f, err := CompileFilter(filterStr)
		if err != nil {
			return nil, fmt.Errorf("fill-from-fasta: -i/-e expression: %w", err)
		}
		p.filter = f
		p.filterExclude = haveExclude
	}

	return out, nil
}

// Process fills the configured target for one record. When a filter is active
// and the record does not pass it, the record is returned unchanged.
func (p *fillFromFastaPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	if p.filter != nil {
		ret := p.filter.Eval(v)
		if p.filterExclude {
			if ret {
				return []*vcf.Variant{v}, nil
			}
		} else if !ret {
			return []*vcf.Variant{v}, nil
		}
	}

	refLen := len(v.Ref)
	if refLen == 0 {
		return []*vcf.Variant{v}, nil
	}
	// faidx_fetch_seq(fai, chr, pos, pos+ref_len-1): 0-based inclusive in C,
	// which is the half-open [pos, pos+ref_len) range here. Pos is 1-based.
	start := int64(v.Pos - 1)
	end := start + int64(refLen)
	fa, err := p.fai.Fetch(v.Chrom, start, end)
	if err != nil {
		return nil, fmt.Errorf("fill-from-fasta: faidx fetch failed at %s:%d: %w", v.Chrom, v.Pos, err)
	}
	// Fetch already uppercases; apply the -N replacement to non-ACGTN bases.
	if p.replaceNonACGT {
		for i, b := range fa {
			if b != 'A' && b != 'C' && b != 'G' && b != 'T' && b != 'N' {
				fa[i] = 'N'
			}
		}
	}

	switch p.anno {
	case fffAnnoRef:
		v.Ref = string(fa)
	case fffAnnoString:
		setInfo(v, p.column, string(fa))
	case fffAnnoInt:
		if refLen == 1 {
			// atoi on the single base, exactly as the C code: a non-digit base
			// parses to 0.
			val := 0
			if fa[0] >= '0' && fa[0] <= '9' {
				val = int(fa[0] - '0')
			}
			setInfo(v, p.column, strconv.Itoa(val))
		}
	}
	return []*vcf.Variant{v}, nil
}

// Destroy releases the FASTA handle.
func (p *fillFromFastaPlugin) Destroy() error {
	if p.fai != nil {
		return p.fai.Close()
	}
	return nil
}
