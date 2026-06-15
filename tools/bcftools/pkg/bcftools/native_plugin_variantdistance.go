// Native port of the upstream `variant-distance` plugin
// (plugins/variant-distance.c). It annotates every site with the distance to
// the nearest variant via a custom INFO tag (default DIST), supporting four
// directionalities: nearest, fwd (next), rev (previous), and both (a Number=2
// tag carrying previous,next). Records sharing a chromosome+position are a
// duplicate block and all receive the same distance. The output VCF/BCF is
// emitted with the new INFO tag; the plugin needs cross-record look-ahead and
// look-back, so it is a serial bufferedPlugin.
package bcftools

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("variant-distance", func() NativePlugin { return &variantDistancePlugin{} })
}

// variant-distance directionalities, matching DIR_* in variant-distance.c.
const (
	vdNearest = iota
	vdFwd
	vdRev
	vdBoth
)

// variantDistancePlugin implements the `variant-distance` plugin.
type variantDistancePlugin struct {
	tag       string
	direction int
}

// Name returns the plugin name.
func (p *variantDistancePlugin) Name() string { return "variant-distance" }

// About returns the one-line description, matching variant-distance.c about().
func (p *variantDistancePlugin) About() string {
	return "Annotate sites with distance to the nearest variant"
}

// RunStyle reports that variant-distance is a run()-style plugin: upstream
// exports a `run` symbol for it, so its options precede the input file and
// there is no `--` separator (e.g. `+variant-distance -d nearest FILE`).
func (p *variantDistancePlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of variant-distance's own flags consumes
// the following CLI token as its value, so the host can separate the lone
// input-file positional from the plugin options. Both -d/--direction and
// -n/--tag-name take a value; the boolean help flags do not.
func (p *variantDistancePlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-d", "--direction", "-n", "--tag-name":
		return true
	}
	return false
}

// Init parses -d/--direction and -n/--tag-name and appends the INFO header
// line describing the chosen direction.
func (p *variantDistancePlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.tag = "DIST"
	p.direction = vdNearest
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-d", "--direction":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("variant-distance: -d requires an argument")
			}
			i++
			switch strings.ToLower(args[i]) {
			case "nearest":
				p.direction = vdNearest
			case "fwd":
				p.direction = vdFwd
			case "rev":
				p.direction = vdRev
			case "both":
				p.direction = vdBoth
			default:
				return nil, fmt.Errorf("variant-distance: unknown argument to --direction: %s", args[i])
			}
		case "-n", "--tag-name":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("variant-distance: -n requires an argument")
			}
			i++
			p.tag = args[i]
		default:
			return nil, fmt.Errorf("variant-distance: unsupported option %q", a)
		}
	}

	nval := 1
	var desc string
	switch p.direction {
	case vdFwd:
		desc = "next"
	case vdRev:
		desc = "previous"
	case vdNearest:
		desc = "nearest"
	default:
		desc = "previous and next"
		nval = 2
	}
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	line := fmt.Sprintf(`##INFO=<ID=%s,Number=%d,Type=Integer,Description="Distance to the %s variant">`, p.tag, nval, desc)
	out.MetaInfo = appendInfoHeader(out.MetaInfo, line)
	return out, nil
}

// ProcessAll annotates the whole ordered stream. Records are grouped into
// duplicate blocks of identical chromosome+position; each block's rev_dist is
// the gap to the preceding block on the same chromosome and its fwd_dist the
// gap to the following block on the same chromosome (0 when there is none, e.g.
// at a chromosome boundary), mirroring the C buffering and flush logic.
func (p *variantDistancePlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	// Identify duplicate-position block boundaries.
	type block struct{ start, end int } // [start,end) over variants
	var blocks []block
	for i := 0; i < len(variants); {
		j := i + 1
		for j < len(variants) && variants[j].Chrom == variants[i].Chrom && variants[j].Pos == variants[i].Pos {
			j++
		}
		blocks = append(blocks, block{i, j})
		i = j
	}

	for bi, b := range blocks {
		var revDist, fwdDist int
		if bi > 0 {
			prev := variants[blocks[bi-1].start]
			cur := variants[b.start]
			if prev.Chrom == cur.Chrom {
				revDist = cur.Pos - prev.Pos
			}
		}
		if bi+1 < len(blocks) {
			next := variants[blocks[bi+1].start]
			cur := variants[b.start]
			if next.Chrom == cur.Chrom {
				fwdDist = next.Pos - cur.Pos
			}
		}

		val := p.distanceValue(revDist, fwdDist)
		if val == "" {
			continue
		}
		for k := b.start; k < b.end; k++ {
			setInfo(variants[k], p.tag, val)
		}
	}
	return variants, nil
}

// distanceValue renders the INFO value for one block from its rev/fwd
// distances, or "" when no tag should be added (matching the nval==0 case in
// the C flush()). For DIR_NEAREST the smaller non-zero distance is chosen.
func (p *variantDistancePlugin) distanceValue(revDist, fwdDist int) string {
	switch p.direction {
	case vdFwd:
		if fwdDist != 0 {
			return strconv.Itoa(fwdDist)
		}
	case vdRev:
		if revDist != 0 {
			return strconv.Itoa(revDist)
		}
	case vdBoth:
		if revDist != 0 || fwdDist != 0 {
			return strconv.Itoa(revDist) + "," + strconv.Itoa(fwdDist)
		}
	case vdNearest:
		if revDist != 0 || fwdDist != 0 {
			var v int
			switch {
			case fwdDist == 0:
				v = revDist
			case revDist == 0:
				v = fwdDist
			case revDist < fwdDist:
				v = revDist
			default:
				v = fwdDist
			}
			return strconv.Itoa(v)
		}
	}
	return ""
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *variantDistancePlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *variantDistancePlugin) Destroy() error { return nil }
