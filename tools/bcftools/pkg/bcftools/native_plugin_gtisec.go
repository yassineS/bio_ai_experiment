// Native port of the upstream `GTisec` plugin (plugins/GTisec.c). It counts, for
// every possible non-empty subset of the input samples, how many genotypes are
// shared by exactly that subset of samples, and prints the counts in banker's
// sequence order. The VCF/BCF output is suppressed; only the textual report is
// written. The plugin is stateful (it accumulates subset counts across every
// record) and therefore runs serially.
package bcftools

import (
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("GTisec", func() NativePlugin { return &gtisecPlugin{} }) }

// GTisec output flags, matching the MISSING/VERBOSE/SMPORDER bits in GTisec.c.
const (
	gtisecMissing  = 1 << 0
	gtisecVerbose  = 1 << 1
	gtisecSmpOrder = 1 << 2
)

// gtisecPlugin implements the `GTisec` plugin.
type gtisecPlugin struct {
	hdr        *vcf.Header
	nsmp       int
	nsmpp2     uint32 // 2^nsmp
	flag       uint8
	bankers    []uint32 // banker's sequence position -> subset bitmask
	quick      []uint64 // n-choose-k memo (the upstream "quick" table)
	missingGTs []uint64 // per-sample missing-genotype counts (when -m)
	smpIs      []uint64 // subset bitmask -> shared-genotype count
	out        io.Writer
	argv       []string
}

// SuppressVCF reports true: GTisec emits only its textual report.
func (p *gtisecPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *gtisecPlugin) SetStdout(w io.Writer) { p.out = w }

// SetArgv records the upstream-equivalent argv for the command-line header line.
func (p *gtisecPlugin) SetArgv(argv []string) { p.argv = argv }

// Name returns the plugin name.
func (p *gtisecPlugin) Name() string { return "GTisec" }

// About returns the one-line description, matching GTisec.c about().
func (p *gtisecPlugin) About() string {
	return "Count genotype intersections across all possible sample subsets in a vcf file.\n"
}

// Parallel reports false: subset counts are accumulated serially.
func (p *gtisecPlugin) Parallel() bool { return false }

// Init parses the -m/-v/-H flags, validates samples and GT, allocates the
// counters, computes the banker's sequence and prints the report header.
func (p *gtisecPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	for i := 0; i < len(args); i++ {
		switch args[i] {
		case "-m", "--missing":
			p.flag |= gtisecMissing
		case "-v", "--verbose":
			p.flag |= gtisecVerbose
		case "-H", "--human-readable":
			p.flag |= gtisecSmpOrder | gtisecVerbose
		default:
			return nil, fmt.Errorf("GTisec: unsupported option %q", args[i])
		}
	}
	p.hdr = hdr
	if len(hdr.Samples) == 0 {
		return nil, fmt.Errorf("GTisec: no samples in input file")
	}
	p.nsmp = len(hdr.Samples)
	if p.nsmp > 32 {
		return nil, fmt.Errorf("GTisec: too many samples; a maximum of 32 is supported")
	}
	if !hasFormatHeader(hdr.MetaInfo, "GT") {
		return nil, fmt.Errorf("GTisec: GT not present in the header")
	}
	p.nsmpp2 = 1 << uint(p.nsmp)
	p.bankers = make([]uint32, p.nsmpp2)
	p.quick = make([]uint64, (p.nsmp*(p.nsmp+1))/4+1)
	if p.flag&gtisecMissing != 0 {
		p.missingGTs = make([]uint64, p.nsmp)
	}
	p.smpIs = make([]uint64, p.nsmpp2)

	p.printHeader()

	for j := uint32(0); j < p.nsmpp2; j++ {
		p.bankers[j] = p.computeBankers(uint64(j))
	}
	return hdr, nil
}

// printHeader prints the multi-line report header, matching GTisec.c init().
func (p *gtisecPlugin) printHeader() {
	if p.out == nil {
		return
	}
	fp := p.out
	// The "produced by" and "command line was" lines are provenance and are
	// stripped from the oracle comparison; emit them so the surrounding
	// structure matches upstream after stripping.
	fmt.Fprint(fp, "# This file was produced by bcftools +GTisec (bio_ai_experiment+htslib-bio_ai_experiment)\n")
	fmt.Fprintf(fp, "# The command line was:\tbcftools +GTisec %s\n", strings.Join(p.argv, " "))
	fmt.Fprint(fp, "# This file can be used as input to the subset plotting tools at:\n"+
		"#   https://github.com/dlaehnemann/bankers2\n")
	fmt.Fprint(fp, "# Genotype intersections across samples:\n")
	fmt.Fprint(fp, "@SMPS")
	for i := p.nsmp - 1; i >= 0; i-- {
		fmt.Fprintf(fp, " %s", p.hdr.Samples[i])
	}
	fmt.Fprintln(fp)
	if p.flag&gtisecMissing != 0 {
		if p.flag&gtisecSmpOrder != 0 {
			fmt.Fprint(fp, "# The first line of each sample contains its count of missing genotypes, with a '-' appended\n"+
				"#   to the sample name.\n")
		} else {
			fmt.Fprintf(fp, "# The first %d lines contain the counts for missing values of each sample in the order provided\n"+
				"#   in the SMPS-line above. Intersection counts only start afterwards.\n", p.nsmp)
		}
	}
	if p.flag&gtisecSmpOrder != 0 {
		fmt.Fprint(fp, "# Human readable output (-H) was requested. Subset intersection counts are therefore sorted by\n"+
			"#   sample and repeated for each contained sample. For each sample, counts are in banker's \n"+
			"#   sequence order regarding all other samples.\n")
	} else {
		fmt.Fprint(fp, "# Subset intersection counts are in global banker's sequence order.\n")
		if p.nsmp > 2 {
			fmt.Fprintf(fp, "#   After exclusive sample counts in order of the SMPS-line, banker's sequence continues with:\n"+
				"#   %s,%s   %s,%s   ...\n",
				p.hdr.Samples[p.nsmp-1], p.hdr.Samples[p.nsmp-2],
				p.hdr.Samples[p.nsmp-1], p.hdr.Samples[p.nsmp-3])
		}
	}
	if p.flag&gtisecVerbose != 0 {
		fmt.Fprint(fp, "# [1] Number of shared non-ref genotypes \t[2] Samples sharing non-ref genotype (GT)\n")
	} else {
		fmt.Fprint(fp, "# [1] Number of shared non-ref genotypes\n")
	}
}

// choose computes the binomial coefficient n-choose-k with the memo table,
// mirroring GTisec.c choose().
func (p *gtisecPlugin) choose(n, k uint) uint64 {
	if n == 0 {
		return 0
	}
	if n == k || k == 0 {
		return 1
	}
	if k > n/2 {
		k = n - k
	}
	i := (n*(n-1))/4 + k - 1
	if p.quick[i] == 0 {
		p.quick[i] = p.choose(n-1, k-1) + p.choose(n-1, k)
	}
	return p.quick[i]
}

// computeBankers returns the subset bitmask at banker's sequence position a,
// mirroring GTisec.c compute_bankers() including the memoised symmetry shortcut.
func (p *gtisecPlugin) computeBankers(a uint64) uint32 {
	if a == 0 {
		return 0
	}
	if p.bankers[a] == 0 {
		if a >= uint64(p.nsmpp2)/2 {
			p.bankers[a] = p.computeBankers(uint64(p.nsmpp2)-(a+1)) ^ (p.nsmpp2 - 1)
			return p.bankers[a]
		}
		var c uint
		n := uint(p.nsmp)
		e := a
		binom := p.choose(n, c)
		for {
			e -= binom
			c++
			binom = p.choose(n, c)
			if binom > e {
				break
			}
		}
		for {
			if e == 0 {
				c--
				p.bankers[a] |= 1
			} else {
				binom = p.choose(n-1, c-1)
				if binom > e {
					c--
					p.bankers[a] |= 1
				} else {
					e -= binom
				}
			}
			n--
			if n == 0 || c == 0 {
				break
			}
			p.bankers[a] <<= 1
		}
		p.bankers[a] <<= n
	}
	return p.bankers[a]
}

// Process accumulates, for one record, the per-subset shared-genotype counts.
// For each genotype value at the site, the samples carrying it form a bitmask,
// and that subset's counter is incremented (one increment per distinct genotype
// observed). Missing genotypes are skipped (and optionally tallied with -m).
func (p *gtisecPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	// Map genotype-key -> sample bitmask. The key is a canonical encoding of the
	// sample's genotype: samples that carry the same genotype map to the same
	// key, and that key's accumulated sample bitmask is incremented once per
	// distinct genotype seen at the site (mirroring GTisec.c's gts2smps hash).
	gts := map[string]uint32{}
	for i := 0; i < p.nsmp; i++ {
		gt, ok := sampleGT(v, i)
		if !ok || len(gt.alleles) == 0 {
			if p.flag&gtisecMissing != 0 {
				p.missingGTs[i]++
			}
			continue
		}
		key, missing := gtisecGTKey(gt)
		if missing {
			// Any missing allele makes the genotype uninformative for sharing,
			// matching upstream's bcf_gt_is_missing checks (generalised here to
			// arbitrary ploidy: upstream only inspects the first two slots
			// because it rejects ploidy>2 outright — see UPSTREAM_BUGS.md).
			if p.flag&gtisecMissing != 0 {
				p.missingGTs[i]++
			}
			continue
		}
		gts[key] |= 1 << uint(i)
	}
	for _, mask := range gts {
		p.smpIs[mask]++
	}
	return nil, nil
}

// gtisecGTKey returns a canonical key identifying a sample's genotype for the
// purpose of GTisec's intersection counting, plus whether the genotype is
// missing (and so should be skipped). Two genotypes share a key iff they are the
// same UNORDERED multiset of alleles at the same ploidy — generalising upstream
// GTisec's bcf_alleles2gt(a,b) (an unordered diploid pair) to arbitrary ploidy.
// A genotype is treated as missing if any of its (non vector-end) alleles is
// missing. For ploidy 1 and 2 the resulting partition is byte-identical to
// upstream's bcf_alleles2gt keys; for higher ploidy it is the natural extension
// (sort the alleles, then join), which upstream rejects outright.
func gtisecGTKey(gt genotype) (key string, missing bool) {
	alleles := make([]int, 0, len(gt.alleles))
	for _, a := range gt.alleles {
		if int64(a) == vectorEndAllele {
			// Defensive: a decoded vector-end pad is not a real allele.
			continue
		}
		if a == missingAllele {
			return "", true
		}
		alleles = append(alleles, a)
	}
	if len(alleles) == 0 {
		return "", true
	}
	// Canonicalise as an unordered multiset: sort ascending, then join. The
	// ploidy is implied by the element count, so [1] (haploid) and [1,1]
	// (diploid homozygous) stay distinct.
	sort.Ints(alleles)
	var b strings.Builder
	for j, a := range alleles {
		if j > 0 {
			b.WriteByte('/')
		}
		b.WriteString(strconv.Itoa(a))
	}
	return b.String(), false
}

// Destroy prints the accumulated subset counts in the requested ordering,
// mirroring GTisec.c destroy().
func (p *gtisecPlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	fp := p.out
	switch {
	case p.flag&gtisecSmpOrder != 0:
		for s := p.nsmp - 1; s >= 0; s-- {
			if p.flag&gtisecMissing != 0 {
				fmt.Fprintf(fp, "%d\t%s-\n", p.missingGTs[s], p.hdr.Samples[s])
			}
			for i := uint32(1); i < p.nsmpp2; i++ {
				if (p.bankers[i]>>uint(s))&1 == 0 {
					continue
				}
				fmt.Fprintf(fp, "%d\t", p.smpIs[p.bankers[i]])
				var b strings.Builder
				b.WriteString(p.hdr.Samples[s])
				for j := p.nsmp - 1; j >= 0; j-- {
					if (p.bankers[i]^(1<<uint(s)))&(1<<uint(j)) != 0 {
						b.WriteByte(',')
						b.WriteString(p.hdr.Samples[j])
					}
				}
				b.WriteByte('\n')
				fp.Write([]byte(b.String()))
			}
		}
	case p.flag&gtisecVerbose != 0:
		if p.flag&gtisecMissing != 0 {
			for s := p.nsmp - 1; s >= 0; s-- {
				fmt.Fprintf(fp, "%d\t%s-\n", p.missingGTs[s], p.hdr.Samples[s])
			}
		}
		for i := uint32(1); i < p.nsmpp2; i++ {
			fmt.Fprintf(fp, "%d\t", p.smpIs[p.bankers[i]])
			n := 0
			for s := p.nsmp - 1; s >= 0; s-- {
				if (p.bankers[i]>>uint(s))&1 != 0 {
					if n != 0 {
						fp.Write([]byte(","))
					}
					fp.Write([]byte(p.hdr.Samples[s]))
					n = 1
				}
			}
			fp.Write([]byte("\n"))
		}
	default:
		if p.flag&gtisecMissing != 0 {
			for s := p.nsmp - 1; s >= 0; s-- {
				fmt.Fprintf(fp, "%d\n", p.missingGTs[s])
			}
		}
		for i := uint32(1); i < p.nsmpp2; i++ {
			fmt.Fprintf(fp, "%d\n", p.smpIs[p.bankers[i]])
		}
	}
	return nil
}

// vectorEndAllele is the decoded value of bcf_int32_vector_end after
// bcf_gt_allele in htslib (= (0x80000001>>1)-1). It marks an unused (padded)
// genotype slot for samples of below-maximum ploidy; gtisecGTKey drops such
// slots so a haploid sample is not conflated with a padded diploid one.
const vectorEndAllele = int64(0x40000000 - 1)
