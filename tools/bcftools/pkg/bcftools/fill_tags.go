// bcftools +fill-tags — native (built-in) port of the upstream dlopen
// plugin reference_code/bcftools/plugins/fill-tags.c.
//
// The default tag set ("all") recomputes, per record, the population
// allele statistics from FORMAT/GT and writes them as INFO fields:
//
//	F_MISSING  fraction of samples with a fully-missing genotype
//	NS         number of samples with data
//	AN         total called alleles
//	AC         per-ALT allele count
//	AF         per-ALT allele frequency
//	MAF        frequency of the second most common allele
//	AC_Hom     per-ALT homozygous allele count
//	AC_Het     per-ALT heterozygous allele count
//	AC_Hemi    per-ALT hemizygous allele count
//	HWE        per-ALT Hardy-Weinberg exact-test p-value (1=good)
//	ExcHet     per-ALT excess-heterozygosity test (1=good)
//
// VAF / VAF1 header lines are declared (they are part of "all") but the
// per-record values are skipped when FORMAT/AD is absent, mirroring
// upstream. END / TYPE are excluded from "all".
//
// The counting model, the per-ALT HWE/ExcHet derivation, and the exact
// HWE test (calc_hwe) follow fill-tags.c line-for-line so the output
// byte-matches upstream (modulo provenance).
package bcftools

import (
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// fillTagsBit is a bitmask over the tags fill-tags can add.
type fillTagsBit uint

const (
	tagAN fillTagsBit = 1 << iota
	tagAC
	tagACHom
	tagACHet
	tagACHemi
	tagAF
	tagNS
	tagMAF
	tagHWE
	tagExcHet
	tagFMissing
	tagVAF
	tagVAF1
)

// tagAll is the default ("all") set: everything except END / TYPE.
const tagAll = tagAN | tagAC | tagACHom | tagACHet | tagACHemi | tagAF |
	tagNS | tagMAF | tagHWE | tagExcHet | tagFMissing | tagVAF | tagVAF1

// fillTagsHeaderLine maps a header-appendable tag to its ##INFO/##FORMAT
// declaration. The slice order is the order upstream appends them.
var fillTagsHeaderLines = []struct {
	bit  fillTagsBit
	id   string
	line string
}{
	{tagFMissing, "F_MISSING", `##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Added by +fill-tags expression F_MISSING:1=F_MISSING">`},
	{tagAN, "AN", `##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">`},
	{tagAC, "AC", `##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">`},
	{tagNS, "NS", `##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">`},
	{tagACHom, "AC_Hom", `##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">`},
	{tagACHet, "AC_Het", `##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">`},
	{tagACHemi, "AC_Hemi", `##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="Allele counts in hemizygous genotypes">`},
	{tagAF, "AF", `##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">`},
	{tagMAF, "MAF", `##INFO=<ID=MAF,Number=1,Type=Float,Description="Frequency of the second most common allele">`},
	{tagHWE, "HWE", `##INFO=<ID=HWE,Number=A,Type=Float,Description="HWE test (PMID:15789306); 1=good, 0=bad">`},
	{tagExcHet, "ExcHet", `##INFO=<ID=ExcHet,Number=A,Type=Float,Description="Test excess heterozygosity; 1=good, 0=bad">`},
	{tagVAF, "VAF", `##FORMAT=<ID=VAF,Number=A,Type=Float,Description="The fraction of reads with alternate allele (nALT/nSumAll)">`},
	{tagVAF1, "VAF1", `##FORMAT=<ID=VAF1,Number=1,Type=Float,Description="The fraction of reads with alternate alleles (nSumALT/nSumAll)">`},
}

// runBuiltinFillTags is the native `+fill-tags` entry point invoked by
// RunPlugin. It parses the upstream plugin flags (-t/--tags is the only
// one that changes the computed set; -d/--drop-missing toggles the
// half-call rule) and the input positional, then fills the tags.
func runBuiltinFillTags(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	valued := map[string]bool{
		"-t": true, "--tags": true,
		"-o": true, "--output": true,
		"-O": true, "--output-type": true,
		"-S": true, "--samples-file": true,
		"-r": true, "--regions": true,
		"-R": true, "--regions-file": true,
	}
	bools := map[string]bool{
		"-d": true, "--drop-missing": true,
	}
	pa, err := parsePluginArgs(opts.Args, valued, bools)
	if err != nil {
		return fmt.Errorf("bcftools +fill-tags: %w", err)
	}
	input := pa.input
	if input == "" {
		input = opts.InputFile
	}
	if input == "" {
		input = "-"
	}
	dropMissing := pa.flags["-d"] || pa.flags["--drop-missing"]

	tags := tagAll
	if t := pa.optVal("-t", "--tags"); t != "" {
		tags, err = parseFillTags(t)
		if err != nil {
			return err
		}
	}

	in, err := iohelper.OpenReader(input)
	if err != nil {
		return fmt.Errorf("bcftools +fill-tags: open %s: %w", input, err)
	}
	defer in.Close()

	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return fmt.Errorf("bcftools +fill-tags: %w", err)
	}

	outHdr := fillTagsHeader(hdr, tags)
	format := opts.OutputFormat
	w, cleanup, err := openOutput(out, ViewOptions{
		OutputFormat:  format,
		CompressLevel: opts.CompressLevel,
	}, outHdr)
	if err != nil {
		return fmt.Errorf("bcftools +fill-tags: %w", err)
	}
	defer cleanup()
	if err := w.WriteHeader(); err != nil {
		return err
	}
	for _, v := range variants {
		fillTagsRecord(v, hdr.Samples, tags, dropMissing)
		if err := w.Write(v); err != nil {
			return err
		}
	}
	return w.Flush()
}

// parseFillTags parses a comma-separated -t/--tags value into a bitmask.
func parseFillTags(s string) (fillTagsBit, error) {
	var flag fillTagsBit
	for _, raw := range strings.Split(s, ",") {
		t := strings.TrimSpace(raw)
		t = strings.TrimPrefix(t, "INFO/")
		t = strings.TrimPrefix(t, "FORMAT/")
		switch strings.ToLower(t) {
		case "all":
			flag |= tagAll
		case "an":
			flag |= tagAN
		case "ac":
			flag |= tagAC
		case "ns":
			flag |= tagNS
		case "ac_hom":
			flag |= tagACHom
		case "ac_het":
			flag |= tagACHet
		case "ac_hemi":
			flag |= tagACHemi
		case "af":
			flag |= tagAF
		case "maf":
			flag |= tagMAF
		case "hwe":
			flag |= tagHWE
		case "exchet":
			flag |= tagExcHet
		case "f_missing":
			flag |= tagFMissing
		case "vaf":
			flag |= tagVAF
		case "vaf1":
			flag |= tagVAF1
		default:
			return 0, fmt.Errorf("bcftools +fill-tags: tag %q not supported", raw)
		}
	}
	return flag, nil
}

// fillTagsHeader returns a copy of hdr with the requested tag ##INFO /
// ##FORMAT lines appended (in upstream order, skipping IDs already
// declared).
func fillTagsHeader(hdr *vcf.Header, tags fillTagsBit) *vcf.Header {
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	existing := map[string]bool{}
	for _, m := range hdr.MetaInfo {
		if id := metaInfoID(m); id != "" {
			existing[id] = true
		}
	}
	for _, h := range fillTagsHeaderLines {
		if tags&h.bit == 0 {
			continue
		}
		if existing[h.id] {
			continue
		}
		out.MetaInfo = append(out.MetaInfo, h.line)
		existing[h.id] = true
	}
	return out
}

// metaInfoID extracts the ID=... value from an ##INFO/##FORMAT line.
func metaInfoID(line string) string {
	if !strings.HasPrefix(line, "##INFO=<") && !strings.HasPrefix(line, "##FORMAT=<") {
		return ""
	}
	i := strings.Index(line, "ID=")
	if i < 0 {
		return ""
	}
	rest := line[i+3:]
	end := strings.IndexAny(rest, ",>")
	if end < 0 {
		return rest
	}
	return rest[:end]
}

// fillCounts accumulates per-allele genotype counts for one record,
// mirroring counts_t in fill-tags.c.
type fillCounts struct {
	nhom  int // homozygous allele count (incremented by 2)
	nhet  int // heterozygous allele count
	nhemi int // hemizygous allele count
	nac   int // half-missing allele count (drop-missing mode)
}

// fillTagsRecord computes and writes the requested INFO tags for one
// record. It re-derives the per-allele counts from FORMAT/GT exactly as
// process_fmt does, then sets each tag in INFO (appending new keys to
// the record's InfoOrder so the writer emits them in upstream order).
func fillTagsRecord(v *vcf.Variant, samples []string, tags fillTagsBit, dropMissing bool) {
	nAllele := 1 + len(v.Alt) // REF + ALTs (Alt always non-empty; "." => monoallelic)
	if len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "") {
		nAllele = 1
	}

	hasGT := false
	for _, f := range v.Format {
		if f == "GT" {
			hasGT = true
			break
		}
	}

	counts := make([]fillCounts, nAllele)
	ns := 0
	nMissingSamples := 0
	nSamples := len(samples)

	if hasGT {
		for i := range v.Samples {
			gt := v.Samples[i].Data["GT"]
			alleles, ploidy := parseGTAlleles(gt)
			if len(alleles) == 0 {
				nMissingSamples++
				continue
			}
			// distinct allele indices
			seen := map[int]bool{}
			for _, a := range alleles {
				if a >= 0 && a < nAllele {
					seen[a] = true
				}
			}
			nals := len(alleles)
			isHom := len(seen) == 1
			var isHalf, isHemi bool
			switch {
			case nals != ploidy:
				if dropMissing {
					isHalf, isHemi = true, false
				} else {
					isHalf, isHemi = false, true
				}
			case nals == 1:
				isHemi = true
			}
			for a := range seen {
				switch {
				case isHalf:
					counts[a].nac++
				case !isHom:
					counts[a].nhet++
				case !isHemi:
					counts[a].nhom += 2
				default:
					counts[a].nhemi++
				}
			}
			ns++
		}
	}

	if tags&tagFMissing != 0 {
		f := 0.0
		if nSamples > 0 {
			f = float64(nMissingSamples) / float64(nSamples)
		}
		setInfo(v, "F_MISSING", formatFloat32G(f))
	}
	if tags&tagNS != 0 {
		setInfo(v, "NS", fmt.Sprintf("%d", ns))
	}

	// Per-allele totals: AC[j] = nhet+nhom+nhemi+nac.
	ac := make([]int, nAllele)
	an := 0
	for j := 0; j < nAllele; j++ {
		ac[j] = counts[j].nhet + counts[j].nhom + counts[j].nhemi + counts[j].nac
		an += ac[j]
	}

	if tags&tagAN != 0 {
		setInfo(v, "AN", fmt.Sprintf("%d", an))
	}
	if tags&tagAC != 0 {
		setInfo(v, "AC", perAlleleInts(ac[1:]))
	}
	if tags&tagAF != 0 {
		af := make([]string, 0, nAllele-1)
		for j := 1; j < nAllele; j++ {
			if an > 0 {
				af = append(af, formatFloat32G(float64(ac[j])/float64(an)))
			} else {
				af = append(af, ".")
			}
		}
		setInfo(v, "AF", strings.Join(af, ","))
	}
	if tags&tagMAF != 0 && nAllele > 1 {
		// Second most common allele frequency.
		freqs := make([]float64, nAllele)
		if an > 0 {
			for j := 0; j < nAllele; j++ {
				freqs[j] = float64(ac[j]) / float64(an)
			}
			sort.Sort(sort.Reverse(sort.Float64Slice(freqs)))
			setInfo(v, "MAF", formatFloat32G(freqs[1]))
		} else {
			setInfo(v, "MAF", ".")
		}
	}
	if tags&tagACHet != 0 {
		vals := make([]int, nAllele-1)
		for j := 1; j < nAllele; j++ {
			vals[j-1] = counts[j].nhet
		}
		setInfo(v, "AC_Het", perAlleleInts(vals))
	}
	if tags&tagACHom != 0 {
		vals := make([]int, nAllele-1)
		for j := 1; j < nAllele; j++ {
			vals[j-1] = counts[j].nhom
		}
		setInfo(v, "AC_Hom", perAlleleInts(vals))
	}
	if tags&tagACHemi != 0 && nAllele > 1 {
		vals := make([]int, nAllele-1)
		for j := 1; j < nAllele; j++ {
			vals[j-1] = counts[j].nhemi
		}
		setInfo(v, "AC_Hemi", perAlleleInts(vals))
	}
	if tags&(tagHWE|tagExcHet) != 0 {
		hwe := make([]string, nAllele-1)
		exc := make([]string, nAllele-1)
		nrefTot := counts[0].nhom
		for j := 0; j < nAllele; j++ {
			nrefTot += counts[j].nhet
		}
		for j := 1; j < nAllele; j++ {
			nref := nrefTot - counts[j].nhet
			nalt := counts[j].nhet + counts[j].nhom
			nhet := counts[j].nhet
			var pHWE, pExc float64 = 1, 1
			if nref > 0 && nalt > 0 {
				pHWE, pExc = calcHWE(nref, nalt, nhet)
			}
			hwe[j-1] = formatFloat32G(pHWE)
			exc[j-1] = formatFloat32G(pExc)
		}
		if tags&tagHWE != 0 && nAllele > 1 {
			setInfo(v, "HWE", strings.Join(hwe, ","))
		}
		if tags&tagExcHet != 0 && nAllele > 1 {
			setInfo(v, "ExcHet", strings.Join(exc, ","))
		}
	}
}

// parseGTAlleles parses a FORMAT/GT string into the list of called
// (non-missing) allele indices and the ploidy (number of allele slots,
// counting missing as a slot). A fully-missing or empty GT returns an
// empty allele slice.
func parseGTAlleles(gt string) (alleles []int, ploidy int) {
	if gt == "" || gt == "." {
		return nil, 0
	}
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	ploidy = len(parts)
	for _, p := range parts {
		if p == "." {
			continue // missing allele, still a slot
		}
		var idx int
		if _, err := fmt.Sscanf(p, "%d", &idx); err != nil {
			return nil, 0
		}
		alleles = append(alleles, idx)
	}
	return alleles, ploidy
}

// perAlleleInts renders a comma-separated list of integers.
func perAlleleInts(vals []int) string {
	if len(vals) == 0 {
		return "."
	}
	parts := make([]string, len(vals))
	for i, v := range vals {
		parts[i] = fmt.Sprintf("%d", v)
	}
	return strings.Join(parts, ",")
}

// calcHWE is the Hardy-Weinberg exact test from fill-tags.c:calc_hwe.
// It returns (HWE p-value, excess-heterozygosity p-value). nref/nalt are
// the reference/alternate allele counts and nhet the heterozygote count.
func calcHWE(nref, nalt, nhet int) (pHWE, pExcHet float64) {
	ngt := (nref + nalt) / 2
	nrare := nref
	if nalt < nref {
		nrare = nalt
	}

	probs := make([]float64, nrare+1)

	mid := int(float64(nrare) * float64(nref+nalt-nrare) / float64(nref+nalt))
	if (nrare&1)^(mid&1) != 0 {
		mid++
	}

	het := mid
	homR := (nrare - mid) / 2
	homC := ngt - het - homR
	probs[mid] = 1.0
	sum := 1.0

	for het = mid; het > 1; het -= 2 {
		probs[het-2] = probs[het] * float64(het) * float64(het-1) /
			(4.0 * float64(homR+1) * float64(homC+1))
		sum += probs[het-2]
		homR++
		homC++
	}

	het = mid
	homR = (nrare - mid) / 2
	homC = ngt - het - homR
	for het = mid; het <= nrare-2; het += 2 {
		probs[het+2] = probs[het] * 4.0 * float64(homR) * float64(homC) /
			(float64(het+2) * float64(het+1))
		sum += probs[het+2]
		homR--
		homC--
	}

	for i := 0; i <= nrare; i++ {
		probs[i] /= sum
	}

	prob := probs[nhet]
	for h := nhet + 1; h <= nrare; h++ {
		prob += probs[h]
	}
	pExcHet = prob

	prob = 0
	for h := 0; h <= nrare; h++ {
		if probs[h] > probs[nhet] {
			continue
		}
		prob += probs[h]
	}
	if prob > 1 {
		prob = 1
	}
	pHWE = prob
	return pHWE, pExcHet
}
