package bcftools

// This file implements the Oxford GEN/sample ("gensample") family of the
// `bcftools convert` subcommand: VCF -> .gen+.sample (-g/--gensample) and the
// reverse .gen+.sample -> VCF (-G/--gensample2vcf). It mirrors the
// vcf_to_gensample / gensample_to_vcf paths of upstream
// reference_code/bcftools/vcfconvert.c plus the probability converters in
// convert.c (process_gt_to_prob3 / process_pl_to_prob3 / process_gp_to_prob3)
// and init_sample2sex.
//
// The .gen format is the IMPUTE2 / SNPTEST genotype layout. Each row carries,
// in order:
//
//	[CHROM]  CHROM:POS_REF_ALT  ID  POS  REF  ALT  <3 probs per sample...>
//
// The leading CHROM column is present only with --3N6 (the "3*N+6" variant);
// otherwise the row has 3*N+5 columns. The ID column holds the VCF ID when
// --vcf-ids is given and a copy of CHROM:POS_REF_ALT otherwise.

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// GenSampleOptions configures the VCF<->GEN/sample conversion modes.
type GenSampleOptions struct {
	// Tag selects the FORMAT field that drives the genotype probabilities
	// written to the .gen file. Supported values are "GT" (the default when
	// empty), "PL", and "GP". It only applies to the VCF -> GEN direction.
	Tag string

	// ThreeN6 selects the 3*N+6 column .gen layout, which prepends a bare
	// CHROM column. It maps to upstream's --3N6 flag.
	ThreeN6 bool

	// VCFIDs writes VCF IDs in the .gen ID column (and, on import, reads them
	// from there) instead of repeating the CHROM:POS_REF_ALT label. It maps
	// to upstream's --vcf-ids flag.
	VCFIDs bool

	// SexFile names a "sample<TAB>[MF]" file. When set, a fourth "sex" column
	// is added to the .sample file. It maps to upstream's --sex flag.
	SexFile string

	// KeepDuplicates retains records that share CHROM/POS with the previous
	// record instead of dropping them (upstream --keep-duplicates).
	KeepDuplicates bool

	// IncludeExpr / ExcludeExpr are the standard bcftools -i / -e filter
	// expressions, applied before a record is written to the .gen file.
	IncludeExpr string
	ExcludeExpr string

	// OutputFormat / CompressLevel forward the -O / -l options to the GEN ->
	// VCF writer.
	OutputFormat  OutputFormat
	CompressLevel int

	// NoVersion suppresses the ##bcftools_convert* provenance lines in the
	// GEN -> VCF header (upstream --no-version). v1 never emits those lines,
	// so this is accepted as a no-op for flag parity.
	NoVersion bool
}

// VCFToGenSampleFile reads the VCF/BCF at inPath and writes the Oxford
// .gen(.gz) and .sample files described by prefix. prefix follows upstream's
// --gensample argument grammar: a single token P is expanded to "P.gen.gz" and
// "P.samples", while "GEN,SAMPLE" names the two files explicitly (a "." or
// empty side skips that file).
func VCFToGenSampleFile(inPath, prefix string, opts GenSampleOptions) error {
	r, err := iohelper.OpenReader(inPath)
	if err != nil {
		return fmt.Errorf("bcftools convert: open %s: %w", inPath, err)
	}
	defer r.Close()
	return VCFToGenSample(r, prefix, opts)
}

// VCFToGenSample streams a VCF/BCF source and writes the .gen / .sample files
// named by prefix. It is the streaming entry point split out so tests can drive
// it without a backing input file.
func VCFToGenSample(in io.Reader, prefix string, opts GenSampleOptions) error {
	genFname, sampleFname, err := parseGenSamplePrefix(prefix, true)
	if err != nil {
		return err
	}

	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return err
	}

	tag, err := normalizeGenTag(opts.Tag)
	if err != nil {
		return err
	}

	if sampleFname != "" {
		if err := writeSampleFile(sampleFname, hdr.Samples, opts.SexFile); err != nil {
			return err
		}
	}
	if genFname == "" {
		return nil
	}

	include, exclude, err := compileExpressions(ViewOptions{
		IncludeExpr: opts.IncludeExpr,
		ExcludeExpr: opts.ExcludeExpr,
	})
	if err != nil {
		return fmt.Errorf("bcftools convert: %w", err)
	}

	w, closeFn, err := openGenWriter(genFname)
	if err != nil {
		return err
	}

	prevChrom, prevPos := "", -1
	havePrev := false
	for _, v := range variants {
		if include != nil && !include.Eval(v) {
			continue
		}
		if exclude != nil && exclude.Eval(v) {
			continue
		}
		// ALT allele required, and only bi-allelic records (REF + 1 ALT).
		if len(v.Alt) == 0 || (len(v.Alt) == 1 && v.Alt[0] == ".") {
			continue
		}
		if len(v.Alt) > 1 {
			continue
		}
		if !opts.KeepDuplicates && havePrev && v.Chrom == prevChrom && v.Pos == prevPos {
			continue
		}
		prevChrom, prevPos = v.Chrom, v.Pos
		havePrev = true

		line, err := genRowForVariant(v, tag, opts)
		if err != nil {
			closeFn()
			return err
		}
		if _, err := io.WriteString(w, line); err != nil {
			closeFn()
			return err
		}
	}
	return closeFn()
}

// genRowForVariant renders one .gen line (including the trailing newline) for
// variant v under the given tag and options.
func genRowForVariant(v *vcf.Variant, tag string, opts GenSampleOptions) (string, error) {
	var b strings.Builder
	label := fmt.Sprintf("%s:%d_%s_%s", v.Chrom, v.Pos, v.Ref, firstAlt(v))
	if opts.ThreeN6 {
		b.WriteString(v.Chrom)
		b.WriteByte(' ')
	}
	b.WriteString(label)
	b.WriteByte(' ')
	if opts.VCFIDs {
		b.WriteString(v.ID)
	} else {
		b.WriteString(label)
	}
	b.WriteByte(' ')
	b.WriteString(strconv.Itoa(v.Pos))
	b.WriteByte(' ')
	b.WriteString(v.Ref)
	b.WriteByte(' ')
	b.WriteString(firstAlt(v))

	probs, err := genProbsForVariant(v, tag)
	if err != nil {
		return "", err
	}
	b.WriteString(probs)
	b.WriteByte('\n')
	return b.String(), nil
}

// firstAlt returns the first ALT allele, or "." when there is none. It mirrors
// upstream's %FIRST_ALT conversion key.
func firstAlt(v *vcf.Variant) string {
	if len(v.Alt) == 0 {
		return "."
	}
	return v.Alt[0]
}

// genProbsForVariant produces the space-prefixed 3-per-sample probability
// columns for one record, dispatching on tag.
func genProbsForVariant(v *vcf.Variant, tag string) (string, error) {
	switch tag {
	case "GT":
		return gtToProb3(v)
	case "PL":
		return plToProb3(v)
	case "GP":
		return gpToProb3(v)
	default:
		// Upstream (vcfconvert.c:vcf_to_gensample) supports only GT, PL and
		// GP for the .gen file and rejects everything else — including GL,
		// which the --help text lists but the code never wires up — with the
		// exact message "todo: --tag %s". We mirror that rejection verbatim so
		// the behaviour is positive parity rather than a divergent stub.
		return "", fmt.Errorf("bcftools convert: todo: --tag %s", tag)
	}
}

// gtToProb3 mirrors process_gt_to_prob3: it derives hard {0,1} probability
// triples from the GT field. Het is "0 1 0", first-ALT hom is "0 0 1", any
// other hom is "1 0 0", and missing diploid is "0.33 0.33 0.33". Haploid calls
// collapse to "1 0 0" / "0 0 1" (missing -> "0.5 0.0 0.5").
func gtToProb3(v *vcf.Variant) (string, error) {
	var b strings.Builder
	for i := range v.Samples {
		gt, ok := v.Samples[i].Data["GT"]
		if !ok || gt == "" {
			return "", fmt.Errorf("bcftools convert: error parsing GT tag at %s:%d", v.Chrom, v.Pos)
		}
		alleles := splitGT(gt)
		switch len(alleles) {
		case 2:
			a0, miss0 := parseAllele(alleles[0])
			a1, _ := parseAllele(alleles[1])
			switch {
			case miss0:
				b.WriteString(" 0.33 0.33 0.33")
			case a0 != a1:
				b.WriteString(" 0 1 0")
			case a0 == 1:
				b.WriteString(" 0 0 1")
			default:
				b.WriteString(" 1 0 0")
			}
		case 1:
			a0, miss0 := parseAllele(alleles[0])
			switch {
			case miss0:
				b.WriteString(" 0.5 0.0 0.5")
			case a0 == 1:
				b.WriteString(" 0 0 1")
			default:
				b.WriteString(" 1 0 0")
			}
		default:
			// Upstream process_gt_to_prob3 (convert.c) only handles haploid
			// (ploidy 1) and diploid (ploidy 2) genotypes; any other ploidy
			// aborts with error("FIXME: not ready for ploidy %d\n", j). We
			// reject the same cases with the same message text so the
			// behaviour is positive parity, not a port-specific extension.
			return "", fmt.Errorf("bcftools convert: FIXME: not ready for ploidy %d at %s:%d", len(alleles), v.Chrom, v.Pos)
		}
	}
	return b.String(), nil
}

// plToProb3 mirrors process_pl_to_prob3: it converts phred-scaled PL values to
// normalised linear probabilities. With n_allele genotype likelihoods present
// the record is treated as haploid (two values, REF/ALT, with a zero het).
func plToProb3(v *vcf.Variant) (string, error) {
	var b strings.Builder
	nAllele := 1 + len(v.Alt)
	for i := range v.Samples {
		raw, ok := v.Samples[i].Data["PL"]
		if !ok || raw == "" || raw == "." {
			return "", fmt.Errorf("bcftools convert: error parsing PL tag at %s:%d", v.Chrom, v.Pos)
		}
		parts := strings.Split(raw, ",")
		vals := make([]float64, 0, len(parts))
		sum := 0.0
		for _, p := range parts {
			if p == "." {
				break
			}
			pl, err := strconv.ParseFloat(p, 64)
			if err != nil {
				return "", fmt.Errorf("bcftools convert: error parsing PL tag at %s:%d", v.Chrom, v.Pos)
			}
			lin := math.Pow(10, -0.1*pl)
			vals = append(vals, lin)
			sum += lin
		}
		if len(vals) == nAllele {
			// haploid
			b.WriteByte(' ')
			b.WriteString(formatProbFloat(vals[0] / sum))
			b.WriteString(" 0 ")
			b.WriteString(formatProbFloat(vals[1] / sum))
		} else {
			// diploid
			b.WriteByte(' ')
			b.WriteString(formatProbFloat(vals[0] / sum))
			b.WriteByte(' ')
			b.WriteString(formatProbFloat(vals[1] / sum))
			b.WriteByte(' ')
			b.WriteString(formatProbFloat(vals[2] / sum))
		}
	}
	return b.String(), nil
}

// gpToProb3 mirrors process_gp_to_prob3: it copies the VCF4.3+ GP genotype
// posterior probabilities through, clamping missing entries to 0 and rejecting
// values outside [0,1]. n_allele values denotes a haploid record.
func gpToProb3(v *vcf.Variant) (string, error) {
	var b strings.Builder
	nAllele := 1 + len(v.Alt)
	for i := range v.Samples {
		raw, ok := v.Samples[i].Data["GP"]
		if !ok || raw == "" {
			return "", fmt.Errorf("bcftools convert: error parsing GP tag at %s:%d", v.Chrom, v.Pos)
		}
		parts := strings.Split(raw, ",")
		vals := make([]float64, 0, len(parts))
		for _, p := range parts {
			if p == "." {
				vals = append(vals, 0)
				continue
			}
			gp, err := strconv.ParseFloat(p, 64)
			if err != nil {
				return "", fmt.Errorf("bcftools convert: error parsing GP tag at %s:%d", v.Chrom, v.Pos)
			}
			if gp < 0 || gp > 1 {
				return "", fmt.Errorf("bcftools convert: [%s:%d:%g] GP value outside range [0,1]", v.Chrom, v.Pos, gp)
			}
			vals = append(vals, gp)
		}
		if len(vals) == nAllele {
			// haploid: ptr[0], 0, ptr[1]
			b.WriteByte(' ')
			b.WriteString(formatGPFloat(vals[0]))
			b.WriteByte(' ')
			b.WriteString(formatGPFloat(0))
			b.WriteByte(' ')
			b.WriteString(formatGPFloat(vals[1]))
		} else {
			b.WriteByte(' ')
			b.WriteString(formatGPFloat(vals[0]))
			b.WriteByte(' ')
			b.WriteString(formatGPFloat(vals[1]))
			b.WriteByte(' ')
			b.WriteString(formatGPFloat(vals[2]))
		}
	}
	return b.String(), nil
}

// formatProbFloat renders a probability with upstream's "%f" (6-decimal) C
// formatting, matching ksprintf in process_pl_to_prob3.
func formatProbFloat(f float64) string {
	return strconv.FormatFloat(f, 'f', 6, 64)
}

// formatGPFloat renders a GP probability with upstream's "%f" (6-decimal) C
// formatting, matching ksprintf in process_gp_to_prob3.
func formatGPFloat(f float64) string {
	return strconv.FormatFloat(f, 'f', 6, 64)
}

// writeSampleFile writes the Oxford .sample file for the given samples. The
// header is two lines ("ID_1 ID_2 missing[ sex]" then a "0 0 0[ 0]" type row),
// followed by one "ID ID 0[ SEX]" row per sample. When sexFile is set the
// optional sex column is appended and every sample must have an entry.
func writeSampleFile(path string, samples []string, sexFile string) error {
	var sample2sex map[string]byte
	if sexFile != "" {
		var err error
		sample2sex, err = initSampleToSex(samples, sexFile)
		if err != nil {
			return err
		}
	}

	w, closeFn, err := openGenWriter(path)
	if err != nil {
		return err
	}
	var b strings.Builder
	if sample2sex != nil {
		b.WriteString("ID_1 ID_2 missing sex\n0 0 0 0\n")
	} else {
		b.WriteString("ID_1 ID_2 missing\n0 0 0\n")
	}
	for _, s := range samples {
		if sample2sex != nil {
			fmt.Fprintf(&b, "%s %s 0 %c\n", s, s, sample2sex[s])
		} else {
			fmt.Fprintf(&b, "%s %s 0\n", s, s)
		}
	}
	if _, err := io.WriteString(w, b.String()); err != nil {
		closeFn()
		return err
	}
	return closeFn()
}

// initSampleToSex parses a "sample<whitespace>[MF]" file into a per-sample sex
// code ('1' for M, '2' for F), mirroring upstream's init_sample2sex. Every
// sample in the header must be present, and only M/F are accepted.
func initSampleToSex(samples []string, sexFile string) (map[string]byte, error) {
	f, err := os.Open(sexFile)
	if err != nil {
		return nil, fmt.Errorf("bcftools convert: could not read %s: %w", sexFile, err)
	}
	defer f.Close()

	known := make(map[string]bool, len(samples))
	for _, s := range samples {
		known[s] = true
	}
	out := make(map[string]byte, len(samples))
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := sc.Text()
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}
		name := fields[0]
		if !known[name] {
			continue
		}
		switch fields[1][0] {
		case 'M':
			out[name] = '1'
		case 'F':
			out[name] = '2'
		default:
			return nil, fmt.Errorf("bcftools convert: could not parse %s: %s", sexFile, line)
		}
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("bcftools convert: read %s: %w", sexFile, err)
	}
	for _, s := range samples {
		if _, ok := out[s]; !ok {
			return nil, fmt.Errorf("bcftools convert: missing sex for sample %s in %s", s, sexFile)
		}
	}
	return out, nil
}

// GenSampleToVCFFile reads the .gen / .sample files named by prefix and writes
// VCF/BCF to out. prefix follows the same grammar as the export direction:
// "P" -> "P.gen.gz"+"P.samples", or "GEN,SAMPLE" naming both explicitly.
func GenSampleToVCFFile(prefix string, out io.Writer, opts GenSampleOptions) error {
	genFname, sampleFname, err := parseGenSamplePrefix(prefix, false)
	if err != nil {
		return err
	}
	genR, err := iohelper.OpenReader(genFname)
	if err != nil {
		return fmt.Errorf("bcftools convert: could not read %s: %w", genFname, err)
	}
	defer genR.Close()
	sampleR, err := iohelper.OpenReader(sampleFname)
	if err != nil {
		return fmt.Errorf("bcftools convert: could not read %s: %w", sampleFname, err)
	}
	defer sampleR.Close()
	return GenSampleToVCF(genR, sampleR, out, opts)
}

// GenSampleToVCF reads an Oxford .gen stream and its .sample stream and writes
// the reconstructed VCF/BCF to out. It mirrors gensample_to_vcf: GT is derived
// from the probability triples (argmax, REF-hom on ties) and GP carries the
// raw triples.
func GenSampleToVCF(genIn, sampleIn io.Reader, out io.Writer, opts GenSampleOptions) error {
	samples, err := readGenSampleNames(sampleIn)
	if err != nil {
		return err
	}

	genScanner := bufio.NewScanner(genIn)
	genScanner.Buffer(make([]byte, 0, 64*1024), 64*1024*1024)
	if !genScanner.Scan() {
		if err := genScanner.Err(); err != nil {
			return err
		}
		return fmt.Errorf("bcftools convert: empty .gen file")
	}
	firstLine := genScanner.Text()

	chrom, err := genChromFromLine(firstLine, opts.ThreeN6)
	if err != nil {
		return err
	}

	hdr := &vcf.Header{
		MetaInfo: []string{
			"##fileformat=VCFv4.2",
			`##FILTER=<ID=PASS,Description="All filters passed">`,
			`##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">`,
			`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`,
			`##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">`,
			fmt.Sprintf("##contig=<ID=%s,length=2147483647>", chrom),
		},
		Samples: samples,
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
	}, hdr)
	if err != nil {
		return fmt.Errorf("bcftools convert: %w", err)
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return err
	}

	line := firstLine
	for {
		v, err := parseGenLine(line, samples, opts)
		if err != nil {
			return err
		}
		if err := w.Write(v); err != nil {
			return err
		}
		if !genScanner.Scan() {
			break
		}
		line = genScanner.Text()
	}
	if err := genScanner.Err(); err != nil {
		return err
	}
	return w.Flush()
}

// readGenSampleNames reads the .sample file and returns the per-individual
// names. Like upstream it skips the two header rows (offset 2) and takes the
// first whitespace-delimited token of each remaining row.
func readGenSampleNames(r io.Reader) ([]string, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	var rows []string
	for sc.Scan() {
		rows = append(rows, sc.Text())
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if len(rows) < 2 {
		return nil, fmt.Errorf("bcftools convert: malformed .sample file")
	}
	var samples []string
	for _, row := range rows[2:] {
		fields := strings.Fields(row)
		if len(fields) == 0 {
			continue
		}
		samples = append(samples, fields[0])
	}
	return samples, nil
}

// genChromFromLine extracts the CHROM name from the first .gen line. With
// threeN6 the CHROM is the bare first column; otherwise it is the part before
// the first ':' of the CHROM:POS_REF_ALT label, searched across the first two
// columns (upstream tolerates the label appearing in either column).
func genChromFromLine(line string, threeN6 bool) (string, error) {
	fields := strings.Fields(line)
	if len(fields) == 0 {
		return "", fmt.Errorf("bcftools convert: could not determine CHROM: %s", line)
	}
	if threeN6 {
		return fields[0], nil
	}
	for _, f := range fields[:genMin(2, len(fields))] {
		if idx := strings.IndexByte(f, ':'); idx >= 0 {
			return f[:idx], nil
		}
	}
	return "", fmt.Errorf("bcftools convert: could not determine CHROM: %s", line)
}

// parseGenLine parses one .gen row into a Variant. It validates the embedded
// POS and REF/ALT against the CHROM:POS_REF_ALT label (allowing the row to swap
// REF/ALT order, in which case the probability triple is reversed, matching
// upstream's tsv_setter_verify_ref_alt + rev_als handling).
func parseGenLine(line string, samples []string, opts GenSampleOptions) (*vcf.Variant, error) {
	fields := strings.Fields(line)
	idx := 0
	if opts.ThreeN6 {
		idx++ // skip bare CHROM column
	}
	// Need: label, id, pos, ref, alt + 3*nsamples probabilities.
	want := idx + 5 + 3*len(samples)
	if len(fields) < want {
		return nil, fmt.Errorf("bcftools convert: wrong number of fields: %s", line)
	}

	col1 := fields[idx]
	col2 := fields[idx+1]
	posStr := fields[idx+2]
	ref := fields[idx+3]
	alt := fields[idx+4]

	// The CHROM:POS_REF_ALT label may appear in the first column, or, when an
	// IMPUTE2 reference panel left the first column as "--", in the second.
	// Mirror upstream's tsv_setter_chrom_pos_ref_alt_or_id /
	// tsv_setter_chrom_pos_ref_alt_id_or_die handshake: try column 1 first and
	// fall back to column 2, and source the VCF ID from whichever column did
	// NOT hold the label (only when --vcf-ids is set).
	id := "."
	chrom, pos, lref, lalt, ok := parseChromPosRefAlt(col1)
	if ok {
		if opts.VCFIDs {
			id = col2
		}
	} else {
		if opts.VCFIDs {
			id = col1
		}
		chrom, pos, lref, lalt, ok = parseChromPosRefAlt(col2)
		if !ok {
			return nil, fmt.Errorf("bcftools convert: could not parse the CHROM:POS_REF_ALT[_END] string: %s", col2)
		}
	}

	gotPos, err := strconv.Atoi(posStr)
	if err != nil {
		return nil, fmt.Errorf("bcftools convert: could not parse POS: %s", posStr)
	}
	if gotPos != pos {
		return nil, fmt.Errorf("bcftools convert: POS mismatch: %s", line)
	}

	// Verify REF/ALT and detect a reversed-allele row.
	revAls := false
	switch {
	case ref == lref && alt == lalt:
		revAls = false
	case ref == lalt && alt == lref:
		revAls = true
	default:
		return nil, fmt.Errorf("bcftools convert: REF/ALT mismatch: [%s][%s]", ref, alt)
	}

	v := &vcf.Variant{
		Chrom:   chrom,
		Pos:     pos,
		ID:      id,
		Ref:     lref,
		Alt:     []string{lalt},
		Qual:    -1,
		Filter:  []string{"."},
		Info:    map[string]string{},
		Format:  []string{"GT", "GP"},
		Samples: make([]vcf.Sample, len(samples)),
	}

	probs := fields[idx+5:]
	for i := range samples {
		aa, err := strconv.ParseFloat(probs[3*i+0], 64)
		if err != nil {
			return nil, fmt.Errorf("bcftools convert: could not parse first value of %d-th sample", i+1)
		}
		ab, err := strconv.ParseFloat(probs[3*i+1], 64)
		if err != nil {
			return nil, fmt.Errorf("bcftools convert: could not parse second value of %d-th sample", i+1)
		}
		bb, err := strconv.ParseFloat(probs[3*i+2], 64)
		if err != nil {
			return nil, fmt.Errorf("bcftools convert: could not parse third value of %d-th sample", i+1)
		}
		if revAls {
			aa, bb = bb, aa
		}
		gt := gtFromProb3(aa, ab, bb)
		gp := formatGenGP(aa) + "," + formatGenGP(ab) + "," + formatGenGP(bb)
		v.Samples[i] = vcf.Sample{
			Name: samples[i],
			Data: map[string]string{"GT": gt, "GP": gp},
		}
	}
	return v, nil
}

// gtFromProb3 derives the GT call from a probability triple, mirroring the
// argmax-with-REF-hom-on-ties logic of upstream's tsv_setter_gt_gp.
func gtFromProb3(aa, ab, bb float64) string {
	if aa >= ab {
		if aa >= bb {
			return "0/0"
		}
		return "1/1"
	}
	if ab >= bb {
		return "0/1"
	}
	return "1/1"
}

// formatGenGP renders a probability triple value for the reconstructed GP field
// using the same compact float formatting htslib applies on output (integers
// print bare; otherwise the shortest float32 round-trip representation).
func formatGenGP(f float64) string {
	f32 := float32(f)
	if f32 == float32(int32(f32)) {
		return strconv.FormatInt(int64(f32), 10)
	}
	return strconv.FormatFloat(float64(f32), 'g', -1, 32)
}

// parseChromPosRefAlt splits a "CHROM:POS_REF_ALT[_END]" label into its parts.
// It returns ok=false when the label does not match that shape.
func parseChromPosRefAlt(label string) (chrom string, pos int, ref, alt string, ok bool) {
	ci := strings.IndexByte(label, ':')
	if ci < 0 {
		return "", 0, "", "", false
	}
	chrom = label[:ci]
	rest := label[ci+1:]
	ui := strings.IndexByte(rest, '_')
	if ui < 0 {
		return "", 0, "", "", false
	}
	posStr := rest[:ui]
	pos, err := strconv.Atoi(posStr)
	if err != nil {
		return "", 0, "", "", false
	}
	rest = rest[ui+1:]
	ui = strings.IndexByte(rest, '_')
	if ui < 0 {
		return "", 0, "", "", false
	}
	ref = rest[:ui]
	alt = rest[ui+1:]
	// Drop an optional _END suffix.
	if ei := strings.IndexByte(alt, '_'); ei >= 0 {
		alt = alt[:ei]
	}
	if chrom == "" || ref == "" || alt == "" {
		return "", 0, "", "", false
	}
	return chrom, pos, ref, alt, true
}

// parseGenSamplePrefix expands the --gensample / --gensample2vcf argument into
// the gen and sample filenames. A single token P expands to P.gen.gz +
// P.samples; "GEN,SAMPLE" names both explicitly, with a "." or empty side
// skipped (export only — import requires both files).
func parseGenSamplePrefix(prefix string, export bool) (genFname, sampleFname string, err error) {
	if prefix == "" {
		return "", "", fmt.Errorf("bcftools convert: error parsing --gensample filenames: %s", prefix)
	}
	if idx := strings.IndexByte(prefix, ','); idx >= 0 {
		gen := prefix[:idx]
		smpl := prefix[idx+1:]
		if export {
			if gen == "." {
				gen = ""
			}
			if smpl == "." {
				smpl = ""
			}
			return gen, smpl, nil
		}
		if gen == "" || smpl == "" {
			return "", "", fmt.Errorf("bcftools convert: error parsing --gensample2vcf filenames: %s", prefix)
		}
		return gen, smpl, nil
	}
	return prefix + ".gen.gz", prefix + ".samples", nil
}

// normalizeGenTag validates the --tag value and returns the canonical tag.
// An empty value defaults to GT.
func normalizeGenTag(tag string) (string, error) {
	switch tag {
	case "", "GT":
		return "GT", nil
	case "PL":
		return "PL", nil
	case "GP":
		return "GP", nil
	default:
		// Matches upstream vcfconvert.c, which accepts only GT/PL/GP for the
		// .gen file and otherwise aborts with "todo: --tag %s".
		return "", fmt.Errorf("bcftools convert: todo: --tag %s", tag)
	}
}

// openGenWriter opens path for writing, bgzf-compressing when it ends in .gz so
// the output matches htslib's bgzf_open("wg"). It returns the writer and a
// close function that flushes and closes both the compressor and the file.
func openGenWriter(path string) (io.Writer, func() error, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, nil, fmt.Errorf("bcftools convert: create %s: %w", path, err)
	}
	if strings.HasSuffix(path, ".gz") {
		gw := bgzf.NewWriter(f)
		closeFn := func() error {
			if err := gw.Close(); err != nil {
				f.Close()
				return err
			}
			return f.Close()
		}
		return gw, closeFn, nil
	}
	return f, f.Close, nil
}

// splitGT splits a GT string into its per-allele tokens on either phase
// separator ('/' or '|').
func splitGT(gt string) []string {
	return strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
}

// parseAllele parses a single GT allele token, returning its numeric index and
// whether it is the missing allele ('.').
func parseAllele(tok string) (allele int, missing bool) {
	if tok == "." || tok == "" {
		return 0, true
	}
	n, err := strconv.Atoi(tok)
	if err != nil {
		return 0, true
	}
	return n, false
}

// min2 returns the smaller of a and b. (A local helper to avoid relying on the
// Go 1.21 builtin min, keeping the file compatible with the CI toolchain.)
func genMin(a, b int) int {
	if a < b {
		return a
	}
	return b
}
