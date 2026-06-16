package bcftools

// TSV2VCF and GVCF2VCF conversion modes for `bcftools convert`, ported
// from upstream's reference_code/bcftools/vcfconvert.c (the tsv_to_vcf /
// gvcf_to_vcf paths) and reference_code/bcftools/tsv2vcf.c (the column
// parser / setters). Both modes require a faidx-indexed reference
// (-f/--fasta-ref): tsv2vcf fills REF (and the contig header) from it, and
// gvcf2vcf fills the per-position REF base when expanding reference blocks.

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// TSV2VCFOptions controls TSV-to-VCF conversion (`bcftools convert
// --tsv2vcf`). It mirrors the subset of upstream args_t fields consulted
// by tsv_to_vcf.
type TSV2VCFOptions struct {
	// FastaRef is the path to the faidx-indexed reference (-f/--fasta-ref).
	// It is mandatory: contigs come from its .fai and REF bases are fetched
	// from it. Upstream errors out when it is absent.
	FastaRef string
	// Columns is the -c/--columns specification (a comma-separated list of
	// column names). When empty, upstream's default "ID,CHROM,POS,AA" is
	// used.
	Columns string
	// Samples / SamplesFile name the genotype columns for the AA (allele
	// array) path. When set, each sample contributes one TSV genotype field
	// per row. SamplesFile takes one name per line.
	Samples     []string
	SamplesFile string
	// KeepDuplicates mirrors upstream's --keep-duplicates flag. Upstream's
	// tsv_to_vcf never consults it (it only affects GEN/SAMPLE export), so it
	// is accepted here for flag parity and has no effect on the output.
	KeepDuplicates bool
	// OutputFormat / CompressLevel select the output container.
	OutputFormat  OutputFormat
	CompressLevel int
	// NoVersion suppresses the ##bcftools_convert{Version,Command} header
	// lines (upstream --no-version). Tests rely on this for deterministic
	// output.
	NoVersion bool
}

// tsvConvert is the shared mutable state threaded through the column
// setters, analogous to upstream's args_t. It is reset per row.
type tsvConvert struct {
	ref      *fasta.RandomAccess
	contigs  map[string]int // contig name → header index (for CHROM validation)
	nSamples int

	// Per-row accumulators.
	ref0  string // REF text for the REF/ALT path
	alt0  string // ALT text for the REF/ALT path
	revAl bool   // REF/ALT verification swapped the alleles
}

// tsvSetter mutates a partially-built variant from one whitespace-delimited
// TSV field. ss..se is the current field. A non-nil error aborts the whole
// run; (skip=true, nil) silently drops the row (matching upstream's ret==-1
// "skipped" path). The remaining-fields slice lets the AA setter consume the
// per-sample genotype columns that trail it.
type tsvSetter func(c *tsvConvert, v *vcf.Variant, field string, rest []string) (skip bool, err error)

// TSV2VCFFile opens the named TSV file and converts it to VCF/BCF on out.
// The reference (opts.FastaRef) must carry a sibling .fai (or be plain so
// one can be built). Returns the number of sites written.
func TSV2VCFFile(path string, out io.Writer, opts TSV2VCFOptions) (int, error) {
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: open %s: %w", path, err)
	}
	defer r.Close()
	return TSV2VCF(r, out, opts)
}

// TSV2VCF streams a TSV source through opts and writes the requested format
// to out. The reference is loaded from opts.FastaRef.
func TSV2VCF(in io.Reader, out io.Writer, opts TSV2VCFOptions) (int, error) {
	if opts.FastaRef == "" {
		return 0, fmt.Errorf("bcftools convert: --tsv2vcf requires the --fasta-ref option")
	}

	ref, err := fasta.OpenRandomAccess(opts.FastaRef)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: could not load the reference %s: %w", opts.FastaRef, err)
	}
	defer ref.Close()

	if opts.SamplesFile != "" {
		names, err := LoadSamplesFile(opts.SamplesFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools convert: %w", err)
		}
		opts.Samples = append(opts.Samples, names...)
	}

	hdr, contigs := buildTSVHeader(ref, opts)

	cv := &tsvConvert{
		ref:      ref,
		contigs:  contigs,
		nSamples: len(opts.Samples),
	}

	cols := opts.Columns
	if cols == "" {
		cols = "ID,CHROM,POS,AA"
	}
	setters, err := buildTSVSetters(cols, opts)
	if err != nil {
		return 0, err
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
	}, hdr)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: %w", err)
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	written := 0
	sc := bufio.NewScanner(in)
	sc.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || line[0] == '#' {
			continue // skip comments / blank lines
		}
		fields := splitTSVFields(line)
		v, skip, err := parseTSVRow(cv, setters, fields, opts)
		if err != nil {
			return written, err
		}
		if skip {
			continue
		}
		if err := w.Write(v); err != nil {
			return written, err
		}
		written++
	}
	if err := sc.Err(); err != nil {
		return written, fmt.Errorf("bcftools convert: read TSV: %w", err)
	}
	return written, w.Flush()
}

// splitTSVFields splits a TSV line on runs of ASCII whitespace, matching
// upstream's isspace_c-based field walk in tsv_parse.
func splitTSVFields(line string) []string {
	return strings.Fields(line)
}

// buildTSVHeader constructs the output header the same way upstream's
// tsv_to_vcf does: fileformat, PASS filter, one ##contig per reference
// sequence, a FORMAT/GT line, an optional version block, and the sample
// list. The returned map keys contig names to their order for CHROM
// validation.
func buildTSVHeader(ref *fasta.RandomAccess, opts TSV2VCFOptions) (*vcf.Header, map[string]int) {
	meta := []string{
		"##fileformat=VCFv4.2",
		`##FILTER=<ID=PASS,Description="All filters passed">`,
	}
	contigs := make(map[string]int)
	for i, e := range ref.Index().Entries() {
		meta = append(meta, fmt.Sprintf("##contig=<ID=%s,length=%d>", e.Name, e.Length))
		contigs[e.Name] = i
	}
	meta = append(meta, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
	if !opts.NoVersion {
		meta = append(meta,
			"##bcftools_convertVersion=bio_ai_experiment",
			"##bcftools_convertCommand=convert --tsv2vcf",
		)
	}
	return &vcf.Header{MetaInfo: meta, Samples: append([]string{}, opts.Samples...)}, contigs
}

// buildTSVSetters resolves the column specification into an ordered list of
// setters. It enforces the same column requirements upstream does in
// tsv_to_vcf: CHROM and POS are mandatory; ID is required unless -c was
// given; and either AA (with -s/-S) or both REF and ALT must be present.
func buildTSVSetters(cols string, opts TSV2VCFOptions) ([]tsvSetter, error) {
	names := strings.Split(cols, ",")
	have := make(map[string]bool)
	setters := make([]tsvSetter, len(names))
	for i, raw := range names {
		name := strings.ToUpper(strings.TrimSpace(raw))
		have[name] = true
		switch name {
		case "CHROM":
			setters[i] = tsvSetChrom
		case "POS":
			setters[i] = tsvSetPos
		case "ID":
			setters[i] = tsvSetID
		case "REF":
			setters[i] = tsvSetRef
		case "ALT":
			setters[i] = tsvSetAlt
		case "AA":
			setters[i] = tsvSetAA
		case "-", "":
			setters[i] = nil // ignored column
		default:
			return nil, fmt.Errorf("bcftools convert: unsupported --columns field %q", raw)
		}
	}

	if !have["CHROM"] {
		return nil, fmt.Errorf("bcftools convert: expected CHROM column")
	}
	if !have["POS"] {
		return nil, fmt.Errorf("bcftools convert: expected POS column")
	}
	if !have["ID"] && opts.Columns == "" {
		return nil, fmt.Errorf("bcftools convert: expected ID column")
	}
	if !have["AA"] {
		if len(opts.Samples) > 0 {
			return nil, fmt.Errorf("bcftools convert: expected AA column with -s/-S")
		}
		if !have["REF"] || !have["ALT"] {
			return nil, fmt.Errorf("bcftools convert: expected REF and ALT columns when AA was not given")
		}
	}
	return setters, nil
}

// parseTSVRow applies every setter to its column and returns the assembled
// variant, or skip=true when a setter requests the row be dropped.
func parseTSVRow(c *tsvConvert, setters []tsvSetter, fields []string, opts TSV2VCFOptions) (*vcf.Variant, bool, error) {
	v := &vcf.Variant{
		ID:     ".",
		Qual:   -1,
		Info:   map[string]string{},
		Filter: nil,
	}
	// Reset per-row REF/ALT accumulators.
	c.ref0, c.alt0, c.revAl = "", "", false

	applied := false
	for i, set := range setters {
		if set == nil {
			continue
		}
		if i >= len(fields) {
			return nil, false, fmt.Errorf("bcftools convert: too few columns in TSV row")
		}
		// The AA setter consumes the per-sample genotype columns that
		// trail it, so it receives the remaining slice.
		rest := fields[i+1:]
		skip, err := set(c, v, fields[i], rest)
		if err != nil {
			return nil, false, err
		}
		if skip {
			return nil, true, nil
		}
		applied = true
	}
	if !applied {
		return nil, true, nil
	}
	return v, false, nil
}

// --- column setters (one per upstream tsv_setter_*) ---

// tsvSetChrom validates the CHROM field against the reference contigs and
// records it. Unknown contigs abort (upstream returns -1 from
// tsv_setter_chrom, which propagates as an error from tsv_parse and ends the
// run).
func tsvSetChrom(c *tsvConvert, v *vcf.Variant, field string, _ []string) (bool, error) {
	if _, ok := c.contigs[field]; !ok {
		return false, fmt.Errorf("bcftools convert: CHROM %q not found in reference", field)
	}
	v.Chrom = field
	return false, nil
}

// tsvSetPos parses the 1-based POS field.
func tsvSetPos(_ *tsvConvert, v *vcf.Variant, field string, _ []string) (bool, error) {
	pos, err := strconv.Atoi(field)
	if err != nil {
		return false, fmt.Errorf("bcftools convert: could not parse POS: %s", field)
	}
	v.Pos = pos
	return false, nil
}

// tsvSetID records the ID column verbatim.
func tsvSetID(_ *tsvConvert, v *vcf.Variant, field string, _ []string) (bool, error) {
	v.ID = field
	return false, nil
}

// tsvSetRef stores the REF text; alleles are finalised once both REF and
// ALT have been seen (mirrors upstream's tsv_setter_ref / _set_ref_alt).
func tsvSetRef(c *tsvConvert, v *vcf.Variant, field string, _ []string) (bool, error) {
	c.ref0 = field
	if c.alt0 != "" {
		setRefAlt(c, v)
	}
	return false, nil
}

// tsvSetAlt stores the ALT text and finalises the alleles if REF is known.
func tsvSetAlt(c *tsvConvert, v *vcf.Variant, field string, _ []string) (bool, error) {
	c.alt0 = field
	if c.ref0 != "" {
		setRefAlt(c, v)
	}
	return false, nil
}

// setRefAlt assembles REF and ALT into the variant. A "." ALT or an ALT
// equal to REF collapses to no alternate allele, matching upstream's
// _set_ref_alt.
func setRefAlt(c *tsvConvert, v *vcf.Variant) {
	v.Ref = c.ref0
	if c.alt0 != "." && c.alt0 != c.ref0 {
		v.Alt = []string{c.alt0}
	} else {
		// No alternate allele: VCF renders ALT as "." (matching htslib's
		// bcf_update_alleles_str with a single REF allele).
		v.Alt = []string{"."}
	}
	c.ref0, c.alt0 = "", ""
}

// acgtTo5 maps a base to 0..3 for A,C,G,T and 4 for anything else, matching
// upstream's acgt_to_5.
func acgtTo5(b byte) int {
	switch b {
	case 'A':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	case 'T':
		return 3
	}
	return 4
}

// tsvSetAA implements upstream's tsv_setter_aa: the field plus the
// per-sample genotype columns are decoded into a per-position allele table
// and the variant's REF/ALT and GT samples are filled from the reference.
// Non-SNP rows (indels, multi-char fields) cause the row to be skipped.
func tsvSetAA(c *tsvConvert, v *vcf.Variant, field string, rest []string) (bool, error) {
	if v.Chrom == "" {
		return false, fmt.Errorf("bcftools convert: AA column requires CHROM before it")
	}
	// Fetch the single reference base at this position (0-based pos-1).
	refBases, err := c.ref.Fetch(v.Chrom, int64(v.Pos-1), int64(v.Pos))
	if err != nil || len(refBases) == 0 {
		return false, fmt.Errorf("bcftools convert: faidx fetch failed at %s:%d", v.Chrom, v.Pos)
	}
	refBase := upperByte(refBases[0])
	iref := acgtTo5(refBase)

	alleles := [5]int{-1, -1, -1, -1, -1}
	alleles[iref] = 0
	nals := 1

	v.Format = []string{"GT"}
	v.Samples = make([]vcf.Sample, c.nSamples)

	for i := 0; i < c.nSamples; i++ {
		var gtField string
		if i == 0 {
			gtField = field
		} else {
			idx := i - 1
			if idx >= len(rest) {
				return false, fmt.Errorf("bcftools convert: too few columns for %d samples at %s:%d", c.nSamples, v.Chrom, v.Pos)
			}
			gtField = rest[idx]
		}
		gt, skip, err := decodeAA1(gtField, &alleles, &nals, iref)
		if err != nil {
			return false, err
		}
		if skip {
			return true, nil // non-SNP: skip whole row (upstream ret==-2)
		}
		v.Samples[i] = vcf.Sample{Name: "", Data: map[string]string{"GT": gt}}
	}

	// Build REF,ALT from the allele table in ACGTN order.
	var sb strings.Builder
	sb.WriteByte(refBase)
	var alt []string
	for i := 0; i < 5; i++ {
		if alleles[i] > 0 {
			alt = append(alt, string("ACGTN"[i]))
		}
	}
	v.Ref = sb.String()
	if len(alt) == 0 {
		// All genotypes were hom-ref: VCF renders ALT as "." (htslib emits a
		// single REF allele).
		alt = []string{"."}
	}
	v.Alt = alt
	return false, nil
}

// decodeAA1 decodes a single sample's two-character allele field into a
// VCF GT string, registering novel alleles in the shared table. It mirrors
// upstream's tsv_setter_aa1: '-'/'.' is missing, 'I'/'D' (indels) request a
// row skip, and >2-char fields are an error.
func decodeAA1(field string, alleles *[5]int, nals *int, iref int) (gt string, skip bool, err error) {
	if len(field) > 2 {
		return "", false, fmt.Errorf("bcftools convert: expected two characters, got %q", field)
	}
	if field == "" {
		return "", false, fmt.Errorf("bcftools convert: empty genotype field")
	}
	switch field[0] {
	case '-', '.':
		return "./.", false, nil
	case 'I', 'D':
		return "", true, nil // insertion/deletion: skip the row
	}
	a0 := acgtTo5(upperByte(field[0]))
	var a1 int
	hasSecond := len(field) == 2
	if hasSecond {
		a1 = acgtTo5(upperByte(field[1]))
	} else {
		a1 = a0
	}
	if alleles[a0] < 0 {
		alleles[a0] = *nals
		*nals++
	}
	if alleles[a1] < 0 {
		alleles[a1] = *nals
		*nals++
	}
	if hasSecond {
		return fmt.Sprintf("%d/%d", alleles[a0], alleles[a1]), false, nil
	}
	// Single-character field: haploid genotype (upstream emits the vector
	// end after the first allele, which renders as a single-allele GT).
	return strconv.Itoa(alleles[a0]), false, nil
}

// GVCFToVCFOptions controls gVCF-block expansion (`bcftools convert
// --gvcf2vcf`).
type GVCFToVCFOptions struct {
	// FastaRef is the faidx-indexed reference (-f/--fasta-ref), mandatory:
	// each expanded position takes its REF base from it.
	FastaRef string
	// IncludeExpr / ExcludeExpr are the standard -i/-e filter expressions.
	// Records that fail the filter are written verbatim (no expansion),
	// matching upstream gvcf_to_vcf.
	IncludeExpr string
	ExcludeExpr string
	// OutputFormat / CompressLevel select the output container.
	OutputFormat  OutputFormat
	CompressLevel int
	// NoVersion suppresses the appended version header line.
	NoVersion bool
}

// GVCFToVCFFile opens the named gVCF input (VCF/BCF, transparently gzipped)
// and writes the expanded VCF to out. Returns the number of records written.
func GVCFToVCFFile(path string, out io.Writer, opts GVCFToVCFOptions) (int, error) {
	if opts.FastaRef == "" {
		return 0, fmt.Errorf("bcftools convert: --gvcf2vcf requires the --fasta-ref option")
	}
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: open %s: %w", path, err)
	}
	defer r.Close()
	return GVCFToVCF(r, out, opts)
}

// GVCFToVCF reads a gVCF stream from in and writes the block-expanded VCF to
// out, fetching per-position REF bases from opts.FastaRef.
func GVCFToVCF(in io.Reader, out io.Writer, opts GVCFToVCFOptions) (int, error) {
	if opts.FastaRef == "" {
		return 0, fmt.Errorf("bcftools convert: --gvcf2vcf requires the --fasta-ref option")
	}
	ref, err := fasta.OpenRandomAccess(opts.FastaRef)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: could not load the fai index for reference %s: %w", opts.FastaRef, err)
	}
	defer ref.Close()

	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: %w", err)
	}

	if !opts.NoVersion {
		hdr.MetaInfo = appendConvertVersion(hdr.MetaInfo)
	}

	include, exclude, err := compileExpressions(ViewOptions{
		IncludeExpr: opts.IncludeExpr,
		ExcludeExpr: opts.ExcludeExpr,
	}, hdr)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: %w", err)
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
	}, hdr)
	if err != nil {
		return 0, fmt.Errorf("bcftools convert: %w", err)
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	written := 0
	for i, v := range variants {
		end1, isBlock := gvcfBlockEnd(v, variants, i)

		// Filter: a record failing the include/exclude test is written
		// verbatim without expansion (upstream gvcf_to_vcf).
		if include != nil || exclude != nil {
			pass := true
			if include != nil && !include.Eval(v) {
				pass = false
			}
			if exclude != nil && exclude.Eval(v) {
				pass = false
			}
			if !pass {
				if err := w.Write(v); err != nil {
					return written, err
				}
				written++
				continue
			}
		}

		if !isBlock {
			if err := w.Write(v); err != nil {
				return written, err
			}
			written++
			continue
		}

		// Expand the block: drop INFO/END, emit one record per position
		// from pos (1-based) to end1 inclusive, REF taken from the
		// reference at each position.
		deleteInfo(v, "END")
		for pos := v.Pos; pos <= end1; pos++ {
			rec := cloneVariant(v)
			rec.Pos = pos
			base, err := ref.Fetch(rec.Chrom, int64(pos-1), int64(pos))
			if err != nil || len(base) == 0 {
				return written, fmt.Errorf("bcftools convert: faidx fetch failed at %s:%d", rec.Chrom, pos)
			}
			rec.Ref = string(upperByte(base[0]))
			if err := w.Write(rec); err != nil {
				return written, err
			}
			written++
		}
	}
	return written, w.Flush()
}

// gvcfBlockEnd reports whether v is a gVCF reference block and, if so, its
// 1-based inclusive end position. A block has a symbolic / absent ALT
// (<*>, <X>, <NON_REF>, or REF-only) plus an INFO/END value. The lookahead
// over the next record clamps END when it would overlap, matching the
// malformed-gVCF handling in upstream gvcf_next_line. Returns (0, false)
// for non-block records.
func gvcfBlockEnd(v *vcf.Variant, all []*vcf.Variant, idx int) (int, bool) {
	if !isGVCFBlockAllele(v) {
		return 0, false
	}
	endStr, ok := v.Info["END"]
	if !ok {
		return 0, false
	}
	end1, err := strconv.Atoi(endStr)
	if err != nil {
		return 0, false
	}
	// Clamp against the next record on the same contig (malformed gVCF).
	// Upstream compares 0-based coordinates: (END-1) >= peek->pos, where
	// peek->pos is 0-based (== peek.Pos-1 here). It then sets the loop bound
	// to peek->pos (0-based), which as a 1-based inclusive end is
	// peek.Pos-1; or to 0 when the block does not start before the peek.
	if idx+1 < len(all) {
		peek := all[idx+1]
		if peek.Chrom == v.Chrom && end1-1 >= peek.Pos-1 {
			if v.Pos < peek.Pos {
				end1 = peek.Pos - 1
			} else {
				end1 = 0
			}
		}
	}
	// When end1 < pos (e.g. clamped to 0 by the overlap guard) the
	// expansion loop runs zero times — the block is dropped entirely,
	// matching upstream which has already removed INFO/END and then writes
	// nothing for pos..end1.
	return end1, true
}

// isGVCFBlockAllele reports whether v's ALT marks a gVCF reference block:
// REF-only (no ALT), or a symbolic ALT of <*>, <X>, or <NON_REF>. Mirrors
// the gallele search in upstream gvcf_next_line.
func isGVCFBlockAllele(v *vcf.Variant) bool {
	// REF-only record (n_allele==1): ALT empty or ".".
	if len(v.Alt) == 0 || (len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "")) {
		return true
	}
	// Upstream only inspects the allele list when the *first* ALT is
	// symbolic (allele[1][0]=='<'), then scans for <*>, <X>, or <NON_REF>.
	if len(v.Alt[0]) == 0 || v.Alt[0][0] != '<' {
		return false
	}
	for _, a := range v.Alt {
		if a == "<*>" || a == "<X>" || a == "<NON_REF>" {
			return true
		}
	}
	return false
}

// deleteInfo removes key from the variant's INFO map and order slice.
func deleteInfo(v *vcf.Variant, key string) {
	if _, ok := v.Info[key]; !ok {
		return
	}
	delete(v.Info, key)
	out := v.InfoOrder[:0]
	for _, k := range v.InfoOrder {
		if k != key {
			out = append(out, k)
		}
	}
	v.InfoOrder = out
}

// appendConvertVersion appends the ##bcftools_convert{Version,Command}
// header lines after the last ## meta line, matching upstream's
// bcf_hdr_append_version placement (just before the #CHROM line).
func appendConvertVersion(meta []string) []string {
	return append(meta,
		"##bcftools_convertVersion=bio_ai_experiment",
		"##bcftools_convertCommand=convert --gvcf2vcf",
	)
}
