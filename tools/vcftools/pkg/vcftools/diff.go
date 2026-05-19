// VCF comparison ("--diff" family) for vcftools.
//
// Supported flags:
//
//	--diff FILE                 second VCF to compare against
//	--diff-site                 emit <prefix>.diff.sites_in_files
//	--diff-indv                 emit <prefix>.diff.indv_in_files
//	--diff-site-discordance     emit <prefix>.diff.sites
//	--diff-indv-discordance     emit <prefix>.diff.indv
//	--diff-indv-map FILE        two-column file that renames file-2 sample
//	                            IDs before matching against file-1 (upstream
//	                            variant_file_diff.cpp:11-34)
//	--diff-discordance-matrix   emit <prefix>.diff.discordance_matrix
//	                            (4x4 genotype-by-genotype counts; upstream
//	                            variant_file_diff.cpp:944)
//
// File-2 (the "--diff" file) is loaded fully into memory keyed by
// (CHR, POS); each variant keeps REF/ALT and its sample genotypes. File-1 is
// streamed through the regular vcftools filter pipeline before being handed
// to this module — that means filters such as --chr, --bed, --minQ apply to
// the file-1 side but the file-2 side is unfiltered (mirroring upstream).
//
// Output formats follow the upstream column layout:
//
//	.diff.sites_in_files       CHROM POS1 POS2 IN_FILE REF1 REF2 ALT1 ALT2
//	                           IN_FILE ∈ {1, 2, B}; "." for absent fields
//
//	.diff.indv_in_files        INDV FILES  (FILES ∈ {1, 2, B})
//
//	.diff.sites                CHROM POS FILES MATCHING_ALLELES
//	                           N_COMMON_CALLED N_DISCORD DISCORDANCE
//	                           Mirrors upstream variant_file_diff.cpp:677
//	                           (the header line). FILES ∈ {1, 2, B} just
//	                           like .diff.sites_in_files; MATCHING_ALLELES
//	                           is 1 when REF and the first ALT match in
//	                           both files (B-rows only) else 0;
//	                           N_COMMON_CALLED / N_DISCORD count common
//	                           samples called in both files; DISCORDANCE
//	                           = N_DISCORD / N_COMMON_CALLED, with
//	                           division by zero rendered as -nan to
//	                           match upstream's C++ printf.
//
//	.diff.indv                 INDV N_COMMON_CALLED N_DISCORD DISCORDANCE
//	                           Lists the *union* of file-1 samples and
//	                           (post-map) file-2 samples in alphabetical
//	                           order, matching upstream's
//	                           combined_individuals std::map iteration
//	                           (variant_file_diff.cpp:619-625). Samples
//	                           not shared between the two files appear
//	                           with 0/0/-nan; sites present in only one
//	                           file contribute 0 to every sample.
//
//	.diff.discordance_matrix   5x5 grid: header row "-\tN_0/0_file1\t...\t
//	                           N_./._file1"; four data rows labelled
//	                           N_<GT>_file2 with four counts each. Counts
//	                           biallelic diploid genotype pairs (ALT must
//	                           match, REF treated as in upstream:
//	                           variant_file_diff.cpp:1072-1083).
//
// Discordance compares unphased, sorted allele indices restricted to the
// first ALT (REF=0, ALT=1, anything else=missing) which is how upstream
// vcftools' diff family treats multi-allelic sites by default.
package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// diffRecord holds the file-2 site information we care about for the diff.
type diffRecord struct {
	ref     string
	alt     string // empty if monomorphic / no ALT
	rawALTs []string
	// genotypes keyed by sample name. Stores the raw GT string.
	genotypes map[string]string
}

// diffData is the loaded second VCF used by --diff.
type diffData struct {
	samples []string
	// sites[chrom][pos] = record
	sites map[string]map[int]*diffRecord
}

// loadDiffVCF reads the second VCF file fully into memory.
func loadDiffVCF(filename string) (*diffData, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, fmt.Errorf("opening --diff file %s: %w", filename, err)
	}
	defer f.Close()

	reader := vcf.NewReader(f)
	hdr, err := reader.ReadHeader()
	if err != nil {
		return nil, fmt.Errorf("reading --diff VCF header: %w", err)
	}
	return loadDiffFromSource(&vcfVariantSource{r: reader, hdr: hdr}, "--diff")
}

// loadDiffBCF reads the --diff-bcf second BCF file fully into memory.
// Mirrors `loadDiffVCF` but uses the shared BCF reader stack
// (BGZF + bcf.Reader → vcf.Variant via Record.ToVariant). The wave-22
// variantSource adapter handles both file ownership and the
// VCF/BCF dispatch uniformly.
func loadDiffBCF(filename string) (*diffData, error) {
	src, err := newBCFVariantSource(filename)
	if err != nil {
		return nil, fmt.Errorf("opening --diff-bcf %s: %w", filename, err)
	}
	defer src.Close()
	return loadDiffFromSource(src, "--diff-bcf")
}

// loadDiffFromSource is the shared body of loadDiffVCF / loadDiffBCF:
// it reads every variant from the source and builds the (CHR,POS)-keyed
// diffData map plus the genotype lookup per sample. flagLabel is the
// user-facing flag name used in error messages.
func loadDiffFromSource(src variantSource, flagLabel string) (*diffData, error) {
	d := &diffData{
		samples: append([]string(nil), src.Header().Samples...),
		sites:   make(map[string]map[int]*diffRecord),
	}
	for {
		v, err := src.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading %s file: %w", flagLabel, err)
		}
		rec := &diffRecord{
			ref:       v.Ref,
			rawALTs:   append([]string(nil), v.Alt...),
			genotypes: make(map[string]string, len(v.Samples)),
		}
		if len(v.Alt) > 0 {
			rec.alt = v.Alt[0]
		}
		for i := range v.Samples {
			if gt, ok := v.Samples[i].Data["GT"]; ok {
				rec.genotypes[v.Samples[i].Name] = gt
			}
		}
		bucket := d.sites[v.Chrom]
		if bucket == nil {
			bucket = make(map[int]*diffRecord)
			d.sites[v.Chrom] = bucket
		}
		// If the same position appears twice in the diff file, the last
		// one wins — matches upstream's "overwrite" behaviour.
		bucket[v.Pos] = rec
	}
	return d, nil
}

// commonPair links a file-1 sample name to its file-2 sample name. When no
// --diff-indv-map is supplied the two names are equal; with a map they may
// differ. The "name" used for per-individual output keying is f1Name (the
// canonical/file-1 ID), matching upstream's combined_individuals map.
type commonPair struct {
	f1Name string
	f2Name string
}

// diffRunner accumulates per-site and per-individual stats while file-1 is
// streamed through Run, then flushes outputs in close().
type diffRunner struct {
	params  *Params
	data    *diffData
	samples []string // file-1 samples (post sample filtering)

	// Set of file-1 samples — quick lookup for the indv_in_files report.
	file1SampleSet map[string]struct{}

	// indvMap stores the parsed --diff-indv-map table: file-2 sample ID →
	// renamed-to ID (typically a file-1 sample name). Empty when the flag
	// isn't supplied. file2RenamedSet records the renamed file-2 IDs so the
	// indv_in_files report classifies them as "B" rather than "2".
	indvMap         map[string]string
	file2RenamedSet map[string]struct{}

	// Set of (chrom,pos) seen in file-1 for the sites_in_files report's
	// "B" classification. Once a site is seen on the file-1 side we mark it
	// here; at close() time we walk file-2's sites and any not in this map
	// are "2"-only.
	seenSites map[string]map[int]struct{}

	// Per-sample discordance accumulators, keyed by the file-1 (canonical)
	// sample name.
	indvCommon  map[string]int
	indvDiscord map[string]int
	// commonPairs is the ordered intersection of file-1 and (mapped) file-2
	// sample names. Iterating it in file-1 order keeps output stable.
	commonPairs []commonPair

	// discMatrix holds the 4x4 genotype-by-genotype counts for
	// --diff-discordance-matrix. Indexed as [file1GT][file2GT] where the
	// genotype code is 0=0/0, 1=0/1, 2=1/1, 3=./. (matching the file2-as-
	// row, file1-as-column layout of the upstream output).
	discMatrix [4][4]int

	// Output writers.
	wSitesInFiles *diffOutFile
	wSites        *diffOutFile

	// switchRun is the per-individual phase-switch tracker for
	// --diff-switch-error. nil when the flag isn't set. See
	// diff_switch.go.
	switchRun *switchRunner
}

// diffOutFile bundles a bufio.Writer with its file handle.
type diffOutFile struct {
	f io.WriteCloser
	w *bufio.Writer
}

func newDiffOutFile(path string, header string) (*diffOutFile, error) {
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", path, err)
	}
	w := bufio.NewWriter(f)
	if header != "" {
		if _, err := w.WriteString(header); err != nil {
			f.Close()
			return nil, fmt.Errorf("writing header to %s: %w", path, err)
		}
	}
	return &diffOutFile{f: f, w: w}, nil
}

func (o *diffOutFile) close() error {
	if o == nil {
		return nil
	}
	var firstErr error
	if o.w != nil {
		if err := o.w.Flush(); err != nil {
			firstErr = err
		}
	}
	if o.f != nil {
		if err := o.f.Close(); err != nil && firstErr == nil {
			firstErr = err
		}
	}
	return firstErr
}

// newDiffRunner prepares per-output-file writers and the intersection sample
// list. It is the caller's responsibility to ensure params.Diff or
// params.DiffBCF is set and at least one --diff-* output flag is set.
func newDiffRunner(params *Params, samples []string) (*diffRunner, error) {
	if params.Diff == "" && params.DiffBCF == "" {
		return nil, nil
	}
	if !params.DiffSite && !params.DiffIndv && !params.DiffSiteDiscordance &&
		!params.DiffIndvDiscordance && !params.DiffDiscordanceMatrix &&
		!params.DiffSwitchError {
		return nil, nil
	}
	var (
		d   *diffData
		err error
	)
	if params.DiffBCF != "" {
		d, err = loadDiffBCF(params.DiffBCF)
	} else {
		d, err = loadDiffVCF(params.Diff)
	}
	if err != nil {
		return nil, err
	}

	r := &diffRunner{
		params:          params,
		data:            d,
		samples:         append([]string(nil), samples...),
		file1SampleSet:  make(map[string]struct{}, len(samples)),
		file2RenamedSet: make(map[string]struct{}),
		seenSites:       make(map[string]map[int]struct{}),
		indvCommon:      make(map[string]int),
		indvDiscord:     make(map[string]int),
	}
	for _, s := range samples {
		r.file1SampleSet[s] = struct{}{}
	}

	// Optional --diff-indv-map: file-2 ID → renamed-to ID. We apply this
	// before forming the intersection so a mapped pair counts as common.
	if params.DiffIndvMap != "" {
		m, err := loadDiffIndvMap(params.DiffIndvMap)
		if err != nil {
			return nil, err
		}
		r.indvMap = m
	}

	// Build the set of effective file-2 sample names (renamed where the map
	// applies). For each effective name we remember the *raw* file-2 ID so
	// we can look up genotypes in the file-2 record.
	f2Effective := make(map[string]string, len(d.samples))
	for _, s2 := range d.samples {
		eff := s2
		if r.indvMap != nil {
			if renamed, ok := r.indvMap[s2]; ok {
				eff = renamed
				r.file2RenamedSet[s2] = struct{}{}
			}
		}
		// If two file-2 IDs collapse to the same effective name the first
		// wins (mirrors upstream's overwrite-into-map semantics: the second
		// occurrence updates combined_individuals[eff].second).
		if _, dup := f2Effective[eff]; !dup {
			f2Effective[eff] = s2
		}
	}

	// Common pairs: file-1 samples that match an effective file-2 name.
	for _, s1 := range samples {
		if raw2, ok := f2Effective[s1]; ok {
			r.commonPairs = append(r.commonPairs, commonPair{f1Name: s1, f2Name: raw2})
			r.indvCommon[s1] = 0
			r.indvDiscord[s1] = 0
		}
	}

	if params.DiffSite {
		path := params.OutPrefix + ".diff.sites_in_files"
		w, err := newDiffOutFile(path, "CHROM\tPOS1\tPOS2\tIN_FILE\tREF1\tREF2\tALT1\tALT2\n")
		if err != nil {
			return nil, err
		}
		r.wSitesInFiles = w
	}
	if params.DiffSiteDiscordance {
		path := params.OutPrefix + ".diff.sites"
		w, err := newDiffOutFile(path, "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n")
		if err != nil {
			return nil, err
		}
		r.wSites = w
	}
	if params.DiffSwitchError {
		sw, err := newSwitchRunner(params.OutPrefix, r.commonPairs)
		if err != nil {
			return nil, err
		}
		r.switchRun = sw
	}
	return r, nil
}

// loadDiffIndvMap parses a two-column whitespace-separated mapping file. Each
// line is "<file-2 ID> <file-1 ID>"; blank lines and lines starting with '#'
// are skipped. Mirrors upstream variant_file_diff.cpp:11-34.
func loadDiffIndvMap(path string) (map[string]string, error) {
	f, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("opening --diff-indv-map %s: %w", path, err)
	}
	defer f.Close()

	m := make(map[string]string)
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		// Skip blank lines and comments. Upstream tests the first byte for
		// '#'; we match that.
		trimmed := strings.TrimLeft(line, " \t")
		if trimmed == "" || strings.HasPrefix(trimmed, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			// Upstream's `map >> indv1 >> indv2` silently leaves indv2 as ""
			// if the second token is missing. We follow suit by skipping —
			// keeps fixtures forgiving of trailing blanks.
			continue
		}
		m[fields[0]] = fields[1]
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("reading --diff-indv-map %s: %w", path, err)
	}
	return m, nil
}

// addVariant is called once per file-1 variant after filtering.
func (r *diffRunner) addVariant(v *vcf.Variant) error {
	if r == nil {
		return nil
	}
	if _, ok := r.seenSites[v.Chrom]; !ok {
		r.seenSites[v.Chrom] = make(map[int]struct{})
	}
	r.seenSites[v.Chrom][v.Pos] = struct{}{}

	var rec *diffRecord
	if bucket, ok := r.data.sites[v.Chrom]; ok {
		rec = bucket[v.Pos]
	}

	if r.wSitesInFiles != nil {
		if err := r.emitSitesInFilesRow(v, rec); err != nil {
			return err
		}
	}
	// --diff-switch-error: only contributes when the site is present in
	// both files (rec != nil). The runner itself enforces het-and-phased.
	r.switchRun.addVariant(v, rec)
	if rec == nil {
		// File-1-only site. .diff.sites still gets a zero row with
		// FILES=1, MATCHING_ALLELES=0, all counts 0, discordance -nan
		// (variant_file_diff.cpp:801-803 emits the row prefix; the
		// "0" MATCHING_ALLELES literal is line 929 + counts at line 933).
		if r.wSites != nil {
			fmt.Fprintf(r.wSites.w, "%s\t%d\t1\t0\t0\t0\t-nan\n", v.Chrom, v.Pos)
		}
		return nil
	}

	// Both files have the site → compute per-site and per-individual
	// discordance, but only over samples present in both files.
	siteCommon, siteDiscord := 0, 0
	// --diff-discordance-matrix counts biallelic loci with matching ALT
	// alleles only (upstream variant_file_diff.cpp:1122-1129). We pre-
	// compute the eligibility flag once per site.
	matrixEligible := r.params.DiffDiscordanceMatrix &&
		isBiallelic(v.Ref, v.Alt) && isBiallelic(rec.ref, rec.rawALTs) &&
		altsMatchFirst(v.Alt, rec.rawALTs)
	for _, pair := range r.commonPairs {
		gt1, ok1 := findSampleGT(v, pair.f1Name)
		gt2, ok2 := rec.genotypes[pair.f2Name]
		if !ok1 || !ok2 {
			continue
		}
		a1, b1, miss1 := canonicalBiallelicGT(gt1)
		a2, b2, miss2 := canonicalBiallelicGT(gt2)

		if matrixEligible {
			r.discMatrix[gtCode(a1, b1, miss1)][gtCode(a2, b2, miss2)]++
		}

		if miss1 || miss2 {
			continue
		}
		siteCommon++
		r.indvCommon[pair.f1Name]++
		if a1 != a2 || b1 != b2 {
			siteDiscord++
			r.indvDiscord[pair.f1Name]++
		}
	}
	if r.wSites != nil {
		// MATCHING_ALLELES is upstream's `alleles_match = (ALT1 == ALT2)
		// && (REF1 == REF2)` (variant_file_diff.cpp:844). For multi-ALT
		// sites we compare the joined ALT strings, which matches the
		// upstream `ALT1 == ALT2` string comparison over the whole list.
		matching := 0
		if v.Ref == rec.ref && strings.Join(v.Alt, ",") == strings.Join(rec.rawALTs, ",") {
			matching = 1
		}
		fmt.Fprintf(r.wSites.w, "%s\t%d\tB\t%d\t%d\t%d\t%s\n",
			v.Chrom, v.Pos, matching, siteCommon, siteDiscord,
			formatDiscordance(siteDiscord, siteCommon))
	}
	return nil
}

// formatDiscordance renders N_DISCORD / N_COMMON_CALLED using the same
// six-significant-digit default precision as upstream's C++ ostream
// (variant_file_diff.cpp:625 + :932). Division by zero produces "-nan"
// to match libstdc++'s `<< nan` output verbatim. Whole values print
// without a decimal point ("1", "0", "0.5") — see `%g` semantics in Go,
// which match the C `%g` format upstream relies on.
func formatDiscordance(numerator, denominator int) string {
	if denominator == 0 {
		return "-nan"
	}
	return fmt.Sprintf("%.6g", float64(numerator)/float64(denominator))
}

// gtCode encodes a canonicalised biallelic diploid genotype as one of:
//
//	0 = 0/0, 1 = 0/1, 2 = 1/1, 3 = ./. (missing)
//
// This matches the indexing the upstream discordance-matrix output uses
// (variant_file_diff.cpp:1160-1175 where N = a+b for non-missing genotypes
// and N=3 for the missing case).
func gtCode(a, b int, missing bool) int {
	if missing {
		return 3
	}
	return a + b
}

// isBiallelic reports whether a site has exactly one REF allele plus one ALT
// allele. The discordance matrix only counts biallelic sites (upstream
// variant_file_diff.cpp:1122).
func isBiallelic(ref string, alt []string) bool {
	if ref == "" {
		return false
	}
	return len(alt) == 1 && alt[0] != ""
}

// altsMatchFirst returns true when both ALT-allele slices start with the
// same first allele. Upstream's discordance-matrix path skips sites with
// mismatching ALT to avoid having to compare full genotype strings
// (variant_file_diff.cpp:1125-1129).
func altsMatchFirst(a, b []string) bool {
	if len(a) == 0 || len(b) == 0 {
		return false
	}
	return a[0] == b[0]
}

// emitSitesInFilesRow writes one row of .diff.sites_in_files for a file-1
// variant. rec is the matching file-2 record (or nil if file-1-only).
func (r *diffRunner) emitSitesInFilesRow(v *vcf.Variant, rec *diffRecord) error {
	pos1 := fmt.Sprintf("%d", v.Pos)
	ref1 := v.Ref
	alt1 := "."
	if len(v.Alt) > 0 {
		alt1 = strings.Join(v.Alt, ",")
	}
	if rec == nil {
		_, err := fmt.Fprintf(r.wSitesInFiles.w, "%s\t%s\t.\t1\t%s\t.\t%s\t.\n",
			v.Chrom, pos1, ref1, alt1)
		return err
	}
	alt2 := "."
	if len(rec.rawALTs) > 0 {
		alt2 = strings.Join(rec.rawALTs, ",")
	}
	_, err := fmt.Fprintf(r.wSitesInFiles.w, "%s\t%s\t%s\tB\t%s\t%s\t%s\t%s\n",
		v.Chrom, pos1, pos1, ref1, rec.ref, alt1, alt2)
	return err
}

// close flushes per-site outputs, walks file-2 to emit file-2-only rows for
// .diff.sites_in_files, and writes the per-individual outputs.
func (r *diffRunner) close() error {
	if r == nil {
		return nil
	}
	// File-2-only sites for sites_in_files and .diff.sites. Both outputs
	// walk the same residual set so we compute it once.
	if r.wSitesInFiles != nil || r.wSites != nil {
		chroms := make([]string, 0, len(r.data.sites))
		for c := range r.data.sites {
			chroms = append(chroms, c)
		}
		sort.Strings(chroms)
		for _, c := range chroms {
			bucket := r.data.sites[c]
			seen := r.seenSites[c]
			positions := make([]int, 0, len(bucket))
			for p := range bucket {
				if _, ok := seen[p]; ok {
					continue
				}
				positions = append(positions, p)
			}
			sort.Ints(positions)
			for _, p := range positions {
				rec := bucket[p]
				if r.wSitesInFiles != nil {
					alt2 := "."
					if len(rec.rawALTs) > 0 {
						alt2 = strings.Join(rec.rawALTs, ",")
					}
					if _, err := fmt.Fprintf(r.wSitesInFiles.w, "%s\t.\t%d\t2\t.\t%s\t.\t%s\n",
						c, p, rec.ref, alt2); err != nil {
						return err
					}
				}
				if r.wSites != nil {
					// File-2-only zero row. Upstream prints "0" for
					// MATCHING_ALLELES on non-B rows and -nan for the
					// 0/0 discordance (variant_file_diff.cpp:929-933).
					if _, err := fmt.Fprintf(r.wSites.w, "%s\t%d\t2\t0\t0\t0\t-nan\n",
						c, p); err != nil {
						return err
					}
				}
			}
		}
	}

	if err := r.wSitesInFiles.close(); err != nil {
		return err
	}
	if err := r.wSites.close(); err != nil {
		return err
	}

	if r.params.DiffIndv {
		if err := r.writeIndvInFiles(); err != nil {
			return err
		}
	}
	if r.params.DiffIndvDiscordance {
		if err := r.writeIndvDiscordance(); err != nil {
			return err
		}
	}
	if r.params.DiffDiscordanceMatrix {
		if err := r.writeDiscordanceMatrix(); err != nil {
			return err
		}
	}
	if err := r.switchRun.close(); err != nil {
		return err
	}
	return nil
}

// writeIndvInFiles writes <prefix>.diff.indv_in_files. With --diff-indv-map
// in effect, file-2 IDs that the map renames appear under their renamed (file-
// 1) ID; unmapped file-2 IDs keep their raw name. This mirrors upstream
// combined_individuals keying (variant_file_diff.cpp:36-57).
func (r *diffRunner) writeIndvInFiles() error {
	path := r.params.OutPrefix + ".diff.indv_in_files"
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return fmt.Errorf("opening %s: %w", path, err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	if _, err := fmt.Fprintln(w, "INDV\tFILES"); err != nil {
		return err
	}

	// f2Effective: set of effective file-2 IDs (post-map). Used for the "in
	// file 2" test below.
	f2Effective := make(map[string]struct{}, len(r.data.samples))
	for _, s := range r.data.samples {
		eff := s
		if r.indvMap != nil {
			if renamed, ok := r.indvMap[s]; ok {
				eff = renamed
			}
		}
		f2Effective[eff] = struct{}{}
	}

	// All names, sorted for stable output. We emit the effective name (the
	// renamed one when the map applies), never the raw file-2 ID for
	// mapped samples.
	all := make(map[string]struct{}, len(r.samples)+len(r.data.samples))
	for _, s := range r.samples {
		all[s] = struct{}{}
	}
	for s := range f2Effective {
		all[s] = struct{}{}
	}
	names := make([]string, 0, len(all))
	for s := range all {
		names = append(names, s)
	}
	sort.Strings(names)

	for _, s := range names {
		_, in1 := r.file1SampleSet[s]
		_, in2 := f2Effective[s]
		var tag string
		switch {
		case in1 && in2:
			tag = "B"
		case in1:
			tag = "1"
		default:
			tag = "2"
		}
		if _, err := fmt.Fprintf(w, "%s\t%s\n", s, tag); err != nil {
			return err
		}
	}
	return nil
}

// writeDiscordanceMatrix writes <prefix>.diff.discordance_matrix in the
// upstream 5x5 layout: a single header row (file-1 genotype labels) followed
// by four data rows (file-2 genotype labels). Cells are
// discMatrix[file1GT][file2GT]; the printed order is file2-as-row,
// file1-as-column (variant_file_diff.cpp:1194-1198).
func (r *diffRunner) writeDiscordanceMatrix() error {
	path := r.params.OutPrefix + ".diff.discordance_matrix"
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return fmt.Errorf("opening %s: %w", path, err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	if _, err := fmt.Fprintln(w, "-\tN_0/0_file1\tN_0/1_file1\tN_1/1_file1\tN_./._file1"); err != nil {
		return err
	}
	rowLabels := [4]string{"N_0/0_file2", "N_0/1_file2", "N_1/1_file2", "N_./._file2"}
	for row := 0; row < 4; row++ {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%d\n",
			rowLabels[row],
			r.discMatrix[0][row], r.discMatrix[1][row],
			r.discMatrix[2][row], r.discMatrix[3][row]); err != nil {
			return err
		}
	}
	return nil
}

// writeIndvDiscordance writes <prefix>.diff.indv. Mirrors upstream's
// output_discordance_by_indv (variant_file_diff.cpp:338-633), which iterates
// `combined_individuals` — the std::map-sorted UNION of file-1 samples and
// (post-map) file-2 samples. Samples that appear in only one file get
// N_COMMON_CALLED=0, N_DISCORD=0, DISCORDANCE=-nan.
func (r *diffRunner) writeIndvDiscordance() error {
	path := r.params.OutPrefix + ".diff.indv"
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return fmt.Errorf("opening %s: %w", path, err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	if _, err := fmt.Fprintln(w, "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE"); err != nil {
		return err
	}
	// Union of all sample names: file-1 names + effective file-2 names
	// (post --diff-indv-map renaming). Upstream's combined_individuals is
	// a std::map<string, ...> so we sort alphabetically to match its
	// iteration order.
	all := make(map[string]struct{}, len(r.samples)+len(r.data.samples))
	for _, s := range r.samples {
		all[s] = struct{}{}
	}
	for _, s2 := range r.data.samples {
		eff := s2
		if r.indvMap != nil {
			if renamed, ok := r.indvMap[s2]; ok {
				eff = renamed
			}
		}
		all[eff] = struct{}{}
	}
	names := make([]string, 0, len(all))
	for s := range all {
		names = append(names, s)
	}
	sort.Strings(names)

	for _, name := range names {
		n := r.indvCommon[name]
		d := r.indvDiscord[name]
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s\n",
			name, n, d, formatDiscordance(d, n)); err != nil {
			return err
		}
	}
	return nil
}

// findSampleGT returns the GT string for the named sample (and whether it
// was found at all).
func findSampleGT(v *vcf.Variant, name string) (string, bool) {
	for i := range v.Samples {
		if v.Samples[i].Name == name {
			gt, ok := v.Samples[i].Data["GT"]
			return gt, ok
		}
	}
	return "", false
}

// canonicalBiallelicGT parses a GT string into a sorted (smaller-first) pair
// of allele indices restricted to {0,1}; anything else (".", "2/2", "./.",
// haploid calls, etc.) is reported as missing.
func canonicalBiallelicGT(gt string) (a, b int, missing bool) {
	if gt == "" || gt == "." || gt == "./." || gt == ".|." {
		return 0, 0, true
	}
	parts := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
	if len(parts) != 2 {
		return 0, 0, true
	}
	parsed := make([]int, 0, 2)
	for _, p := range parts {
		if p == "." {
			return 0, 0, true
		}
		// We restrict to "0" or "1"; any larger index is treated as missing.
		switch p {
		case "0":
			parsed = append(parsed, 0)
		case "1":
			parsed = append(parsed, 1)
		default:
			return 0, 0, true
		}
	}
	a, b = parsed[0], parsed[1]
	if a > b {
		a, b = b, a
	}
	return a, b, false
}

// envDiffWarn writes a one-time stderr warning that the diff file lacks
// samples present in file-1 — kept here to centralise the message.
func envDiffWarn(missing []string) {
	if len(missing) == 0 {
		return
	}
	fmt.Fprintf(os.Stderr, "warning: %d sample(s) from --diff file are absent from file-1\n", len(missing))
}
