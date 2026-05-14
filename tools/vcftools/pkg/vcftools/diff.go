// VCF comparison ("--diff" family) for vcftools.
//
// Supported flags in this v1:
//
//	--diff FILE                 second VCF to compare against
//	--diff-site                 emit <prefix>.diff.sites_in_files
//	--diff-indv                 emit <prefix>.diff.indv_in_files
//	--diff-site-discordance     emit <prefix>.diff.sites
//	--diff-indv-discordance     emit <prefix>.diff.indv
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
//	.diff.sites                CHROM POS N_COMMON_CALLED N_DISCORD
//	                           (sites present in both files; counts only
//	                           samples called in both)
//
//	.diff.indv                 INDV N_COMMON_CALLED N_DISCORD
//	                           (samples present in both files; counts only
//	                           sites called in both that are present in
//	                           both files)
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

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
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

	d := &diffData{
		samples: append([]string(nil), hdr.Samples...),
		sites:   make(map[string]map[int]*diffRecord),
	}

	for {
		v, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading --diff VCF: %w", err)
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
		// If the same position appears twice in file-2, the last one wins —
		// matches upstream's "overwrite" behaviour.
		bucket[v.Pos] = rec
	}
	return d, nil
}

// diffRunner accumulates per-site and per-individual stats while file-1 is
// streamed through Run, then flushes outputs in close().
type diffRunner struct {
	params  *Params
	data    *diffData
	samples []string // file-1 samples (post sample filtering)

	// Set of file-1 samples — quick lookup for the indv_in_files report.
	file1SampleSet map[string]struct{}

	// Set of (chrom,pos) seen in file-1 for the sites_in_files report's
	// "B" classification. Once a site is seen on the file-1 side we mark it
	// here; at close() time we walk file-2's sites and any not in this map
	// are "2"-only.
	seenSites map[string]map[int]struct{}

	// Per-sample discordance accumulators (samples shared by both files).
	indvCommon   map[string]int
	indvDiscord  map[string]int
	commonSample []string // sorted, intersection of file1 and file2 samples

	// Output writers.
	wSitesInFiles *diffOutFile
	wSites        *diffOutFile
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
// list. It is the caller's responsibility to ensure params.Diff != "" and at
// least one --diff-* output flag is set.
func newDiffRunner(params *Params, samples []string) (*diffRunner, error) {
	if params.Diff == "" {
		return nil, nil
	}
	if !params.DiffSite && !params.DiffIndv && !params.DiffSiteDiscordance && !params.DiffIndvDiscordance {
		return nil, nil
	}
	d, err := loadDiffVCF(params.Diff)
	if err != nil {
		return nil, err
	}

	r := &diffRunner{
		params:         params,
		data:           d,
		samples:        append([]string(nil), samples...),
		file1SampleSet: make(map[string]struct{}, len(samples)),
		seenSites:      make(map[string]map[int]struct{}),
		indvCommon:     make(map[string]int),
		indvDiscord:    make(map[string]int),
	}
	for _, s := range samples {
		r.file1SampleSet[s] = struct{}{}
	}
	// Common samples: those in both file1 and file2, preserving file-1 order.
	f2Set := make(map[string]struct{}, len(d.samples))
	for _, s := range d.samples {
		f2Set[s] = struct{}{}
	}
	for _, s := range samples {
		if _, ok := f2Set[s]; ok {
			r.commonSample = append(r.commonSample, s)
			r.indvCommon[s] = 0
			r.indvDiscord[s] = 0
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
		w, err := newDiffOutFile(path, "CHROM\tPOS\tN_COMMON_CALLED\tN_DISCORD\n")
		if err != nil {
			return nil, err
		}
		r.wSites = w
	}
	return r, nil
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
	if rec == nil {
		return nil
	}

	// Both files have the site → compute per-site and per-individual
	// discordance, but only over samples present in both files.
	siteCommon, siteDiscord := 0, 0
	for _, name := range r.commonSample {
		gt1, ok1 := findSampleGT(v, name)
		gt2, ok2 := rec.genotypes[name]
		if !ok1 || !ok2 {
			continue
		}
		a1, b1, miss1 := canonicalBiallelicGT(gt1)
		a2, b2, miss2 := canonicalBiallelicGT(gt2)
		if miss1 || miss2 {
			continue
		}
		siteCommon++
		r.indvCommon[name]++
		if a1 != a2 || b1 != b2 {
			siteDiscord++
			r.indvDiscord[name]++
		}
	}
	if r.wSites != nil {
		fmt.Fprintf(r.wSites.w, "%s\t%d\t%d\t%d\n", v.Chrom, v.Pos, siteCommon, siteDiscord)
	}
	return nil
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
	// File-2-only sites for sites_in_files.
	if r.wSitesInFiles != nil {
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
				alt2 := "."
				if len(rec.rawALTs) > 0 {
					alt2 = strings.Join(rec.rawALTs, ",")
				}
				if _, err := fmt.Fprintf(r.wSitesInFiles.w, "%s\t.\t%d\t2\t.\t%s\t.\t%s\n",
					c, p, rec.ref, alt2); err != nil {
					return err
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
	return nil
}

// writeIndvInFiles writes <prefix>.diff.indv_in_files.
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

	f2Set := make(map[string]struct{}, len(r.data.samples))
	for _, s := range r.data.samples {
		f2Set[s] = struct{}{}
	}

	// All names, sorted for stable output.
	all := make(map[string]struct{}, len(r.samples)+len(r.data.samples))
	for _, s := range r.samples {
		all[s] = struct{}{}
	}
	for _, s := range r.data.samples {
		all[s] = struct{}{}
	}
	names := make([]string, 0, len(all))
	for s := range all {
		names = append(names, s)
	}
	sort.Strings(names)

	for _, s := range names {
		_, in1 := r.file1SampleSet[s]
		_, in2 := f2Set[s]
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

// writeIndvDiscordance writes <prefix>.diff.indv.
func (r *diffRunner) writeIndvDiscordance() error {
	path := r.params.OutPrefix + ".diff.indv"
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return fmt.Errorf("opening %s: %w", path, err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	if _, err := fmt.Fprintln(w, "INDV\tN_COMMON_CALLED\tN_DISCORD"); err != nil {
		return err
	}
	// Emit in file-1 sample order for the intersection.
	for _, s := range r.commonSample {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\n", s, r.indvCommon[s], r.indvDiscord[s]); err != nil {
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
