// Read-group → sample dispatch for `bcftools mpileup -G/--read-groups`.
//
// This is a faithful port of upstream's bam_sample.c. The -G file maps
// read-group IDs to output sample names; reads are routed to output
// sample columns BY THEIR RG TAG rather than by input file. One BAM
// whose reads carry multiple RGs can therefore populate multiple output
// columns, and conversely several RGs (even across files) can collapse
// onto a single sample column.
//
// The default single-column-per-BAM path (deriveSample) is untouched;
// this model is only built when opts.ReadGroups (or --ignore-RG) asks
// for it.
package bcftools

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// bamSmplFile mirrors upstream's per-file file_t (bam_sample.c:37-43): a
// map from read-group ID to the resolved output sample column index, plus
// a defaultIdx used when every read in the file maps to one sample.
type bamSmplFile struct {
	fname      string
	rg2idx     map[string]int // RG-ID → output sample column index
	defaultIdx int            // column index when all RGs collapse to one sample, else -1
}

// bamSmpl is the Go port of bam_sample.c's _bam_smpl_t. It accumulates
// the output sample column list as BAMs are added, resolving each file's
// @RG lines against the -G map (rgList / rgLogic) and the -s/-S sample
// map (which mpileup applies separately; only rg handling lives here).
type bamSmpl struct {
	files    []bamSmplFile
	smpl     []string       // output sample names, in column order
	name2idx map[string]int // sample name → column index

	ignoreRG bool

	// rgList is the -G map: RG-ID → sample-rename target. An empty
	// string value means "keep the @RG SM name" (upstream stores "\t").
	// rgLogic mirrors upstream: true=include (default), false=exclude
	// (the `^` prefix).
	rgList  map[string]string
	rgLogic bool
}

// rgKeep is the upstream "\t" sentinel: the RG is listed in -G but with no
// rename, so the @RG SM tag is kept. We use a distinct sentinel rather than
// the empty string so a genuinely empty rename can never be confused.
const rgKeepSentinel = "\t"

// newBamSmpl creates an empty sample model with the name→index hash
// initialised (bam_sample.c:58 bam_smpl_init).
func newBamSmpl() *bamSmpl {
	return &bamSmpl{name2idx: map[string]int{}}
}

// parseReadGroupsFile loads upstream's -G file (bam_sample.c:324
// bam_smpl_add_readgroups). A leading `^` on the whole list flips to
// exclude mode. Each line is whitespace-separated into up to three
// fields:
//
//	RG-ID                 → include with the @RG SM name (no rename)
//	RG-ID  SAMPLE         → include, renaming the sample to SAMPLE
//	RG-ID  FILE  SAMPLE   → file-qualified form (key "RG-ID\tFILE")
//
// Backslash escaping of whitespace is honoured, matching upstream's
// escaped-char handling. The argument is the option value; `-G -`
// (literal dash) is handled by the caller before this is reached.
func parseReadGroupsFile(path string) (rgList map[string]string, rgLogic bool, err error) {
	logic := true
	if strings.HasPrefix(path, "^") {
		logic = false
		path = path[1:]
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, false, fmt.Errorf("could not read read-groups file: %w", err)
	}
	defer f.Close()

	out := map[string]string{}
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := sc.Text()
		if line == "" {
			continue
		}
		fld1, fld2, fld3 := splitRGFields(line)
		if fld1 == "" {
			continue
		}
		// ID FILE SAMPLE form: the key becomes "ID\tFILE" and the value
		// is SAMPLE (bam_sample.c:375-382).
		key := fld1
		val := fld2
		if fld3 != "" {
			key = fld1 + "\t" + fld2
			val = fld3
		}
		if val == "" {
			val = rgKeepSentinel
		}
		if prev, ok := out[key]; ok {
			if prev != val {
				return nil, false, fmt.Errorf("the read group %q was assigned to two different samples: %q and %q",
					key, prev, val)
			}
			continue
		}
		out[key] = val
	}
	if err := sc.Err(); err != nil {
		return nil, false, err
	}
	return out, logic, nil
}

// splitRGFields splits one -G line into up to three whitespace-separated
// fields with backslash escaping, matching the kstring loops in
// bam_smpl_add_readgroups (bam_sample.c:342-374).
func splitRGFields(line string) (f1, f2, f3 string) {
	fields := make([]string, 0, 3)
	var cur strings.Builder
	escaped := false
	inField := false
	flush := func() {
		if inField {
			fields = append(fields, cur.String())
			cur.Reset()
			inField = false
		}
	}
	for i := 0; i < len(line) && len(fields) < 3; i++ {
		c := line[i]
		if c == '\\' && !escaped {
			escaped = true
			inField = true
			continue
		}
		if (c == ' ' || c == '\t') && !escaped {
			flush()
			continue
		}
		inField = true
		cur.WriteByte(c)
		escaped = false
	}
	flush()
	switch len(fields) {
	case 0:
		return "", "", ""
	case 1:
		return fields[0], "", ""
	case 2:
		return fields[0], fields[1], ""
	default:
		return fields[0], fields[1], fields[2]
	}
}

// addReadgroup registers (rgID → smplName) for file (bam_sample.c:90
// bsmpl_add_readgroup). A nil-equivalent smplName ("" with seen=false)
// records that the RG was seen but excluded. rgID "*" sets the file's
// default column. New sample names extend smpl[] in encounter order,
// which is what fixes the output column ordering.
func (bs *bamSmpl) addReadgroup(file *bamSmplFile, rgID, smplName string, hasSample bool) {
	ismpl := -1
	if hasSample {
		if idx, ok := bs.name2idx[smplName]; ok {
			ismpl = idx
		} else {
			ismpl = len(bs.smpl)
			bs.smpl = append(bs.smpl, smplName)
			bs.name2idx[smplName] = ismpl
		}
	}
	if rgID == "*" {
		file.defaultIdx = ismpl
		return
	}
	if file.rg2idx == nil {
		file.rg2idx = map[string]int{}
	}
	if _, ok := file.rg2idx[rgID]; ok {
		return // duplicate @RG ID
	}
	file.rg2idx[rgID] = ismpl
}

// keepReadgroup is the Go port of bsmpl_keep_readgroup (bam_sample.c:114).
// It looks up rgID in the -G map (trying the bare ID, then "ID\tFILE",
// then "*\tFILE"), applies include/exclude logic, and rewrites smplName
// with the rename target when one is present. Returns whether the RG is
// kept.
func (bs *bamSmpl) keepReadgroup(file *bamSmplFile, rgID string, smplName *string) bool {
	rgSmpl, ok := bs.rgList[rgID]
	if !ok {
		rgSmpl, ok = bs.rgList[rgID+"\t"+file.fname]
	}
	if !ok {
		rgSmpl, ok = bs.rgList["*\t"+file.fname]
	}
	if !ok && bs.rgLogic {
		return false
	}
	if ok && !bs.rgLogic {
		return false
	}
	// A rename target that is not the keep-sentinel renames the sample.
	if ok && rgSmpl != rgKeepSentinel && rgSmpl != "" && rgSmpl[0] != '\t' {
		*smplName = rgSmpl
	}
	return true
}

// addBAM registers one input file's @RG lines, building its rg2idx /
// defaultIdx. It is the Go port of bam_smpl_add_bam (bam_sample.c:151),
// restricted to the -G / --ignore-RG paths mpileup uses here (the -s/-S
// sample_list is applied by MpileupFile separately). Returns the file
// index, or -1 if the file has no usable read group (upstream drops it).
//
// rgs are the parsed @RG entries in header order; the (id, sm) pairs are
// processed in that order so the output column order matches upstream.
func (bs *bamSmpl) addBAM(fname string, rgs []sam.ReadGroup) int {
	idx := len(bs.files)
	bs.files = append(bs.files, bamSmplFile{fname: fname, defaultIdx: -1})
	file := &bs.files[idx]

	// --ignore-RG or no @RG: the whole file is one sample named after
	// the file (bam_sample.c:160-165).
	if bs.ignoreRG || len(rgs) == 0 {
		bs.addReadgroup(file, "*", fname, true)
		return idx
	}

	firstSmpl := -1
	nskipped := 0
	bamSmpls := map[string]struct{}{}
	for _, rg := range rgs {
		id := rg.ID
		sm := rgSM(rg)
		if id == "" || sm == "" {
			continue
		}
		// restrict / rename based on -G (bam_sample.c:205-216).
		acceptRG := true
		r := sm
		if bs.rgList != nil {
			acceptRG = bs.keepReadgroup(file, id, &r)
		}
		if acceptRG {
			bs.addReadgroup(file, id, r, true)
		} else {
			bs.addReadgroup(file, id, "", false) // seen but excluded
			nskipped++
		}
		if firstSmpl < 0 {
			if v, ok := bs.name2idx[r]; ok {
				firstSmpl = v
			}
		}
		bamSmpls[r] = struct{}{}
	}
	nsmpls := len(bamSmpls)

	// Decide how reads with no/unknown RG ("?") are handled
	// (bam_sample.c:232-254).
	smplName := ""
	hasSmplName := false
	acceptNullRG := true
	if bs.rgList != nil {
		var rn string
		if !bs.keepReadgroup(file, "?", &rn) {
			acceptNullRG = false
		} else if rn != "" {
			smplName, hasSmplName = rn, true
		}
	}
	if !acceptNullRG && firstSmpl == -1 {
		// No usable read group: drop the file (bam_sample.c:237-244).
		bs.files = bs.files[:idx]
		return -1
	}
	if !acceptNullRG {
		return idx
	}
	if nsmpls == 1 && nskipped == 0 {
		file.defaultIdx = firstSmpl
		return idx
	}
	if !hasSmplName {
		if firstSmpl == -1 {
			smplName = fname
		} else {
			smplName = bs.smpl[firstSmpl]
		}
	}
	bs.addReadgroup(file, "?", smplName, true)
	return idx
}

// rgSM extracts the SM field of an @RG entry, "" when absent.
func rgSM(rg sam.ReadGroup) string {
	for _, f := range rg.Extra {
		if f.Tag == "SM" {
			return f.Value
		}
	}
	return ""
}

// sampleID resolves a read to its output sample column, the Go port of
// bam_smpl_get_sample_id (bam_sample.c:263). A negative result means the
// read is dropped (its RG is excluded). fileIdx is the index returned by
// addBAM.
func (bs *bamSmpl) sampleID(fileIdx int, rec *sam.Record) int {
	file := &bs.files[fileIdx]
	if file.defaultIdx >= 0 {
		return file.defaultIdx
	}
	auxRG := "?"
	if a, ok := rec.GetAux("RG"); ok {
		if s, ok := a.String(); ok {
			auxRG = s
		}
	}
	if id, ok := file.rg2idx[auxRG]; ok {
		return id
	}
	if id, ok := file.rg2idx["?"]; ok {
		return id
	}
	return -1
}

// samples returns the resolved output sample column names in order.
func (bs *bamSmpl) samples() []string { return bs.smpl }

// dispatchByReadGroup builds the bam_sample.c model from the -G file (or
// --ignore-RG) and re-partitions every input's reads into per-output-
// column record buckets. headers/paths/perInputRecs are parallel slices
// over the opened input BAMs (one entry per file). The returned samples
// and colRecs are parallel slices over the OUTPUT sample columns: column
// i carries the reads whose RG resolved to sample i, regardless of which
// file they came from.
//
// This is the data-model lift that lets one BAM populate several columns
// and several RGs collapse onto one column (mpileup.c:1652 →
// bam_smpl_add_bam / bam_smpl_get_sample_id).
func dispatchByReadGroup(opts MpileupOptions, headers []*sam.Header, paths []string,
	perInputRecs []map[string][]*sam.Record) (samples []string, colRecs []map[string][]*sam.Record, err error) {

	bs := newBamSmpl()
	bs.ignoreRG = opts.IgnoreRG
	if opts.ReadGroups != "" {
		bs.rgList, bs.rgLogic, err = parseReadGroupsFile(opts.ReadGroups)
		if err != nil {
			return nil, nil, fmt.Errorf("bcftools mpileup: %w", err)
		}
	}

	// fileCol[f] is the addBAM index for input file f, or -1 if upstream
	// dropped the file (no usable read group).
	fileCol := make([]int, len(headers))
	for f := range headers {
		fileCol[f] = bs.addBAM(paths[f], headers[f].ReadGroups)
	}

	nCol := len(bs.smpl)
	colRecs = make([]map[string][]*sam.Record, nCol)
	for c := range colRecs {
		colRecs[c] = map[string][]*sam.Record{}
	}

	// Dispatch every kept read to its resolved sample column. Records
	// are already chrom-bucketed and position-sorted per file; appending
	// in file order then re-sorting per (column,chrom) preserves the
	// coordinate order the pileup engine expects.
	for f := range perInputRecs {
		fc := fileCol[f]
		if fc < 0 {
			continue // file dropped: no usable read group
		}
		for chrom, recs := range perInputRecs[f] {
			for _, rec := range recs {
				col := bs.sampleID(fc, rec)
				if col < 0 {
					continue // RG excluded in include mode (or unknown)
				}
				colRecs[col][chrom] = append(colRecs[col][chrom], rec)
			}
		}
	}
	for c := range colRecs {
		for chrom := range colRecs[c] {
			recs := colRecs[c][chrom]
			sort.SliceStable(recs, func(i, j int) bool { return recs[i].Pos < recs[j].Pos })
		}
	}

	return bs.samples(), colRecs, nil
}
