// Mendelian inconsistency detection for vcftools.
//
// CLI flag: --mendel <PED>.
//
// Upstream source: reference_code/vcftools/src/cpp/variant_file_output.cpp:
// 5332-5470 (variant_file::output_mendel_inconsistencies). Registered in
// parameters.cpp:293:
//
//	else if (in_str == "--mendel") {
//	    mendel_ped_file = get_arg(i+1); i++; num_outputs++;
//	}
//
// PED format (whitespace-separated, four columns parsed):
//
//	family  child  father  mother
//
// The first line of the file is *always* skipped (upstream calls
// PED.ignore(..., '\n') unconditionally). Lines starting with '#' or empty
// lines are also skipped. Rows where any of child/father/mother is the
// literal "0" are dropped (upstream lines 5356-5357), as is any trio whose
// child/mother/father isn't in the VCF sample list (lines 5370-5381).
//
// Output: <prefix>.mendel — header
//
//	CHR\tPOS\tREF\tALT\tFAMILY\tCHILD\tFATHER\tMOTHER
//
// Where FAMILY (the upstream `family_ids` entry) is the string
// `<child>_<father>_<mother>` (upstream line 5380); CHILD / FATHER / MOTHER
// are the three "a/b" allele-index pairs ('a' separated from 'b' by '/')
// at the offending site.
//
// A site is reported as an inconsistency when **neither** the child's
// (allele1, allele2) pair nor its swapped (allele2, allele1) pair is in
// the set of four possible Mendelian children genotypes
// (mother×father pairings) — see upstream lines 5450-5458.
//
// Missing alleles (".") in any of the three trio members → site skipped
// for that trio. The implementation here operates per-trio per-site;
// independent trios are independent.

package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// mendelTrio captures one resolved (child, father, mother) tuple, plus the
// per-trio family-id string used in the output column.
type mendelTrio struct {
	family string // family_ids entry: "<child>_<father>_<mother>"
	child  string
	father string
	mother string
	// Sample-index pointers into the filtered VCF sample list. Stored at
	// initialisation time so we don't have to scan v.Samples per row.
	childIdx  int
	fatherIdx int
	motherIdx int
}

// mendelRunner buffers an output writer and the trio list; trios are matched
// to the filtered sample list once at construction time so the per-row hot
// path is two map lookups.
type mendelRunner struct {
	prefix string
	trios  []mendelTrio

	out *mendelOutFile
}

type mendelOutFile struct {
	f io.WriteCloser
	w *bufio.Writer
}

// newMendelRunner reads the PED file, intersects it with `samples`, opens the
// output file, and returns a ready-to-use runner. Returns (nil, nil) when no
// trios overlap with the VCF — upstream errors out with "No PED individuals
// found in VCF." but we surface the same error to the caller for clarity.
func newMendelRunner(prefix, pedPath string, samples []string) (*mendelRunner, error) {
	trios, err := loadMendelPED(pedPath, samples)
	if err != nil {
		return nil, err
	}
	if len(trios) == 0 {
		return nil, fmt.Errorf("--mendel: no PED individuals found in VCF")
	}
	// Upstream prints "Found N trios in the VCF file.\n" to stderr; mirror
	// that for log-equivalent behaviour. Tests that rely on byte-for-byte
	// stderr would have to mask this line just like they do upstream's
	// other LOG.printLOG calls.
	fmt.Fprintf(os.Stderr, "Found %d trios in the VCF file.\n", len(trios))

	path := prefix + ".mendel"
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", path, err)
	}
	w := bufio.NewWriter(f)
	if _, err := w.WriteString("CHR\tPOS\tREF\tALT\tFAMILY\tCHILD\tFATHER\tMOTHER\n"); err != nil {
		f.Close()
		return nil, fmt.Errorf("writing header to %s: %w", path, err)
	}
	return &mendelRunner{
		prefix: prefix,
		trios:  trios,
		out:    &mendelOutFile{f: f, w: w},
	}, nil
}

// addVariant emits one row per trio whose genotype at v is Mendel-inconsistent.
// v must be the post-filter variant (i.e. sample slice matches the filtered
// header used to construct the runner).
func (r *mendelRunner) addVariant(v *vcf.Variant) error {
	if r == nil {
		return nil
	}
	for _, trio := range r.trios {
		if trio.childIdx < 0 || trio.fatherIdx < 0 || trio.motherIdx < 0 {
			continue
		}
		if trio.childIdx >= len(v.Samples) ||
			trio.fatherIdx >= len(v.Samples) ||
			trio.motherIdx >= len(v.Samples) {
			continue
		}
		childGT := v.Samples[trio.childIdx].Data["GT"]
		fatherGT := v.Samples[trio.fatherIdx].Data["GT"]
		motherGT := v.Samples[trio.motherIdx].Data["GT"]
		c1, c2, cMiss := parseDiploidGT(childGT)
		f1, f2, fMiss := parseDiploidGT(fatherGT)
		m1, m2, mMiss := parseDiploidGT(motherGT)
		if cMiss || fMiss || mMiss {
			continue
		}

		// Build the four possible child genotypes from mother × father.
		// Upstream uses a `set<pair<int,int>>` which collapses
		// duplicates; we check membership via direct comparison.
		if isMendelConsistent(c1, c2, f1, f2, m1, m2) {
			continue
		}
		// Mendel error: emit one row. CHROM / POS / REF / ALT come
		// from v; FAMILY/CHILD/FATHER/MOTHER follow upstream's
		// formatting (no zero-padding, alleles slash-separated).
		altOut := strings.Join(v.Alt, ",")
		if altOut == "" {
			altOut = "."
		}
		if _, err := fmt.Fprintf(r.out.w,
			"%s\t%d\t%s\t%s\t%s\t%d/%d\t%d/%d\t%d/%d\n",
			v.Chrom, v.Pos, v.Ref, altOut,
			trio.family,
			c1, c2, f1, f2, m1, m2); err != nil {
			return fmt.Errorf("writing %s.mendel: %w", r.prefix, err)
		}
	}
	return nil
}

// close flushes the .mendel output.
func (r *mendelRunner) close() error {
	if r == nil {
		return nil
	}
	var firstErr error
	if err := r.out.w.Flush(); err != nil {
		firstErr = err
	}
	if err := r.out.f.Close(); err != nil && firstErr == nil {
		firstErr = err
	}
	return firstErr
}

// loadMendelPED reads the four-column PED file and intersects it with the
// filtered sample list, returning the resolved per-trio sample indices.
// Rows that don't fully resolve in `samples` are dropped (upstream behaviour).
// The first line is always treated as a header (skipped). Comments are lines
// whose first byte is '#'.
func loadMendelPED(path string, samples []string) ([]mendelTrio, error) {
	f, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("opening --mendel PED file %s: %w", path, err)
	}
	defer f.Close()

	// Build a name → index map for the VCF samples.
	idx := make(map[string]int, len(samples))
	for i, s := range samples {
		idx[s] = i
	}

	var trios []mendelTrio
	scanner := bufio.NewScanner(f)
	// Upstream unconditionally skips the first line.
	first := true
	for scanner.Scan() {
		if first {
			first = false
			continue
		}
		line := scanner.Text()
		if line == "" {
			continue
		}
		if line[0] == '#' {
			continue
		}
		// PED is whitespace-separated; upstream uses `ss >> family >>
		// child >> father >> mother;` so we mirror that and accept any
		// run of whitespace (tab or space) as a separator.
		fields := strings.Fields(line)
		if len(fields) < 4 {
			continue
		}
		family, child, father, mother := fields[0], fields[1], fields[2], fields[3]
		if child == "0" || father == "0" || mother == "0" {
			continue
		}
		ci, cOK := idx[child]
		fi, fOK := idx[father]
		mi, mOK := idx[mother]
		if !cOK || !fOK || !mOK {
			continue
		}
		// family_ids upstream is "<child>_<father>_<mother>"; family
		// itself is read but not emitted. We keep the raw column to
		// note we parsed it.
		_ = family
		trios = append(trios, mendelTrio{
			family:    family,
			child:     child,
			father:    father,
			mother:    mother,
			childIdx:  ci,
			fatherIdx: fi,
			motherIdx: mi,
		})
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("reading --mendel PED %s: %w", path, err)
	}
	// Family-id string (upstream lines 5380):
	for i := range trios {
		trios[i].family = trios[i].child + "_" + trios[i].father + "_" + trios[i].mother
	}
	return trios, nil
}

// parseDiploidGT parses "a/b" or "a|b" into the two integer allele indices.
// Returns (0, 0, true) for any missing slot or malformed GT. Indices >= 0 are
// preserved verbatim so a tri-allelic site (index 2) still participates in the
// Mendel check; upstream's get_indv_GENOTYPE_ids returns the same numeric
// indices.
func parseDiploidGT(gt string) (a, b int, missing bool) {
	if gt == "" || gt == "." {
		return 0, 0, true
	}
	sep := -1
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			sep = i
			break
		}
	}
	if sep < 0 {
		// Haploid call — upstream treats unknown / haploid as missing
		// here (get_indv_GENOTYPE_ids returns -1 for the second slot).
		return 0, 0, true
	}
	left := gt[:sep]
	right := gt[sep+1:]
	if left == "." || right == "." || left == "" || right == "" {
		return 0, 0, true
	}
	la, lOK := parseLDhatAllele(left)
	rb, rOK := parseLDhatAllele(right)
	if !lOK || !rOK || la < 0 || rb < 0 {
		return 0, 0, true
	}
	return la, rb, false
}

// isMendelConsistent returns true when the child genotype (c1,c2) can be
// formed by drawing one allele from the mother and one from the father.
// Upstream builds a set of four candidate genotypes (mother × father), then
// checks both (c1,c2) and (c2,c1) against it.
func isMendelConsistent(c1, c2, f1, f2, m1, m2 int) bool {
	// Four possible orderings of (mother_allele, father_allele).
	candidates := [4][2]int{
		{m1, f1},
		{m1, f2},
		{m2, f1},
		{m2, f2},
	}
	for _, cand := range candidates {
		if cand[0] == c1 && cand[1] == c2 {
			return true
		}
		if cand[0] == c2 && cand[1] == c1 {
			return true
		}
	}
	return false
}
