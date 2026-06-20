// This file adds small helpers shared by the conformance and edge-case test
// batteries: locating fixture corpora under reference_code/ and normalising
// SAM text so that two functionally identical alignment files compare equal
// despite cosmetic differences (notably @PG provenance lines, whose CL: and
// VN: fields legitimately differ between our binaries and upstream).

package upstream

import (
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// HtslibTestDir returns the absolute path to reference_code/htslib/test and
// whether it exists (i.e. the htslib submodule is initialised). Callers should
// SKIP gracefully when ok is false.
func HtslibTestDir() (dir string, ok bool) {
	root, err := RepoRoot()
	if err != nil {
		return "", false
	}
	dir = filepath.Join(root, "reference_code", "htslib", "test")
	if st, err := os.Stat(dir); err == nil && st.IsDir() {
		return dir, true
	}
	return dir, false
}

// HtscodecsTestDir returns the absolute path to
// reference_code/htscodecs/tests and whether it exists (the htscodecs
// submodule corpus). Callers should SKIP gracefully when ok is false.
func HtscodecsTestDir() (dir string, ok bool) {
	root, err := RepoRoot()
	if err != nil {
		return "", false
	}
	dir = filepath.Join(root, "reference_code", "htscodecs", "tests")
	if st, err := os.Stat(dir); err == nil && st.IsDir() {
		return dir, true
	}
	return dir, false
}

// NormalizeSAM canonicalises SAM text for comparison between two producers.
// It drops @PG header lines entirely (their CL:/VN:/PP: fields differ between
// our binaries and upstream and carry no record-level information), strips the
// volatile reference-provenance subfields M5: and UR: from @SQ lines (UR: is a
// machine-specific absolute path that can never match across checkouts; M5: is
// the reference MD5 that CRAM auto-annotates), and trims a trailing newline.
// Record lines and all other header lines/fields are preserved verbatim. The
// optional sortRecords flag additionally sorts the alignment records, which is
// useful when comparing producers that may legitimately differ in record order
// (e.g. CRAM slice repacking) but must preserve the set of records and their
// fields.
func NormalizeSAM(sam string, sortRecords bool) string {
	lines := strings.Split(sam, "\n")
	var hd string // the single @HD line, if any (must stay first per spec)
	var header []string
	var records []string
	for _, ln := range lines {
		if ln == "" {
			continue
		}
		if strings.HasPrefix(ln, "@") {
			if strings.HasPrefix(ln, "@PG") {
				continue
			}
			if strings.HasPrefix(ln, "@HD") {
				hd = ln
				continue
			}
			if strings.HasPrefix(ln, "@SQ") {
				ln = stripSQProvenance(ln)
			}
			header = append(header, ln)
			continue
		}
		records = append(records, ln)
	}
	// Sort the non-@HD header lines so that producers which differ only in the
	// relative ordering of @SQ vs @CO lines (a cosmetic, spec-allowed
	// difference — only @HD is positionally significant) still compare equal.
	// Header *content* parity is still enforced; only line order is relaxed.
	sort.Strings(header)
	if sortRecords {
		sort.Strings(records)
	}
	var out []string
	if hd != "" {
		out = append(out, hd)
	}
	out = append(out, header...)
	out = append(out, records...)
	return strings.Join(out, "\n")
}

// stripSQProvenance removes the M5: and UR: tab-delimited subfields from an
// @SQ header line, leaving SN:/LN: and any other tags intact.
func stripSQProvenance(line string) string {
	parts := strings.Split(line, "\t")
	kept := parts[:0]
	for _, p := range parts {
		if strings.HasPrefix(p, "M5:") || strings.HasPrefix(p, "UR:") {
			continue
		}
		kept = append(kept, p)
	}
	return strings.Join(kept, "\t")
}
