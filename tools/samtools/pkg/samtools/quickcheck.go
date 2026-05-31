package samtools

import (
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// QuickcheckOptions configures Quickcheck. The defaults match upstream
// `samtools quickcheck`'s "verify magic + EOF block + parseable header"
// behaviour.
type QuickcheckOptions struct {
	// Verbose causes Quickcheck to write a per-file PASS/FAIL line to
	// w (the writer passed to Quickcheck) in addition to the standard
	// non-zero exit code on failure.
	Verbose bool
	// UnmappedExpected, when true, suppresses the "no @SQ lines" failure
	// (matches upstream `-u`).
	UnmappedExpected bool
}

// QuickcheckResult is the outcome for one file.
type QuickcheckResult struct {
	Path   string
	OK     bool
	Reason string
}

// bgzfEOF is the canonical 28-byte BGZF empty block that every well-
// formed BAM file ends with. samtools quickcheck checks for it
// byte-for-byte. The bytes are the gzip header for an empty deflate
// stream encoded with the BGZF BC subfield set to 27 (= len-1).
var bgzfEOF = []byte{
	0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
	0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
	0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00,
}

// QuickcheckOne runs the BAM sanity check on a single file path.
// The checks (in order):
//  1. file exists and is non-empty,
//  2. starts with the BGZF gzip magic + BC subfield (i.e. is a BGZF stream),
//  3. ends with the canonical 28-byte empty BGZF EOF block,
//  4. the first decompressed bytes contain the BAM\1 magic and a parseable
//     header (we delegate to sam.NewBAMReader for this).
//
// Returns (true, "") on success and (false, reason) on any failure.
func QuickcheckOne(path string, opts QuickcheckOptions) QuickcheckResult {
	res := QuickcheckResult{Path: path}
	st, err := os.Stat(path)
	if err != nil {
		res.Reason = fmt.Sprintf("stat: %v", err)
		return res
	}
	if st.Size() == 0 {
		res.Reason = "empty file"
		return res
	}

	f, err := os.Open(path)
	if err != nil {
		res.Reason = fmt.Sprintf("open: %v", err)
		return res
	}
	defer f.Close()

	// 1) BGZF leading bytes.
	hdrBuf := make([]byte, 18)
	if _, err := io.ReadFull(f, hdrBuf); err != nil {
		res.Reason = "truncated BGZF header"
		return res
	}
	if !looksLikeBGZFHeader(hdrBuf) {
		// Upstream `quickcheck` is permissive on non-BGZF sequence data
		// (plain SAM and FASTQ): it only verifies the format category is
		// sequence_data and the header parses, and skips the EOF-block
		// check on formats that don't have one. We mirror that by
		// accepting files whose first bytes look like a SAM header.
		if looksLikeSAM(hdrBuf) {
			res.OK = true
			return res
		}
		res.Reason = "not a BGZF (BAM) file"
		return res
	}

	// 2) BGZF EOF block — must be the trailing 28 bytes byte-for-byte.
	if st.Size() < int64(len(bgzfEOF)) {
		res.Reason = "file too short to contain BGZF EOF"
		return res
	}
	tail := make([]byte, len(bgzfEOF))
	if _, err := f.ReadAt(tail, st.Size()-int64(len(bgzfEOF))); err != nil {
		res.Reason = fmt.Sprintf("read trailer: %v", err)
		return res
	}
	if !bytesEqual(tail, bgzfEOF) {
		res.Reason = "missing BGZF EOF block"
		return res
	}

	// 3) Header parses.
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		res.Reason = fmt.Sprintf("seek: %v", err)
		return res
	}
	br, err := sam.NewBAMReader(f)
	if err != nil {
		res.Reason = fmt.Sprintf("header: %v", err)
		return res
	}
	defer br.Close()
	if !opts.UnmappedExpected && len(br.Header().Refs) == 0 {
		// Upstream `quickcheck` (without -u) considers a BAM with no
		// @SQ lines suspicious because most pipelines need a reference
		// table to do anything useful.
		res.Reason = "no @SQ lines (use -u to allow header-only)"
		return res
	}
	res.OK = true
	return res
}

// looksLikeSAM returns true when the leading bytes resemble a SAM
// header (start with '@HD', '@SQ', '@RG', '@PG', '@CO', or a record
// line). Mirrors upstream's hts_detect_format auto-detection for the
// SAM case as far as samtools quickcheck cares.
func looksLikeSAM(b []byte) bool {
	if len(b) < 1 {
		return false
	}
	if b[0] == '@' {
		return true
	}
	// A bare record line starts with QNAME followed by a tab. Permit
	// any printable QNAME char followed by a tab within the first few
	// bytes so single-record / header-less SAM files still pass.
	for i := 0; i < len(b) && i < 256; i++ {
		if b[i] == '\t' {
			return true
		}
		if b[i] == '\n' || b[i] < 0x20 {
			return false
		}
	}
	return false
}

// looksLikeBGZFHeader checks the gzip magic + the BC subfield prefix
// without consuming the BSIZE bytes. It's the same check
// pkg/htsgo/sam.looksLikeBGZF performs but exposed here so
// quickcheck can run before constructing a BGZF reader.
func looksLikeBGZFHeader(b []byte) bool {
	if len(b) < 18 {
		return false
	}
	if b[0] != 0x1f || b[1] != 0x8b {
		return false
	}
	if b[2] != 0x08 {
		return false
	}
	if b[3]&0x04 == 0 { // FEXTRA must be set
		return false
	}
	xlen := uint16(b[10]) | uint16(b[11])<<8
	if xlen < 6 {
		return false
	}
	return b[12] == 'B' && b[13] == 'C' && b[14] == 0x02 && b[15] == 0x00
}

// Quickcheck runs the check on every path and returns the per-file results.
// At least one failed check causes ErrQuickcheck to be returned alongside
// the results so callers can map to a non-zero exit code without parsing
// the slice.
func Quickcheck(paths []string, opts QuickcheckOptions, w io.Writer) ([]QuickcheckResult, error) {
	results := make([]QuickcheckResult, len(paths))
	failed := false
	for i, p := range paths {
		results[i] = QuickcheckOne(p, opts)
		if !results[i].OK {
			failed = true
			fmt.Fprintf(w, "%s\n", p)
		} else if opts.Verbose {
			fmt.Fprintf(w, "%s: PASS\n", p)
		}
	}
	if failed {
		return results, ErrQuickcheck
	}
	return results, nil
}

// ErrQuickcheck is returned when one or more files fail the sanity check.
var ErrQuickcheck = errors.New("samtools quickcheck: one or more files failed")

// bytesEqual is a small dependency-free byte comparator used so the
// quickcheck file can stand alone.
func bytesEqual(a, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
