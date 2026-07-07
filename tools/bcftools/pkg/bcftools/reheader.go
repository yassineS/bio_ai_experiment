// Package bcftools — `bcftools reheader` subcommand.
//
// `bcftools reheader` rewrites the header of a VCF/BCF in place. The body
// records are passed through unchanged (subject to the new header's contig
// order for BCF outputs). Three input modes are supported:
//
//   - `-h FILE` — replace the entire header with the contents of FILE.
//   - `-s FILE` — rename samples; FILE is either a one-name-per-line list
//     (positional rename) or a tab-separated `OLD\tNEW` map (by name).
//   - `-f FAI`  — rebuild `##contig` lines from a samtools FAI sidecar.
//
// Modes can be combined: the order of operations is (header replace) →
// (FAI contigs) → (sample rename), matching upstream behaviour.
//
// For BCF input, upstream htslib requires a careful header rewrite so the
// dictionary indices are preserved; this v1 port reads the BCF, applies the
// edits to the in-memory header, and re-emits the records via the BCF
// writer. The records are conceptually unchanged but are re-encoded.
package bcftools

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// ReheaderOptions controls the `bcftools reheader` subcommand.
type ReheaderOptions struct {
	// HeaderFile, when non-empty, is read verbatim and used as the new
	// header. Lines starting with `#` are kept; the file is expected to
	// terminate with the `#CHROM` header line.
	HeaderFile string
	// SamplesFile, when non-empty, drives the sample-rename step. The
	// format is auto-detected: if any non-comment line contains a tab, it
	// is treated as a `OLD\tNEW` mapping; otherwise the file is a flat
	// list of new names applied positionally.
	SamplesFile string
	// FaiFile, when non-empty, is a samtools FAI index whose first two
	// columns (NAME, LENGTH) are used to rebuild the `##contig` lines.
	FaiFile string
	// OutputFormat selects the output encoding. Defaults to OutputVCF.
	OutputFormat OutputFormat
	// OutputFormatExplicit reports whether the caller set OutputFormat
	// deliberately (our -O/--output-type extension). When false, ReheaderFile
	// mirrors the input's compression: upstream `bcftools reheader` has no -O
	// flag and emits BGZF for a BGZF (.vcf.gz) input and plain text for a plain
	// input (reheader.c main: bgzf input -> reheader_vcf_gz, plain ->
	// reheader_vcf).
	OutputFormatExplicit bool
	// CompressLevel is the gzip level for -O z output (negative means
	// gzip's default).
	CompressLevel int
	// Threads is upstream's -@/--threads. When greater than 1 it enables
	// parallel BGZF compression of -O z and -O b output via bgzf.MultiWriter;
	// the framed result decodes byte-identically regardless of thread count.
	Threads int
}

// ReheaderFile is the file-aware entry point for `bcftools reheader`. It
// opens path through iohelper, applies the requested header edits, and
// writes the records (with the new header) to out.
func ReheaderFile(path string, out io.Writer, opts ReheaderOptions) (int, error) {
	// Mirror the input's compression unless the caller forced a format via -O:
	// upstream reheader re-emits a BGZF (.vcf.gz) input as BGZF and a plain
	// input as plain text. We default OutputFormat to plain VCF, so only the
	// compressed case needs an upgrade.
	if !opts.OutputFormatExplicit && opts.OutputFormat == OutputVCF {
		gz, err := pathIsGzip(path)
		if err != nil {
			return 0, fmt.Errorf("bcftools reheader: open %s: %w", path, err)
		}
		if gz {
			opts.OutputFormat = OutputVCFGz
			// Fast path: for a BGZF (.vcf.gz) input carrying a text-VCF payload,
			// mirror upstream reheader_vcf_gz and raw-copy the compressed body
			// blocks verbatim rather than inflate+re-deflate the whole file. This
			// is eligible only when the output stays BGZF (no -O change) and the
			// single-threaded writer is in play; anything else, or a plain-gzip
			// (non-BGZF) .gz, falls through to the decode path below. A negative
			// return with an incompatibility is signalled by rawEligible=false.
			if opts.Threads <= 1 {
				n, done, err := reheaderVCFGzRaw(path, out, opts)
				if err != nil {
					return n, err
				}
				if done {
					return n, nil
				}
			}
		}
	}
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools reheader: open %s: %w", path, err)
	}
	defer in.Close()
	return Reheader(in, out, opts)
}

// reheaderVCFGzRaw implements upstream reheader_vcf_gz: for a BGZF (.vcf.gz)
// input carrying a text-VCF payload it rewrites only the header and raw-copies
// the compressed body blocks verbatim, without inflating and re-deflating the
// record body.
//
// It reports (records, done, err). done is true when the raw path handled the
// whole file; when done is false and err is nil the input was not raw-eligible
// (e.g. a plain-gzip .gz with no BGZF BC subfield, or a BCF-in-gz payload) and
// the caller must fall back to the decode/stream path. records is always 0 on
// this path (records are never materialised).
//
// The subtlety mirrored from upstream's skip_until: the header rarely ends on a
// BGZF block boundary. We therefore re-compress the block(s) that contain the
// header plus the start of the body — the new header, the body bytes we had to
// inflate while scanning for the header end, and the decompressed remainder of
// the block that straddles the header/body split are all written through a fresh
// BGZF writer and flushed. Only the SUBSEQUENT whole compressed blocks are
// raw-copied. No BGZF EOF marker is emitted until the very end.
func reheaderVCFGzRaw(path string, out io.Writer, opts ReheaderOptions) (records int, done bool, err error) {
	f, err := os.Open(path)
	if err != nil {
		return 0, false, fmt.Errorf("bcftools reheader: open %s: %w", path, err)
	}
	defer f.Close()

	br, err := bgzip.NewReader(f)
	if err != nil {
		// Not a BGZF stream (e.g. plain gzip). Not raw-eligible — fall back.
		return 0, false, nil
	}
	defer br.Close()

	// Read the leading VCF header directly from the BGZF reader, tracking how
	// many decoded body bytes we had to pull past the header end. Reading via
	// the reader's own Read keeps the underlying stream positioned at a block
	// boundary (nextBlock reads whole blocks and never reads ahead), so after
	// the header scan DecompressedRemainder + RawRemaining are valid.
	headerBytes, leadingBody, isText, err := readVCFHeaderFromBGZF(br)
	if err != nil {
		// The stream parsed a BGZF header (NewReader accepted it) but decoding
		// the first block failed — most commonly a plain-gzip .gz whose single
		// member lacks the BC subfield, which NewReader tolerates but the first
		// Read rejects. Rather than fail, fall back to the transparent-gzip
		// decode path, which handles plain gzip and reports genuine corruption.
		return 0, false, nil
	}
	if !isText {
		// Payload is not text VCF (e.g. a BCF stream framed in BGZF). The raw
		// header-text swap does not apply; fall back to the decode path.
		return 0, false, nil
	}

	newHeader, err := editVCFHeaderText(headerBytes, opts)
	if err != nil {
		return 0, false, err
	}

	bw := bgzip.NewWriter(out)
	// New (edited) header first.
	if _, err := bw.Write(newHeader); err != nil {
		bw.Close()
		return 0, false, err
	}
	// The body bytes we inflated while scanning for the header end.
	if len(leadingBody) > 0 {
		if _, err := bw.Write(leadingBody); err != nil {
			bw.Close()
			return 0, false, err
		}
	}
	// The decompressed remainder of the block that straddles the header/body
	// split — these body bytes live in the same BGZF block as (part of) the
	// header, so they must be re-compressed rather than raw-copied.
	if rem := br.DecompressedRemainder(); len(rem) > 0 {
		if _, err := bw.Write(rem); err != nil {
			bw.Close()
			return 0, false, err
		}
	}
	// Flush the header/straddle block(s) so the raw compressed blocks that
	// follow are appended directly after them in the output stream.
	if err := bw.Flush(); err != nil {
		bw.Close()
		return 0, false, err
	}
	// Raw-copy the remaining whole compressed blocks verbatim, trimming the
	// input's trailing EOF marker. When the header consumed the whole input up
	// to (and including) the input EOF block — e.g. a header-only VCF.gz — there
	// are no trailing blocks left and CopyRawTrimEOF reports ErrTruncated; that
	// is expected here (we supply our own single EOF marker below), so tolerate
	// it while still surfacing any genuine copy/IO error.
	if _, err := bgzip.CopyRawTrimEOF(out, br.RawRemaining()); err != nil && !errors.Is(err, bgzip.ErrTruncated) {
		bw.Close()
		return 0, false, fmt.Errorf("bcftools reheader: %w", err)
	}
	// Close the writer WITHOUT letting it emit its own trailing EOF after the
	// raw blocks — instead we emit exactly one EOF marker to terminate the
	// concatenated stream. bw has already been fully flushed above, so closing
	// it here only appends its EOF block; we must not double it. Emit our own.
	if _, err := out.Write(bgzip.EOFBlock); err != nil {
		return 0, false, err
	}
	return 0, true, nil
}

// readVCFHeaderFromBGZF reads the contiguous leading '#'-prefixed VCF header
// directly from a BGZF reader. It returns the header bytes, any body bytes that
// had to be inflated past the header end (the first non-'#' bytes already
// delivered by the reader), and whether the payload looks like text VCF (its
// first byte is '#'). The reader is left positioned so that DecompressedRemainder
// yields the rest of the header/body-straddling block and RawRemaining yields the
// subsequent whole compressed blocks.
func readVCFHeaderFromBGZF(br *bgzip.Reader) (header, leadingBody []byte, isText bool, err error) {
	// Read one byte at a time so we stop pulling decompressed bytes as soon as
	// the header ends — minimising the body prefix we must re-compress. Headers
	// are small (a few KB), so per-byte reads here are cheap relative to the
	// raw-copy of the (large) body that follows.
	var (
		buf     []byte
		atBOL   = true // at beginning of a line
		one     [1]byte
		hdrLen  int // number of bytes belonging to the header
		started bool
	)
	for {
		n, rerr := br.Read(one[:])
		if n > 0 {
			c := one[0]
			if !started {
				started = true
				// The first byte decides text-VCF vs. non-text payload.
				if c != '#' {
					return nil, nil, false, nil
				}
			}
			if atBOL && c != '#' {
				// First byte of a line that is not '#': the header ended at the
				// previous newline. This byte is the first body byte.
				buf = append(buf, c)
				header = buf[:hdrLen:hdrLen]
				leadingBody = buf[hdrLen:]
				return header, leadingBody, true, nil
			}
			buf = append(buf, c)
			hdrLen = len(buf)
			atBOL = c == '\n'
		}
		if rerr == io.EOF {
			// Whole input was header (no body records).
			return buf, nil, started, nil
		}
		if rerr != nil {
			return nil, nil, false, rerr
		}
	}
}

// pathIsGzip reports whether path begins with the gzip/BGZF magic bytes
// (0x1f 0x8b). A "-" (stdin) path is treated as not compressed: upstream
// detects stdin's format via hts_open, but the reheader matrix and CLI use a
// real file path, and a plain stdin stays plain. Errors other than a short
// read are returned.
func pathIsGzip(path string) (bool, error) {
	if path == "" || path == "-" {
		return false, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return false, err
	}
	defer f.Close()
	var magic [2]byte
	n, err := io.ReadFull(f, magic[:])
	if err == io.EOF || err == io.ErrUnexpectedEOF {
		return false, nil // too short to be gzip
	}
	if err != nil {
		return false, err
	}
	return n == 2 && magic[0] == 0x1f && magic[1] == 0x8b, nil
}

// Reheader edits the header of in per opts and writes the result to out.
//
// For text VCF and VCF.gz (BGZF/gzip) input the body records are never decoded:
// the leading `#`-prefixed header is split off and swapped, and the compressed
// or plain body bytes are copied through verbatim (mirroring upstream
// `bcftools reheader`, which edits the header text and passes the body through).
// BCF input keeps the decode + re-encode path because the record dictionary
// indices are tied to the header ordering.
func Reheader(in io.Reader, out io.Writer, opts ReheaderOptions) (int, error) {
	// Peek the leading bytes to distinguish BCF (needs full re-encode) from
	// text VCF (stream-copy). iohelper has already transparently decompressed
	// gzip/BGZF input by the time we get here, so a VCF.gz body arrives as the
	// underlying text.
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return 0, fmt.Errorf("bcftools reheader: %w", err)
	}
	if len(head) >= 3 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return reheaderBCF(br, out, opts)
	}
	return reheaderStream(br, out, opts)
}

// reheaderStream implements the streaming header-swap path for text VCF /
// VCF.gz input. It reads the leading `#` header without decoding any record,
// applies the requested header edits to that text, and streams the body bytes
// through verbatim (never buffering the whole body in memory). Output is wrapped
// in a BGZF writer when the selected format is VCF.gz. The returned count is
// always 0 (records are never materialised on this path).
func reheaderStream(in *bufio.Reader, out io.Writer, opts ReheaderOptions) (int, error) {
	// Read only the leading `#`-prefixed header lines; leave the body in the
	// reader so it can be streamed straight to the output.
	headerBytes, err := readVCFHeaderLines(in)
	if err != nil {
		return 0, fmt.Errorf("bcftools reheader: %w", err)
	}

	newHeader, err := editVCFHeaderText(headerBytes, opts)
	if err != nil {
		return 0, err
	}

	// Emit: new header, then the body streamed verbatim. Wrap in BGZF for
	// VCF.gz output. in still holds every byte after the header.
	if opts.OutputFormat == OutputVCFGz {
		var bw io.WriteCloser
		if opts.Threads > 1 {
			mw, werr := newBGZFOutput(out, opts.CompressLevel, opts.Threads)
			if werr != nil {
				return 0, werr
			}
			bw = mw
		} else {
			bw = bgzip.NewWriter(out)
		}
		if _, err := bw.Write(newHeader); err != nil {
			bw.Close()
			return 0, err
		}
		if _, err := io.Copy(bw, in); err != nil {
			bw.Close()
			return 0, err
		}
		return 0, bw.Close()
	}

	// Plain VCF output.
	if _, err := out.Write(newHeader); err != nil {
		return 0, err
	}
	if _, err := io.Copy(out, in); err != nil {
		return 0, err
	}
	return 0, nil
}

// editVCFHeaderText applies the requested header edits (-h replace, -f FAI
// contigs, -s sample rename) to the original VCF header text and returns the
// new header text (always '\n'-terminated). It is shared by the plain/BGZF
// stream path and the raw compressed-block passthrough path so both produce an
// identical header.
func editVCFHeaderText(headerBytes []byte, opts ReheaderOptions) ([]byte, error) {
	newHeader := headerBytes
	if opts.HeaderFile != "" {
		fileHdr, err := os.ReadFile(opts.HeaderFile)
		if err != nil {
			return nil, fmt.Errorf("bcftools reheader: load -h %s: %w", opts.HeaderFile, err)
		}
		// Preserve the original #CHROM line (and thus the old sample list)
		// when the replacement header omits it — matches upstream's
		// "keep old samples" fallback.
		if !headerTextHasChrom(fileHdr) {
			fileHdr = appendChromLine(fileHdr, oldChromLine(headerBytes))
		}
		newHeader = ensureTrailingNewline(fileHdr)
	}

	if opts.FaiFile != "" || opts.SamplesFile != "" {
		// -f/-s operate on the vcf.Header; rebuild the header text from the
		// current (post -h) header so their edits compose correctly.
		hdr := parseHeaderText(newHeader)
		var err error
		if opts.FaiFile != "" {
			if hdr, err = applyFaiContigs(hdr, opts.FaiFile); err != nil {
				return nil, fmt.Errorf("bcftools reheader: -f %s: %w", opts.FaiFile, err)
			}
		}
		if opts.SamplesFile != "" {
			mapping, names, err := loadSamplesRename(opts.SamplesFile)
			if err != nil {
				return nil, fmt.Errorf("bcftools reheader: -s %s: %w", opts.SamplesFile, err)
			}
			hdr = renameHeaderSamples(hdr, mapping, names)
		}
		newHeader = serialiseHeaderText(hdr)
	}
	return newHeader, nil
}

// readVCFHeaderLines consumes the leading VCF header from br — the contiguous
// run of lines whose first byte is '#', including the terminating newline of the
// last header line — and returns those bytes. The reader is left positioned at
// the first body byte, so the body can be streamed without buffering. When the
// input has no header line the result is empty and br is left untouched.
func readVCFHeaderLines(br *bufio.Reader) ([]byte, error) {
	var header []byte
	for {
		b, err := br.Peek(1)
		if err == io.EOF {
			return header, nil
		}
		if err != nil {
			return header, err
		}
		if b[0] != '#' {
			return header, nil
		}
		line, err := br.ReadBytes('\n')
		header = append(header, line...)
		if err == io.EOF {
			return header, nil
		}
		if err != nil {
			return header, err
		}
	}
}

// reheaderBCF keeps the decode + re-encode path for BCF input: the record
// dictionary indices depend on the header ordering, so the body cannot be
// copied through verbatim.
func reheaderBCF(in io.Reader, out io.Writer, opts ReheaderOptions) (int, error) {
	hdr, recs, err := readAllBCF(in)
	if err != nil {
		return 0, fmt.Errorf("bcftools reheader: %w", err)
	}

	if opts.HeaderFile != "" {
		newHdr, err := loadHeaderFromFile(opts.HeaderFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools reheader: load -h %s: %w", opts.HeaderFile, err)
		}
		if len(newHdr.Samples) == 0 {
			newHdr.Samples = hdr.Samples
		}
		hdr = newHdr
	}

	if opts.FaiFile != "" {
		if hdr, err = applyFaiContigs(hdr, opts.FaiFile); err != nil {
			return 0, fmt.Errorf("bcftools reheader: -f %s: %w", opts.FaiFile, err)
		}
	}

	if opts.SamplesFile != "" {
		mapping, names, err := loadSamplesRename(opts.SamplesFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools reheader: -s %s: %w", opts.SamplesFile, err)
		}
		hdr = renameHeaderSamples(hdr, mapping, names)
		for _, v := range recs {
			renameVariantSamples(v, mapping, names, hdr.Samples)
		}
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}
	for _, v := range recs {
		if err := w.Write(v); err != nil {
			return 0, err
		}
	}
	return len(recs), w.Flush()
}

// splitVCFHeader splits data into its leading VCF header (the contiguous run of
// lines whose first byte is '#', including the terminating newline of the last
// header line) and the remaining body bytes. When data has no header line the
// header is empty and the whole input is the body.
func splitVCFHeader(data []byte) (header, body []byte) {
	pos := 0
	for pos < len(data) && data[pos] == '#' {
		nl := bytes.IndexByte(data[pos:], '\n')
		if nl < 0 {
			// Header line runs to EOF with no newline: the whole file is header.
			return data, nil
		}
		pos += nl + 1
	}
	return data[:pos], data[pos:]
}

// parseHeaderText parses a VCF header text block into a vcf.Header, keeping the
// meta lines verbatim and extracting the sample list from the #CHROM line.
func parseHeaderText(header []byte) *vcf.Header {
	out := &vcf.Header{}
	sc := bufio.NewScanner(bytes.NewReader(header))
	sc.Buffer(make([]byte, 0, 64<<10), 64<<20)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "##") {
			out.MetaInfo = append(out.MetaInfo, line)
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			fields := strings.Split(line, "\t")
			if len(fields) > 9 {
				out.Samples = append([]string{}, fields[9:]...)
			}
		}
	}
	return out
}

// serialiseHeaderText renders a vcf.Header to VCF header text in the same form
// the vcf.Writer emits (meta lines each on their own line, then the #CHROM line,
// each line terminated by '\n').
func serialiseHeaderText(hdr *vcf.Header) []byte {
	var b bytes.Buffer
	for _, meta := range hdr.MetaInfo {
		b.WriteString(meta)
		b.WriteByte('\n')
	}
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	if len(hdr.Samples) > 0 {
		b.WriteString("\tFORMAT\t")
		b.WriteString(strings.Join(hdr.Samples, "\t"))
	}
	b.WriteByte('\n')
	return b.Bytes()
}

// headerTextHasChrom reports whether the header text contains a #CHROM line.
func headerTextHasChrom(header []byte) bool {
	for _, line := range bytes.Split(header, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("#CHROM")) {
			return true
		}
	}
	return false
}

// oldChromLine returns the #CHROM line (without a trailing newline) from a VCF
// header text block, or an empty slice if none is present.
func oldChromLine(header []byte) []byte {
	for _, line := range bytes.Split(header, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("#CHROM")) {
			return line
		}
	}
	return nil
}

// appendChromLine appends chrom (a #CHROM line) to header, ensuring header ends
// with a newline first so the two never merge.
func appendChromLine(header, chrom []byte) []byte {
	out := ensureTrailingNewline(append([]byte{}, header...))
	out = append(out, chrom...)
	out = append(out, '\n')
	return out
}

// ensureTrailingNewline returns b with a single terminating '\n' appended when
// it does not already end in one (and b is non-empty).
func ensureTrailingNewline(b []byte) []byte {
	if len(b) > 0 && b[len(b)-1] != '\n' {
		return append(b, '\n')
	}
	return b
}

// loadHeaderFromFile reads a plain VCF header (lines beginning with `#`)
// from path. Trailing data lines are ignored — they're conceptually
// disallowed in `-h` files anyway.
func loadHeaderFromFile(path string) (*vcf.Header, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	out := &vcf.Header{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 10<<20)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "##") {
			out.MetaInfo = append(out.MetaInfo, line)
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			fields := strings.Split(line, "\t")
			if len(fields) > 9 {
				out.Samples = append([]string{}, fields[9:]...)
			}
			break
		}
	}
	return out, sc.Err()
}

// applyFaiContigs rewrites every `##contig=<ID=...>` line in hdr to match the
// (NAME, LENGTH) entries of a samtools FAI sidecar. New contigs that are
// present in the FAI but missing from the header are appended in FAI order;
// contigs in the header that are absent from the FAI are dropped (this
// matches upstream's `-f` behaviour).
func applyFaiContigs(hdr *vcf.Header, faiPath string) (*vcf.Header, error) {
	f, err := os.Open(faiPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	type entry struct {
		name string
		len  int
	}
	var entries []entry
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 2 {
			return nil, fmt.Errorf("invalid FAI line %q", line)
		}
		ln, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid FAI length %q: %w", fields[1], err)
		}
		entries = append(entries, entry{name: fields[0], len: ln})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	out := &vcf.Header{Samples: hdr.Samples}
	// Insertion strategy: keep every non-contig meta line in its original
	// position; replace the *first* contig block with the FAI-derived
	// lines.
	contigInserted := false
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##contig=") {
			if !contigInserted {
				for _, e := range entries {
					out.MetaInfo = append(out.MetaInfo,
						fmt.Sprintf("##contig=<ID=%s,length=%d>", e.name, e.len))
				}
				contigInserted = true
			}
			continue
		}
		out.MetaInfo = append(out.MetaInfo, m)
	}
	if !contigInserted {
		// Header had no contigs at all — append the FAI block at the end.
		for _, e := range entries {
			out.MetaInfo = append(out.MetaInfo,
				fmt.Sprintf("##contig=<ID=%s,length=%d>", e.name, e.len))
		}
	}
	return out, nil
}

// loadSamplesRename parses the `-s/--samples` rename file. It returns a
// `OLD->NEW` map (non-nil only when the file is in tab-separated form) and a
// flat list of new names (non-nil only when the file is a plain list).
// Comment lines starting with `#` and blank lines are skipped.
func loadSamplesRename(path string) (map[string]string, []string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()
	var names []string
	mapping := map[string]string{}
	tabsSeen := false
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		if strings.IndexByte(line, '\t') >= 0 {
			tabsSeen = true
			fields := strings.SplitN(line, "\t", 2)
			mapping[fields[0]] = strings.TrimSpace(fields[1])
		} else {
			names = append(names, strings.TrimSpace(line))
		}
	}
	if err := sc.Err(); err != nil {
		return nil, nil, err
	}
	if tabsSeen {
		return mapping, nil, nil
	}
	return nil, names, nil
}

// renameHeaderSamples returns a copy of hdr with samples renamed per the
// chosen rule. Exactly one of `mapping` or `names` is set; the other is nil.
func renameHeaderSamples(hdr *vcf.Header, mapping map[string]string, names []string) *vcf.Header {
	out := &vcf.Header{MetaInfo: append([]string{}, hdr.MetaInfo...)}
	if mapping != nil {
		for _, s := range hdr.Samples {
			if n, ok := mapping[s]; ok && n != "" {
				out.Samples = append(out.Samples, n)
			} else {
				out.Samples = append(out.Samples, s)
			}
		}
		return out
	}
	if names != nil {
		// Positional rename — clip to whichever is shorter so we never
		// drop or invent samples.
		out.Samples = make([]string, len(hdr.Samples))
		for i, s := range hdr.Samples {
			if i < len(names) && names[i] != "" {
				out.Samples[i] = names[i]
			} else {
				out.Samples[i] = s
			}
		}
		return out
	}
	out.Samples = append([]string{}, hdr.Samples...)
	return out
}

// renameVariantSamples re-keys each sample on v.Samples so the downstream
// writer sees the new names. We have to keep positions consistent with the
// rewritten header.
func renameVariantSamples(v *vcf.Variant, mapping map[string]string, names []string, newNames []string) {
	if len(v.Samples) == 0 {
		return
	}
	for i, s := range v.Samples {
		if mapping != nil {
			if n, ok := mapping[s.Name]; ok && n != "" {
				v.Samples[i].Name = n
			}
			continue
		}
		if names != nil {
			if i < len(names) && names[i] != "" {
				v.Samples[i].Name = names[i]
			}
			continue
		}
	}
	// As a safety net, ensure positions match the new header sample order.
	if len(newNames) == len(v.Samples) {
		for i := range v.Samples {
			v.Samples[i].Name = newNames[i]
		}
	}
}
