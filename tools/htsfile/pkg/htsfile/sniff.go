// Package htsfile implements format detection for the bioinformatics
// file types this repo can handle: SAM, BAM, CRAM, VCF, BCF, FASTA,
// FASTQ, BED, GFF, plus the BGZF / plain-gzip / raw-text containers
// they ride in. It mirrors htslib's `htsfile` binary at a behavioural
// level — same one-line summary form, same set of formats — without
// linking against libhts.
//
// The public surface is intentionally small: callers either feed a
// path to `Identify` or a pre-opened `io.ReaderAt` to `IdentifyReader`,
// and get back a `*Format` describing the detected container +
// payload. The sniff peeks at a small prefix (≤64 KiB) and never
// consumes more than one BGZF block to make the decision.
package htsfile

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// Compression labels the on-disk wrapping of a file's payload.
type Compression int

const (
	// CompressionPlain is uncompressed text.
	CompressionPlain Compression = iota
	// CompressionGzip is a plain gzip stream (no BGZF block structure).
	CompressionGzip
	// CompressionBGZF is the htslib block-gzip variant — gzip with
	// the BSIZE/"BC" subfield in the extra field of every block.
	CompressionBGZF
	// CompressionUnknown means the bytes were unreadable or the
	// container couldn't be sniffed.
	CompressionUnknown
)

func (c Compression) String() string {
	switch c {
	case CompressionPlain:
		return "plain"
	case CompressionGzip:
		return "gzip-compressed"
	case CompressionBGZF:
		return "BGZF-compressed"
	default:
		return "unknown-compression"
	}
}

// Payload identifies the bioinformatics format carried by a file.
type Payload int

const (
	// PayloadUnknown means the sniff couldn't decide.
	PayloadUnknown Payload = iota
	PayloadSAM
	PayloadBAM
	PayloadCRAM
	PayloadVCF
	PayloadBCF
	PayloadFASTA
	PayloadFASTQ
	PayloadBED
	PayloadGFF
	// PayloadText means a plain ASCII/UTF-8 text file we couldn't
	// classify more specifically.
	PayloadText
	// PayloadBinary means non-text bytes that don't match any
	// known format.
	PayloadBinary
)

func (p Payload) String() string {
	switch p {
	case PayloadSAM:
		return "SAM"
	case PayloadBAM:
		return "BAM"
	case PayloadCRAM:
		return "CRAM"
	case PayloadVCF:
		return "VCF"
	case PayloadBCF:
		return "BCF"
	case PayloadFASTA:
		return "FASTA"
	case PayloadFASTQ:
		return "FASTQ"
	case PayloadBED:
		return "BED"
	case PayloadGFF:
		return "GFF"
	case PayloadText:
		return "text"
	case PayloadBinary:
		return "binary"
	default:
		return "unknown"
	}
}

// Format is the structured result of an htsfile identification.
type Format struct {
	Compression Compression
	Payload     Payload
	// Version, when known, is the format version string the file
	// declares in its header (e.g. "1.6" for SAM, "4.2" for VCF,
	// "3.0" for CRAM). Empty when the format doesn't carry a
	// version field or the field wasn't reached by the sniff.
	Version string
}

// Describe returns the one-line description htslib's htsfile would
// print for this format. The shape is:
//
//	<format-name> [version <X.Y> ]<compression> <kind>
//
// where <kind> is "sequence data" (SAM/BAM/CRAM/FASTA/FASTQ),
// "variant calling data" (VCF/BCF), "genomic interval data" (BED/GFF),
// or "text" / "binary" for the catch-alls.
func (f *Format) Describe() string {
	if f == nil {
		return "unknown"
	}
	kind := payloadKind(f.Payload)
	var b strings.Builder
	b.WriteString(f.Payload.String())
	if f.Version != "" {
		b.WriteString(" version ")
		b.WriteString(f.Version)
	}
	b.WriteByte(' ')
	b.WriteString(f.Compression.String())
	b.WriteByte(' ')
	b.WriteString(kind)
	return b.String()
}

func payloadKind(p Payload) string {
	switch p {
	case PayloadSAM, PayloadBAM, PayloadCRAM, PayloadFASTA, PayloadFASTQ:
		return "sequence data"
	case PayloadVCF, PayloadBCF:
		return "variant calling data"
	case PayloadBED, PayloadGFF:
		return "genomic interval data"
	case PayloadText:
		return "text"
	default:
		return "data"
	}
}

// Identify opens the file at path and sniffs its format.
func Identify(path string) (*Format, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("htsfile: open %s: %w", path, err)
	}
	defer f.Close()
	return IdentifyReader(f)
}

// IdentifyReader sniffs the format of a reader. The reader must
// support enough buffered reads for the sniff (typically the first
// ~256 bytes uncompressed, more if the container is bgzipped).
func IdentifyReader(r io.Reader) (*Format, error) {
	br := bufio.NewReaderSize(r, 1<<16)
	prefix, err := br.Peek(18)
	if err != nil && err != io.EOF && err != bufio.ErrBufferFull {
		// Even short files might be valid (e.g. empty FASTA); tolerate
		// io.ErrUnexpectedEOF by falling through to the binary fallback.
		if len(prefix) == 0 {
			return nil, fmt.Errorf("htsfile: read prefix: %w", err)
		}
	}
	c := detectCompression(prefix)
	switch c {
	case CompressionPlain:
		return identifyPlain(br)
	case CompressionGzip, CompressionBGZF:
		return identifyCompressed(br, c)
	}
	// Unknown: try the text path anyway — a short file with a hint
	// in its first byte (e.g. ">", "@", "##fileformat=") still gets
	// classified, with Compression marked unknown.
	return identifyPlain(br)
}

// detectCompression inspects the first ~18 bytes for the gzip magic
// and the BGZF "BC" extra subfield. Plain text is everything else.
func detectCompression(p []byte) Compression {
	if len(p) >= 2 && p[0] == 0x1f && p[1] == 0x8b {
		// Gzip family. Look for the BGZF BC subfield in the extra
		// field: the gzip header byte 3 carries the FEXTRA flag (0x04);
		// bytes 12-13 begin the extra fields with the BC tag.
		if len(p) >= 18 && (p[3]&0x04) != 0 && p[12] == 'B' && p[13] == 'C' {
			return CompressionBGZF
		}
		return CompressionGzip
	}
	if len(p) == 0 {
		return CompressionUnknown
	}
	return CompressionPlain
}

func identifyCompressed(br *bufio.Reader, c Compression) (*Format, error) {
	// gzip.NewReader handles both plain gzip and BGZF (the BC subfield
	// is opaque to compress/gzip). Decode just enough of the first
	// block to look at the payload prefix.
	zr, err := gzip.NewReader(br)
	if err != nil {
		return nil, fmt.Errorf("htsfile: gzip reader: %w", err)
	}
	defer zr.Close()
	prefix := make([]byte, 4096)
	n, err := io.ReadFull(zr, prefix)
	if err != nil && err != io.EOF && err != io.ErrUnexpectedEOF {
		return nil, fmt.Errorf("htsfile: gzip read: %w", err)
	}
	f := classifyPayload(prefix[:n])
	f.Compression = c
	return f, nil
}

func identifyPlain(br *bufio.Reader) (*Format, error) {
	prefix, _ := br.Peek(4096)
	f := classifyPayload(prefix)
	if f.Compression == CompressionUnknown {
		f.Compression = CompressionPlain
	}
	return f, nil
}

// classifyPayload inspects the (already-decompressed) prefix bytes and
// fills in Payload + Version. Compression is left as its zero value
// (CompressionPlain) for the caller to override.
func classifyPayload(prefix []byte) *Format {
	if len(prefix) == 0 {
		return &Format{Payload: PayloadUnknown}
	}

	// Binary magic-byte formats first.
	if bytes.HasPrefix(prefix, []byte("BAM\x01")) {
		return &Format{Payload: PayloadBAM}
	}
	if bytes.HasPrefix(prefix, []byte("BCF\x02\x02")) {
		return &Format{Payload: PayloadBCF, Version: "2.2"}
	}
	if bytes.HasPrefix(prefix, []byte("BCF\x02\x01")) {
		return &Format{Payload: PayloadBCF, Version: "2.1"}
	}
	if bytes.HasPrefix(prefix, []byte("CRAM")) && len(prefix) >= 6 {
		// CRAM magic is "CRAM" + uint8 major + uint8 minor.
		ver := fmt.Sprintf("%d.%d", prefix[4], prefix[5])
		return &Format{Payload: PayloadCRAM, Version: ver}
	}

	// Text formats: sniff the first non-blank line.
	line := firstLine(prefix)
	if line == "" {
		// Empty file — call it text rather than guess.
		return &Format{Payload: PayloadText}
	}

	switch {
	case strings.HasPrefix(line, "##fileformat=VCFv"):
		return &Format{Payload: PayloadVCF, Version: strings.TrimPrefix(line, "##fileformat=VCFv")}
	case strings.HasPrefix(line, "@HD") || strings.HasPrefix(line, "@SQ") || strings.HasPrefix(line, "@RG") || strings.HasPrefix(line, "@PG") || strings.HasPrefix(line, "@CO"):
		return &Format{Payload: PayloadSAM, Version: samVersion(prefix)}
	case strings.HasPrefix(line, "##gff-version"):
		v := strings.TrimSpace(strings.TrimPrefix(line, "##gff-version"))
		return &Format{Payload: PayloadGFF, Version: v}
	case strings.HasPrefix(line, ">"):
		return &Format{Payload: PayloadFASTA}
	case strings.HasPrefix(line, "@") && looksLikeFASTQ(prefix):
		return &Format{Payload: PayloadFASTQ}
	}

	// BED heuristic: first non-comment line is 3+ whitespace-separated
	// columns where columns 2 and 3 parse as non-negative integers.
	if looksLikeBED(prefix) {
		return &Format{Payload: PayloadBED}
	}

	// Fall through to text if it's ASCII-looking, else binary.
	if isMostlyText(prefix) {
		return &Format{Payload: PayloadText}
	}
	return &Format{Payload: PayloadBinary}
}

// firstLine returns the first non-blank line in p (without the
// trailing newline) or empty if there isn't one within `prefix`.
func firstLine(p []byte) string {
	for len(p) > 0 {
		idx := bytes.IndexByte(p, '\n')
		var line []byte
		if idx < 0 {
			line = p
			p = nil
		} else {
			line = p[:idx]
			p = p[idx+1:]
		}
		// Strip CR for CRLF inputs.
		if len(line) > 0 && line[len(line)-1] == '\r' {
			line = line[:len(line)-1]
		}
		if len(bytes.TrimSpace(line)) > 0 {
			return string(line)
		}
	}
	return ""
}

// samVersion returns the SAM major.minor from the @HD VN: field, if
// present in the prefix. Returns "" when no @HD line is visible.
func samVersion(p []byte) string {
	for _, ln := range bytes.SplitN(p, []byte("\n"), 4) {
		if !bytes.HasPrefix(ln, []byte("@HD")) {
			continue
		}
		for _, field := range bytes.Split(ln, []byte("\t")) {
			if bytes.HasPrefix(field, []byte("VN:")) {
				return string(field[3:])
			}
		}
	}
	return ""
}

// looksLikeFASTQ checks for the @id / sequence / + / quality 4-line
// pattern in the prefix. The third line must start with '+' or be a
// FASTQ separator; the second and fourth lines must be the same
// length.
func looksLikeFASTQ(p []byte) bool {
	lines := bytes.SplitN(p, []byte("\n"), 5)
	if len(lines) < 4 {
		return false
	}
	if len(lines[0]) == 0 || lines[0][0] != '@' {
		return false
	}
	if len(lines[2]) == 0 || lines[2][0] != '+' {
		return false
	}
	if len(lines[1]) != len(lines[3]) {
		return false
	}
	return true
}

// looksLikeBED checks the first non-comment, non-blank line for the
// chrom/start/end shape.
func looksLikeBED(p []byte) bool {
	for _, ln := range bytes.SplitN(p, []byte("\n"), 8) {
		ln = bytes.TrimRight(ln, "\r")
		if len(ln) == 0 || ln[0] == '#' || bytes.HasPrefix(ln, []byte("track")) || bytes.HasPrefix(ln, []byte("browser")) {
			continue
		}
		fields := bytes.Fields(ln)
		if len(fields) < 3 {
			return false
		}
		if !isNonNegInt(fields[1]) || !isNonNegInt(fields[2]) {
			return false
		}
		return true
	}
	return false
}

func isNonNegInt(b []byte) bool {
	if len(b) == 0 {
		return false
	}
	for _, c := range b {
		if c < '0' || c > '9' {
			return false
		}
	}
	return true
}

func isMostlyText(p []byte) bool {
	if len(p) == 0 {
		return false
	}
	textish := 0
	for _, c := range p {
		if c == 0 {
			return false
		}
		if c == '\t' || c == '\n' || c == '\r' || (c >= 0x20 && c <= 0x7e) {
			textish++
		}
	}
	return textish*100/len(p) >= 95
}
