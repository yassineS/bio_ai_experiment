package bcf

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// Writer encodes BCF records to an io.Writer. The first call to WriteHeader (or
// the implicit one inside Write) emits the BCF magic and the length-prefixed
// text header; subsequent calls to Write emit one (l_shared, l_indiv, body)
// frame per record.
//
// Writer does not own the underlying io.Writer. Call Flush before discarding
// the Writer; pkg/bioformats/iohelper and tools/bgzip handle compression at
// the layer above us.
type Writer struct {
	bw      *bufio.Writer
	header  *Header
	wroteHd bool
	// chromIndex maps contig name → dictionary index. Built once from the
	// header so Write can convert vcf.Variant.Chrom (a string) into the
	// int32 the wire wants.
	chromIndex map[string]int32
	// infoIndex / fmtIndex do the same for INFO/FILTER keys and FORMAT keys.
	infoIndex map[string]int32
	fmtIndex  map[string]int32
}

// NewWriter returns a Writer that emits BCF for the given header. The header
// is copied by reference; mutating it after construction may produce invalid
// output.
func NewWriter(w io.Writer, header *Header) *Writer {
	wr := &Writer{
		bw:         bufio.NewWriter(w),
		header:     header,
		chromIndex: make(map[string]int32, len(header.Contigs)),
		infoIndex:  make(map[string]int32, len(header.InfoTags)),
		fmtIndex:   make(map[string]int32, len(header.FmtTags)),
	}
	for i, c := range header.Contigs {
		wr.chromIndex[c.ID] = int32(i)
	}
	for i, t := range header.InfoTags {
		wr.infoIndex[t.ID] = int32(i)
	}
	for i, t := range header.FmtTags {
		wr.fmtIndex[t.ID] = int32(i)
	}
	return wr
}

// NewWriterFromVCFHeader builds a BCF Header from a vcf.Header (parsed text
// metadata + sample list) and returns a Writer pre-configured for it. The
// returned Header is also accessible via Writer.Header() for callers that need
// it (e.g. to look up dictionary indices manually).
func NewWriterFromVCFHeader(w io.Writer, vh *vcf.Header) (*Writer, error) {
	if vh == nil {
		return nil, fmt.Errorf("bcf: nil vcf header")
	}
	text := buildBCFTextHeader(vh)
	hdr, err := parseTextHeader(text)
	if err != nil {
		return nil, err
	}
	return NewWriter(w, hdr), nil
}

// Header returns the Header the Writer is bound to.
func (w *Writer) Header() *Header { return w.header }

// WriteHeader emits the BCF magic and the length-prefixed VCF-style text
// header. Calling it twice is a no-op.
func (w *Writer) WriteHeader() error {
	if w.wroteHd {
		return nil
	}
	if _, err := w.bw.Write(Magic[:]); err != nil {
		return err
	}
	text := w.header.Text
	if text == "" {
		// Synthesize a minimal text header from VCF metadata if Text was not
		// preserved (e.g. the Header was hand-built rather than parsed).
		text = buildBCFTextHeader(w.header.VCF)
	}
	// htslib NUL-terminates the text portion so consumers can use C strings.
	withNul := text + "\x00"
	if err := binary.Write(w.bw, binary.LittleEndian, uint32(len(withNul))); err != nil {
		return err
	}
	if _, err := w.bw.WriteString(withNul); err != nil {
		return err
	}
	w.wroteHd = true
	return nil
}

// Write encodes one VCF variant as a BCF record and writes it to the stream.
// The variant is converted to wire form (1-based POS → 0-based, allele text →
// length-prefixed typed char vectors, INFO/FORMAT names → dictionary indices)
// and emitted as the (l_shared, l_indiv, body) triple defined by the spec.
func (w *Writer) Write(v *vcf.Variant) error {
	if !w.wroteHd {
		if err := w.WriteHeader(); err != nil {
			return err
		}
	}
	shared, indiv, err := w.encodeRecord(v)
	if err != nil {
		return err
	}
	if err := binary.Write(w.bw, binary.LittleEndian, uint32(len(shared))); err != nil {
		return err
	}
	if err := binary.Write(w.bw, binary.LittleEndian, uint32(len(indiv))); err != nil {
		return err
	}
	if _, err := w.bw.Write(shared); err != nil {
		return err
	}
	if _, err := w.bw.Write(indiv); err != nil {
		return err
	}
	return nil
}

// WriteRecord encodes an already-decoded BCF Record. This is the round-trip
// path: it bypasses the vcf.Variant translation and re-emits the dictionary
// indices verbatim.
func (w *Writer) WriteRecord(r *Record) error {
	if !w.wroteHd {
		if err := w.WriteHeader(); err != nil {
			return err
		}
	}
	shared, indiv, err := encodeRecordRaw(r)
	if err != nil {
		return err
	}
	if err := binary.Write(w.bw, binary.LittleEndian, uint32(len(shared))); err != nil {
		return err
	}
	if err := binary.Write(w.bw, binary.LittleEndian, uint32(len(indiv))); err != nil {
		return err
	}
	if _, err := w.bw.Write(shared); err != nil {
		return err
	}
	if _, err := w.bw.Write(indiv); err != nil {
		return err
	}
	return nil
}

// Flush drains the buffered output.
func (w *Writer) Flush() error { return w.bw.Flush() }

// buildBCFTextHeader returns the VCF-style header text that a freshly built
// BCF file should carry. We use the MetaInfo lines verbatim and reconstruct
// the #CHROM line from the sample list.
func buildBCFTextHeader(vh *vcf.Header) string {
	if vh == nil {
		return "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	}
	var sb strings.Builder
	sawFileformat := false
	for _, m := range vh.MetaInfo {
		if strings.HasPrefix(m, "##fileformat=") {
			sawFileformat = true
		}
		sb.WriteString(m)
		sb.WriteByte('\n')
	}
	if !sawFileformat {
		// htslib refuses to parse a BCF whose text lacks fileformat=; emit one.
		sb.WriteString("##fileformat=VCFv4.2\n")
	}
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	if len(vh.Samples) > 0 {
		sb.WriteString("\tFORMAT")
		for _, s := range vh.Samples {
			sb.WriteByte('\t')
			sb.WriteString(s)
		}
	}
	sb.WriteByte('\n')
	return sb.String()
}
