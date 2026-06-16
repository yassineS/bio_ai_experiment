// Binary alignment output for `bedtools intersect` when the query (-a) is a BAM
// or CRAM file and -bed is NOT given. Upstream writes the intersecting
// ALIGNMENTS back out as binary by default — the surviving original alignment
// records, in input order, under the alignment-level flags that make sense for
// binary output (-u / -v / -wa / default). This path preserves the original
// SAM/BAM/CRAM records verbatim (re-encoded by pkg/htsgo/sam's BAMWriter or
// pkg/htsgo/alnio's CRAM writer) and emits the original header, rather than the
// BED12 projection used by the text-output path.
//
// Output framing mirrors upstream's RecordOutputMgr exactly:
//
//   - A BAM query writes BAM.
//   - A CRAM query writes CRAM, but ONLY when a CRAM reference is available
//     (upstream gates CRAM output on --cram-ref / CRAM_REFERENCE: its writer
//     opens "wc" when a reference is set and "wb" otherwise). Without a
//     reference a CRAM query writes BAM, exactly like upstream.
//   - -ubam forces uncompressed BAM regardless of the query format.
//
// Flag gating mirrors upstream ContextIntersect::isValidState exactly:
//
//   - -c / -C (writeCount): error "writeCount option is not valid with BAM
//     query input, unless bed output is specified with -bed option."
//   - -wo / -wao (writeOverlap / writeAllOverlap): error "writeAllOverlap
//     option is not valid with BAM query input, unless bed output is specified
//     with -bed option."
//   - -wb / -loj: a stderr WARNING; the flags are ignored and output proceeds
//     as BAM.
//   - -header: a stderr WARNING; ignored.
//   - -u / -v / -wa / default: produce BAM output.
//
// These gating decisions are made in the CLI (cmd/bedintersect/main.go) before
// this function runs; this function implements only the surviving-record
// selection and BAM encoding.
package bedintersect

import (
	"bufio"
	"bytes"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// QueryFormat classifies the query (-a) stream's alignment format, which
// upstream `bedtools intersect` uses to decide the output type. A text query
// (BED/GFF/VCF) writes text; a BAM or CRAM query writes binary alignments by
// default (BAM, or CRAM when a CRAM reference is available — see the CLI's
// output-format selection).
type QueryFormat int

const (
	// QueryText is a BED/GFF/VCF (or BGZF-compressed text) query: text output.
	QueryText QueryFormat = iota
	// QueryBAM is a BAM query: binary BAM output by default.
	QueryBAM
	// QueryCRAM is a CRAM query: CRAM output when a reference is available,
	// otherwise BAM (mirroring upstream, whose CRAM output is gated on a
	// CRAM reference being set).
	QueryCRAM
)

// ClassifyQueryInput sniffs r and reports its alignment QueryFormat, returning a
// replacement reader that re-yields the bytes consumed while sniffing (so the
// returned reader must be used in place of r). It is the query-format probe the
// CLI uses to decide the output type (text vs BAM vs CRAM), distinguishing BAM
// from CRAM so a CRAM query can be re-emitted as CRAM the way upstream does.
//
// The probe deliberately does NOT classify "any BGZF stream" as BAM: a
// BGZF-compressed BED/GFF/VCF (e.g. a `.bed.gz` piped to `-a -`, which iohelper
// does not transparently decompress for stdin) is still text and must produce
// text output. So a gzip/BGZF stream is QueryBAM only when its decompressed
// prefix is the "BAM\1" magic, mirroring readInRecords' own detection. CRAM is
// recognised by its four-byte `CRAM` magic; a raw (already-decompressed) BAM
// body by the "BAM\1" magic. Everything else (SAM text, plain BED/GFF/VCF) is
// QueryText.
func ClassifyQueryInput(r io.Reader) (QueryFormat, io.Reader, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	magic, _ := br.Peek(4)
	switch {
	case len(magic) >= 4 && string(magic) == "CRAM":
		return QueryCRAM, br, nil
	case len(magic) >= 4 && string(magic) == "BAM\x01":
		return QueryBAM, br, nil
	case len(magic) >= 4 && magic[0] == 0x1f && magic[1] == 0x8b:
		// gzip/BGZF: BAM is BGZF-wrapped, so a BAM stream's decompressed prefix is
		// "BAM\1"; a plain gzipped/BGZF text file is not. Buffer the whole stream
		// (so the returned reader can replay it) and probe the decompressed prefix.
		buf, err := io.ReadAll(br)
		if err != nil {
			return QueryText, bytes.NewReader(buf), err
		}
		if isGzippedBAM(buf) {
			return QueryBAM, bytes.NewReader(buf), nil
		}
		return QueryText, bytes.NewReader(buf), nil
	default:
		return QueryText, br, nil
	}
}

// IsBAMOrCRAMInput reports whether r is a BAM or CRAM alignment stream, returning
// a replacement reader that re-yields the bytes consumed while sniffing (so the
// returned reader must be used in place of r). It is a thin wrapper over
// ClassifyQueryInput preserved for callers that only need the binary/text
// distinction: it is true for both BAM and CRAM queries.
func IsBAMOrCRAMInput(r io.Reader) (bool, io.Reader, error) {
	format, rr, err := ClassifyQueryInput(r)
	return format != QueryText, rr, err
}

// AlnOutputFormat selects the binary framing used to re-emit the surviving
// alignments of a BAM/CRAM query (-a) without -bed.
type AlnOutputFormat int

const (
	// OutputBAM writes BGZF-compressed BAM. It is the default for a BAM query,
	// and the fallback for a CRAM query with no reference available.
	OutputBAM AlnOutputFormat = iota
	// OutputCRAM writes CRAM. Upstream selects it for a CRAM query only when a
	// CRAM reference is available (--cram-ref / CRAM_REFERENCE).
	OutputCRAM
)

// AlnOutputOptions configures the binary alignment-output path.
type AlnOutputOptions struct {
	// Format selects BAM or CRAM framing for the surviving alignments.
	Format AlnOutputFormat
	// Uncompressed selects uncompressed (level-0) BAM output, the effect of
	// upstream `bedtools intersect`'s -ubam flag. The BAM is still BGZF-framed,
	// just with stored DEFLATE blocks. It is a no-op for CRAM output (whose
	// framing -ubam does not affect upstream either).
	Uncompressed bool
	// ReferenceFASTA names the CRAM decode reference (from --cram-ref or the
	// CRAM_REFERENCE environment variable). It reconstructs reference-backed
	// bases when reading a reference-compressed CRAM query; an empty value
	// decodes reference-derived bases as 'N' (the documented fallback). It does
	// not affect a BAM query, which carries its sequence inline.
	ReferenceFASTA string
}

// IntersectBinaryOutput intersects a BAM/CRAM query (-a) against one or more B
// files and writes the surviving alignment records back out as binary BAM or
// CRAM (per out.Format), matching upstream `bedtools intersect`'s default
// behaviour. readerA must be a BAM or CRAM stream (the caller has already
// classified it); readersB may be any supported B input (BED/GFF/VCF/BAM/CRAM).
// The original alignment header and the surviving alignment records are written
// verbatim, in input order.
//
// Only the alignment-level output modes are honoured here (selected by opts):
// the default and -wa/-u write each A alignment that has at least one overlap
// once; -v writes each A alignment that has no overlap. -wb/-loj are ignored
// (the caller has already emitted the upstream warning); -c/-C/-wo/-wao are
// rejected by the caller before reaching this path. It returns the number of
// alignment records written.
func IntersectBinaryOutput(readerA io.Reader, readersB []io.Reader, writer io.Writer, opts IntersectOptions, out AlnOutputOptions) (int, error) {
	bRecords, _, _, err := readAllB(readersB)
	if err != nil {
		return 0, err
	}
	finder := newFinder(bRecords, opts)

	aReader, err := alnio.NewReaderWithReference(readerA, out.ReferenceFASTA)
	if err != nil {
		return 0, fmt.Errorf("error opening BAM/CRAM input: %w", err)
	}

	hdr := aReader.Header()
	w, err := newAlnWriter(writer, out.Format, out.Uncompressed)
	if err != nil {
		return 0, fmt.Errorf("error opening %s output: %w", out.Format, err)
	}
	if err := w.WriteHeader(hdr); err != nil {
		return 0, fmt.Errorf("error writing %s header: %w", out.Format, err)
	}

	count := 0
	for {
		rec, err := aReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, fmt.Errorf("error reading BAM/CRAM input: %w", err)
		}
		a := aRecordFromAlignment(rec)
		var hits []rawHit
		if !a.unmapped {
			hits = finder.overlaps(a, opts)
		}
		if !bamRecordSurvives(len(hits), opts) {
			continue
		}
		if err := w.Write(rec); err != nil {
			return count, fmt.Errorf("error writing %s record: %w", out.Format, err)
		}
		count++
	}
	if err := w.Close(); err != nil {
		return count, fmt.Errorf("error finalising %s output: %w", out.Format, err)
	}
	return count, nil
}

// newAlnWriter constructs the sam.Writer matching format: a BGZF BAM writer for
// OutputBAM (level-0 stored blocks when uncompressed is set, matching upstream
// -ubam), or alnio's CRAM v3.0 writer for OutputCRAM. Both satisfy the
// WriteHeader/Write/Close shape this path drives. It returns an error only if
// the BAM writer cannot be constructed.
func newAlnWriter(w io.Writer, format AlnOutputFormat, uncompressed bool) (sam.Writer, error) {
	if format == OutputCRAM {
		return alnio.NewCRAMWriter(w), nil
	}
	return sam.NewBAMWriterOptions(w, sam.BAMWriterOptions{Uncompressed: uncompressed})
}

// String renders the output format for diagnostics.
func (f AlnOutputFormat) String() string {
	if f == OutputCRAM {
		return "CRAM"
	}
	return "BAM"
}

// IntersectBAMOutput is IntersectBinaryOutput pinned to BAM output. It is kept
// for callers (and tests) that only need the BAM path; new callers should use
// IntersectBinaryOutput so CRAM output can be selected.
func IntersectBAMOutput(readerA io.Reader, readersB []io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	return IntersectBinaryOutput(readerA, readersB, writer, opts, AlnOutputOptions{Format: OutputBAM})
}

// aRecordFromAlignment builds the interval model used for overlap detection from
// a single alignment, reusing the BED12 projection helpers so the overlap math
// (CIGAR blocks, strand, -split block sums) is identical to the text-output
// path. Unmapped alignments yield a record flagged unmapped (never overlaps),
// matching upstream's printUnmapped handling.
func aRecordFromAlignment(rec *sam.Record) *inRecord {
	if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" {
		return unmappedBAMRecord(rec)
	}
	return bamToBED12(rec)
}

// bamRecordSurvives reports whether an A alignment with the given overlap count
// is written under the active BAM-output mode, mirroring upstream
// RecordOutputMgr's selection for BAM output:
//
//   - -v (noHit): written when it has no overlap.
//   - default / -wa / -u: written when it has at least one overlap (written
//     once regardless of the number of hits).
//
// -wb and -loj are ignored for BAM output (they fall through to the default
// "any overlap" rule), and -c/-C/-wo/-wao never reach this path.
func bamRecordSurvives(numHits int, opts IntersectOptions) bool {
	if opts.NoOverlap {
		return numHits == 0
	}
	return numHits > 0
}
