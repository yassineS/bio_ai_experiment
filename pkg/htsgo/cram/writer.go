package cram

import (
	"bytes"
	"crypto/md5"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"hash/crc32"
	"io"
	"os"
	"runtime"
	"sync"

	kgzip "github.com/klauspost/compress/gzip"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram/codec"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// Version selects the CRAM format version a RecordWriter emits. Both
// versions share the v3 container, slice and record layout and the v3
// CRC32 fields; they differ only in the per-block compression codecs the
// writer is allowed to use.
type Version int

const (
	// VersionV30 is CRAM v3.0: the writer compresses each block with raw,
	// gzip or bzip2. It is the default and the format every existing
	// caller gets.
	VersionV30 Version = iota
	// VersionV31 is CRAM v3.1: in addition to raw, gzip and bzip2 the
	// writer may compress a block with the rANS 4x16 codec (block method
	// 5), the distinguishing capability of v3.1.
	VersionV31
	// VersionV40 is CRAM v4.0 (the draft format). It keeps the v3
	// container/block/slice tree and CRC32 fields but encodes every
	// variable-length integer as a uint7 LEB128 varint rather than ITF-8 /
	// LTF-8, widens the alignment coordinates to 64-bit, stores its integer
	// data series through the VARINT_UNSIGNED / VARINT_SIGNED codecs (rather
	// than EXTERNAL-of-ITF8), and terminates the file with the distinct
	// 31-byte v4 EOF marker. Block compression is the v3.0 set
	// (raw/gzip/bzip2): the v4 framing change is orthogonal to which block
	// codec is chosen.
	VersionV40
)

// fileDefinition returns the on-disk file-definition version for v. CRAM
// v3.0 and v3.1 carry major version 3 (only the minor differs); v4.0
// carries major version 4.
func (v Version) fileDefinition() FileDefinition {
	switch v {
	case VersionV31:
		return FileDefinition{Major: 3, Minor: 1}
	case VersionV40:
		return FileDefinition{Major: 4, Minor: 0}
	default:
		return FileDefinition{Major: 3, Minor: 0}
	}
}

// usesUint7 reports whether this writer version frames its
// variable-length integers as uint7 varints (CRAM v4.0) rather than
// ITF-8 / LTF-8 (v2.x / v3.x). It is the single predicate that steers the
// container, block, slice-header, compression-header and data-series
// framing onto the v4 path, mirroring FileDefinition.usesUint7 on the
// decode side.
func (v Version) usesUint7() bool { return v.fileDefinition().Major >= 4 }

// String returns the version as a "major.minor" string.
func (v Version) String() string {
	return v.fileDefinition().VersionString()
}

// defaultRecordsPerSlice caps how many records the writer packs into one
// slice (and, since the writer emits one slice per container, one
// container). A new container is started once this many records have
// been buffered. The value is a balance between container-header
// overhead and per-slice memory; it has no effect on correctness.
const defaultRecordsPerSlice = 10000

// refIDUnset marks RecordWriter.bufRefID as carrying no value: the slice
// buffer is empty, so the next record adopts its reference id. It differs
// from -1 (the genuine unmapped/"*" reference id) and -2 (multi-reference)
// so an empty buffer is never mistaken for an unmapped run.
const refIDUnset = -3

// gzipBlockThreshold is the smallest external-block payload the writer
// bothers to gzip-compress. Below it the deflate header and trailer cost
// more than they save, so the block is written raw (method 0).
const gzipBlockThreshold = 32

// content-id constants. Every external data series the writer emits is
// keyed by a stable small integer; the compression header's encoding map
// points each two-character series at its block by this id. The values
// are an implementation detail — any distinct positive integers would
// do — but keeping them fixed makes the on-disk layout predictable and
// the writer easy to read against the reader.
const (
	cidSAMHeader = 0  // the SAM file-header block.
	cidBF        = 1  // bit flags (SAM FLAG).
	cidCF        = 2  // CRAM per-record flags.
	cidRI        = 3  // reference id (multi-reference slices).
	cidRL        = 4  // read length.
	cidAP        = 5  // alignment position.
	cidRG        = 6  // read group (always -1; RG travels as a tag).
	cidRN        = 7  // read names.
	cidMF        = 8  // mate flags.
	cidNS        = 9  // mate reference id.
	cidNP        = 10 // mate alignment position.
	cidTS        = 11 // template size.
	cidTL        = 12 // tag-line (tag-combination index).
	cidMQ        = 13 // mapping quality.
	cidFN        = 14 // read-feature count.
	cidFC        = 15 // read-feature codes.
	cidFP        = 16 // read-feature positions.
	cidBB        = 17 // base-stretch bytes (feature 'b').
	cidBBLen     = 18 // base-stretch lengths.
	cidBA        = 19 // single bases (unmapped reads).
	cidQS        = 20 // quality scores.
	cidIN        = 21 // inserted bases (feature 'I').
	cidINLen     = 22 // inserted-base lengths.
	cidSC        = 23 // soft-clipped bases (feature 'S').
	cidSCLen     = 24 // soft-clipped lengths.
	cidDL        = 25 // deletion lengths (feature 'D').
	cidRS        = 26 // reference-skip lengths (feature 'N').
	cidPD        = 27 // padding lengths (feature 'P').
	cidHC        = 28 // hard-clip lengths (feature 'H').
	cidBS        = 29 // base-substitution codes (feature 'X', reference-based encoding).
	cidTagBase   = 64 // first content id handed out to auxiliary tag series.
)

// tagContentIDs returns the two content ids an auxiliary tag's
// BYTE_ARRAY_LEN series uses: a length block and a value block. The
// i-th distinct tag gets the pair (cidTagBase+2i, cidTagBase+2i+1).
func tagContentIDs(i int) (lenID, valID int32) {
	base := cidTagBase + int32(i)*2
	return base, base + 1
}

// RecordWriter encodes a stream of sam.Record values as a CRAM v3.0
// file. It mirrors RecordReader: construct one with NewRecordWriter,
// feed records with Write, and finalise the file with Close.
//
// The writer produces reference-free CRAM — a mapped read's bases are
// stored literally rather than diffed against an external reference — so
// the file is fully self-contained and decodes without a reference
// FASTA. Records are buffered and flushed one slice (and one container)
// at a time; Close flushes the final partial slice and appends the CRAM
// EOF marker.
//
// A RecordWriter is not safe for concurrent use.
type RecordWriter struct {
	w      io.Writer
	closer io.Closer
	header *sam.Header

	// version is the CRAM format this writer emits. It is a per-writer
	// field — not package-level state — so two writers targeting
	// different versions can run concurrently. It governs both the
	// file-definition minor version and the per-block codec set.
	version Version

	// binning is the lossy quality-binning scheme applied to every
	// record's QUAL before it reaches the QS data series. BinningNone
	// (the zero value) leaves quality untouched, so the default writer is
	// losslessly exact and existing callers are unaffected.
	binning QualityBinning

	// refProvider supplies reference bases per slice (see
	// WriterOptions.ReferenceProvider). The writer fetches only each slice's
	// coordinate window — never a whole contig — so peak memory is bounded by
	// one slice's span. Accessed only from the single-threaded Write/flush
	// path; each encode job carries its own immutable window slice.
	refProvider ReferenceProvider

	// referencePath is the reference FASTA path passed with -T/--reference,
	// emitted verbatim as the @SQ UR: tag (see WriterOptions.ReferencePath).
	// Empty suppresses UR injection.
	referencePath string

	// precomputedM5 holds the @SQ M5 (reference MD5) tags computed up front in
	// headerForEncode, keyed by contig name. The whole-genome M5 hash is the
	// writer's dominant header-write cost (every reference-present contig is
	// hashed in full), so it is fanned across a worker pool — each worker on its
	// own independent FASTA handle — instead of being hashed serially inside the
	// @SQ loop. A value of "" records a contig whose bases could not be fetched
	// (no M5 emitted, matching upstream). Nil when the precompute did not run
	// (no reference, or the serial-from-augment fallback), in which case
	// augmentSQLine falls back to contigMD5. The map's bytes are identical to
	// the serial result: the worker/handle/goroutine choice never changes a
	// contig's MD5.
	precomputedM5 map[string]string

	// disableM5, when true, suppresses computed @SQ M5 injection for the ENTIRE
	// header: upstream htslib disables M5 for all @SQ (embed_ref=2 fallback) as
	// soon as any single M5-less @SQ contig is unresolvable in the reference. It
	// is computed once in headerForEncode and read by augmentSQLine. Pre-existing
	// M5 tags in the input header are unaffected.
	disableM5 bool

	// refIndex maps a reference name to its zero-based @SQ position, so a
	// record's RName / RNext can be turned into the integer ids the CRAM
	// data series store.
	refIndex map[string]int32

	// buf holds the records of the slice currently being assembled; it is
	// flushed to a container once it reaches recordsPerSlice entries.
	buf []*sam.Record
	// bufRefID is the reference id shared by every record currently in buf,
	// or refIDUnset when buf is empty. The writer flushes the buffer before
	// appending a record whose reference id differs, so every container — and
	// thus every slice — is single-reference. A multi-reference slice
	// (RefSeqID -2) cannot be reconstructed against an external reference (the
	// decoder has no per-record reference span to copy match runs from), so a
	// slice straddling a contig boundary would emit unreconstructable SEQ.
	// Splitting at the boundary also matches htslib, which starts a fresh
	// container at every new reference.
	bufRefID int32
	// recordsPerSlice is the slice-size cap; defaultRecordsPerSlice
	// unless overridden for testing.
	recordsPerSlice int
	// recordCounter is the running total of records written to previous
	// containers, stored in each container and slice header.
	recordCounter int64

	// wroteHeader records whether the file-definition and SAM-header
	// container have been emitted yet.
	wroteHeader bool
	// closed records whether Close has run, so a double Close is a no-op.
	closed bool
	// err latches the first write error; every later call returns it.
	err error

	// encodeThreads is the number of containers encoded concurrently. The
	// encode of a container — dominated by per-block DEFLATE — is the writer's
	// hot path, and containers are independent, so they are farmed to a worker
	// pool and written back in submission order. The emitted bytes are
	// identical for any thread count. 1 keeps the legacy synchronous path.
	encodeThreads int
	// Pipeline state, lazily started on the first asynchronous flush:
	pipeStarted bool
	jobs        chan *encodeJob        // submitted to workers, in order
	ordered     chan chan encodeResult // result channels, in submission order
	pipeDone    chan struct{}          // closed when the writer goroutine exits
	pipeWG      sync.WaitGroup         // tracks encode workers
	pipeErr     error                  // first encode/write error from the pipeline
}

// ReferenceProvider supplies reference bases on demand so a reference-based
// CRAM writer need not hold the whole genome resident. The writer fetches only
// the coordinate window each slice actually covers (a few records' span),
// exactly as upstream htslib streams the reference per slice, so peak memory is
// bounded by one slice's span rather than a whole chromosome. *fasta.RandomAccess
// satisfies it directly. See WriterOptions.ReferenceProvider for the contract.
type ReferenceProvider interface {
	// Length returns the named contig's base count, or 0 if it is absent.
	Length(name string) int64
	// Fetch returns the contig's bases over the half-open, 0-based range
	// [start, end), clamped to the contig. The returned slice is read-only.
	Fetch(name string, start, end int64) ([]byte, error)
}

// mapReferenceProvider adapts an in-memory whole-genome map to the
// ReferenceProvider interface, so the windowing encode path is the single code
// path whether the reference is supplied lazily (the production path) or as a
// resident map (small references and tests).
type mapReferenceProvider map[string][]byte

func (m mapReferenceProvider) Length(name string) int64 { return int64(len(m[name])) }

func (m mapReferenceProvider) Fetch(name string, start, end int64) ([]byte, error) {
	b, ok := m[name]
	if !ok {
		return nil, nil
	}
	if start < 0 {
		start = 0
	}
	if end > int64(len(b)) {
		end = int64(len(b))
	}
	if start >= end {
		return nil, nil
	}
	return b[start:end], nil
}

// encodeJob is one container's worth of work handed to an encode worker. Its
// reference window (refWindow, the bases the slice spans, based at refWindowStart
// 0-based; hasRef marks the contig as present in the reference even when the
// window is empty) is immutable for the job's lifetime, so workers read it
// without synchronising against the writer.
type encodeJob struct {
	records        []*sam.Record
	refWindow      []byte
	refWindowStart int32
	hasRef         bool
	counter        int64
	out            chan encodeResult
}

// encodeResult carries an encoded container (or the error that produced none).
type encodeResult struct {
	data []byte
	err  error
}

// NewRecordWriter returns a RecordWriter that encodes records to w as a
// CRAM v3.0 file. The SAM header is written immediately as the first
// container's single block, so h must be complete before the call. It
// returns an error only if the initial header write to w fails.
//
// To target CRAM v3.1 instead, use NewRecordWriterVersion.
func NewRecordWriter(w io.Writer, h *sam.Header) (*RecordWriter, error) {
	return NewRecordWriterVersion(w, h, VersionV30)
}

// NewRecordWriterVersion returns a RecordWriter that encodes records to w
// as a CRAM file of the requested version. VersionV30 produces a v3.0
// file (raw/gzip block compression); VersionV31 produces a v3.1 file,
// which may additionally compress blocks with the rANS 4x16 codec. The
// SAM header is written immediately, so h must be complete before the
// call. It returns an error only if the initial header write to w fails.
func NewRecordWriterVersion(w io.Writer, h *sam.Header, version Version) (*RecordWriter, error) {
	return NewRecordWriterOpts(w, h, WriterOptions{Version: version})
}

// WriterOptions configures a RecordWriter at construction time. Its zero
// value is the default writer: CRAM v3.0 with no quality binning, exactly
// what NewRecordWriter produces.
//
// Construction-time configuration is required for any option that
// influences the CRAM file header — the header container is the file's
// first container and is written immediately by the constructor, so it
// must be settable before that write. Quality binning records its
// provenance in the SAM header, so it lives here rather than in a
// post-construction setter.
type WriterOptions struct {
	// Version selects the CRAM format version (v3.0 or v3.1). The zero
	// value, VersionV30, is the default.
	Version Version
	// Binning selects the lossy quality-binning scheme. The zero value,
	// BinningNone, disables binning and keeps the writer losslessly
	// exact. When a real scheme is set, the writer records a @CO
	// provenance line in the embedded SAM header.
	Binning QualityBinning
	// Reference maps a contig name to its full reference bases (upper-case
	// ACGTN). When provided, a mapped read on a known contig is encoded
	// reference-based — only its mismatches are stored, as substitution
	// features — exactly as upstream CRAM does, so the file is far smaller
	// and faster to encode but requires the same reference to decode. When
	// nil (the zero value, what NewRecordWriter produces) the writer stays
	// reference-free: every read's bases are carried literally, so the file
	// is self-contained and decodes without a FASTA.
	//
	// Holding the whole genome resident costs gigabytes for a human
	// reference; ReferenceProvider is the lazy alternative and is preferred
	// for large references.
	Reference map[string][]byte
	// ReferenceProvider lazily supplies one contig's reference bases on
	// demand, so a reference-based writer never holds more than the contig
	// (or two, across a boundary) currently being encoded — matching
	// upstream htslib, which streams the reference per contig rather than
	// loading the whole FASTA. When set it takes precedence over Reference.
	// The provider is invoked only from the writer's single-threaded Write
	// path, at most once per contig change; it must return the contig's full
	// upper-cased bases, or (nil, nil) for a contig absent from the reference
	// (which puts that contig on the reference-free path). The returned slice
	// is treated as read-only.
	ReferenceProvider ReferenceProvider
	// EncodeThreads sets how many containers are encoded concurrently. The
	// emitted file is byte-identical for any value. 0 (the default) auto-sizes
	// to the machine's CPU count (capped); 1 forces the synchronous path.
	EncodeThreads int
	// ReferencePath is the reference FASTA path the file was encoded against
	// (the -T/--reference argument). When non-empty it is written verbatim as
	// the UR: tag on every @SQ line that lacks one, mirroring upstream htslib
	// (cram_write_SAM_hdr's full_path(ref_fn)). The caller passes the path
	// exactly as samtools would emit it — already absolute for the fixtures —
	// so the embedded header byte-matches upstream. Empty suppresses UR
	// injection. It is independent of Reference: Reference supplies the bases
	// for M5 and reference-based encoding, ReferencePath supplies the UR text.
	ReferencePath string
}

// NewRecordWriterOpts returns a RecordWriter that encodes records to w as
// a CRAM file configured by opts. It is the most general RecordWriter
// constructor; NewRecordWriter, NewRecordWriterVersion and the CreateCRAM
// family are thin wrappers over it.
//
// When opts.Binning selects a real binning scheme, a @CO line documenting
// the lossy quality transform is appended to a copy of h before the SAM
// header is written, so a downstream reader can tell the qualities were
// binned. The caller's *sam.Header is not modified.
//
// The SAM header is written immediately, so h must be complete before the
// call. It returns an error only if opts is invalid or the initial header
// write to w fails.
func NewRecordWriterOpts(w io.Writer, h *sam.Header, opts WriterOptions) (*RecordWriter, error) {
	if !opts.Binning.valid() {
		return nil, fmt.Errorf("cram: unknown quality-binning scheme %d", int(opts.Binning))
	}
	if h == nil {
		h = &sam.Header{}
	}
	if opts.Binning != BinningNone {
		// Record the lossy transform in the embedded header. A copy is
		// taken so the caller's header is left untouched.
		h = headerWithBinningProvenance(h, opts.Binning)
	}
	rw := &RecordWriter{
		w:               w,
		header:          h,
		version:         opts.Version,
		binning:         opts.Binning,
		refProvider:     resolveReferenceProvider(opts),
		referencePath:   opts.ReferencePath,
		encodeThreads:   resolveEncodeThreads(opts.EncodeThreads),
		refIndex:        make(map[string]int32, len(h.Refs)),
		recordsPerSlice: defaultRecordsPerSlice,
		bufRefID:        refIDUnset,
	}
	for i, ref := range h.Refs {
		rw.refIndex[ref.Name] = int32(i)
	}
	if err := rw.writeFileHeader(); err != nil {
		return nil, err
	}
	return rw, nil
}

// headerWithBinningProvenance returns a shallow copy of h with one extra
// @CO comment line documenting that lossy quality-score binning was
// applied. The copy shares the slice backing arrays of the original
// except for Lines and Comments, which are freshly allocated, so the
// caller's header is not mutated.
func headerWithBinningProvenance(h *sam.Header, b QualityBinning) *sam.Header {
	note := "samtools/htsgo: lossy quality-score binning applied (scheme " + b.String() + ")"
	cp := *h
	cp.Lines = make([]sam.HeaderLine, len(h.Lines), len(h.Lines)+1)
	copy(cp.Lines, h.Lines)
	cp.Lines = append(cp.Lines, sam.HeaderLine{
		Tag:    "CO",
		Fields: []sam.HeaderField{{Value: note}},
	})
	cp.Comments = make([]string, len(h.Comments), len(h.Comments)+1)
	copy(cp.Comments, h.Comments)
	cp.Comments = append(cp.Comments, note)
	return &cp
}

// CreateCRAM creates the named file and returns a RecordWriter over it
// targeting CRAM v3.0. The caller must call Close, which flushes the
// final container and releases the file handle.
func CreateCRAM(path string, h *sam.Header) (*RecordWriter, error) {
	return CreateCRAMVersion(path, h, VersionV30)
}

// CreateCRAMVersion creates the named file and returns a RecordWriter
// over it targeting the requested CRAM version. The caller must call
// Close, which flushes the final container and releases the file handle.
func CreateCRAMVersion(path string, h *sam.Header, version Version) (*RecordWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	rw, err := NewRecordWriterVersion(f, h, version)
	if err != nil {
		f.Close()
		os.Remove(path)
		return nil, err
	}
	rw.closer = f
	return rw, nil
}

// WriteCRAM encodes header and every record in records to w as a
// complete CRAM v3.0 file, including the trailing EOF marker. It is a
// convenience wrapper over NewRecordWriter / Write / Close.
func WriteCRAM(w io.Writer, h *sam.Header, records []*sam.Record) error {
	return WriteCRAMVersion(w, h, records, VersionV30)
}

// WriteCRAMVersion encodes header and every record in records to w as a
// complete CRAM file of the requested version, including the trailing
// EOF marker. It is a convenience wrapper over NewRecordWriterVersion /
// Write / Close.
func WriteCRAMVersion(w io.Writer, h *sam.Header, records []*sam.Record, version Version) error {
	rw, err := NewRecordWriterVersion(w, h, version)
	if err != nil {
		return err
	}
	for i, rec := range records {
		if err := rw.Write(rec); err != nil {
			return fmt.Errorf("cram: writing record %d: %w", i, err)
		}
	}
	return rw.Close()
}

// Write buffers one record for encoding. Records accumulate until a
// slice is full, at which point a container is assembled and flushed to
// the underlying writer. A record shape the simple writer cannot encode
// (see the package documentation) is rejected here with a clear error
// and nothing is written.
func (rw *RecordWriter) Write(rec *sam.Record) error {
	if rw.err != nil {
		return rw.err
	}
	if rw.closed {
		return fmt.Errorf("cram: Write after Close")
	}
	if rec == nil {
		return fmt.Errorf("cram: cannot write a nil record")
	}
	// Defensive lazy guard: the CRAM→BAM decode passthrough may hand the writer
	// a record carrying RawAux instead of a decoded Aux slice (it should not on
	// the gated path, but the writer buffers the pointer and reads .Aux below /
	// at encode time, so materialise the aux fields here so correctness never
	// depends on the gate). A no-op when Aux is already set.
	rec.MaterialiseAux()
	if err := rw.checkRecord(rec); err != nil {
		return err
	}
	// Keep every slice single-reference: flush the buffer before adding a
	// record that maps to a different reference than the records already
	// buffered. A slice mixing references would be encoded as multi-reference
	// (RefSeqID -2), which the reference-based decode path cannot reconstruct
	// — it has no per-record reference span to copy match runs from, so the
	// matched bases come back as 'N'. htslib likewise starts a fresh container
	// at every reference change.
	id := recordRefID(rec, rw.refIndex)
	if len(rw.buf) > 0 && id != rw.bufRefID {
		if err := rw.flushContainer(); err != nil {
			return err
		}
	}
	if len(rw.buf) == 0 {
		rw.bufRefID = id
	}
	rw.buf = append(rw.buf, rec)
	if len(rw.buf) >= rw.recordsPerSlice {
		return rw.flushContainer()
	}
	return nil
}

// Close flushes the final partial slice, appends the version-appropriate
// CRAM EOF marker (the 38-byte v3 sentinel, or the 31-byte v4 sentinel for
// CRAM v4.0), and releases the file handle when the writer was created by
// CreateCRAM. It is safe to call Close more than once.
func (rw *RecordWriter) Close() error {
	if rw.closed {
		return rw.err
	}
	rw.closed = true
	if rw.err == nil && len(rw.buf) > 0 {
		rw.flushContainer()
	}
	// Drain the async encode pipeline (if any) so every container is written
	// before the EOF marker, and surface the first error it latched.
	if perr := rw.finishPipeline(); perr != nil && rw.err == nil {
		rw.err = perr
	}
	if rw.err == nil {
		marker := eofMarkerV3
		if rw.version.usesUint7() {
			marker = eofMarkerV4
		}
		if _, err := rw.w.Write(marker); err != nil {
			rw.err = fmt.Errorf("cram: writing EOF marker: %w", err)
		}
	}
	if rw.closer != nil {
		if cerr := rw.closer.Close(); cerr != nil && rw.err == nil {
			rw.err = cerr
		}
	}
	return rw.err
}

// checkRecord rejects a record the simple reference-free writer cannot
// faithfully encode. The writer aims for correctness over coverage: a
// shape it cannot round-trip is an explicit error, never silent data
// loss.
func (rw *RecordWriter) checkRecord(rec *sam.Record) error {
	if rec.RName != "" && rec.RName != "*" {
		if _, ok := rw.refIndex[rec.RName]; !ok {
			return fmt.Errorf("cram: record %q references %q, absent from the SAM header @SQ lines",
				rec.QName, rec.RName)
		}
	}
	if rec.RNext != "" && rec.RNext != "*" && rec.RNext != "=" {
		if _, ok := rw.refIndex[rec.RNext]; !ok {
			return fmt.Errorf("cram: record %q mate references %q, absent from the SAM header @SQ lines",
				rec.QName, rec.RNext)
		}
	}
	if rec.Flag&sam.FlagUnmapped == 0 {
		// A mapped record's CIGAR must consume exactly the SEQ length, or
		// the read-feature reconstruction cannot recover the sequence.
		if len(rec.Cigar) > 0 && rec.Seq != "" && rec.Seq != "*" {
			if rec.Cigar.QueryLength() != len(rec.Seq) {
				return fmt.Errorf("cram: record %q CIGAR query length %d does not match SEQ length %d",
					rec.QName, rec.Cigar.QueryLength(), len(rec.Seq))
			}
		}
	}
	for _, a := range rec.Aux {
		if _, err := encodeTagValue(a); err != nil {
			return fmt.Errorf("cram: record %q tag %q: %w", rec.QName, a.Tag, err)
		}
	}
	return nil
}

// writeFileHeader emits the 26-byte CRAM file definition and the
// file-header container that carries the SAM text header as its single
// block.
func (rw *RecordWriter) writeFileHeader() error {
	if rw.wroteHeader {
		return nil
	}
	fd := rw.version.fileDefinition()
	var def [fileDefSize]byte
	copy(def[0:4], fileDefMagic[:])
	def[4] = fd.Major
	def[5] = fd.Minor
	// FileID is left NUL — the conventional originating-file-name slot is
	// optional and a reader trims trailing NULs.
	if _, err := rw.w.Write(def[:]); err != nil {
		return fmt.Errorf("cram: writing file definition: %w", err)
	}

	// The SAM header block payload is a 4-byte little-endian text length
	// followed by the header text, matching what readSAMHeader expects. The
	// header is serialised in htslib's emission order (verbatim input order
	// with @HD hoisted first) so the embedded SAM header byte-matches what
	// upstream samtools writes into a CRAM. When a reference is supplied the
	// @SQ lines are first augmented with the M5 (reference MD5) and UR
	// (reference path) tags upstream injects (cram_write_SAM_hdr).
	text := []byte(rw.headerForEncode().TextCanonical())
	payload := make([]byte, 4+len(text))
	binary.LittleEndian.PutUint32(payload[:4], uint32(len(text)))
	copy(payload[4:], text)
	headerBlock := encodeBlock(rw.version, ContentFileHeader, cidSAMHeader, payload)

	// The file-header container holds exactly the header block. Its
	// reference fields are the unmapped/no-data sentinels.
	body := headerBlock
	hdr := containerHeaderBytes(rw.version, containerFields{
		length:        int32(len(body)),
		refSeqID:      0,
		startPos:      0,
		alignmentSpan: 0,
		numRecords:    0,
		recordCounter: 0,
		numBases:      0,
		numBlocks:     1,
		landmarks:     nil,
	})
	if _, err := rw.w.Write(hdr); err != nil {
		return fmt.Errorf("cram: writing file-header container header: %w", err)
	}
	if _, err := rw.w.Write(body); err != nil {
		return fmt.Errorf("cram: writing file-header container body: %w", err)
	}
	rw.wroteHeader = true
	return nil
}

// headerForEncode returns the SAM header to embed in the CRAM file. When the
// writer has neither a reference map nor a reference path it returns the
// original header unchanged (the reference-free case is left byte-identical).
// Otherwise it returns a copy whose @SQ lines carry the M5 (reference MD5) and
// UR (reference path) tags upstream htslib injects in cram_write_SAM_hdr:
//
//   - M5 is the lower-case hex MD5 of the reference sequence for that contig,
//     hashed exactly as htslib does — over the upper-cased, whitespace-stripped
//     bases. The bases in rw.reference are already supplied upper-cased with
//     newlines removed (the FASTA random-access reader canonicalises them), so
//     md5.Sum over them reproduces the htslib M5 byte-for-byte. M5 is only
//     filled when the contig's bases are available; a contig absent from the
//     reference map is left without M5 (matching upstream, which cannot hash a
//     reference it could not load).
//   - UR is rw.referencePath verbatim — the -T/--reference argument, which the
//     caller passes already as upstream's full_path() would emit it. It is
//     injected onto every @SQ that ends up with an M5 (computed OR pre-existing);
//     an @SQ with no M5 at all (a contig absent from the reference and without a
//     pre-existing M5, encoded reference-free) is left bare, matching upstream.
//
// An @SQ that already carries an M5 or UR is left intact: upstream never
// overwrites an existing tag, only fills the absent one.
func (rw *RecordWriter) headerForEncode() *sam.Header {
	if rw.refProvider == nil && rw.referencePath == "" {
		return rw.header
	}
	// All-or-nothing M5: upstream htslib (cram_io.c cram_write_SAM_hdr) computes
	// an M5 for every @SQ that lacks one, but if ANY M5-less @SQ contig cannot be
	// resolved in the reference (absent, or its bases cannot be loaded) it falls
	// back to embed_ref=2 and emits NO computed M5 tags at all — M5 injection is
	// disabled for the WHOLE header, not per-contig. This matters for a partial
	// reference (e.g. a ~2580-contig BAM against a chr20-only -T FASTA): chr20 is
	// resolvable and would otherwise get an M5 the other 2579 lack, diverging from
	// upstream, which leaves them all bare. Pre-existing M5 tags in the input
	// header are left intact either way (upstream never overwrites them).
	rw.disableM5 = rw.anyUnresolvableSQ()

	// Compute every absent @SQ M5 up front, in parallel where possible; the
	// per-contig hash is the writer's dominant header-write cost. augmentSQLine
	// then reads the result by name (falling back to a serial hash for any
	// contig the precompute did not cover), so the emitted bytes are unchanged.
	// Skipped entirely when M5 is disabled — no contig gets a computed M5.
	if !rw.disableM5 {
		rw.precomputedM5 = rw.precomputeContigM5s()
	}
	cp := *rw.header
	cp.Lines = make([]sam.HeaderLine, len(rw.header.Lines))
	for i, line := range rw.header.Lines {
		if line.Tag != "SQ" {
			cp.Lines[i] = line
			continue
		}
		cp.Lines[i] = rw.augmentSQLine(line)
	}
	return &cp
}

// anyUnresolvableSQ reports whether any @SQ line WITHOUT a pre-existing M5 has a
// contig that cannot be resolved in the reference (absent, or Length <= 0). It
// is the all-or-nothing trigger that mirrors upstream htslib: one unresolvable
// M5-less @SQ disables computed M5 injection for the entire header. An @SQ that
// already carries its own M5 never triggers the fallback (upstream keeps it and
// does not need to hash it). With no reference provider configured, no contig is
// resolvable, so any M5-less @SQ trivially triggers the disable — which is
// correct: without a reference the writer computes no M5 anyway.
func (rw *RecordWriter) anyUnresolvableSQ() bool {
	for _, line := range rw.header.Lines {
		if line.Tag != "SQ" {
			continue
		}
		var name string
		haveM5 := false
		for _, f := range line.Fields {
			switch f.Tag {
			case "SN":
				name = f.Value
			case "M5":
				haveM5 = true
			}
		}
		if haveM5 {
			continue // pre-existing M5: never needs computing, never triggers.
		}
		resolvable := rw.refProvider != nil && name != "" && rw.refProvider.Length(name) > 0
		if !resolvable {
			return true
		}
	}
	return false
}

// augmentSQLine returns a copy of an @SQ header line with the M5 and UR tags
// appended when they are absent and computable. The original line (and its
// Fields slice) is never mutated, so the caller's header is untouched.
func (rw *RecordWriter) augmentSQLine(line sam.HeaderLine) sam.HeaderLine {
	var name string
	haveM5, haveUR := false, false
	for _, f := range line.Fields {
		switch f.Tag {
		case "SN":
			name = f.Value
		case "M5":
			haveM5 = true
		case "UR":
			haveUR = true
		}
	}

	// M5 is computed only for a contig whose bases were loaded from the
	// reference. UR (the -T path) is then injected onto every @SQ that ENDS UP
	// WITH an M5 — whether one was just computed or was already present in the
	// input header. An @SQ with NO M5 (a contig absent from the reference and
	// carrying no pre-existing M5, e.g. a name-mismatched chr-prefixed contig)
	// gets NEITHER tag and is encoded reference-free, exactly as upstream htslib
	// does (verified against samtools on a name-mismatched GIAB BAM).
	// Load the contig's bases (lazily, when a provider is configured) just to
	// hash; refBasesFor caches at most one contig, so M5'ing the whole header
	// streams the reference rather than holding it. A load error is treated as
	// "contig absent" — upstream likewise cannot hash a reference it could not
	// load, and the reference-based encode path surfaces a genuine I/O failure.
	contigInRef := rw.refProvider != nil && rw.refProvider.Length(name) > 0
	var m5 string
	if !haveM5 && contigInRef && !rw.disableM5 {
		// Prefer the value precomputed (in parallel) by precomputeContigM5s;
		// fall back to a serial hash for any name the precompute did not cover
		// (e.g. it ran serial-from-augment, or a name was somehow missed) so no
		// contig silently loses its M5. The two paths produce identical bytes.
		if precomputed, ok := rw.precomputedM5[name]; ok {
			m5 = precomputed
		} else {
			m5 = rw.contigMD5(name)
		}
	}
	addUR := !haveUR && rw.referencePath != "" && (haveM5 || m5 != "")

	if m5 == "" && !addUR {
		return line // nothing to add; leave the line exactly as it was.
	}

	out := line
	out.Fields = make([]sam.HeaderField, len(line.Fields), len(line.Fields)+2)
	copy(out.Fields, line.Fields)
	if m5 != "" {
		out.Fields = append(out.Fields, sam.HeaderField{Tag: "M5", Value: m5})
	}
	if addUR {
		out.Fields = append(out.Fields, sam.HeaderField{Tag: "UR", Value: rw.referencePath})
	}
	return out
}

// contigMD5 returns the lower-case hex MD5 of a contig's reference bases for the
// @SQ M5 tag, hashed in fixed-size chunks so the whole contig is never held
// resident — M5'ing a human chromosome up front would otherwise spike peak RSS
// by ~one chromosome, which dominated the writer's memory. The MD5 of the
// in-order chunk concatenation equals the MD5 of the whole contig. It returns ""
// when a chunk cannot be fetched, matching upstream's "cannot hash a reference
// it could not load" (no M5 emitted).
func (rw *RecordWriter) contigMD5(name string) string {
	n := rw.refProvider.Length(name)
	if n <= 0 {
		return ""
	}
	return hashContigBases(func(off, end int64) ([]byte, error) {
		return rw.refProvider.Fetch(name, off, end)
	}, n)
}

// hashContigBases returns the lower-case hex MD5 of a contig's n bases, fetched
// in 1 MiB chunks via fetch. It returns "" if any chunk cannot be fetched
// (matching upstream's "cannot hash a reference it could not load"). The MD5 of
// the in-order chunk concatenation equals the MD5 of the whole contig, so the
// chunk size, the handle fetch reads through, and the goroutine that drives it
// never change the result — this is what lets contigMD5 (shared provider) and
// precomputeContigM5s (per-worker independent handles) agree byte-for-byte. The
// 1 MiB chunk bounds the transient hash buffer so peak RSS stays bounded even
// across N parallel workers (N × ~1 MiB).
func hashContigBases(fetch func(off, end int64) ([]byte, error), n int64) string {
	if n <= 0 {
		return ""
	}
	const chunk = 1 << 20 // 1 MiB: bounds the transient hash buffer.
	h := md5.New()
	for off := int64(0); off < n; off += chunk {
		end := off + chunk
		if end > n {
			end = n
		}
		bases, err := fetch(off, end)
		if err != nil || bases == nil {
			return ""
		}
		h.Write(bases)
	}
	return hex.EncodeToString(h.Sum(nil))
}

// precomputeContigM5s computes, before the @SQ loop runs, the M5 tag for every
// reference-present @SQ contig that lacks one — exactly the set augmentSQLine
// would otherwise hash serially. For a -T FASTA/BGZF reference (referencePath
// set) it fans the per-contig hashing across a worker pool: the writer's
// dominant header-write cost on a whole-genome reference is hashing all the
// contigs (~3.1 Gbp for GIAB), and the contigs are independent. Each worker
// opens its OWN fasta.RandomAccess handle (its own file descriptor and its own
// seek/inflate state), so no seek-stateful handle is shared across goroutines.
// Length lookups go through the shared provider (a pure index read, safe for
// concurrent use). The result is byte-identical to the serial path: hashing
// order, handle, and goroutine never change a contig's MD5. A "" entry records
// an unfetchable contig (no M5 emitted, exactly as today). It returns nil when
// no contig needs an M5.
func (rw *RecordWriter) precomputeContigM5s() map[string]string {
	if rw.refProvider == nil {
		return nil
	}
	// Collect, in header order, the @SQ names that LACK an M5 and are present
	// in the reference — the precise set augmentSQLine hashes today.
	var names []string
	seen := make(map[string]struct{})
	for _, line := range rw.header.Lines {
		if line.Tag != "SQ" {
			continue
		}
		var name string
		haveM5 := false
		for _, f := range line.Fields {
			switch f.Tag {
			case "SN":
				name = f.Value
			case "M5":
				haveM5 = true
			}
		}
		if haveM5 || name == "" {
			continue
		}
		if _, dup := seen[name]; dup {
			continue
		}
		if rw.refProvider.Length(name) <= 0 {
			continue
		}
		seen[name] = struct{}{}
		names = append(names, name)
	}
	if len(names) == 0 {
		return nil
	}

	// Without a reference path we cannot open independent file handles (the
	// provider may be an in-memory map or a custom provider whose only handle
	// is the shared, seek-stateful one), so hash serially through the shared
	// provider — unchanged behaviour. The map provider is already fast.
	if rw.referencePath == "" {
		out := make(map[string]string, len(names))
		for _, name := range names {
			out[name] = rw.contigMD5(name)
		}
		return out
	}

	// referencePath is set: parallelise across worker-private handles.
	workers := rw.encodeThreads
	if workers < 1 {
		workers = 1
	}
	if workers > len(names) {
		workers = len(names)
	}

	jobs := make(chan string, len(names))
	for _, name := range names {
		jobs <- name
	}
	close(jobs)

	out := make(map[string]string, len(names))
	var mu sync.Mutex
	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			// Each worker opens its OWN handle: an independent file
			// descriptor with its own .fai/.gzi and inflate state, so no
			// seek-stateful state is shared between goroutines. If the open
			// fails, this worker degrades to the shared serial provider
			// rather than crashing; the bytes are identical either way.
			ra, err := fasta.OpenRandomAccess(rw.referencePath)
			if err != nil {
				ra = nil
			} else {
				defer ra.Close()
			}
			for name := range jobs {
				var m5 string
				if ra != nil {
					n := rw.refProvider.Length(name)
					nm := name
					m5 = hashContigBases(func(off, end int64) ([]byte, error) {
						return ra.Fetch(nm, off, end)
					}, n)
				} else {
					// Fallback: the shared provider is seek-stateful, so
					// serialise these reads behind the mutex.
					mu.Lock()
					m5 = rw.contigMD5(name)
					mu.Unlock()
				}
				mu.Lock()
				out[name] = m5
				mu.Unlock()
			}
		}()
	}
	wg.Wait()
	return out
}

// resolveReferenceProvider chooses the writer's reference source: an explicit
// ReferenceProvider when given (the lazy, per-slice production path), otherwise
// the resident Reference map wrapped so the single windowing encode path serves
// both. Nil when neither is set (reference-free).
func resolveReferenceProvider(opts WriterOptions) ReferenceProvider {
	if opts.ReferenceProvider != nil {
		return opts.ReferenceProvider
	}
	if opts.Reference != nil {
		return mapReferenceProvider(opts.Reference)
	}
	return nil
}

// bufContigName is the contig name shared by the records currently buffered,
// or "" when they are unmapped (bufRefID < 0). Every buffered slice is
// single-reference (Write flushes on a reference change), so one name covers
// the whole container.
func (rw *RecordWriter) bufContigName() string {
	if rw.bufRefID < 0 || int(rw.bufRefID) >= len(rw.header.Refs) {
		return ""
	}
	return rw.header.Refs[rw.bufRefID].Name
}

// containerWindow resolves the reference window for the container about to be
// flushed: only the coordinate span the buffered records actually cover is
// fetched, so the writer never holds more than one slice's worth of reference
// resident (this is what bounds peak RSS, matching upstream's per-slice
// reference streaming). It returns the window bases, the window's 0-based start
// position, and whether the contig exists in the reference at all (true even
// when the window is empty, e.g. an all-no-SEQ slice, so the container is still
// flagged reference-using). The returned slice is owned by the encode job and
// read unsynchronised by workers.
func (rw *RecordWriter) containerWindow() (window []byte, winStart int32, hasRef bool, err error) {
	name := rw.bufContigName()
	if name == "" || rw.refProvider == nil {
		return nil, 0, false, nil
	}
	contigLen := rw.refProvider.Length(name)
	if contigLen <= 0 {
		return nil, 0, false, nil // contig absent → reference-free
	}
	start, span := sliceSpan(rw.buf)
	if span <= 0 {
		return nil, 0, true, nil // mapped contig but no positioned bases (all no-SEQ)
	}
	// sliceSpan returns a 1-based start and the span; the window is the
	// half-open 0-based range [start-1, start-1+span), clamped to the contig.
	lo := int64(start) - 1
	if lo < 0 {
		lo = 0
	}
	hi := int64(start) - 1 + int64(span)
	if hi > contigLen {
		hi = contigLen
	}
	bases, ferr := rw.refProvider.Fetch(name, lo, hi)
	if ferr != nil {
		return nil, 0, false, ferr
	}
	return bases, int32(lo), true, nil
}

// resolveEncodeThreads turns the WriterOptions value into a concrete worker
// count: a positive value is honoured; 0 auto-sizes to the CPU count, capped at
// 8 to bound the number of in-flight containers (and their memory).
func resolveEncodeThreads(n int) int {
	if n > 0 {
		return n
	}
	c := runtime.NumCPU()
	if c > 8 {
		c = 8
	}
	if c < 1 {
		c = 1
	}
	return c
}

// flushContainer encodes every buffered record into one data container (a
// compression-header block plus a single slice) and writes it to the
// underlying writer. With encodeThreads == 1 it does so synchronously; with
// more it hands the container off to the worker pool (see startPipeline),
// which encodes containers concurrently and writes them in submission order.
// Either way the buffer ownership is released and the record counter advanced.
func (rw *RecordWriter) flushContainer() error {
	if rw.err != nil {
		return rw.err
	}
	if len(rw.buf) == 0 {
		return nil
	}
	window, winStart, hasRef, err := rw.containerWindow()
	if err != nil {
		rw.err = err
		return err
	}
	if rw.encodeThreads <= 1 {
		container, err := encodeContainer(rw.version, rw.binning, rw.buf, rw.refIndex, window, winStart, hasRef, rw.recordCounter)
		if err != nil {
			rw.err = err
			return err
		}
		if _, err := rw.w.Write(container); err != nil {
			rw.err = fmt.Errorf("cram: writing container: %w", err)
			return rw.err
		}
		rw.recordCounter += int64(len(rw.buf))
		rw.buf = rw.buf[:0]
		rw.bufRefID = refIDUnset
		return nil
	}

	if !rw.pipeStarted {
		rw.startPipeline()
	}
	// Hand the buffered records off to the pool and take a fresh buffer; the
	// records must remain stable until encoded, which holds under the same
	// caller contract as the synchronous path (records are already retained
	// across a whole slice). Order is preserved by pushing the result channel
	// onto rw.ordered before the writer can reach it.
	job := &encodeJob{records: rw.buf, refWindow: window, refWindowStart: winStart, hasRef: hasRef, counter: rw.recordCounter, out: make(chan encodeResult, 1)}
	rw.ordered <- job.out
	rw.jobs <- job
	rw.recordCounter += int64(len(rw.buf))
	rw.buf = make([]*sam.Record, 0, rw.recordsPerSlice)
	rw.bufRefID = refIDUnset
	return rw.err
}

// startPipeline spins up the encode worker pool and the single ordered-writer
// goroutine. Workers pull jobs and encode containers concurrently; the writer
// consumes result channels in submission order and writes each container,
// latching the first error into pipeErr.
func (rw *RecordWriter) startPipeline() {
	rw.pipeStarted = true
	rw.jobs = make(chan *encodeJob, rw.encodeThreads)
	rw.ordered = make(chan chan encodeResult, 2*rw.encodeThreads)
	rw.pipeDone = make(chan struct{})

	for i := 0; i < rw.encodeThreads; i++ {
		rw.pipeWG.Add(1)
		go func() {
			defer rw.pipeWG.Done()
			for job := range rw.jobs {
				data, err := encodeContainer(rw.version, rw.binning, job.records, rw.refIndex, job.refWindow, job.refWindowStart, job.hasRef, job.counter)
				job.out <- encodeResult{data: data, err: err}
			}
		}()
	}
	go func() {
		defer close(rw.pipeDone)
		for out := range rw.ordered {
			res := <-out
			if rw.pipeErr != nil {
				continue // drain remaining results after an error
			}
			if res.err != nil {
				rw.pipeErr = res.err
				continue
			}
			if _, err := rw.w.Write(res.data); err != nil {
				rw.pipeErr = fmt.Errorf("cram: writing container: %w", err)
			}
		}
	}()
}

// finishPipeline drains the encode pipeline: it stops the workers, waits for
// the ordered writer to flush every remaining container, and returns the first
// error any of them latched. It is a no-op when the pipeline never started.
func (rw *RecordWriter) finishPipeline() error {
	if !rw.pipeStarted {
		return nil
	}
	close(rw.jobs)    // workers exit once the in-flight jobs drain
	rw.pipeWG.Wait()  // all results have now been sent
	close(rw.ordered) // writer exits once it has consumed them in order
	<-rw.pipeDone
	rw.pipeStarted = false
	return rw.pipeErr
}

// encodeBlock assembles a complete on-disk CRAM block from a content
// type, content id and uncompressed payload. The payload is compressed
// with whichever method chooseBlockCompression picks for the writer's
// version — never larger than raw — and the trailing IEEE CRC32 over the
// whole block is appended. The content id and the two size fields are
// framed as ITF-8 for CRAM v2/v3 and as uint7 varints for v4, the inverse
// of readBlock.
func encodeBlock(version Version, ct BlockContentType, contentID int32, payload []byte) []byte {
	method, stored := chooseBlockCompression(version, payload)
	iw := newIntWriter(version)
	var b []byte
	b = append(b, byte(method), byte(ct))
	b = iw.u32(b, contentID)
	b = iw.u32(b, int32(len(stored)))
	b = iw.u32(b, int32(len(payload)))
	b = append(b, stored...)
	crc := crc32.Checksum(b, crc32.IEEETable)
	var crcBuf [4]byte
	binary.LittleEndian.PutUint32(crcBuf[:], crc)
	return append(b, crcBuf[:]...)
}

// chooseBlockCompression picks the block compression method for a
// payload and returns it together with the bytes to store. The candidate
// set depends on the CRAM version: every version considers raw (method
// 0) and gzip (method 1); v3.1 additionally considers rANS 4x16 (method
// 5), its distinguishing codec. The smallest candidate wins and raw is
// always in the running, so the stored payload is never larger than the
// input — the block stays decodable even if a codec misbehaves.
//
// bzip2 is deliberately NOT in the default candidate set. The in-tree
// pure-Go bzip2 encoder is ~10x slower than gzip and, brute-forced on
// every block, dominated encode time — a 48 MB BAM -> CRAM took 97 s with
// it versus 10 s without — while shrinking the output by only ~0.5%.
// Upstream htslib likewise omits bzip2 from its default CRAM profile
// (using it only in the opt-in "archive"/"small" profiles). The decoder
// still reads bzip2 blocks (block.go's Decompress) for files that use
// them, and codec.Bzip2Encode stays in-tree for an explicit archive
// profile; it is simply never auto-selected by the default writer.
func chooseBlockCompression(version Version, payload []byte) (CompressionMethod, []byte) {
	method := CompRaw
	stored := payload
	// Below the threshold the per-codec framing overhead outweighs any
	// saving, so a tiny block is always stored raw.
	if len(payload) < gzipBlockThreshold {
		return method, stored
	}
	if gz := gzipCompress(payload); len(gz) < len(stored) {
		method, stored = CompGzip, gz
	}
	// rANS: the entropy coder upstream leans on for CRAM. v3.0 uses rANS 4x8,
	// v3.1 the newer rANS 4x16; offer the version-appropriate codec at both
	// orders and keep whichever is smallest. order 0 (frequency model) is
	// trusted — it is corpus-validated and round-trips every input including the
	// empty/degenerate blocks. order 1 (context model) usually wins big on
	// quality scores and base calls, but it does NOT round-trip every degenerate
	// input, so each order-1 candidate is decoded back and only kept if it
	// reproduces the payload exactly.
	tryRANS := func(m CompressionMethod, enc func([]byte, int) ([]byte, error), dec func([]byte) ([]byte, error)) {
		if r, err := enc(payload, 0); err == nil && len(r) < len(stored) {
			method, stored = m, r
		}
		if r, err := enc(payload, 1); err == nil && len(r) < len(stored) {
			if back, derr := dec(r); derr == nil && bytes.Equal(back, payload) {
				method, stored = m, r
			}
		}
	}
	if version == VersionV31 {
		tryRANS(CompRANS4x16, codec.RANS4x16Encode, codec.RANS4x16Decode)
	} else {
		tryRANS(CompRANS4x8, codec.RANS4x8Encode, codec.RANS4x8Decode)
	}
	return method, stored
}

// cramGzipLevel is the klauspost deflate level used for CRAM gzip blocks.
// klauspost's level scale differs from the stdlib's: its level 7 reproduces the
// ratio of stdlib gzip's default (level 6) — which the writer targeted before —
// whereas klauspost level 6 trades a large chunk of ratio for speed. Level 7
// keeps the CRAM size on par with upstream while decoding/encoding faster.
const cramGzipLevel = 7

// gzipWriterPool recycles gzip writers across blocks. A CRAM container holds
// many series blocks and chooseBlockCompression gzips every one, so without
// pooling the deflate state (the bulk of the per-block allocation) would be
// re-allocated thousands of times per file.
var gzipWriterPool = sync.Pool{
	New: func() any {
		zw, _ := kgzip.NewWriterLevel(io.Discard, cramGzipLevel)
		return zw
	},
}

// gzipCompress returns the gzip (RFC 1952) compression of in. It is the
// encode-side counterpart of block.go's gunzip. It uses the klauspost/compress
// gzip backend — the same faster, pure-Go DEFLATE implementation BGZF I/O uses
// — emitting a standard gzip stream that the stdlib reader and upstream htslib
// both decode, and recycles the writer through gzipWriterPool.
func gzipCompress(in []byte) []byte {
	buf := bytes.NewBuffer(make([]byte, 0, len(in)/2+64))
	zw := gzipWriterPool.Get().(*kgzip.Writer)
	zw.Reset(buf)
	if _, err := zw.Write(in); err != nil {
		gzipWriterPool.Put(zw)
		return in
	}
	if err := zw.Close(); err != nil {
		gzipWriterPool.Put(zw)
		return in
	}
	gzipWriterPool.Put(zw)
	return buf.Bytes()
}

// containerFields collects the variable-length fields of a CRAM
// container header so containerHeaderBytes can serialise them.
type containerFields struct {
	length        int32
	refSeqID      int32
	startPos      int32
	alignmentSpan int32
	numRecords    int32
	recordCounter int64
	numBases      int64
	numBlocks     int32
	landmarks     []int32
}

// containerHeaderBytes serialises a CRAM container header and the trailing
// IEEE CRC32 over it. It is the writer-side inverse of
// parseContainerHeader.
//
// For CRAM v2/v3 the length is a fixed 4-byte little-endian int32 (matching
// htslib's on-disk reality, which back-patches it) and the remaining fields
// are ITF-8 / LTF-8. For CRAM v4.0 every field — the length included — is a
// uint7 varint: ref_seq_id is signed (zig-zag), the alignment start and
// span are 64-bit varints, and the counts and landmarks are unsigned
// varints (cram_io.c cram_write_container, major>=4 branch).
func containerHeaderBytes(version Version, f containerFields) []byte {
	iw := newIntWriter(version)
	var b []byte
	if version.usesUint7() {
		// v4: the length is a uint7 varint, not a fixed 4-byte int32.
		b = iw.u32(b, f.length)
		// ref_seq_id is signed; the alignment start and span widen to 64-bit.
		b = iw.s32(b, f.refSeqID)
		b = iw.u64(b, int64(f.startPos))
		b = iw.u64(b, int64(f.alignmentSpan))
		b = iw.u32(b, f.numRecords)
		b = iw.u64(b, f.recordCounter)
		b = iw.u64(b, f.numBases)
		b = iw.u32(b, f.numBlocks)
		b = iw.u32(b, int32(len(f.landmarks)))
		for _, lm := range f.landmarks {
			b = iw.u32(b, lm)
		}
	} else {
		var lenBuf [4]byte
		binary.LittleEndian.PutUint32(lenBuf[:], uint32(f.length))
		b = append(b, lenBuf[:]...)
		b = appendITF8(b, f.refSeqID)
		b = appendITF8(b, f.startPos)
		b = appendITF8(b, f.alignmentSpan)
		b = appendITF8(b, f.numRecords)
		b = appendLTF8(b, f.recordCounter)
		b = appendLTF8(b, f.numBases)
		b = appendITF8(b, f.numBlocks)
		b = appendITF8(b, int32(len(f.landmarks)))
		for _, lm := range f.landmarks {
			b = appendITF8(b, lm)
		}
	}
	crc := crc32.Checksum(b, crc32.IEEETable)
	var crcBuf [4]byte
	binary.LittleEndian.PutUint32(crcBuf[:], crc)
	return append(b, crcBuf[:]...)
}
