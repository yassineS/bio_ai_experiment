// BAM binary output for `bedtools intersect` when the query (-a) is a BAM (or
// CRAM) file and -bed is NOT given. Upstream writes the intersecting ALIGNMENTS
// back out as BAM by default — the surviving original alignment records, in
// input order, under the alignment-level flags that make sense for BAM output
// (-u / -v / -wa / default). This path preserves the original SAM/BAM records
// verbatim (re-encoded by pkg/htsgo/sam's BAMWriter) and emits the original
// header, rather than the BED12 projection used by the text-output path.
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

// IsBAMOrCRAMInput reports whether r is a BAM or CRAM alignment stream, returning
// a replacement reader that re-yields the bytes consumed while sniffing (so the
// returned reader must be used in place of r). It is the query-format probe the
// CLI uses to decide whether BAM-A output gating applies: upstream determines
// the output type (BAM vs BED) from the query file's type.
//
// The probe deliberately does NOT classify "any BGZF stream" as BAM: a
// BGZF-compressed BED/GFF/VCF (e.g. a `.bed.gz` piped to `-a -`, which iohelper
// does not transparently decompress for stdin) is still text and must produce
// text output. So a gzip/BGZF stream is only BAM when its decompressed prefix is
// the "BAM\1" magic, mirroring readInRecords' own detection. CRAM is recognised
// by its four-byte `CRAM` magic; a raw (already-decompressed) BAM body by the
// "BAM\1" magic. Everything else (SAM text, plain BED/GFF/VCF) returns false.
func IsBAMOrCRAMInput(r io.Reader) (bool, io.Reader, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	magic, _ := br.Peek(4)
	switch {
	case len(magic) >= 4 && string(magic) == "CRAM":
		return true, br, nil
	case len(magic) >= 4 && string(magic) == "BAM\x01":
		return true, br, nil
	case len(magic) >= 4 && magic[0] == 0x1f && magic[1] == 0x8b:
		// gzip/BGZF: BAM is BGZF-wrapped, so a BAM stream's decompressed prefix is
		// "BAM\1"; a plain gzipped/BGZF text file is not. Buffer the whole stream
		// (so the returned reader can replay it) and probe the decompressed prefix.
		buf, err := io.ReadAll(br)
		if err != nil {
			return false, bytes.NewReader(buf), err
		}
		return isGzippedBAM(buf), bytes.NewReader(buf), nil
	default:
		return false, br, nil
	}
}

// IntersectBAMOutput intersects a BAM/CRAM query (-a) against one or more B
// files and writes the surviving alignment records back out as BAM, matching
// upstream `bedtools intersect`'s default BAM-output behaviour. readerA must be
// a BAM or CRAM stream (the caller has already classified it); readersB may be
// any supported B input (BED/GFF/VCF/BAM/CRAM). The original BAM header and the
// surviving alignment records are written verbatim, in input order.
//
// Only the alignment-level output modes are honoured here (selected by opts):
// the default and -wa/-u write each A alignment that has at least one overlap
// once; -v writes each A alignment that has no overlap. -wb/-loj are ignored
// (the caller has already emitted the upstream warning); -c/-C/-wo/-wao are
// rejected by the caller before reaching this path. It returns the number of
// alignment records written.
func IntersectBAMOutput(readerA io.Reader, readersB []io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	bRecords, _, _, err := readAllB(readersB)
	if err != nil {
		return 0, err
	}
	finder := newFinder(bRecords, opts)

	aReader, err := alnio.NewReaderWithReference(readerA, "")
	if err != nil {
		return 0, fmt.Errorf("error opening BAM/CRAM input: %w", err)
	}

	hdr := aReader.Header()
	bw := sam.NewBAMWriter(writer)
	if err := bw.WriteHeader(hdr); err != nil {
		return 0, fmt.Errorf("error writing BAM header: %w", err)
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
		if err := bw.Write(rec); err != nil {
			return count, fmt.Errorf("error writing BAM record: %w", err)
		}
		count++
	}
	if err := bw.Close(); err != nil {
		return count, fmt.Errorf("error finalising BAM output: %w", err)
	}
	return count, nil
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
