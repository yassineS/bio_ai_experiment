package cram

import (
	"bufio"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// RegionReader answers .crai-indexed CRAM region queries by seeking to the
// containers a query overlaps instead of streaming the whole file. It reads
// the file definition and SAM-header container once (from offset 0) to gather
// the per-file context record reconstruction needs — ref names, read groups,
// and the file's CRAM version — then, for each query, seeks the underlying
// io.ReadSeeker to each matching container's byte offset and decodes only
// those containers.
//
// A RegionReader reuses the sequential RecordReader's header-parse and
// container-decode paths verbatim: a CRAM container is self-contained (it
// carries its own compression-header block), so it can be decoded at an
// arbitrary offset given the offset-0 context the RecordReader already
// gathered. Reference-backed decode is honoured exactly as for RecordReader;
// attach a reference with SetReference / SetReferenceFASTA / SetRefCache /
// UseRefCacheFromEnv before the first Query.
//
// A RegionReader is not safe for concurrent use: each Query advances the
// shared io.ReadSeeker.
type RegionReader struct {
	rs  io.ReadSeeker
	rr  *RecordReader
	idx *CRAIIndex
	def FileDefinition
}

// NewRegionReader builds a seek-based region reader over rs using the parsed
// .crai index idx. It reads the CRAM file definition and SAM-header container
// from offset 0, so rs must be positioned anywhere — NewRegionReader seeks to
// the start itself. It returns an error if rs is not a CRAM stream or its
// embedded header cannot be parsed.
func NewRegionReader(rs io.ReadSeeker, idx *CRAIIndex) (*RegionReader, error) {
	if _, err := rs.Seek(0, io.SeekStart); err != nil {
		return nil, wrapf(err, "seeking to the start of the CRAM stream")
	}
	rd, err := NewReader(rs)
	if err != nil {
		return nil, err
	}
	rr := &RecordReader{rd: rd}
	if err := rr.readSAMHeader(); err != nil {
		return nil, err
	}
	return &RegionReader{rs: rs, rr: rr, idx: idx, def: rd.FileDefinition()}, nil
}

// Header returns the SAM header parsed from the CRAM file's first container.
// It is available immediately after NewRegionReader and is the header a
// caller resolves region chrom names against.
func (rg *RegionReader) Header() *sam.Header { return rg.rr.Header() }

// SetReference attaches an external reference source so reference-backed
// mapped reads reconstruct their bases instead of being filled with 'N'. It
// mirrors RecordReader.SetReference and must be called before the first
// Query.
func (rg *RegionReader) SetReference(src ReferenceSource) { rg.rr.SetReference(src) }

// SetReferenceFASTA opens the named FASTA file as the decode reference and
// attaches it. It mirrors RecordReader.SetReferenceFASTA; the FASTA's file
// handle is released by Close.
func (rg *RegionReader) SetReferenceFASTA(path string) error { return rg.rr.SetReferenceFASTA(path) }

// SetRefCache attaches the htslib REF_CACHE directory as a decode reference,
// looked up by each slice header's MD5. It mirrors RecordReader.SetRefCache.
func (rg *RegionReader) SetRefCache(dir string) { rg.rr.SetRefCache(dir) }

// UseRefCacheFromEnv attaches the REF_CACHE directory named by the REF_CACHE
// environment variable, if set, as a decode reference. It reports whether a
// cache was attached and mirrors RecordReader.UseRefCacheFromEnv.
func (rg *RegionReader) UseRefCacheFromEnv() bool { return rg.rr.UseRefCacheFromEnv() }

// UseRefPathFromEnv attaches the network REF_PATH URL-fetch source named by the
// REF_PATH environment variable, if set, mirroring
// RecordReader.UseRefPathFromEnv.
func (rg *RegionReader) UseRefPathFromEnv() bool { return rg.rr.UseRefPathFromEnv() }

// Close releases the underlying reference FASTA handle, if any was attached.
// It does not close the io.ReadSeeker, whose lifetime the caller owns.
func (rg *RegionReader) Close() error {
	if rg.rr.refResolver != nil && rg.rr.refResolver.fasta != nil {
		return rg.rr.refResolver.fasta.Close()
	}
	return nil
}

// Query returns every reconstructed record overlapping the resolved region,
// in file order. It uses the .crai index to find the containers the region
// touches, seeks to each distinct container exactly once, decodes its
// records, and keeps only the records whose reference span overlaps the
// region.
//
// Overlap matches upstream samtools' coordinate-region semantics exactly — the
// same rule the BAM indexed path applies: a record overlaps chr:beg-end iff it
// is placed on the region's reference (its RName matches) and its
// [POS, POS+refLen) footprint intersects the region. A flag-unmapped but
// mate-placed read is included at its recorded POS (as upstream's coordinate
// iterator does); only a truly unplaced read (no reference name) is excluded.
func (rg *RegionReader) Query(reg region.ResolvedRegion) ([]*sam.Record, error) {
	hits := rg.idx.QueryRegion(reg)
	if len(hits) == 0 {
		return nil, nil
	}
	regionName := rg.rr.refNameForRegion(reg.RefID)
	var out []*sam.Record
	var lastOffset int64 = -1
	for _, e := range hits {
		// QueryRegion returns entries sorted by container offset; a single
		// container can hold several overlapping slices, so decode each
		// distinct container exactly once.
		if e.ContainerOffset == lastOffset {
			continue
		}
		lastOffset = e.ContainerOffset
		c, err := rg.readContainerAt(e.ContainerOffset)
		if err != nil {
			return out, err
		}
		var recs []*sam.Record
		if err := rg.rr.decodeContainerInto(c, &recs); err != nil {
			return out, err
		}
		for _, rec := range recs {
			if regionOverlap(rec, reg, regionName) {
				out = append(out, rec)
			}
		}
	}
	return out, nil
}

// readContainerAt seeks rs to off and reads the single container that begins
// there. The container is self-contained, so only the file definition (its
// CRAM version, which fixes the structural widths) is needed beyond the
// stream itself.
func (rg *RegionReader) readContainerAt(off int64) (*Container, error) {
	if _, err := rg.rs.Seek(off, io.SeekStart); err != nil {
		return nil, wrapf(err, "seeking to container at offset %d", off)
	}
	rd := &Reader{def: rg.def, r: bufio.NewReader(rg.rs)}
	c, err := rd.Next()
	if err != nil {
		return nil, wrapf(err, "reading container at offset %d", off)
	}
	return c, nil
}

// refNameForRegion returns the @SQ name for a region's reference id, or ""
// when the id is out of range. The empty string never matches a mapped
// record's RName, so an out-of-range region simply yields no hits.
func (rr *RecordReader) refNameForRegion(refID int) string {
	if refID < 0 || refID >= len(rr.refNames) {
		return ""
	}
	return rr.refNames[refID]
}

// regionOverlap reports whether rec overlaps the resolved region. It mirrors
// upstream samtools' coordinate-region rule exactly — the same rule the BAM
// indexed path's buildRegionFilter applies: the record must be placed on the
// region's reference (its RName matches, regardless of the unmapped flag —
// htslib's coordinate iterator returns a flag-unmapped but mate-placed read
// at its recorded POS), and its [POS, POS+refLen) footprint must intersect
// the region's half-open 0-based [Beg0, End0).
//
// A record whose CIGAR consumes no reference (a pure-insertion/clip CIGAR, or
// the "*" CIGAR a placed-unmapped read carries) is treated as a length-1
// footprint at its POS, matching buildRegionFilter's refLen<=0 -> 1 rule and
// upstream's inclusion of such reads in a region containing their POS.
func regionOverlap(rec *sam.Record, reg region.ResolvedRegion, regionName string) bool {
	if regionName == "" || rec.RName == "" || rec.RName == "*" {
		return false
	}
	if rec.RName != regionName {
		return false
	}
	refLen := rec.Cigar.ReferenceLength()
	if refLen <= 0 {
		refLen = 1
	}
	recBeg0 := int(rec.Pos) - 1
	if recBeg0 < 0 {
		recBeg0 = 0
	}
	recEnd0 := recBeg0 + refLen
	return recBeg0 < reg.End0 && recEnd0 > reg.Beg0
}
