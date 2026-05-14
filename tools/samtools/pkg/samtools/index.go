package samtools

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
)

// IndexOptions configures the index builder. The zero value selects the
// .bai default; CSI is accepted via SelectCSI but rejected with a clear
// error in v1 so callers can wire the flag now and a later slice can flip
// the bit without touching the CLI.
type IndexOptions struct {
	// SelectCSI requests a .csi index (>512Mb chromosomes). Not yet
	// implemented; setting it surfaces ErrCSIUnsupported.
	SelectCSI bool
	// CSIMinShift is honoured only when SelectCSI is true; accepted for
	// flag-compatibility with upstream samtools.
	CSIMinShift int
	// Threads is accepted but ignored — the v1 index pipeline is
	// single-threaded.
	Threads int
}

// ErrCSIUnsupported is returned when an index build asks for CSI output.
// CSI is deferred to a later slice; v1 only emits BAI.
var ErrCSIUnsupported = errors.New("samtools index: CSI output (-c/--csi) is not yet implemented; v1 emits BAI only")

// Index reads a coordinate-sorted BAM stream from in, builds a BAI index in
// memory, and writes it to out. The input must be a coordinate-sorted BAM
// (the @HD SO:coordinate header is not strictly required, but the records
// must arrive in increasing (refID, pos) order — for our use case the
// `samtools sort` output guarantees this).
func Index(in io.Reader, out io.Writer, opts IndexOptions) error {
	if opts.SelectCSI {
		return ErrCSIUnsupported
	}
	br, err := sam.NewBAMReader(in)
	if err != nil {
		return err
	}
	defer br.Close()

	hdr := br.Header()
	idx, err := BuildBAI(br, len(hdr.Refs))
	if err != nil {
		return err
	}
	return WriteBAI(out, idx)
}

// BuildBAI streams every record from a BAMReader and returns the assembled
// BAI index. It does not validate sort order — callers must pass a
// coordinate-sorted BAM (where coordinates means (refID, 0-based pos)).
func BuildBAI(br *sam.BAMReader, numRefs int) (*BAIIndex, error) {
	bld := NewBAIBuilder(numRefs)
	for {
		vBeg := br.VirtualOffset()
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		vEnd := br.VirtualOffset()

		// Resolve refID from RName via header order.
		refID := -1
		if rec.RName != "" && rec.RName != "*" {
			refID = br.Header().RefIndex(rec.RName)
			if refID < 0 {
				return nil, fmt.Errorf("samtools index: record references unknown @SQ %q", rec.RName)
			}
		}
		mapped := !rec.IsUnmapped()
		// Unmapped records that nevertheless carry a refID + pos still go
		// into the regular bin/linear index per the SAM spec: htslib treats
		// them as "placed but unmapped". Records with refID == -1 are the
		// truly unplaced ones that bump n_no_coor.
		beg := int(rec.Pos) - 1
		if beg < 0 {
			beg = 0
		}
		end := beg + rec.Cigar.ReferenceLength()
		if err := bld.AddRecord(refID, beg, end, vBeg, vEnd, mapped); err != nil {
			return nil, err
		}
	}
	return bld.Finish(), nil
}

// IndexFile reads the BAM at inPath, builds a BAI index, and writes it to
// outPath. If outPath is empty, the index is written to <inPath>.bai.
func IndexFile(inPath, outPath string, opts IndexOptions) error {
	if opts.SelectCSI {
		return ErrCSIUnsupported
	}
	if outPath == "" {
		outPath = inPath + ".bai"
	}
	in, err := os.Open(inPath)
	if err != nil {
		return err
	}
	defer in.Close()

	// Materialise to a sibling tmp file then rename so a half-written index
	// never replaces a real one.
	tmp, err := os.CreateTemp(filepath.Dir(outPath), ".bai.tmp.")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	cleanup := func() { _ = os.Remove(tmpName) }
	if err := Index(in, tmp, opts); err != nil {
		_ = tmp.Close()
		cleanup()
		return err
	}
	if err := tmp.Close(); err != nil {
		cleanup()
		return err
	}
	return os.Rename(tmpName, outPath)
}
