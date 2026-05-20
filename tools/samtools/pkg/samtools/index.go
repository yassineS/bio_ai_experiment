package samtools

import (
	"errors"
	"io"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
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
	return bam.WriteBAI(out, idx)
}

// BuildBAI is preserved as a thin delegation to pkg/htsgo/bam.BuildBAI
// for backwards compatibility with any external caller importing the
// samtools package directly. The real implementation moved to the
// htsgo/bam package in PR-D.
func BuildBAI(br *sam.BAMReader, numRefs int) (*bam.BAIIndex, error) {
	return bam.BuildBAI(br, numRefs)
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
