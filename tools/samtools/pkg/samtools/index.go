package samtools

import (
	"errors"
	"io"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
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

// IndexFile reads the alignment file at inPath, builds the appropriate
// index, and writes it to outPath. The index kind is chosen from the
// input format: a BAM file gets a BAI written to <inPath>.bai, while a
// CRAM file gets a CRAI written to <inPath>.crai. When outPath is
// non-empty it overrides the default destination.
func IndexFile(inPath, outPath string, opts IndexOptions) error {
	// CSI is a BAM-only index kind and is not yet implemented; reject it
	// up front, before the input file is even opened, so the deferral is
	// surfaced regardless of the input format.
	if opts.SelectCSI {
		return ErrCSIUnsupported
	}
	isCRAM, err := inputIsCRAM(inPath)
	if err != nil {
		return err
	}
	if isCRAM {
		if outPath == "" {
			outPath = inPath + ".crai"
		}
		return indexCRAMFile(inPath, outPath)
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

// inputIsCRAM reports whether the file at inPath is a CRAM stream. It
// sniffs the leading bytes through iohelper's format detector; a SAM or
// BAM file reports false.
func inputIsCRAM(inPath string) (bool, error) {
	f, err := os.Open(inPath)
	if err != nil {
		return false, err
	}
	defer f.Close()
	format, _, err := iohelper.DetectFormat(f)
	if err != nil {
		return false, err
	}
	return format == iohelper.FormatCRAM, nil
}

// indexCRAMFile builds a .crai index for the CRAM at inPath and writes it
// to outPath. The index is materialised to a sibling temp file and then
// renamed so a half-written .crai never replaces a real one.
func indexCRAMFile(inPath, outPath string) error {
	entries, err := cram.BuildCRAIFile(inPath)
	if err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(outPath), ".crai.tmp.")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	cleanup := func() { _ = os.Remove(tmpName) }
	if err := cram.WriteCRAI(tmp, entries); err != nil {
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
