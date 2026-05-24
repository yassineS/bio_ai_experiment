package samtools

import (
	"io"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// IndexOptions configures the index builder. The zero value selects the
// .bai default; setting SelectCSI requests a coordinate-sorted index
// (.csi), which is required for reference sequences longer than the BAI
// 2^29 bp ceiling.
type IndexOptions struct {
	// SelectCSI requests a .csi index (>512Mb chromosomes).
	SelectCSI bool
	// CSIMinShift is the CSI bin-hierarchy min_shift; it is honoured only
	// when SelectCSI is true. Zero (or any non-positive value) selects the
	// htslib default of 14.
	CSIMinShift int
	// Threads is accepted for upstream samtools index -@/--threads parity.
	// The BAI/CSI build is intrinsically sequential — one entry per ref
	// bin per BGZF virtual offset, walked in coordinate order — so
	// pipelining the bgzf decode side is the only available win and is
	// not yet implemented. This field is therefore currently a no-op;
	// the CLI accepts the flag without silently changing behaviour.
	Threads int
}

// Index reads a coordinate-sorted BAM stream from in, builds an index in
// memory, and writes it to out. The input must be a coordinate-sorted BAM
// (the @HD SO:coordinate header is not strictly required, but the records
// must arrive in increasing (refID, pos) order — for our use case the
// `samtools sort` output guarantees this). When opts.SelectCSI is set a
// BGZF-compressed .csi index is written; otherwise a plain .bai index.
func Index(in io.Reader, out io.Writer, opts IndexOptions) error {
	br, err := sam.NewBAMReader(in)
	if err != nil {
		return err
	}
	defer br.Close()

	hdr := br.Header()
	if opts.SelectCSI {
		idx, err := bam.BuildCSI(br, len(hdr.Refs), int32(opts.CSIMinShift), bam.DefaultCSIDepth)
		if err != nil {
			return err
		}
		return bam.WriteCSI(out, idx)
	}
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
// input format and options: a BAM file gets a BAI written to
// <inPath>.bai by default, or a .csi written to <inPath>.csi when
// opts.SelectCSI is set, while a CRAM file gets a CRAI written to
// <inPath>.crai. When outPath is non-empty it overrides the default
// destination.
func IndexFile(inPath, outPath string, opts IndexOptions) error {
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
		if opts.SelectCSI {
			outPath = inPath + ".csi"
		} else {
			outPath = inPath + ".bai"
		}
	}
	in, err := os.Open(inPath)
	if err != nil {
		return err
	}
	defer in.Close()

	// Materialise to a sibling tmp file then rename so a half-written index
	// never replaces a real one.
	tmp, err := os.CreateTemp(filepath.Dir(outPath), ".idx.tmp.")
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
