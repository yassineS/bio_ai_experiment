// Package bedigv implements `bedtools igv`: it reads BED intervals and emits
// an IGV batch-mode script that takes one snapshot per interval.
//
// Upstream reference: reference_code/bedtools/src/bedToIgv/bedToIgv.cpp.
//
// Output shape (matches the upstream binary byte-for-byte for the supported
// option subset):
//
//	snapshotDirectory <path>
//	[load <session>]
//	goto <chrom>:<start-slop>-<end+slop>
//	[sort <sortType>]
//	[collapse]
//	snapshot <filename>.<imageType>
//	... (repeat per interval)
//
// Notes:
//   - Snapshot filename defaults to `<chrom>_<start>_<end>` (with `_<name>`
//     appended when `-name` is set and the BED record's name column is
//     non-empty, and `_slop<N>` appended when `-slop` > 0). Coordinates in
//     the filename are the *original* BED coordinates (slop is applied only
//     to the `goto` locus), mirroring upstream.
//   - `-clps` emits a single `collapse` line per interval, identical to
//     upstream (the user-facing "collapse single-record output into a single
//     goto" reading in the task brief does not match the upstream code, which
//     always emits one goto + one snapshot per record).
//   - `-name` in upstream is a boolean (use BED column 4 as part of the
//     filename) rather than the `{NUM|POS}` value flag described in the task
//     brief. We follow upstream so the byte-for-byte output matches.
package bedigv

import (
	"bufio"
	"fmt"
	"io"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// ImageType is the snapshot extension passed to IGV's `snapshot` command.
// Upstream accepts any string but documents `png`, `eps`, `svg`. We
// validate against that set when configured via the CLI; the library
// accepts any non-empty string.
type ImageType string

// Common IGV snapshot extensions.
const (
	ImagePNG ImageType = "png"
	ImageEPS ImageType = "eps"
	ImageSVG ImageType = "svg"
	ImageJPG ImageType = "jpg"
)

// SortType is the optional BAM-sort directive emitted between each `goto`
// and its `snapshot`. The empty string disables sort emission. Upstream
// accepts the six values below (any other value is rejected).
type SortType string

// Valid SortType values per upstream bedToIgv.cpp.
const (
	SortNone      SortType = ""
	SortBase      SortType = "base"
	SortPosition  SortType = "position"
	SortStrand    SortType = "strand"
	SortQuality   SortType = "quality"
	SortSample    SortType = "sample"
	SortReadGroup SortType = "readGroup"
)

// IsValidSort returns true if s is a recognised upstream sort directive
// (including the empty/no-sort default).
func IsValidSort(s SortType) bool {
	switch s {
	case SortNone, SortBase, SortPosition, SortStrand, SortQuality, SortSample, SortReadGroup:
		return true
	}
	return false
}

// Options configures Run. Zero-valued fields take the upstream defaults
// (no session, no sort, no collapse, no name, slop=0, image=png, path=./).
type Options struct {
	// Path is the directory IGV will write snapshots into. Emitted as
	// `snapshotDirectory <Path>`. Defaults to "./" when empty.
	Path string

	// Session is the path to an existing IGV session file to load before
	// taking snapshots. When empty, no `load` line is emitted.
	Session string

	// Sort selects an optional BAM-sort directive emitted between each
	// interval's `goto` and `snapshot`. SortNone (empty) disables it.
	Sort SortType

	// Collapse causes a `collapse` line to be emitted after each `goto`
	// (and after `sort` if both are set).
	Collapse bool

	// UseNames, when true, appends `_<name>` to each snapshot filename,
	// where `<name>` is the BED record's name column (column 4). If the
	// record's name is empty, Run returns an error mirroring upstream.
	UseNames bool

	// Slop extends each interval by Slop bp on both sides for the `goto`
	// locus. Must be >= 0. Does NOT change the filename's coordinates
	// (matching upstream).
	Slop int

	// ImageType is the snapshot extension. Defaults to ImagePNG when
	// empty.
	ImageType ImageType
}

// Run reads BED records from r and writes the IGV batch script to w. The
// returned count is the number of `snapshot` lines emitted (one per valid
// BED record).
func Run(r io.Reader, w io.Writer, opts Options) (int, error) {
	if opts.Slop < 0 {
		return 0, fmt.Errorf("slop must be >= 0 (got %d)", opts.Slop)
	}
	if opts.Sort != SortNone && !IsValidSort(opts.Sort) {
		return 0, fmt.Errorf("invalid sort type %q", opts.Sort)
	}

	path := opts.Path
	if path == "" {
		path = "./"
	}
	imgType := opts.ImageType
	if imgType == "" {
		imgType = ImagePNG
	}

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	if _, err := fmt.Fprintf(bw, "snapshotDirectory %s\n", path); err != nil {
		return 0, err
	}
	if opts.Session != "" {
		if _, err := fmt.Fprintf(bw, "load %s\n", opts.Session); err != nil {
			return 0, err
		}
	}

	br := bed.NewReader(r)
	n := 0
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return n, err
		}

		startStr := strconv.Itoa(rec.ChromStart)
		endStr := strconv.Itoa(rec.ChromEnd)
		filename := rec.Chrom + "_" + startStr + "_" + endStr

		// `goto` locus uses slop-expanded coordinates. Upstream emits
		// the raw arithmetic result, including negative values when
		// start-slop underflows; we preserve that behaviour.
		locus := rec.Chrom + ":" + strconv.Itoa(rec.ChromStart-opts.Slop) + "-" + strconv.Itoa(rec.ChromEnd+opts.Slop)

		if opts.UseNames {
			if rec.Name == "" {
				return n, fmt.Errorf("requested -name but BED record %d has empty name column", n+1)
			}
			filename = filename + "_" + rec.Name
		}
		if opts.Slop > 0 {
			filename = filename + "_slop" + strconv.Itoa(opts.Slop)
		}

		if _, err := fmt.Fprintf(bw, "goto %s\n", locus); err != nil {
			return n, err
		}
		if opts.Sort != SortNone {
			if _, err := fmt.Fprintf(bw, "sort %s\n", opts.Sort); err != nil {
				return n, err
			}
		}
		if opts.Collapse {
			if _, err := fmt.Fprintln(bw, "collapse"); err != nil {
				return n, err
			}
		}
		if _, err := fmt.Fprintf(bw, "snapshot %s.%s\n", filename, imgType); err != nil {
			return n, err
		}
		n++
	}
	return n, nil
}
