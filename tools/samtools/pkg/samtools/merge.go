package samtools

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// MergeOptions configures Merge. Defaults to coordinate-sorted input
// (the upstream default); set ByName to merge name-sorted inputs.
type MergeOptions struct {
	// ByName (-n) treats inputs as name-sorted (lexicographic).
	ByName bool
	// ByNameNatural (-n with natural ordering) — exposed separately so the
	// CLI can route -n to the natural variant when SS:queryname:natural is
	// detected on the @HD line.
	ByNameNatural bool
	// HeaderOverride is a path to a SAM header text file; mirrors `-h FILE`.
	HeaderOverride string
	// ForceRGLine, when non-empty, is added to the merged @RG table and
	// applied to every record without an existing RG aux. Mirrors `-r RG`.
	ForceRGLine string
	// CompressLevel picks the BGZF deflate level for the output. -1 keeps
	// the bgzip default.
	CompressLevel int
	// CollapsePG (-c) collapses identical @PG chains to a single line.
	CollapsePG bool
	// PreservePG (-p) keeps every @PG even when collapsable.
	PreservePG bool
	// NoPG suppresses the @PG line injection. By default, merge appends
	// a @PG line documenting the command via the shared InjectPG helper.
	NoPG bool
	// PGCommand is the raw command-line stored under @PG:CL when NoPG
	// is false. The CLI populates this with os.Args.
	PGCommand string
	// Threads is accepted for upstream-CLI compatibility; ignored.
	Threads int
}

// Merge does a streaming k-way merge of inputs into out. Inputs must
// already be sorted in the same order (coordinate by default, by-name
// when ByName is set). The output header is built from the union of
// every input's @SQ table (in first-input order; subsequent inputs must
// match) and every input's @RG and @PG records.
func Merge(inputs []io.Reader, out io.Writer, opts MergeOptions) error {
	if len(inputs) == 0 {
		return errors.New("samtools merge: no input files")
	}

	readers := make([]*sam.BAMReader, 0, len(inputs))
	for i, in := range inputs {
		br, err := sam.NewBAMReader(in)
		if err != nil {
			return fmt.Errorf("samtools merge: input %d: %w", i, err)
		}
		readers = append(readers, br)
	}

	// Header construction.
	hdr := readers[0].Header()
	if opts.HeaderOverride != "" {
		f, err := os.Open(opts.HeaderOverride)
		if err != nil {
			return err
		}
		raw, err := io.ReadAll(f)
		_ = f.Close()
		if err != nil {
			return err
		}
		nh, err := sam.ParseHeaderText(string(raw))
		if err != nil {
			return err
		}
		hdr = nh
	} else {
		// Cross-check @SQ ordering across every input.
		for i := 1; i < len(readers); i++ {
			if !sameRefTable(hdr.Refs, readers[i].Header().Refs) {
				return fmt.Errorf("samtools merge: input %d has a different @SQ table", i)
			}
		}
		// Union @RG entries.
		seenRG := make(map[string]bool, len(hdr.ReadGroups))
		for _, rg := range hdr.ReadGroups {
			seenRG[rg.ID] = true
		}
		for i := 1; i < len(readers); i++ {
			for _, rg := range readers[i].Header().ReadGroups {
				if seenRG[rg.ID] {
					continue
				}
				seenRG[rg.ID] = true
				hdr.ReadGroups = append(hdr.ReadGroups, rg)
				hdr.Lines = append(hdr.Lines, sam.HeaderLine{
					Tag:    "RG",
					Fields: append([]sam.HeaderField{{Tag: "ID", Value: rg.ID}}, rg.Extra...),
				})
			}
		}
		// Union @PG entries.
		seenPG := make(map[string]bool, len(hdr.Programs))
		for _, pg := range hdr.Programs {
			seenPG[pg.ID] = true
		}
		for i := 1; i < len(readers); i++ {
			for _, pg := range readers[i].Header().Programs {
				if seenPG[pg.ID] && (opts.CollapsePG || !opts.PreservePG) {
					continue
				}
				seenPG[pg.ID] = true
				hdr.Programs = append(hdr.Programs, pg)
				hdr.Lines = append(hdr.Lines, sam.HeaderLine{
					Tag:    "PG",
					Fields: append([]sam.HeaderField{{Tag: "ID", Value: pg.ID}}, pg.Extra...),
				})
			}
		}
	}

	// Force-RG handling: parse the line, ensure it lands in the header.
	forcedRGID := ""
	if opts.ForceRGLine != "" {
		newRG, err := parseRGLine(opts.ForceRGLine)
		if err != nil {
			return err
		}
		forcedRGID = newRG.ID
		if findRG(hdr, forcedRGID) < 0 {
			hdr.ReadGroups = append(hdr.ReadGroups, newRG)
			hdr.Lines = append(hdr.Lines, sam.HeaderLine{
				Tag:    "RG",
				Fields: append([]sam.HeaderField{{Tag: "ID", Value: forcedRGID}}, newRG.Extra...),
			})
		}
	}

	// Build the comparator.
	refIndex := make(map[string]int, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refIndex[ref.Name] = i
	}
	var less func(a, b *sam.Record) bool
	switch {
	case opts.ByNameNatural:
		less = func(a, b *sam.Record) bool { return naturalLess(a.QName, b.QName) }
	case opts.ByName:
		less = func(a, b *sam.Record) bool { return a.QName < b.QName }
	default:
		less = func(a, b *sam.Record) bool { return coordLess(a, b, refIndex) }
	}

	bw := sam.NewBAMWriter(out)
	hdr = InjectPG(hdr, "samtools", "samtools", "0.1.0", opts.PGCommand, opts.NoPG)
	if err := bw.WriteHeader(hdr); err != nil {
		return err
	}
	it := newMergeIterator(readers, less)
	for {
		rec, err := it.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if forcedRGID != "" {
			setRecordRG(rec, forcedRGID, AddReplaceRGOverwriteAll)
		}
		if err := bw.Write(rec); err != nil {
			return err
		}
	}
	return bw.Close()
}

// MergeFiles is the high-level CLI entry: opens each input path, runs
// Merge, and closes every file. The output writer is owned by caller.
func MergeFiles(paths []string, out io.Writer, opts MergeOptions) error {
	readers := make([]io.Reader, 0, len(paths))
	closers := make([]io.Closer, 0, len(paths))
	for _, p := range paths {
		f, err := os.Open(p)
		if err != nil {
			for _, c := range closers {
				_ = c.Close()
			}
			return err
		}
		readers = append(readers, f)
		closers = append(closers, f)
	}
	defer func() {
		for _, c := range closers {
			_ = c.Close()
		}
	}()
	return Merge(readers, out, opts)
}

// LoadFOFN reads a "file of file names" — one BAM path per line, with
// '#' comments and blank lines ignored. Used by `samtools merge -b
// FILE-LIST`.
func LoadFOFN(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var out []string
	scan := bufio.NewScanner(f)
	for scan.Scan() {
		line := strings.TrimSpace(scan.Text())
		if line == "" || line[0] == '#' {
			continue
		}
		out = append(out, line)
	}
	return out, scan.Err()
}
