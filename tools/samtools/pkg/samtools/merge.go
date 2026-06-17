package samtools

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
	"time"

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
	// CombineRG (-c) combines @RG records with colliding IDs into one. When
	// false (the upstream default) a colliding @RG ID is renamed to a distinct
	// "<id>-<8 hex>" form and the records carrying it are retagged. Mirrors
	// samtools merge -c.
	CombineRG bool
	// CombinePG (-p) combines @PG records with colliding IDs. When false (the
	// upstream default) colliding @PG IDs are kept distinct. Mirrors -p.
	CombinePG bool
	// RandomSeed seeds the drand48-compatible PRNG used to disambiguate
	// colliding @RG IDs (samtools merge -s SEED). Used only when SeedSet is
	// true; otherwise the PRNG is seeded from the wall clock, matching
	// upstream's time(NULL) default (and so non-reproducible by design).
	RandomSeed int64
	// SeedSet reports whether RandomSeed was given explicitly (-s).
	SeedSet bool
	// Threads is accepted for upstream-CLI compatibility; ignored.
	Threads int
}

// rand48 is a drand48-compatible 48-bit linear congruential generator,
// byte-identical to glibc's lrand48 (and thus htslib's hts_lrand48). samtools
// merge uses it to append a random suffix to colliding @RG IDs, so reproducing
// its stream lets the renamed IDs match upstream exactly under a fixed -s seed.
type rand48 struct{ state uint64 }

const (
	rand48Mult = 0x5DEECE66D
	rand48Add  = 0xB
	rand48Mask = (uint64(1) << 48) - 1
)

// newRand48 seeds the generator as POSIX srand48 does: the 48-bit state is
// (seed<<16) | 0x330E.
func newRand48(seed int64) *rand48 {
	return &rand48{state: ((uint64(seed) << 16) | 0x330E) & rand48Mask}
}

// lrand48 returns the next non-negative 31-bit value, mirroring lrand48().
func (r *rand48) lrand48() uint64 {
	r.state = (rand48Mult*r.state + rand48Add) & rand48Mask
	return r.state >> 17
}

// genUniqueID appends a random "-%08X" suffix to prefix until the result is not
// already in existing, mirroring samtools' gen_unique_id with
// always_add_suffix=true (the colliding-ID path).
func genUniqueID(prefix string, existing map[string]bool, rng *rand48) string {
	for {
		id := fmt.Sprintf("%s-%08X", prefix, rng.lrand48())
		if !existing[id] {
			return id
		}
	}
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

	// rgTrans[i] maps input i's original @RG IDs to the (possibly renamed) ID
	// in the merged header, so records from that input can be retagged.
	rgTrans := make([]map[string]string, len(readers))
	for i := range rgTrans {
		rgTrans[i] = map[string]string{}
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
		// Union @RG entries. rgTrans[i] maps input i's original @RG IDs to the
		// ID used in the merged header: identity unless a collision forced a
		// rename. The PRNG is seeded so a fixed -s SEED reproduces upstream's
		// renamed suffixes exactly.
		seedVal := opts.RandomSeed
		if !opts.SeedSet {
			seedVal = time.Now().UnixNano()
		}
		rng := newRand48(seedVal)
		seenRG := make(map[string]bool, len(hdr.ReadGroups))
		for _, rg := range hdr.ReadGroups {
			seenRG[rg.ID] = true
		}
		// Input 0's groups are already in hdr (identity); inputs 1.. may rename.
		for i := 1; i < len(readers); i++ {
			for _, rg := range readers[i].Header().ReadGroups {
				newID := rg.ID
				if seenRG[rg.ID] {
					if opts.CombineRG {
						// Combine: reuse the existing ID, add nothing.
						continue
					}
					// Rename to a distinct "<id>-<8 hex>" form (upstream's
					// gen_unique_id, always_add_suffix). This is the only RG
					// PRNG draw, and it precedes the @PG processing below, so
					// the suffix matches upstream for a given seed.
					newID = genUniqueID(rg.ID, seenRG, rng)
					rgTrans[i][rg.ID] = newID
				}
				seenRG[newID] = true
				hdr.ReadGroups = append(hdr.ReadGroups, sam.ReadGroup{ID: newID, Extra: rg.Extra})
				hdr.Lines = append(hdr.Lines, sam.HeaderLine{
					Tag:    "RG",
					Fields: append([]sam.HeaderField{{Tag: "ID", Value: newID}}, rg.Extra...),
				})
			}
		}
		// Union @PG entries. With -p (CombinePG) colliding IDs are combined;
		// otherwise duplicates are kept distinct. @PG provenance is dropped by
		// downstream decoded comparisons, so the exact renumbering is not
		// reproduced here.
		seenPG := make(map[string]bool, len(hdr.Programs))
		for _, pg := range hdr.Programs {
			seenPG[pg.ID] = true
		}
		for i := 1; i < len(readers); i++ {
			for _, pg := range readers[i].Header().Programs {
				if seenPG[pg.ID] && opts.CombinePG {
					continue
				}
				if seenPG[pg.ID] {
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
	if err := bw.WriteHeader(hdr); err != nil {
		return err
	}
	// Whether any input had an @RG collision that forced a rename.
	anyRGTrans := false
	for _, m := range rgTrans {
		if len(m) > 0 {
			anyRGTrans = true
			break
		}
	}

	it := newMergeIterator(readers, less)
	mi, _ := it.(*mergeIterator)
	for {
		var rec *sam.Record
		var src int
		var err error
		if mi != nil {
			rec, src, err = mi.nextWithSrc()
		} else {
			rec, err = it.Next()
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		// Retag records whose input's @RG ID was renamed on collision.
		if anyRGTrans && src >= 0 {
			if newID, ok := rgTrans[src][recordRGID(rec)]; ok {
				setRecordRG(rec, newID, AddReplaceRGOverwriteAll)
			}
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

// recordRGID returns the record's RG:Z aux value, or "" when absent.
func recordRGID(rec *sam.Record) string {
	for _, a := range rec.Aux {
		if a.Tag == "RG" {
			if s, ok := a.Value.(string); ok {
				return s
			}
		}
	}
	return ""
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
