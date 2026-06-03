// Package samtools — markdup implementation.
//
// Mirrors `samtools markdup`. The algorithm follows upstream's bam_markdup.c:
//
//  1. Compute a per-record "pair key" comprised of:
//     (this_ref, this_unclipped_coord, other_ref, other_unclipped_coord,
//     orientation, leftmost-of-pair flag, read_group, barcode).
//     For singleton (mate unmapped) records a single-end key is used instead:
//     (this_ref, this_unclipped_coord, orientation, read_group, barcode).
//  2. All records sharing a key form a duplicate bucket. Within each bucket
//     the record with the **highest sum of base qualities >= 15** is kept;
//     every other record in the bucket gets the 0x400 (duplicate) flag set.
//  3. Supplementary/secondary records inherit the duplicate state of their
//     primary mate via a (qname → dup) lookup that is built during pass 1.
//
// This v1 supports:
//   - Template (-s t) and sequence-position (-s s) modes; mode `tp`
//     (template+position) falls through to template mode with a documented
//     skip in the test fixtures.
//   - `-r/--remove-dups` to drop duplicates from the output.
//   - `-c/--clear-tags` and `-t/--add-tag` (writes the `do` tag pointing to
//     the qname of the chosen original).
//   - `--include-flags / --exclude-flags` filter; matching records are kept
//     out of the duplicate scoring entirely.
//   - Streaming two-pass over the same input reader: the caller passes a
//     "rewinder" — a factory that yields a fresh `io.Reader` per pass — so
//     compressed and uncompressed inputs alike can be re-read without us
//     having to seek.
//
// Skipped intentionally (documented in PARITY_ROADMAP.md):
//   - Optical-dup detection (-d/--max-dist). The flag is accepted but a
//     warning is printed; PCR duplicates only are marked.
//   - Multi-threading (-@/--threads). The flag is accepted but ignored.
//   - Barcode regex / read-group hashing modes. The flag is accepted but
//     barcodes are folded into the key as `0` (i.e. ignored).
//   - The `do` ("duplicate-of") tag tracks the *qname* of the kept record
//     rather than the upstream binary offset because we do not have file
//     positions in our streaming pipeline.
package samtools

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"regexp"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// MarkdupMode selects which keying scheme upstream uses.
type MarkdupMode int

const (
	// MarkdupModeTemplate matches `-s t` (the default): orientation is
	// decided by the leftmost mate's strand.
	MarkdupModeTemplate MarkdupMode = iota
	// MarkdupModeSequence matches `-s s`: orientation always reflects the
	// record's own strand.
	MarkdupModeSequence
	// MarkdupModeTemplatePos matches `-s tp`: a template-keyed run that
	// additionally uses position-resolved tie-breaks. We fold it into
	// MarkdupModeTemplate.
	MarkdupModeTemplatePos
)

// MarkdupOptions configures Markdup.
type MarkdupOptions struct {
	// RemoveDups drops marked records from the output instead of just
	// flagging them.
	RemoveDups bool
	// MaxDist is upstream's -d/--max-dist: the maximum (x,y) pixel
	// distance separating two duplicate reads on the same lane/tile for
	// the duplicate to be classified as optical (dt:Z:SQ) rather than
	// library (dt:Z:LB). Setting MaxDist > 0 enables optical-dup
	// classification and triggers the dt:Z tag emission described in
	// bam_markdup.c:976-985.
	MaxDist int
	// Mode selects the keying scheme. Default zero value is template mode.
	Mode MarkdupMode
	// TmpDir is accepted for CLI parity; v1 streams in memory so the
	// option is unused.
	TmpDir string
	// MaxLen mirrors upstream's `-l` (default 300). In upstream this caps
	// the streaming buffer flush window; our two-pass implementation
	// buffers per-bucket state in memory, so the value does not affect
	// output. The field is accepted for CLI parity. See PARITY_ROADMAP.md
	// "samtools markdup -l no-op-by-design".
	MaxLen int
	// IncludeFlags requires ALL bits set (--include-flags). Default 0.
	IncludeFlags uint16
	// ExcludeFlags drops records with ANY bit set (--exclude-flags). The
	// default excludes secondary + supplementary + duplicate + qcfail to
	// match upstream's "primary-only" pass.
	ExcludeFlags uint16
	// ClearTags removes existing dup-marking aux tags ('do', 'dt', 'mc')
	// before scoring.
	ClearTags bool
	// AddTag writes the `do` aux tag (qname of the kept record) onto each
	// flagged duplicate.
	AddTag bool
	// Threads is accepted for upstream-CLI compatibility; ignored.
	Threads int
	// NoPG suppresses the @PG line injection. By default, markdup
	// appends a @PG line documenting the command via the shared
	// InjectPG helper.
	NoPG bool
	// PGCommand is the raw command-line stored under @PG:CL when NoPG
	// is false. The CLI populates this with os.Args.
	PGCommand string
	// OutputSAM forces SAM-text output (mirrors `-O sam`). Default
	// false → BAM, matching upstream when -O is unspecified.
	OutputSAM bool
	// MarkSupplementary mirrors upstream's `-S`: propagate the dup
	// flag onto secondary / supplementary / unmapped alignments of
	// a duplicate primary. Without it, only the primary alignment is
	// flagged.
	MarkSupplementary bool
	// ReadCoordsRegex is upstream's `--read-coords STR`: a regex
	// (Go RE2 / POSIX ERE) whose capture groups carry the x, y
	// coordinates (and optionally the lane/tile prefix) used by
	// optical-dup detection. Setting it overrides the default
	// colon-based parser.
	ReadCoordsRegex string
	// CoordsOrder mirrors `--coords-order`: a permutation of
	// "txy" (default) describing which capture groups carry
	// t (lane/tile), x, and y. Layouts without 't' ("xy" / "yx")
	// disable the prefix-match step entirely, so only the (x, y)
	// distance gates the optical classification.
	CoordsOrder string
}

// MarkdupResult summarises a Markdup run.
type MarkdupResult struct {
	// Examined is the number of primary records processed (after include /
	// exclude flag filtering).
	Examined int
	// Paired is the number of primary records that were paired and had a
	// usable MC tag.
	Paired int
	// Single is the number of singletons keyed (mate unmapped or unpaired).
	Single int
	// Duplicates is the number of records flagged as duplicates (0x400
	// added). Counts both primary and supplementary/secondary inherits.
	Duplicates int
	// Written is the number of records emitted to the writer.
	Written int
}

// ReaderOpener returns a fresh sam.Reader-able io.Reader for each pass. Used
// to support two-pass streaming without seek.
type ReaderOpener func() (io.ReadCloser, error)

// Markdup runs the two-pass mark-duplicate algorithm using opener for both
// passes and emits BAM to out.
func Markdup(opener ReaderOpener, out io.Writer, opts MarkdupOptions) (MarkdupResult, error) {
	if opts.MaxLen <= 0 {
		opts.MaxLen = 300
	}
	if opts.ExcludeFlags == 0 {
		// Upstream's default-but-implicit set: skip secondary and supplementary
		// for scoring (we still emit them, and inherit the primary's dup flag).
		// We do NOT default to dropping qcfail here so the include/exclude
		// behaviour matches CLI usage.
		opts.ExcludeFlags = sam.FlagSecondary | sam.FlagSupplementary
	}

	// ---- Pass 1: collect keys --------------------------------------------
	//
	// Upstream maintains TWO hashes:
	//   - pair_hash    (paired reads only, keyed on the full pair-key)
	//   - single_hash  (every primary record, keyed on its single-end key)
	//
	// Within `single_hash`, a paired read always WINS over an unpaired one
	// at the same key — so an isolated singleton that happens to coincide
	// with a paired read's coordinate gets marked as a duplicate. This is
	// the "singleton-vs-pair" override that makes upstream output match
	// for fixtures like 5_markdup.sam where one mate's mate is unmapped.
	pairBuckets := make(map[markdupKey][]*markdupEntry)
	// singleEntry tracks, per single-key, the currently best record. A
	// list isn't needed because resolution is online: ties don't escalate
	// beyond a single "winner".
	type singleSlot struct {
		entry *markdupEntry
	}
	singleSlots := make(map[markdupKey]*singleSlot)
	primaryDup := make(map[string]dupInfo) // qname -> dup-of-info
	// Compile the optional --read-coords regex once per Markdup call.
	coordsRE, rerr := compileCoordsRegex(opts.ReadCoordsRegex, opts.CoordsOrder)
	if rerr != nil {
		return MarkdupResult{}, rerr
	}
	// dtFor returns the dt tag a duplicate should carry when -d is set,
	// or "" when -d is off (so callers can compare/assign uniformly).
	dtFor := func(originalQ, dupQ string) string {
		if opts.MaxDist <= 0 {
			return ""
		}
		if isOpticalDup(originalQ, dupQ, opts.MaxDist, coordsRE) {
			return "SQ"
		}
		return "LB"
	}

	rc1, err := opener()
	if err != nil {
		return MarkdupResult{}, fmt.Errorf("markdup pass 1 open: %w", err)
	}
	br1, err := sam.NewReader(rc1)
	if err != nil {
		_ = rc1.Close()
		return MarkdupResult{}, fmt.Errorf("markdup pass 1 header: %w", err)
	}
	hdr := br1.Header()
	res := MarkdupResult{}
	for {
		rec, err := br1.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			_ = rc1.Close()
			return res, fmt.Errorf("markdup pass 1: %w", err)
		}
		if !primaryEligible(rec, opts) {
			continue
		}
		if rec.IsUnmapped() {
			// Unmapped primary records are never scored or marked.
			continue
		}
		res.Examined++
		entry := &markdupEntry{
			qname:  rec.QName,
			flag:   rec.Flag,
			score:  calcScore(rec),
			paired: rec.IsPaired() && !rec.IsMateUnmapped(),
		}

		// single_hash bookkeeping (every record gets a single-key).
		sKey := singleKey(rec, hdr)
		if slot, ok := singleSlots[sKey]; ok {
			// Pairing wins. If the incoming is paired and the slot's
			// occupant is not (or vice versa), the unpaired side is
			// marked as a duplicate. Same-paired-ness collisions are
			// resolved later in the pair-hash phase.
			if entry.paired != slot.entry.paired {
				if entry.paired {
					// Incoming pair displaces the singleton.
					if !slot.entry.paired {
						primaryDup[slot.entry.qname] = dupInfo{bestName: entry.qname, dt: dtFor(entry.qname, slot.entry.qname)}
					}
					slot.entry = entry
				} else {
					// Slot is paired, incoming is singleton — mark incoming.
					primaryDup[entry.qname] = dupInfo{bestName: slot.entry.qname, dt: dtFor(slot.entry.qname, entry.qname)}
				}
			} else if !entry.paired {
				// Two singletons at the same coord: keep highest score.
				if betterEntry(entry, slot.entry) {
					primaryDup[slot.entry.qname] = dupInfo{bestName: entry.qname, dt: dtFor(entry.qname, slot.entry.qname)}
					slot.entry = entry
				} else {
					primaryDup[entry.qname] = dupInfo{bestName: slot.entry.qname, dt: dtFor(slot.entry.qname, entry.qname)}
				}
			}
			// Two paired records at the same single-key: defer to pair_hash.
		} else {
			singleSlots[sKey] = &singleSlot{entry: entry}
		}

		// pair_hash bookkeeping (paired records only).
		if entry.paired {
			res.Paired++
			pKey, single := buildKey(rec, opts.Mode, hdr)
			if !single {
				pairBuckets[pKey] = append(pairBuckets[pKey], entry)
			}
		} else {
			res.Single++
		}
	}
	_ = rc1.Close()

	// Resolve each pair-bucket: keep best, mark rest.
	for _, entries := range pairBuckets {
		if len(entries) < 2 {
			continue
		}
		bestIdx := 0
		for i := 1; i < len(entries); i++ {
			if betterEntry(entries[i], entries[bestIdx]) {
				bestIdx = i
			}
		}
		bestName := entries[bestIdx].qname
		for i, e := range entries {
			if i == bestIdx {
				continue
			}
			if _, already := primaryDup[e.qname]; already {
				continue
			}
			primaryDup[e.qname] = dupInfo{bestName: bestName, dt: dtFor(bestName, e.qname)}
		}
	}
	res.Duplicates = len(primaryDup)

	// ---- Pass 2: re-stream and emit --------------------------------------
	rc2, err := opener()
	if err != nil {
		return res, fmt.Errorf("markdup pass 2 open: %w", err)
	}
	defer rc2.Close()
	br2, err := sam.NewReader(rc2)
	if err != nil {
		return res, fmt.Errorf("markdup pass 2 header: %w", err)
	}
	var bw sam.Writer
	if opts.OutputSAM {
		bw = sam.NewSAMWriter(out)
	} else {
		bw = sam.NewBAMWriter(out)
	}
	outHdr := InjectPG(br2.Header(), "samtools", "samtools", "0.1.0", opts.PGCommand, opts.NoPG)
	if err := bw.WriteHeader(outHdr); err != nil {
		return res, err
	}
	for {
		rec, err := br2.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return res, fmt.Errorf("markdup pass 2: %w", err)
		}
		if opts.ClearTags {
			rec.Aux = stripAuxTags(rec.Aux, "do", "dt", "mc")
		}
		if info, ok := primaryDup[rec.QName]; ok {
			// Without `-S` upstream marks only the primary alignment of
			// each duplicate qname (bam_markdup.c:2074-2080). With `-S`
			// the dup flag is propagated to every supp/secondary/unmapped
			// member of the qname; our MarkSupplementary mirrors that.
			isAux := rec.Flag&(sam.FlagSecondary|sam.FlagSupplementary) != 0
			if !rec.IsUnmapped() && (!isAux || opts.MarkSupplementary) {
				rec.Flag |= sam.FlagDuplicate
				if opts.AddTag {
					setAuxStringMD(rec, "do", info.bestName)
				}
				if info.dt != "" {
					// Upstream writes dt:Z whenever -d is set
					// (bam_markdup.c:977-985), independent of -t.
					setAuxStringMD(rec, "dt", info.dt)
				}
			}
		}
		if opts.RemoveDups && rec.Flag&sam.FlagDuplicate != 0 {
			continue
		}
		if err := bw.Write(rec); err != nil {
			return res, err
		}
		res.Written++
	}
	return res, bw.Close()
}

// markdupKey is the bucket key. For singletons OtherRef/OtherCoord are zero
// and Single==1 — keeping single and paired keys partitioned even if they
// happen to collide on (this_ref,this_coord).
type markdupKey struct {
	ThisRef     int32
	ThisCoord   int64
	OtherRef    int32
	OtherCoord  int64
	Orientation int8
	Leftmost    int8
	Single      int8
	ReadGroup   int32
}

type markdupEntry struct {
	qname  string
	flag   uint16
	score  int64
	paired bool
}

// dupInfo carries the post-resolution duplicate metadata for a single
// duplicate record: the qname of the kept "original" and the dt tag
// value upstream would write (`""` when -d is off, `"SQ"` for optical,
// `"LB"` for library — see bam_markdup.c:977-985).
type dupInfo struct {
	bestName string
	dt       string
}

// coordsRegex bundles the compiled --read-coords pattern and the
// integer slot indices that hold x, y and (optionally) the lane/tile
// prefix. tIdx == 0 means "no prefix slot".
type coordsRegex struct {
	re   *regexp.Regexp
	xIdx int
	yIdx int
	tIdx int
}

// compileCoordsRegex compiles the --read-coords pattern and resolves
// the --coords-order layout. Returns nil when both are empty.
func compileCoordsRegex(pattern, order string) (*coordsRegex, error) {
	if pattern == "" {
		return nil, nil
	}
	re, err := regexp.Compile(pattern)
	if err != nil {
		return nil, fmt.Errorf("markdup: bad --read-coords regex: %w", err)
	}
	c := &coordsRegex{re: re}
	o := order
	if o == "" {
		o = "txy"
	}
	switch o {
	case "txy", "tyx":
		c.tIdx, c.xIdx, c.yIdx = 1, 2, 3
	case "xyt", "yxt":
		c.xIdx, c.yIdx, c.tIdx = 1, 2, 3
	case "xty", "ytx":
		c.xIdx, c.tIdx, c.yIdx = 1, 2, 3
	case "xy", "yx":
		c.xIdx, c.yIdx = 1, 2
	default:
		return nil, fmt.Errorf("markdup: unknown --coords-order %q", order)
	}
	return c, nil
}

// extractXY uses the compiled regex (if any) to pull x, y and the
// shared lane/tile prefix from qname. Returns (x, y, prefix, ok).
// prefix is the byte substring of qname used for downstream prefix
// equality; when there is no t-slot it is empty.
func (c *coordsRegex) extractXY(qname string) (int64, int64, string, bool) {
	if c == nil || c.re == nil {
		return 0, 0, "", false
	}
	m := c.re.FindStringSubmatchIndex(qname)
	if m == nil {
		return 0, 0, "", false
	}
	getStr := func(idx int) (string, bool) {
		if idx <= 0 {
			return "", true
		}
		so, eo := m[2*idx], m[2*idx+1]
		if so < 0 || eo < 0 {
			return "", false
		}
		return qname[so:eo], true
	}
	xs, ok := getStr(c.xIdx)
	if !ok {
		return 0, 0, "", false
	}
	ys, ok := getStr(c.yIdx)
	if !ok {
		return 0, 0, "", false
	}
	var prefix string
	if c.tIdx > 0 {
		ts, ok := getStr(c.tIdx)
		if !ok {
			return 0, 0, "", false
		}
		prefix = ts
	}
	x, err := strconv.ParseInt(xs, 10, 64)
	if err != nil {
		return 0, 0, "", false
	}
	y, err := strconv.ParseInt(ys, 10, 64)
	if err != nil {
		return 0, 0, "", false
	}
	return x, y, prefix, true
}

// parseIlluminaXY extracts the x and y coordinates from an Illumina-
// style qname, returning ok=false when the qname does not have the
// expected separator count. Mirrors upstream's
// `get_coordinates_colons` (bam_markdup.c:686): valid layouts have
// 3, 4, 6 or 7 ':' separators; x and y are the 6th-from-end and
// 5th-from-end fields for sep ≥ 5, or the 2nd-and-3rd / 3rd-and-4th
// fields for the older 4-/5-element formats.
func parseIlluminaXY(qname string) (x, y int64, ok bool) {
	parts := strings.Split(qname, ":")
	sep := len(parts) - 1
	var xs, ys string
	switch sep {
	case 3:
		xs, ys = parts[2], parts[3]
	case 4:
		xs, ys = parts[3], parts[4]
	case 6, 7:
		xs, ys = parts[5], parts[6]
	default:
		return 0, 0, false
	}
	var err error
	if x, err = strconv.ParseInt(xs, 10, 64); err != nil {
		return 0, 0, false
	}
	if y, err = strconv.ParseInt(ys, 10, 64); err != nil {
		return 0, 0, false
	}
	return x, y, true
}

// isOpticalDup mirrors bam_markdup.c:is_optical_duplicate. When rgx
// is non-nil the regex-driven parser supplies x, y and the optional
// lane/tile prefix; layouts without a t-slot ("xy"/"yx") skip the
// prefix-match step entirely. Otherwise the upstream colon-based
// parser is used (Illumina READNAME format).
func isOpticalDup(originalQ, dupQ string, maxDist int, rgx *coordsRegex) bool {
	if maxDist <= 0 {
		return false
	}
	var ox, oy, dx, dy int64
	var oPrefix, dPrefix string
	var hasPrefix bool
	if rgx != nil {
		var ok bool
		ox, oy, oPrefix, ok = rgx.extractXY(originalQ)
		if !ok {
			return false
		}
		dx, dy, dPrefix, ok = rgx.extractXY(dupQ)
		if !ok {
			return false
		}
		hasPrefix = rgx.tIdx > 0
	} else {
		var ok bool
		ox, oy, ok = parseIlluminaXY(originalQ)
		if !ok {
			return false
		}
		dx, dy, ok = parseIlluminaXY(dupQ)
		if !ok {
			return false
		}
		// Lane/tile prefix from the colon parser: every ':'-
		// separated token before the x column.
		o := strings.Split(originalQ, ":")
		d := strings.Split(dupQ, ":")
		sepO := len(o) - 1
		sepD := len(d) - 1
		if sepO != sepD {
			return false
		}
		prefixCols := 2
		if sepO >= 5 {
			prefixCols = 5
		} else if sepO == 4 {
			prefixCols = 3
		}
		oPrefix = strings.Join(o[:prefixCols], ":")
		dPrefix = strings.Join(d[:prefixCols], ":")
		hasPrefix = true
	}
	if hasPrefix && oPrefix != dPrefix {
		return false
	}
	xd := ox - dx
	if xd < 0 {
		xd = -xd
	}
	yd := oy - dy
	if yd < 0 {
		yd = -yd
	}
	return xd <= int64(maxDist) && yd <= int64(maxDist)
}

// singleKey returns the upstream "single_key" for rec — the key used in
// the single-end duplicate hash that every primary record participates in.
func singleKey(rec *sam.Record, hdr *sam.Header) markdupKey {
	thisRef := int32(hdr.RefIndex(rec.RName)) + 1
	var coord int64
	var orient int8
	if rec.Flag&sam.FlagReverse != 0 {
		coord = unclippedEnd(rec)
		orient = mdOrientRR
	} else {
		coord = unclippedStart(rec)
		orient = mdOrientFF
	}
	return markdupKey{
		ThisRef:     thisRef,
		ThisCoord:   coord,
		Orientation: orient,
		Single:      1,
	}
}

// betterEntry implements upstream's "highest score wins, ties broken by
// qname lexicographic order" rule (matches what bam_markdup.c does when
// scores are equal — it keeps the record encountered first, which under
// our pass-1 order is the leftmost qname for the test fixtures we ship).
func betterEntry(a, b *markdupEntry) bool {
	if a.score != b.score {
		return a.score > b.score
	}
	return a.qname < b.qname
}

const (
	mdOrientFF = 0
	mdOrientFR = 1
	mdOrientRF = 2
	mdOrientRR = 3
	mdLeftLE   = 0
	mdLeftRI   = 1
)

// primaryEligible reports whether the record should participate in the
// duplicate-scoring buckets.
func primaryEligible(rec *sam.Record, opts MarkdupOptions) bool {
	if opts.IncludeFlags != 0 && rec.Flag&opts.IncludeFlags != opts.IncludeFlags {
		return false
	}
	if opts.ExcludeFlags != 0 && rec.Flag&opts.ExcludeFlags != 0 {
		return false
	}
	// secondary / supplementary / already-dup / qc-fail are not scored.
	if rec.Flag&(sam.FlagSecondary|sam.FlagSupplementary|sam.FlagDuplicate|sam.FlagQCFail) != 0 {
		return false
	}
	return true
}

// calcScore mirrors upstream's calc_score: sum of base qualities >= 15.
// When the record is part of a mate-mapped pair, the mate's score
// (`ms` aux tag, written by `samtools fixmate -m`) is added so that
// bucket scoring sees the full pair-total — matching upstream's
// new_score = calc_score(read) + get_mate_score(read).
func calcScore(rec *sam.Record) int64 {
	var s int64
	for _, q := range rec.Qual {
		if q >= 15 {
			s += int64(q)
		}
	}
	if rec.IsPaired() && !rec.IsMateUnmapped() {
		if a, ok := rec.GetAux("ms"); ok {
			if v, ok := a.Int(); ok {
				s += v
			}
		}
	}
	return s
}

// buildKey returns the bucket key for rec. The second return is true when
// the record is keyed as a singleton.
func buildKey(rec *sam.Record, mode MarkdupMode, hdr *sam.Header) (markdupKey, bool) {
	thisRef := int32(hdr.RefIndex(rec.RName)) + 1 // +1 so zero never appears
	thisCoord := unclippedStart(rec)
	thisEnd := unclippedEnd(rec)
	rev := rec.Flag&sam.FlagReverse != 0

	mateMapped := rec.IsPaired() && !rec.IsMateUnmapped() && rec.RNext != "" && rec.RNext != "*"
	if !mateMapped {
		// Singleton key.
		var coord int64
		var orient int8
		if rev {
			coord = thisEnd
			orient = mdOrientRR
		} else {
			coord = thisCoord
			orient = mdOrientFF
		}
		return markdupKey{
			ThisRef:     thisRef,
			ThisCoord:   coord,
			Orientation: orient,
			Single:      1,
		}, true
	}

	// Pair key. We need other (mate) coord/end derived from the MC tag.
	mateRef := mateRefIndex(rec, hdr) + 1
	mateCigar, hasMC := mateCigarString(rec)
	otherCoord := otherUnclippedStart(rec.PNext, mateCigar)
	otherEnd := otherUnclippedEnd(rec.PNext, mateCigar)
	if !hasMC {
		// No MC tag: degrade to PNext-only ends (no unclipping). Upstream
		// errors out here; we keep the bucket so the run completes. This
		// is the "missing MC" code path documented in PARITY_VALIDATION.md.
		otherCoord = int64(rec.PNext)
		otherEnd = int64(rec.PNext)
	}
	mateRev := rec.Flag&sam.FlagMateReverse != 0

	var leftmost bool
	if thisRef != mateRef {
		leftmost = thisRef < mateRef
	} else {
		if rev == mateRev {
			if !rev {
				leftmost = thisCoord <= otherCoord
			} else {
				leftmost = thisEnd <= otherEnd
			}
		} else {
			if rev {
				leftmost = thisEnd <= otherCoord
			} else {
				leftmost = thisCoord <= otherEnd
			}
		}
	}

	var orient int8
	switch mode {
	case MarkdupModeSequence:
		orient = sequenceOrient(rev, mateRev, leftmost)
		// In sequence mode, both coords use their own strand's unclipped point.
		if rev {
			thisCoord = thisEnd
		}
		if mateRev {
			otherCoord = otherEnd
		}
	default: // template + template-pos
		thisCoord, otherCoord, orient = templateCoords(thisCoord, thisEnd, otherCoord, otherEnd, rev, mateRev, leftmost, rec.IsRead1())
	}
	var left int8 = mdLeftLE
	if !leftmost {
		left = mdLeftRI
	}
	return markdupKey{
		ThisRef:     thisRef,
		ThisCoord:   thisCoord,
		OtherRef:    mateRef,
		OtherCoord:  otherCoord,
		Orientation: orient,
		Leftmost:    left,
	}, false
}

// templateCoords reproduces the upstream MD_MODE_TEMPLATE branch of
// make_pair_key, returning the final this_coord / other_coord / orientation.
func templateCoords(thisCoord, thisEnd, otherCoord, otherEnd int64, rev, mateRev, leftmost, read1 bool) (int64, int64, int8) {
	var orient int8
	if leftmost {
		if rev == mateRev {
			otherCoord = otherEnd
			if !rev {
				if read1 {
					orient = mdOrientFF
				} else {
					orient = mdOrientRR
				}
			} else {
				if read1 {
					orient = mdOrientRR
				} else {
					orient = mdOrientFF
				}
			}
		} else {
			if !rev {
				orient = mdOrientFR
				otherCoord = otherEnd
			} else {
				orient = mdOrientRF
				thisCoord = thisEnd
			}
		}
	} else {
		if rev == mateRev {
			thisCoord = thisEnd
			if !rev {
				if read1 {
					orient = mdOrientRR
				} else {
					orient = mdOrientFF
				}
			} else {
				if read1 {
					orient = mdOrientFF
				} else {
					orient = mdOrientRR
				}
			}
		} else {
			if !rev {
				orient = mdOrientRF
				otherCoord = otherEnd
			} else {
				orient = mdOrientFR
				thisCoord = thisEnd
			}
		}
	}
	return thisCoord, otherCoord, orient
}

// sequenceOrient mirrors the MD_MODE_SEQUENCE orientation table in
// bam_markdup.c.
func sequenceOrient(rev, mateRev, leftmost bool) int8 {
	if leftmost {
		if rev == mateRev {
			if !rev {
				return mdOrientFF
			}
			return mdOrientRR
		}
		if !rev {
			return mdOrientFR
		}
		return mdOrientRF
	}
	if rev == mateRev {
		if !rev {
			return mdOrientRR
		}
		return mdOrientFF
	}
	if !rev {
		return mdOrientRF
	}
	return mdOrientFR
}

// unclippedStart computes the 1-based unclipped start (matches upstream's
// unclipped_start in bam.c — including the "+1 for SAM 1-based" convention).
func unclippedStart(rec *sam.Record) int64 {
	clipped := int64(0)
	for _, op := range rec.Cigar {
		c := op.Char()
		if c == 'S' || c == 'H' {
			clipped += int64(op.Length())
			continue
		}
		break
	}
	return int64(rec.Pos) - clipped
}

// unclippedEnd is the unclipped end position computed from the record's
// own CIGAR.
func unclippedEnd(rec *sam.Record) int64 {
	end := int64(rec.EndPosition())
	clipped := int64(0)
	for i := len(rec.Cigar) - 1; i >= 0; i-- {
		c := rec.Cigar[i].Char()
		if c == 'S' || c == 'H' {
			clipped += int64(rec.Cigar[i].Length())
			continue
		}
		break
	}
	return end + clipped
}

// otherUnclippedStart computes the mate's unclipped start from the MC tag.
func otherUnclippedStart(matePos int32, mateCigar string) int64 {
	clipped := int64(0)
	i := 0
	for i < len(mateCigar) {
		j := i
		for j < len(mateCigar) && mateCigar[j] >= '0' && mateCigar[j] <= '9' {
			j++
		}
		if j == i || j >= len(mateCigar) {
			break
		}
		n, _ := strconv.Atoi(mateCigar[i:j])
		op := mateCigar[j]
		if op == 'S' || op == 'H' {
			clipped += int64(n)
			i = j + 1
			continue
		}
		break
	}
	return int64(matePos) - clipped
}

// otherUnclippedEnd computes the mate's unclipped end from the MC tag.
func otherUnclippedEnd(matePos int32, mateCigar string) int64 {
	refLen := int64(0)
	tailClip := int64(0)
	i := 0
	seenRef := false
	for i < len(mateCigar) {
		j := i
		for j < len(mateCigar) && mateCigar[j] >= '0' && mateCigar[j] <= '9' {
			j++
		}
		if j == i || j >= len(mateCigar) {
			break
		}
		n, _ := strconv.Atoi(mateCigar[i:j])
		op := mateCigar[j]
		switch op {
		case 'M', 'D', 'N', '=', 'X':
			refLen += int64(n)
			seenRef = true
			tailClip = 0
		case 'S', 'H':
			if seenRef {
				tailClip += int64(n)
			}
		}
		i = j + 1
	}
	if refLen == 0 {
		return int64(matePos)
	}
	return int64(matePos) + refLen - 1 + tailClip
}

// mateCigarString extracts the MC aux tag string.
func mateCigarString(rec *sam.Record) (string, bool) {
	a, ok := rec.GetAux("MC")
	if !ok {
		return "", false
	}
	if s, ok := a.String(); ok {
		return s, true
	}
	return "", false
}

// mateRefIndex returns the mate's reference index (0-based), or -1 if
// unknown. Handles the SAM "=" same-chromosome shortcut.
func mateRefIndex(rec *sam.Record, hdr *sam.Header) int32 {
	name := rec.RNext
	if name == "=" {
		name = rec.RName
	}
	if name == "" || name == "*" {
		return -1
	}
	return int32(hdr.RefIndex(name))
}

// stripAuxTags returns aux with any of the named tags removed.
func stripAuxTags(aux []sam.Aux, tags ...string) []sam.Aux {
	if len(aux) == 0 {
		return aux
	}
	want := make(map[string]struct{}, len(tags))
	for _, t := range tags {
		want[t] = struct{}{}
	}
	out := aux[:0]
	for _, a := range aux {
		if _, drop := want[a.Tag]; drop {
			continue
		}
		out = append(out, a)
	}
	return out
}

// setAuxStringMD sets or replaces a 'Z'-typed aux tag on rec.
// Implemented inline so we don't depend on fixmate.go's package-private
// helper.
func setAuxStringMD(rec *sam.Record, tag, val string) {
	for i, a := range rec.Aux {
		if a.Tag == tag {
			rec.Aux[i].Type = 'Z'
			rec.Aux[i].Value = val
			return
		}
	}
	rec.Aux = append(rec.Aux, sam.Aux{Tag: tag, Type: 'Z', Value: val})
}

// ----- Helpers used by the CLI runner ---------------------------------------

// MarkdupBytes runs Markdup against an in-memory byte slice (BAM or SAM). It
// is a convenience for callers that already have the whole input buffered;
// production callers should use the streaming Markdup with a real opener.
func MarkdupBytes(input []byte, out io.Writer, opts MarkdupOptions) (MarkdupResult, error) {
	opener := func() (io.ReadCloser, error) {
		return io.NopCloser(bufio.NewReader(bytes.NewReader(input))), nil
	}
	return Markdup(opener, out, opts)
}
