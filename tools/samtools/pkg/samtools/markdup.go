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
//  3. By default supplementary/secondary records are emitted untouched. Under
//     `-S` they inherit the duplicate state of their primary mate via a
//     (qname → dup) lookup built during pass 1.
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
//   - `-S` (mark supplementary/secondary duplicates) and the upstream default
//     of leaving non-primary records untouched.
//   - `-s/--stats` and `-f FILE` duplicate-statistics reporting, including the
//     Picard-style ESTIMATED_LIBRARY_SIZE solve.
//   - `-d` optical-duplicate detection: each flagged duplicate is compared to
//     the chosen original by Illumina read-name tile coordinates and tagged
//     dt:Z:SQ (optical) or dt:Z:LB (library) accordingly.
//   - Streaming two-pass over the same input reader: the caller passes a
//     "rewinder" — a factory that yields a fresh `io.Reader` per pass — so
//     compressed and uncompressed inputs alike can be re-read without us
//     having to seek.
//
// Skipped intentionally (documented in docs/PARITY_ROADMAP.md):
//   - Barcode regex / read-group hashing modes. The flag is accepted but
//     barcodes are folded into the key as `0` (i.e. ignored).
//   - The `do` ("duplicate-of") tag tracks the *qname* of the kept record
//     rather than the upstream binary offset because we do not have file
//     positions in our streaming pipeline.
//   - Optical detection uses the pairwise original-vs-duplicate comparison;
//     upstream's `--read-coords` regex and the check_chain whole-cluster
//     re-expansion (which can promote a duplicate-of-a-duplicate to optical)
//     are not reproduced. Read names must use the colon-delimited Illumina
//     layout for `-d` to take effect.
package samtools

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"

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
	// MaxDist is upstream's `-d` optical-duplicate distance. When greater
	// than zero, each flagged duplicate is compared to the chosen original by
	// Illumina read-name tile coordinates; if both axes are within MaxDist the
	// duplicate is counted as optical and tagged dt:Z:SQ, otherwise dt:Z:LB.
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
	// Threads sets the BGZF compression worker count from `-@/--threads`. A
	// value above 1 compresses the marked BAM output in parallel; the decoded
	// records are byte-identical to the single-threaded path.
	Threads int
	// Supp mirrors upstream's `-S`: when set, supplementary, secondary and
	// mate-unmapped non-primary alignments whose primary was flagged as a
	// duplicate are themselves flagged (and counted under "DUPLICATE NON
	// PRIMARY"). When clear — the upstream default — non-primary records are
	// never marked, even if a same-qname primary is a duplicate.
	Supp bool
	// Stats requests the summary statistics report (upstream `-s`). When
	// StatsFile is empty the report goes to StatsWriter (typically stderr);
	// otherwise it is written to StatsFile.
	Stats bool
	// StatsFile, when non-empty, writes the stats report to the named file
	// (upstream `-f FILE`, which also implies Stats).
	StatsFile string
	// StatsWriter receives the text stats report when Stats is set and
	// StatsFile is empty. A nil value silences the report.
	StatsWriter io.Writer
	// NoPG suppresses @PG injection. v1 never injects @PG so this is a
	// no-op kept for flag-compat.
	NoPG bool
}

// MarkdupResult summarises a Markdup run. The fields mirror the counters
// upstream's bam_markdup.c prints under `-s`.
type MarkdupResult struct {
	// Reading is the total number of records read from the input.
	Reading int
	// Writing is the total number of records emitted to the writer.
	Writing int
	// Excluded is the number of records skipped from scoring because they hit
	// the exclude mask (secondary, supplementary, unmapped or qcfail).
	Excluded int
	// Examined is the number of primary records processed (after include /
	// exclude flag filtering).
	Examined int
	// Paired is the number of primary records that were paired and had a
	// usable MC tag.
	Paired int
	// Single is the number of singletons keyed (mate unmapped or unpaired).
	Single int
	// DuplicatePair is the number of paired primary records flagged duplicate.
	DuplicatePair int
	// DuplicateSingle is the number of single primary records flagged duplicate.
	DuplicateSingle int
	// DuplicatePairOptical / DuplicateSingleOptical count the optical subset of
	// the paired / single duplicates (nonzero only when MaxDist > 0).
	DuplicatePairOptical   int
	DuplicateSingleOptical int
	// DuplicateNonPrimary is the number of non-primary (supplementary /
	// secondary) records flagged duplicate (nonzero only with Supp).
	DuplicateNonPrimary int
	// DuplicateNonPrimaryOptical is the optical subset of DuplicateNonPrimary.
	DuplicateNonPrimaryOptical int
	// Duplicates is the total number of records flagged as duplicates (0x400
	// added), kept for backward compatibility with existing callers/tests.
	Duplicates int
	// Written is the number of records emitted to the writer (alias of
	// Writing, retained for backward compatibility).
	Written int
}

// ReaderOpener returns a fresh sam.Reader-able io.Reader for each pass. Used
// to support two-pass streaming without seek.
type ReaderOpener func() (io.ReadCloser, error)

// dupInfo records, for a qname flagged as a duplicate, the qname of the
// chosen original and whether the duplicate was keyed as a pair or a single.
type dupInfo struct {
	dupOf  string
	paired bool
	// suppLink mirrors the chosen-duplicate's suppLink: only when set is the
	// qname eligible for non-primary (-S) duplicate marking, matching
	// upstream's add_duplicate gate.
	suppLink bool
	// optical is set during pass 2 once the read names have been compared.
	optical bool
}

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

	// ---- Pass 1: collect keys (coordinate-windowed, online resolution) ---
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
	//
	// Memory: upstream keeps only a MOVING WINDOW of reads — once a buffered
	// read's coordinate + max_length is behind the current read's coordinate
	// (same tid), it can never match another read and is purged from both
	// hashes. We reproduce that here: each hash slot holds only the CURRENT
	// winner (resolution is fully online), and slots whose coordinate has
	// fallen behind the window are evicted, so peak memory is O(active
	// window) rather than O(whole file). The set of records flagged duplicate
	// is identical to the accumulate-then-resolve approach because a beaten
	// winner is re-marked in primaryDup the moment it loses.
	pairSlots := make(map[markdupKey]*markdupSlot)
	singleSlots := make(map[markdupKey]*markdupSlot)
	primaryDup := make(map[string]*dupInfo) // qname -> dup classification

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

	// Sweep window state. curTid/curCoord track the current read's reference
	// and (raw, unclipped-single-key) coordinate; slots on curTid whose coord
	// falls more than MaxLen behind curCoord are purged. maxLen mirrors
	// upstream's param->max_length (-l, default 300).
	maxLen := int64(opts.MaxLen)
	curTid := int32(-2)
	curCoord := int64(0)

	// evictBehind removes hash slots that can no longer match any future read.
	// A slot on the current reference is dead once its coordinate is more than
	// maxLen behind the sweep front; every slot on a previous reference is
	// dead (a new reference resets the coordinate frame). Marking is already
	// resolved online, so purging a slot loses no information.
	evictBehind := func() {
		for k, s := range pairSlots {
			if s.tid != curTid || s.coord+maxLen <= curCoord {
				delete(pairSlots, k)
			}
		}
		for k, s := range singleSlots {
			if s.tid != curTid || s.coord+maxLen <= curCoord {
				delete(singleSlots, k)
			}
		}
	}

	// markPairDup records e as a pair duplicate of dupOf. When e was itself a
	// former winner already recorded as an original, its accumulated supp link
	// is OR-ed in (upstream registers a qname for -S marking when ANY of its
	// reads carries SA/XA or an unmapped mate).
	markPairDup := func(e *markdupEntry, dupOf string) {
		if d, already := primaryDup[e.qname]; already && d.paired {
			d.suppLink = d.suppLink || e.suppLink
			return
		}
		primaryDup[e.qname] = &dupInfo{dupOf: dupOf, paired: true, suppLink: e.suppLink}
	}
	markSingleDup := func(e *markdupEntry, dupOf string) {
		primaryDup[e.qname] = &dupInfo{dupOf: dupOf, paired: false, suppLink: e.suppLink}
	}

	for {
		rec, err := br1.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			_ = rc1.Close()
			return res, fmt.Errorf("markdup pass 1: %w", err)
		}
		res.Reading++
		// Upstream's exclude mask is SECONDARY|SUPPLEMENTARY|UNMAP|QCFAIL
		// (bam_markdup.c:1692-1696). Records hitting it — including unmapped
		// reads — are counted as EXCLUDED and skipped; everything else,
		// including records that already carry the duplicate flag, is EXAMINED
		// and (re-)scored.
		if !primaryEligible(rec, opts) {
			res.Excluded++
			continue
		}
		res.Examined++
		entry := &markdupEntry{
			qname:    rec.QName,
			flag:     rec.Flag,
			score:    calcScore(rec),
			paired:   hasMate(rec),
			suppLink: hasSuppLink(rec),
		}

		// Advance the sweep front and purge slots that fell behind the window.
		// The coordinate frame is the single-key's this_coord (the same value
		// upstream stores in in_read->pos for the window test). A new reference
		// resets the frame and drops every prior-reference slot.
		sKey := singleKey(rec, hdr)
		recTid := sKey.ThisRef
		recCoord := sKey.ThisCoord
		if recTid != curTid {
			curTid = recTid
			curCoord = recCoord
		} else if recCoord > curCoord {
			curCoord = recCoord
		}
		evictBehind()

		// single_hash bookkeeping (every record gets a single-key).
		if slot, ok := singleSlots[sKey]; ok {
			cur := slot.entry
			// Pairing wins. If the incoming is paired and the slot's occupant is
			// not (or vice versa), the unpaired side is marked as a duplicate.
			if entry.paired != cur.paired {
				if entry.paired {
					// Incoming pair displaces the singleton occupant.
					if !cur.paired {
						markSingleDup(cur, entry.qname)
					}
					slot.entry = entry
				} else {
					// Slot is paired, incoming is singleton — mark incoming.
					markSingleDup(entry, cur.qname)
				}
			} else if !entry.paired {
				// Two singletons at the same coord: keep highest score. Upstream's
				// single-hash swap is `if (new_score > old_score)` (bam_markdup.c
				// ~1901) — a STRICT comparison with NO qname tie-break — so on an
				// exact score tie the incumbent (the record encountered first in
				// coordinate order) stays and the new arrival is the duplicate.
				if betterSingle(entry, cur) {
					markSingleDup(cur, entry.qname)
					slot.entry = entry
				} else {
					markSingleDup(entry, cur.qname)
				}
			}
			// Two paired records at the same single-key: resolved in pair_hash.
			slot.tid, slot.coord = recTid, recCoord
		} else {
			singleSlots[sKey] = &markdupSlot{entry: entry, tid: recTid, coord: recCoord}
		}

		// pair_hash bookkeeping (paired records only), resolved ONLINE against
		// the current winner so no whole-file bucket list is retained.
		if entry.paired {
			res.Paired++
			pKey, single := buildKey(rec, opts.Mode, hdr)
			if !single {
				if slot, ok := pairSlots[pKey]; ok {
					// Collision: the higher score (pair tie-break: smaller qname)
					// stays as the winner; the loser is flagged a pair duplicate.
					if betterPair(entry, slot.entry) {
						markPairDup(slot.entry, entry.qname)
						slot.entry = entry
					} else {
						markPairDup(entry, slot.entry.qname)
					}
					slot.tid, slot.coord = recTid, recCoord
				} else {
					pairSlots[pKey] = &markdupSlot{entry: entry, tid: recTid, coord: recCoord}
				}
			}
		} else {
			res.Single++
		}
	}
	_ = rc1.Close()

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
	bw := sam.NewBAMWriterThreads(out, opts.Threads)
	// closeBW drains the parallel BGZF back end's worker goroutines on any
	// early error return. Close is idempotent, so the success path's explicit
	// Close still reports the flush/EOF error.
	closed := false
	closeBW := func() {
		if !closed {
			_ = bw.Close()
		}
	}
	if err := bw.WriteHeader(br2.Header()); err != nil {
		closeBW()
		return res, err
	}
	// opticalCache memoises the per-qname-pair optical classification so the
	// (potentially repeated) read-name parsing runs at most once per pair.
	opticalCache := make(map[string]bool)
	isOptical := func(d *dupInfo, dupQName string) bool {
		if opts.MaxDist <= 0 {
			return false
		}
		key := d.dupOf + "\x00" + dupQName
		if v, ok := opticalCache[key]; ok {
			return v
		}
		v := opticalDuplicate(d.dupOf, dupQName, opts.MaxDist)
		opticalCache[key] = v
		return v
	}

	for {
		rec, err := br2.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			closeBW()
			return res, fmt.Errorf("markdup pass 2: %w", err)
		}
		if opts.ClearTags {
			rec.Aux = stripAuxTags(rec.Aux, "do", "dt", "mc")
		}
		nonPrimary := rec.Flag&(sam.FlagSecondary|sam.FlagSupplementary) != 0 || rec.IsUnmapped()
		if d, ok := primaryDup[rec.QName]; ok && !rec.IsUnmapped() {
			if !nonPrimary {
				// Primary duplicate: always marked. Count per record so a
				// duplicate pair contributes 2 to DUPLICATE PAIR, matching
				// upstream's per-pair-key accounting.
				rec.Flag |= sam.FlagDuplicate
				optical := isOptical(d, rec.QName)
				if d.paired {
					res.DuplicatePair++
					if optical {
						res.DuplicatePairOptical++
					}
				} else {
					res.DuplicateSingle++
					if optical {
						res.DuplicateSingleOptical++
					}
				}
				// Upstream mark_duplicates writes 'do' before 'dt'.
				if opts.AddTag {
					setAuxStringMD(rec, "do", d.dupOf)
				}
				if opts.MaxDist > 0 {
					if optical {
						setAuxStringMD(rec, "dt", "SQ")
					} else {
						setAuxStringMD(rec, "dt", "LB")
					}
				}
			} else if opts.Supp && d.suppLink {
				// Non-primary record whose primary is a duplicate: only marked
				// under -S, and only when the chosen duplicate primary carried
				// an SA/XA tag or an unmapped mate (upstream's add_duplicate
				// gate populates the dup_hash only in that case).
				rec.Flag |= sam.FlagDuplicate
				res.DuplicateNonPrimary++
				optical := isOptical(d, rec.QName)
				// Upstream's supplementary pass writes 'do' before 'dt'.
				if opts.AddTag {
					setAuxStringMD(rec, "do", d.dupOf)
				}
				if opts.MaxDist > 0 {
					if optical {
						setAuxStringMD(rec, "dt", "SQ")
						res.DuplicateNonPrimaryOptical++
					} else {
						setAuxStringMD(rec, "dt", "LB")
					}
				}
			}
		}
		if opts.RemoveDups && rec.Flag&sam.FlagDuplicate != 0 {
			continue
		}
		if err := bw.Write(rec); err != nil {
			closeBW()
			return res, err
		}
		res.Writing++
	}
	res.Written = res.Writing
	res.Duplicates = res.DuplicatePair + res.DuplicateSingle + res.DuplicateNonPrimary
	if werr := writeMarkdupStats(opts, &res); werr != nil {
		closeBW()
		return res, werr
	}
	closed = true
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

// markdupSlot is the current winner held in a coordinate-windowed hash. It
// carries the winning entry plus the reference id and coordinate used to purge
// the slot once the sweep front has moved more than MaxLen past it.
type markdupSlot struct {
	entry *markdupEntry
	tid   int32
	coord int64
}

type markdupEntry struct {
	qname  string
	flag   uint16
	score  int64
	paired bool
	// suppLink records whether this record carries the links upstream uses to
	// gate supplementary-duplicate marking (an SA or XA aux tag, or a flagged
	// unmapped mate). Only when the *chosen duplicate* primary has suppLink set
	// does upstream add the qname to its dup_hash, so a non-primary alignment
	// is marked under -S only in that case.
	suppLink bool
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

// betterPair reports whether a should replace incumbent b in a PAIR bucket.
// It mirrors upstream's pair-hash swap (bam_markdup.c ~1795): the higher score
// wins, and an exact score tie is broken by lexicographically smaller qname
// (tie_add = +1 when qname(a) < qname(b), so a swaps in). a is the newer
// arrival, b the incumbent.
func betterPair(a, b *markdupEntry) bool {
	if a.score != b.score {
		return a.score > b.score
	}
	return a.qname < b.qname
}

// betterSingle reports whether a should replace incumbent b in the SINGLE hash.
// Upstream's single-hash swap is a STRICT `new_score > old_score`
// (bam_markdup.c ~1901) with NO qname tie-break, so an exact score tie leaves
// the incumbent in place and marks the newcomer as the duplicate.
func betterSingle(a, b *markdupEntry) bool {
	return a.score > b.score
}

const (
	mdOrientFF = 0
	mdOrientFR = 1
	mdOrientRF = 2
	mdOrientRR = 3
	mdLeftLE   = 0
	mdLeftRI   = 1
)

// primaryEligible reports whether the record is EXAMINED (eligible for
// duplicate scoring) rather than EXCLUDED. It mirrors upstream's exclude mask
// SECONDARY|SUPPLEMENTARY|UNMAP|QCFAIL (bam_markdup.c:1692-1699): unmapped
// reads are excluded, but records that merely already carry the duplicate
// flag are NOT — they are re-examined (and, under -c, would have had FDUP
// cleared first).
func primaryEligible(rec *sam.Record, opts MarkdupOptions) bool {
	if opts.IncludeFlags != 0 && rec.Flag&opts.IncludeFlags != opts.IncludeFlags {
		return false
	}
	if opts.ExcludeFlags != 0 && rec.Flag&opts.ExcludeFlags != 0 {
		return false
	}
	if rec.Flag&(sam.FlagSecondary|sam.FlagSupplementary|sam.FlagUnmapped|sam.FlagQCFail) != 0 {
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
	if hasMate(rec) {
		if a, ok := rec.GetAux("ms"); ok {
			if v, ok := a.Int(); ok {
				s += v
			}
		}
	}
	return s
}

// hasMate mirrors upstream bam_markdup.c has_mate: a record counts as paired
// (participates in the pair hash) only when it is PAIRED, its mate is NOT
// flagged unmapped, and its mate has a real coordinate — i.e. NOT the
// (mtid==-1, mpos==-1) sentinel a fixmate singleton fixup leaves behind. In
// our record model mtid==-1 is RNext "*"/"" and mpos==-1 is PNext 0.
func hasMate(rec *sam.Record) bool {
	if !rec.IsPaired() || rec.IsMateUnmapped() {
		return false
	}
	mateNoRef := rec.RNext == "" || rec.RNext == "*"
	if mateNoRef && rec.PNext == 0 {
		return false
	}
	return true
}

// buildKey returns the bucket key for rec. The second return is true when
// the record is keyed as a singleton.
func buildKey(rec *sam.Record, mode MarkdupMode, hdr *sam.Header) (markdupKey, bool) {
	thisRef := int32(hdr.RefIndex(rec.RName)) + 1 // +1 so zero never appears
	thisCoord := unclippedStart(rec)
	thisEnd := unclippedEnd(rec)
	rev := rec.Flag&sam.FlagReverse != 0

	mateMapped := hasMate(rec)
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
func otherUnclippedStart(matePos int64, mateCigar string) int64 {
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
func otherUnclippedEnd(matePos int64, mateCigar string) int64 {
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

// hasSuppLink reports whether rec carries an SA or XA aux tag or has an
// unmapped mate. Upstream's add_duplicate only registers a duplicate's qname
// in the dup_hash (used to mark non-primary alignments under -S) when one of
// these holds, so only such duplicates can propagate to their supplementary /
// secondary alignments.
func hasSuppLink(rec *sam.Record) bool {
	if rec.IsPaired() && rec.IsMateUnmapped() {
		return true
	}
	if _, ok := rec.GetAux("SA"); ok {
		return true
	}
	if _, ok := rec.GetAux("XA"); ok {
		return true
	}
	return false
}

// ----- Optical-duplicate detection (-d) -------------------------------------

// markdupCoords holds the parsed Illumina read-name coordinates: the x/y
// tile coordinates and the [0,tEnd) prefix (machine/run/flowcell/lane/tile)
// that must match between two read names for them to be on the same tile.
type markdupCoords struct {
	x, y int64
	tEnd int
}

// parseMarkdupCoords reproduces upstream bam_markdup.c get_coordinates_colons:
// it counts the colon separators in qname and, for the recognised Illumina
// layouts (3, 4, 6 or 7 colons), extracts the x and y coordinates and the
// prefix length used as the same-tile test. It returns ok=false when the name
// does not match a known layout or a coordinate cannot be parsed.
func parseMarkdupCoords(qname string) (markdupCoords, bool) {
	sep := 0
	xpos, ypos := 0, 0
	for pos := 0; pos < len(qname); pos++ {
		if qname[pos] != ':' {
			continue
		}
		sep++
		switch sep {
		case 2:
			xpos = pos + 1
		case 3:
			ypos = pos + 1
		case 4: // HiSeq style names
			xpos = ypos
			ypos = pos + 1
		case 5: // Newer Illumina format
			xpos = pos + 1
		case 6:
			ypos = pos + 1
		}
	}
	if !(sep == 3 || sep == 4 || sep == 6 || sep == 7) {
		return markdupCoords{}, false
	}
	x, okx := parseLeadingInt(qname[xpos:])
	if !okx {
		return markdupCoords{}, false
	}
	y, oky := parseLeadingInt(qname[ypos:])
	if !oky {
		return markdupCoords{}, false
	}
	return markdupCoords{x: x, y: y, tEnd: xpos}, true
}

// parseLeadingInt mirrors C's strtol(s, &end, 10) used for the x/y
// coordinates: it parses an optional sign and leading decimal digits and
// reports ok=false when no digits were consumed (end == s).
func parseLeadingInt(s string) (int64, bool) {
	i := 0
	neg := false
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		neg = s[i] == '-'
		i++
	}
	start := i
	var v int64
	overflow := false
	const maxInt64 = int64(^uint64(0) >> 1)
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		d := int64(s[i] - '0')
		// Saturate at int64 max on overflow, matching C strtol clamping to
		// LONG_MAX (so an out-of-range coordinate stays "far away" rather
		// than wrapping to a small value).
		if overflow || v > (maxInt64-d)/10 {
			overflow = true
			v = maxInt64
		} else {
			v = v*10 + d
		}
		i++
	}
	if i == start {
		return 0, false
	}
	if neg {
		v = -v
	}
	return v, true
}

// opticalDuplicate reports whether dup is an optical duplicate of ori: the
// two read names must share the same tile prefix and their x and y tile
// coordinates must each be within maxDist. Mirrors upstream's
// is_optical_duplicate.
func opticalDuplicate(ori, dup string, maxDist int) bool {
	oc, ok := parseMarkdupCoords(ori)
	if !ok {
		return false
	}
	dc, ok := parseMarkdupCoords(dup)
	if !ok {
		return false
	}
	if oc.tEnd != dc.tEnd || ori[:oc.tEnd] != dup[:dc.tEnd] {
		return false
	}
	xdiff := oc.x - dc.x
	if xdiff < 0 {
		xdiff = -xdiff
	}
	if xdiff > int64(maxDist) {
		return false
	}
	ydiff := oc.y - dc.y
	if ydiff < 0 {
		ydiff = -ydiff
	}
	return ydiff <= int64(maxDist)
}

// ----- Stats report (-s / -f) -----------------------------------------------

// writeMarkdupStats emits the upstream-format duplicate statistics summary
// when requested. The report goes to opts.StatsFile when set, otherwise to
// opts.StatsWriter. The leading "COMMAND:" line that upstream prints is
// intentionally omitted: it echoes the verbatim argv, which the library layer
// does not have. The CLI runner prepends an equivalent line.
func writeMarkdupStats(opts MarkdupOptions, res *MarkdupResult) error {
	if !opts.Stats && opts.StatsFile == "" {
		return nil
	}
	var w io.Writer
	if opts.StatsFile != "" {
		f, err := os.Create(opts.StatsFile)
		if err != nil {
			return fmt.Errorf("markdup: cannot write stats to %s: %w", opts.StatsFile, err)
		}
		defer f.Close()
		w = f
	} else {
		if opts.StatsWriter == nil {
			return nil
		}
		w = opts.StatsWriter
	}
	FormatMarkdupStats(w, res)
	return nil
}

// FormatMarkdupStats writes the duplicate-statistics block (everything from
// the READ line onward) to w, matching upstream bam_markdup.c write_stats.
func FormatMarkdupStats(w io.Writer, res *MarkdupResult) {
	els := estimateLibrarySize(res.Paired, res.DuplicatePair, res.DuplicatePairOptical)
	fmt.Fprintf(w,
		"READ: %d\n"+
			"WRITTEN: %d\n"+
			"EXCLUDED: %d\n"+
			"EXAMINED: %d\n"+
			"PAIRED: %d\n"+
			"SINGLE: %d\n"+
			"DUPLICATE PAIR: %d\n"+
			"DUPLICATE SINGLE: %d\n"+
			"DUPLICATE PAIR OPTICAL: %d\n"+
			"DUPLICATE SINGLE OPTICAL: %d\n"+
			"DUPLICATE NON PRIMARY: %d\n"+
			"DUPLICATE NON PRIMARY OPTICAL: %d\n"+
			"DUPLICATE PRIMARY TOTAL: %d\n"+
			"DUPLICATE TOTAL: %d\n"+
			"ESTIMATED_LIBRARY_SIZE: %d\n",
		res.Reading, res.Writing, res.Excluded, res.Examined, res.Paired, res.Single,
		res.DuplicatePair, res.DuplicateSingle, res.DuplicatePairOptical, res.DuplicateSingleOptical,
		res.DuplicateNonPrimary, res.DuplicateNonPrimaryOptical,
		res.DuplicateSingle+res.DuplicatePair,
		res.DuplicateSingle+res.DuplicatePair+res.DuplicateNonPrimary, els)
}

// coverageEquation is the rearranged Lander/Waterman coverage equation solved
// by estimateLibrarySize. Mirrors upstream bam_markdup.c coverage_equation.
func coverageEquation(x, c, n float64) float64 {
	return c/x - 1 + math.Exp(-n/x)
}

// estimateLibrarySize ports upstream bam_markdup.c estimate_library_size: a
// bisection solve of the coverage equation over the paired-read counts. It
// returns 0 when the inputs cannot yield an estimate (matching upstream's
// guarded early returns; the accompanying warnings are not reproduced).
func estimateLibrarySize(pairedReads, pairedDuplicateReads, optical int) int {
	nonOpticalPairs := (pairedReads - optical) / 2
	uniquePairs := (pairedReads - pairedDuplicateReads) / 2
	duplicatePairs := (pairedDuplicateReads - optical) / 2

	if !(nonOpticalPairs != 0 && duplicatePairs != 0 && uniquePairs != 0 && nonOpticalPairs > duplicatePairs) {
		return 0
	}
	m := 1.0
	M := 100.0
	if coverageEquation(m*float64(uniquePairs), float64(uniquePairs), float64(nonOpticalPairs)) < 0 {
		return 0
	}
	for coverageEquation(M*float64(uniquePairs), float64(uniquePairs), float64(nonOpticalPairs)) > 0 {
		M *= 10
	}
	for i := 0; i < 40; i++ {
		r := (m + M) / 2
		u := coverageEquation(r*float64(uniquePairs), float64(uniquePairs), float64(nonOpticalPairs))
		if u > 0 {
			m = r
		} else if u < 0 {
			M = r
		} else {
			break
		}
	}
	return int(float64(uniquePairs) * (m + M) / 2)
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
