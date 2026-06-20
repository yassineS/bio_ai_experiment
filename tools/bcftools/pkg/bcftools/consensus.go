package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// HaplotypeSelector identifies which allele the consensus engine should
// pull from each record's FORMAT/GT field. The upstream codes are
// case-insensitive; we accept them all up-front.
//
//	HapAuto    — apply every ALT (no sample restriction). Upstream default.
//	HapIndex   — apply the N-th allele (1-based). HapIndexValue is N.
//	HapRef     — REF allele in het genotypes.
//	HapAlt     — ALT allele in het genotypes.
//	HapIUPAC   — IUPAC ambiguity code for all genotypes.
//	HapLongRef — longer allele, breaking ties with REF.
//	HapLongAlt — longer allele, breaking ties with ALT.
//	HapShortRef — shorter allele, breaking ties with REF.
//	HapShortAlt — shorter allele, breaking ties with ALT.
//	HapPhasedIUPAC — the N-th allele for phased genotypes, IUPAC code for
//	  unphased ones (upstream's "NpIu" form). HapIndexValue is N.
type HaplotypeSelector int

// HaplotypeSelector enumeration; see the type doc for the upstream code
// mapping.
const (
	HapAuto HaplotypeSelector = iota
	HapIndex
	HapRef
	HapAlt
	HapIUPAC
	HapLongRef
	HapLongAlt
	HapShortRef
	HapShortAlt
	HapPhasedIUPAC
)

// MarkCase encodes the "uc" / "lc" / single-char highlight modes used by
// --mark-ins, --mark-del, --mark-snv, and --mask-with. The zero value is
// MarkNone (no highlight).
type MarkCase int

// MarkCase enumeration.
const (
	MarkNone  MarkCase = iota
	MarkUpper          // "uc"
	MarkLower          // "lc"
	MarkChar           // single character (Char field on MarkSpec)
)

// MarkSpec is the parsed form of a `--mark-*` flag.
type MarkSpec struct {
	Mode MarkCase
	Char byte
}

// ParseMarkSpec parses upstream's `uc|lc|CHAR` syntax shared by
// --mark-ins, --mark-snv, and --mask-with. --mark-del only accepts a
// single character; use the Mode==MarkChar branch.
func ParseMarkSpec(s string) (MarkSpec, error) {
	switch strings.ToLower(s) {
	case "":
		return MarkSpec{Mode: MarkNone}, nil
	case "uc":
		return MarkSpec{Mode: MarkUpper}, nil
	case "lc":
		return MarkSpec{Mode: MarkLower}, nil
	}
	if len(s) == 1 {
		return MarkSpec{Mode: MarkChar, Char: s[0]}, nil
	}
	return MarkSpec{}, fmt.Errorf("expected uc|lc|<single-char>, got %q", s)
}

// ParseHaplotypeSelector parses upstream's -H/--haplotype argument. The
// codes are case-insensitive (per the upstream usage block).
//
//	N (1,2,...)  -> HapIndex with HapIndexValue=N
//	R            -> HapRef
//	A            -> HapAlt
//	I            -> HapIUPAC
//	LR / LA      -> HapLongRef / HapLongAlt
//	SR / SA      -> HapShortRef / HapShortAlt
func ParseHaplotypeSelector(s string) (HaplotypeSelector, int, error) {
	switch strings.ToUpper(s) {
	case "":
		return HapAuto, 0, nil
	case "R":
		return HapRef, 0, nil
	case "A":
		return HapAlt, 0, nil
	case "I":
		return HapIUPAC, 0, nil
	case "LR", "L":
		// Upstream consensus.c:1312 sets PICK_LONG|PICK_REF for
		// both "L" and "LR".
		return HapLongRef, 0, nil
	case "LA":
		return HapLongAlt, 0, nil
	case "SR", "S":
		// Upstream consensus.c:1313 sets PICK_SHORT|PICK_REF for
		// both "S" and "SR".
		return HapShortRef, 0, nil
	case "SA":
		return HapShortAlt, 0, nil
	}
	// The "NpIu" form: index N for phased genotypes, IUPAC for unphased.
	// Upstream parses the leading integer then accepts a trailing "pIu"
	// (case-insensitive). "1pIu" / "2pIu" are the documented examples.
	digits := s
	suffix := false
	if i := indexNonDigit(s); i >= 0 {
		digits = s[:i]
		if !strings.EqualFold(s[i:], "pIu") {
			return 0, 0, fmt.Errorf("could not parse haplotype %q, expected a number optionally followed by \"pIu\"", s)
		}
		suffix = true
	}
	n, err := strconv.Atoi(digits)
	if err != nil {
		return 0, 0, fmt.Errorf("unrecognised haplotype selector %q", s)
	}
	if n < 1 {
		return 0, 0, fmt.Errorf("haplotype index %d must be >= 1", n)
	}
	if suffix {
		return HapPhasedIUPAC, n, nil
	}
	return HapIndex, n, nil
}

// indexNonDigit returns the index of the first non-ASCII-digit byte in s,
// or -1 if s is all digits (or empty).
func indexNonDigit(s string) int {
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			return i
		}
	}
	return -1
}

// ConsensusOptions controls the behaviour of Consensus / ConsensusFile.
//
// The v1 port covers the common-case "apply SNPs and simple indels from
// a VCF to a reference FASTA, optionally restricted to one sample's
// genotype" pipeline, plus the --chain liftover file (ChainFile) and the
// --mask BED replacement (Mask / MaskWith) — both byte-for-byte matched
// against upstream. Remaining upstream-only surface (BCF/regions input via
// the synced reader) is tracked in docs/PARITY_ROADMAP.md.
type ConsensusOptions struct {
	// Reference holds the parsed FASTA records (in input order).
	Reference []*fasta.Record

	// Sample limits the consensus to a single sample's GT field. The
	// empty string means "apply every ALT" (upstream's no-samples mode).
	Sample string

	// Haplotype selects which allele to apply when Sample is set. The
	// zero value (HapAuto) means: pull the 1st ALT for hom-alt sites,
	// REF for hom-ref, and (deterministically) the 1st allele of het
	// sites. HapIndexValue is the 1-based ALT index for HapIndex.
	Haplotype      HaplotypeSelector
	HaplotypeIndex int

	// Mask is the parsed BED of regions to replace with MaskWith. v1
	// stores the parsed regions; the CLI is responsible for loading.
	Mask     []MaskRegion
	MaskWith MarkSpec
	// MaskBED stores the (unread) BED file path so the CLI can report it.
	MaskBED string

	// MarkIns / MarkSnv highlight inserted / substituted bases.
	// MarkDel inserts a literal character per deleted base.
	MarkIns MarkSpec
	MarkSnv MarkSpec
	MarkDel MarkSpec // only Mode==MarkChar is meaningful here

	// Missing is the upstream -M/--missing CHAR: a single character to
	// emit instead of skipping a missing genotype. The zero value (0)
	// means "skip missing GTs" (upstream default).
	Missing byte

	// Absent is the upstream -a/--absent CHAR: a single character to
	// insert at positions that don't appear in any VCF record. The zero
	// value means "leave the reference base as-is".
	Absent byte

	// Prefix is added to each output sequence name (-p/--prefix).
	Prefix string

	// IUPACCodes asks for IUPAC encoding of hetero sites (-I).
	IUPACCodes bool

	// IncludeExpr / ExcludeExpr re-use the Filter type from view.go.
	IncludeExpr string
	ExcludeExpr string

	// LineWidth controls the wrapped FASTA output width. Zero means use
	// the upstream default of 60.
	LineWidth int

	// ChainFile, when non-empty, is the path of a liftover chain file to
	// write alongside the consensus FASTA. It mirrors upstream's
	// -c/--chain. The chain maps reference coordinates to the modified
	// consensus coordinates in UCSC chain format.
	ChainFile string
}

// MaskRegion is one half-open BED-style mask range.
type MaskRegion struct {
	Chrom string
	Beg   int // 1-based inclusive
	End   int // 1-based inclusive
}

// ConsensusFile is the file-aware entry point. It opens path (transparent
// gzip / BCF), reads the FASTA referenced by fastaPath, and applies the
// variants to out as wrapped FASTA.
func ConsensusFile(path, fastaPath string, out io.Writer, opts ConsensusOptions) (int, error) {
	if fastaPath == "" {
		return 0, fmt.Errorf("bcftools consensus: -f/--fasta-ref is required")
	}
	fa, err := readFasta(fastaPath)
	if err != nil {
		return 0, fmt.Errorf("bcftools consensus: read fasta: %w", err)
	}
	opts.Reference = fa
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools consensus: open %s: %w", path, err)
	}
	defer r.Close()
	return Consensus(r, out, opts)
}

// Consensus streams the VCF in in and writes the modified FASTA to out.
// The reference sequences in opts.Reference must already be populated.
func Consensus(in io.Reader, out io.Writer, opts ConsensusOptions) (int, error) {
	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return 0, fmt.Errorf("bcftools consensus: %w", err)
	}
	include, exclude, err := compileExpressions(ViewOptions{
		IncludeExpr: opts.IncludeExpr,
		ExcludeExpr: opts.ExcludeExpr,
	}, hdr)
	if err != nil {
		return 0, fmt.Errorf("bcftools consensus: %w", err)
	}

	sampleIdx := -1
	if opts.Sample != "" {
		for i, s := range hdr.Samples {
			if s == opts.Sample {
				sampleIdx = i
				break
			}
		}
		if sampleIdx < 0 {
			return 0, fmt.Errorf("bcftools consensus: sample %q not found in header", opts.Sample)
		}
	}

	// Determine upstream's iupac_GTs mode (consensus.c init_data): when the
	// VCF carries samples and neither -H (a haplotype/allele pick) nor an
	// explicit allele selector is given, the consensus applies IUPAC ambiguity
	// codes derived from FORMAT/GT — across the chosen sample (-s) or across
	// ALL samples when -s is absent. Sites-only VCFs keep the "apply 1st ALT"
	// behaviour. iupacSamples is the list of sample indices to OR together;
	// nil means the iupac_GTs mode is off.
	var iupacSamples []int
	if len(hdr.Samples) > 0 && opts.Haplotype == HapAuto && !opts.IUPACCodes {
		if sampleIdx >= 0 {
			iupacSamples = []int{sampleIdx}
		} else {
			iupacSamples = make([]int, len(hdr.Samples))
			for i := range hdr.Samples {
				iupacSamples[i] = i
			}
		}
	}

	// Group variants by chromosome (in stable input order). Use a map
	// keyed by chrom to allow us to skip whole records that don't match
	// any reference sequence.
	byChrom := make(map[string][]*vcf.Variant)
	for _, v := range variants {
		if include != nil && !include.Eval(v) {
			continue
		}
		if exclude != nil && exclude.Eval(v) {
			continue
		}
		byChrom[v.Chrom] = append(byChrom[v.Chrom], v)
	}

	lineWidth := opts.LineWidth
	if lineWidth <= 0 {
		lineWidth = 60
	}

	// Optional liftover chain output. The chain writer is shared across
	// all sequences so the running chain identifier auto-increments.
	var cw *chainWriter
	if opts.ChainFile != "" {
		f, err := openChainFile(opts.ChainFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools consensus: %w", err)
		}
		cw = &chainWriter{bw: bufio.NewWriter(f), closer: f}
		defer func() {
			cw.flush()
			cw.closer.Close()
		}()
	}

	bw := bufio.NewWriter(out)
	applied := 0
	for _, rec := range opts.Reference {
		seq := append([]byte(nil), rec.Sequence...)
		seq = applyMask(seq, rec.ID, opts.Mask, opts.MaskWith)
		if opts.Absent != 0 {
			// Upstream's -a CHAR replaces positions NOT covered by any
			// variant. v1 implementation: mark every position not later
			// rewritten as opts.Absent. We approximate by zeroing the
			// sequence first and writing back ref/alt as we go; this is
			// the standard "absent" semantics used by upstream when the
			// VCF has at most one variant per position.
			for i := range seq {
				seq[i] = opts.Absent
			}
		}

		vars := byChrom[rec.ID]
		// Apply in left-to-right order. Indels shift downstream positions,
		// so we maintain an offset relative to the original reference.
		sort.SliceStable(vars, func(i, j int) bool { return vars[i].Pos < vars[j].Pos })
		offset := 0
		// frzPos mirrors upstream consensus.c's args->fa_frz_pos: the last
		// 0-based reference position consumed by an applied variant
		// (rec->pos + rec->rlen - 1). It starts at -1 (nothing consumed).
		// prevIsInsert mirrors args->prev_is_insert: whether the most
		// recently applied variant was a net insertion. Together they drive
		// the overlap/skip decision so the port reproduces upstream's
		// "The site CHR:POS overlaps with another variant, skipping..."
		// behaviour byte-for-byte: a record at pos <= frzPos is skipped
		// unless it is a clean insertion (trim_beg, non-zero length delta)
		// landing exactly on frzPos and not following another insertion.
		frzPos := -1
		prevIsInsert := false
		// prevBasePos mirrors upstream args->prev_base_pos: the 0-based
		// genomic position of the last reference base that a previous variant
		// occupied (its last consumed base for substitutions/deletions, or its
		// anchor base for insertions). It drives the ibeg trim below, which
		// keeps an insertion from overwriting the bases an earlier
		// same-position variant already edited (consensus.c lines 1029-1037).
		prevBasePos := -1
		// chain accumulates the liftover blocks for this sequence when
		// -c/--chain is requested. The reference origin is 0 because the
		// port emits whole contigs (no -r region windowing).
		var chain *Chain
		if cw != nil {
			chain = NewChain(0)
		}
		for _, v := range vars {
			alt, ok := selectAllele(v, sampleIdx, iupacSamples, opts)
			if !ok {
				continue
			}
			ref := v.Ref
			// Upstream returns early when the reference allele is selected
			// (ialt==0): it neither edits the sequence nor advances the
			// freeze position (unless --absent, which we honour by restoring
			// the REF bases below). A selected allele identical to REF is the
			// ref-allele case, so it must not freeze the position — otherwise
			// a following insertion at the same coordinate would be wrongly
			// skipped as overlapping.
			isRefAllele := alt == ref
			// Overlap handling mirrors apply_variant() in consensus.c. The
			// length delta of the applied event (positive=insertion,
			// negative=deletion, zero=substitution) and trim_beg (whether the
			// indel carries the shared anchor base) decide whether a record
			// landing within the frozen region is dropped.
			lenDelta := len(alt) - len(ref)
			isIndel := lenDelta != 0
			trimBeg := isIndel && len(ref) > 0 && len(alt) > 0 &&
				lowerByte(ref[0]) == lowerByte(alt[0])
			origStartCheck := v.Pos - 1
			if !isRefAllele && origStartCheck <= frzPos {
				// pos <= frz: only a clean insertion exactly on frz, not
				// following another insertion and not landing on a base a real
				// edit already changed, may proceed; everything else overlaps
				// and is skipped (upstream prints a warning here).
				overlap := origStartCheck < frzPos || !trimBeg || lenDelta == 0 || prevIsInsert
				if overlap {
					fmt.Fprintf(os.Stderr, "The site %s:%d overlaps with another variant, skipping...\n", v.Chrom, v.Pos)
					continue
				}
			}
			// 1-based VCF POS -> 0-based ref index, adjusted for prior
			// indel-induced shifts. We also use the un-shifted POS to
			// drive opts.Absent's "fill the original ref" logic when
			// absent is set: in that path we restore the REF bases.
			origStart := v.Pos - 1
			start := origStart + offset
			end := start + len(ref)
			if opts.Absent != 0 {
				// Restore the REF bases (zeroed above) before we modify.
				copy(seq[start:end], ref)
			}
			if start < 0 || end > len(seq) {
				continue
			}
			if isRefAllele {
				// Reference allele selected (upstream ialt==0): the sequence
				// is left unchanged, but the freeze position still advances to
				// the last reference base the record spans (verified against
				// upstream: a 0/0 record at a position causes a subsequent
				// same-position substitution to be dropped as overlapping,
				// while a clean insertion there is still applied). A ref
				// allele is never a net insertion, so prevIsInsert is false.
				if newFrz := origStart + len(ref) - 1; newFrz > frzPos {
					frzPos = newFrz
					prevIsInsert = false
				}
				continue
			}
			// Apply highlight casing if requested.
			alt = applyMarks(ref, alt, opts)
			emitted := len(alt)
			// ibeg is the count of leading ALT bases that coincide with the
			// reference anchor of a previous same-position variant and must be
			// preserved rather than rewritten. For a net insertion that lands
			// on bases an earlier variant already edited (e.g. a substitution
			// at the same POS turning the anchor into an IUPAC code), upstream
			// skips writing those shared leading bases so the prior edit
			// survives and only the truly inserted bases are added. For all
			// other cases ibeg is 0 and the full ALT replaces the REF run.
			ibeg := 0
			if lenDelta > 0 {
				for ibeg < len(alt) && ibeg < len(ref) &&
					lowerByte(ref[ibeg]) == lowerByte(alt[ibeg]) &&
					origStart+ibeg <= prevBasePos {
					ibeg++
				}
			}
			newSeq := make([]byte, 0, len(seq)+len(alt)-len(ref))
			newSeq = append(newSeq, seq[:start+ibeg]...)
			if opts.MarkDel.Mode == MarkChar && len(alt) < len(ref) {
				// Pad the deletion with a literal character so the
				// downstream coordinates stay aligned. The total
				// emitted length equals len(ref).
				newSeq = append(newSeq, alt[ibeg:]...)
				for i := 0; i < len(ref)-len(alt); i++ {
					newSeq = append(newSeq, opts.MarkDel.Char)
				}
				emitted = len(ref)
			} else {
				newSeq = append(newSeq, alt[ibeg:]...)
			}
			newSeq = append(newSeq, seq[end:]...)
			seq = newSeq
			// Record the indel as a chain gap before advancing the offset.
			// Upstream pushes only when the emitted length differs from the
			// reference run (len_diff != 0); mark-del padding keeps the run
			// equal and therefore emits no gap. offset is fa_mod_off here:
			// the cumulative shift *before* this variant.
			if chain != nil {
				lenDiff := emitted - len(ref)
				if lenDiff != 0 {
					if len(ref) > 0 && len(alt) > 0 && lowerByte(ref[0]) == lowerByte(alt[0]) {
						// Indels usually carry the base before the event:
						// extend the ungapped block by one base.
						chain.pushGap(origStart+1, len(ref)-1, origStart+1+offset, emitted-1)
					} else {
						chain.pushGap(origStart, len(ref), origStart+offset, emitted)
					}
				}
			}
			// Track the coordinate shift by the ACTUAL emitted run
			// vs. ref consumed. With mark-del padding, emitted==len(ref)
			// so offset is unchanged (downstream coordinates align).
			offset += emitted - len(ref)
			// Advance the freeze position to the last reference base this
			// record consumed (upstream args->fa_frz_pos = rec->pos +
			// rec->rlen - 1) and record whether it was a net insertion.
			frzPos = origStart + len(ref) - 1
			prevIsInsert = lenDelta > 0
			if lenDelta > 0 {
				// Insertion: the anchor base is the only reference position it
				// occupies (upstream prev_base_pos = rec->pos).
				prevBasePos = origStart
			} else {
				prevBasePos = origStart + len(ref) - 1
			}
			applied++
		}
		if chain != nil {
			if err := cw.writeChain(chain, rec.ID, len(rec.Sequence)); err != nil {
				return applied, err
			}
		}

		name := opts.Prefix + rec.ID
		if _, err := fmt.Fprintf(bw, ">%s\n", name); err != nil {
			return applied, err
		}
		for i := 0; i < len(seq); i += lineWidth {
			j := i + lineWidth
			if j > len(seq) {
				j = len(seq)
			}
			if _, err := bw.Write(seq[i:j]); err != nil {
				return applied, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return applied, err
			}
		}
	}
	return applied, bw.Flush()
}

// selectAllele picks the ALT (or REF) string to inject per opts. Returns
// (allele, true) when the record should be applied, or ("", false) when
// it should be skipped (missing GT, ./. with no -M, mismatched REF, ...).
func selectAllele(v *vcf.Variant, sampleIdx int, iupacSamples []int, opts ConsensusOptions) (string, bool) {
	if len(v.Alt) == 0 {
		return "", false
	}
	// iupac_GTs mode (upstream default when the VCF has samples and no -H/-a):
	// OR the alleles present in every selected sample's GT into one IUPAC code.
	if len(iupacSamples) > 0 {
		return iupacAlleleAcrossSamples(v, iupacSamples, opts)
	}
	if sampleIdx < 0 {
		// No -s flag and no samples (sites-only VCF): apply the 1st ALT
		// (upstream "all variants" default).
		return v.Alt[0], true
	}
	if sampleIdx >= len(v.Samples) {
		return "", false
	}
	gt, ok := v.Samples[sampleIdx].Data["GT"]
	if !ok {
		return "", false
	}
	parts := splitGT(gt)
	if len(parts) == 0 {
		return "", false
	}
	phased := gtIsPhased(gt)
	// Missing GT handling.
	allMissing := true
	for _, p := range parts {
		if p != "." {
			allMissing = false
			break
		}
	}
	if allMissing {
		if opts.Missing != 0 {
			return string(opts.Missing), true
		}
		return "", false
	}
	// Phased-index / unphased-IUPAC (upstream "NpIu"). Mirrors
	// apply_variant: when allele==PICK_IUPAC and a haplotype index is set,
	// the IUPAC branch is taken only for unphased genotypes; phased
	// genotypes fall through to plain haplotype-index selection.
	if opts.Haplotype == HapPhasedIUPAC {
		if phased {
			return phasedIndexAllele(parts, opts.HaplotypeIndex, v, opts)
		}
		return iupacAllele(parts, v, opts)
	}
	// IUPAC mode (-I or -H I). Upstream applies the IUPAC ambiguity code
	// over every allele present in the genotype.
	if opts.IUPACCodes || opts.Haplotype == HapIUPAC {
		return iupacAllele(parts, v, opts)
	}
	switch opts.Haplotype {
	case HapIndex:
		// -H N selects the N-th haplotype slot of the genotype (1-based),
		// regardless of phasing, then resolves it to an allele. This
		// mirrors upstream's use_hap branch (ialt = GT[haplotype-1]).
		return phasedIndexAllele(parts, opts.HaplotypeIndex, v, opts)
	case HapRef:
		if isHet(parts) {
			return v.Ref, true
		}
		return alleleString(parts[0], v), true
	case HapAlt:
		if isHet(parts) {
			return v.Alt[0], true
		}
		return alleleString(parts[0], v), true
	case HapLongRef, HapLongAlt, HapShortRef, HapShortAlt:
		// Pick longer/shorter allele, breaking ties per the suffix.
		la := alleleString(parts[0], v)
		ra := la
		if len(parts) > 1 {
			ra = alleleString(parts[1], v)
		}
		switch {
		case opts.Haplotype == HapLongRef || opts.Haplotype == HapLongAlt:
			if len(la) > len(ra) {
				return la, true
			}
			if len(ra) > len(la) {
				return ra, true
			}
			// tie -> REF or ALT
			if opts.Haplotype == HapLongRef {
				return v.Ref, true
			}
			return v.Alt[0], true
		case opts.Haplotype == HapShortRef || opts.Haplotype == HapShortAlt:
			if len(la) < len(ra) {
				return la, true
			}
			if len(ra) < len(la) {
				return ra, true
			}
			if opts.Haplotype == HapShortRef {
				return v.Ref, true
			}
			return v.Alt[0], true
		}
	}
	// HapAuto / default: take the 1st allele slot.
	return alleleString(parts[0], v), true
}

// gtIsPhased reports whether a raw GT string represents a phased genotype.
// Upstream treats a genotype as phased when its alleles carry the phase
// bit, which in VCF text corresponds to a '|' separator. A haploid call is
// trivially phased.
func gtIsPhased(gt string) bool {
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' {
			return false
		}
	}
	// No '/' separator: either a haploid call (trivially phased) or '|'
	// separators throughout (phased).
	return true
}

// phasedIndexAllele returns the allele indexed by the 1-based GT slot idx,
// mirroring upstream's use_hap branch (ialt = GT[haplotype-1]).
//
// When idx exceeds the genotype's ploidy, upstream emits the -M/--missing
// character only if the genotype's first or last slot is itself missing,
// otherwise it warns and skips the record. When the addressed slot names a
// missing allele, the missing character is emitted if -M is set, otherwise
// the record is skipped.
func phasedIndexAllele(parts []string, idx int, v *vcf.Variant, opts ConsensusOptions) (string, bool) {
	if idx < 1 {
		idx = 1
	}
	if idx > len(parts) {
		// Out-of-range haplotype index: only emit missing when the
		// genotype edges are themselves missing (upstream's check on
		// ptr[0] / ptr[fmt->n-1]); otherwise skip.
		if opts.Missing != 0 && len(parts) > 0 && (parts[0] == "." || parts[len(parts)-1] == ".") {
			return string(opts.Missing), true
		}
		return "", false
	}
	slot := parts[idx-1]
	if slot == "." || slot == "" {
		if opts.Missing != 0 {
			return string(opts.Missing), true
		}
		return "", false
	}
	return alleleString(slot, v), true
}

// iupacAlleleAcrossSamples implements upstream's iupac_GTs apply path
// (consensus.c apply_variant): it collects the allele indices present in the
// FORMAT/GT of every sample in samples (iupac_add_gt) and folds them into one
// IUPAC ambiguity code via iupac_set_allele. When no sample carries a usable
// allele the record is skipped (or the -M missing character is emitted).
func iupacAlleleAcrossSamples(v *vcf.Variant, samples []int, opts ConsensusOptions) (string, bool) {
	var parts []string
	for _, si := range samples {
		if si < 0 || si >= len(v.Samples) {
			continue
		}
		gt, ok := v.Samples[si].Data["GT"]
		if !ok {
			continue
		}
		parts = append(parts, splitGT(gt)...)
	}
	// Keep only present (non-missing) allele tokens; if none are set the
	// record is skipped unless -M provides a missing character.
	anySet := false
	for _, p := range parts {
		if p != "." && p != "" {
			anySet = true
			break
		}
	}
	if !anySet {
		if opts.Missing != 0 {
			return string(opts.Missing), true
		}
		return "", false
	}
	return iupacAllele(parts, v, opts)
}

// iupacAllele encodes the alleles present in a genotype as an IUPAC
// ambiguity string, mirroring upstream's iupac_set_allele: it OR-s the
// per-position IUPAC bitmasks across every (non-symbolic) allele in the
// genotype and emits the result on the longest allele. When no allele
// admits an IUPAC encoding it falls back to the first allele present.
func iupacAllele(parts []string, v *vcf.Variant, opts ConsensusOptions) (string, bool) {
	// Collect the distinct allele indices present in the genotype.
	seen := map[int]bool{}
	order := []int{}
	for _, p := range parts {
		if p == "." || p == "" {
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			continue
		}
		if !seen[n] {
			seen[n] = true
			order = append(order, n)
		}
	}
	if len(order) == 0 {
		if opts.Missing != 0 {
			return string(opts.Missing), true
		}
		return "", false
	}

	var fallback string
	haveFallback := false
	maxLen := 0
	var bitmask []byte
	target := ""
	targetLen := 0
	for _, n := range order {
		al := alleleIndexString(n, v)
		if !haveFallback {
			fallback = al
			haveFallback = true
		}
		mask := make([]byte, len(al))
		ok := true
		for j := 0; j < len(al); j++ {
			b := iupac2bitmask(al[j])
			if b < 0 {
				ok = false
				break
			}
			mask[j] = byte(b)
		}
		if !ok {
			continue // symbolic allele / invalid character
		}
		if len(al) > maxLen {
			grown := make([]byte, len(al))
			copy(grown, bitmask)
			bitmask = grown
			maxLen = len(al)
		}
		for j := 0; j < len(al); j++ {
			bitmask[j] |= mask[j]
		}
		// Upstream tracks the longest *ALT* allele (index > 0) as target.
		if n > 0 && len(al) > targetLen {
			targetLen = len(al)
			target = al
		}
	}
	if targetLen == 0 {
		// No usable ALT allele; emit the fallback (first allele present).
		return fallback, true
	}
	out := []byte(target)
	for j := 0; j < targetLen; j++ {
		out[j] = bitmask2iupac(bitmask[j])
	}
	return string(out), true
}

// alleleIndexString returns the allele string for a numeric GT allele
// index (0 == REF, N == ALT[N-1]).
func alleleIndexString(n int, v *vcf.Variant) string {
	if n == 0 {
		return v.Ref
	}
	if n-1 < len(v.Alt) {
		return v.Alt[n-1]
	}
	return v.Ref
}

// iupac2bitmask maps a single nucleotide or IUPAC ambiguity character to a
// 4-bit mask (A=1, C=2, G=4, T=8). It returns -1 for any character that is
// not a defined IUPAC code (symbolic alleles, gaps, etc.). The mapping is
// case-insensitive and matches htslib/bcftools' iupac2bitmask table.
func iupac2bitmask(c byte) int {
	if c >= 'a' && c <= 'z' {
		c -= 32
	}
	switch c {
	case 'A':
		return 1
	case 'C':
		return 2
	case 'G':
		return 4
	case 'T':
		return 8
	case 'M':
		return 1 | 2
	case 'R':
		return 1 | 4
	case 'W':
		return 1 | 8
	case 'S':
		return 2 | 4
	case 'Y':
		return 2 | 8
	case 'K':
		return 4 | 8
	case 'V':
		return 1 | 2 | 4
	case 'H':
		return 1 | 2 | 8
	case 'D':
		return 1 | 4 | 8
	case 'B':
		return 2 | 4 | 8
	case 'N':
		return 1 | 2 | 4 | 8
	}
	return -1
}

// bitmask2iupac is the inverse of iupac2bitmask: it maps a 4-bit nucleotide
// mask (1..15) to its IUPAC ambiguity character, returning 0 for an
// out-of-range mask. It matches htslib/bcftools' bitmask2iupac table.
func bitmask2iupac(mask byte) byte {
	table := [16]byte{'.', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}
	if mask == 0 || mask > 15 {
		return 0
	}
	return table[mask]
}

// alleleString returns the full allele string (REF or ALT[N-1]) for one
// slot of a GT field. Missing or malformed slots return v.Ref.
func alleleString(slot string, v *vcf.Variant) string {
	if slot == "." || slot == "" {
		return v.Ref
	}
	n, err := strconv.Atoi(slot)
	if err != nil {
		return v.Ref
	}
	if n == 0 {
		return v.Ref
	}
	if n-1 < len(v.Alt) {
		return v.Alt[n-1]
	}
	return v.Ref
}

// isHet returns true when the parsed GT slots represent a heterozygous
// genotype (at least two distinct allele indices).
func isHet(slots []string) bool {
	if len(slots) < 2 {
		return false
	}
	first := slots[0]
	for _, s := range slots[1:] {
		if s != first {
			return true
		}
	}
	return false
}

// iupacCode returns the IUPAC ambiguity code for two unphased bases, or
// 0 if there is no defined code.
func iupacCode(a, b byte) byte {
	if a > b {
		a, b = b, a
	}
	switch string([]byte{a, b}) {
	case "AA":
		return 'A'
	case "CC":
		return 'C'
	case "GG":
		return 'G'
	case "TT":
		return 'T'
	case "AC":
		return 'M'
	case "AG":
		return 'R'
	case "AT":
		return 'W'
	case "CG":
		return 'S'
	case "CT":
		return 'Y'
	case "GT":
		return 'K'
	}
	return 0
}

// applyMarks rewrites alt according to MarkSnv / MarkIns. A SNV (equal-
// length REF/ALT, len 1) uses MarkSnv; a longer ALT (insertion) uses
// MarkIns for the inserted suffix.
func applyMarks(ref, alt string, opts ConsensusOptions) string {
	if len(ref) == 1 && len(alt) == 1 && ref != alt {
		return applyCase(alt, opts.MarkSnv)
	}
	if len(alt) > len(ref) && strings.HasPrefix(alt, ref) {
		// Simple insertion: ref is the prefix.
		head := alt[:len(ref)]
		ins := alt[len(ref):]
		return head + applyCase(ins, opts.MarkIns)
	}
	return alt
}

// applyCase folds a string per a MarkSpec.
func applyCase(s string, spec MarkSpec) string {
	switch spec.Mode {
	case MarkUpper:
		return strings.ToUpper(s)
	case MarkLower:
		return strings.ToLower(s)
	case MarkChar:
		b := make([]byte, len(s))
		for i := range b {
			b[i] = spec.Char
		}
		return string(b)
	}
	return s
}

// applyMask rewrites bases inside any Mask region according to MaskWith.
// A MarkUpper/MarkLower spec simply changes case; MarkChar overwrites
// every base with the literal character; MarkNone defaults to 'N'.
func applyMask(seq []byte, chrom string, masks []MaskRegion, mw MarkSpec) []byte {
	if len(masks) == 0 {
		return seq
	}
	for _, m := range masks {
		if m.Chrom != chrom {
			continue
		}
		beg := m.Beg - 1
		end := m.End
		if beg < 0 {
			beg = 0
		}
		if end > len(seq) {
			end = len(seq)
		}
		if beg >= end {
			continue
		}
		switch mw.Mode {
		case MarkUpper:
			for i := beg; i < end; i++ {
				seq[i] = upperByte(seq[i])
			}
		case MarkLower:
			for i := beg; i < end; i++ {
				seq[i] = lowerByte(seq[i])
			}
		case MarkChar:
			for i := beg; i < end; i++ {
				seq[i] = mw.Char
			}
		default:
			// Upstream default is 'N'.
			for i := beg; i < end; i++ {
				seq[i] = 'N'
			}
		}
	}
	return seq
}

func upperByte(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

func lowerByte(b byte) byte {
	if b >= 'A' && b <= 'Z' {
		return b + 32
	}
	return b
}

// readFasta is a thin wrapper around fasta.Reader.ReadAll that opens the
// named file through iohelper (transparent gzip / bgzf).
func readFasta(path string) ([]*fasta.Record, error) {
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer r.Close()
	return fasta.NewReader(r).ReadAll()
}

// LoadMaskBED parses a BED file (one CHROM <tab> BEG <tab> END per line)
// into MaskRegion slices using 1-based inclusive coordinates. Lines
// starting with '#' or blank lines are skipped.
func LoadMaskBED(path string) ([]MaskRegion, error) {
	regs, err := LoadRegionsFile(path)
	if err != nil {
		return nil, err
	}
	out := make([]MaskRegion, 0, len(regs))
	for _, r := range regs {
		// LoadRegionsFile emits "chrom:beg-end" or "chrom" lines.
		colon := strings.IndexByte(r, ':')
		if colon < 0 {
			// Whole chromosome.
			out = append(out, MaskRegion{Chrom: r, Beg: 1, End: 1 << 30})
			continue
		}
		chrom := r[:colon]
		rest := r[colon+1:]
		dash := strings.IndexByte(rest, '-')
		if dash < 0 {
			n, err := strconv.Atoi(rest)
			if err != nil {
				return nil, fmt.Errorf("bcftools consensus: bad mask %q", r)
			}
			out = append(out, MaskRegion{Chrom: chrom, Beg: n, End: n})
			continue
		}
		beg, err := strconv.Atoi(rest[:dash])
		if err != nil {
			return nil, fmt.Errorf("bcftools consensus: bad mask %q", r)
		}
		end, err := strconv.Atoi(rest[dash+1:])
		if err != nil {
			return nil, fmt.Errorf("bcftools consensus: bad mask %q", r)
		}
		out = append(out, MaskRegion{Chrom: chrom, Beg: beg, End: end})
	}
	return out, nil
}
