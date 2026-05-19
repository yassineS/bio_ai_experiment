package bcftools

import (
	"bufio"
	"fmt"
	"io"
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
//	HapIUPAC   — IUPAC ambiguity code in het genotypes.
//	HapLongRef — longer allele, breaking ties with REF.
//	HapLongAlt — longer allele, breaking ties with ALT.
//	HapShortRef — shorter allele, breaking ties with REF.
//	HapShortAlt — shorter allele, breaking ties with ALT.
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
	// Try a plain integer (the "N" form). Upstream also recognises a
	// trailing "pIu" qualifier ("phased index / IUPAC for unphased")
	// which v1 rejects with a roadmap pointer in the CLI layer.
	n, err := strconv.Atoi(s)
	if err != nil {
		return 0, 0, fmt.Errorf("unrecognised haplotype selector %q", s)
	}
	if n < 1 {
		return 0, 0, fmt.Errorf("haplotype index %d must be >= 1", n)
	}
	return HapIndex, n, nil
}

// ConsensusOptions controls the behaviour of Consensus / ConsensusFile.
//
// The v1 port covers the common-case "apply SNPs and simple indels from
// a VCF to a reference FASTA, optionally restricted to one sample's
// genotype" pipeline. Upstream-only features (the --chain liftover file,
// --mask BED replacement, BCF/regions input via the synced reader) are
// tracked in docs/PARITY_ROADMAP.md and gated in the CLI.
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
	})
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
		lastEnd := 0
		for _, v := range vars {
			// Skip overlapping variants (upstream emits a warning; we
			// follow the "first wins" rule for v1 to keep coordinates
			// well-defined).
			if v.Pos <= lastEnd {
				continue
			}
			alt, ok := selectAllele(v, sampleIdx, opts)
			if !ok {
				continue
			}
			ref := v.Ref
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
			// Apply highlight casing if requested.
			alt = applyMarks(ref, alt, opts)
			newSeq := make([]byte, 0, len(seq)+len(alt)-len(ref))
			newSeq = append(newSeq, seq[:start]...)
			emitted := len(alt)
			if opts.MarkDel.Mode == MarkChar && len(alt) < len(ref) {
				// Pad the deletion with a literal character so the
				// downstream coordinates stay aligned. The total
				// emitted length equals len(ref).
				newSeq = append(newSeq, alt...)
				for i := 0; i < len(ref)-len(alt); i++ {
					newSeq = append(newSeq, opts.MarkDel.Char)
				}
				emitted = len(ref)
			} else {
				newSeq = append(newSeq, alt...)
			}
			newSeq = append(newSeq, seq[end:]...)
			seq = newSeq
			// Track the coordinate shift by the ACTUAL emitted run
			// vs. ref consumed. With mark-del padding, emitted==len(ref)
			// so offset is unchanged (downstream coordinates align).
			offset += emitted - len(ref)
			lastEnd = origStart + len(ref)
			applied++
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
func selectAllele(v *vcf.Variant, sampleIdx int, opts ConsensusOptions) (string, bool) {
	if len(v.Alt) == 0 {
		return "", false
	}
	if sampleIdx < 0 {
		// No -s flag: apply the 1st ALT (upstream "all variants" default).
		return v.Alt[0], true
	}
	if sampleIdx >= len(v.Samples) {
		return "", false
	}
	gt, ok := v.Samples[sampleIdx].Data["GT"]
	if !ok {
		return "", false
	}
	// Normalize phase.
	parts := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
	if len(parts) == 0 {
		return "", false
	}
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
	// IUPAC mode (-I or HapIUPAC). Only meaningful for biallelic SNPs;
	// fall back to ALT for other cases.
	if opts.IUPACCodes || opts.Haplotype == HapIUPAC {
		if len(v.Ref) == 1 && len(v.Alt[0]) == 1 {
			a1 := alleleBase(parts[0], v)
			a2 := a1
			if len(parts) > 1 {
				a2 = alleleBase(parts[1], v)
			}
			if iupac := iupacCode(a1, a2); iupac != 0 {
				return string([]byte{iupac}), true
			}
		}
		return v.Alt[0], true
	}
	switch opts.Haplotype {
	case HapIndex:
		// 1-based ALT index from -H N. ALT[0] corresponds to N=1; N==0
		// means REF.
		if opts.HaplotypeIndex == 0 {
			return v.Ref, true
		}
		if opts.HaplotypeIndex-1 < len(v.Alt) {
			return v.Alt[opts.HaplotypeIndex-1], true
		}
		return "", false
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

// alleleBase returns the first base of allele N for IUPAC handling.
func alleleBase(slot string, v *vcf.Variant) byte {
	s := alleleString(slot, v)
	if len(s) == 0 {
		return 'N'
	}
	return s[0]
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
