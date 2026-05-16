// bcftools csq — predict variant consequences against a GFF annotation.
//
// V1 SIMPLIFICATION: this port implements ONLY the protein-coding SNP
// consequence classifier. The full upstream csq.c performs haplotype-
// aware phasing, indel handling, splice-site and stop-gain detection,
// and compound-het bookkeeping; none of those are implemented here.
// See docs/PARITY_ROADMAP.md#bcftools for the deferred tail.
//
// What v1 DOES do:
//   - Load a GFF3 file and build a per-gene CDS index keyed by transcript ID.
//   - Load the reference FASTA into memory.
//   - For each input VCF record:
//     * Skip non-SNPs (REF length != 1 OR ALT length != 1).
//     * Locate the transcript(s) whose CDS exons cover POS.
//     * Compute the codon and amino-acid change.
//     * Emit INFO/BCSQ as "consequence|gene|transcript|biotype|strand|aa_change|dna_change"
//       (one entry per matching transcript, comma-separated).
//   - Pass non-coding sites through unchanged (no BCSQ).
//
// The CLI accepts the full upstream getopt_long surface (see csq.c at
// `static struct option loptions[]`). Flags we don't yet honour either
// no-op or hard-reject with a roadmap pointer.

package bcftools

import (
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/gff"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// CSQOptions controls CSQ / CSQFile.
type CSQOptions struct {
	// FastaRef is the reference FASTA path (-f/--fasta-ref). Required.
	FastaRef string
	// GFFAnnot is the GFF3 annotation path (-g/--gff-annot). Required.
	GFFAnnot string

	// CustomTag is the INFO tag to write under, default "BCSQ".
	CustomTag string

	// LocalCSQ matches upstream's -l/--local-csq (predict per-record
	// rather than per-haplotype). v1 always operates per-record so this
	// flag is effectively a no-op but accepted for parity.
	LocalCSQ bool

	// Phase is upstream's -p/--phase {a|m|r|R|s}. v1 stores it for
	// future use; the per-record SNP classifier does not depend on
	// haplotype phasing.
	Phase byte

	// NCSQ is upstream's -n/--ncsq cap on per-haplotype consequences.
	// v1 reports every matching transcript without a cap; the field is
	// retained for parity.
	NCSQ int
	// TrimProteinSeq is upstream's -B/--trim-protein-seq. v1 accepts
	// and stores it; no trimming is performed.
	TrimProteinSeq int
	// GeneticCode is upstream's -C/--genetic-code. v1 supports only
	// the standard table (0). Non-zero values are rejected.
	GeneticCode int
	// Verbosity is upstream's -v/--verbose / --verbosity.
	Verbosity int
	// Quiet is upstream's -q/--quiet (deprecated upstream; accepted).
	Quiet bool
	// Force is upstream's --force (skip sanity checks).
	Force bool
	// NoVersion is upstream's --no-version.
	NoVersion bool

	// Statistics is upstream's -x/--statistics path.
	Statistics string

	// IncludeExpr / ExcludeExpr are upstream's -i / -e expressions.
	// v1 stores them but does not evaluate (the SNP classifier runs on
	// every input record).
	IncludeExpr string
	ExcludeExpr string

	// Sample selection (upstream -s / -S). v1 stores them; no
	// per-sample subsetting is done in v1 (consequences are
	// position-driven, not sample-driven).
	Samples     []string
	SamplesFile string

	// Regions / Targets — post-filters.
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string

	// OutputFormat — upstream supports v/z/u/b/t. v1 emits only `v`
	// (uncompressed VCF text).
	OutputFormat OutputFormat

	// DumpGFF is upstream's --dump-gff. Accepted; v1 does not produce
	// debug dumps.
	DumpGFF string
	// UnifyChrNames is upstream's --unify-chr-names. Accepted; v1
	// honours the special "0" value (no rewriting); other specs are
	// stored but unused.
	UnifyChrNames string
}

// CSQFile streams the VCF at vcfPath through the CSQ pipeline and
// writes the annotated VCF to w. It returns the number of records
// emitted.
func CSQFile(vcfPath string, w io.Writer, opts CSQOptions) (int, error) {
	if opts.FastaRef == "" {
		return 0, fmt.Errorf("bcftools csq: -f/--fasta-ref is required")
	}
	if opts.GFFAnnot == "" {
		return 0, fmt.Errorf("bcftools csq: -g/--gff-annot is required")
	}
	idx, err := loadCSQIndex(opts.FastaRef, opts.GFFAnnot)
	if err != nil {
		return 0, err
	}
	r, err := iohelper.OpenReader(vcfPath)
	if err != nil {
		return 0, fmt.Errorf("open %q: %w", vcfPath, err)
	}
	defer r.Close()
	return CSQ(r, w, idx, opts)
}

// CSQ reads VCF records from r, annotates them with the BCSQ INFO tag,
// and writes them to w.
func CSQ(r io.Reader, w io.Writer, idx *CSQIndex, opts CSQOptions) (int, error) {
	if opts.CustomTag == "" {
		opts.CustomTag = "BCSQ"
	}

	vr := vcf.NewReader(r)
	hdr, err := vr.ReadHeader()
	if err != nil {
		return 0, fmt.Errorf("read VCF header: %w", err)
	}

	// Inject INFO meta-line for the consequence tag.
	tag := opts.CustomTag
	hdr.MetaInfo = ensureCSQInfoLine(hdr.MetaInfo, tag)

	out := vcf.NewWriter(w, hdr)
	if err := out.WriteHeader(); err != nil {
		return 0, fmt.Errorf("write VCF header: %w", err)
	}

	regionSpecs := append([]string(nil), opts.Regions...)
	regionSpecs = append(regionSpecs, opts.Targets...)

	written := 0
	for {
		v, err := vr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return written, fmt.Errorf("read variant: %w", err)
		}
		if len(regionSpecs) > 0 && !regionMatches(v, regionSpecs) {
			continue
		}

		// Annotate per-transcript consequences and merge into INFO.
		entries := classifyCSQVariant(v, idx)
		if len(entries) > 0 {
			merged := strings.Join(entries, ",")
			setInfoField(v, tag, merged)
		}
		if err := out.Write(v); err != nil {
			return written, fmt.Errorf("write variant: %w", err)
		}
		written++
	}
	if err := out.Flush(); err != nil {
		return written, err
	}
	return written, nil
}

// ensureCSQInfoLine inserts (or replaces) the ##INFO=<ID=...> header
// for the consequence tag. Existing entries are kept untouched if
// they already declare the same ID.
func ensureCSQInfoLine(meta []string, tag string) []string {
	prefix := fmt.Sprintf("##INFO=<ID=%s,", tag)
	for _, m := range meta {
		if strings.HasPrefix(m, prefix) {
			return meta
		}
	}
	line := fmt.Sprintf("##INFO=<ID=%s,Number=.,Type=String,Description=\"Consequence prediction (bcftools csq v1)\">", tag)
	return append(meta, line)
}

// setInfoField inserts or replaces a single INFO key in v.
func setInfoField(v *vcf.Variant, key, val string) {
	if v.Info == nil {
		v.Info = make(map[string]string)
	}
	if _, ok := v.Info[key]; !ok {
		v.InfoOrder = append(v.InfoOrder, key)
	}
	v.Info[key] = val
}

// CSQIndex holds the parsed GFF (gene/mRNA/CDS hierarchy) plus the
// per-contig reference sequence used by the SNP classifier.
type CSQIndex struct {
	// Refs holds the per-contig reference sequence, keyed by FASTA
	// record ID (column 1 of the GFF must match).
	Refs map[string][]byte
	// Transcripts is keyed by transcript ID and lists the CDS exons
	// in genomic order, plus the parent gene name and biotype.
	Transcripts map[string]*CSQTranscript
	// ByChrom indexes transcript IDs by contig for fast position
	// lookup.
	ByChrom map[string][]*CSQTranscript
}

// CSQTranscript is a single transcript's CDS structure.
type CSQTranscript struct {
	ID       string
	Gene     string
	Biotype  string
	Chrom    string
	Strand   gff.Strand
	CDSExons []CSQExon // sorted by genomic Start
}

// CSQExon is one CDS exon (genomic coordinates, 1-based inclusive).
type CSQExon struct {
	Start int
	End   int
	Phase int
}

// loadCSQIndex reads the FASTA + GFF and constructs the cross-reference.
func loadCSQIndex(fastaPath, gffPath string) (*CSQIndex, error) {
	idx := &CSQIndex{
		Refs:        make(map[string][]byte),
		Transcripts: make(map[string]*CSQTranscript),
		ByChrom:     make(map[string][]*CSQTranscript),
	}

	// Load FASTA.
	fr, err := iohelper.OpenReader(fastaPath)
	if err != nil {
		return nil, fmt.Errorf("open fasta %q: %w", fastaPath, err)
	}
	defer fr.Close()
	rec, err := fasta.NewReader(fr).ReadAll()
	if err != nil {
		return nil, fmt.Errorf("read fasta: %w", err)
	}
	for _, r := range rec {
		idx.Refs[r.ID] = r.Sequence
	}

	// Load GFF.
	gr, err := iohelper.OpenReader(gffPath)
	if err != nil {
		return nil, fmt.Errorf("open gff %q: %w", gffPath, err)
	}
	defer gr.Close()
	feats, err := gff.NewReader(gr).ReadAll()
	if err != nil {
		return nil, fmt.Errorf("read gff: %w", err)
	}

	// First pass: map gene ID -> gene name + biotype.
	geneInfo := map[string]struct{ name, biotype string }{}
	for _, f := range feats {
		if f.Type != "gene" {
			continue
		}
		id := f.ID()
		if id == "" {
			continue
		}
		name := f.Attributes["Name"]
		if name == "" {
			name = f.Attributes["gene_name"]
		}
		if name == "" {
			name = id
		}
		biotype := f.Attributes["biotype"]
		if biotype == "" {
			biotype = f.Attributes["gene_biotype"]
		}
		if biotype == "" {
			biotype = "protein_coding"
		}
		geneInfo[id] = struct{ name, biotype string }{name, biotype}
	}

	// Second pass: collect transcripts.
	for _, f := range feats {
		if f.Type != "mRNA" && f.Type != "transcript" {
			continue
		}
		tid := f.ID()
		if tid == "" {
			continue
		}
		parent := f.Parent()
		gi := geneInfo[parent]
		if gi.name == "" {
			gi.name = parent
		}
		if gi.biotype == "" {
			gi.biotype = "protein_coding"
		}
		idx.Transcripts[tid] = &CSQTranscript{
			ID:      tid,
			Gene:    gi.name,
			Biotype: gi.biotype,
			Chrom:   f.Seqid,
			Strand:  f.Strand,
		}
	}

	// Third pass: collect CDS exons under each transcript.
	for _, f := range feats {
		if f.Type != "CDS" {
			continue
		}
		parent := f.Parent()
		// Parent may be a comma-list per GFF3; we accept the first.
		if i := strings.IndexByte(parent, ','); i >= 0 {
			parent = parent[:i]
		}
		t, ok := idx.Transcripts[parent]
		if !ok {
			continue
		}
		t.CDSExons = append(t.CDSExons, CSQExon{
			Start: f.Start,
			End:   f.End,
			Phase: f.Phase,
		})
	}

	// Sort each transcript's exons by genomic start, then build the
	// per-contig index.
	for _, t := range idx.Transcripts {
		sort.Slice(t.CDSExons, func(i, j int) bool { return t.CDSExons[i].Start < t.CDSExons[j].Start })
		idx.ByChrom[t.Chrom] = append(idx.ByChrom[t.Chrom], t)
	}
	return idx, nil
}

// classifyCSQVariant returns BCSQ entries (one per matching transcript)
// for the given variant. Empty slice means "no coding consequence".
func classifyCSQVariant(v *vcf.Variant, idx *CSQIndex) []string {
	if v == nil || idx == nil {
		return nil
	}
	// SNP-only: REF and (at least one) ALT must be a single base.
	if len(v.Ref) != 1 {
		return nil
	}
	alt := ""
	for _, a := range v.Alt {
		if len(a) == 1 && a != "." && a != v.Ref {
			alt = a
			break
		}
	}
	if alt == "" {
		return nil
	}
	pos := v.Pos
	out := []string{}
	for _, t := range idx.ByChrom[v.Chrom] {
		if !cdsCovers(t, pos) {
			continue
		}
		entry, ok := classifyForTranscript(t, idx.Refs[v.Chrom], pos, v.Ref[0], alt[0])
		if !ok {
			continue
		}
		out = append(out, entry)
	}
	return out
}

// cdsCovers reports whether t's CDS exons include pos.
func cdsCovers(t *CSQTranscript, pos int) bool {
	for _, e := range t.CDSExons {
		if pos >= e.Start && pos <= e.End {
			return true
		}
	}
	return false
}

// classifyForTranscript returns the BCSQ entry for one (transcript, SNP)
// pair. The reference base must match v.Ref; if it does not we return
// false so the variant is left unannotated for that transcript.
func classifyForTranscript(t *CSQTranscript, refSeq []byte, pos int, refBase, altBase byte) (string, bool) {
	if len(refSeq) == 0 {
		return "", false
	}
	// Sanity: position must be within sequence and match the reference.
	if pos < 1 || pos > len(refSeq) {
		return "", false
	}
	if refSeq[pos-1] != refBase && upper(refSeq[pos-1]) != upper(refBase) {
		return "", false
	}

	// Compute coding offset (0-based codon position within the CDS).
	codingOff, ok := cdsOffset(t, pos)
	if !ok {
		return "", false
	}

	// Identify the codon: which CDS-coordinate triplet contains
	// codingOff? Then fetch the 3 codon bases from the genomic
	// sequence using cdsOffsetToGenomic for the *other* two positions.
	codonIdx := codingOff / 3
	withinCodon := codingOff % 3

	// Triplet covers cds positions [codonIdx*3, codonIdx*3+1, codonIdx*3+2].
	codonStart := codonIdx * 3
	codon := [3]byte{}
	for i := 0; i < 3; i++ {
		g, ok := cdsToGenomic(t, codonStart+i)
		if !ok {
			return "", false
		}
		codon[i] = upper(refSeq[g-1])
	}
	// On the reverse strand the codon bases come from the reverse
	// complement of the genomic sequence; flip both the codon and the
	// alt base.
	if t.Strand == gff.StrandReverse {
		codon = revcompCodon(codon)
		altBase = complementBase(altBase)
		withinCodon = 2 - withinCodon
	}

	refCodon := codon
	mutCodon := codon
	mutCodon[withinCodon] = upper(altBase)

	refAA := translateCodon(refCodon)
	altAA := translateCodon(mutCodon)

	consequence := "synonymous"
	if refAA == '*' && altAA != '*' {
		consequence = "stop_lost"
	} else if altAA == '*' && refAA != '*' {
		consequence = "stop_gained"
	} else if codingOff == 0 && refAA == 'M' && altAA != 'M' {
		// First codon of CDS — start-loss.
		consequence = "start_lost"
	} else if refAA != altAA {
		consequence = "missense"
	}

	// dna_change: 1-based CDS coordinate of the changed position.
	dnaCoord := codingOff + 1
	dnaChange := fmt.Sprintf("%d%c>%c", dnaCoord, codon[withinCodon], upper(altBase))
	aaChange := fmt.Sprintf("%d%c>%c", codonIdx+1, refAA, altAA)

	// Format: consequence|gene|transcript|biotype|strand|aa_change|dna_change
	strand := "+"
	if t.Strand == gff.StrandReverse {
		strand = "-"
	}
	return fmt.Sprintf("%s|%s|%s|%s|%s|%s|%s",
		consequence, t.Gene, t.ID, t.Biotype, strand, aaChange, dnaChange,
	), true
}

// cdsOffset returns the 0-based offset of pos within the transcript's
// CDS, walking exons in strand order.
func cdsOffset(t *CSQTranscript, pos int) (int, bool) {
	if t.Strand == gff.StrandReverse {
		// Iterate exons from highest genomic coord downward.
		off := 0
		for i := len(t.CDSExons) - 1; i >= 0; i-- {
			e := t.CDSExons[i]
			if pos >= e.Start && pos <= e.End {
				return off + (e.End - pos), true
			}
			off += e.End - e.Start + 1
		}
		return 0, false
	}
	// Forward strand.
	off := 0
	for _, e := range t.CDSExons {
		if pos >= e.Start && pos <= e.End {
			return off + (pos - e.Start), true
		}
		off += e.End - e.Start + 1
	}
	return 0, false
}

// cdsToGenomic converts a CDS-coordinate offset back to a 1-based
// genomic position. Returns false if the offset is past the end of
// the CDS.
func cdsToGenomic(t *CSQTranscript, off int) (int, bool) {
	if t.Strand == gff.StrandReverse {
		for i := len(t.CDSExons) - 1; i >= 0; i-- {
			e := t.CDSExons[i]
			n := e.End - e.Start + 1
			if off < n {
				return e.End - off, true
			}
			off -= n
		}
		return 0, false
	}
	for _, e := range t.CDSExons {
		n := e.End - e.Start + 1
		if off < n {
			return e.Start + off, true
		}
		off -= n
	}
	return 0, false
}

// translateCodon returns the standard-table single-letter amino acid
// for the given 3-base codon. Unknown bases yield 'X'.
func translateCodon(codon [3]byte) byte {
	const table = "" +
		// T              C              A              G
		"FFLLSSSSYY**CC*W" + // TTx TCx TAx TGx
		"LLLLPPPPHHQQRRRR" + // CTx CCx CAx CGx
		"IIIMTTTTNNKKSSRR" + // ATx ACx AAx AGx
		"VVVVAAAADDEEGGGG" //   GTx GCx GAx GGx
	idx := codonIndex(codon)
	if idx < 0 {
		return 'X'
	}
	return table[idx]
}

func codonIndex(c [3]byte) int {
	i0 := baseIndex(c[0])
	i1 := baseIndex(c[1])
	i2 := baseIndex(c[2])
	if i0 < 0 || i1 < 0 || i2 < 0 {
		return -1
	}
	return i0*16 + i1*4 + i2
}

func baseIndex(b byte) int {
	switch b {
	case 'T', 't', 'U', 'u':
		return 0
	case 'C', 'c':
		return 1
	case 'A', 'a':
		return 2
	case 'G', 'g':
		return 3
	}
	return -1
}

func complementBase(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'T':
		return 'A'
	case 'U':
		return 'A'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'a':
		return 't'
	case 't':
		return 'a'
	case 'u':
		return 'a'
	case 'c':
		return 'g'
	case 'g':
		return 'c'
	}
	return 'N'
}

func revcompCodon(c [3]byte) [3]byte {
	return [3]byte{
		complementBase(c[2]),
		complementBase(c[1]),
		complementBase(c[0]),
	}
}

func upper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

// ParseCSQPhase parses upstream's -p/--phase {a|m|r|R|s}.
func ParseCSQPhase(s string) (byte, error) {
	if s == "" {
		return 'r', nil
	}
	switch s[0] {
	case 'a', 'm', 'r', 'R', 's':
		return s[0], nil
	}
	return 0, fmt.Errorf("bcftools csq: bad -p/--phase %q (want a|m|r|R|s)", s)
}
