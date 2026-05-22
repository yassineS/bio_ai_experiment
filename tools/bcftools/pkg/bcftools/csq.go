// bcftools csq — predict variant consequences against a GFF annotation.
//
// SCOPE: this port implements the full PER-RECORD (single-variant,
// non-haplotype-phased) consequence classifier. The haplotype engine
// (compound consequences, the hap_node_t tree, -p/--phase modes,
// -n/--ncsq) is a separate, later slice; see docs/PARITY_ROADMAP.md
// "csq full-parity slicing plan" for the ordered roadmap.
//
// What this port DOES do:
//   - Load a GFF3 file and build per-transcript CDS / exon / UTR /
//     transcript-span indexes keyed by transcript ID.
//   - Load the reference FASTA into memory.
//   - For each input VCF record, classify each ALT allele against every
//     overlapping transcript and emit INFO/BCSQ as
//     "consequence|gene|transcript|biotype|strand" (one entry per
//     matching transcript, comma-separated).
//   - Cover the full per-record SO-term set: synonymous, missense,
//     stop_gained, stop_lost, start_lost, splice_donor/acceptor/region,
//     5_prime_utr, 3_prime_utr, intron, non_coding, and the indel
//     consequences (inframe insertion/deletion, frameshift,
//     feature_elongation/truncation). The dispatch and splice machinery
//     are ported from csq.c's test_cds/test_utr/test_splice/
//     test_tscript and the splice_csq family — see csq_classify.go.
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

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
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
		// classifyCSQRecord covers the full per-record SO-term set
		// (CDS / UTR / intron / non-coding / splice, SNPs + indels).
		entries := classifyCSQRecord(v, idx)
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

// CSQTranscript is a single transcript's structure: its CDS exons,
// full exons, derived UTR regions and transcript span. These mirror
// upstream csq.c's idx_cds / idx_exon / idx_utr / idx_tscript indexes.
type CSQTranscript struct {
	ID       string
	Gene     string
	Biotype  string
	Chrom    string
	Strand   gff.Strand
	CDSExons []CSQExon // CDS exons, sorted by genomic Start
	Exons    []CSQExon // full exons (CDS + UTR), sorted by genomic Start
	UTRs     []CSQUTR  // 5'/3' UTR regions, sorted by genomic Start

	// Beg and End are the transcript's genomic span (1-based,
	// inclusive). When no transcript/mRNA feature was seen they are
	// derived as the min/max over Exons (or CDSExons).
	Beg int
	End int

	// Coding reports whether the transcript has any CDS exon. Upstream
	// distinguishes coding transcripts (intron consequence) from
	// non-coding ones (non_coding consequence) via GF_is_coding.
	Coding bool

	// Trim5 / Trim3 mark an incomplete CDS at the 5' / 3' end of the
	// transcript. Upstream's TRIM_5PRIME / TRIM_3PRIME suppress
	// start_lost / stop_lost on incomplete annotations.
	Trim5 bool
	Trim3 bool
}

// CSQExon is one exon (genomic coordinates, 1-based inclusive). For a
// CDS exon Phase holds the GFF3 reading-frame phase; for a full exon
// Phase is unused.
type CSQExon struct {
	Start int
	End   int
	Phase int
}

// CSQUTR is one untranslated-region span. Prime5 is true for a
// 5'-UTR, false for a 3'-UTR.
type CSQUTR struct {
	Start  int
	End    int
	Prime5 bool
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

	// Second pass: collect transcripts. A transcript is any feature
	// whose type is not gene/CDS/exon/UTR but which has an ID; csq
	// recognises mRNA, transcript and the various non-coding RNA
	// biotypes (lnc_RNA, miRNA, ...). We treat anything that is a
	// declared Parent of a CDS or exon as a transcript.
	for _, f := range feats {
		if !isTranscriptType(f.Type) {
			continue
		}
		tid := f.ID()
		if tid == "" {
			continue
		}
		parent := firstParent(f.Parent())
		gi := geneInfo[parent]
		if gi.name == "" {
			gi.name = parent
		}
		if gi.biotype == "" {
			gi.biotype = "protein_coding"
		}
		biotype := f.Attributes["biotype"]
		if biotype == "" {
			biotype = f.Attributes["transcript_biotype"]
		}
		if biotype == "" {
			biotype = gi.biotype
		}
		idx.Transcripts[tid] = &CSQTranscript{
			ID:      tid,
			Gene:    gi.name,
			Biotype: biotype,
			Chrom:   f.Seqid,
			Strand:  f.Strand,
			Beg:     f.Start,
			End:     f.End,
		}
	}

	// Third pass: collect CDS and full exons under each transcript.
	for _, f := range feats {
		switch f.Type {
		case "CDS":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.CDSExons = append(t.CDSExons, CSQExon{Start: f.Start, End: f.End, Phase: f.Phase})
			}
		case "exon":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.Exons = append(t.Exons, CSQExon{Start: f.Start, End: f.End})
			}
		case "five_prime_UTR", "5UTR", "five_prime_utr":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.UTRs = append(t.UTRs, CSQUTR{Start: f.Start, End: f.End, Prime5: true})
			}
		case "three_prime_UTR", "3UTR", "three_prime_utr":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.UTRs = append(t.UTRs, CSQUTR{Start: f.Start, End: f.End, Prime5: false})
			}
		}
	}

	// Finalise each transcript: sort exons, derive missing UTRs and
	// transcript span, set Coding / Trim flags, then index by contig.
	for _, t := range idx.Transcripts {
		sort.Slice(t.CDSExons, func(i, j int) bool { return t.CDSExons[i].Start < t.CDSExons[j].Start })
		sort.Slice(t.Exons, func(i, j int) bool { return t.Exons[i].Start < t.Exons[j].Start })
		finalizeTranscript(t, idx.Refs[t.Chrom])
		sort.Slice(t.UTRs, func(i, j int) bool { return t.UTRs[i].Start < t.UTRs[j].Start })
		idx.ByChrom[t.Chrom] = append(idx.ByChrom[t.Chrom], t)
	}
	return idx, nil
}

// isTranscriptType reports whether a GFF3 feature type denotes a
// transcript (the second level of the gene/transcript/exon hierarchy).
func isTranscriptType(t string) bool {
	switch t {
	case "mRNA", "transcript", "lnc_RNA", "lncRNA", "ncRNA", "miRNA",
		"snRNA", "snoRNA", "rRNA", "tRNA", "scRNA", "pseudogenic_transcript",
		"processed_transcript", "nc_primary_transcript", "antisense_RNA":
		return true
	}
	return false
}

// firstParent returns the first ID in a possibly comma-separated GFF3
// Parent attribute.
func firstParent(p string) string {
	if i := strings.IndexByte(p, ','); i >= 0 {
		return p[:i]
	}
	return p
}

// finalizeTranscript fills in derived fields after the GFF passes:
// Coding, the transcript span, the incomplete-CDS trim flags, and any
// UTR regions not explicitly present in the GFF (derived as the parts
// of full exons that fall outside the CDS span).
func finalizeTranscript(t *CSQTranscript, ref []byte) {
	t.Coding = len(t.CDSExons) > 0

	// Transcript span: prefer the mRNA/transcript feature's own span;
	// otherwise derive from exons or CDS exons.
	if t.Beg == 0 || t.End == 0 {
		spanFrom := t.Exons
		if len(spanFrom) == 0 {
			spanFrom = t.CDSExons
		}
		for i, e := range spanFrom {
			if i == 0 || e.Start < t.Beg {
				t.Beg = e.Start
			}
			if i == 0 || e.End > t.End {
				t.End = e.End
			}
		}
	}

	// Derive UTRs from exon minus CDS when none were given explicitly.
	if len(t.UTRs) == 0 && t.Coding && len(t.Exons) > 0 {
		cdsBeg, cdsEnd := t.CDSExons[0].Start, t.CDSExons[len(t.CDSExons)-1].End
		for _, e := range t.Exons {
			if e.Start < cdsBeg {
				end := e.End
				if end >= cdsBeg {
					end = cdsBeg - 1
				}
				t.UTRs = append(t.UTRs, CSQUTR{Start: e.Start, End: end, Prime5: t.Strand != gff.StrandReverse})
			}
			if e.End > cdsEnd {
				start := e.Start
				if start <= cdsEnd {
					start = cdsEnd + 1
				}
				t.UTRs = append(t.UTRs, CSQUTR{Start: start, End: e.End, Prime5: t.Strand == gff.StrandReverse})
			}
		}
	}

	// Mark an incomplete CDS: a coding transcript whose first codon
	// (in transcript orientation) is not ATG. Upstream uses this to
	// suppress spurious start_lost / stop_lost calls.
	if t.Coding && len(ref) > 0 {
		if c, ok := codingFirstCodon(t, ref); ok && c != "ATG" {
			t.Trim5 = true
		}
	}

	// Mark an incomplete 3' CDS: a coding transcript whose total
	// coding length is not a multiple of three. Upstream (gff.c:844
	// `if (len%3 != 0) tr->trim |= TRIM_3PRIME`) computes len over the
	// phase-trimmed CDS exons, so the equivalent here is the sum of CDS
	// exon lengths minus the 5' reading-frame phase. Without this flag
	// an incomplete-3' transcript still runs the stop check in
	// classifyCDS and can emit a spurious stop_lost on its last,
	// non-stop codon. csq.c:1646-1650 gates check_stop on !TRIM_3PRIME.
	if t.Coding {
		codingLen := -transcriptFirstPhase(t)
		for _, e := range t.CDSExons {
			codingLen += e.End - e.Start + 1
		}
		if codingLen%3 != 0 {
			t.Trim3 = true
		}
	}
}

// codingFirstCodon returns the transcript's first CDS codon (ATG for a
// complete annotation), reading in transcript orientation.
func codingFirstCodon(t *CSQTranscript, ref []byte) (string, bool) {
	var b [3]byte
	for i := 0; i < 3; i++ {
		g, ok := cdsToGenomic(t, i)
		if !ok || g < 1 || g > len(ref) {
			return "", false
		}
		b[i] = upper(ref[g-1])
	}
	if t.Strand == gff.StrandReverse {
		b = revcompCodon(b)
	}
	return string(b[:]), true
}

// cdsCovers reports whether t's CDS exons include pos. The full
// per-record dispatch lives in classifyCSQRecord (csq_classify.go);
// cdsCovers is retained as a small reusable predicate.
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
	} else if codonIdx == 0 && refAA == 'M' && altAA != 'M' {
		// First codon of CDS — start-loss. Any of the 3 ATG
		// positions count, not just the leading base.
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
	// Honour GFF3 Phase on the first CDS exon (in transcript order).
	// Per GFF3 spec, Phase = number of bases of the 5' end of the
	// exon that are leftover from the previous codon, so the reading
	// frame is offset by `phase`. Most Ensembl/GENCODE transcripts
	// whose CDS starts mid-exon (e.g. truncated annotations) have
	// non-zero phase here; ignoring it puts every codon out of frame.
	framePhase := transcriptFirstPhase(t)
	if t.Strand == gff.StrandReverse {
		// Iterate exons from highest genomic coord downward.
		off := 0
		for i := len(t.CDSExons) - 1; i >= 0; i-- {
			e := t.CDSExons[i]
			if pos >= e.Start && pos <= e.End {
				return off + (e.End - pos) - framePhase, true
			}
			off += e.End - e.Start + 1
		}
		return 0, false
	}
	// Forward strand.
	off := 0
	for _, e := range t.CDSExons {
		if pos >= e.Start && pos <= e.End {
			return off + (pos - e.Start) - framePhase, true
		}
		off += e.End - e.Start + 1
	}
	return 0, false
}

// transcriptFirstPhase returns the GFF3 phase of the first CDS exon
// in transcript order (highest-Start for reverse strand, lowest-Start
// for forward). Defaults to 0 when the parser hasn't recorded a phase.
func transcriptFirstPhase(t *CSQTranscript) int {
	if len(t.CDSExons) == 0 {
		return 0
	}
	if t.Strand == gff.StrandReverse {
		p := t.CDSExons[len(t.CDSExons)-1].Phase
		if p < 0 || p > 2 {
			return 0
		}
		return p
	}
	p := t.CDSExons[0].Phase
	if p < 0 || p > 2 {
		return 0
	}
	return p
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
