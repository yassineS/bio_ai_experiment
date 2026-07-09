// bcftools csq — haplotype-aware variant consequence caller against a
// GFF annotation.
//
// SCOPE: this port implements the full haplotype-aware engine. Variants
// overlapping a coding region are phased onto per-sample haplotypes and
// walked together (the hap_node_t tree), so compound consequences are
// called jointly and the reference vs haplotype-altered CDS are
// translated and diffed. The INFO/BCSQ output matches upstream
// byte-for-byte on the targeted goldens (see csq_golden_test.go).
//
// File roles:
//   - csq.go            — GFF/FASTA index build (CSQIndex/CSQTranscript)
//                         plus the CSQ / CSQFile entry points.
//   - csq_classify.go   — the shared CSQ_* SO-term bit constants and the
//                         csqStrings table consumed by the engine.
//   - csq_hap.go        — engine data model: hapEngine, hapTranscript,
//                         hapNode, the genetic-code tables.
//   - csq_splice.go     — splice_csq with set_refalt / splice_build_hap.
//   - csq_engine.go     — hap_init, cds_translate, hap_finalize,
//                         hap_add_csq, csq_push, kput_vcsq.
//   - csq_process.go    — the vbuf/pos2vbuf buffer and the
//                         test_cds/utr/splice/tscript dispatch.
//
// Slice 4 added the GFF/output tail: FORMAT/TBCSQ per-haplotype text
// expansion (the query path, see query.go expandTBCSQ),
// --unify-chr-names contig reconciliation (parseUnifyChrNames /
// unifyChrName), --dump-gff model dumping (csq_dump.go, byte-exact vs
// upstream gff_dump), and non-text -O b|u|z output via the in-tree
// BCF/BGZF writers (openCSQOutput). -l/--local-csq (the per-record,
// non-haplotype-aware caller, upstream test_cds_local) is ported in
// csq_local.go and selected when CSQOptions.LocalCSQ is set. The CLI
// accepts the full upstream getopt_long surface; flags we don't yet
// honour either no-op or hard-reject with a roadmap pointer.

package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
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

	// LocalCSQ matches upstream's -l/--local-csq: predict consequences
	// per-record (non-haplotype-aware) via the ported test_cds_local
	// (see csq_local.go). Each record's coding consequence is derived
	// from its own ref/alt against the spliced reference, so compound
	// consequences spanning several records are not joined.
	LocalCSQ bool

	// Phase is upstream's -p/--phase {a|m|r|R|s}. v1 stores it for
	// future use; the per-record SNP classifier does not depend on
	// haplotype phasing.
	Phase byte

	// NCSQ is upstream's -n/--ncsq cap on per-haplotype consequences.
	// v1 reports every matching transcript without a cap; the field is
	// retained for parity.
	NCSQ int
	// TrimProteinSeq is upstream's -B/--trim-protein-seq (and the alias
	// -b/--brief-predictions, which sets it to 1). When >0 each
	// amino-acid prediction in INFO/BCSQ is abbreviated to its first N
	// residues followed by "..<index>". Mirrors args->brief_predictions.
	TrimProteinSeq int
	// GeneticCode is upstream's -C/--genetic-code: the NCBI translation
	// table id used for codon->amino-acid translation. 0 is the standard
	// simplified table; see gencodeTables for the supported ids.
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

	// OutputFormat — upstream supports v/z/u/b/t. This port emits v
	// (VCF text), z (BGZF VCF), b (BCF) and u (uncompressed BCF) via the
	// in-tree writers. The streaming-text `t` form (FT_TAB_TEXT) is
	// selected by TextOutput instead, and ignores OutputFormat.
	OutputFormat OutputFormat
	// TextOutput selects upstream's `-O t` streaming-text output
	// (FT_TAB_TEXT): instead of a VCF/BCF stream, csq emits a per-
	// consequence "CSQ\t<sample>\t<haplotype>\t<chrom>\t<pos>\t<csq>"
	// text stream (with a leading provenance/column-header block),
	// matching upstream csq.c's text_print_vcsq path byte-for-byte.
	TextOutput bool
	// TextArgv is the original CLI argv (e.g. ["-f", "ref.fa", "-g",
	// "ann.gff", "-O", "t", "in.vcf"]) used to reproduce upstream's
	// "# The command line was:" provenance line in text mode. It is
	// ignored unless TextOutput is set.
	TextArgv []string
	// CompressLevel is the gzip level for -O z output (negative selects
	// the package default).
	CompressLevel int
	// Threads is upstream's -@/--threads. When greater than 1 it enables
	// parallel BGZF compression of -O z and -O b output via bgzf.MultiWriter;
	// the framed result decodes byte-identically regardless of thread count.
	Threads int

	// DumpGFF is upstream's --dump-gff FILE. When set, CSQFile writes the
	// parsed GFF model (genes/transcripts/CDS/UTR/exons) to FILE as a
	// BGZF-compressed trimmed GFF3, byte-exact with upstream gff_dump on
	// position-ordered inputs (see DumpGFF / csq_dump.go).
	DumpGFF string
	// UnifyChrNames is upstream's --unify-chr-names VCF,GFF,FAI. The
	// special value "0" (or empty) disables rewriting; otherwise the
	// three comma-separated prefixes reconcile VCF/GFF/FASTA contig
	// namespaces (see parseUnifyChrNames / unifyChrName).
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
	prefixVCF, prefixGFF, prefixFAI, err := parseUnifyChrNames(opts.UnifyChrNames)
	if err != nil {
		return 0, err
	}
	idx, err := loadCSQIndexUnified(opts.FastaRef, opts.GFFAnnot, prefixVCF, prefixGFF, prefixFAI,
		csqPhaseCheck{Force: opts.Force, Verbosity: opts.Verbosity})
	if err != nil {
		return 0, err
	}
	if opts.DumpGFF != "" {
		if err := dumpCSQGFFToFile(opts.DumpGFF, idx); err != nil {
			return 0, err
		}
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

	// Inject INFO meta-line for the consequence tag. In streaming-text
	// mode (-O t) upstream never writes a VCF header, but the INFO/FORMAT
	// description text is still required so the engine's consequence
	// rendering matches; the lines are simply not emitted.
	tag := opts.CustomTag
	hdr.MetaInfo = ensureCSQInfoLine(hdr.MetaInfo, tag)
	if len(hdr.Samples) > 0 {
		hdr.MetaInfo = ensureCSQFormatLine(hdr.MetaInfo, tag)
	}

	regionSpecs := append([]string(nil), opts.Regions...)
	regionSpecs = append(regionSpecs, opts.Targets...)

	// The haplotype-aware engine buffers variants until the enclosing
	// transcript boundary is crossed, so compound consequences across
	// several records can be called jointly. Records are emitted in
	// input order once their annotations are final.
	opts.CustomTag = tag

	if opts.TextOutput {
		return csqText(w, vr, idx, opts, hdr, regionSpecs)
	}

	out, cleanup, err := openCSQOutput(w, opts.OutputFormat, opts.CompressLevel, opts.Threads, hdr)
	if err != nil {
		return 0, err
	}
	defer cleanup()
	if err := out.WriteHeader(); err != nil {
		return 0, fmt.Errorf("write VCF header: %w", err)
	}

	eng := newHapEngine(idx, opts, hdr)
	written := 0
	emit := func() error {
		for _, v := range eng.out {
			if err := out.Write(v); err != nil {
				return fmt.Errorf("write variant: %w", err)
			}
			written++
		}
		eng.out = eng.out[:0]
		return nil
	}
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
		if err := eng.process(v); err != nil {
			return written, err
		}
		if err := emit(); err != nil {
			return written, err
		}
	}
	eng.finish()
	if err := emit(); err != nil {
		return written, err
	}
	if err := out.Flush(); err != nil {
		return written, err
	}
	return written, nil
}

// csqText drives the engine in -O t streaming-text mode. It writes the
// upstream provenance/column-header block, then emits each finalised
// "CSQ\t..." text line as records flush. It returns the number of text
// lines written (the CSQ rows, excluding the header block).
func csqText(w io.Writer, vr *vcf.Reader, idx *CSQIndex, opts CSQOptions, hdr *vcf.Header, regionSpecs []string) (int, error) {
	bw := bufio.NewWriter(w)
	if err := writeCSQTextHeader(bw, opts.TextArgv); err != nil {
		return 0, err
	}

	eng := newHapEngine(idx, opts, hdr)
	written := 0
	emit := func() error {
		for _, line := range eng.outText {
			if _, err := bw.WriteString(line); err != nil {
				return err
			}
			written++
		}
		eng.outText = eng.outText[:0]
		return nil
	}
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
		if err := eng.process(v); err != nil {
			return written, err
		}
		if err := emit(); err != nil {
			return written, err
		}
	}
	eng.finish()
	if err := emit(); err != nil {
		return written, err
	}
	return written, bw.Flush()
}

// writeCSQTextHeader writes the four-line header block upstream emits at
// the top of -O t output: a version provenance line, the command line,
// a "# LOG" line, and the "# CSQ" column header. The first two lines are
// version/command provenance (parity tests strip them); the latter two
// are stable.
func writeCSQTextHeader(w io.Writer, argv []string) error {
	if _, err := fmt.Fprintf(w, "# This file was produced by: bcftools +csq(%s+htslib-%s)\n",
		csqProducedVersion(), csqHtslibVersion()); err != nil {
		return err
	}
	var sb strings.Builder
	sb.WriteString("# The command line was:\tbcftools +csq")
	for _, a := range argv {
		sb.WriteByte(' ')
		sb.WriteString(a)
	}
	sb.WriteByte('\n')
	if _, err := io.WriteString(w, sb.String()); err != nil {
		return err
	}
	if _, err := io.WriteString(w, "# LOG\t[2]Message\n"); err != nil {
		return err
	}
	_, err := io.WriteString(w,
		"# CSQ\t[2]Sample\t[3]Haplotype\t[4]Chromosome\t[5]Position\t[6]Consequence\n")
	return err
}

// csqProducedVersion / csqHtslibVersion return the version strings
// embedded in the -O t provenance line. They are placeholders for
// parity: upstream emits its own build version here, which always
// differs, so the parity tests strip this line. Keeping a stable,
// non-empty value documents the format without claiming byte parity on
// the version itself.
func csqProducedVersion() string { return "bio_ai_experiment" }
func csqHtslibVersion() string   { return "bio_ai_experiment" }

// openCSQOutput wraps the destination writer in a variantWriter for the
// requested -O format. Supports VCF text, VCF.gz, compressed BCF and
// uncompressed BCF (shared with `bcftools view` so the BCF writer logic
// stays in one place).
func openCSQOutput(w io.Writer, format OutputFormat, level, threads int, hdr *vcf.Header) (variantWriter, func(), error) {
	return openOutput(w, ViewOptions{OutputFormat: format, CompressLevel: level, Threads: threads}, hdr)
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
	line := fmt.Sprintf("##INFO=<ID=%s,Number=.,Type=String,Description=\"Haplotype-aware consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change\">", tag)
	return append(meta, line)
}

// ensureCSQFormatLine inserts the ##FORMAT=<ID=...> header for the
// per-sample bitmask emitted alongside INFO/BCSQ. The description
// matches upstream csq.c:800 so `bcftools query -f'[%TBCSQ]'` can
// detect and expand the bitmask.
func ensureCSQFormatLine(meta []string, tag string) []string {
	prefix := fmt.Sprintf("##FORMAT=<ID=%s,", tag)
	for _, m := range meta {
		if strings.HasPrefix(m, prefix) {
			return meta
		}
	}
	line := fmt.Sprintf("##FORMAT=<ID=%s,Number=.,Type=Integer,Description=\"Bitmask of indexes to INFO/BCSQ, with interleaved first/second haplotype. Use \\\"bcftools query -f'[%%CHROM\\t%%POS\\t%%SAMPLE\\t%%TBCSQ\\n]'\\\" to translate.\">", tag)
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

	// SpliceExons is the flattened per-contig exon list of all coding
	// transcripts, ordered the way htslib's idx_exon returns overlaps (by
	// start ascending, end descending). test_splice iterates idx_exon, not the
	// per-transcript exon lists, so two overlapping transcripts contribute
	// their splice_region consequences in exon-position order — which differs
	// from transcript-start order and decides the INFO/BCSQ entry order.
	SpliceExons map[string][]csqExonRef

	// UTRRegions is the flattened per-contig UTR list of all transcripts,
	// ordered the way htslib's idx_utr returns overlaps (start ascending,
	// end descending). test_utr iterates idx_utr, not the per-transcript
	// UTR lists, so two overlapping transcripts contribute their
	// 5'/3'_prime_utr consequences in UTR-position order — which differs
	// from transcript-start order and decides the INFO/BCSQ entry order.
	UTRRegions map[string][]csqUTRRef

	// Genes holds the parsed gene features in GFF appearance order.
	// Retained only for --dump-gff (mirrors upstream's gid2gene table).
	Genes []*CSQGene
}

// CSQGene is one gene feature retained for the --dump-gff output. It
// mirrors the upstream gff.c gf_gene_t fields that gff_dump emits.
type CSQGene struct {
	// ID is the bare gene accession (Ensembl "gene:" prefix stripped).
	ID string
	// Name is the gene Name/gene_name attribute (falls back to ID).
	Name string
	// Chrom is the contig (rewritten under --unify-chr-names).
	Chrom string
	// Strand is the gene strand.
	Strand gff.Strand
	// Beg and End are the 1-based inclusive gene span.
	Beg int
	End int
	// Used reports whether at least one of the gene's transcripts was
	// linked to a child CDS/exon/UTR feature, matching upstream's
	// gene->used flag set in gff_parse.
	Used bool
}

// CSQTranscript is a single transcript's structure: its CDS exons,
// full exons, derived UTR regions and transcript span. These mirror
// upstream csq.c's idx_cds / idx_exon / idx_utr / idx_tscript indexes.
type CSQTranscript struct {
	ID       string
	Gene     string
	GeneID   string // bare parent gene accession, for --dump-gff Parent=
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

	// Coding reports whether the transcript's biotype is protein-coding,
	// derived from the biotype string via codingBiotype (mirroring
	// upstream gff.c gff_parse_biotype + GF_is_coding). It is NOT a
	// function of CDS presence: a transcript may carry CDS exons yet be
	// non-coding by biotype (skipped by test_cds), and the intron vs
	// non_coding consequence at the transcript level is chosen from this
	// flag. CDS-derived processing (phase trim, UTR derivation, 3' trim,
	// the splice index) instead keys off actual CDS-exon presence.
	Coding bool

	// Trim5 / Trim3 mark an incomplete CDS at the 5' / 3' end of the
	// transcript. Upstream's TRIM_5PRIME / TRIM_3PRIME suppress
	// start_lost / stop_lost on incomplete annotations.
	Trim5 bool
	Trim3 bool

	// Used reports whether the transcript was linked to at least one
	// child CDS/exon/UTR feature. Mirrors upstream's tr->used flag and
	// is emitted by --dump-gff.
	Used bool
}

// csqExonRef is one exon paired with its parent transcript, used to
// build the flattened idx_exon-ordered splice index (CSQIndex.SpliceExons).
type csqExonRef struct {
	Start int
	End   int
	Tr    *CSQTranscript
}

// csqUTRRef is one UTR span paired with its parent transcript, used to
// build the flattened idx_utr-ordered UTR index (CSQIndex.UTRRegions).
type csqUTRRef struct {
	Start  int
	End    int
	Prime5 bool
	Tr     *CSQTranscript
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

// parseUnifyChrNames decodes the upstream --unify-chr-names spec
// `VCF,GFF,FAI` into three prefixes (empty for "-" or "0", which
// disables the unifier). Returns an error for malformed input.
func parseUnifyChrNames(spec string) (vcfPfx, gffPfx, faiPfx string, err error) {
	if spec == "" || spec == "0" {
		return "", "", "", nil
	}
	parts := strings.Split(spec, ",")
	if len(parts) != 3 {
		return "", "", "", fmt.Errorf("bcftools csq: --unify-chr-names: expected three comma-separated prefixes, got %q", spec)
	}
	for i, p := range parts {
		if p == "-" {
			parts[i] = ""
		}
	}
	return parts[0], parts[1], parts[2], nil
}

// unifyChrName ports csq.c's unify_chr_name: strip the source prefix,
// then prepend the destination prefix. Empty prefixes pass the name
// through unchanged.
func unifyChrName(chr, srcPfx, dstPfx string) string {
	if srcPfx == "" && dstPfx == "" {
		return chr
	}
	if srcPfx != "" && strings.HasPrefix(chr, srcPfx) {
		chr = chr[len(srcPfx):]
	}
	if dstPfx != "" {
		return dstPfx + chr
	}
	return chr
}

// csqPhaseCheck carries the --force / --verbosity settings that govern
// the per-CDS reading-frame phase-consistency check performed while
// finalising each transcript (mirroring upstream gff.c
// tscript_init_cds). Force turns the hard "inconsistent phase" error
// into a warn-and-skip; Verbosity gates whether the warning prints.
type csqPhaseCheck struct {
	Force     bool
	Verbosity int
}

// defaultCSQPhaseCheck matches upstream's csq defaults (force off,
// verbosity 1) so the phase-consistency check errors out on a bad GFF
// unless the caller opts into --force. Test helpers that don't care
// about phase validation use this via loadCSQIndex.
func defaultCSQPhaseCheck() csqPhaseCheck { return csqPhaseCheck{Force: false, Verbosity: 1} }

// loadCSQIndex reads the FASTA + GFF and constructs the cross-reference.
func loadCSQIndex(fastaPath, gffPath string) (*CSQIndex, error) {
	return loadCSQIndexUnified(fastaPath, gffPath, "", "", "", defaultCSQPhaseCheck())
}

// loadCSQIndexUnified is loadCSQIndex with --unify-chr-names prefixes.
// All GFF / FASTA contig keys are rewritten into VCF-prefix form so
// the engine's per-record lookups can use rec.Chrom directly. pc carries
// the --force / --verbosity settings for the phase-consistency check.
func loadCSQIndexUnified(fastaPath, gffPath, prefixVCF, prefixGFF, prefixFAI string, pc csqPhaseCheck) (*CSQIndex, error) {
	idx := &CSQIndex{
		Refs:        make(map[string][]byte),
		Transcripts: make(map[string]*CSQTranscript),
		ByChrom:     make(map[string][]*CSQTranscript),
		SpliceExons: make(map[string][]csqExonRef),
		UTRRegions:  make(map[string][]csqUTRRef),
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
		// Rewrite FASTA contig names into VCF-prefix form when
		// --unify-chr-names is in effect, so the engine's
		// idx.Refs[t.Chrom] lookups (with t.Chrom also rewritten
		// below) succeed transparently.
		idx.Refs[unifyChrName(r.ID, prefixFAI, prefixVCF)] = r.Sequence
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

	// First pass: map gene ID -> gene name + biotype, and retain the
	// gene span/strand for --dump-gff. genesByID points into idx.Genes
	// (preserving GFF appearance order) so the transcript pass can mark
	// the parent gene Used.
	geneInfo := map[string]struct{ name, biotype string }{}
	genesByID := map[string]*CSQGene{}
	for _, f := range feats {
		if !isGeneType(f.Type) {
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
		// Key by the bare accession (the GFF "gene:" prefix stripped), so
		// the transcript pass — which records t.GeneID via
		// stripGFFPrefix(parent) — can find and mark the parent gene Used.
		bareID := stripGFFPrefix(id)
		if _, seen := genesByID[bareID]; !seen {
			g := &CSQGene{
				ID:     bareID,
				Name:   name,
				Chrom:  unifyChrName(f.Seqid, prefixGFF, prefixVCF),
				Strand: f.Strand,
				Beg:    f.Start,
				End:    f.End,
			}
			genesByID[bareID] = g
			idx.Genes = append(idx.Genes, g)
		}
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
			ID:      stripGFFPrefix(tid),
			Gene:    gi.name,
			GeneID:  stripGFFPrefix(parent),
			Biotype: biotype,
			Chrom:   unifyChrName(f.Seqid, prefixGFF, prefixVCF),
			Strand:  f.Strand,
			Beg:     f.Start,
			End:     f.End,
		}
	}

	// Third pass: collect CDS and full exons under each transcript. A
	// transcript (and, transitively, its gene) is marked Used the moment
	// it is linked to a child feature, mirroring upstream gff_parse.
	for _, f := range feats {
		switch f.Type {
		case "CDS":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.CDSExons = append(t.CDSExons, CSQExon{Start: f.Start, End: f.End, Phase: f.Phase})
				t.Used = true
			}
		case "exon":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.Exons = append(t.Exons, CSQExon{Start: f.Start, End: f.End})
				t.Used = true
			}
		case "five_prime_UTR", "5UTR", "five_prime_utr":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.UTRs = append(t.UTRs, CSQUTR{Start: f.Start, End: f.End, Prime5: true})
				t.Used = true
			}
		case "three_prime_UTR", "3UTR", "three_prime_utr":
			if t := idx.Transcripts[firstParent(f.Parent())]; t != nil {
				t.UTRs = append(t.UTRs, CSQUTR{Start: f.Start, End: f.End, Prime5: false})
				t.Used = true
			}
		}
	}

	// Propagate transcript Used flags up to their parent genes for the
	// --dump-gff `used=` column.
	for _, t := range idx.Transcripts {
		if t.Used {
			if g := genesByID[t.GeneID]; g != nil {
				g.Used = true
			}
		}
	}

	// Finalise each transcript: sort exons, derive missing UTRs and
	// transcript span, set Coding / Trim flags, then index by contig.
	//
	// Iterate in a deterministic transcript-ID order. Upstream walks the
	// khash id2tr in unspecified order, but the phase-consistency check
	// below can hard-error naming a specific transcript, so a stable order
	// makes our error message reproducible run-to-run.
	tids := make([]string, 0, len(idx.Transcripts))
	for tid := range idx.Transcripts {
		tids = append(tids, tid)
	}
	sort.Strings(tids)
	// wrongPhaseWarned mirrors gff->warned.wrong_phase: the warn-and-skip
	// message prints once at verbosity 1 and every time at verbosity > 1.
	wrongPhaseWarned := 0
	for _, tid := range tids {
		t := idx.Transcripts[tid]
		sort.Slice(t.CDSExons, func(i, j int) bool { return t.CDSExons[i].Start < t.CDSExons[j].Start })
		sort.Slice(t.Exons, func(i, j int) bool { return t.Exons[i].Start < t.Exons[j].Start })
		// Phase-consistency check (upstream gff.c tscript_init_cds): verify
		// that every CDS exon's GFF3 phase agrees with the reading frame
		// implied by the cumulative CDS length. Without --force a mismatch
		// is a hard error; with --force the transcript is skipped (its CDS
		// exons dropped so it degrades to an intron container).
		if bad, badPos, phase, plen := transcriptPhaseInconsistency(t); bad {
			if !pc.Force {
				return nil, fmt.Errorf("Error: GFF3 assumption failed for transcript %s, CDS=%d: phase!=len%%3 (phase=%d, len=%d). Use the --force option to proceed anyway (at your own risk).", t.ID, badPos, phase, plen)
			}
			if pc.Verbosity > 0 {
				if wrongPhaseWarned == 0 || pc.Verbosity > 1 {
					fmt.Fprintf(os.Stderr, "Warning: The GFF has inconsistent phase column in transcript %s, skipping. CDS pos=%d: phase!=len%%3 (phase=%d, len=%d)\n", t.ID, badPos, phase, plen)
				}
				wrongPhaseWarned++
			}
			// Drop the CDS so the transcript degrades to a bare intron
			// container (matching upstream, which keeps it in idx_tscript
			// but skips CDS indexing). Keep it in ByChrom so region
			// lookups still classify overlapping variants as intron.
			t.CDSExons = nil
		}
		finalizeTranscript(t, idx.Refs[t.Chrom])
		sort.Slice(t.UTRs, func(i, j int) bool { return t.UTRs[i].Start < t.UTRs[j].Start })
		idx.ByChrom[t.Chrom] = append(idx.ByChrom[t.Chrom], t)
	}
	// Order each contig's transcript list the way htslib's regidx returns
	// overlaps: by start ascending, then end DESCENDING (longer spans first —
	// cmp_regs in regidx.c), with a deterministic ID tie-break. csq.c iterates
	// idx_tscript / idx_cds / idx_utr in this order, and csq_push keeps the
	// FIRST-staged consequence on a dedup (e.g. an intron shared by two
	// overlapping genes collapses onto the earlier-starting transcript), so the
	// scan order decides the winner. Building ByChrom from a Go map left it in
	// randomised iteration order; matching regidx makes that selection — and the
	// BCSQ consequence ordering it feeds — deterministic and upstream-faithful.
	for chrom := range idx.ByChrom {
		list := idx.ByChrom[chrom]
		sort.Slice(list, func(i, j int) bool {
			if list[i].Beg != list[j].Beg {
				return list[i].Beg < list[j].Beg
			}
			if list[i].End != list[j].End {
				return list[i].End > list[j].End
			}
			return list[i].ID < list[j].ID
		})
	}
	// Flatten transcripts-with-CDS exons into the idx_exon-ordered splice index
	// (start ascending, end descending, transcript-ID tie-break). test_splice
	// only considers transcripts that have CDS exons (the `!tr->ncds` skip in
	// csq.c — keyed off actual CDS presence, NOT the biotype coding flag), so
	// transcripts without any CDS are excluded here.
	for chrom, list := range idx.ByChrom {
		for _, t := range list {
			if len(t.CDSExons) == 0 {
				continue
			}
			for _, ex := range t.Exons {
				idx.SpliceExons[chrom] = append(idx.SpliceExons[chrom], csqExonRef{Start: ex.Start, End: ex.End, Tr: t})
			}
		}
	}
	for chrom := range idx.SpliceExons {
		ex := idx.SpliceExons[chrom]
		sort.Slice(ex, func(i, j int) bool {
			if ex[i].Start != ex[j].Start {
				return ex[i].Start < ex[j].Start
			}
			if ex[i].End != ex[j].End {
				return ex[i].End > ex[j].End
			}
			return ex[i].Tr.ID < ex[j].Tr.ID
		})
	}
	// Flatten every transcript's UTRs into the idx_utr-ordered UTR index
	// (start ascending, end descending, transcript-ID tie-break). test_utr
	// iterates idx_utr, so overlapping transcripts contribute their UTR
	// consequences in UTR-position order rather than transcript-start
	// order — matching htslib regidx cmp_regs (beg asc, end desc) and the
	// resulting INFO/BCSQ entry order.
	for chrom, list := range idx.ByChrom {
		for _, t := range list {
			for _, u := range t.UTRs {
				idx.UTRRegions[chrom] = append(idx.UTRRegions[chrom],
					csqUTRRef{Start: u.Start, End: u.End, Prime5: u.Prime5, Tr: t})
			}
		}
	}
	for chrom := range idx.UTRRegions {
		ur := idx.UTRRegions[chrom]
		sort.Slice(ur, func(i, j int) bool {
			if ur[i].Start != ur[j].Start {
				return ur[i].Start < ur[j].Start
			}
			if ur[i].End != ur[j].End {
				return ur[i].End > ur[j].End
			}
			// Tie-break by transcript ID for a deterministic order. Upstream
			// returns equal-span UTRs in htslib regidx pointer-address order
			// (its khash id2tr push order), which is not reproducible without
			// replicating that hash; the ID tie-break is the closest stable
			// approximation and matches upstream on the position-ordered
			// majority.
			return ur[i].Tr.ID < ur[j].Tr.ID
		})
	}
	return idx, nil
}

// isGeneType reports whether a GFF3 feature type denotes a gene (the
// top level of the gene/transcript/exon hierarchy). Ensembl uses
// biotype-qualified names such as "lincRNA_gene" alongside plain
// "gene", so any type ending in "gene" is accepted.
func isGeneType(t string) bool {
	return t == "gene" || strings.HasSuffix(t, "_gene")
}

// isTranscriptType reports whether a GFF3 feature type denotes a
// transcript (the second level of the gene/transcript/exon hierarchy).
func isTranscriptType(t string) bool {
	switch t {
	case "mRNA", "transcript", "lnc_RNA", "lncRNA", "lincRNA", "ncRNA", "miRNA",
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

// stripGFFPrefix removes the Ensembl-style "transcript:" / "gene:"
// namespace prefix from a GFF3 ID, matching upstream gff_id2string
// which stores and reports the bare accession.
func stripGFFPrefix(id string) string {
	if i := strings.IndexByte(id, ':'); i >= 0 {
		return id[i+1:]
	}
	return id
}

// classifyBiotype maps a raw GFF3 biotype string onto its upstream
// canonical form and coding flag, mirroring bcftools gff.c
// gff_parse_biotype (raw string -> GF_* enum) followed by
// gf_type2gff_string (enum -> canonical string) and GF_is_coding (enum
// -> coding bit). The canonical string is what upstream emits in the
// BCSQ biotype field and what its csq_push dedup compares on (upstream
// stores the enum, not the raw string), so e.g. both
// "protein_coding_CDS_not_defined" and "protein_coding" canonicalise to
// "protein_coding" and thus collapse to a single BCSQ term.
//
// The prefix rules follow upstream exactly: "protein_coding" is a
// 14-char prefix match (so protein_coding_CDS_not_defined etc. still map
// to GF_PROTEIN_CODING); "mRNA"/"MRNA" map case-insensitively to
// GF_PROTEIN_CODING; the *_pseudogene branches are tested before the
// bare coding IG_*/TR_* branches. When the biotype is unrecognised the
// raw string is returned unchanged with coding=false, matching the
// engine's earlier lenient behaviour for biotypes upstream would ignore.
func classifyBiotype(biotype string) (canonical string, coding bool) {
	if biotype == "" {
		return biotype, false
	}
	hasPrefix := func(p string) bool { return strings.HasPrefix(biotype, p) }
	hasPrefixFold := func(p string) bool {
		return len(biotype) >= len(p) && strings.EqualFold(biotype[:len(p)], p)
	}
	switch biotype[0] {
	case 'p':
		if hasPrefix("protein_coding") {
			return "protein_coding", true
		}
		if hasPrefix("polymorphic_pseudogene") {
			return "polymorphic_pseudogene", true
		}
		if hasPrefix("processed_transcript") {
			return "processed_transcript", false
		}
		if hasPrefix("processed_pseudogene") {
			return "processed_pseudogene", false
		}
		if hasPrefix("pseudogene") {
			return "pseudogene", false
		}
	case 'a':
		if hasPrefix("artifact") {
			return "artifact", false
		}
		if hasPrefix("antisense") {
			return "antisense", false
		}
		if hasPrefix("ambiguous_orf") {
			return "ambiguous_orf", false
		}
	case 'I':
		if hasPrefix("IG_pseudogene") {
			return "IG_pseudogene", false
		}
		if hasPrefix("IG_C_pseudogene") {
			return "IG_C_pseudogene", false
		}
		if hasPrefix("IG_J_pseudogene") {
			return "IG_J_pseudogene", false
		}
		if hasPrefix("IG_V_pseudogene") {
			return "IG_V_pseudogene", false
		}
		if hasPrefix("IG_C") {
			return "IG_C", true
		}
		if hasPrefix("IG_D") {
			return "IG_D", true
		}
		if hasPrefix("IG_J") {
			return "IG_J", true
		}
		if hasPrefix("IG_V") {
			return "IG_V", true
		}
		if hasPrefix("IG_LV") {
			return "IG_LV", true
		}
	case 'T':
		if hasPrefix("TR_V_pseudogene") {
			return "TR_V_pseudogene", false
		}
		if hasPrefix("TR_J_pseudogene") {
			return "TR_J_pseudogene", false
		}
		if hasPrefix("TR_C") {
			return "TR_C", true
		}
		if hasPrefix("TR_D") {
			return "TR_D", true
		}
		if hasPrefix("TR_J") {
			return "TR_J", true
		}
		if hasPrefix("TR_V") {
			return "TR_V", true
		}
	case 'M':
		if hasPrefix("Mt_tRNA_pseudogene") {
			return "MT_tRNA_pseudogene", false
		}
		if hasPrefixFold("Mt_tRNA") {
			return "MT_tRNA", false
		}
		if hasPrefixFold("Mt_rRNA") {
			return "MT_rRNA", false
		}
		if hasPrefixFold("MRNA") {
			return "protein_coding", true
		}
	case 'l':
		if hasPrefix("lincRNA") {
			return "lincRNA", false
		}
		if hasPrefix("lncRNA") {
			return "lncRNA", false
		}
	case 'm':
		if hasPrefix("macro_lncRNA") {
			return "macro_lncRNA", false
		}
		if hasPrefix("misc_RNA_pseudogene") {
			return "misc_RNA_pseudogene", false
		}
		if hasPrefix("miRNA_pseudogene") {
			return "miRNA_pseudogene", false
		}
		if hasPrefix("miRNA") {
			return "miRNA", false
		}
		if hasPrefix("misc_RNA") {
			return "misc_RNA", false
		}
		if hasPrefixFold("mRNA") {
			return "protein_coding", true
		}
	case 'r':
		if hasPrefix("rRNA") {
			return "rRNA", false
		}
		if hasPrefix("ribozyme") {
			return "ribozyme", false
		}
		if hasPrefix("retained_intron") {
			return "retained_intron", false
		}
		if hasPrefix("retrotransposed") {
			return "retrotransposed", false
		}
	case 's':
		if hasPrefix("snRNA") {
			return "snRNA", false
		}
		if hasPrefix("sRNA") {
			return "sRNA", false
		}
		if hasPrefix("scRNA") {
			return "scRNA", false
		}
		if hasPrefix("scaRNA") {
			return "scaRNA", false
		}
		if hasPrefix("snoRNA") {
			return "snoRNA", false
		}
		if hasPrefix("sense_intronic") {
			return "sense_intronic", false
		}
		if hasPrefix("sense_overlapping") {
			return "sense_overlapping", false
		}
	case 't':
		if hasPrefix("tRNA_pseudogene") {
			return "Trna_pseudogene", false
		}
		if hasPrefix("transcribed_processed_pseudogene") {
			return "transcribed_processed_pseudogene", false
		}
		if hasPrefix("transcribed_unprocessed_pseudogene") {
			return "transcribed_unprocessed_pseudogene", false
		}
		if hasPrefix("transcribed_unitary_pseudogene") {
			return "transcribed_unitary_pseudogene", false
		}
		if hasPrefix("translated_unprocessed_pseudogene") {
			return "translated_unprocessed_pseudogene", false
		}
		if hasPrefix("translated_processed_pseudogene") {
			return "translated_processed_pseudogene", false
		}
	case 'n':
		if hasPrefix("nonsense_mediated_decay") {
			return "NMD", true
		}
		if hasPrefix("non_stop_decay") {
			return "non_stop_decay", true
		}
	case 'N':
		if hasPrefix("NMD") {
			return "NMD", true
		}
	case 'k':
		if hasPrefix("known_ncrna") {
			return "known_ncRNA", false
		}
	case 'u':
		if hasPrefix("unitary_pseudogene") {
			return "unitary_pseudogene", false
		}
		if hasPrefix("unprocessed_pseudogene") {
			return "unprocessed_pseudogene", false
		}
	case 'L':
		if hasPrefix("LRG_gene") {
			return "LRG_gene", false
		}
	case '3':
		if hasPrefixFold("3prime_overlapping_ncRNA") || hasPrefixFold("3_prime_overlapping_ncRNA") {
			return "3_prime_overlapping_ncRNA", false
		}
	case 'd':
		if hasPrefix("disrupted_domain") {
			return "disrupted_domain", false
		}
	case 'v':
		if hasPrefix("vaultRNA") {
			return "vaultRNA", false
		}
	case 'b':
		if hasPrefix("bidirectional_promoter_lncRNA") {
			return "bidirectional_promoter_lncRNA", false
		}
	}
	return biotype, false
}

// transcriptPhaseInconsistency replicates upstream gff.c tscript_init_cds's
// per-CDS reading-frame sanity check. It first applies the same 5' phase
// trim upstream does (single-exon for forward strand; the multi-exon
// walk-back for reverse strand — see finalizeTranscript), then walks the
// CDS exons in transcript order verifying that each exon's GFF3 phase
// agrees with the frame implied by the cumulative CDS length. On the first
// mismatch it returns bad=true with the offending CDS's 1-based genomic
// start (badPos), the computed phase (3-phase, or 0), and the cumulative
// length so far (plen) — the exact values upstream prints. CDS exons whose
// phase is unknown (the parser records -1 / out-of-range) are skipped, as
// upstream skips CDS_PHASE_UNKN. The transcript's stored exons are not
// mutated; the check operates on local copies.
func transcriptPhaseInconsistency(t *CSQTranscript) (bad bool, badPos, phase, plen int) {
	n := len(t.CDSExons)
	if n == 0 {
		return false, 0, 0, 0
	}
	// Local (beg, len, phase) copies in genomic order (CDSExons is sorted
	// by Start ascending, so index 0 is the lowest coordinate).
	type cds struct {
		beg, ln, ph int
	}
	cs := make([]cds, n)
	for i, e := range t.CDSExons {
		ph := e.Phase
		if ph < 0 || ph > 2 {
			ph = csqPhaseUnknown
		}
		cs[i] = cds{beg: e.Start, ln: e.End - e.Start + 1, ph: ph}
	}

	if t.Strand == gff.StrandReverse {
		// 5' trim + multi-exon walk-back (upstream STRAND_REV branch).
		if cs[n-1].ph != csqPhaseUnknown {
			i := n - 1
			ph := cs[i].ph
			for i >= 0 && ph > cs[i].ln {
				ph -= cs[i].ln
				cs[i].ph = 0
				cs[i].ln = 0
				i--
			}
			cs[i].ln -= cs[i].ph
			cs[i].ph = 0
		}
		// Sanity check phase, walking 3'->5' in genomic terms (highest
		// coordinate first is the 5' end on the reverse strand).
		length := 0
		for i := n - 1; i >= 0; i-- {
			if cs[i].ph == csqPhaseUnknown {
				length += cs[i].ln
				continue
			}
			p := 0
			if cs[i].ph != 0 {
				p = 3 - cs[i].ph
			}
			if p != length%3 {
				return true, cs[i].beg, p, length
			}
			length += cs[i].ln
		}
		return false, 0, 0, 0
	}

	// Forward strand (upstream STRAND_FWD branch): single 5' exon trim.
	if cs[0].ph != csqPhaseUnknown {
		cs[0].beg += cs[0].ph
		cs[0].ln -= cs[0].ph
		cs[0].ph = 0
	}
	length := 0
	for i := 0; i < n; i++ {
		if cs[i].ph == csqPhaseUnknown {
			length += cs[i].ln
			continue
		}
		p := 0
		if cs[i].ph != 0 {
			p = 3 - cs[i].ph
		}
		if p != length%3 {
			return true, cs[i].beg, p, length
		}
		length += cs[i].ln
	}
	return false, 0, 0, 0
}

// csqPhaseUnknown marks a CDS exon whose GFF3 phase column was '.'; it
// mirrors upstream's CDS_PHASE_UNKN sentinel.
const csqPhaseUnknown = 3

// finalizeTranscript fills in derived fields after the GFF passes:
// Coding, the transcript span, the incomplete-CDS trim flags, and any
// UTR regions not explicitly present in the GFF (derived as the parts
// of full exons that fall outside the CDS span). The Coding flag is
// derived from the biotype string (upstream GF_is_coding), NOT from CDS
// presence; the CDS-derived blocks below (phase trim, UTR derivation,
// 3' trim) instead key off the actual presence of CDS exons, matching
// upstream gff.c which operates on tr->ncds there.
func finalizeTranscript(t *CSQTranscript, ref []byte) {
	// Canonicalise the biotype string to its upstream form and derive the
	// coding flag from it (upstream gff.c gff_parse_biotype +
	// gf_type2gff_string + GF_is_coding). Both the BCSQ biotype field and
	// the csq_push dedup then operate on the canonical string, so distinct
	// raw biotypes that share a canonical form (e.g.
	// protein_coding_CDS_not_defined and protein_coding) collapse exactly
	// as upstream collapses their shared GF_* enum.
	t.Biotype, t.Coding = classifyBiotype(t.Biotype)
	hasCDS := len(t.CDSExons) > 0

	// Trim the GFF3 reading-frame phase off the 5' CDS exon, mirroring
	// gff.c (the "trim non-coding start" block). After trimming, every
	// CDS exon begins exactly on a codon boundary and Phase is 0, so the
	// haplotype engine and the per-record classifier both treat the
	// spliced CDS as frame-aligned. A non-zero leading phase also marks
	// the CDS as 5' incomplete (TRIM_5PRIME).
	if hasCDS {
		if t.Strand == gff.StrandReverse {
			// The 5' end on the reverse strand is the highest-coordinate CDS
			// exon (last after the Start-ascending sort). Upstream gff.c
			// checks that the phase is not bigger than that exon's length and,
			// when it is (a 5' incomplete CDS whose leading phase spills past
			// the first exon), walks the trim back across MULTIPLE exons,
			// fully consuming each until the residual phase fits. This is the
			// reverse-strand incomplete-CDS case documented in upstream's
			// test/csq/ENST00000520868; GENCODE chr20 never triggers it.
			last := len(t.CDSExons) - 1
			ph := t.CDSExons[last].Phase
			if ph >= 1 && ph <= 2 {
				t.Trim5 = true
			}
			if ph < 0 || ph > 2 {
				ph = 0
			}
			i := last
			for i >= 0 && ph > t.CDSExons[i].End-t.CDSExons[i].Start+1 {
				ph -= t.CDSExons[i].End - t.CDSExons[i].Start + 1
				// Fully consume this exon (len -> 0) and drop out of the CDS.
				t.CDSExons[i].End = t.CDSExons[i].Start - 1
				t.CDSExons[i].Phase = 0
				i--
			}
			if i >= 0 {
				// Trim the residual exon by its OWN GFF phase (gff.c:780
				// `cds[i]->len -= cds[i]->phase`), matching the sibling check
				// in transcriptPhaseInconsistency, NOT the loop's residual ph.
				t.CDSExons[i].End -= t.CDSExons[i].Phase
				t.CDSExons[i].Phase = 0
			}
			// Drop any exons that the walk-back fully consumed (len 0).
			trimmed := t.CDSExons[:0]
			for _, e := range t.CDSExons {
				if e.End >= e.Start {
					trimmed = append(trimmed, e)
				}
			}
			t.CDSExons = trimmed
		} else {
			ph := t.CDSExons[0].Phase
			if ph >= 1 && ph <= 2 {
				t.CDSExons[0].Start += ph
				t.Trim5 = true
			}
			t.CDSExons[0].Phase = 0
		}
	}

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
	if len(t.UTRs) == 0 && hasCDS && len(t.Exons) > 0 {
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

	// Note: upstream's TRIM_5PRIME comes solely from a non-zero GFF3
	// CDS phase (handled by the phase-trim block above). A first codon
	// that is not ATG does NOT set TRIM_5PRIME; instead the haplotype
	// engine suppresses the start check locally per variant (see
	// hapInit's checkStart handling), so no extra flag is needed here.
	_ = ref

	// Mark an incomplete 3' CDS: a coding transcript whose total
	// coding length is not a multiple of three. Upstream (gff.c:844
	// `if (len%3 != 0) tr->trim |= TRIM_3PRIME`) computes len over the
	// phase-trimmed CDS exons, so the equivalent here is the sum of CDS
	// exon lengths minus the 5' reading-frame phase. Without this flag
	// an incomplete-3' transcript still runs the stop check in the
	// haplotype engine and can emit a spurious stop_lost on its last,
	// non-stop codon. csq.c:1646-1650 gates check_stop on !TRIM_3PRIME;
	// the equivalent gate lives in hapInit (csq_engine.go).
	if hasCDS {
		codingLen := -transcriptFirstPhase(t)
		for _, e := range t.CDSExons {
			codingLen += e.End - e.Start + 1
		}
		if codingLen%3 != 0 {
			t.Trim3 = true
		}
	}
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
