package bcftools

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
)

// dumpCSQGFFToFile opens path and writes the BGZF GFF dump of idx to it.
// It is the file-path wrapper used by CSQFile for --dump-gff.
func dumpCSQGFFToFile(path string, idx *CSQIndex) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("bcftools csq: --dump-gff %q: %w", path, err)
	}
	if err := DumpGFF(f, idx); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// DumpGFF writes the parsed GFF model held in idx to w as a BGZF-
// compressed, trimmed GFF3 dump, mirroring upstream gff.c's gff_dump
// (csq --dump-gff FILE). The output has five sections in upstream's
// order — genes, transcripts, CDS, UTRs, exons — with these per-line
// formats:
//
//	SEQ \t . \t gene             \t BEG \t END \t . \t STRAND \t . \t ID=GID;Name=NAME;used=N
//	SEQ \t . \t TYPE             \t BEG \t END \t . \t STRAND \t . \t ID=TID;Parent=GID;biotype=BIO;used=N
//	SEQ \t . \t CDS              \t BEG \t END \t . \t STRAND \t PHASE \t Parent=TID
//	SEQ \t . \t {five,three}_prime_UTR \t BEG \t END \t . \t STRAND \t . \t Parent=TID
//	SEQ \t . \t exon             \t BEG \t END \t . \t STRAND \t . \t Parent=TID
//
// The transcript TYPE column is "mRNA" for protein-coding transcripts
// and the transcript biotype otherwise, matching upstream's
// gf_type2gff_string mapping.
//
// Ordering note: upstream emits the gene section in hash-table order,
// which is not deterministic. This port instead sorts every section by
// (contig, begin, end, id) for reproducibility; for inputs whose genes
// are already position-ordered the output is byte-identical to upstream.
// The biotype column reproduces the input biotype verbatim rather than
// re-deriving upstream's internal biotype enum, so non-standard biotype
// spellings may differ in normalisation.
func DumpGFF(w io.Writer, idx *CSQIndex) error {
	bw := bgzf.NewWriter(w)
	var sb strings.Builder

	// Section 1: genes.
	genes := append([]*CSQGene(nil), idx.Genes...)
	sort.SliceStable(genes, func(i, j int) bool { return geneLess(genes[i], genes[j]) })
	for _, g := range genes {
		fmt.Fprintf(&sb, "%s\t.\tgene\t%d\t%d\t.\t%c\t.\tID=%s;Name=%s;used=%d\n",
			g.Chrom, g.Beg, g.End, strandChar(g.Strand), g.ID, g.Name, b2i(g.Used))
	}

	// Section 2: transcripts.
	trs := make([]*CSQTranscript, 0, len(idx.Transcripts))
	for _, t := range idx.Transcripts {
		trs = append(trs, t)
	}
	sort.SliceStable(trs, func(i, j int) bool { return tscriptLess(trs[i], trs[j]) })
	for _, t := range trs {
		fmt.Fprintf(&sb, "%s\t.\t%s\t%d\t%d\t.\t%c\t.\tID=%s;Parent=%s;biotype=%s;used=%d\n",
			t.Chrom, dumpTranscriptType(t), t.Beg, t.End, strandChar(t.Strand),
			t.ID, t.GeneID, t.Biotype, b2i(t.Used))
	}

	// Section 3: CDS.
	type cdsRow struct {
		t *CSQTranscript
		e CSQExon
	}
	var cdss []cdsRow
	for _, t := range trs {
		for _, e := range dumpCDSExons(t) {
			cdss = append(cdss, cdsRow{t, e})
		}
	}
	sort.SliceStable(cdss, func(i, j int) bool {
		return rowLess(cdss[i].t.Chrom, cdss[i].e.Start, cdss[i].e.End, cdss[i].t.ID, cdss[j].t.Chrom, cdss[j].e.Start, cdss[j].e.End, cdss[j].t.ID)
	})
	for _, r := range cdss {
		fmt.Fprintf(&sb, "%s\t.\tCDS\t%d\t%d\t.\t%c\t%s\tParent=%s\n",
			r.t.Chrom, r.e.Start, r.e.End, strandChar(r.t.Strand), phaseChar(r.e.Phase), r.t.ID)
	}

	// Section 4: UTRs.
	type utrRow struct {
		t *CSQTranscript
		u CSQUTR
	}
	var utrs []utrRow
	for _, t := range trs {
		for _, u := range t.UTRs {
			utrs = append(utrs, utrRow{t, u})
		}
	}
	sort.SliceStable(utrs, func(i, j int) bool {
		return rowLess(utrs[i].t.Chrom, utrs[i].u.Start, utrs[i].u.End, utrs[i].t.ID, utrs[j].t.Chrom, utrs[j].u.Start, utrs[j].u.End, utrs[j].t.ID)
	})
	for _, r := range utrs {
		kind := "three_prime_UTR"
		if r.u.Prime5 {
			kind = "five_prime_UTR"
		}
		fmt.Fprintf(&sb, "%s\t.\t%s\t%d\t%d\t.\t%c\t.\tParent=%s\n",
			r.t.Chrom, kind, r.u.Start, r.u.End, strandChar(r.t.Strand), r.t.ID)
	}

	// Section 5: exons.
	type exonRow struct {
		t *CSQTranscript
		e CSQExon
	}
	var exons []exonRow
	for _, t := range trs {
		for _, e := range t.Exons {
			exons = append(exons, exonRow{t, e})
		}
	}
	sort.SliceStable(exons, func(i, j int) bool {
		return rowLess(exons[i].t.Chrom, exons[i].e.Start, exons[i].e.End, exons[i].t.ID, exons[j].t.Chrom, exons[j].e.Start, exons[j].e.End, exons[j].t.ID)
	})
	for _, r := range exons {
		fmt.Fprintf(&sb, "%s\t.\texon\t%d\t%d\t.\t%c\t.\tParent=%s\n",
			r.t.Chrom, r.e.Start, r.e.End, strandChar(r.t.Strand), r.t.ID)
	}

	if _, err := io.WriteString(bw, sb.String()); err != nil {
		bw.Close()
		return fmt.Errorf("bcftools csq: --dump-gff write: %w", err)
	}
	return bw.Close()
}

// dumpCDSExons returns a transcript's CDS exons with upstream's
// TRIM_3PRIME adjustment applied (the TRIM_5PRIME phase trim is already
// baked into t.CDSExons by finalizeTranscript). When the total coding
// length is not a multiple of three, upstream gff.c (tscript_init_cds)
// trims the surplus 1-2 bases off the 3' end: the last exon (highest
// coord) on the forward strand, the first exon (lowest coord) on the
// reverse strand, spilling into adjacent exons if a CDS exon is shorter
// than the surplus. The engine keeps the untrimmed coordinates and
// gates on t.Trim3 instead, so this trim is reproduced here only for the
// dump's coordinate columns.
func dumpCDSExons(t *CSQTranscript) []CSQExon {
	out := make([]CSQExon, len(t.CDSExons))
	copy(out, t.CDSExons)
	if !t.Coding {
		return out
	}
	total := 0
	for _, e := range out {
		total += e.End - e.Start + 1
	}
	rem := total % 3
	if rem == 0 {
		return out
	}
	if t.Strand == gff.StrandReverse {
		// Trim the low-coordinate (3') end forward.
		for i := 0; i < len(out) && rem > 0; i++ {
			cdsLen := out[i].End - out[i].Start + 1
			d := rem
			if cdsLen < d {
				d = cdsLen
			}
			out[i].Start += d
			rem -= d
		}
	} else {
		// Trim the high-coordinate (3') end backward.
		for i := len(out) - 1; i >= 0 && rem > 0; i-- {
			cdsLen := out[i].End - out[i].Start + 1
			d := rem
			if cdsLen < d {
				d = cdsLen
			}
			out[i].End -= d
			rem -= d
		}
	}
	return out
}

// dumpTranscriptType returns the GFF3 type column for a transcript in
// the --dump-gff output: "mRNA" for protein-coding transcripts, the
// transcript biotype otherwise. This mirrors upstream gff_dump's
// `tr->type==GF_PROTEIN_CODING ? "mRNA" : gf_type2gff_string(tr->type)`.
func dumpTranscriptType(t *CSQTranscript) string {
	if t.Biotype == "protein_coding" {
		return "mRNA"
	}
	return t.Biotype
}

// geneLess orders genes by (contig, begin, end, id).
func geneLess(a, b *CSQGene) bool {
	return rowLess(a.Chrom, a.Beg, a.End, a.ID, b.Chrom, b.Beg, b.End, b.ID)
}

// tscriptLess orders transcripts by (contig, begin, end, id).
func tscriptLess(a, b *CSQTranscript) bool {
	return rowLess(a.Chrom, a.Beg, a.End, a.ID, b.Chrom, b.Beg, b.End, b.ID)
}

// rowLess is the shared (contig, begin, end, id) comparator used to sort
// each --dump-gff section deterministically.
func rowLess(ca string, ba, ea int, ia string, cb string, bb, eb int, ib string) bool {
	if ca != cb {
		return ca < cb
	}
	if ba != bb {
		return ba < bb
	}
	if ea != eb {
		return ea < eb
	}
	return ia < ib
}

// strandChar renders a GFF strand as the dump's single-character column.
func strandChar(s gff.Strand) byte {
	switch s {
	case gff.StrandForward:
		return '+'
	case gff.StrandReverse:
		return '-'
	default:
		return '.'
	}
}

// phaseChar renders a CDS phase as upstream gff_dump does: '.' for
// unknown (encoded as 3 / out of range), otherwise the digit '0'..'2'.
func phaseChar(p int) string {
	if p < 0 || p > 2 {
		return "."
	}
	return string(rune('0' + p))
}

// b2i maps a bool to upstream's used= integer (1 or 0).
func b2i(b bool) int {
	if b {
		return 1
	}
	return 0
}
