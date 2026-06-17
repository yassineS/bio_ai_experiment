package fixtures

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"math/rand"
	"os"
	"sort"
	"strings"
)

// illuminaAdapter is the canonical Illumina TruSeq 3' adapter prefix that the
// trimming tools (sickle, skewer, fastp, prinseq) know how to remove. It is
// appended (truncated to fit) to a fraction of reads so adapter-trimming has
// real work to do.
const illuminaAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

// phred33 converts a Phred quality (0..40) to its Sanger/Illumina-1.8 ASCII
// character (offset 33).
func phred33(q int) byte {
	if q < 0 {
		q = 0
	}
	if q > 40 {
		q = 40
	}
	return byte('!' + q)
}

// randRead synthesises one realistic read: a body of random bases at high
// quality, an optional adapter contaminating the 3' end (on adapterFrac of
// reads), and a low-quality / N tail (on tailFrac of reads). The returned
// length varies around p.FastqReadLen so length filters have something to do.
// It draws every value from rng so the whole fixture is reproducible.
func randRead(p Params, rng *rand.Rand) (seq, qual []byte) {
	mean := p.FastqReadLen
	// Vary length by +/- 10% of the mean (at least 1 base of jitter).
	jitter := mean / 10
	if jitter < 1 {
		jitter = 1
	}
	rl := mean - jitter + rng.Intn(2*jitter+1)
	if rl < 20 {
		rl = 20
	}
	seq = make([]byte, rl)
	qual = make([]byte, rl)
	for i := 0; i < rl; i++ {
		seq[i] = bases[rng.Intn(4)]
		// High base quality (Q28..Q40) for the bulk of the read.
		qual[i] = phred33(28 + rng.Intn(13))
	}
	// 25% of reads carry adapter contamination starting somewhere past the
	// halfway point.
	if rng.Intn(4) == 0 {
		start := rl/2 + rng.Intn(rl/4+1)
		for i := 0; start+i < rl && i < len(illuminaAdapter); i++ {
			seq[start+i] = illuminaAdapter[i]
		}
	}
	// 30% of reads carry a degraded 3' tail: falling quality and the
	// occasional N, mirroring real sequencer end-of-read behaviour so
	// quality-trim and N-handling flags exercise.
	if rng.Intn(10) < 3 {
		tail := rng.Intn(rl/4 + 1)
		for i := rl - tail; i < rl; i++ {
			qual[i] = phred33(rng.Intn(12)) // Q0..Q11
			if rng.Intn(5) == 0 {
				seq[i] = 'N'
			}
		}
	}
	return seq, qual
}

// writeFastqSE writes a single-end FASTQ to plain and (when gzPath != "")
// gzip-compressed files, returning after both are flushed. Records are named
// read000001/1 .. and drawn deterministically from rng.
func writeFastqSE(plainPath, gzPath string, p Params, rng *rand.Rand) error {
	pf, err := os.Create(plainPath)
	if err != nil {
		return err
	}
	defer pf.Close()
	pw := bufio.NewWriterSize(pf, 1<<20)

	var gf *os.File
	var gz *gzip.Writer
	var gw *bufio.Writer
	if gzPath != "" {
		gf, err = os.Create(gzPath)
		if err != nil {
			return err
		}
		defer gf.Close()
		gz = gzip.NewWriter(gf)
		gw = bufio.NewWriterSize(gz, 1<<20)
	}

	for i := 0; i < p.FastqReads; i++ {
		seq, qual := randRead(p, rng)
		rec := fmt.Sprintf("@read%06d/1\n%s\n+\n%s\n", i, seq, qual)
		if _, err := pw.WriteString(rec); err != nil {
			return err
		}
		if gw != nil {
			if _, err := gw.WriteString(rec); err != nil {
				return err
			}
		}
	}
	if err := pw.Flush(); err != nil {
		return err
	}
	if gw != nil {
		if err := gw.Flush(); err != nil {
			return err
		}
		if err := gz.Close(); err != nil {
			return err
		}
	}
	return nil
}

// writeFastqPE writes a paired-end FASTQ set: <r1Path> and <r2Path> with
// matching record names (readNNNNNN/1 and /2). Mate 2 sequences are independent
// reads (not reverse complements) — sufficient for the trimming tools, which
// process the two files in lockstep. Records are deterministic from rng.
func writeFastqPE(r1Path, r2Path string, p Params, rng *rand.Rand) error {
	f1, err := os.Create(r1Path)
	if err != nil {
		return err
	}
	defer f1.Close()
	f2, err := os.Create(r2Path)
	if err != nil {
		return err
	}
	defer f2.Close()
	w1 := bufio.NewWriterSize(f1, 1<<20)
	w2 := bufio.NewWriterSize(f2, 1<<20)

	for i := 0; i < p.FastqReads; i++ {
		s1, q1 := randRead(p, rng)
		s2, q2 := randRead(p, rng)
		if _, err := fmt.Fprintf(w1, "@read%06d/1\n%s\n+\n%s\n", i, s1, q1); err != nil {
			return err
		}
		if _, err := fmt.Fprintf(w2, "@read%06d/2\n%s\n+\n%s\n", i, s2, q2); err != nil {
			return err
		}
	}
	if err := w1.Flush(); err != nil {
		return err
	}
	return w2.Flush()
}

// writeGFF writes a valid GFF3 over the same contigs as the rest of the fixture
// set. Each gene locus expands into a gene, an mRNA, and a few exon/CDS child
// rows with proper Parent= attributes, sorted by (contig order, start) so it is
// valid input to interval tools and bcftools csq. Coordinates are 1-based,
// inclusive, matching the GFF3 spec.
func writeGFF(path string, contigs []contig, p Params, rng *rand.Rand) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriterSize(f, 1<<20)

	if _, err := fmt.Fprintln(w, "##gff-version 3"); err != nil {
		return err
	}
	for _, c := range contigs {
		if _, err := fmt.Fprintf(w, "##sequence-region %s 1 %d\n", c.Name, c.Len); err != nil {
			return err
		}
	}

	type gene struct {
		ci, start, end int
		strand         byte
	}
	genes := make([]gene, 0, p.Genes)
	for i := 0; i < p.Genes; i++ {
		ci := rng.Intn(len(contigs))
		span := 500 + rng.Intn(4000)
		if span >= contigs[ci].Len {
			span = contigs[ci].Len / 2
		}
		start := 1 + rng.Intn(contigs[ci].Len-span)
		strand := byte('+')
		if rng.Intn(2) == 0 {
			strand = '-'
		}
		genes = append(genes, gene{ci, start, start + span, strand})
	}
	sort.Slice(genes, func(i, j int) bool {
		if genes[i].ci != genes[j].ci {
			return genes[i].ci < genes[j].ci
		}
		return genes[i].start < genes[j].start
	})

	for gi, g := range genes {
		c := contigs[g.ci]
		gid := fmt.Sprintf("gene%05d", gi)
		tid := fmt.Sprintf("mRNA%05d", gi)
		fmt.Fprintf(w, "%s\tparity\tgene\t%d\t%d\t.\t%c\t.\tID=%s;Name=%s\n",
			c.Name, g.start, g.end, g.strand, gid, gid)
		fmt.Fprintf(w, "%s\tparity\tmRNA\t%d\t%d\t.\t%c\t.\tID=%s;Parent=%s\n",
			c.Name, g.start, g.end, g.strand, tid, gid)
		// 1..3 exons inside the gene span, each with a CDS row.
		nExon := 1 + rng.Intn(3)
		total := g.end - g.start
		seg := total / (nExon + 1)
		if seg < 1 {
			seg = 1
		}
		for e := 0; e < nExon; e++ {
			es := g.start + e*seg
			ee := es + seg - 1
			if ee >= g.end {
				ee = g.end
			}
			phase := (e * 2) % 3
			fmt.Fprintf(w, "%s\tparity\texon\t%d\t%d\t.\t%c\t.\tID=%s.exon%d;Parent=%s\n",
				c.Name, es, ee, g.strand, tid, e+1, tid)
			fmt.Fprintf(w, "%s\tparity\tCDS\t%d\t%d\t.\t%c\t%d\tID=%s.cds%d;Parent=%s\n",
				c.Name, es, ee, g.strand, phase, tid, e+1, tid)
		}
	}
	return w.Flush()
}

// writeMultiSampleVCF writes a valid VCF 4.2 with p.MultiSamples sample columns
// of biallelic SNPs, so vcftools modes that need more than one sample
// (relatedness, het, LD) have meaningful input. It reuses the contig sequences
// for REF bases and is sorted in coordinate order so it can be bgzipped and
// tabixed. Output is deterministic from rng.
func writeMultiSampleVCF(path string, contigs []contig, p Params, rng *rand.Rand) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriterSize(f, 1<<20)

	fmt.Fprintln(w, "##fileformat=VCFv4.2")
	fmt.Fprintln(w, "##source=parity-pipeline-fixtures")
	for _, c := range contigs {
		fmt.Fprintf(w, "##contig=<ID=%s,length=%d>\n", c.Name, c.Len)
	}
	fmt.Fprintln(w, `##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">`)
	fmt.Fprintln(w, `##FILTER=<ID=PASS,Description="All filters passed">`)
	fmt.Fprintln(w, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
	fmt.Fprintln(w, `##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">`)
	nS := p.MultiSamples
	if nS < 2 {
		nS = 2
	}
	var hdr strings.Builder
	hdr.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for s := 0; s < nS; s++ {
		fmt.Fprintf(&hdr, "\tsample%d", s+1)
	}
	fmt.Fprintln(w, hdr.String())

	type rec struct{ ci, pos int }
	// Use the same number of variants as the single-sample VCF for size parity.
	recs := make([]rec, 0, p.Variants)
	for i := 0; i < p.Variants; i++ {
		ci := rng.Intn(len(contigs))
		pos := 1 + rng.Intn(contigs[ci].Len-10)
		recs = append(recs, rec{ci, pos})
	}
	sort.Slice(recs, func(i, j int) bool {
		if recs[i].ci != recs[j].ci {
			return recs[i].ci < recs[j].ci
		}
		return recs[i].pos < recs[j].pos
	})

	gts := []string{"0/0", "0/1", "1/1"}
	for i, r := range recs {
		c := contigs[r.ci]
		ref := string(c.Seq[r.pos-1])
		alt := string(bases[(strings.IndexByte(bases, c.Seq[r.pos-1])+1+rng.Intn(3))%4])
		qual := 5 + rng.Intn(95)
		var line strings.Builder
		fmt.Fprintf(&line, "%s\t%d\trs%d\t%s\t%s\t%d\tPASS\tDP=%d\tGT:DP",
			c.Name, r.pos, i+1, ref, alt, qual, 10+rng.Intn(90))
		for s := 0; s < nS; s++ {
			fmt.Fprintf(&line, "\t%s:%d", gts[rng.Intn(len(gts))], 5+rng.Intn(60))
		}
		fmt.Fprintln(w, line.String())
	}
	return w.Flush()
}
