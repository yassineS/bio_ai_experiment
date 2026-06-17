package fixtures

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
)

const bases = "ACGT"

// writeFasta writes a multi-contig reference with 60-base lines, filling each
// contig's sequence with seeded random bases. The sequences are stored back on
// the contig structs so SAM/VCF generation can reference real reference bases.
func writeFasta(path string, contigs []contig, rng *rand.Rand) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriterSize(f, 1<<20)
	defer w.Flush()

	for ci := range contigs {
		c := &contigs[ci]
		c.Seq = make([]byte, c.Len)
		for i := 0; i < c.Len; i++ {
			c.Seq[i] = bases[rng.Intn(4)]
		}
		if _, err := fmt.Fprintf(w, ">%s\n", c.Name); err != nil {
			return err
		}
		for i := 0; i < c.Len; i += 60 {
			end := i + 60
			if end > c.Len {
				end = c.Len
			}
			if _, err := w.Write(c.Seq[i:end]); err != nil {
				return err
			}
			if err := w.WriteByte('\n'); err != nil {
				return err
			}
		}
	}
	return nil
}

// writeSAM writes an unsorted SAM with a header derived from the contigs and
// Reads single-end records of length ReadLen drawn from the reference (with a
// low seeded mutation rate so MD/NM would be meaningful and quality strings are
// realistic). The records are intentionally emitted out of coordinate order so
// `samtools sort` has real work to do.
func writeSAM(path string, contigs []contig, p Params, rng *rand.Rand) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriterSize(f, 1<<20)
	defer w.Flush()

	if _, err := fmt.Fprintln(w, "@HD\tVN:1.6\tSO:unsorted"); err != nil {
		return err
	}
	for _, c := range contigs {
		if _, err := fmt.Fprintf(w, "@SQ\tSN:%s\tLN:%d\n", c.Name, c.Len); err != nil {
			return err
		}
	}
	if _, err := fmt.Fprintln(w, "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA"); err != nil {
		return err
	}

	rl := p.ReadLen
	seq := make([]byte, rl)
	qual := make([]byte, rl)
	for i := 0; i < p.Reads; i++ {
		c := contigs[rng.Intn(len(contigs))]
		maxStart := c.Len - rl
		if maxStart <= 0 {
			continue
		}
		pos0 := rng.Intn(maxStart) // 0-based
		copy(seq, c.Seq[pos0:pos0+rl])
		// Sprinkle a few mismatches.
		for m := 0; m < rl/40; m++ {
			j := rng.Intn(rl)
			seq[j] = bases[rng.Intn(4)]
		}
		for j := 0; j < rl; j++ {
			qual[j] = byte('!' + rng.Intn(40)) // Phred 0..39
		}
		flag := 0
		if rng.Intn(2) == 0 {
			flag |= 0x10 // reverse strand
		}
		mapq := 20 + rng.Intn(40)
		name := fmt.Sprintf("read%07d", i)
		// QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
		if _, err := fmt.Fprintf(w, "%s\t%d\t%s\t%d\t%d\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n",
			name, flag, c.Name, pos0+1, mapq, rl, seq, qual); err != nil {
			return err
		}
	}
	return nil
}

// writeVCF writes a valid VCF 4.2 with sorted records across the contigs, a
// single sample column, and biallelic SNPs/indels with realistic INFO/FORMAT
// fields. Records are emitted in coordinate order (contig order matches the
// header) so tabix indexing succeeds.
func writeVCF(path string, contigs []contig, p Params, rng *rand.Rand) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriterSize(f, 1<<20)
	defer w.Flush()

	fmt.Fprintln(w, "##fileformat=VCFv4.2")
	fmt.Fprintln(w, "##source=parity-pipeline-fixtures")
	for _, c := range contigs {
		fmt.Fprintf(w, "##contig=<ID=%s,length=%d>\n", c.Name, c.Len)
	}
	fmt.Fprintln(w, `##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">`)
	fmt.Fprintln(w, `##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">`)
	fmt.Fprintln(w, `##FILTER=<ID=PASS,Description="All filters passed">`)
	fmt.Fprintln(w, `##FILTER=<ID=q10,Description="Quality below 10">`)
	fmt.Fprintln(w, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
	fmt.Fprintln(w, `##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">`)
	fmt.Fprintln(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1")

	// Distribute variants across contigs, then sort per contig by position.
	type rec struct {
		ci, pos int
	}
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

	// All generated sites are biallelic (one REF, one ALT), so only genotypes
	// that index allele 0 or 1 are valid. A "1/2" genotype would reference a
	// non-existent second ALT allele and is undefined behaviour in downstream
	// tools (vcftools mis-counts N_CHR / emits -nan), so it is deliberately
	// excluded.
	gts := []string{"0/1", "1/1", "0/0"}
	for i, r := range recs {
		c := contigs[r.ci]
		ref := string(c.Seq[r.pos-1])
		alt := string(bases[(strings.IndexByte(bases, c.Seq[r.pos-1])+1+rng.Intn(3))%4])
		// Occasionally make it an insertion.
		if rng.Intn(8) == 0 {
			alt = ref + string(bases[rng.Intn(4)])
		}
		qual := 5 + rng.Intn(95)
		filt := "PASS"
		if qual < 10 {
			filt = "q10"
		}
		dp := 10 + rng.Intn(90)
		af := float64(1+rng.Intn(99)) / 100.0
		gt := gts[rng.Intn(len(gts))]
		// Format AF with the shortest representation so it round-trips
		// byte-identically through tools (like bcftools) that re-print floats
		// using shortest-decimal formatting; a fixed %.3f would emit trailing
		// zeros that upstream normalises away, producing spurious divergence.
		afStr := strconv.FormatFloat(af, 'g', -1, 64)
		fmt.Fprintf(w, "%s\t%d\trs%d\t%s\t%s\t%d\t%s\tDP=%d;AF=%s\tGT:DP\t%s:%d\n",
			c.Name, r.pos, i+1, ref, alt, qual, filt, dp, afStr, gt, dp)
	}
	return nil
}

// writeBED writes a plain BED6 file, a BED12 file, and a genome (chrom sizes)
// file, all over the same coordinate space as the reference. Records are
// sorted by (chrom-order, start) so they are valid input to interval tools
// that expect sorted input.
func writeBED(bedPath, bed12Path, genomePath string, contigs []contig, p Params, rng *rand.Rand) error {
	// genome file
	gf, err := os.Create(genomePath)
	if err != nil {
		return err
	}
	for _, c := range contigs {
		fmt.Fprintf(gf, "%s\t%d\n", c.Name, c.Len)
	}
	gf.Close()

	type ivl struct {
		ci, start, end int
		strand         byte
	}
	ivls := make([]ivl, 0, p.Intervals)
	for i := 0; i < p.Intervals; i++ {
		ci := rng.Intn(len(contigs))
		span := 100 + rng.Intn(2000)
		if span >= contigs[ci].Len {
			span = contigs[ci].Len / 2
		}
		start := rng.Intn(contigs[ci].Len - span)
		strand := byte('+')
		if rng.Intn(2) == 0 {
			strand = '-'
		}
		ivls = append(ivls, ivl{ci, start, start + span, strand})
	}
	sort.Slice(ivls, func(i, j int) bool {
		if ivls[i].ci != ivls[j].ci {
			return ivls[i].ci < ivls[j].ci
		}
		return ivls[i].start < ivls[j].start
	})

	bf, err := os.Create(bedPath)
	if err != nil {
		return err
	}
	w := bufio.NewWriterSize(bf, 1<<20)
	b12, err := os.Create(bed12Path)
	if err != nil {
		bf.Close()
		return err
	}
	w12 := bufio.NewWriterSize(b12, 1<<20)

	for i, v := range ivls {
		name := fmt.Sprintf("feat%d", i)
		score := rng.Intn(1000)
		c := contigs[v.ci]
		fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%d\t%c\n", c.Name, v.start, v.end, name, score, v.strand)

		// BED12: a couple of blocks within the interval.
		blockCount := 1 + rng.Intn(2)
		var sizes, starts []string
		if blockCount == 1 {
			sizes = []string{fmt.Sprintf("%d", v.end-v.start)}
			starts = []string{"0"}
		} else {
			half := (v.end - v.start) / 3
			sizes = []string{fmt.Sprintf("%d", half), fmt.Sprintf("%d", half)}
			starts = []string{"0", fmt.Sprintf("%d", (v.end-v.start)-half)}
		}
		fmt.Fprintf(w12, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t0,0,0\t%d\t%s\t%s\n",
			c.Name, v.start, v.end, name, score, v.strand, v.start, v.end,
			blockCount, strings.Join(sizes, ",")+",", strings.Join(starts, ",")+",")
	}
	if err := w.Flush(); err != nil {
		bf.Close()
		b12.Close()
		return err
	}
	if err := w12.Flush(); err != nil {
		bf.Close()
		b12.Close()
		return err
	}
	bf.Close()
	b12.Close()
	return nil
}
