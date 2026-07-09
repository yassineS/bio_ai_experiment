package realbench

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// deriveInputs synthesises the extra bed* inputs that a plain BED3 (the NIST
// high-confidence intervals) cannot satisfy, and returns a copy of in with the
// derived paths filled. The synthetics are written into dir (the run's scratch
// dir) so they exist on the real run too, and they are derived DETERMINISTICALLY
// from the real inputs so both ours and upstream see byte-identical inputs:
//
//   - BED4: the BED3 with a synthetic name column (region_N) appended. Feeds
//     the bed* subcommands that require a name field (tobam -> QNAME).
//   - Window: consecutive BED records paired into an 8-field
//     "chrom start end name" x2 line (the shape `bedtools window` emits),
//     for `bedoverlap -cols 2,3,6,7`.
//   - BEDPE: consecutive BED records paired into a 10-field BEDPE line
//     (chrom1 s1 e1 chrom2 s2 e2 name score strand1 strand2), for
//     pairtopair / pairtobed.
//   - BedGraph: a 4-column BedGraph (chrom start end value) derived from the
//     BED3, for `bedtools unionbedg` (which SIGABRTs on a bare BED3).
//   - SampleRename: a one-line sample-rename file for `bcftools reheader -s`.
//
// The samtoolsBin (a resolved samtools binary; upstream preferred, else ours)
// additionally drives two PREREQUISITE BAM transforms, produced ONCE and fed to
// BOTH sides identically so the comparison stays fair:
//
//   - NameBAM: `samtools sort -n` of the BAM — the name-collated input that
//     upstream `samtools fixmate` requires (it errors on coord-sorted input).
//   - FixmateBAM: `sort -n | fixmate -m | sort` — the markdup-ready input that
//     upstream `samtools markdup` requires (needs ms + MC then coord order).
//
// When a required real input is empty (or samtoolsBin is empty), the
// corresponding derivation is skipped and the dependent cells SKIP. BED-shape
// synthesis errors are returned so the caller can surface them; the BAM
// transforms are best-effort (a failure leaves the derived path empty).
func deriveInputs(in Inputs, dir, samtoolsBin string) (Inputs, error) {
	// Plain-FASTQ derivation: decompress Fastq1 into the scratch dir so the
	// prinseq cells (both ours and upstream) read a plain FASTQ. prinseq-lite.pl
	// 0.20.4 cannot read gzip, so feeding it the bgzipped R1 yields no output on
	// either side. This is independent of the BED synthesis below; a failure is
	// non-fatal (the dependent cells SKIP because FastqPlain stays empty).
	if in.Fastq1 != "" {
		plain := filepath.Join(dir, "derived.plain.fastq")
		if err := decompressFastq(in.Fastq1, plain); err == nil {
			in.FastqPlain = plain
		}
	}

	// Prerequisite BAM transforms for the samtools fixmate/markdup cells. These
	// are best-effort: a failure (or an unresolved samtools) leaves the derived
	// path empty and the dependent cell SKIPs rather than aborting the run. The
	// SAME derived file feeds both ours and upstream, so the comparison is fair.
	if in.BAM != "" && samtoolsBin != "" {
		nameBAM := filepath.Join(dir, "derived.namecollated.bam")
		if err := deriveNameCollatedBAM(samtoolsBin, in.BAM, nameBAM, dir); err == nil {
			in.NameBAM = nameBAM
		}
		fixmateBAM := filepath.Join(dir, "derived.fixmate.bam")
		if err := deriveMarkdupReadyBAM(samtoolsBin, in.BAM, fixmateBAM, dir); err == nil {
			in.FixmateBAM = fixmateBAM
		}
	}

	// One-line sample-rename map for `bcftools reheader -s`. Deterministic and
	// independent of the BED synthesis below; only needs a VCF to be relevant.
	if in.VCF != "" {
		rename := filepath.Join(dir, "derived.rename.txt")
		if err := os.WriteFile(rename, []byte("RB_SAMPLE\n"), 0o644); err == nil {
			in.SampleRename = rename
		}
	}

	// csq-normalised GFF: inject the `biotype=` attribute upstream `bcftools
	// csq` needs (upstream exits 255 without it on a bare GENCODE GFF3). The
	// same plain-text normalised GFF feeds both ours and upstream so the csq
	// cell is a fair byte-exact parity comparison. Best-effort: a failure
	// leaves NormGFF empty and the csq cell SKIPs (NeedNormGFF).
	if in.GFF != "" {
		normGFF := filepath.Join(dir, "derived.csqnorm.gff")
		if err := normaliseGFFForCsq(in.GFF, normGFF); err == nil {
			in.NormGFF = normGFF
		}
	}

	if in.BED == "" {
		return in, nil
	}
	recs, err := readBED3(in.BED)
	if err != nil {
		return in, fmt.Errorf("reading BED %s for synthesis: %w", in.BED, err)
	}
	if len(recs) == 0 {
		return in, nil
	}

	bed4 := filepath.Join(dir, "derived.bed4")
	if err := writeBED4(bed4, recs); err != nil {
		return in, err
	}
	in.BED4 = bed4

	bedgraph := filepath.Join(dir, "derived.bedgraph")
	if err := writeBedGraph(bedgraph, recs); err != nil {
		return in, err
	}
	in.BedGraph = bedgraph

	win := filepath.Join(dir, "derived.window.bed")
	if err := writeWindow(win, recs); err != nil {
		return in, err
	}
	in.Window = win

	bedpe := filepath.Join(dir, "derived.bedpe")
	if err := writeBEDPE(bedpe, recs); err != nil {
		return in, err
	}
	in.BEDPE = bedpe

	return in, nil
}

// decompressFastq reads src (transparently gzip/bgzip-decoded via iohelper) and
// writes the plain payload to dst. It is used to give the prinseq cells a plain
// FASTQ that prinseq-lite.pl 0.20.4 can actually read. When src is already
// plain, iohelper passes the bytes through unchanged, so dst is a verbatim copy.
func decompressFastq(src, dst string) error {
	r, err := iohelper.OpenReader(src)
	if err != nil {
		return err
	}
	defer r.Close()
	f, err := os.Create(dst)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	if _, err := io.Copy(w, r); err != nil {
		return err
	}
	return w.Flush()
}

// normaliseGFFForCsq reads src (transparently gzip/bgzip-decoded) and writes a
// plain GFF3 to dst with the `biotype=` attribute injected on gene / transcript
// features that lack it, mirroring reference_code/bcftools/misc/gff2gff. This is
// the fix upstream `bcftools csq` requires to parse a bare GENCODE GFF3 (it
// derives the transcript biotype from `transcript_type`, falling back to the
// parent gene's `gene_type`; genes fall back to `gene_type`). Only the
// `biotype=` injection is reproduced — the ID/Parent/Name back-fill in gff2gff
// is unnecessary for GENCODE, which already carries those. `bcftools csq`
// accepts a plain (uncompressed) GFF, so no bgzip/tabix step is needed.
func normaliseGFFForCsq(src, dst string) error {
	r, err := iohelper.OpenReader(src)
	if err != nil {
		return err
	}
	defer r.Close()

	// First pass over a buffered copy is avoided by remembering gene biotypes
	// as we go: GENCODE emits a gene before its transcripts, so a transcript's
	// parent gene_type is already known. We still fall back to the
	// transcript's own gene_type attribute (GENCODE carries it) when the parent
	// map misses, so a single streaming pass suffices.
	geneBiotype := map[string]string{}

	f, err := os.Create(dst)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 64*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "#") {
			if _, err := fmt.Fprintln(w, line); err != nil {
				return err
			}
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 9 {
			if _, err := fmt.Fprintln(w, line); err != nil {
				return err
			}
			continue
		}
		attrs := cols[8]
		typ := cols[2]
		isGene := typ == "gene"
		isTranscript := typ == "transcript" ||
			(gffAttr(attrs, "Parent") != "" && gffAttr(attrs, "biotype") == "")
		if !isGene && !isTranscript {
			if _, err := fmt.Fprintln(w, line); err != nil {
				return err
			}
			continue
		}
		if isGene {
			bt := gffAttr(attrs, "biotype")
			if bt == "" {
				bt = gffAttr(attrs, "gene_type")
			}
			if bt == "" {
				bt = gffAttr(attrs, "gene_biotype")
			}
			if id := gffAttr(attrs, "ID"); id != "" && bt != "" {
				geneBiotype[id] = bt
			}
			if gffAttr(attrs, "biotype") == "" && bt != "" {
				attrs += ";biotype=" + bt
			}
		} else { // transcript
			bt := gffAttr(attrs, "biotype")
			if bt == "" {
				bt = gffAttr(attrs, "transcript_type")
			}
			if bt == "" {
				bt = gffAttr(attrs, "transcript_biotype")
			}
			if bt == "" {
				if p := gffAttr(attrs, "Parent"); p != "" {
					bt = geneBiotype[p]
				}
			}
			if bt == "" {
				bt = gffAttr(attrs, "gene_type")
			}
			if gffAttr(attrs, "biotype") == "" && bt != "" {
				attrs += ";biotype=" + bt
			}
		}
		cols[8] = attrs
		if _, err := fmt.Fprintln(w, strings.Join(cols, "\t")); err != nil {
			return err
		}
	}
	if err := sc.Err(); err != nil {
		return err
	}
	return w.Flush()
}

// gffAttr extracts the value of a `key=value` attribute from a GFF3 column-9
// attribute string (semicolon-separated). It returns "" when the key is absent.
func gffAttr(attrs, key string) string {
	for _, kv := range strings.Split(attrs, ";") {
		kv = strings.TrimSpace(kv)
		if strings.HasPrefix(kv, key+"=") {
			return kv[len(key)+1:]
		}
	}
	return ""
}

// bedRec is a single parsed BED3 interval.
type bedRec struct {
	chrom      string
	start, end string // kept as strings so the coordinates round-trip verbatim
}

// readBED3 reads the first three columns of every non-comment, non-track line of
// a (optionally gzipped) BED file. Only chrom/start/end are retained; extra
// columns are ignored so the synthesis works off a strict BED3 view.
func readBED3(path string) ([]bedRec, error) {
	f, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var recs []bedRec
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || strings.HasPrefix(line, "#") ||
			strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") {
			continue
		}
		f := strings.Split(line, "\t")
		if len(f) < 3 {
			continue
		}
		recs = append(recs, bedRec{chrom: f[0], start: f[1], end: f[2]})
	}
	return recs, sc.Err()
}

// deriveNameCollatedBAM produces a name-collated (queryname-grouped) BAM from
// src via `samtools sort -n`, the input upstream `samtools fixmate` requires.
// tmp is the scratch dir for samtool's intermediate spill files. The output is
// deterministic given src, so both ours and upstream digest the same fixmate
// input. Returns an error (leaving no output) when the transform fails.
func deriveNameCollatedBAM(samtoolsBin, src, dst, tmp string) error {
	prefix := filepath.Join(tmp, "rb-namesort")
	_, err := runOnce(samtoolsBin,
		[]string{"sort", "-n", "-T", prefix, "-o", dst, src}, "", "", nil, nil)
	return err
}

// deriveMarkdupReadyBAM produces the input upstream `samtools markdup` requires:
// name-sort, `fixmate -m` (adds the ms mate-score and MC tags markdup consumes),
// then coordinate re-sort. It is produced ONCE with the resolved samtools and
// fed to both sides, so ours-vs-upstream markdup compares like for like. tmp
// holds the intermediate name-sorted/fixmate'd BAMs and samtools' spill files.
func deriveMarkdupReadyBAM(samtoolsBin, src, dst, tmp string) error {
	nameSorted := filepath.Join(tmp, "rb-md-namesort.bam")
	if _, err := runOnce(samtoolsBin,
		[]string{"sort", "-n", "-T", filepath.Join(tmp, "rb-md-ns"), "-o", nameSorted, src},
		"", "", nil, nil); err != nil {
		return err
	}
	defer os.Remove(nameSorted)

	fixmated := filepath.Join(tmp, "rb-md-fixmate.bam")
	if _, err := runOnce(samtoolsBin,
		[]string{"fixmate", "-m", nameSorted, fixmated}, "", "", nil, nil); err != nil {
		return err
	}
	defer os.Remove(fixmated)

	_, err := runOnce(samtoolsBin,
		[]string{"sort", "-T", filepath.Join(tmp, "rb-md-cs"), "-o", dst, fixmated},
		"", "", nil, nil)
	return err
}

// writeBedGraph writes a 4-column BedGraph (chrom start end value) from BED3
// records: `bedtools unionbedg` needs a value column and SIGABRTs on a bare
// BED3. The value is the interval length (end-start), a deterministic integer
// derived from the coordinates; a non-numeric coordinate falls back to the
// 1-based record index so the column is always a valid integer.
func writeBedGraph(path string, recs []bedRec) error {
	return writeLines(path, func(w *bufio.Writer) error {
		for i, r := range recs {
			val := i + 1
			if s, err := strconv.Atoi(r.start); err == nil {
				if e, err := strconv.Atoi(r.end); err == nil && e >= s {
					val = e - s
				}
			}
			if _, err := fmt.Fprintf(w, "%s\t%s\t%s\t%d\n", r.chrom, r.start, r.end, val); err != nil {
				return err
			}
		}
		return nil
	})
}

// writeBED4 writes chrom/start/end plus a deterministic name column.
func writeBED4(path string, recs []bedRec) error {
	return writeLines(path, func(w *bufio.Writer) error {
		for i, r := range recs {
			if _, err := fmt.Fprintf(w, "%s\t%s\t%s\tregion_%d\n", r.chrom, r.start, r.end, i+1); err != nil {
				return err
			}
		}
		return nil
	})
}

// writeWindow pairs each record with the next one into an 8-field line
// (chrom start end name)x2 — the shape `bedtools window` produces and the shape
// `bedoverlap -cols 2,3,6,7` consumes. The last odd record is dropped so every
// emitted line is well-formed.
func writeWindow(path string, recs []bedRec) error {
	return writeLines(path, func(w *bufio.Writer) error {
		for i := 0; i+1 < len(recs); i += 2 {
			a, b := recs[i], recs[i+1]
			if _, err := fmt.Fprintf(w, "%s\t%s\t%s\tregion_%d\t%s\t%s\t%s\tregion_%d\n",
				a.chrom, a.start, a.end, i+1, b.chrom, b.start, b.end, i+2); err != nil {
				return err
			}
		}
		return nil
	})
}

// writeBEDPE pairs each record with the next into a 10-field BEDPE line
// (chrom1 s1 e1 chrom2 s2 e2 name score strand1 strand2). The last odd record is
// dropped so every line has the full pair.
func writeBEDPE(path string, recs []bedRec) error {
	return writeLines(path, func(w *bufio.Writer) error {
		for i := 0; i+1 < len(recs); i += 2 {
			a, b := recs[i], recs[i+1]
			if _, err := fmt.Fprintf(w, "%s\t%s\t%s\t%s\t%s\t%s\tpair_%d\t0\t+\t+\n",
				a.chrom, a.start, a.end, b.chrom, b.start, b.end, i/2+1); err != nil {
				return err
			}
		}
		return nil
	})
}

// writeLines creates path, runs fn against a buffered writer, and flushes.
func writeLines(path string, fn func(*bufio.Writer) error) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	if err := fn(w); err != nil {
		return err
	}
	return w.Flush()
}
