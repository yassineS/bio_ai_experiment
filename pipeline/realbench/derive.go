package realbench

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// deriveInputs synthesises the extra bed* inputs that a plain BED3 (the NIST
// high-confidence intervals) cannot satisfy, and returns a copy of in with the
// derived paths filled. The synthetics are written into dir (the run's scratch
// dir) so they exist on the real run too, and they are derived DETERMINISTICALLY
// from the real BED so both ours and upstream see byte-identical inputs:
//
//   - BED4: the BED3 with a synthetic name column (region_N) appended. Feeds
//     the bed* subcommands that require a name field (tobam -> QNAME).
//   - Window: consecutive BED records paired into an 8-field
//     "chrom start end name" x2 line (the shape `bedtools window` emits),
//     for `bedoverlap -cols 2,3,6,7`.
//   - BEDPE: consecutive BED records paired into a 10-field BEDPE line
//     (chrom1 s1 e1 chrom2 s2 e2 name score strand1 strand2), for
//     pairtopair / pairtobed.
//
// When in.BED is empty nothing is derived (the dependent cells SKIP). A
// synthesis error is returned so the caller can surface it rather than silently
// running invalid cells.
func deriveInputs(in Inputs, dir string) (Inputs, error) {
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
