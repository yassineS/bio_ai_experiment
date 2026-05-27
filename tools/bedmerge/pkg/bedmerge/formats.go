// Format adapters for non-BED inputs. Upstream `bedtools merge` accepts BED,
// GFF and VCF inputs; this file provides thin streaming converters so the
// rest of the package can keep treating its input as BED text.
//
// The adapters emit `chrom\tstart\tend\n` lines only (the upstream `merge`
// output for VCF/GFF is BED3 — see test-merge.sh cases t13/t14), which is
// what the merge layer ultimately needs.
package bedmerge

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// NewGFFToBEDReader wraps an io.Reader containing GFF/GTF text and emits
// BED3 lines (chrom, start-1, end) suitable for feeding into Merge. Header
// lines ("##", "browser", "track") and blank lines are dropped. Lines with
// fewer than 5 tab-separated columns or with unparseable coordinates are
// skipped (matching upstream's permissive behaviour on auxiliary header
// text).
func NewGFFToBEDReader(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		sc := bufio.NewScanner(r)
		sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
		for sc.Scan() {
			line := sc.Text()
			trim := strings.TrimSpace(line)
			if trim == "" || strings.HasPrefix(trim, "#") ||
				strings.HasPrefix(trim, "browser") || strings.HasPrefix(trim, "track") {
				continue
			}
			fields := strings.Split(line, "\t")
			if len(fields) < 5 {
				continue
			}
			s, err := strconv.Atoi(fields[3])
			if err != nil || s < 1 {
				continue
			}
			e, err := strconv.Atoi(fields[4])
			if err != nil || e < s {
				continue
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\n", fields[0], s-1, e); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
		if err := sc.Err(); err != nil {
			_ = pw.CloseWithError(err)
		}
	}()
	return pr
}

// NewVCFToBEDReader wraps an io.Reader containing VCF text and emits BED3
// lines (chrom, POS-1, POS-1+len(REF)) suitable for feeding into Merge.
// VCF header lines (those starting with "#" — including "##" metadata and
// the "#CHROM" column header) are dropped. The interval span is derived
// from the REF allele length so deletions and complex variants span the
// reference bases they consume; matches upstream `bedtools merge -i x.vcf`.
func NewVCFToBEDReader(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		sc := bufio.NewScanner(r)
		sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
		for sc.Scan() {
			line := sc.Text()
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			fields := strings.Split(line, "\t")
			if len(fields) < 4 {
				continue
			}
			pos, err := strconv.Atoi(fields[1])
			if err != nil || pos < 1 {
				continue
			}
			ref := fields[3]
			refLen := len(ref)
			if refLen == 0 {
				refLen = 1
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\n", fields[0], pos-1, pos-1+refLen); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
		if err := sc.Err(); err != nil {
			_ = pw.CloseWithError(err)
		}
	}()
	return pr
}
