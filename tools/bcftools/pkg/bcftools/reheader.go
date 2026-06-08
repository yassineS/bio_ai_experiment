// Package bcftools — `bcftools reheader` subcommand.
//
// `bcftools reheader` rewrites the header of a VCF/BCF in place. The body
// records are passed through unchanged (subject to the new header's contig
// order for BCF outputs). Three input modes are supported:
//
//   - `-h FILE` — replace the entire header with the contents of FILE.
//   - `-s FILE` — rename samples; FILE is either a one-name-per-line list
//     (positional rename) or a tab-separated `OLD\tNEW` map (by name).
//   - `-f FAI`  — rebuild `##contig` lines from a samtools FAI sidecar.
//
// Modes can be combined: the order of operations is (header replace) →
// (FAI contigs) → (sample rename), matching upstream behaviour.
//
// For BCF input, upstream htslib requires a careful header rewrite so the
// dictionary indices are preserved; this v1 port reads the BCF, applies the
// edits to the in-memory header, and re-emits the records via the BCF
// writer. The records are conceptually unchanged but are re-encoded.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// ReheaderOptions controls the `bcftools reheader` subcommand.
type ReheaderOptions struct {
	// HeaderFile, when non-empty, is read verbatim and used as the new
	// header. Lines starting with `#` are kept; the file is expected to
	// terminate with the `#CHROM` header line.
	HeaderFile string
	// SamplesFile, when non-empty, drives the sample-rename step. The
	// format is auto-detected: if any non-comment line contains a tab, it
	// is treated as a `OLD\tNEW` mapping; otherwise the file is a flat
	// list of new names applied positionally.
	SamplesFile string
	// SamplesList is upstream `-n/--samples-list LIST` — a single
	// comma-separated list of new sample names. When set the names are
	// applied positionally (same semantics as SamplesFile's flat-list
	// branch). Both SamplesList and SamplesFile being set is an error.
	SamplesList []string
	// FaiFile, when non-empty, is a samtools FAI index whose first two
	// columns (NAME, LENGTH) are used to rebuild the `##contig` lines.
	FaiFile string
	// OutputFormat selects the output encoding. Defaults to OutputVCF.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output (negative means
	// gzip's default).
	CompressLevel int
}

// ReheaderFile is the file-aware entry point for `bcftools reheader`. It
// opens path through iohelper, applies the requested header edits, and
// writes the records (with the new header) to out.
func ReheaderFile(path string, out io.Writer, opts ReheaderOptions) (int, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools reheader: open %s: %w", path, err)
	}
	defer in.Close()
	return Reheader(in, out, opts)
}

// Reheader reads every record from in, edits the header per opts, and writes
// the records (with the new header) to out.
func Reheader(in io.Reader, out io.Writer, opts ReheaderOptions) (int, error) {
	hdr, recs, err := readAllVariants(in)
	if err != nil {
		return 0, fmt.Errorf("bcftools reheader: %w", err)
	}

	if opts.HeaderFile != "" {
		newHdr, err := loadHeaderFromFile(opts.HeaderFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools reheader: load -h %s: %w", opts.HeaderFile, err)
		}
		// Preserve the original sample list when the new header has none —
		// matches upstream's behaviour of complaining about #sample mismatch
		// otherwise. Here we keep the old samples if the new header omitted
		// the `#CHROM` line.
		if len(newHdr.Samples) == 0 {
			newHdr.Samples = hdr.Samples
		}
		hdr = newHdr
	}

	if opts.FaiFile != "" {
		if hdr, err = applyFaiContigs(hdr, opts.FaiFile); err != nil {
			return 0, fmt.Errorf("bcftools reheader: -f %s: %w", opts.FaiFile, err)
		}
	}

	if opts.SamplesFile != "" && len(opts.SamplesList) > 0 {
		return 0, fmt.Errorf("bcftools reheader: --samples-file and --samples-list are mutually exclusive")
	}
	if opts.SamplesFile != "" {
		mapping, names, err := loadSamplesRename(opts.SamplesFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools reheader: -N %s: %w", opts.SamplesFile, err)
		}
		hdr = renameHeaderSamples(hdr, mapping, names)
		for _, v := range recs {
			renameVariantSamples(v, mapping, names, hdr.Samples)
		}
	} else if len(opts.SamplesList) > 0 {
		hdr = renameHeaderSamples(hdr, nil, opts.SamplesList)
		for _, v := range recs {
			renameVariantSamples(v, nil, opts.SamplesList, hdr.Samples)
		}
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:      opts.OutputFormat,
		CompressLevel:     opts.CompressLevel,
		SkipPASSInjection: true,
	}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}
	for _, v := range recs {
		if err := w.Write(v); err != nil {
			return 0, err
		}
	}
	return len(recs), w.Flush()
}

// loadHeaderFromFile reads a plain VCF header (lines beginning with `#`)
// from path. Trailing data lines are ignored — they're conceptually
// disallowed in `-h` files anyway.
func loadHeaderFromFile(path string) (*vcf.Header, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	out := &vcf.Header{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 10<<20)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "##") {
			out.MetaInfo = append(out.MetaInfo, line)
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			fields := strings.Split(line, "\t")
			if len(fields) > 9 {
				out.Samples = append([]string{}, fields[9:]...)
			}
			break
		}
	}
	return out, sc.Err()
}

// applyFaiContigs rewrites every `##contig=<ID=...>` line in hdr to match the
// (NAME, LENGTH) entries of a samtools FAI sidecar. New contigs that are
// present in the FAI but missing from the header are appended in FAI order;
// contigs in the header that are absent from the FAI are dropped (this
// matches upstream's `-f` behaviour).
func applyFaiContigs(hdr *vcf.Header, faiPath string) (*vcf.Header, error) {
	f, err := os.Open(faiPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	type entry struct {
		name string
		len  int
	}
	var entries []entry
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 2 {
			return nil, fmt.Errorf("invalid FAI line %q", line)
		}
		ln, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid FAI length %q: %w", fields[1], err)
		}
		entries = append(entries, entry{name: fields[0], len: ln})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	out := &vcf.Header{Samples: hdr.Samples}
	// Insertion strategy: keep every non-contig meta line in its original
	// position; replace the *first* contig block with the FAI-derived
	// lines.
	contigInserted := false
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##contig=") {
			if !contigInserted {
				for _, e := range entries {
					out.MetaInfo = append(out.MetaInfo,
						fmt.Sprintf("##contig=<ID=%s,length=%d>", e.name, e.len))
				}
				contigInserted = true
			}
			continue
		}
		out.MetaInfo = append(out.MetaInfo, m)
	}
	if !contigInserted {
		// Header had no contigs at all — append the FAI block at the end.
		for _, e := range entries {
			out.MetaInfo = append(out.MetaInfo,
				fmt.Sprintf("##contig=<ID=%s,length=%d>", e.name, e.len))
		}
	}
	return out, nil
}

// loadSamplesRename parses the `-s/--samples` rename file. It returns a
// `OLD->NEW` map (non-nil only when the file is in tab-separated form) and a
// flat list of new names (non-nil only when the file is a plain list).
// Comment lines starting with `#` and blank lines are skipped.
func loadSamplesRename(path string) (map[string]string, []string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()
	var names []string
	mapping := map[string]string{}
	tabsSeen := false
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		if strings.IndexByte(line, '\t') >= 0 {
			tabsSeen = true
			fields := strings.SplitN(line, "\t", 2)
			mapping[fields[0]] = strings.TrimSpace(fields[1])
		} else {
			names = append(names, strings.TrimSpace(line))
		}
	}
	if err := sc.Err(); err != nil {
		return nil, nil, err
	}
	if tabsSeen {
		return mapping, nil, nil
	}
	return nil, names, nil
}

// renameHeaderSamples returns a copy of hdr with samples renamed per the
// chosen rule. Exactly one of `mapping` or `names` is set; the other is nil.
func renameHeaderSamples(hdr *vcf.Header, mapping map[string]string, names []string) *vcf.Header {
	out := &vcf.Header{MetaInfo: append([]string{}, hdr.MetaInfo...)}
	if mapping != nil {
		for _, s := range hdr.Samples {
			if n, ok := mapping[s]; ok && n != "" {
				out.Samples = append(out.Samples, n)
			} else {
				out.Samples = append(out.Samples, s)
			}
		}
		return out
	}
	if names != nil {
		// Positional rename — clip to whichever is shorter so we never
		// drop or invent samples.
		out.Samples = make([]string, len(hdr.Samples))
		for i, s := range hdr.Samples {
			if i < len(names) && names[i] != "" {
				out.Samples[i] = names[i]
			} else {
				out.Samples[i] = s
			}
		}
		return out
	}
	out.Samples = append([]string{}, hdr.Samples...)
	return out
}

// renameVariantSamples re-keys each sample on v.Samples so the downstream
// writer sees the new names. We have to keep positions consistent with the
// rewritten header.
func renameVariantSamples(v *vcf.Variant, mapping map[string]string, names []string, newNames []string) {
	if len(v.Samples) == 0 {
		return
	}
	for i, s := range v.Samples {
		if mapping != nil {
			if n, ok := mapping[s.Name]; ok && n != "" {
				v.Samples[i].Name = n
			}
			continue
		}
		if names != nil {
			if i < len(names) && names[i] != "" {
				v.Samples[i].Name = names[i]
			}
			continue
		}
	}
	// As a safety net, ensure positions match the new header sample order.
	if len(newNames) == len(v.Samples) {
		for i := range v.Samples {
			v.Samples[i].Name = newNames[i]
		}
	}
}
