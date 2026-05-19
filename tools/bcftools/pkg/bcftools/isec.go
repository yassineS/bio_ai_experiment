// Package bcftools — `bcftools isec` subcommand.
//
// `bcftools isec` performs set operations on N sorted VCF/BCF inputs. The
// canonical use is "find variants present in at least N of the inputs", but
// the flags also support exact membership counts and explicit inclusion
// bitmasks. Outputs can be either:
//
//   - the per-input projections written to `<prefix>/000<i>.vcf` (one file
//     per input, holding the records selected from that input — the
//     `-p`/`--prefix` mode), or
//   - one or more inputs written to stdout (the `-w`/`--write` mode).
//
// The two modes can be combined: `-p` always writes the per-input projections,
// and `-w` additionally streams the chosen inputs to stdout.
//
// Records collapse on (CHROM, POS, REF, ALT) by default. The `-c/--collapse`
// flag widens the comparison: `none` keeps the strict tuple, `snps` collapses
// any pair of SNP-records at the same site regardless of allele, `indels`
// likewise for indels, `both` is the union of `snps`+`indels`, `all`
// collapses everything at the same (CHROM, POS), `id` matches on the ID
// column, and `some` is upstream's hybrid (REF must match, any ALT in common
// is enough). The default is `none`.
package bcftools

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// CollapseMode controls how `isec` decides two records "match".
type CollapseMode int

const (
	// CollapseNone requires identical (CHROM, POS, REF, ALT). Default.
	CollapseNone CollapseMode = iota
	// CollapseSNPs collapses any two SNP records at the same site.
	CollapseSNPs
	// CollapseIndels collapses any two indel records at the same site.
	CollapseIndels
	// CollapseBoth is the union of CollapseSNPs and CollapseIndels.
	CollapseBoth
	// CollapseAll collapses every record at the same (CHROM, POS).
	CollapseAll
	// CollapseSome requires REF match and at least one ALT in common.
	CollapseSome
	// CollapseID matches on the ID column.
	CollapseID
)

// ParseCollapseMode parses the -c/--collapse flag.
func ParseCollapseMode(s string) (CollapseMode, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "none":
		return CollapseNone, nil
	case "snps":
		return CollapseSNPs, nil
	case "indels":
		return CollapseIndels, nil
	case "both":
		return CollapseBoth, nil
	case "all":
		return CollapseAll, nil
	case "some":
		return CollapseSome, nil
	case "id":
		return CollapseID, nil
	}
	return CollapseNone, fmt.Errorf("bcftools isec: unknown -c mode %q (expect none|snps|indels|both|all|some|id)", s)
}

// NfilesSpec encodes the `-n/--nfiles` argument.
type NfilesSpec struct {
	// Mode: '~' bitmask, '+' at-least-N, '=' exactly-N, 0 unset.
	Mode byte
	// N is the threshold for `+N` and `=N`. For `~bits` it is unused.
	N int
	// Bits is the inclusion mask for `~bits` mode (one bit per input,
	// bit 0 = first input).
	Bits []bool
}

// ParseNfilesSpec parses the -n/--nfiles argument. Empty s returns the
// "no constraint" spec (Mode == 0).
func ParseNfilesSpec(s string) (NfilesSpec, error) {
	out := NfilesSpec{}
	s = strings.TrimSpace(s)
	if s == "" {
		return out, nil
	}
	if len(s) >= 2 && (s[0] == '+' || s[0] == '=') {
		var n int
		_, err := fmt.Sscanf(s[1:], "%d", &n)
		if err != nil {
			return out, fmt.Errorf("bcftools isec: bad -n %q", s)
		}
		out.Mode = s[0]
		out.N = n
		return out, nil
	}
	if s[0] == '~' {
		bits := make([]bool, 0, len(s)-1)
		for _, r := range s[1:] {
			switch r {
			case '0':
				bits = append(bits, false)
			case '1':
				bits = append(bits, true)
			default:
				return out, fmt.Errorf("bcftools isec: bad -n bitmask %q (expect 0/1 chars)", s)
			}
		}
		out.Mode = '~'
		out.Bits = bits
		return out, nil
	}
	// Bare integer == "+N".
	var n int
	_, err := fmt.Sscanf(s, "%d", &n)
	if err != nil {
		return out, fmt.Errorf("bcftools isec: bad -n %q", s)
	}
	out.Mode = '+'
	out.N = n
	return out, nil
}

// IsecOptions controls the `bcftools isec` subcommand.
type IsecOptions struct {
	// Nfiles is the membership constraint. When unset (Mode == 0) all
	// records that appear in any input are emitted.
	Nfiles NfilesSpec
	// Collapse controls how records collapse. Defaults to CollapseNone.
	Collapse CollapseMode
	// Prefix is the per-input output directory: when non-empty, one file
	// `<prefix>/000<i>.vcf[.gz|.bcf]` is written per input.
	Prefix string
	// Write is the 1-based input indices to dump to stdout (the `-w` flag).
	// When empty and Prefix is also empty, every selected record is dumped
	// from input 1.
	Write []int
	// OutputFormat selects the encoding for both the prefix files and the
	// stdout stream.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
}

// IsecFiles is the file-aware entry point. It opens each path through
// iohelper, fully reads them into memory, and runs Isec.
func IsecFiles(paths []string, out io.Writer, opts IsecOptions) (int, error) {
	if len(paths) < 2 {
		return 0, fmt.Errorf("bcftools isec: need at least two input files (got %d)", len(paths))
	}
	headers := make([]*vcf.Header, len(paths))
	groups := make([][]*vcf.Variant, len(paths))
	for i, p := range paths {
		in, err := iohelper.OpenReader(p)
		if err != nil {
			return 0, fmt.Errorf("bcftools isec: open %s: %w", p, err)
		}
		hdr, recs, err := readAllVariants(in)
		_ = in.Close()
		if err != nil {
			return 0, fmt.Errorf("bcftools isec: %s: %w", p, err)
		}
		headers[i] = hdr
		groups[i] = recs
	}
	return Isec(headers, groups, out, opts)
}

// Isec performs the set operation on pre-loaded headers/groups, writing the
// per-input projections to `<prefix>/000<i>.vcf*` (when Prefix is set) and
// the stdout stream defined by `-w` (when Write is set, defaulting to input
// 1 when neither is set).
func Isec(headers []*vcf.Header, groups [][]*vcf.Variant, stdout io.Writer, opts IsecOptions) (int, error) {
	if len(headers) != len(groups) {
		return 0, fmt.Errorf("bcftools isec: header / variant count mismatch")
	}
	n := len(headers)
	if n == 0 {
		return 0, fmt.Errorf("bcftools isec: no inputs")
	}
	if opts.Nfiles.Mode == '~' && len(opts.Nfiles.Bits) != n {
		return 0, fmt.Errorf("bcftools isec: -n ~bits length %d != number of inputs %d", len(opts.Nfiles.Bits), n)
	}

	// Build per-input keyed maps and the union key list with membership.
	keyOf := func(v *vcf.Variant) string {
		return collapseKey(v, opts.Collapse)
	}
	type cell struct {
		variant   *vcf.Variant
		groupIdx  int
		variantIx int // index within its group
	}
	keyMembership := map[string]*[]bool{}
	keyOrder := []string{}
	keyVariants := map[string][]cell{}
	for gi, g := range groups {
		for vi, v := range g {
			k := keyOf(v)
			if _, ok := keyMembership[k]; !ok {
				m := make([]bool, n)
				keyMembership[k] = &m
				keyOrder = append(keyOrder, k)
			}
			(*keyMembership[k])[gi] = true
			keyVariants[k] = append(keyVariants[k], cell{variant: v, groupIdx: gi, variantIx: vi})
		}
	}

	// Sort the union key list using the first-input contig order as the
	// primary signal. Ties fall back to lexicographic.
	primaryOrder := contigOrder(headers[0])
	sort.SliceStable(keyOrder, func(i, j int) bool {
		ai, bi := keyVariants[keyOrder[i]][0].variant, keyVariants[keyOrder[j]][0].variant
		ka := keyFor(ai, primaryOrder)
		kb := keyFor(bi, primaryOrder)
		if !ka.equal(kb) {
			return ka.less(kb)
		}
		return keyOrder[i] < keyOrder[j]
	})

	// Decide which keys pass the -n constraint.
	keep := map[string]bool{}
	for _, k := range keyOrder {
		mem := *keyMembership[k]
		if !nfilesPasses(mem, opts.Nfiles) {
			continue
		}
		keep[k] = true
	}

	// Open per-input output files when -p is set.
	var perInputW []variantWriter
	var perInputClose []func()
	if opts.Prefix != "" {
		if err := os.MkdirAll(opts.Prefix, 0755); err != nil {
			return 0, fmt.Errorf("bcftools isec: mkdir %s: %w", opts.Prefix, err)
		}
		// Drop a small README in the prefix (matches upstream convention).
		readmePath := filepath.Join(opts.Prefix, "README.txt")
		_ = os.WriteFile(readmePath, []byte(isecReadme(headers, opts)), 0644)
		// One file per input, named 0000.vcf, 0001.vcf, ... matching upstream.
		ext := outputExt(opts.OutputFormat)
		for i := range groups {
			path := filepath.Join(opts.Prefix, fmt.Sprintf("%04d.vcf%s", i, ext))
			f, err := os.Create(path)
			if err != nil {
				return 0, fmt.Errorf("bcftools isec: create %s: %w", path, err)
			}
			w, finish, err := openOutput(f, ViewOptions{
				OutputFormat:  opts.OutputFormat,
				CompressLevel: opts.CompressLevel,
			}, headers[i])
			if err != nil {
				_ = f.Close()
				return 0, err
			}
			closeF := f
			perInputW = append(perInputW, w)
			perInputClose = append(perInputClose, func() {
				finish()
				_ = closeF.Close()
			})
			if err := w.WriteHeader(); err != nil {
				return 0, err
			}
		}
		// Also write a sites.txt listing kept positions and membership bits
		// — the upstream shape is `CHROM POS REF ALT N\tBITS`.
		sitesPath := filepath.Join(opts.Prefix, "sites.txt")
		sitesF, err := os.Create(sitesPath)
		if err != nil {
			return 0, fmt.Errorf("bcftools isec: create %s: %w", sitesPath, err)
		}
		defer sitesF.Close()
		for _, k := range keyOrder {
			if !keep[k] {
				continue
			}
			cells := keyVariants[k]
			v := cells[0].variant
			mem := *keyMembership[k]
			bits := make([]byte, n)
			for i, b := range mem {
				if b {
					bits[i] = '1'
				} else {
					bits[i] = '0'
				}
			}
			fmt.Fprintf(sitesF, "%s\t%d\t%s\t%s\t%s\n", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","), string(bits))
		}
	}

	// Open the stdout writer if -w is set or both -p and -w are empty
	// (default = dump from the first input).
	writeFromInputs := opts.Write
	openStdout := opts.Prefix == "" || len(writeFromInputs) > 0
	var stdoutW variantWriter
	var stdoutFinish func()
	if openStdout {
		// Header source: union (first input) when multiple write-targets,
		// or the chosen input when exactly one.
		hdrIdx := 0
		if len(writeFromInputs) == 1 {
			if writeFromInputs[0]-1 >= 0 && writeFromInputs[0]-1 < n {
				hdrIdx = writeFromInputs[0] - 1
			}
		}
		w, finish, err := openOutput(stdout, ViewOptions{
			OutputFormat:  opts.OutputFormat,
			CompressLevel: opts.CompressLevel,
		}, headers[hdrIdx])
		if err != nil {
			return 0, err
		}
		stdoutW = w
		stdoutFinish = finish
		if err := stdoutW.WriteHeader(); err != nil {
			return 0, err
		}
		if len(writeFromInputs) == 0 {
			writeFromInputs = []int{1}
		}
	}

	// Walk the keys and emit.
	totalKept := 0
	for _, k := range keyOrder {
		if !keep[k] {
			continue
		}
		totalKept++
		cells := keyVariants[k]
		// Per-input projection: write each cell to its respective per-input
		// file.
		if opts.Prefix != "" {
			for _, c := range cells {
				if err := perInputW[c.groupIdx].Write(c.variant); err != nil {
					return totalKept, err
				}
			}
		}
		// stdout dump: pick the first cell whose group is in writeFromInputs.
		if openStdout {
			for _, want := range writeFromInputs {
				wi := want - 1
				if wi < 0 || wi >= n {
					continue
				}
				for _, c := range cells {
					if c.groupIdx != wi {
						continue
					}
					if err := stdoutW.Write(c.variant); err != nil {
						return totalKept, err
					}
					goto wrote
				}
			}
		wrote:
		}
	}
	if openStdout {
		_ = stdoutW.Flush()
		stdoutFinish()
	}
	for _, w := range perInputW {
		_ = w.Flush()
	}
	for _, c := range perInputClose {
		c()
	}
	return totalKept, nil
}

// collapseKey returns the canonical membership key for a variant under the
// given collapse mode.
func collapseKey(v *vcf.Variant, mode CollapseMode) string {
	switch mode {
	case CollapseNone:
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseAll:
		return fmt.Sprintf("%s\t%d", v.Chrom, v.Pos)
	case CollapseID:
		return v.ID
	case CollapseSNPs:
		if isSNPRecord(v) {
			return fmt.Sprintf("%s\t%d\tSNP", v.Chrom, v.Pos)
		}
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseIndels:
		if isIndelRecord(v) {
			return fmt.Sprintf("%s\t%d\tINDEL", v.Chrom, v.Pos)
		}
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseBoth:
		switch {
		case isSNPRecord(v):
			return fmt.Sprintf("%s\t%d\tSNP", v.Chrom, v.Pos)
		case isIndelRecord(v):
			return fmt.Sprintf("%s\t%d\tINDEL", v.Chrom, v.Pos)
		}
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseSome:
		return fmt.Sprintf("%s\t%d\t%s", v.Chrom, v.Pos, v.Ref)
	}
	return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
}

// nfilesPasses reports whether mem (the per-input membership bitmap) is
// admissible under spec.
func nfilesPasses(mem []bool, spec NfilesSpec) bool {
	if spec.Mode == 0 {
		// No constraint — keep every key that appears in at least one
		// input (which is always the case).
		return true
	}
	count := 0
	for _, b := range mem {
		if b {
			count++
		}
	}
	switch spec.Mode {
	case '+':
		return count >= spec.N
	case '=':
		return count == spec.N
	case '~':
		if len(spec.Bits) != len(mem) {
			return false
		}
		for i := range mem {
			if mem[i] != spec.Bits[i] {
				return false
			}
		}
		return true
	}
	return true
}

// outputExt returns the filename extension matching the chosen output
// format.
func outputExt(f OutputFormat) string {
	switch f {
	case OutputVCFGz:
		return ".gz"
	case OutputBCF:
		return "" // .vcf.bcf doesn't make sense — upstream uses .bcf, see below
	case OutputBCFUncompressed:
		return ""
	}
	return ""
}

// isecReadme is a small explainer dropped into the prefix dir so users know
// what the per-input files contain. Matches upstream's habit of writing a
// README.txt next to the outputs.
func isecReadme(headers []*vcf.Header, opts IsecOptions) string {
	var b strings.Builder
	b.WriteString("This directory was produced by `bcftools isec` (Go port).\n\n")
	b.WriteString("Per-input projections:\n")
	for i := range headers {
		b.WriteString(fmt.Sprintf("  000%d.vcf - records from input %d that pass the membership constraint\n", i, i+1))
	}
	b.WriteString("\nMembership constraint:\n")
	switch opts.Nfiles.Mode {
	case 0:
		b.WriteString("  (none — all keys present in any input)\n")
	case '+':
		b.WriteString(fmt.Sprintf("  +%d (present in at least %d input files)\n", opts.Nfiles.N, opts.Nfiles.N))
	case '=':
		b.WriteString(fmt.Sprintf("  =%d (present in exactly %d input files)\n", opts.Nfiles.N, opts.Nfiles.N))
	case '~':
		bits := make([]byte, len(opts.Nfiles.Bits))
		for i, b2 := range opts.Nfiles.Bits {
			if b2 {
				bits[i] = '1'
			} else {
				bits[i] = '0'
			}
		}
		b.WriteString(fmt.Sprintf("  ~%s (bitmask)\n", string(bits)))
	}
	return b.String()
}
