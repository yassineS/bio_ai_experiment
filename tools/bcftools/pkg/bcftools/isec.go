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
	// Threads is upstream's -@/--threads. When greater than 1 it enables
	// parallel BGZF compression of -O z and -O b output via bgzf.MultiWriter;
	// the framed result decodes byte-identically regardless of thread count.
	Threads int
	// InputPaths are the input file names, used only for the README.txt legend
	// (which names each private/shared output file after its source). Set by
	// IsecFiles; may be empty when Isec is called directly.
	InputPaths []string
	// CmdArgv is the `isec ...` argument vector echoed into README.txt's "The
	// command line was:" line. Set by the CLI; may be empty.
	CmdArgv []string
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
	if opts.InputPaths == nil {
		opts.InputPaths = paths
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

	// Sort the union key list by (contig order, POS) only, stably. Upstream's
	// synced reader presents records at the same position in file order, so the
	// tie-break must preserve first-appearance order (keyOrder is built in file
	// order) rather than re-sort by REF/ALT — otherwise two records at one POS
	// (e.g. A>T then A>C) would be emitted A>C-first.
	primaryOrder := contigOrder(headers[0])
	sort.SliceStable(keyOrder, func(i, j int) bool {
		ai, bi := keyVariants[keyOrder[i]][0].variant, keyVariants[keyOrder[j]][0].variant
		return keyFor(ai, primaryOrder).less(keyFor(bi, primaryOrder))
	})

	// The -n constraint is applied per matched site in the emit loop below (a
	// site is one paired occurrence across inputs), mirroring upstream's
	// per-record synced-reader test rather than a per-position one.

	// VENN mode: with exactly two inputs and no membership constraint, upstream
	// (vcfisec.c: nfiles==2 && !isec_op -> OP_VENN) writes the four-way Venn
	// decomposition — 0000 records private to input 1, 0001 private to input 2,
	// 0002 the input-1 copy of shared records, 0003 the input-2 copy — rather
	// than one file per input.
	venn := n == 2 && opts.Nfiles.Mode == 0

	// Open per-input output files when -p is set. In VENN mode there are four
	// files; file i draws its header from input vennReader[i].
	var perInputW []variantWriter
	var perInputClose []func()
	var sitesF *os.File
	vennReader := []int{0, 1, 0, 1}
	if opts.Prefix != "" {
		if err := os.MkdirAll(opts.Prefix, 0755); err != nil {
			return 0, fmt.Errorf("bcftools isec: mkdir %s: %w", opts.Prefix, err)
		}
		ext := outputExt(opts.OutputFormat)
		nOut := n
		hdrFor := func(i int) *vcf.Header { return headers[i] }
		if venn {
			nOut = 4
			hdrFor = func(i int) *vcf.Header { return headers[vennReader[i]] }
		}
		// README.txt names each output file after its source; build it with the
		// real output paths so the legend matches upstream.
		outNames := make([]string, nOut)
		for i := 0; i < nOut; i++ {
			outNames[i] = filepath.Join(opts.Prefix, fmt.Sprintf("%04d.vcf%s", i, ext))
		}
		readmePath := filepath.Join(opts.Prefix, "README.txt")
		_ = os.WriteFile(readmePath, []byte(isecReadme(opts, outNames, venn)), 0644)
		for i := 0; i < nOut; i++ {
			f, err := os.Create(outNames[i])
			if err != nil {
				return 0, fmt.Errorf("bcftools isec: create %s: %w", outNames[i], err)
			}
			w, finish, err := openOutput(f, ViewOptions{
				OutputFormat:  opts.OutputFormat,
				CompressLevel: opts.CompressLevel,
				Threads:       opts.Threads,
			}, hdrFor(i))
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
		// sites.txt (one line per matched site, written in the per-occurrence
		// emit loop below) lists every site's CHROM POS REF ALT and membership
		// bits, matching upstream's `CHROM\tPOS\tREF\tALT\tBITS` shape.
		sitesPath := filepath.Join(opts.Prefix, "sites.txt")
		sf, err := os.Create(sitesPath)
		if err != nil {
			return 0, fmt.Errorf("bcftools isec: create %s: %w", sitesPath, err)
		}
		defer sf.Close()
		sitesF = sf
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
			Threads:       opts.Threads,
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

	// Walk the keys and emit one logical site per matched occurrence. When a
	// key holds several records from one input (intra-position duplicates),
	// upstream's synced reader pairs the k-th occurrence across inputs, so we
	// iterate occurrences rather than collapse them — each occurrence has its
	// own membership and its own sites.txt line.
	totalKept := 0
	bits := make([]byte, n)
	for _, k := range keyOrder {
		// Group this key's records by input, preserving file order.
		byFile := make([][]*vcf.Variant, n)
		for _, c := range keyVariants[k] {
			byFile[c.groupIdx] = append(byFile[c.groupIdx], c.variant)
		}
		maxOcc := 0
		for gi := range byFile {
			if len(byFile[gi]) > maxOcc {
				maxOcc = len(byFile[gi])
			}
		}
		for occ := 0; occ < maxOcc; occ++ {
			present := make([]bool, n)
			firstGi := -1
			for gi := range byFile {
				if occ < len(byFile[gi]) {
					present[gi] = true
					if firstGi < 0 {
						firstGi = gi
					}
				}
			}
			if !nfilesPasses(present, opts.Nfiles) {
				continue
			}
			totalKept++
			rep := byFile[firstGi][occ]
			// sites.txt: CHROM POS REF ALT BITS (first present input's record).
			if sitesF != nil {
				for gi := 0; gi < n; gi++ {
					if present[gi] {
						bits[gi] = '1'
					} else {
						bits[gi] = '0'
					}
				}
				ref := rep.Ref
				if ref == "" {
					ref = "."
				}
				alt := "."
				if len(rep.Alt) > 0 {
					alt = strings.Join(rep.Alt, ",")
				}
				fmt.Fprintf(sitesF, "%s\t%d\t%s\t%s\t%s\n", rep.Chrom, rep.Pos, ref, alt, string(bits))
			}
			// Per-input projection.
			if opts.Prefix != "" {
				if venn && present[0] && present[1] {
					// Shared: input-1 copy -> 0002, input-2 copy -> 0003.
					if err := perInputW[2].Write(byFile[0][occ]); err != nil {
						return totalKept, err
					}
					if err := perInputW[3].Write(byFile[1][occ]); err != nil {
						return totalKept, err
					}
				} else {
					for gi := range byFile {
						if present[gi] {
							if err := perInputW[gi].Write(byFile[gi][occ]); err != nil {
								return totalKept, err
							}
						}
					}
				}
			}
			// stdout dump: first requested input that is present at this site.
			if openStdout {
				for _, want := range writeFromInputs {
					wi := want - 1
					if wi < 0 || wi >= n || !present[wi] {
						continue
					}
					if err := stdoutW.Write(byFile[wi][occ]); err != nil {
						return totalKept, err
					}
					break
				}
			}
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

// isecReadme renders the README.txt legend in upstream vcfisec's exact shape:
// a provenance preamble, then one "<outfile>\tfor records private to\t<input>"
// (and, in VENN mode, "<outfile>\tfor records from <input> shared by both\t..."
// ) line per output file. outNames are the full output paths; venn selects the
// four-file Venn legend over the one-file-per-input legend.
func isecReadme(opts IsecOptions, outNames []string, venn bool) string {
	var b strings.Builder
	b.WriteString("This file was produced by vcfisec.\n")
	b.WriteString("The command line was:\tbcftools")
	if len(opts.CmdArgv) > 0 {
		b.WriteString(" ")
		b.WriteString(opts.CmdArgv[0])
		b.WriteString(" ")
		for _, a := range opts.CmdArgv[1:] {
			b.WriteString(" ")
			b.WriteString(a)
		}
	}
	b.WriteString("\n\nUsing the following file names:\n")
	in := opts.InputPaths
	get := func(i int) string {
		if i < len(in) {
			return in[i]
		}
		return fmt.Sprintf("input%d", i+1)
	}
	if venn && len(outNames) == 4 {
		fmt.Fprintf(&b, "%s\tfor records private to\t%s\n", outNames[0], get(0))
		fmt.Fprintf(&b, "%s\tfor records private to\t%s\n", outNames[1], get(1))
		fmt.Fprintf(&b, "%s\tfor records from %s shared by both\t%s %s\n", outNames[2], get(0), get(0), get(1))
		fmt.Fprintf(&b, "%s\tfor records from %s shared by both\t%s %s\n", outNames[3], get(1), get(0), get(1))
	} else {
		for i, name := range outNames {
			fmt.Fprintf(&b, "%s\tfor stripped\t%s\n", name, get(i))
		}
	}
	return b.String()
}
