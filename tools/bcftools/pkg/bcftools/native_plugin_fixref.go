// Native port of the upstream `fixref` plugin (plugins/fixref.c). It determines
// and fixes REF-allele strand orientation against a FASTA reference. The
// supported modes mirror the C tool: stats (collect+print stats, no VCF
// output), ref-alt, swap, flip, flip-all, and top (Illumina TOP -> fwd with
// ambiguous-pair sequence walking). Each converted record is annotated with an
// INFO/FIXREF (or -t named) string recording the change (none/flip/swap/GT/...).
//
// The `id`/`--use-id` mode (MODE_USE_ID) is implemented in the companion file
// native_plugin_fixref_id.go: it determines the REF allele from a separate
// dbSNP VCF keyed by the ID column and swaps REF/ALT (and the genotypes)
// accordingly.
package bcftools

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("fixref", func() NativePlugin { return &fixrefPlugin{} })
}

// fixref modes, mirroring the MODE_* constants in fixref.c.
const (
	fixrefModeStats      = iota + 1 // collect and print stats, no VCF output
	fixrefModeTop2Fwd               // Illumina TOP strand -> fwd
	fixrefModeFlip2Fwd              // flip/swap non-ambiguous SNPs
	fixrefModeUseID                 // determine REF from a dbSNP VCF keyed by ID
	fixrefModeRefAlt                // swap/flip to match REF, leave GTs
	fixrefModeFlipAll               // flip/swap all SNPs including ambiguous
	fixrefModeSwapRefAlt            // only swap to match REF, leave GTs
)

// fixref FIX_* dirty bits, mirroring fixref.c. The annotation strings are
// emitted in bit order err,skip,none,flip,swap,GT.
const (
	fixFixErr  = 1 << 0
	fixFixSkip = 1 << 1
	fixFixNone = 1 << 2
	fixFixFlip = 1 << 3
	fixFixSwap = 1 << 4
	fixFixGT   = 1 << 5
)

var fixrefInfoAnnots = []string{"err", "skip", "none", "flip", "swap", "GT"}

// fixrefPlugin implements the `fixref` plugin. Conversions only read from the
// FASTA, so Process is per-record and parallel. The aggregate stats printed by
// Destroy are accumulated under a mutex-free per-record path because parallel
// counting would race; we therefore declare Parallel()==false to keep counts
// and the (stderr) stats deterministic and correct.
type fixrefPlugin struct {
	mode    int
	discard bool
	infoTag string
	fai     *fasta.RandomAccess
	stderr  io.Writer

	skipRID string // sequence name to skip after a faidx miss

	// --use-id (MODE_USE_ID) state. dbsnpFname is the dbSNP VCF; dbsnpMap is the
	// per-chromosome rsID -> {pos, ref} map rebuilt on a chromosome change
	// (dbsnpRID tracks the chromosome the map was built for). dbsnpPrevPos and
	// dbsnpUnsorted reproduce the one-shot "unsorted VCF" warning.
	dbsnpFname    string
	dbsnpMap      map[string]fixrefMarker
	dbsnpRID      string
	dbsnpPrevPos  int
	dbsnpUnsorted bool

	// statistics (mirrors the args_t counters)
	nsite, nok, nflip, nunresolved, nswap, nflipSwap    uint32
	nonSNP, nonACGT, nonbiallelic, nerr, dirty, npesErr uint32
	count                                               [4][4]uint32
}

// Name returns the plugin name.
func (p *fixrefPlugin) Name() string { return "fixref" }

// About returns the one-line description, matching fixref.c about().
func (p *fixrefPlugin) About() string {
	return "Fix reference strand orientation, e.g. from Illumina/TOP to fwd."
}

// Parallel reports false: Destroy prints aggregate stats accumulated across
// records, so Process must run serially in input order.
func (p *fixrefPlugin) Parallel() bool { return false }

// SetStderr wires the host stderr for the end-of-run stats block.
func (p *fixrefPlugin) SetStderr(w io.Writer) { p.stderr = w }

// SuppressVCF reports whether the VCF/BCF output is suppressed. In stats mode
// upstream init() returns 1 to emit only the stderr stats; mirror that here.
func (p *fixrefPlugin) SuppressVCF() bool { return p.mode == fixrefModeStats }

// SetStdout is part of outputSuppressor; fixref writes its report to stderr,
// not stdout, so the stdout writer is ignored.
func (p *fixrefPlugin) SetStdout(io.Writer) {}

// Init parses the plugin options, loads the FASTA, and appends the INFO/FIXREF
// header line.
func (p *fixrefPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.mode = fixrefModeStats
	p.infoTag = "FIXREF"
	var refFname string
	for i := 0; i < len(args); i++ {
		a := args[i]
		needVal := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("fixref: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-m", "--mode":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			switch strings.ToLower(v) {
			case "top":
				p.mode = fixrefModeTop2Fwd
			case "flip":
				p.mode = fixrefModeFlip2Fwd
			case "flip-all":
				p.mode = fixrefModeFlipAll
			case "id":
				p.mode = fixrefModeUseID
			case "ref-alt":
				p.mode = fixrefModeRefAlt
			case "swap":
				p.mode = fixrefModeSwapRefAlt
			case "stats":
				p.mode = fixrefModeStats
			default:
				return nil, fmt.Errorf("fixref: the source strand convention not recognised: %s", v)
			}
		case "-i", "--use-id":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			p.dbsnpFname = v
			p.mode = fixrefModeUseID
		case "-d", "--discard":
			p.discard = true
		case "-f", "--fasta-ref":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			refFname = v
		case "-t", "--tag-name":
			v, err := needVal()
			if err != nil {
				return nil, err
			}
			p.infoTag = v
		default:
			return nil, fmt.Errorf("fixref: unsupported option %q", a)
		}
	}

	if p.mode == fixrefModeUseID && p.dbsnpFname == "" {
		return nil, fmt.Errorf("fixref: No ID file specified, use -i/--use-id")
	}
	if refFname == "" {
		return nil, fmt.Errorf("fixref: expected the -f option")
	}
	fai, err := fasta.OpenRandomAccess(refFname)
	if err != nil {
		return nil, fmt.Errorf("fixref: failed to load the fai index: %w", err)
	}
	p.fai = fai

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	line := fmt.Sprintf(`##INFO=<ID=%s,Number=.,Type=String,Description="The change made by bcftools/fixref">`, p.infoTag)
	out.MetaInfo = appendInfoHeader(out.MetaInfo, line)
	return out, nil
}

// fixrefNt2int maps a base to 0/1/2/3 for A/C/G/T (case-insensitive), or -1.
func fixrefNt2int(nt byte) int {
	switch nt {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	}
	return -1
}

// fixrefInt2nt maps 0/1/2/3 back to A/C/G/T.
func fixrefInt2nt(x int) byte { return "ACGT"[x] }

// fixrefRevint complements an A/C/G/T index (A<->T, C<->G), mirroring revint().
func fixrefRevint(x int) int { return int("3210"[x] - '0') }

// fetchRefBase returns the 0/1/2/3 reference-base index at the record position,
// -1 for a non-ACGT base, or -2 when the contig is unknown (then the contig is
// flagged to be skipped, matching the C fetch_ref behaviour).
func (p *fixrefPlugin) fetchRefBase(v *vcf.Variant) int {
	start := int64(v.Pos - 1)
	fa, err := p.fai.Fetch(v.Chrom, start, start+1)
	if err != nil {
		if p.fai.Length(v.Chrom) < 0 {
			if p.stderr != nil {
				fmt.Fprintf(p.stderr, "Ignoring sequence \"%s\"\n", v.Chrom)
			}
			p.skipRID = v.Chrom
			return -2
		}
		return -3 // hard error sentinel
	}
	return fixrefNt2int(fa[0])
}

// setRefAlt rewrites REF/ALT to the single-base ref/alt and, when swapGT is
// set, swaps allele 0<->1 in every sample GT, mirroring set_ref_alt().
func (p *fixrefPlugin) setRefAlt(v *vcf.Variant, ref, alt byte, swapGT bool) {
	v.Ref = string([]byte{ref})
	if len(v.Alt) > 0 {
		v.Alt[0] = string([]byte{alt})
	}
	if !swapGT {
		return
	}
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		changed := false
		for j, a := range gt.alleles {
			if a == 0 {
				gt.alleles[j] = 1
				changed = true
			} else if a == 1 {
				gt.alleles[j] = 0
				changed = true
			}
		}
		if changed {
			v.Samples[i].Data["GT"] = gt.String()
		}
	}
}

// Process applies the configured mode to one record and (for non-stats modes)
// annotates it with the FIXREF change. It returns zero records when the record
// is discarded.
func (p *fixrefPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	p.dirty = 0
	keep, err := p.processRecord(v)
	if err != nil {
		return nil, err
	}
	if keep && p.dirty != 0 && p.mode != fixrefModeStats {
		var b strings.Builder
		for i := 0; i < 6; i++ {
			if p.dirty&(1<<i) == 0 {
				continue
			}
			if b.Len() > 0 {
				b.WriteByte(',')
			}
			b.WriteString(fixrefInfoAnnots[i])
		}
		setInfo(v, p.infoTag, b.String())
	}
	if !keep {
		return nil, nil
	}
	if p.mode == fixrefModeStats {
		return nil, nil
	}
	return []*vcf.Variant{v}, nil
}

// emit folds the discard flag for a non-converted record: it returns false
// (drop) when discarding, otherwise keep.
func (p *fixrefPlugin) emit() bool { return !p.discard }

// processRecord is the core dispatch mirroring fixref.c process_record. It
// returns whether to keep the record (false drops it), setting p.dirty as a
// side effect.
func (p *fixrefPlugin) processRecord(v *vcf.Variant) (bool, error) {
	if p.skipRID != "" && v.Chrom == p.skipRID {
		return false, nil
	}

	p.nsite++
	p.dirty = fixFixSkip

	// Skip non-SNPs (must be exactly a SNP record, any ALT non-SNP excludes it).
	if variantTypeMask(v) != vtSNP {
		p.nonSNP++
		return p.emit(), nil
	}

	ir := p.fetchRefBase(v)
	if ir == -3 {
		return false, fmt.Errorf("fixref: faidx fetch failed at %s:%d", v.Chrom, v.Pos)
	}
	if ir == -2 {
		return false, nil
	}
	if ir == -1 {
		p.nonACGT++
		return p.emit(), nil
	}

	if len(v.Alt)+1 != 2 {
		p.nonbiallelic++
		return p.emit(), nil
	}

	ia := -1
	if len(v.Ref) > 0 {
		ia = fixrefNt2int(v.Ref[0])
	}
	if ia < 0 {
		p.nonACGT++
		return p.emit(), nil
	}
	ib := -1
	if len(v.Alt[0]) > 0 {
		ib = fixrefNt2int(v.Alt[0][0])
	}
	if ib < 0 {
		p.nonACGT++
		return p.emit(), nil
	}
	if ia == ib {
		p.nonSNP++
		return p.emit(), nil
	}
	p.count[ia][ib]++
	if ir == ia {
		p.nok++
	}

	switch p.mode {
	case fixrefModeUseID:
		return p.applyUseID(v, ir, ia, ib)
	case fixrefModeRefAlt:
		return p.applyRefAlt(v, ir, ia, ib), nil
	case fixrefModeSwapRefAlt:
		return p.applySwapOnly(v, ir, ia, ib), nil
	case fixrefModeFlip2Fwd, fixrefModeFlipAll:
		return p.applyFlip(v, ir, ia, ib), nil
	case fixrefModeTop2Fwd:
		return p.applyTop(v, ir, ia, ib)
	}
	return p.emit(), nil
}

// applyRefAlt is MODE_REF_ALT: swap or flip REF/ALT to match the reference,
// leaving genotypes unchanged.
func (p *fixrefPlugin) applyRefAlt(v *vcf.Variant, ir, ia, ib int) bool {
	switch {
	case ir == ia:
		p.dirty = fixFixNone
	case ir == ib:
		p.dirty = fixFixSwap
		p.nswap++
		p.setRefAlt(v, fixrefInt2nt(ib), fixrefInt2nt(ia), false)
	case ir == fixrefRevint(ia):
		p.dirty = fixFixFlip
		p.nflip++
		p.setRefAlt(v, fixrefInt2nt(fixrefRevint(ia)), fixrefInt2nt(fixrefRevint(ib)), false)
	case ir == fixrefRevint(ib):
		p.dirty = fixFixFlip | fixFixSwap
		p.nflipSwap++
		p.setRefAlt(v, fixrefInt2nt(fixrefRevint(ib)), fixrefInt2nt(fixrefRevint(ia)), false)
	default:
		p.dirty = fixFixErr
		p.nerr++
	}
	return true
}

// applySwapOnly is MODE_SWAP_REF_ALT: only swap REF/ALT, never flip.
func (p *fixrefPlugin) applySwapOnly(v *vcf.Variant, ir, ia, ib int) bool {
	switch {
	case ir == ia:
		p.dirty = fixFixNone
	case ir == ib:
		p.dirty = fixFixSwap
		p.nswap++
		p.setRefAlt(v, fixrefInt2nt(ib), fixrefInt2nt(ia), false)
	default:
		p.dirty = fixFixErr
		p.nerr++
	}
	return true
}

// applyFlip is MODE_FLIP2FWD / MODE_FLIP_ALL: swap or flip REF/ALT and GTs for
// non-ambiguous SNPs (flip2fwd) or all SNPs (flip-all).
func (p *fixrefPlugin) applyFlip(v *vcf.Variant, ir, ia, ib int) bool {
	pair := 1<<ia | 1<<ib
	if p.mode == fixrefModeFlip2Fwd && (pair == 0x9 || pair == 0x6) {
		p.nunresolved++
		if p.discard {
			return false
		}
		// dirty stays fixFixSkip from processRecord -> annotated "skip"? No:
		// upstream leaves dirty as set by the caller. Here dirty==FIX_SKIP was
		// set before; but the C code never touched dirty for this branch, so it
		// remains FIX_SKIP. Match that.
		p.dirty = fixFixSkip
		return true
	}
	switch {
	case ir == ia:
		p.dirty = fixFixNone
	case ir == ib:
		p.dirty = fixFixSwap | fixFixGT
		p.nswap++
		p.setRefAlt(v, fixrefInt2nt(ib), fixrefInt2nt(ia), true)
	case ir == fixrefRevint(ia):
		p.dirty = fixFixFlip
		p.nflip++
		p.setRefAlt(v, fixrefInt2nt(fixrefRevint(ia)), fixrefInt2nt(fixrefRevint(ib)), false)
	case ir == fixrefRevint(ib):
		p.dirty = fixFixFlip | fixFixSwap | fixFixGT
		p.nflipSwap++
		p.setRefAlt(v, fixrefInt2nt(fixrefRevint(ib)), fixrefInt2nt(fixrefRevint(ia)), true)
	default:
		p.dirty = fixFixErr
		p.nerr++
	}
	return true
}

// applyTop is MODE_TOP2FWD: convert from Illumina TOP strand to fwd, performing
// sequence walking for ambiguous A/T and C/G pairs.
func (p *fixrefPlugin) applyTop(v *vcf.Variant, ir, ia, ib int) (bool, error) {
	pair := 1<<ia | 1<<ib
	if pair != 0x9 && pair != 0x6 {
		// unambiguous pair: A/C or A/G
		if ir == ia {
			return true, nil
		}
		iaRev := fixrefRevint(ia)
		if ir == iaRev {
			p.dirty = fixFixFlip
			p.nflip++
			p.setRefAlt(v, fixrefInt2nt(iaRev), fixrefInt2nt(fixrefRevint(ib)), false)
			return true, nil
		}
		if ir == ib {
			p.dirty = fixFixSwap | fixFixGT
			p.nswap++
			p.setRefAlt(v, fixrefInt2nt(ib), fixrefInt2nt(ia), true)
			return true, nil
		}
		if ib != fixrefRevint(ir) {
			p.dirty = fixFixErr
			p.nerr++
			return true, nil
		}
		p.dirty = fixFixFlip | fixFixSwap | fixFixGT
		p.nflipSwap++
		p.setRefAlt(v, fixrefInt2nt(fixrefRevint(ib)), fixrefInt2nt(fixrefRevint(ia)), true)
		return true, nil
	}

	// ambiguous pair: sequence walking around the position.
	win := v.Pos - 1
	if win > 100 {
		win = 100
	}
	beg := int64(v.Pos - 1 - win)
	end := int64(v.Pos-1+win) + 1 // half-open
	ref, err := p.fai.Fetch(v.Chrom, beg, end)
	if err != nil {
		return false, fmt.Errorf("fixref: faidx fetch failed at %s:%d: %w", v.Chrom, v.Pos, err)
	}
	mid := win // index of the record position within ref
	strand := 0
	for i := 1; i <= win; i++ {
		ra := fixrefNt2int(ref[mid-i])
		rb := fixrefNt2int(ref[mid+i])
		if ra < 0 || rb < 0 || ra == rb {
			continue
		}
		pr := 1<<ra | 1<<rb
		if pr == 0x9 || pr == 0x6 {
			continue
		}
		if (1<<ra)&0x9 != 0 {
			strand = 1
		} else {
			strand = -1
		}
		break
	}

	if strand == 1 {
		if ir == ia {
			p.dirty = fixFixNone
			return true, nil
		}
		if ir == ib {
			p.dirty = fixFixSwap | fixFixGT
			p.nswap++
			p.setRefAlt(v, fixrefInt2nt(ib), fixrefInt2nt(ia), true)
			return true, nil
		}
	} else if strand == -1 {
		iaRev := fixrefRevint(ia)
		ibRev := fixrefRevint(ib)
		if ir == iaRev {
			p.dirty = fixFixFlip
			p.nflip++
			p.setRefAlt(v, fixrefInt2nt(iaRev), fixrefInt2nt(ibRev), false)
			return true, nil
		}
		if ir == ibRev {
			p.dirty = fixFixFlip | fixFixSwap | fixFixGT
			p.nflipSwap++
			p.setRefAlt(v, fixrefInt2nt(ibRev), fixrefInt2nt(iaRev), true)
			return true, nil
		}
	}

	p.nunresolved++
	if p.discard {
		return false, nil
	}
	p.dirty = fixFixSkip
	return true, nil
}

// Destroy prints the upstream stats block to stderr. The stats are not part of
// the parity-checked stdout, but they are emitted to match upstream behaviour.
func (p *fixrefPlugin) Destroy() error {
	if p.fai != nil {
		_ = p.fai.Close()
	}
	if p.stderr == nil {
		return nil
	}
	w := p.stderr

	var topMask = [4][4]int{
		{0, 1, 1, 1},
		{0, 0, 1, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
	}
	var botMask = [4][4]int{
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 1, 0, 0},
		{1, 1, 1, 0},
	}
	var tot uint32
	var topErr, botErr uint32
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			tot += p.count[i][j]
			if topMask[i][j] == 0 && p.count[i][j] != 0 {
				topErr++
			}
			if botMask[i][j] == 0 && p.count[i][j] != 0 {
				botErr++
			}
		}
	}
	nskip := p.nonACGT + p.nonSNP + p.nonbiallelic
	ncmp := p.nsite - nskip

	pct := func(num, den uint32) float64 {
		if den == 0 {
			return 0
		}
		return 100 * float64(num) / float64(den)
	}

	fmt.Fprintf(w, "# SC, guessed strand convention\n")
	scTop, scBot := 1, 1
	if topErr != 0 {
		scTop = 0
	}
	if botErr != 0 {
		scBot = 0
	}
	fmt.Fprintf(w, "SC\tTOP-compatible\t%d\n", scTop)
	fmt.Fprintf(w, "SC\tBOT-compatible\t%d\n", scBot)

	fmt.Fprintf(w, "# ST, substitution types\n")
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if i == j {
				continue
			}
			var p100 float64
			if tot != 0 {
				p100 = float64(p.count[i][j]) * 100 / float64(tot)
			}
			fmt.Fprintf(w, "ST\t%c>%c\t%d\t%.1f%%\n", fixrefInt2nt(i), fixrefInt2nt(j), p.count[i][j], p100)
		}
	}
	fmt.Fprintf(w, "# NS, Number of sites:\n")
	fmt.Fprintf(w, "NS\ttotal        \t%d\n", p.nsite)
	fmt.Fprintf(w, "NS\tref match    \t%d\t%.1f%%\n", p.nok, pct(p.nok, ncmp))
	fmt.Fprintf(w, "NS\tref mismatch \t%d\t%.1f%%\n", ncmp-p.nok, pct(ncmp-p.nok, ncmp))
	if p.mode != fixrefModeStats {
		den := p.nsite - nskip
		fmt.Fprintf(w, "NS\tflipped      \t%d\t%.1f%%\n", p.nflip, pct(p.nflip, den))
		fmt.Fprintf(w, "NS\tswapped      \t%d\t%.1f%%\n", p.nswap, pct(p.nswap, den))
		fmt.Fprintf(w, "NS\tflip+swap    \t%d\t%.1f%%\n", p.nflipSwap, pct(p.nflipSwap, den))
		fmt.Fprintf(w, "NS\tunresolved   \t%d\t%.1f%%\n", p.nunresolved, pct(p.nunresolved, den))
		fmt.Fprintf(w, "NS\tfixed pos    \t%d\t%.1f%%\n", p.npesErr, pct(p.npesErr, den))
	}
	fmt.Fprintf(w, "NS\terrors       \t%d\n", p.nerr)
	fmt.Fprintf(w, "NS\tskipped      \t%d\n", nskip)
	fmt.Fprintf(w, "NS\tnon-ACGT     \t%d\n", p.nonACGT)
	fmt.Fprintf(w, "NS\tnon-SNP      \t%d\n", p.nonSNP)
	fmt.Fprintf(w, "NS\tnon-biallelic\t%d\n", p.nonbiallelic)
	return nil
}
