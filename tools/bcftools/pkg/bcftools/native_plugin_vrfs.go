// Native port of the upstream `vrfs` plugin (plugins/vrfs.c). vrfs (variant
// read frequency score) assesses site noisiness by piling up a large set of
// BAM/CRAM alignments against a FASTA reference: for every indexed site it
// walks the reads, bins each sample's variant-allele fraction (VAF) into a
// per-site histogram, then derives a per-site score from the across-sample
// variance and emits the SITE/MEAN/VAR2 profile.
//
// Pileup semantics replicated (the exact mpileup2 LEGACY_MODE configuration
// vrfs.c sets via mpileup_set):
//
//   - MIN_MQ 0, SKIP_ANY_SET {UNMAP,SECONDARY,QCFAIL,DUP}: a read is dropped
//     only if unmapped, secondary, QC-fail or duplicate, or has tid<0. No
//     mapping-quality floor, no proper-pair / orphan filtering, no
//     supplementary drop (vrfs does NOT set BAM_FSUPPLEMENTARY).
//   - LEGACY_MODE 1: crucially, upstream's legacy_mplp_func has the BAQ /
//     realignment step stubbed out ("// realign" in mpileup2/mpileup.c). The
//     legacy pileup therefore applies NO BAQ and NO base-quality adjustment,
//     and the vrfs counting loop reads bases with no base-quality floor
//     (MIN_BQ is irrelevant — it never consults the qual). So a pure-SNP
//     pileup is byte-exact reproducible without the baq engine. The MAX_BQ /
//     DELTA_BQ / MIN_REALN_* / MAX_DP_PER_SAMPLE knobs vrfs sets only affect
//     the (unused) new-mode path; in legacy mode they are no-ops for the
//     count, which is why this port reaches byte-for-byte parity on SNV sites
//     without porting them.
//   - At each site, per sample: total = reads whose base is ref or a
//     non-ref base or an indel; nalt[b] counts reads supporting alt base b
//     (A/C/G/T -> 0..3) or a generic indel (-> 4). A read with an indel
//     immediately after the site position counts as an indel alt; otherwise
//     the read base is compared to the reference base. Reads with total <
//     min-depth are skipped. nn2bin maps (nref,nalt) -> VAF bin.
//
// vrfs reads no VCF/BCF input — it is registered as a fullPlugin and owns the
// whole invocation (reading the alignment list, the sites file and the FASTA,
// then writing the text profile). The implementation reuses pkg/htsgo (sam,
// alnio, fasta) for BAM/CRAM/FASTA access.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("vrfs", func() NativePlugin { return &vrfsPlugin{} }) }

// vrfsDefaultNBins is the default number of VAF bins (vrfs.c N_BINS).
const vrfsDefaultNBins = 20

// vrfsDefaultMinDepth is the default minimum read depth (vrfs.c args->min_dp).
const vrfsDefaultMinDepth = 10

// vrfsHardcodedVar2 is the pre-computed per-bin variance profile from vrfs.c
// (the var2[N_BINS] table). Used by the default "hc" recalc mode and rescaled
// when --nbins differs from 20.
var vrfsHardcodedVar2 = []float64{
	1, 1.441327e-03, 4.382657e-05, 6.160600e-07, 1.414270e-08, 2.828540e-09,
	2.357117e-09, 2.020386e-09, 1.767838e-09, 1.571411e-09, 1.414270e-09,
	1.285700e-09, 1.178558e-09, 1.087900e-09, 1.010193e-09, 9.428468e-10,
	8.839188e-10, 8.319236e-10, 7.857056e-10, 7.443527e-10,
}

// vrfsPlugin is the native entry for `+vrfs`.
type vrfsPlugin struct{}

// Name returns the plugin name.
func (p *vrfsPlugin) Name() string { return "vrfs" }

// About returns the one-line description, matching vrfs.c about().
func (p *vrfsPlugin) About() string {
	return "Localised assessment of sequencing artefacts, estimate site noisiness (variant read frequency score)\n"
}

// RunStyle reports that vrfs is a run()-style plugin: upstream exports a `run`
// symbol, so its options precede any input with no `--` separator.
func (p *vrfsPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of vrfs's own flags consumes the following
// CLI token as its value. -i/--use-index is a boolean toggle; -v/--verbosity is
// optional-argument (its value, when present, is attached as -v3, not a
// separate token), so neither consumes the next token.
func (p *vrfsPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-a", "--alns", "-f", "--fasta-ref", "-s", "--sites",
		"-d", "--min-depth", "-n", "--nbins", "-r", "--recalc",
		"-b", "--batch", "-m", "--merge-batches", "-M", "--merge-files",
		"-o", "--output", "-O", "--output-type":
		return true
	}
	return false
}

// RunFull executes vrfs end to end: it parses vrfs's own argv, then either
// builds a profile from alignments (-a) or merges batch files (-m/-M), writing
// the text profile to out.
func (p *vrfsPlugin) RunFull(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	cfg, err := parseVrfsArgs(opts.Args)
	if err != nil {
		return err
	}

	// --batch k=N just reports the number of batches and exits, before reading
	// anything but the alignment list.
	if cfg.batch != "" && strings.HasPrefix(cfg.batch, "k=") {
		return runVrfsBatchCount(cfg, out)
	}

	doMerge := cfg.batchFile != "" || len(cfg.mergeFiles) > 0
	doProfile := cfg.alnFile != ""
	if doMerge && doProfile {
		return fmt.Errorf("vrfs: cannot both build a profile (-a) and merge (-m/-M)")
	}
	if !doMerge && !doProfile {
		return fmt.Errorf("vrfs: missing required options (need -a for profiling or -m/-M for merging)")
	}

	if doProfile {
		return runVrfsProfile(cfg, out, stderr)
	}
	return runVrfsMerge(cfg, out)
}

// Init is required by the NativePlugin interface but unused: RunFull owns the
// whole invocation, so vrfs never enters the per-record pipeline.
func (p *vrfsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return hdr, nil
}

// Process is never reached (RunFull owns the invocation).
func (p *vrfsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *vrfsPlugin) Destroy() error { return nil }

// vrfsConfig holds the parsed vrfs options.
type vrfsConfig struct {
	alnFile    string
	fastaFile  string
	sitesFile  string
	output     string
	outputGz   bool
	clevel     int
	minDepth   int
	nbins      int
	recalc     string // "hc", "data", or "file:/path"
	useBamIdx  bool
	batch      string
	batchFile  string   // -m
	mergeFiles []string // -M and positionals
	verbose    int
}

// parseVrfsArgs parses vrfs's getopt-style argv into a vrfsConfig, mirroring
// vrfs.c run()'s option handling and defaults.
func parseVrfsArgs(args []string) (vrfsConfig, error) {
	cfg := vrfsConfig{
		output:   "-",
		clevel:   -1,
		minDepth: vrfsDefaultMinDepth,
		nbins:    vrfsDefaultNBins,
		recalc:   "hc",
	}
	outputTypeSet := false
	want := func(i int) (string, error) {
		if i+1 >= len(args) {
			return "", fmt.Errorf("vrfs: missing argument after %q", args[i])
		}
		return args[i+1], nil
	}
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch {
		case a == "-d" || a == "--min-depth":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			n, perr := strconv.Atoi(v)
			if perr != nil || n < 0 {
				return cfg, fmt.Errorf("vrfs: could not parse argument: --min-depth %s", v)
			}
			cfg.minDepth = n
		case a == "-n" || a == "--nbins":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			n, perr := strconv.Atoi(v)
			if perr != nil || n < 10 {
				return cfg, fmt.Errorf("vrfs: could not parse argument: --nbins %s; the minimum value is 10", v)
			}
			cfg.nbins = n
		case a == "-i" || a == "--use-index":
			cfg.useBamIdx = true
		case a == "-r" || a == "--recalc":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.recalc = v
		case a == "-f" || a == "--fasta-ref":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.fastaFile = v
		case a == "-o" || a == "--output":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.output = v
		case a == "-O" || a == "--output-type":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			outputTypeSet = true
			if err := parseVrfsOutputType(&cfg, v); err != nil {
				return cfg, err
			}
		case a == "-s" || a == "--sites":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.sitesFile = v
		case a == "-a" || a == "--alns":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.alnFile = v
		case a == "-b" || a == "--batch":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.batch = v
		case a == "-m" || a == "--merge-batches":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.batchFile = v
		case a == "-M" || a == "--merge-files":
			v, err := want(i)
			if err != nil {
				return cfg, err
			}
			i++
			cfg.mergeFiles = append(cfg.mergeFiles, v)
		case a == "-v" || a == "--verbose" || a == "--verbosity":
			cfg.verbose++
		case strings.HasPrefix(a, "-v") && len(a) > 2:
			n, perr := strconv.Atoi(a[2:])
			if perr != nil || n < 0 {
				return cfg, fmt.Errorf("vrfs: could not parse argument: --verbosity %s", a[2:])
			}
			cfg.verbose = n
		case a == "--":
			cfg.mergeFiles = append(cfg.mergeFiles, args[i+1:]...)
			i = len(args)
		case len(a) > 0 && a[0] == '-' && a != "-":
			return cfg, fmt.Errorf("vrfs: unrecognized option %q", a)
		default:
			// trailing positional = additional merge file (-M ... FILE)
			cfg.mergeFiles = append(cfg.mergeFiles, a)
		}
	}
	// Infer output type from filename when -O was not given (.gz -> compressed).
	if !outputTypeSet {
		if len(cfg.output) > 3 && strings.EqualFold(cfg.output[len(cfg.output)-3:], ".gz") {
			cfg.outputGz = true
		}
	}
	return cfg, nil
}

// parseVrfsOutputType parses the -O t|z[0-9] argument.
func parseVrfsOutputType(cfg *vrfsConfig, v string) error {
	if v == "" {
		return fmt.Errorf("vrfs: the output type %q not recognised", v)
	}
	switch v[0] {
	case 't':
		cfg.outputGz = false
	case 'z':
		cfg.outputGz = true
	default:
		n, err := strconv.Atoi(v)
		if err != nil || n < 0 || n > 9 {
			return fmt.Errorf("vrfs: the output type %q not recognised", v)
		}
		cfg.clevel = n
		cfg.outputGz = true
	}
	if len(v) > 1 {
		n, err := strconv.Atoi(v[1:])
		if err != nil || n < 0 || n > 9 {
			return fmt.Errorf("vrfs: could not parse argument: --compression-level %s", v[1:])
		}
		cfg.clevel = n
	}
	return nil
}

// vrfsSite is one indexed site: a (chrom,pos,ref,alt) tuple with its VAF
// histogram. Indels are preserved verbatim in ref/alt but treated as a single
// generic "indel" type internally; SNVs keep only the first ref/alt base.
type vrfsSite struct {
	chrom string
	pos0  int // 0-based position
	ref   string
	alt   string
	dist  []uint32 // VAF histogram, length nbins
	nval  int      // number of values contributing
}

// isIndel reports whether the site is an indel (ref or alt has length > 1
// after the SNV normalisation parseVrfsSites applies).
func (s *vrfsSite) isIndel() bool { return len(s.ref) > 1 || len(s.alt) > 1 }

// altClass returns the 0..4 alt-allele index for the site: A/C/G/T -> 0..3,
// indel -> 4. It mirrors vrfs.c's seq_nt16 encoding of the alt base.
func (s *vrfsSite) altClass() int {
	if s.isIndel() {
		return 4
	}
	return baseClass(s.alt[0])
}

// baseClass maps a base byte to 0..3 (A/C/G/T), or -1 for anything else.
// Matches htslib's seq_nt16_int[seq_nt16_table[base]] for A/C/G/T (case
// insensitive).
func baseClass(b byte) int {
	switch b {
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

// parseVrfsSites parses the tab/space-delimited sites file ("chr pos ref alt",
// 1-based pos) into vrfsSite records, mirroring vrfs.c parse_sites. Blank lines
// and lines starting with '#' are skipped. When ref and alt have equal length
// (a putative SNV/MNV), only the first ref/alt base is retained (vrfs treats it
// as an SNV). Each site's dist is allocated with nbins slots.
func parseVrfsSites(r io.Reader, nbins int) ([]vrfsSite, error) {
	var sites []vrfsSite
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 4<<20)
	for sc.Scan() {
		line := sc.Text()
		trimmed := strings.TrimLeft(line, " \t")
		if trimmed == "" || trimmed[0] == '#' {
			continue
		}
		fields := strings.Fields(trimmed)
		if len(fields) < 2 {
			return nil, fmt.Errorf("vrfs: could not parse the CHR/POS part of the line: %s", line)
		}
		chrom := fields[0]
		pos1, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("vrfs: could not parse the POS part of the line: %s", line)
		}
		if len(fields) < 3 || fields[2] == "" {
			return nil, fmt.Errorf("vrfs: could not parse the REF part of the line: %s", line)
		}
		if len(fields) < 4 || fields[3] == "" {
			return nil, fmt.Errorf("vrfs: could not parse the ALT part of the line: %s", line)
		}
		ref := fields[2]
		alt := fields[3]
		if len(ref) == len(alt) {
			// treat as SNV: keep only the first base of each
			ref = ref[:1]
			alt = alt[:1]
		}
		sites = append(sites, vrfsSite{
			chrom: chrom,
			pos0:  pos1 - 1,
			ref:   ref,
			alt:   alt,
			dist:  make([]uint32, nbins),
		})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return sites, nil
}

// parseAlnList reads a one-path-per-line alignment list (vrfs.c hts_readlist).
// Lines are trimmed of surrounding whitespace; blank lines are skipped.
func parseAlnList(r io.Reader) ([]string, error) {
	var out []string
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 4<<20)
	for sc.Scan() {
		ln := strings.TrimSpace(sc.Text())
		if ln == "" {
			continue
		}
		out = append(out, ln)
	}
	return out, sc.Err()
}

// nn2bin maps a (nref, nalt) read count to a VAF bin in [0, nbins), mirroring
// vrfs.c nn2bin. -1 means no reads (skip). bin 0 is "no alt"; otherwise the bin
// is floor((nbins-1) * nalt / (nref+nalt)).
func nn2bin(nbins, nref, nalt int) int {
	if nalt == 0 && nref == 0 {
		return -1
	}
	if nalt == 0 {
		return 0
	}
	return int(float64(nbins-1) * float64(nalt) / float64(nref+nalt))
}

// vrfsSiteIndex indexes sites by (chrom, pos0) for fast overlap lookup during
// the pileup, preserving the sorted iteration order vrfs.c's regidx produces
// for output.
type vrfsSiteIndex struct {
	sites []*vrfsSite
	byPos map[string]map[int][]*vrfsSite // chrom -> pos0 -> sites
	nbins int
}

// buildVrfsIndex builds a position index over sites, returning the sites sorted
// by (chrom, pos0) — the order vrfs.c emits SITE lines (regidx sorted order).
func buildVrfsIndex(sites []vrfsSite, nbins int) *vrfsSiteIndex {
	idx := &vrfsSiteIndex{
		byPos: map[string]map[int][]*vrfsSite{},
		nbins: nbins,
	}
	idx.sites = make([]*vrfsSite, len(sites))
	for i := range sites {
		idx.sites[i] = &sites[i]
	}
	// Stable sort by chrom (first-seen order) then position. vrfs.c's regidx
	// groups by sequence in insertion order and sorts positions within a
	// sequence; replicate that by tracking chrom first-appearance order.
	chromOrder := map[string]int{}
	next := 0
	for _, s := range idx.sites {
		if _, ok := chromOrder[s.chrom]; !ok {
			chromOrder[s.chrom] = next
			next++
		}
	}
	sort.SliceStable(idx.sites, func(i, j int) bool {
		ci, cj := chromOrder[idx.sites[i].chrom], chromOrder[idx.sites[j].chrom]
		if ci != cj {
			return ci < cj
		}
		return idx.sites[i].pos0 < idx.sites[j].pos0
	})
	for _, s := range idx.sites {
		m := idx.byPos[s.chrom]
		if m == nil {
			m = map[int][]*vrfsSite{}
			idx.byPos[s.chrom] = m
		}
		m[s.pos0] = append(m[s.pos0], s)
	}
	return idx
}

// runVrfsBatchCount handles --batch k=N: read the aln list and print the number
// of batches needed, matching vrfs.c.
func runVrfsBatchCount(cfg vrfsConfig, out io.Writer) error {
	f, err := os.Open(cfg.alnFile)
	if err != nil {
		return fmt.Errorf("vrfs: open aln list %q: %w", cfg.alnFile, err)
	}
	defer f.Close()
	bams, err := parseAlnList(f)
	if err != nil {
		return err
	}
	k, perr := strconv.Atoi(strings.TrimPrefix(cfg.batch, "k="))
	if perr != nil || k <= 0 {
		return fmt.Errorf("vrfs: could not parse --batch %s", cfg.batch)
	}
	fmt.Fprintf(out, "# Number of required batches with %d files total and max %d files per batch:\n", len(bams), k)
	fmt.Fprintf(out, "%.0f\n", math.Ceil(float64(len(bams))/float64(k)))
	return nil
}

// selectBatch applies --batch I/N to the bam list, returning the I-th slice.
func selectBatch(bams []string, batch string) ([]string, error) {
	parts := strings.SplitN(batch, "/", 2)
	if len(parts) != 2 {
		return nil, fmt.Errorf("vrfs: could not parse --batch %s", batch)
	}
	ith, e1 := strconv.Atoi(parts[0])
	nbatches, e2 := strconv.Atoi(parts[1])
	if e1 != nil || ith <= 0 {
		return nil, fmt.Errorf("vrfs: could not parse --batch %s", batch)
	}
	if e2 != nil || nbatches <= 0 {
		return nil, fmt.Errorf("vrfs: could not parse --batch %s", batch)
	}
	if ith > len(bams) {
		return nil, fmt.Errorf("vrfs: asked for %d-th batch in a list of %d files", ith, len(bams))
	}
	if ith > nbatches {
		return nil, fmt.Errorf("vrfs: the batch index is outside the permitted range [1,%d]", nbatches)
	}
	if nbatches > len(bams) {
		return nil, fmt.Errorf("vrfs: cannot create %d batches from a list of %d files", nbatches, len(bams))
	}
	nper := int(math.Ceil(float64(len(bams)) / float64(nbatches)))
	isrc := (ith - 1) * nper
	if isrc >= len(bams) {
		return nil, nil // empty batch
	}
	end := isrc + nper
	if end > len(bams) {
		end = len(bams)
	}
	return bams[isrc:end], nil
}

// runVrfsProfile builds the VAF profile from the alignment list and writes it.
func runVrfsProfile(cfg vrfsConfig, out io.Writer, stderr io.Writer) error {
	if cfg.fastaFile == "" || cfg.sitesFile == "" {
		return fmt.Errorf("vrfs: -f/--fasta-ref and -s/--sites are required for profiling")
	}

	af, err := os.Open(cfg.alnFile)
	if err != nil {
		return fmt.Errorf("vrfs: open aln list %q: %w", cfg.alnFile, err)
	}
	bams, err := parseAlnList(af)
	af.Close()
	if err != nil {
		return err
	}

	if cfg.batch != "" {
		bams, err = selectBatch(bams, cfg.batch)
		if err != nil {
			return err
		}
	}

	sf, err := os.Open(cfg.sitesFile)
	if err != nil {
		return fmt.Errorf("vrfs: open sites %q: %w", cfg.sitesFile, err)
	}
	sites, err := parseVrfsSites(sf, cfg.nbins)
	sf.Close()
	if err != nil {
		return err
	}
	idx := buildVrfsIndex(sites, cfg.nbins)

	ref, err := fasta.OpenRandomAccess(cfg.fastaFile)
	if err != nil {
		return fmt.Errorf("vrfs: open reference %q: %w", cfg.fastaFile, err)
	}
	defer ref.Close()

	for _, bam := range bams {
		if err := vrfsPileupBam(bam, idx, ref, cfg.minDepth); err != nil {
			return err
		}
	}

	prof := computeVrfsProfile(idx.sites, cfg)
	return writeVrfsProfile(out, idx.sites, prof, cfg)
}

// vrfsPileupBam streams one BAM/CRAM and accumulates per-sample VAF counts into
// the indexed sites. It mirrors vrfs.c's batch_profile_run1 over the legacy
// pileup: read-level flag filtering (drop unmapped/secondary/qcfail/dup), no
// MAPQ floor, no BAQ, then per-site per-sample ref/alt counting.
func vrfsPileupBam(path string, idx *vrfsSiteIndex, ref *fasta.RandomAccess, minDepth int) error {
	f, err := os.Open(path)
	if err != nil {
		return fmt.Errorf("vrfs: open alignment %q: %w", path, err)
	}
	defer f.Close()
	rd, err := alnio.NewReaderWithReference(f, "")
	if err != nil {
		return fmt.Errorf("vrfs: read %q: %w", path, err)
	}
	hdr := rd.Header()

	// Determine the sample for every read group, defaulting to the bam path
	// when a read has no RG / the bam has no @RG (mirrors bam_smpl_add_bam).
	rg2smpl, defaultSmpl := vrfsSampleMap(hdr, path)
	acc := newVrfsAccumulator()

	refCache := map[string][]byte{}
	refBaseAt := func(chrom string, pos0 int) (byte, bool) {
		seq, ok := refCache[chrom]
		if !ok {
			length := ref.Length(chrom)
			if length <= 0 {
				refCache[chrom] = nil
				return 0, false
			}
			b, ferr := ref.Fetch(chrom, 0, length)
			if ferr != nil {
				refCache[chrom] = nil
				return 0, false
			}
			seq = b
			refCache[chrom] = seq
		}
		if seq == nil || pos0 < 0 || pos0 >= len(seq) {
			return 0, false
		}
		return seq[pos0], true
	}

	for {
		rec, rerr := rd.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return fmt.Errorf("vrfs: read %q: %w", path, rerr)
		}
		if !vrfsKeepRead(rec) {
			continue
		}
		m := idx.byPos[rec.RName]
		if m == nil {
			continue
		}
		smpl := vrfsReadSample(rec, rg2smpl, defaultSmpl)
		accumulateReadAtSites(rec, m, smpl, refBaseAt, acc)
	}

	acc.flush(idx.nbins, minDepth)
	return nil
}

// vrfsKeepRead applies vrfs's read-level filter: drop tid<0 (unmapped RName),
// unmapped, secondary, QC-fail, duplicate. There is no MAPQ floor (MIN_MQ 0)
// and supplementary reads are NOT dropped (vrfs.c does not set
// BAM_FSUPPLEMENTARY in SKIP_ANY_SET).
func vrfsKeepRead(rec *sam.Record) bool {
	if rec.Pos <= 0 || rec.RName == "" {
		return false
	}
	if rec.Flag&(sam.FlagUnmapped|sam.FlagSecondary|sam.FlagQCFail|sam.FlagDuplicate) != 0 {
		return false
	}
	return true
}

// vrfsSampleMap builds a read-group-id -> sample-name map from the @RG lines,
// plus the default sample used for reads with no RG (or when the bam has no
// usable read groups): the bam file name, matching htslib bam_smpl_add_bam.
func vrfsSampleMap(hdr *sam.Header, path string) (map[string]string, string) {
	rg2smpl := map[string]string{}
	for _, rg := range hdr.ReadGroups {
		sm := path
		for _, f := range rg.Extra {
			if f.Tag == "SM" {
				sm = f.Value
				break
			}
		}
		rg2smpl[rg.ID] = sm
	}
	return rg2smpl, path
}

// vrfsReadSample resolves a read's sample name from its RG:Z aux tag.
func vrfsReadSample(rec *sam.Record, rg2smpl map[string]string, defaultSmpl string) string {
	if len(rg2smpl) == 0 {
		return defaultSmpl
	}
	if a, ok := rec.GetAux("RG"); ok {
		if id, ok := a.String(); ok {
			if sm, ok := rg2smpl[id]; ok {
				return sm
			}
		}
	}
	return defaultSmpl
}

// vrfsPosCounts holds, for one (sample, site-position) cell, the per-base and
// indel supporting-read counts plus the running total. It mirrors the ntot /
// nalt[5] accounting in vrfs.c's per-sample loop.
type vrfsPosCounts struct {
	ntot int
	nalt [5]int // A,C,G,T,indel
	site [5]*vrfsSite
}

// vrfsAccumulator gathers per-(sample, position) read counts. vrfs counts
// contributions per sample, then bins once per (sample, site). The key is
// "sample\x00chrom\x00pos0".
type vrfsAccumulator struct {
	cells map[string]*vrfsPosCounts
}

func newVrfsAccumulator() *vrfsAccumulator {
	return &vrfsAccumulator{cells: map[string]*vrfsPosCounts{}}
}

// accumulateReadAtSites adds one read's contribution to every site position it
// covers, classifying it as ref / a non-ref base / an indel exactly as vrfs.c's
// inner loop over the legacy pileup column does. refBaseAt yields the reference
// base for the no-ref-match guard (a read column is skipped when the reference
// base differs from the site's recorded ref, matching vrfs.c).
func accumulateReadAtSites(rec *sam.Record, sitesByPos map[int][]*vrfsSite, smpl string,
	refBaseAt func(string, int) (byte, bool), acc *vrfsAccumulator) {
	for _, col := range vrfsReadColumns(rec, sitesByPos) {
		sites := sitesByPos[col.pos0]
		if len(sites) == 0 {
			continue
		}
		refb, ok := refBaseAt(rec.RName, col.pos0)
		if !ok {
			continue
		}
		// vrfs skips the whole column when the reference base does not match
		// the site's recorded REF (No ref match).
		matched := false
		for _, s := range sites {
			if len(s.ref) > 0 && s.ref[0] == refb {
				matched = true
				break
			}
		}
		if !matched {
			continue
		}

		key := smpl + "\x00" + rec.RName + "\x00" + strconv.Itoa(col.pos0)
		c := acc.cells[key]
		if c == nil {
			c = &vrfsPosCounts{}
			// keep one representative site per alt class (vrfs ialt2site)
			for _, s := range sites {
				ia := s.altClass()
				if ia >= 0 && c.site[ia] == nil {
					c.site[ia] = s
				}
			}
			acc.cells[key] = c
		}

		if col.indel {
			c.nalt[4]++
			c.ntot++
			continue
		}
		if col.base == refb {
			c.ntot++
			continue
		}
		ia := baseClass(col.base)
		if ia >= 0 {
			c.nalt[ia]++
			c.ntot++
		}
		// A base that is neither ref nor a valid A/C/G/T (e.g. N) is ignored
		// for ntot, matching vrfs's branch where only a real alt base or a
		// ref match increments ntot.
	}
}

// vrfsReadColumn is one read base at a covered site position.
type vrfsReadColumn struct {
	pos0  int
	base  byte
	indel bool // an insertion or deletion begins immediately after this base
}

// vrfsReadColumns walks rec's CIGAR and returns, for every reference position
// that coincides with a site in sitesByPos, the read base htslib's pileup
// engine would report there and whether an indel begins right after it
// (mirroring bam_pileup1_t's qpos/indel/is_del).
//
//   - M/=/X positions report the aligned read base. The last base of an M-run
//     immediately followed by an I or D op carries indel=true (htslib sets the
//     prior column's `indel` field), so vrfs classifies the read as a generic
//     indel alt there.
//   - A position consumed by a D (deletion) op is an is_del column. htslib
//     leaves qpos pointing at the query base immediately AFTER the deletion
//     (the next aligned base) and `indel`=0; vrfs reads that base with
//     bam_seqi and counts it as ref or an alt SNV base accordingly. So a
//     deletion read contributes its post-deletion base to every reference
//     position the deletion spans. We replicate this exactly. (A deletion
//     reaching the read's end leaves no next base; htslib reports the last
//     query base — handled by clamping the query index.)
//   - N (ref-skip) columns are treated the same way (rare; vrfs reads the
//     base after the skip), keeping behaviour faithful when one occurs.
func vrfsReadColumns(rec *sam.Record, sitesByPos map[int][]*vrfsSite) []vrfsReadColumn {
	var cols []vrfsReadColumn
	refPos := int(rec.Pos) - 1
	queryPos := 0
	ops := rec.Cigar
	for oi := 0; oi < len(ops); oi++ {
		op := ops[oi]
		l := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			for k := 0; k < l; k++ {
				p := refPos + k
				if _, want := sitesByPos[p]; want {
					q := queryPos + k
					var base byte = 'N'
					if q < len(rec.Seq) {
						base = rec.Seq[q]
					}
					indel := false
					// An indel begins right after this base if it is the last
					// base of this M-run and the next op is I or D.
					if k == l-1 && oi+1 < len(ops) {
						switch ops[oi+1].Op() {
						case sam.CigarInsertion, sam.CigarDeletion:
							indel = true
						}
					}
					cols = append(cols, vrfsReadColumn{pos0: p, base: base, indel: indel})
				}
			}
			refPos += l
			queryPos += l
		case sam.CigarInsertion:
			queryPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			// is_del / is_refskip columns: htslib's qpos points at the query
			// base just after the gap (queryPos here, since query has not yet
			// advanced past the M run preceding the gap). vrfs reads that base.
			q := queryPos
			if q >= len(rec.Seq) {
				q = len(rec.Seq) - 1
			}
			var base byte = 'N'
			if q >= 0 && q < len(rec.Seq) {
				base = rec.Seq[q]
			}
			for k := 0; k < l; k++ {
				p := refPos + k
				if _, want := sitesByPos[p]; want {
					cols = append(cols, vrfsReadColumn{pos0: p, base: base, indel: false})
				}
			}
			refPos += l
		case sam.CigarSoftClip:
			queryPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// consumes neither
		}
	}
	return cols
}

// flush bins every accumulated (sample, position) cell into the site
// histograms, applying the min-depth filter and nn2bin, mirroring vrfs.c's
// "increment site counters" block. Order is irrelevant: each cell increments
// independent site dist counters by integers, and addition commutes.
func (acc *vrfsAccumulator) flush(nbins, minDepth int) {
	for _, c := range acc.cells {
		if c.ntot < minDepth {
			continue
		}
		for j := 0; j < 5; j++ {
			s := c.site[j]
			if s == nil {
				continue
			}
			ifreq := nn2bin(nbins, c.ntot-c.nalt[j], c.nalt[j])
			if ifreq < 0 {
				continue
			}
			s.nval++
			s.dist[ifreq]++
		}
	}
}

// vrfsProfile holds the computed mean / var2 arrays for the SITE/MEAN/VAR2
// output.
type vrfsProfile struct {
	mean []float64
	var2 []float64
	nval int
}

// computeVrfsProfile computes the per-bin mean and var2 across all sites with
// data, mirroring vrfs.c batch_profile_set_mean_var2 + init_var2.
func computeVrfsProfile(sites []*vrfsSite, cfg vrfsConfig) vrfsProfile {
	nbins := cfg.nbins
	prof := vrfsProfile{
		mean: make([]float64, nbins),
	}

	// init_var2: returns a non-nil accumulator slice only in "data" mode.
	var var2Acc []float64
	recalc := cfg.recalc
	switch {
	case strings.HasPrefix(recalc, "file:"):
		ori := vrfsReadVar2File(recalc[len("file:"):])
		prof.var2 = vrfsRescaleVar2(ori, nbins)
	case recalc == "hc":
		prof.var2 = vrfsRescaleVar2(vrfsHardcodedVar2, nbins)
	case recalc == "data":
		prof.var2 = make([]float64, nbins)
		var2Acc = prof.var2
	default:
		// unknown -> behave like hc (defensive; CLI parsing accepts any string)
		prof.var2 = vrfsRescaleVar2(vrfsHardcodedVar2, nbins)
	}

	nval := 0
	for _, s := range sites {
		if s.nval == 0 {
			continue
		}
		maxVal := float64(s.dist[0])
		for i := 0; i < nbins; i++ {
			if v := float64(s.dist[i]); v > maxVal {
				maxVal = v
			}
		}
		for i := 0; i < nbins; i++ {
			val := float64(s.dist[i]) / maxVal
			prof.mean[i] += val
			if var2Acc != nil {
				var2Acc[i] += val * val
			}
		}
		nval++
	}
	for i := 0; i < nbins; i++ {
		prof.mean[i] = prof.mean[i] / float64(nval)
	}
	prof.nval = nval

	if var2Acc != nil && nval != 0 {
		minNonzero := 1.0
		for i := 0; i < nbins; i++ {
			var2Acc[i] = var2Acc[i]/float64(nval) - prof.mean[i]*prof.mean[i]
			if var2Acc[i] > 0 && var2Acc[i] < minNonzero {
				minNonzero = var2Acc[i]
			}
		}
		for i := 0; i < nbins; i++ {
			if var2Acc[i] == 0 {
				div := float64(i + 1)
				if i == 0 {
					div = 1
				}
				var2Acc[i] = minNonzero / div
			}
		}
	}
	return prof
}

// vrfsReadVar2File reads a one-value-per-line variance file (recalc file:PATH).
func vrfsReadVar2File(path string) []float64 {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	var out []float64
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		ln := strings.TrimSpace(sc.Text())
		if ln == "" {
			continue
		}
		v, perr := strconv.ParseFloat(strings.Fields(ln)[0], 64)
		if perr != nil {
			continue
		}
		out = append(out, v)
	}
	return out
}

// vrfsRescaleVar2 interpolates an ori-length variance profile onto nbins bins,
// mirroring vrfs.c init_var2's linear intrapolation. When the lengths match it
// returns a copy.
func vrfsRescaleVar2(ori []float64, nbins int) []float64 {
	nori := len(ori)
	if nori == 0 {
		return make([]float64, nbins)
	}
	if nori == nbins {
		out := make([]float64, nbins)
		copy(out, ori)
		return out
	}
	out := make([]float64, nbins)
	out[0] = ori[0]
	out[nbins-1] = ori[nori-1]
	dx := 1.0 / float64(nbins-1)
	DX := 1.0 / float64(nori-1)
	J := 1
	for i := 1; i < nbins-1; i++ {
		x := float64(i) * dx
		for float64(J)*DX < x {
			J++
		}
		X := float64(J-1) * DX
		if x == X {
			out[i] = ori[J-1]
		} else {
			out[i] = ori[J-1] + (ori[J]-ori[J-1])*(x-X)/DX
		}
	}
	return out
}

// scoreVrfsSite computes a site's raw score from its histogram and the variance
// profile, mirroring vrfs.c score_site.
func scoreVrfsSite(s *vrfsSite, var2 []float64, nbins int) float64 {
	maxVal := float64(s.dist[0])
	for i := 0; i < nbins; i++ {
		if v := float64(s.dist[i]); v > maxVal {
			maxVal = v
		}
	}
	if maxVal == 0 {
		maxVal = 1
	}
	score := 0.0
	for i := 1; i < nbins; i++ {
		tmp := float64(s.dist[i]) / maxVal
		score += tmp * tmp / var2[i]
	}
	return 10 * math.Log(1+score)
}

// writeVrfsProfile emits the SITE / MEAN / VAR2 lines, mirroring vrfs.c
// write_batch + print_buffered_sites. Sites sharing a position are buffered and
// re-scored (the lower score is boosted by 75% of its gap to the max; indels
// share the first indel's distribution and score).
func writeVrfsProfile(out io.Writer, sites []*vrfsSite, prof vrfsProfile, cfg vrfsConfig) error {
	w, closeFn, err := vrfsOpenOutput(out, cfg)
	if err != nil {
		return err
	}
	bw := bufio.NewWriter(w)
	nbins := cfg.nbins

	score := make(map[*vrfsSite]float64, len(sites))
	for _, s := range sites {
		score[s] = scoreVrfsSite(s, prof.var2, nbins)
	}

	i := 0
	for i < len(sites) {
		j := i + 1
		for j < len(sites) && sites[j].chrom == sites[i].chrom && sites[j].pos0 == sites[i].pos0 {
			j++
		}
		writeVrfsBufferedSites(bw, sites[i:j], score, nbins)
		i = j
	}

	bw.WriteString("MEAN\t")
	for i := 0; i < nbins; i++ {
		bw.WriteByte(' ')
		bw.WriteString(vrfsFormatE(prof.mean[i]))
	}
	bw.WriteByte('\n')
	bw.WriteString("VAR2\t")
	for i := 0; i < nbins; i++ {
		bw.WriteByte(' ')
		bw.WriteString(vrfsFormatE(prof.var2[i]))
	}
	bw.WriteByte('\n')

	if err := bw.Flush(); err != nil {
		return err
	}
	return closeFn()
}

// writeVrfsBufferedSites emits one group of co-located sites with the score
// rescaling and indel-sharing rules from vrfs.c print_buffered_sites.
func writeVrfsBufferedSites(bw *bufio.Writer, buf []*vrfsSite, score map[*vrfsSite]float64, nbins int) {
	maxScore := 0.0
	for _, s := range buf {
		if score[s] > maxScore {
			maxScore = score[s]
		}
	}
	var fstIndel *vrfsSite
	for _, s := range buf {
		sc := (maxScore-score[s])*0.75 + score[s]
		dist := s.dist
		if s.isIndel() {
			if fstIndel == nil {
				fstIndel = s
			}
			dist = fstIndel.dist
			sc = score[fstIndel]
		}
		bw.WriteString("SITE\t")
		bw.WriteString(s.chrom)
		bw.WriteByte('\t')
		bw.WriteString(strconv.Itoa(s.pos0 + 1))
		bw.WriteByte('\t')
		bw.WriteString(s.ref)
		bw.WriteByte('\t')
		bw.WriteString(s.alt)
		bw.WriteByte('\t')
		bw.WriteString(vrfsFormatE(sc))
		bw.WriteByte('\t')
		for k := 0; k < nbins; k++ {
			if k != 0 {
				bw.WriteByte('-')
			}
			bw.WriteString(strconv.FormatUint(uint64(dist[k]), 10))
		}
		bw.WriteByte('\n')
	}
}

// vrfsFormatE formats a float in C printf "%e" form (six fractional digits,
// two-digit exponent), with NaN rendered as "nan"/"-nan" according to its sign
// bit and +/-Inf as "inf"/"-inf", matching glibc printf (which vrfs.c relies on
// for the empty-profile MEAN line). The empty-profile mean is a positive
// 0.0/0.0 NaN, which glibc prints as "nan".
func vrfsFormatE(v float64) string {
	switch {
	case math.IsNaN(v):
		// glibc prints the NaN's sign bit: "-nan" when set, "nan" when
		// clear. (This used to hard-code "-nan" to the x86-64 result;
		// a positive 0.0/0.0 NaN such as the empty-profile mean is "nan".)
		if math.Signbit(v) {
			return "-nan"
		}
		return "nan"
	case math.IsInf(v, 1):
		return "inf"
	case math.IsInf(v, -1):
		return "-inf"
	}
	return strconv.FormatFloat(v, 'e', 6, 64)
}

// vrfsOpenOutput resolves vrfs's -o/-O output destination. When cfg.output is
// "-" (the default) the profile is written to the host-supplied out writer;
// otherwise the named file is created. -O z (or a .gz output name) wraps the
// stream in a BGZF writer, matching vrfs.c's bgzf_open("wg") path. The returned
// close function flushes/closes any BGZF and file handle (never the host out).
func vrfsOpenOutput(out io.Writer, cfg vrfsConfig) (io.Writer, func() error, error) {
	var sink io.Writer = out
	var file *os.File
	if cfg.output != "" && cfg.output != "-" {
		f, err := os.Create(cfg.output)
		if err != nil {
			return nil, nil, fmt.Errorf("vrfs: create output %q: %w", cfg.output, err)
		}
		file = f
		sink = f
	}
	if !cfg.outputGz {
		return sink, func() error {
			if file != nil {
				return file.Close()
			}
			return nil
		}, nil
	}
	level := cfg.clevel
	var bw *bgzf.Writer
	var err error
	if level < 0 {
		bw = bgzf.NewWriter(sink)
	} else {
		bw, err = bgzf.NewWriterLevel(sink, level)
		if err != nil {
			if file != nil {
				file.Close()
			}
			return nil, nil, fmt.Errorf("vrfs: bgzf writer: %w", err)
		}
	}
	return bw, func() error {
		if cerr := bw.Close(); cerr != nil {
			if file != nil {
				file.Close()
			}
			return cerr
		}
		if file != nil {
			return file.Close()
		}
		return nil
	}, nil
}

// runVrfsMerge merges batch files (-m FILE listing batch files, or -M FILE...)
// into a single profile and writes it, mirroring vrfs.c merge().
func runVrfsMerge(cfg vrfsConfig, out io.Writer) error {
	var fnames []string
	if len(cfg.mergeFiles) > 0 {
		fnames = cfg.mergeFiles
	} else {
		f, err := os.Open(cfg.batchFile)
		if err != nil {
			return fmt.Errorf("vrfs: open %q: %w", cfg.batchFile, err)
		}
		fnames, err = parseAlnList(f)
		f.Close()
		if err != nil {
			return err
		}
	}

	var merged []vrfsSite
	mergedIdx := map[string]int{} // "chrom\x00pos0\x00ref\x00alt" -> index
	nbins := 0
	for _, fn := range fnames {
		sites, fileBins, err := readVrfsBatch(fn)
		if err != nil {
			return err
		}
		if fileBins == 0 {
			continue
		}
		if nbins == 0 {
			nbins = fileBins
		} else if nbins != fileBins {
			return fmt.Errorf("vrfs: different bin size in %q", fn)
		}
		for _, s := range sites {
			key := s.chrom + "\x00" + strconv.Itoa(s.pos0) + "\x00" + s.ref + "\x00" + s.alt
			if mi, ok := mergedIdx[key]; ok {
				for k := 0; k < nbins; k++ {
					merged[mi].dist[k] += s.dist[k]
				}
				merged[mi].nval += s.nval
			} else {
				mergedIdx[key] = len(merged)
				merged = append(merged, s)
			}
		}
	}
	if nbins == 0 {
		return fmt.Errorf("vrfs: failed to merge the files, no usable data found")
	}
	idx := buildVrfsIndex(merged, nbins)
	cfg.nbins = nbins
	prof := computeVrfsProfile(idx.sites, cfg)
	return writeVrfsProfile(out, idx.sites, prof, cfg)
}

// readVrfsBatch parses a batch file (the SITE/MEAN output of a profiling run)
// back into sites, mirroring vrfs.c parse_batch. It returns the sites and the
// number of bins inferred from the histogram column.
func readVrfsBatch(path string) ([]vrfsSite, int, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, 0, fmt.Errorf("vrfs: open batch %q: %w", path, err)
	}
	defer f.Close()
	var sites []vrfsSite
	nbins := 0
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 8<<20)
	for sc.Scan() {
		line := sc.Text()
		if !strings.HasPrefix(line, "SITE\t") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 7 {
			return nil, 0, fmt.Errorf("vrfs: malformed SITE line: %s", line)
		}
		pos1, perr := strconv.Atoi(fields[2])
		if perr != nil {
			return nil, 0, fmt.Errorf("vrfs: could not parse POS: %s", line)
		}
		binStrs := strings.Split(fields[6], "-")
		dist := make([]uint32, len(binStrs))
		nval := 0
		for i, bs := range binStrs {
			n, berr := strconv.ParseUint(bs, 10, 32)
			if berr != nil {
				return nil, 0, fmt.Errorf("vrfs: could not parse DIST: %s", line)
			}
			dist[i] = uint32(n)
			nval += int(n)
		}
		if nbins == 0 {
			nbins = len(dist)
		} else if nbins != len(dist) {
			return nil, 0, fmt.Errorf("vrfs: different number of bins in %q", path)
		}
		sites = append(sites, vrfsSite{
			chrom: fields[1],
			pos0:  pos1 - 1,
			ref:   fields[3],
			alt:   fields[4],
			dist:  dist,
			nval:  nval,
		})
	}
	if err := sc.Err(); err != nil {
		return nil, 0, err
	}
	return sites, nbins, nil
}
