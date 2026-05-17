// Package prinseq: graph-data (.gd) collection and emission.
//
// This file ports the `--graph_data` path of upstream `prinseq-lite.pl`
// (vendored at `reference_code/prinseq/prinseq-lite.pl`, version 0.20.4).
// The upstream emitter assembles a single-line JSON-shaped string and
// writes it to a `.gd` file, prefixed by two `#`-comment lines. The
// per-read stat collectors are spread across:
//
//   - getSeqStats          (prinseq-lite.pl:4564-4744)
//   - getQualStats         (prinseq-lite.pl:4755-4833)
//   - generateStatsType    (prinseq-lite.pl:4087-4214)
//   - dinucOdds            (prinseq-lite.pl:3977-4023)
//   - checkForDupl         (prinseq-lite.pl:4217-4400)
//   - getTagFrequency      (prinseq-lite.pl:4403-4529)
//   - getBinVal            (prinseq-lite.pl:4835-4849)
//   - the final JSON-string emit block (prinseq-lite.pl:2050-2287)
//
// INTENTIONAL DIVERGENCE FROM UPSTREAM: Perl 5.18+ randomises hash
// iteration on every interpreter start. Upstream's emitter therefore
// produces a key-permuted file on every run — byte parity with
// upstream is mathematically impossible. We sort every map key
// lexicographically before emitting. This is a deliberate
// improvement, not a bug — see `docs/UPSTREAM_BUGS.md > prinseq` and
// the bottom-line note in `docs/PARITY_ROADMAP.md > prinseq-lite`.
// Validation is performed via a JSON-normalised semantic diff against
// the upstream `.gd` (see `tools/PARITY_VALIDATION.md`).

package prinseq

import (
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"
)

// Upstream constants (prinseq-lite.pl lines 38-79).
const (
	gdWindowSize             = 64
	gdWindowStep             = 32
	gdTagLength              = 20
	gdMidCheckLength         = 15
	gdGraphDataSeqMaxLength  = 1000
	gdReverseLogBase         = 62 // ONEOVERLOG62 base
	gdComplexityDustScaling  = 31 // upstream scales dust to 100/31
	gdComplexityEntropyScale = 100
)

// gdMIDS holds the 26 MID/barcode sequences upstream scans for in
// getSeqStats (prinseq-lite.pl:52-77). Order does not affect output
// because MID detection picks the first hit per read.
var gdMIDS = []string{
	"ACGAGTGCGT", "ACGCTCGACA", "AGACGCACTC", "AGCACTGTAG",
	"ATCAGACACG", "ATATCGCGAG", "CGTGTCTCTA", "CTCGCGTGTC",
	"TAGTATCAGC", "TCTCTATGCG", "TGATACGTCT", "TACTGAGCTA",
	"CATAGTAGTG", "CGAGAGATAC", "ACACGACGACT", "ACACGTAGTAT",
	"ACACTACTCGT", "ACGACACGTAT", "ACGAGTAGACT", "ACGCGTCTAGT",
	"ACGTACACACT", "ACGTACTGTGT", "ACGTAGATCGT", "ACTACGTCTCT",
	"ACTATACGAGT", "ACTCGCGTCGT",
}

// gdDinucKeys is the canonical 16-entry dinucleotide alphabet
// (prinseq-lite.pl:79). dinucOdds initialises this with zero counts
// before splitting on runs of N.
var gdDinucKeys = []string{
	"AA", "AC", "AG", "AT",
	"CA", "CC", "CG", "CT",
	"GA", "GC", "GG", "GT",
	"TA", "TC", "TG", "TT",
}

// GraphDataOptions controls graph-data collection. The flag names
// match the upstream CLI knobs that gate which sub-tables get
// emitted. `Enabled` reflects `--graph_data`; the remaining bools
// mirror `--graph_stats` (upstream's `%graphstats`, line 999).
type GraphDataOptions struct {
	// Enabled toggles the entire collector. When false, none of the
	// graph-data paths run (zero-cost).
	Enabled bool

	// Stats selectors. Each bool maps to an upstream `graph_stats`
	// shorthand (line 80, %GRAPH_OPTIONS). When --graph_stats is
	// not supplied upstream sets all of these to 1; we follow suit
	// in DefaultGraphDataOptions().
	LD bool // length distribution
	GC bool // GC distribution
	QD bool // quality distribution
	NS bool // N-content distribution
	PT bool // poly-A/T tail counts
	TS bool // tag-sequence / base-frequency end maps
	AQ bool // (upstream defines `aq` in %GRAPH_OPTIONS but never
	//       references it; we accept the flag for parity but it
	//       has no effect — see UPSTREAM_BUGS.md)
	DE bool // exact duplicates
	DA bool // all duplicate flavours
	SC bool // sequence complexity (dust + entropy)
	DN bool // dinucleotide odds

	// Filename1 / Filename2 override the embedded filename1 /
	// filename2 values in the emitted JSON (upstream looks at
	// $params{filename1} / $params{filename2}, line 2198). Both
	// strings are emitted unchanged when supplied; when empty, the
	// emitter writes the hex-encoded basename of the input path
	// (upstream `convertStringToInt` + `getFileName`, lines
	// 4851-4861).
	Filename1 string
	Filename2 string

	// Phred64 selects offset 64 for ASCII-to-int qual decoding
	// (upstream `--phred64`, lines 230-232). Sanger / Phred+33 is
	// the default.
	Phred64 bool

	// QualNoscale disables the relative (100-bin) quality-by-
	// position table. Upstream's $scale defaults to 1; setting
	// --qual_noscale flips it to 0 (lines 989-993). The absolute
	// `quala` table is always populated regardless.
	QualNoscale bool

	// ExactOnly mirrors upstream's $exactonly switch (referenced
	// at line 2198 only; it controls which duplicate variant the
	// downstream `prinseq-graphs.pl` chooses to draw). Defaults to
	// 0, written verbatim into the JSON.
	ExactOnly int

	// IsFasta makes the emitted "format1" field equal "fasta"
	// instead of the default "fastq". Upstream picks "fasta" when
	// $params{fasta} is defined (line 2198).
	IsFasta bool
}

// DefaultGraphDataOptions returns the upstream default selection,
// where every stat in %GRAPH_OPTIONS is enabled (line 998-999).
func DefaultGraphDataOptions() GraphDataOptions {
	return GraphDataOptions{
		Enabled: true,
		LD:      true,
		GC:      true,
		QD:      true,
		NS:      true,
		PT:      true,
		TS:      true,
		AQ:      true,
		DE:      true,
		DA:      true,
		SC:      true,
		DN:      true,
	}
}

// ParseGraphStatsCSV applies the upstream `--graph_stats` syntax: a
// comma-separated list of two-letter codes (e.g. `gc,qd,ns`). When
// non-empty, only the listed codes are enabled. An unknown code
// returns a non-nil error mirroring upstream's `printError`
// behaviour (line 1013).
func ParseGraphStatsCSV(csv string, opts *GraphDataOptions) error {
	if csv == "" {
		return nil
	}
	// Reset all to false (upstream lines 1004-1007).
	opts.LD = false
	opts.GC = false
	opts.QD = false
	opts.NS = false
	opts.PT = false
	opts.TS = false
	opts.AQ = false
	opts.DE = false
	opts.DA = false
	opts.SC = false
	opts.DN = false
	for _, code := range strings.Split(csv, ",") {
		code = strings.TrimSpace(strings.ToLower(code))
		switch code {
		case "ld":
			opts.LD = true
		case "gc":
			opts.GC = true
		case "qd":
			opts.QD = true
		case "ns":
			opts.NS = true
		case "pt":
			opts.PT = true
		case "ts":
			opts.TS = true
		case "aq":
			opts.AQ = true
		case "de":
			opts.DE = true
		case "da":
			opts.DA = true
		case "sc":
			opts.SC = true
		case "dn":
			opts.DN = true
		default:
			return fmt.Errorf(`unknown option "%s" for -graph_stats`, code)
		}
	}
	return nil
}

// GraphData is the in-memory representation of upstream's
// `%graphdata` hash. Fields use exported names so a future
// consumer (e.g. an HTML report) can introspect them, but the
// emitter (EmitGD) is the canonical user.
type GraphData struct {
	// Per-file totals (lines 2197-2198).
	NumSeqs   int
	NumBases  int
	NumSeqs2  int
	NumBases2 int
	Pairs     int

	// Sequence-length max across all reads, used for binval and
	// the relative-quality stretch / shrink path.
	MaxLength int

	// `counts` sub-tables (length, gc, ns, tail5, tail3). The
	// upstream layout is `counts -> kind -> value -> count`.
	Counts  map[string]map[int]int
	Counts2 map[string]map[int]int // populated only when paired-end

	// `quals` (relative, 100 bins) and `quala` (absolute, per-pos)
	// quality-by-position tables: pos -> value -> count.
	Quals  map[int]map[int]int
	Quala  map[int]map[int]int
	Quals2 map[int]map[int]int
	Quala2 map[int]map[int]int

	// QualsMean: per-read mean qual histogram (line 4773).
	QualsMean  map[int]int
	QualsMean2 map[int]int

	// `freqs` (TS): site (5 / 3) -> position (0..TAG_LENGTH-1) ->
	// base (A/C/G/T/N) -> count. Upstream rescales to int percent
	// in lines 2107-2117.
	Freqs  map[int]map[int]map[byte]int
	Freqs2 map[int]map[int]map[byte]int

	// `kmers` (TS): site -> 5-mer -> count, used to drive the
	// getTagFrequency reduction at lines 2132-2169.
	Kmers  map[int]map[string]int
	Kmers2 map[int]map[string]int

	// `mids`: 5'-tag MID-barcode hits (line 4636).
	Mids map[string]int

	// Sequence complexity scalars (dust + entropy), each
	// stored as the 0..100-scaled mean across windows.
	ComplDust    map[int]int // value -> count
	ComplEntropy map[int]int

	// ComplVals captures min/max sequences for dust/entropy (lines
	// 4700-4718). The score is integer; the sequence is preserved
	// up to GRAPH_DATA_SEQ_MAX_LENGTH chars (then "..." suffix).
	ComplVals map[string]map[string]any // "dust"/"entropy" -> {minval,maxval,minseq,maxseq}

	// DinucOdds accumulates per-read sums; the final emit divides
	// by (numseqs+numseqs2) at line 2175.
	DinucOdds map[string]float64

	// All-sequences buffer for checkForDupl. Each entry is
	// {uppercaseSeq, originalIndex, length}. Stored only when DE
	// or DA is enabled.
	allSeqs []dupSeq

	opts GraphDataOptions
}

// dupSeq is the upstream `[$seq, $idx, $length]` triplet used by
// checkForDupl (line 4218 comment). We keep the upper-case copy.
type dupSeq struct {
	seq string
	idx int
	ln  int
}

// NewGraphData allocates an empty GraphData primed for collection.
func NewGraphData(opts GraphDataOptions) *GraphData {
	return &GraphData{
		Counts:       map[string]map[int]int{},
		Counts2:      map[string]map[int]int{},
		Quals:        map[int]map[int]int{},
		Quala:        map[int]map[int]int{},
		Quals2:       map[int]map[int]int{},
		Quala2:       map[int]map[int]int{},
		QualsMean:    map[int]int{},
		QualsMean2:   map[int]int{},
		Freqs:        map[int]map[int]map[byte]int{},
		Freqs2:       map[int]map[int]map[byte]int{},
		Kmers:        map[int]map[string]int{},
		Kmers2:       map[int]map[string]int{},
		Mids:         map[string]int{},
		ComplDust:    map[int]int{},
		ComplEntropy: map[int]int{},
		ComplVals:    map[string]map[string]any{},
		DinucOdds:    map[string]float64{},
		opts:         opts,
	}
}

// AddSeq updates the per-sequence graph-data tables, matching the
// upstream `getSeqStats` subroutine (lines 4564-4744). `seq` MUST be
// upper-case; callers should upper-case before invocation (upstream
// does this once before dispatch).
func (g *GraphData) AddSeq(seq string) {
	length := len(seq)
	if length == 0 {
		return
	}
	if length > g.MaxLength {
		g.MaxLength = length
	}
	g.NumSeqs++
	g.NumBases += length

	bylength := 100.0 / float64(length)

	// Tail 5'/3' detection (PT). Match upstream's "min 5 char repeats"
	// at lines 4591-4612.
	begin, end := 0, 0
	var str5, str3 string
	if length >= 5 {
		str5 = seq[:5]
		str3 = seq[length-5:]
	}
	if g.opts.PT && length >= 5 {
		if str5 == "AAAAA" || str5 == "TTTTT" {
			c := str5[0]
			begin = 5
			for i := 5; i < length; i++ {
				if seq[i] != c && seq[i] != 'N' {
					break
				}
				begin++
			}
		}
		if str3 == "AAAAA" || str3 == "TTTTT" {
			c := str3[0]
			end = 5
			for i := length - 6; i >= 0; i-- {
				if seq[i] != c && seq[i] != 'N' {
					break
				}
				end++
			}
		}
	}

	// GC content. Upstream uses `tr/GC//` and scales to integer
	// percent via int($gc*$bylength).
	var gcInt int
	if g.opts.GC {
		gc := 0
		for i := 0; i < length; i++ {
			if seq[i] == 'G' || seq[i] == 'C' {
				gc++
			}
		}
		gcInt = int(float64(gc) * bylength)
	}

	// N content. Upstream bumps any sub-percent value to 1 (line 4587).
	var nsInt int
	if g.opts.NS {
		ns := 0
		for i := 0; i < length; i++ {
			if seq[i] == 'N' {
				ns++
			}
		}
		scaled := float64(ns) * bylength
		if ns > 0 && scaled < 1 {
			nsInt = 1
		} else {
			nsInt = int(scaled)
		}
	}

	// Base frequencies + 5-mer kmers + MID hits (TS).
	if g.opts.TS {
		if length >= gdTagLength {
			if g.Freqs[5] == nil {
				g.Freqs[5] = map[int]map[byte]int{}
			}
			if g.Freqs[3] == nil {
				g.Freqs[3] = map[int]map[byte]int{}
			}
			for i := 0; i < gdTagLength; i++ {
				bs5 := seq[i]
				bs3 := seq[length-gdTagLength+i]
				if g.Freqs[5][i] == nil {
					g.Freqs[5][i] = map[byte]int{}
				}
				if g.Freqs[3][i] == nil {
					g.Freqs[3][i] = map[byte]int{}
				}
				g.Freqs[5][i][bs5]++
				g.Freqs[3][i][bs3]++
			}
		}
		if length >= 5 {
			if !(begin > 0) && str5 != "CCCCC" && str5 != "GGGGG" && str5 != "NNNNN" {
				if g.Kmers[5] == nil {
					g.Kmers[5] = map[string]int{}
				}
				g.Kmers[5][str5]++
			}
			if !(end > 0) && str3 != "CCCCC" && str3 != "GGGGG" && str3 != "NNNNN" {
				if g.Kmers[3] == nil {
					g.Kmers[3] = map[string]int{}
				}
				g.Kmers[3][str3]++
			}
		}
		if length >= gdMidCheckLength {
			head := seq[:gdMidCheckLength]
			for _, mid := range gdMIDS {
				if strings.Contains(head, mid) {
					g.Mids[mid]++
					break
				}
			}
		}
	}

	// Sequence complexity (SC): dust + entropy across sliding windows.
	if g.opts.SC {
		dustMean, entropyMean := complexityScores(seq)
		dustVal := int(dustMean * 100.0 / float64(gdComplexityDustScaling))
		entropyVal := int(entropyMean * float64(gdComplexityEntropyScale))
		g.ComplDust[dustVal]++
		g.ComplEntropy[entropyVal]++
		g.updateComplVals("dust", dustVal, seq)
		g.updateComplVals("entropy", entropyVal, seq)
	}

	// Dinucleotide odds ratio (DN).
	if g.opts.DN {
		dinucOdds(seq, g.DinucOdds)
	}

	// Stash all-uppercase seqs for the later checkForDupl call.
	if g.opts.DE || g.opts.DA {
		g.allSeqs = append(g.allSeqs, dupSeq{seq: seq, idx: len(g.allSeqs), ln: length})
	}

	// counts table.
	if g.opts.LD {
		g.ensureCount("length")[length]++
	}
	if begin > 0 {
		g.ensureCount("tail5")[begin]++
	}
	if end > 0 {
		g.ensureCount("tail3")[end]++
	}
	if g.opts.GC {
		g.ensureCount("gc")[gcInt]++
	}
	if nsInt > 0 {
		g.ensureCount("ns")[nsInt]++
	}
}

func (g *GraphData) ensureCount(kind string) map[int]int {
	m, ok := g.Counts[kind]
	if !ok {
		m = map[int]int{}
		g.Counts[kind] = m
	}
	return m
}

func (g *GraphData) updateComplVals(kind string, val int, seq string) {
	store, ok := g.ComplVals[kind]
	if !ok {
		store = map[string]any{}
		g.ComplVals[kind] = store
	}
	if v, ok := store["minval"]; !ok || val < v.(int) {
		store["minval"] = val
		store["minseq"] = clipSeqForGraph(seq)
	}
	if v, ok := store["maxval"]; !ok || val > v.(int) {
		store["maxval"] = val
		store["maxseq"] = clipSeqForGraph(seq)
	}
}

func clipSeqForGraph(seq string) string {
	if len(seq) > gdGraphDataSeqMaxLength {
		return seq[:gdGraphDataSeqMaxLength] + "..."
	}
	return seq
}

// AddQual records the per-position quality table for one read,
// porting upstream's `getQualStats` (lines 4755-4833). `quals` must
// already be in integer-Phred form; the caller is responsible for
// the Phred+33 vs Phred+64 decoding (matching upstream's branching
// at lines 4763-4770).
func (g *GraphData) AddQual(quals []int) {
	if !g.opts.QD {
		return
	}
	length := len(quals)
	if length == 0 {
		return
	}
	// Per-read mean (line 4773).
	sum := 0
	for _, v := range quals {
		sum += v
	}
	mean := sum / length // upstream uses getArrayMean -> int()
	g.QualsMean[mean]++

	// Relative table (line 4777-4828). Only populated when
	// !qual_noscale.
	if !g.opts.QualNoscale {
		switch {
		case length == 100:
			for i := 0; i < 100; i++ {
				g.ensureRelQual(i)[quals[i]]++
			}
		case length < 100:
			factor := 100.0 / float64(length)
			for i := 0; i < length; i++ {
				v := quals[i]
				jStart := int(float64(i) * factor)
				jEnd := int(float64(i+1)*factor) - 1
				for j := jStart; j <= jEnd; j++ {
					g.ensureRelQual(j)[v]++
				}
			}
		default: // length > 100
			factor := float64(length) / 100.0
			for i := 0; i < 100; i++ {
				tmp, count := 0, 0
				jStart := int(float64(i) * factor)
				jEnd := int(float64(i+1)*factor) - 1
				for j := jStart; j <= jEnd; j++ {
					tmp += quals[j]
					count++
				}
				if count > 0 {
					g.ensureRelQual(i)[tmp/count]++
				}
			}
		}
	}

	// Absolute (per-base) table — always emitted (line 4830-4832).
	for i := 0; i < length; i++ {
		if g.Quala[i] == nil {
			g.Quala[i] = map[int]int{}
		}
		g.Quala[i][quals[i]]++
	}
}

func (g *GraphData) ensureRelQual(pos int) map[int]int {
	m, ok := g.Quals[pos]
	if !ok {
		m = map[int]int{}
		g.Quals[pos] = m
	}
	return m
}

// complexityScores computes the dust and entropy means across the
// upstream sliding-window scheme (lines 4646-4695). Returns the
// pre-scaling values; the caller applies the *100/31 and *100
// scalings.
func complexityScores(seq string) (dustMean, entropyMean float64) {
	length := len(seq)
	var steps, rest int
	if length <= gdWindowSize {
		rest = length
		steps = 0
	} else {
		steps = (length-gdWindowSize)/gdWindowStep + 1
		rest = length - steps*gdWindowStep
		if !(rest > gdWindowStep) {
			rest += gdWindowStep
			steps--
		}
	}
	num := gdWindowSize - 2
	bynum := 1.0 / float64(num)
	num--
	var dust, entropy []float64
	counts := make(map[string]int, 64)
	for i := 0; i < steps; i++ {
		window := seq[i*gdWindowStep : i*gdWindowStep+gdWindowSize]
		for k := range counts {
			delete(counts, k)
		}
		for j := 0; j < gdWindowSize-2; j++ {
			counts[window[j:j+3]]++
		}
		dustScore := 0.0
		entropyVal := 0.0
		for _, c := range counts {
			fc := float64(c)
			dustScore += fc * (fc - 1) * 0.5
			p := fc * bynum
			entropyVal -= p * math.Log(p)
		}
		dust = append(dust, dustScore*bynum)
		entropy = append(entropy, entropyVal/math.Log(float64(gdReverseLogBase)))
	}
	if rest > 5 {
		window := seq[steps*gdWindowStep : steps*gdWindowStep+rest]
		for k := range counts {
			delete(counts, k)
		}
		num = rest - 2
		bynum := 1.0 / float64(num)
		for j := 0; j < num; j++ {
			counts[window[j:j+3]]++
		}
		dustScore := 0.0
		entropyVal := 0.0
		for _, c := range counts {
			fc := float64(c)
			dustScore += fc * (fc - 1) * 0.5
			p := fc * bynum
			entropyVal -= p * math.Log(p)
		}
		// Upstream line 4690-4691: special final-step scaling.
		// dustScore /= (num-1), then * ((WINDOWSIZE-2) / num).
		dust = append(dust, (dustScore/float64(num-1))*(float64(gdWindowSize-2)/float64(num)))
		entropy = append(entropy, entropyVal/math.Log(float64(num)))
	} else {
		// Line 4693-4694: assign maximum dust score, zero entropy.
		dust = append(dust, 31)
		entropy = append(entropy, 0)
	}
	return mean64(dust), mean64(entropy)
}

func mean64(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	s := 0.0
	for _, v := range xs {
		s += v
	}
	return s / float64(len(xs))
}

// dinucOdds accumulates the per-read odds-ratio contribution into
// the shared map. Mirrors upstream lines 3977-4023.
func dinucOdds(seq string, odds map[string]float64) {
	// Split on runs of N. The %DN_DI hash is reset per-read in
	// upstream; we mirror that with a fresh map seeded with the
	// 16 canonical dinucleotide keys.
	di := make(map[string]int, 16)
	for _, k := range gdDinucKeys {
		di[k] = 0
	}
	mono := map[string]int{"AT": 0, "GC": 0}
	for _, sub := range splitNRuns(seq) {
		ltmp := len(sub) - 1
		if ltmp <= 0 {
			continue
		}
		for i := 0; i < len(sub); i++ {
			switch sub[i] {
			case 'A', 'T':
				mono["AT"]++
			case 'G', 'C':
				mono["GC"]++
			}
		}
		for i := 0; i < ltmp; i++ {
			di[sub[i:i+2]]++
		}
	}
	dinum := 0
	for _, v := range di {
		dinum += v
	}
	if dinum == 0 {
		return
	}
	mononum := mono["AT"] + mono["GC"]
	factor := 2.0 * float64(mononum) * float64(mononum) / float64(dinum)
	AT := mono["AT"]
	GC := mono["GC"]
	if AT > 0 {
		AT2 := factor / float64(AT*AT)
		odds["AATT"] += float64(di["AA"]+di["TT"]) * AT2
		odds["AT"] += 2 * float64(di["AT"]) * AT2
		odds["TA"] += 2 * float64(di["TA"]) * AT2
		if GC > 0 {
			ATGC := factor / float64(AT*GC)
			odds["ACGT"] += float64(di["AC"]+di["GT"]) * ATGC
			odds["AGCT"] += float64(di["AG"]+di["CT"]) * ATGC
			odds["CATG"] += float64(di["CA"]+di["TG"]) * ATGC
			odds["GATC"] += float64(di["GA"]+di["TC"]) * ATGC
			GC2 := factor / float64(GC*GC)
			odds["CCGG"] += float64(di["CC"]+di["GG"]) * GC2
			odds["CG"] += 2 * float64(di["CG"]) * GC2
			odds["GC"] += 2 * float64(di["GC"]) * GC2
		}
	} else if GC > 0 {
		GC2 := factor / float64(GC*GC)
		odds["CCGG"] += float64(di["CC"]+di["GG"]) * GC2
		odds["CG"] += 2 * float64(di["CG"]) * GC2
		odds["GC"] += 2 * float64(di["GC"]) * GC2
	}
}

// splitNRuns matches Perl's `split(/N+/, $seq)`. Empty substrings
// before the first N or after the last N are preserved, matching
// Perl's default behaviour for trailing-empty drop (we drop them by
// emitting only non-empty substrings, which is what upstream
// effectively gets after the `length-1 > 0` filter).
func splitNRuns(seq string) []string {
	var out []string
	start := 0
	for i := 0; i < len(seq); i++ {
		if seq[i] == 'N' {
			if start < i {
				out = append(out, seq[start:i])
			}
			// skip the run of N
			for i < len(seq) && seq[i] == 'N' {
				i++
			}
			start = i
			// i++ will happen at loop bottom
			i--
		}
	}
	if start < len(seq) {
		out = append(out, seq[start:])
	}
	return out
}

// statsType is the upstream `generateStatsType` output for a single
// `kind`: a fixed set of summary statistics derived from a
// histogram. Each value is stored as the upstream-formatted scalar
// so we can re-emit it verbatim (mean/std are sprintf("%.2f"); the
// rest are ints, written as plain JSON numbers).
type statsType struct {
	Min     int
	Max     int
	Range   int
	Modeval int
	Mode    int
	Median  int
	P25     int
	P75     int
	Mean    string // formatted "%.2f"
	Std     string // formatted "%.2f"
}

// generateStatsType ports the upstream subroutine of the same name
// (lines 4087-4214). Returns a map kind -> statsType.
func generateStatsType(counts map[int]int) statsType {
	// Sort histogram by value ascending.
	keys := make([]int, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	var (
		st          statsType
		minVal      = -1
		maxVal, num int
		modeval     int
		mode        int
		mean        float64
		std         float64
	)
	type cv struct {
		count, val int
	}
	vals := make([]cv, 0, len(keys))
	for _, x1 := range keys {
		c := counts[x1]
		if minVal == -1 {
			minVal = x1
		}
		if maxVal < x1 {
			maxVal = x1
		}
		if modeval < c {
			modeval = c
			mode = x1
		}
		mean += float64(x1 * c)
		num += c
		vals = append(vals, cv{c, x1})
	}
	if num == 0 {
		return statsType{Mean: "0.00", Std: "0.00"}
	}
	mean /= float64(num)
	for x, c := range counts {
		dx := float64(x) - mean
		std += float64(c) * dx * dx
	}

	var median, p25, p75 float64
	switch {
	case num == 1:
		median = float64(vals[0].val)
		p25 = median
		p75 = median
	case num == 2:
		if vals[0].count == 1 {
			p25 = float64(vals[0].val)
			p75 = float64(vals[1].val)
			median = (float64(vals[0].val) + float64(vals[1].val)) / 2
		} else {
			p25 = float64(vals[0].val)
			p75 = p25
			median = p25
		}
	case num > 2:
		var numq int
		if num%2 == 1 {
			i, j := 0, 0
			limit := (num - 1) / 2
			for i <= limit {
				median = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			numq = (num + 1) / 2
		} else {
			i, j := 0, 0
			var median1, median2 float64
			limit := num/2 - 1
			for i <= limit {
				median1 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			median2 = median1
			for i <= num/2 {
				median2 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			median = (median1 + median2) / 2
			numq = num / 2
		}
		if numq%2 == 1 {
			i, j := 0, 0
			limit := (numq - 1) / 2
			for i <= limit {
				p25 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			p75 = p25
			limit2 := num - (numq-1)/2 - 1
			for i <= limit2 {
				p75 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
		} else {
			i, j := 0, 0
			var p251, p252, p751, p752 float64
			limit := numq/2 - 1
			for i <= limit {
				p251 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			p252 = p251
			for i <= numq/2 {
				p252 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			p751 = p252
			limit2 := num - numq/2 - 1
			for i <= limit2 {
				p751 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			p752 = p751
			for i <= num-numq/2 {
				p752 = float64(vals[j].val)
				i += vals[j].count
				j++
			}
			p25 = (p251 + p252) / 2
			p75 = (p751 + p752) / 2
		}
	}

	st.Min = minVal
	st.Max = maxVal
	st.Range = maxVal - minVal + 1
	st.Modeval = modeval
	st.Mode = mode
	st.Median = int(median)
	st.P25 = int(p25)
	st.P75 = int(p75)
	st.Mean = formatFloat2(mean)
	st.Std = formatFloat2(math.Sqrt(std / float64(num)))
	return st
}

// formatFloat2 mirrors Perl's sprintf("%.2f", ...) behaviour. Go's
// strconv.FormatFloat with prec=2 + 'f' matches by default; we
// avoid -0.00 by clamping near-zero values to plain "0.00".
func formatFloat2(f float64) string {
	if f == 0 || (f > -0.005 && f < 0.005) {
		return "0.00"
	}
	return strconv.FormatFloat(f, 'f', 2, 64)
}

// formatFloat9 mirrors sprintf("%.9f", ...). Used only for
// dinucOdds output (line 2175).
func formatFloat9(f float64) string {
	if f == 0 || (f > -0.0000000005 && f < 0.0000000005) {
		return "0.000000000"
	}
	return strconv.FormatFloat(f, 'f', 9, 64)
}

// getBinVal ports upstream lines 4835-4849 verbatim.
func getBinVal(val int) int {
	switch {
	case val <= 0 || val <= 100:
		return 1
	case val < 10000:
		bin := val / 100
		if val%100 != 0 {
			bin++
		}
		return bin
	case val < 100000:
		return 1000
	default:
		step := 1000000
		var xmax int
		if val%step != 0 {
			xmax = (val/step + 1) * step
		} else {
			xmax = val
		}
		return xmax / 100
	}
}

// checkForDuplResult mirrors upstream's three-value return tuple
// (counts, lens, dups).
type checkForDuplResult struct {
	Counts map[int]map[int]int // precount -> type -> count
	Lens   map[int]map[int]int // length  -> type -> count
}

// checkForDupl ports upstream lines 4217-4400. `types` is the set
// of duplicate flavours to detect: 0=exact, 1=prefix, 2=suffix,
// 3=revcomp-exact, 4=revcomp-prefix/suffix.
func checkForDupl(seqs []dupSeq, types map[int]bool) checkForDuplResult {
	out := checkForDuplResult{
		Counts: map[int]map[int]int{},
		Lens:   map[int]map[int]int{},
	}
	if len(seqs) == 0 {
		return out
	}
	dupls := make(map[int]int, len(seqs))

	bumpCount := func(precount, pretype int) {
		if precount <= 0 {
			return
		}
		m, ok := out.Counts[precount]
		if !ok {
			m = map[int]int{}
			out.Counts[precount] = m
		}
		m[pretype]++
	}
	bumpLens := func(length, t int) {
		m, ok := out.Lens[length]
		if !ok {
			m = map[int]int{}
			out.Lens[length] = m
		}
		m[t]++
	}

	// Phase 1: exact + prefix.
	if types[0] || types[1] || types[2] {
		sorted := append([]dupSeq(nil), seqs...)
		sort.Slice(sorted, func(i, j int) bool { return sorted[i].seq < sorted[j].seq })
		pretype, precount := -1, 0
		for i := 0; i < len(sorted)-1; i++ {
			a, b := sorted[i], sorted[i+1]
			if types[0] && a.ln == b.ln && a.seq == b.seq {
				dupls[a.idx] = 0
				bumpLens(a.ln, 0)
				if pretype == 0 {
					precount++
				} else {
					if pretype == 1 && precount > 0 {
						bumpCount(precount, pretype)
					}
					pretype = 0
					precount = 1
				}
			} else if types[1] && a.ln < b.ln && a.seq == b.seq[:a.ln] {
				dupls[a.idx] = 1
				bumpLens(a.ln, 1)
				if pretype == 1 {
					precount++
				} else {
					if pretype == 0 && precount > 0 {
						bumpCount(precount, pretype)
					}
					pretype = 1
					precount = 1
				}
			} else {
				if precount > 0 {
					bumpCount(precount, pretype)
					precount = 0
				}
				pretype = -1
			}
		}
		if precount > 0 {
			bumpCount(precount, pretype)
		}
	}
	// Phase 2: suffix duplicates (reverse the string, exact-prefix).
	if types[2] {
		rev := make([]dupSeq, 0, len(seqs))
		for _, s := range seqs {
			if _, ok := dupls[s.idx]; ok {
				continue
			}
			rev = append(rev, dupSeq{seq: reverseStr(s.seq), idx: s.idx, ln: s.ln})
		}
		if len(rev) > 1 {
			sort.Slice(rev, func(i, j int) bool { return rev[i].seq < rev[j].seq })
			pretype, precount := -1, 0
			for i := 0; i < len(rev)-1; i++ {
				a, b := rev[i], rev[i+1]
				if a.ln < b.ln && a.seq == b.seq[:a.ln] {
					dupls[a.idx] = 2
					bumpLens(a.ln, 2)
					if pretype == 2 {
						precount++
					} else {
						pretype = 2
						precount = 1
					}
				} else {
					if precount > 0 {
						bumpCount(precount, pretype)
						precount = 0
					}
					pretype = -1
				}
			}
			if precount > 0 {
				bumpCount(precount, pretype)
			}
		}
	}
	// Phase 3: revcomp exact + revcomp prefix.
	if types[3] || types[4] {
		type rcSeq struct {
			seq string
			idx int
			ln  int
			rc  int // 0=original, 1=revcomp
		}
		expanded := make([]rcSeq, 0, len(seqs)*2)
		for _, s := range seqs {
			if _, ok := dupls[s.idx]; ok {
				continue
			}
			expanded = append(expanded, rcSeq{seq: s.seq, idx: s.idx, ln: s.ln, rc: 0})
			expanded = append(expanded, rcSeq{seq: revcompUC(s.seq), idx: s.idx, ln: s.ln, rc: 1})
		}
		if len(expanded) > 1 {
			sort.Slice(expanded, func(i, j int) bool { return expanded[i].seq < expanded[j].seq })
			pretype, precount := -1, 0
			for i := 0; i < len(expanded)-1; i++ {
				a, b := expanded[i], expanded[i+1]
				if a.rc == b.rc || a.idx == b.idx {
					continue
				}
				if _, dup := dupls[a.idx]; dup {
					continue
				}
				if types[3] && a.ln == b.ln && a.seq == b.seq {
					dupls[a.idx] = 3
					bumpLens(a.ln, 3)
					if pretype == 3 {
						precount++
					} else {
						if pretype == 4 && precount > 0 {
							bumpCount(precount, pretype)
						}
						pretype = 3
						precount = 1
					}
				} else if types[4] && a.ln < b.ln && a.seq == b.seq[:a.ln] {
					dupls[a.idx] = 4
					bumpLens(a.ln, 4)
					if pretype == 4 {
						precount++
					} else {
						if pretype == 3 && precount > 0 {
							bumpCount(precount, pretype)
						}
						pretype = 4
						precount = 1
					}
				} else {
					if precount > 0 {
						bumpCount(precount, pretype)
						precount = 0
					}
					pretype = -1
				}
			}
			if precount > 0 {
				bumpCount(precount, pretype)
			}
		}
	}
	return out
}

func reverseStr(s string) string {
	b := []byte(s)
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
	return string(b)
}

func revcompUC(s string) string {
	b := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		c := s[len(s)-1-i]
		switch c {
		case 'A':
			c = 'T'
		case 'T':
			c = 'A'
		case 'G':
			c = 'C'
		case 'C':
			c = 'G'
		}
		b[i] = c
	}
	return string(b)
}

// getTagFrequency ports upstream lines 4403-4529. Returns a map
// site -> sum of retained kmer counts. The input `kmers` is
// mutated in place (entries removed) to match upstream semantics.
func getTagFrequency(kmers map[int]map[string]int, numseqs int) map[int]int {
	type maxInfo struct {
		max int
		ten int
		one int
	}
	percentOne := float64(numseqs) / 100.0
	percentTen := float64(numseqs) / 10.0
	most := map[int]*maxInfo{}
	for sp, table := range kmers {
		mi := &maxInfo{}
		for _, c := range table {
			fc := float64(c)
			if fc >= percentTen {
				mi.ten++
			} else if fc >= percentOne {
				mi.one++
			}
			if mi.max < c {
				mi.max = c
			}
		}
		most[sp] = mi
	}
	// Filter kmers in place.
	for sp, table := range kmers {
		mi := most[sp]
		for k, v := range table {
			fv := float64(v)
			if mi.ten > 0 {
				if fv < percentTen {
					delete(table, k)
				}
			} else if mi.one > 0 {
				if fv < percentOne {
					delete(table, k)
				}
			} else {
				if v != mi.max {
					delete(table, k)
				}
			}
		}
	}
	// Collapse to per-site sums. Upstream walks 5' before 3' but
	// the only observable effect is the alignment pass below; we
	// match by iterating in descending site order.
	keys := make([]int, 0, len(kmers))
	for sp := range kmers {
		keys = append(keys, sp)
	}
	sort.Sort(sort.Reverse(sort.IntSlice(keys)))
	out := map[int]int{}
	for _, sp := range keys {
		table := kmers[sp]
		switch len(table) {
		case 0:
			// nothing to record
		case 1:
			for _, v := range table {
				out[sp] += v
			}
		default:
			out[sp] += tagFrequencyMulti(table, sp)
		}
	}
	return out
}

// tagFrequencyMulti reproduces the alignment-based kmer roll-up
// from upstream lines 4444-4520. We compute the per-kmer pairwise
// shift matrix, walk row 0, and count how many sister kmers are
// reachable via a chain of shifts <= 2; the corresponding kmer
// counts are then summed.
func tagFrequencyMulti(table map[string]int, site int) int {
	// Sort kmers by descending count for deterministic order.
	keys := make([]string, 0, len(table))
	for k := range table {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		if table[keys[i]] != table[keys[j]] {
			return table[keys[i]] > table[keys[j]]
		}
		return keys[i] < keys[j]
	})
	n := len(keys)
	if n < 2 {
		sum := 0
		for _, v := range table {
			sum += v
		}
		return sum
	}
	matrix := make([][]cell2, n)
	for i := 0; i < n-1; i++ {
		matrix[i] = make([]cell2, n-i-1)
		for j := i + 1; j < n; j++ {
			a, _ := align2seqs(keys[j], keys[i])
			matrix[i][j-i-1] = a
		}
	}
	countgood := 0
	for i := 0; i < n-1; i++ {
		// matrix[0][i] is the alignment of keys[i+1] against keys[0].
		if matrix[0][i].set == 0 && i > 0 {
			// Try to derive via chain.
			count := 0
			foundShift := false
			for j := 1; j < n-1; j++ {
				if i-j < 0 {
					break
				}
				count++
				if j < len(matrix) && i-j < len(matrix[j]) && matrix[j][i-j].set != 0 {
					foundShift = true
					break
				}
				_ = foundShift
			}
			if count < n-1 {
				sum := 0
				signSet := false
				signNeg := false
				for j := 0; j <= count && j < len(matrix); j++ {
					if i-1 < 0 || i-1 >= len(matrix[j]) {
						continue
					}
					c := matrix[j][i-1]
					if c.set == 0 {
						continue
					}
					if signSet {
						if signNeg {
							if c.set&1 != 0 && c.v0 < 0 {
								sum += c.v0
							} else if c.set&2 != 0 && c.v1 < 0 {
								sum += c.v1
							}
						} else {
							if c.set&1 != 0 && c.v0 > 0 {
								sum += c.v0
							} else if c.set&2 != 0 && c.v1 > 0 {
								sum += c.v1
							}
						}
					} else if c.set&1 != 0 {
						sum += c.v0
					}
					if c.set&1 != 0 {
						signSet = true
						signNeg = c.v0 < 0
					} else if c.set&2 != 0 {
						signSet = true
						signNeg = c.v1 < 0
					}
				}
				if signSet {
					matrix[0][i] = cell2{v0: sum, set: 1}
				}
			}
		}
		if matrix[0][i].set == 0 {
			break
		}
		countgood++
	}
	_ = site
	if countgood == 0 {
		return table[keys[0]]
	}
	sum := table[keys[0]]
	for i := 0; i < countgood; i++ {
		sum += table[keys[i+1]]
	}
	return sum
}

// align2seqs reproduces upstream lines 4531-4548. Returns the first
// and second shift values (if any) packed into a cell.
func align2seqs(seq1, seq2 string) (a cell2, b struct{}) {
	if len(seq1) < 5 || len(seq2) < 5 {
		return cell2{}, struct{}{}
	}
	// shift right by 1, 2
	if seq1[:4] == seq2[1:5] {
		a.v0 = 1
		a.set |= 1
	} else if seq1[:3] == seq2[2:5] {
		a.v0 = 2
		a.set |= 1
	}
	if seq1[1:5] == seq2[:4] {
		a.v1 = -1
		a.set |= 2
	} else if seq1[2:5] == seq2[:3] {
		a.v1 = -2
		a.set |= 2
	}
	return a, b
}

type cell2 struct {
	v0, v1 int
	set    int
}

// EmitGD writes the upstream .gd file structure to `w`. It runs
// the final-pass transformations (binning, percent rescaling,
// stats aggregation) before serialising, mirroring lines
// 2050-2287 of prinseq-lite.pl.
//
// `header` controls the two leading `#`-prefixed comment lines.
// Pass an empty header to omit them entirely (handy for tests).
func (g *GraphData) EmitGD(w io.Writer, header GDHeader) error {
	body := g.buildJSONBody()
	// Write header.
	if header.Format != "" || header.Command != "" || header.Timestamp != "" {
		if _, err := fmt.Fprintln(w, "#Graph data"); err != nil {
			return err
		}
		if _, err := fmt.Fprintf(w, "#[prinseq-lite-%s] [%s] Command: \"%s\"\n",
			defaultIfEmpty(header.Version, "0.20.4"),
			defaultIfEmpty(header.Timestamp, "<timestamp>"),
			defaultIfEmpty(header.Command, "<argv>"),
		); err != nil {
			return err
		}
	}
	if _, err := io.WriteString(w, body); err != nil {
		return err
	}
	return nil
}

// GDHeader configures the optional comment lines preceding the
// JSON body. Leave any field empty to omit the comments.
type GDHeader struct {
	Version   string // typically "0.20.4"
	Timestamp string // formatted upstream-style: MM/DD/YYYY HH:MM:SS
	Command   string // verbatim command-line echo
	Format    string // reserved (currently unused) — set to force header
}

func defaultIfEmpty(s, def string) string {
	if s == "" {
		return def
	}
	return s
}

// buildJSONBody assembles the upstream-shaped JSON-string. Map
// keys are emitted in sort order (lexicographic by string
// representation) for deterministic output. Numeric values are
// emitted unquoted; sprintf-formatted floats (mean/std/dinucodds)
// and explicit string fields are quoted, matching upstream's
// `$v =~ /^\d+$/` test at lines 2205, 2221.
func (g *GraphData) buildJSONBody() string {
	var b strings.Builder
	maxLength := g.MaxLength
	binval := getBinVal(maxLength)

	// Run final-pass transforms.
	statsCounts := transformCounts(g.Counts, g.opts.NS, g.opts.GC, g.opts.LD)
	gdStats := genStatsMap(statsCounts) // stats per kind

	// Compute quals/qualsbin.
	var qualsStats map[int]statsType
	if g.opts.QD && !g.opts.QualNoscale && len(g.Quals) > 0 {
		qualsStats = genStatsMapInt(g.Quals)
	}
	var qualsbinStats map[int]statsType
	if g.opts.QD && len(g.Quala) > 0 {
		// Bin the absolute table. Upstream computes
		// `int(($pos-1)/$binval)` (line 2065), which puts pos=0
		// into bin 0 even though (-1/binval) would be slightly
		// negative — Perl's int() truncates toward zero, so the
		// quotient is 0. We reproduce that with `(pos-1)/binval`
		// when pos > 0, and explicitly clamp pos=0 to bin 0.
		qualsbin := map[int]map[int]int{}
		for pos, byval := range g.Quala {
			var binPos int
			if pos == 0 {
				binPos = 0
			} else {
				binPos = (pos - 1) / binval
			}
			dst, ok := qualsbin[binPos]
			if !ok {
				dst = map[int]int{}
				qualsbin[binPos] = dst
			}
			for v, c := range byval {
				dst[v] += c
			}
		}
		qualsbinStats = genStatsMapInt(qualsbin)
	}

	// Compute freqs (percent rescaling).
	freqsPct := computeFreqsPct(g.Freqs, g.NumSeqs)

	// Tag-probability.
	var tagprob map[int]int
	var tagmidseq string
	tagmidnum := 0
	if len(g.Kmers) > 0 {
		// Clone kmers (getTagFrequency mutates).
		kClone := cloneKmers(g.Kmers)
		kmersum := getTagFrequency(kClone, g.NumSeqs)
		// MID detection.
		midsum := 0
		var midseqs []string
		midKeys := make([]string, 0, len(g.Mids))
		for k := range g.Mids {
			midKeys = append(midKeys, k)
		}
		sort.Strings(midKeys)
		threshold := g.NumSeqs / 34
		for _, mid := range midKeys {
			c := g.Mids[mid]
			midsum += c
			if c > threshold {
				tagmidnum++
				midseqs = append(midseqs, mid)
			}
		}
		if tagmidnum > 0 {
			tagmidseq = strings.Join(midseqs, ",")
		}
		if midsum > kmersum[5] {
			kmersum[5] = midsum
		}
		tagprob = map[int]int{}
		if g.NumSeqs > 0 {
			for kmer, sum := range kmersum {
				tagprob[kmer] = int(100.0 / float64(g.NumSeqs) * float64(sum))
			}
		}
	}

	// dinucodds (mean per read).
	dinucoddsStr := map[string]string{}
	if g.opts.DN {
		divisor := float64(g.NumSeqs + g.NumSeqs2)
		if divisor > 0 {
			keys := make([]string, 0, len(g.DinucOdds))
			for k := range g.DinucOdds {
				keys = append(keys, k)
			}
			sort.Strings(keys)
			for _, k := range keys {
				dinucoddsStr[k] = formatFloat9(g.DinucOdds[k] / divisor)
			}
		}
	}

	// Duplicate detection.
	var dubscounts map[int]map[int]int
	var dubslength map[int]map[int]int
	if g.opts.DE || g.opts.DA {
		types := map[int]bool{}
		if g.opts.DA {
			types[0] = true
			types[1] = true
			types[2] = true
			types[3] = true
			types[4] = true
		} else {
			types[0] = true
		}
		res := checkForDupl(g.allSeqs, types)
		dubscounts = res.Counts
		dubslength = res.Lens
	}

	// Header scalars.
	b.WriteString(`{"numseqs":`)
	b.WriteString(strconv.Itoa(g.NumSeqs))
	b.WriteString(`,"numbases":`)
	b.WriteString(strconv.Itoa(g.NumBases))
	b.WriteString(`,"pairedend":0,"maxlength":`)
	b.WriteString(strconv.Itoa(maxLength))
	b.WriteString(`,"binval":`)
	b.WriteString(strconv.Itoa(binval))
	b.WriteString(`,"exactonly":`)
	b.WriteString(strconv.Itoa(g.opts.ExactOnly))
	b.WriteString(`,"tagmidnum":`)
	b.WriteString(strconv.Itoa(tagmidnum))
	b.WriteString(`,"scale":`)
	if g.opts.QualNoscale {
		b.WriteString("0")
	} else {
		b.WriteString("1")
	}
	b.WriteString(`,"filename1":`)
	b.WriteString(quoteJSON(filenameOrHex(g.opts.Filename1, "")))
	b.WriteString(`,"format1":`)
	b.WriteString(quoteJSON(g.opts.format1Value()))

	// counts.
	if len(statsCounts) > 0 {
		writeIntIntIntBlock(&b, ",\"counts\":", statsCounts)
	}
	// stats (from counts).
	if len(gdStats) > 0 {
		writeStatsTypeBlock(&b, ",\"stats\":", gdStats)
	}
	// quals (relative).
	if len(qualsStats) > 0 {
		writeStatsTypeIntKeyBlock(&b, ",\"quals\":", qualsStats)
	}
	// qualsbin.
	if len(qualsbinStats) > 0 {
		writeStatsTypeIntKeyBlock(&b, ",\"qualsbin\":", qualsbinStats)
	}
	// complvals.
	if g.opts.SC && len(g.ComplVals) > 0 {
		writeComplValsBlock(&b, ",\"complvals\":", g.ComplVals)
	}
	// dubscounts / dubslength.
	if len(dubscounts) > 0 {
		writeIntIntIntBlock(&b, ",\"dubscounts\":", mapIntKey(dubscounts))
	}
	if len(dubslength) > 0 {
		writeIntIntIntBlock(&b, ",\"dubslength\":", mapIntKey(dubslength))
	}
	// flat int->int maps: qualsmean, tagprob, compldust, complentropy.
	if g.opts.QD && len(g.QualsMean) > 0 {
		writeIntIntFlatBlock(&b, ",\"qualsmean\":", g.QualsMean)
	}
	if len(tagprob) > 0 {
		writeIntIntFlatBlock(&b, ",\"tagprob\":", tagprob)
	}
	if g.opts.SC && len(g.ComplDust) > 0 {
		writeIntIntFlatBlock(&b, ",\"compldust\":", g.ComplDust)
	}
	if g.opts.SC && len(g.ComplEntropy) > 0 {
		writeIntIntFlatBlock(&b, ",\"complentropy\":", g.ComplEntropy)
	}
	if g.opts.DN && len(dinucoddsStr) > 0 {
		writeStringStringFlatBlock(&b, ",\"dinucodds\":", dinucoddsStr)
	}

	// tail flag (line 2230).
	hasTail := 0
	if _, ok := g.Counts["tail5"]; ok {
		hasTail = 1
	}
	if _, ok := g.Counts["tail3"]; ok {
		hasTail = 1
	}
	b.WriteString(`,"tail":`)
	b.WriteString(strconv.Itoa(hasTail))

	// freqs (TS).
	if g.opts.TS && len(freqsPct) > 0 {
		writeFreqsBlock(&b, ",\"freqs\":", freqsPct)
	}

	// tagmidseq (string, only when present).
	if tagmidseq != "" {
		b.WriteString(`,"tagmidseq":`)
		b.WriteString(quoteJSON(tagmidseq))
	}

	b.WriteString("}")
	return b.String()
}

// format1Value returns the upstream-shaped "format1" string for
// the emitter, taking IsFasta into account.
func (o GraphDataOptions) format1Value() string {
	if o.IsFasta {
		return "fasta"
	}
	return "fastq"
}

func filenameOrHex(supplied, fallback string) string {
	if supplied != "" {
		return supplied
	}
	if fallback == "" {
		return "stdin"
	}
	// Upstream uses convertStringToInt: hex-encode each byte.
	var b strings.Builder
	for i := 0; i < len(fallback); i++ {
		fmt.Fprintf(&b, "%x", fallback[i])
	}
	return b.String()
}

func quoteJSON(s string) string {
	// Minimal JSON string escaping: this is the same subset upstream
	// uses (no embedded quotes/backslashes/newlines expected in any
	// of our fields).
	var b strings.Builder
	b.WriteByte('"')
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch c {
		case '"':
			b.WriteString(`\"`)
		case '\\':
			b.WriteString(`\\`)
		case '\n':
			b.WriteString(`\n`)
		case '\t':
			b.WriteString(`\t`)
		case '\r':
			b.WriteString(`\r`)
		default:
			if c < 0x20 {
				fmt.Fprintf(&b, `\u%04x`, c)
			} else {
				b.WriteByte(c)
			}
		}
	}
	b.WriteByte('"')
	return b.String()
}

// transformCounts applies upstream lines 2099-2105 (ensure ns key
// exists when NS is on) and returns the int-keyed counts. It does
// not mutate the receiver.
func transformCounts(in map[string]map[int]int, ns, gc, ld bool) map[string]map[int]int {
	out := map[string]map[int]int{}
	for k, v := range in {
		out[k] = v
	}
	if ns {
		if cur, ok := out["ns"]; !ok || len(cur) == 0 {
			out["ns"] = map[int]int{0: 0}
		}
	}
	return out
}

func cloneKmers(in map[int]map[string]int) map[int]map[string]int {
	out := map[int]map[string]int{}
	for k, v := range in {
		c := map[string]int{}
		for kk, vv := range v {
			c[kk] = vv
		}
		out[k] = c
	}
	return out
}

func computeFreqsPct(in map[int]map[int]map[byte]int, numseqs int) map[int]map[int]map[byte]int {
	out := map[int]map[int]map[byte]int{}
	if numseqs == 0 {
		return out
	}
	for site, posMap := range in {
		dst := map[int]map[byte]int{}
		for pos := 0; pos < gdTagLength; pos++ {
			bs := map[byte]int{
				'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
			}
			if cur, ok := posMap[pos]; ok {
				for _, base := range []byte{'A', 'C', 'G', 'T', 'N'} {
					if c, ok := cur[base]; ok {
						bs[base] = c * 100 / numseqs
					}
				}
			}
			dst[pos] = bs
		}
		out[site] = dst
	}
	return out
}

func genStatsMap(in map[string]map[int]int) map[string]statsType {
	out := map[string]statsType{}
	for k, v := range in {
		out[k] = generateStatsType(v)
	}
	return out
}

func genStatsMapInt(in map[int]map[int]int) map[int]statsType {
	out := map[int]statsType{}
	for k, v := range in {
		out[k] = generateStatsType(v)
	}
	return out
}

func mapIntKey(in map[int]map[int]int) map[string]map[int]int {
	out := map[string]map[int]int{}
	for k, v := range in {
		out[strconv.Itoa(k)] = v
	}
	return out
}

// JSON-writer helpers (deterministic by sort order).

func writeIntIntIntBlock(b *strings.Builder, prefix string, m map[string]map[int]int) {
	b.WriteString(prefix)
	b.WriteByte('{')
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(quoteJSON(k))
		b.WriteByte(':')
		writeIntIntFlatBlock(b, "", m[k])
	}
	b.WriteByte('}')
}

func writeIntIntFlatBlock(b *strings.Builder, prefix string, m map[int]int) {
	b.WriteString(prefix)
	b.WriteByte('{')
	keys := make([]int, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		// Lexicographic string sort to match Perl hash-key iteration shape.
		return strconv.Itoa(keys[i]) < strconv.Itoa(keys[j])
	})
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(quoteJSON(strconv.Itoa(k)))
		b.WriteByte(':')
		b.WriteString(strconv.Itoa(m[k]))
	}
	b.WriteByte('}')
}

func writeStringStringFlatBlock(b *strings.Builder, prefix string, m map[string]string) {
	b.WriteString(prefix)
	b.WriteByte('{')
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(quoteJSON(k))
		b.WriteByte(':')
		b.WriteString(quoteJSON(m[k]))
	}
	b.WriteByte('}')
}

func writeStatsTypeBlock(b *strings.Builder, prefix string, m map[string]statsType) {
	b.WriteString(prefix)
	b.WriteByte('{')
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(quoteJSON(k))
		b.WriteByte(':')
		writeStatsTypeObj(b, m[k])
	}
	b.WriteByte('}')
}

func writeStatsTypeIntKeyBlock(b *strings.Builder, prefix string, m map[int]statsType) {
	b.WriteString(prefix)
	b.WriteByte('{')
	keys := make([]int, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		return strconv.Itoa(keys[i]) < strconv.Itoa(keys[j])
	})
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(quoteJSON(strconv.Itoa(k)))
		b.WriteByte(':')
		writeStatsTypeObj(b, m[k])
	}
	b.WriteByte('}')
}

func writeStatsTypeObj(b *strings.Builder, s statsType) {
	// Keys emitted in alphabetic order: max, mean, median, min,
	// mode, modeval, p25, p75, range, std.
	b.WriteString(`{"max":`)
	b.WriteString(strconv.Itoa(s.Max))
	b.WriteString(`,"mean":`)
	b.WriteString(quoteJSON(s.Mean))
	b.WriteString(`,"median":`)
	b.WriteString(strconv.Itoa(s.Median))
	b.WriteString(`,"min":`)
	b.WriteString(strconv.Itoa(s.Min))
	b.WriteString(`,"mode":`)
	b.WriteString(strconv.Itoa(s.Mode))
	b.WriteString(`,"modeval":`)
	b.WriteString(strconv.Itoa(s.Modeval))
	b.WriteString(`,"p25":`)
	b.WriteString(strconv.Itoa(s.P25))
	b.WriteString(`,"p75":`)
	b.WriteString(strconv.Itoa(s.P75))
	b.WriteString(`,"range":`)
	b.WriteString(strconv.Itoa(s.Range))
	b.WriteString(`,"std":`)
	b.WriteString(quoteJSON(s.Std))
	b.WriteByte('}')
}

func writeComplValsBlock(b *strings.Builder, prefix string, m map[string]map[string]any) {
	b.WriteString(prefix)
	b.WriteByte('{')
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for i, k := range keys {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(quoteJSON(k))
		b.WriteString(`:{`)
		inner := m[k]
		ikeys := make([]string, 0, len(inner))
		for ik := range inner {
			ikeys = append(ikeys, ik)
		}
		sort.Strings(ikeys)
		for j, ik := range ikeys {
			if j > 0 {
				b.WriteByte(',')
			}
			b.WriteString(quoteJSON(ik))
			b.WriteByte(':')
			switch v := inner[ik].(type) {
			case int:
				b.WriteString(strconv.Itoa(v))
			case string:
				b.WriteString(quoteJSON(v))
			default:
				b.WriteString(quoteJSON(fmt.Sprintf("%v", v)))
			}
		}
		b.WriteByte('}')
	}
	b.WriteByte('}')
}

func writeFreqsBlock(b *strings.Builder, prefix string, m map[int]map[int]map[byte]int) {
	b.WriteString(prefix)
	b.WriteByte('{')
	siteKeys := make([]int, 0, len(m))
	for k := range m {
		siteKeys = append(siteKeys, k)
	}
	sort.Slice(siteKeys, func(i, j int) bool {
		return strconv.Itoa(siteKeys[i]) < strconv.Itoa(siteKeys[j])
	})
	for i, site := range siteKeys {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(quoteJSON(strconv.Itoa(site)))
		b.WriteString(`:{`)
		posMap := m[site]
		posKeys := make([]int, 0, len(posMap))
		for k := range posMap {
			posKeys = append(posKeys, k)
		}
		sort.Slice(posKeys, func(i, j int) bool {
			return strconv.Itoa(posKeys[i]) < strconv.Itoa(posKeys[j])
		})
		for j, pos := range posKeys {
			if j > 0 {
				b.WriteByte(',')
			}
			b.WriteString(quoteJSON(strconv.Itoa(pos)))
			b.WriteString(`:{`)
			bases := posMap[pos]
			bkeys := make([]byte, 0, len(bases))
			for k := range bases {
				bkeys = append(bkeys, k)
			}
			sort.Slice(bkeys, func(i, j int) bool { return bkeys[i] < bkeys[j] })
			for k, base := range bkeys {
				if k > 0 {
					b.WriteByte(',')
				}
				b.WriteString(quoteJSON(string(base)))
				b.WriteByte(':')
				b.WriteString(strconv.Itoa(bases[base]))
			}
			b.WriteByte('}')
		}
		b.WriteByte('}')
	}
	b.WriteByte('}')
}
