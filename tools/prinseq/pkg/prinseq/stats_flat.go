package prinseq

// This file ports upstream prinseq-lite.pl's `-stats_*` reporting path
// (the flat "stats_<group>\t<key>\t<value>" TSV written to STDOUT), as
// distinct from the `--graph_data` `.gd` path in graphdata.go. The two
// share several primitives (dinucOdds, getTagFrequency, checkForDupl,
// the length-histogram statistics), but the stats path collects a
// slightly different kmer set and computes N-content and assembly (Nx)
// figures the graph-data path does not surface. The collector below
// mirrors `calcSeqStats` (prinseq-lite.pl:3918-3974) and the summary
// block at lines 1944-2048 exactly.

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"
)

// StatsGroups selects which upstream `-stats_*` groups to emit. It maps
// one-to-one to the upstream bool flags. StatsAll enables every group
// (upstream `-stats_all`, lines 676-685).
type StatsGroups struct {
	Info     bool
	Len      bool
	Dupl     bool
	Tag      bool
	Dinuc    bool
	Ns       bool
	Assembly bool
}

// StatsGroupsAll returns the group selection upstream applies for
// -stats_all (single-end: every group; -stats_dupl is skipped for paired
// input, which this port does not handle here). See prinseq-lite.pl
// lines 676-685.
func StatsGroupsAll() StatsGroups {
	return StatsGroups{Info: true, Len: true, Dupl: true, Tag: true, Dinuc: true, Ns: true, Assembly: true}
}

// Any reports whether at least one group is selected.
func (s StatsGroups) Any() bool {
	return s.Info || s.Len || s.Dupl || s.Tag || s.Dinuc || s.Ns || s.Assembly
}

// statsCollector accumulates the per-read counters used by the flat
// stats path. Field semantics mirror the upstream lexical variables
// referenced by calcSeqStats and the summary block.
type statsCollector struct {
	groups StatsGroups

	numSeqs  int
	numBases int

	lengthCounts map[int]int // stats_len / stats_assembly

	odds map[string]float64 // stats_dinuc

	kmers map[int]map[string]int // stats_tag (5' / 3')
	mids  map[string]int         // stats_tag MID hits

	seqswithn int // stats_ns
	maxn      int
	maxp      int

	allSeqs []dupSeq // stats_dupl
}

// newStatsCollector primes a collector for the requested groups.
func newStatsCollector(groups StatsGroups) *statsCollector {
	return &statsCollector{
		groups:       groups,
		lengthCounts: map[int]int{},
		odds:         map[string]float64{},
		kmers:        map[int]map[string]int{},
		mids:         map[string]int{},
	}
}

// addSeq folds one upper-cased sequence into the collector, mirroring
// calcSeqStats (prinseq-lite.pl:3918-3974) plus the numseqs/numbases
// accumulation the main loop performs before dispatch.
func (c *statsCollector) addSeq(seq string) {
	length := len(seq)
	if length == 0 {
		return
	}
	c.numSeqs++
	c.numBases += length

	if c.groups.Len || c.groups.Assembly {
		c.lengthCounts[length]++
	}

	if c.groups.Dinuc {
		dinucOdds(seq, c.odds)
	}

	if c.groups.Tag {
		if length >= 5 {
			str5 := seq[:5]
			str3 := seq[length-5:]
			if !isHomo5(str5) {
				if c.kmers[5] == nil {
					c.kmers[5] = map[string]int{}
				}
				c.kmers[5][str5]++
			}
			if !isHomo5(str3) {
				if c.kmers[3] == nil {
					c.kmers[3] = map[string]int{}
				}
				c.kmers[3][str3]++
			}
		}
		if length >= gdMidCheckLength {
			head := seq[:gdMidCheckLength]
			for _, mid := range gdMIDS {
				if strings.Contains(head, mid) {
					c.mids[mid]++
					break
				}
			}
		}
	}

	if c.groups.Ns {
		ns := 0
		for i := 0; i < length; i++ {
			if seq[i] == 'N' {
				ns++
			}
		}
		if ns > 0 {
			c.seqswithn++
		}
		if ns > c.maxn {
			c.maxn = ns
		}
		// Upstream: $ns = ($ns > 0 && $ns*bylength < 1 ? 1 : int($ns*bylength)).
		bylength := 100.0 / float64(length)
		var nsp int
		scaled := float64(ns) * bylength
		if ns > 0 && scaled < 1 {
			nsp = 1
		} else {
			nsp = int(scaled)
		}
		if nsp > c.maxp {
			c.maxp = nsp
		}
	}

	if c.groups.Dupl {
		c.allSeqs = append(c.allSeqs, dupSeq{seq: seq, idx: len(c.allSeqs), ln: length})
	}
}

// isHomo5 reports whether a 5-char window is one of the homopolymers
// upstream excludes from the kmer tables (AAAAA/TTTTT/CCCCC/GGGGG/NNNNN;
// prinseq-lite.pl:3939-3944).
func isHomo5(s string) bool {
	switch s {
	case "AAAAA", "TTTTT", "CCCCC", "GGGGG", "NNNNN":
		return true
	}
	return false
}

// statsLine is one emitted "stats_<group>\t<key>\t<value>" row.
type statsLine struct {
	group string
	key   string
	value string
}

// buildLines produces the fully-ordered set of output lines. Upstream
// sorts by group name then key name (both string sorts) before printing
// (prinseq-lite.pl:2043-2046).
func (c *statsCollector) buildLines() []statsLine {
	// Assemble a group -> (key -> value) map, then flatten in sorted order.
	groups := map[string]map[string]string{}
	put := func(g, k, v string) {
		if groups[g] == nil {
			groups[g] = map[string]string{}
		}
		groups[g][k] = v
	}

	if c.groups.Info {
		put("stats_info", "reads", strconv.Itoa(c.numSeqs))
		put("stats_info", "bases", strconv.Itoa(c.numBases))
	}

	if c.groups.Len {
		st := generateStatsLen(c.lengthCounts)
		put("stats_len", "min", strconv.Itoa(st.min))
		put("stats_len", "max", strconv.Itoa(st.max))
		put("stats_len", "range", strconv.Itoa(st.rangeVal))
		put("stats_len", "modeval", strconv.Itoa(st.modeval))
		put("stats_len", "mode", strconv.Itoa(st.mode))
		put("stats_len", "mean", st.mean)
		put("stats_len", "stddev", st.stddev)
		put("stats_len", "median", strconv.Itoa(st.median))
	}

	if c.groups.Dinuc {
		divisor := float64(c.numSeqs)
		for k, v := range c.odds {
			var val string
			if divisor > 0 {
				val = formatFloat9(v / divisor)
			} else {
				val = formatFloat9(0)
			}
			put("stats_dinuc", strings.ToLower(k), val)
		}
	}

	if c.groups.Tag {
		// Clone kmers (getTagFrequency mutates its input).
		kClone := cloneKmers(c.kmers)
		kmersum := getTagFrequency(kClone, c.numSeqs)
		// MID detection (prinseq-lite.pl:1966-1983).
		midsum := 0
		midcount := 0
		var midseqs []string
		midKeys := make([]string, 0, len(c.mids))
		for k := range c.mids {
			midKeys = append(midKeys, k)
		}
		sort.Strings(midKeys)
		threshold := float64(c.numSeqs) / 34.0
		for _, mid := range midKeys {
			cnt := c.mids[mid]
			midsum += cnt
			if float64(cnt) > threshold {
				midcount++
				midseqs = append(midseqs, mid)
			}
		}
		put("stats_tag", "midnum", strconv.Itoa(midcount))
		if midcount > 0 {
			put("stats_tag", "midseq", strings.Join(midseqs, ","))
		}
		if midsum > kmersum[5] {
			kmersum[5] = midsum
		}
		if c.numSeqs > 0 {
			for kmer, sum := range kmersum {
				prob := int(100.0 / float64(c.numSeqs) * float64(sum))
				put("stats_tag", "prob"+strconv.Itoa(kmer), strconv.Itoa(prob))
			}
		}
	}

	// Upstream only autovivifies the stats_ns keys when at least one read
	// contained an N (all three assignments are gated on ns > 0). A file with
	// no ambiguous bases therefore emits no stats_ns lines at all
	// (prinseq-lite.pl:3958-3973).
	if c.groups.Ns && c.seqswithn > 0 {
		put("stats_ns", "seqswithn", strconv.Itoa(c.seqswithn))
		put("stats_ns", "maxn", strconv.Itoa(c.maxn))
		put("stats_ns", "maxp", strconv.Itoa(c.maxp))
	}

	if c.groups.Assembly {
		nx := assemblyNx(c.lengthCounts, c.numBases)
		put("stats_assembly", "N50", nx[50])
		put("stats_assembly", "N75", nx[75])
		put("stats_assembly", "N90", nx[90])
		put("stats_assembly", "N95", nx[95])
	}

	if c.groups.Dupl {
		types := map[int]bool{0: true, 1: true, 2: true, 3: true, 4: true}
		res := checkForDupl(c.allSeqs, types)
		// Aggregate precount -> type counts into per-type totals and
		// maxd, plus the grand total (prinseq-lite.pl:2028-2041).
		typeNames := map[int]string{0: "exact", 1: "5", 2: "3", 3: "exactrevcomp", 4: "revcomp"}
		totals := map[string]int{}
		maxd := map[string]int{}
		for _, name := range typeNames {
			totals[name] = 0
			maxd[name] = 0
		}
		total := 0
		for n, byType := range res.Counts {
			for t, cnt := range byType {
				name := typeNames[t]
				totals[name] += cnt * n
				if n > maxd[name] {
					maxd[name] = n
				}
				total += cnt * n
			}
		}
		for _, name := range typeNames {
			put("stats_dupl", name, strconv.Itoa(totals[name]))
			put("stats_dupl", name+"maxd", strconv.Itoa(maxd[name]))
		}
		// Upstream pre-initialises the per-type counts to 0 but only
		// autovivifies `total` inside the duplicate loop, so a file with no
		// duplicates emits every per-type line but omits `total`
		// (prinseq-lite.pl:2031-2041).
		if total > 0 {
			put("stats_dupl", "total", strconv.Itoa(total))
		}
	}

	// Flatten in sorted group/key order.
	groupKeys := make([]string, 0, len(groups))
	for g := range groups {
		groupKeys = append(groupKeys, g)
	}
	sort.Strings(groupKeys)
	var out []statsLine
	for _, g := range groupKeys {
		keys := make([]string, 0, len(groups[g]))
		for k := range groups[g] {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		for _, k := range keys {
			out = append(out, statsLine{group: g, key: k, value: groups[g][k]})
		}
	}
	return out
}

// statsLenResult holds the generateStats() output for the length
// histogram (prinseq-lite.pl:4026-4085). Unlike generateStatsType it
// has no p25/p75.
type statsLenResult struct {
	min, max, rangeVal, modeval, mode, median int
	mean, stddev                              string
}

// generateStatsLen ports the upstream generateStats subroutine used for
// stats_len (prinseq-lite.pl:4026-4085).
func generateStatsLen(counts map[int]int) statsLenResult {
	var (
		min          = -1
		max, modeval int
		mode, count  int
		mean         float64
		std          float64
	)
	keys := make([]int, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	num := 0
	for _, x := range keys {
		c := counts[x]
		if min == -1 {
			min = x
		} else if min > x {
			min = x
		}
		if max < x {
			max = x
		}
		if modeval < c {
			modeval = c
			mode = x
		}
		mean += float64(x * c)
		count += c
		num += c
	}
	if count == 0 {
		return statsLenResult{min: 0, mean: "0.00", stddev: "0.00"}
	}
	mean /= float64(count)
	for x, c := range counts {
		dx := float64(x) - mean
		std += float64(c) * dx * dx
	}

	// median: upstream expands @vals (sorted) then indexes.
	var median int
	switch {
	case num == 1:
		median = keys[0]
	case num == 2:
		// vals holds two entries; upstream averages them (integer keys
		// so the /2 may be fractional but is assigned to $median and
		// printed verbatim). Reconstruct the two values.
		vals := expandVals(keys, counts, 2)
		median = (vals[0] + vals[1]) / 2
	default:
		if num%2 == 1 {
			median = nthVal(keys, counts, (num-1)/2)
		} else {
			a := nthVal(keys, counts, num/2)
			b := nthVal(keys, counts, num/2-1)
			median = (a + b) / 2
		}
	}

	return statsLenResult{
		min:      min,
		max:      max,
		rangeVal: max - min + 1,
		modeval:  modeval,
		mode:     mode,
		mean:     formatFloat2(mean),
		stddev:   formatFloat2(math.Sqrt(std / float64(count))),
		median:   median,
	}
}

// expandVals rebuilds the first n entries of the sorted value list from
// a histogram (ascending keys, each repeated by its count).
func expandVals(sortedKeys []int, counts map[int]int, n int) []int {
	out := make([]int, 0, n)
	for _, k := range sortedKeys {
		for i := 0; i < counts[k] && len(out) < n; i++ {
			out = append(out, k)
		}
		if len(out) >= n {
			break
		}
	}
	return out
}

// nthVal returns the value at index idx (0-based) in the ascending,
// count-expanded value list without materialising the whole slice.
func nthVal(sortedKeys []int, counts map[int]int, idx int) int {
	acc := 0
	for _, k := range sortedKeys {
		acc += counts[k]
		if idx < acc {
			return k
		}
	}
	if len(sortedKeys) == 0 {
		return 0
	}
	return sortedKeys[len(sortedKeys)-1]
}

// assemblyNx computes the N50/N75/N90/N95 assembly figures from the
// length histogram, mirroring prinseq-lite.pl:1994-2022. Missing values
// are reported as "-".
func assemblyNx(lengthCounts map[int]int, numBases int) map[int]string {
	out := map[int]string{50: "-", 75: "-", 90: "-", 95: "-"}
	// Sort lengths descending.
	lens := make([]int, 0, len(lengthCounts))
	for l := range lengthCounts {
		lens = append(lens, l)
	}
	sort.Sort(sort.Reverse(sort.IntSlice(lens)))
	n50 := float64(numBases) * 0.5
	n75 := float64(numBases) * 0.75
	n90 := float64(numBases) * 0.9
	n95 := float64(numBases) * 0.95
	curlen := 0.0
	var got50, got75, got90, got95 bool
	for _, l := range lens {
		for i := 0; i < lengthCounts[l]; i++ {
			curlen += float64(l)
			if curlen >= n50 && !got50 {
				out[50] = strconv.Itoa(l)
				got50 = true
			} else if curlen >= n75 && !got75 {
				out[75] = strconv.Itoa(l)
				got75 = true
			} else if curlen >= n90 && !got90 {
				out[90] = strconv.Itoa(l)
				got90 = true
			} else if curlen >= n95 && !got95 {
				out[95] = strconv.Itoa(l)
				got95 = true
			}
		}
	}
	return out
}

// CollectFlatStats reads a FASTA/FASTQ stream and returns the ordered
// list of "stats_<group>\t<key>\t<value>" rows for the requested groups,
// mirroring upstream prinseq-lite.pl's `-stats_*` output. The rows are
// already sorted by group then key. `isFastq` selects the parser.
func CollectFlatStats(reader io.Reader, isFastq bool, groups StatsGroups) ([]string, error) {
	c := newStatsCollector(groups)
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	if isFastq {
		if err := collectFlatFastq(scanner, c); err != nil {
			return nil, err
		}
	} else {
		if err := collectFlatFasta(scanner, c); err != nil {
			return nil, err
		}
	}
	lines := c.buildLines()
	out := make([]string, 0, len(lines))
	for _, l := range lines {
		out = append(out, fmt.Sprintf("%s\t%s\t%s", l.group, l.key, l.value))
	}
	return out, nil
}

func collectFlatFasta(scanner *bufio.Scanner, c *statsCollector) error {
	var seq strings.Builder
	first := true
	flush := func() {
		if seq.Len() == 0 {
			return
		}
		c.addSeq(strings.ToUpper(seq.String()))
		seq.Reset()
	}
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if !first {
				flush()
			}
			first = false
			continue
		}
		seq.WriteString(line)
	}
	flush()
	return scanner.Err()
}

func collectFlatFastq(scanner *bufio.Scanner, c *statsCollector) error {
	state := 0
	var seq string
	for scanner.Scan() {
		line := scanner.Text()
		switch state {
		case 0:
			if !strings.HasPrefix(line, "@") {
				return fmt.Errorf("expected '@' at start of FASTQ header, got: %q", line)
			}
			state = 1
		case 1:
			seq = strings.ToUpper(line)
			state = 2
		case 2:
			if !strings.HasPrefix(line, "+") {
				return fmt.Errorf("expected '+' separator, got: %q", line)
			}
			state = 3
		case 3:
			c.addSeq(seq)
			seq = ""
			state = 0
		}
	}
	return scanner.Err()
}
