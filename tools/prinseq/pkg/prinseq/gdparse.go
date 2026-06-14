package prinseq

// gdparse.go reads back a prinseq-lite `.gd` graph-data file (the
// JSON-shaped single line emitted by EmitGD / upstream
// prinseq-lite.pl) into a typed in-memory model. This is the input
// side of the PNG rendering flow ported from upstream
// `prinseq-graphs.pl` (vendored at reference_code/prinseq/).
//
// The .gd payload is one JSON object on a single line, optionally
// preceded by `#`-comment lines. Numeric histogram keys are JSON
// object keys (hence strings); we decode them back to ints. Mean /
// std / dinucodds values are sprintf-formatted decimal strings.

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// GDStat mirrors a single `stats`/`quals`/`qualsbin` entry produced
// by upstream generateStatsType. Mean and Std arrive as formatted
// strings in the .gd; the rest are integers.
type GDStat struct {
	Min     int
	Max     int
	Range   int
	Modeval int
	Mode    int
	Median  int
	P25     int
	P75     int
	Mean    float64
	Std     float64
}

// GDData is the parsed representation of a `.gd` file. Only the
// fields consumed by the PNG renderer are populated; unknown keys
// are ignored. Histogram maps use int value -> count, matching the
// upstream `counts->kind->value->count` layout after key decoding.
type GDData struct {
	NumSeqs   int
	NumBases  int
	MaxLength int
	BinVal    int
	Scale     int
	Filename1 string
	Format1   string

	// counts sub-tables: length, gc, ns, tail5, tail3.
	Counts map[string]map[int]int
	// stats sub-tables keyed by the same kind names.
	Stats map[string]GDStat

	// Tail flag (1 when poly-A/T tail tables exist).
	Tail int

	// quals: relative (100-bin) per-position quality boxplot stats.
	Quals map[int]GDStat
	// qualsbin: absolute (binned) per-position quality boxplot stats.
	QualsBin map[int]GDStat
	// qualsmean: per-read mean-quality histogram.
	QualsMean map[int]int

	// Sequence-complexity histograms.
	ComplDust    map[int]int
	ComplEntropy map[int]int

	// dinucodds: dinucleotide odds-ratio means (key -> value).
	DinucOdds map[string]float64

	// Duplicate tables: precount/length -> dup-type -> count.
	DubsCounts map[int]map[int]int
	DubsLength map[int]map[int]int
}

// rawGD is the JSON shape used purely for decoding. Histogram keys
// are strings in JSON; we re-key them to ints after unmarshalling.
type rawGD struct {
	NumSeqs   int    `json:"numseqs"`
	NumBases  int    `json:"numbases"`
	MaxLength int    `json:"maxlength"`
	BinVal    int    `json:"binval"`
	Scale     int    `json:"scale"`
	Filename1 string `json:"filename1"`
	Format1   string `json:"format1"`
	Tail      int    `json:"tail"`

	Counts       map[string]map[string]int `json:"counts"`
	Stats        map[string]map[string]any `json:"stats"`
	Quals        map[string]map[string]any `json:"quals"`
	QualsBin     map[string]map[string]any `json:"qualsbin"`
	QualsMean    map[string]int            `json:"qualsmean"`
	ComplDust    map[string]int            `json:"compldust"`
	ComplEntropy map[string]int            `json:"complentropy"`
	DinucOdds    map[string]string         `json:"dinucodds"`
	DubsCounts   map[string]map[string]int `json:"dubscounts"`
	DubsLength   map[string]map[string]int `json:"dubslength"`
}

// ParseGD reads a `.gd` payload from r, skipping leading `#`-comment
// lines, and decodes the JSON body into a *GDData. It returns an
// error if no JSON object is found or the body is malformed.
func ParseGD(r io.Reader) (*GDData, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 64*1024*1024)
	var body string
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		if strings.TrimSpace(line) == "" {
			continue
		}
		// The JSON object lives on a single line; upstream emits one
		// object per file. Use the first non-comment, non-blank line.
		body = line
		break
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	if body == "" {
		return nil, fmt.Errorf("no graph-data JSON object found in input")
	}

	var raw rawGD
	if err := json.Unmarshal([]byte(body), &raw); err != nil {
		return nil, fmt.Errorf("parsing graph-data JSON: %w", err)
	}

	d := &GDData{
		NumSeqs:      raw.NumSeqs,
		NumBases:     raw.NumBases,
		MaxLength:    raw.MaxLength,
		BinVal:       raw.BinVal,
		Scale:        raw.Scale,
		Filename1:    raw.Filename1,
		Format1:      raw.Format1,
		Tail:         raw.Tail,
		Counts:       map[string]map[int]int{},
		Stats:        map[string]GDStat{},
		Quals:        map[int]GDStat{},
		QualsBin:     map[int]GDStat{},
		QualsMean:    rekeyIntInt(raw.QualsMean),
		ComplDust:    rekeyIntInt(raw.ComplDust),
		ComplEntropy: rekeyIntInt(raw.ComplEntropy),
		DinucOdds:    map[string]float64{},
		DubsCounts:   rekeyIntIntInt(raw.DubsCounts),
		DubsLength:   rekeyIntIntInt(raw.DubsLength),
	}

	for kind, hist := range raw.Counts {
		d.Counts[kind] = rekeyIntInt(hist)
	}
	for kind, st := range raw.Stats {
		d.Stats[kind] = decodeStat(st)
	}
	for k, st := range raw.Quals {
		if pos, err := strconv.Atoi(k); err == nil {
			d.Quals[pos] = decodeStat(st)
		}
	}
	for k, st := range raw.QualsBin {
		if pos, err := strconv.Atoi(k); err == nil {
			d.QualsBin[pos] = decodeStat(st)
		}
	}
	for k, v := range raw.DinucOdds {
		f, err := strconv.ParseFloat(v, 64)
		if err != nil {
			return nil, fmt.Errorf("parsing dinucodds value %q: %w", v, err)
		}
		d.DinucOdds[k] = f
	}

	return d, nil
}

// rekeyIntInt converts a string-keyed histogram into an int-keyed
// one, dropping entries whose key is not an integer.
func rekeyIntInt(in map[string]int) map[int]int {
	out := make(map[int]int, len(in))
	for k, v := range in {
		if ik, err := strconv.Atoi(k); err == nil {
			out[ik] = v
		}
	}
	return out
}

// rekeyIntIntInt converts a doubly string-keyed table into an
// int->int->count table.
func rekeyIntIntInt(in map[string]map[string]int) map[int]map[int]int {
	out := make(map[int]map[int]int, len(in))
	for k, inner := range in {
		ik, err := strconv.Atoi(k)
		if err != nil {
			continue
		}
		out[ik] = rekeyIntInt(inner)
	}
	return out
}

// decodeStat converts a JSON-decoded stats map (mixed number / string
// values) into a typed GDStat. Mean and std are formatted decimal
// strings in the .gd; the rest decode as JSON numbers (float64).
func decodeStat(m map[string]any) GDStat {
	var st GDStat
	st.Min = anyInt(m["min"])
	st.Max = anyInt(m["max"])
	st.Range = anyInt(m["range"])
	st.Modeval = anyInt(m["modeval"])
	st.Mode = anyInt(m["mode"])
	st.Median = anyInt(m["median"])
	st.P25 = anyInt(m["p25"])
	st.P75 = anyInt(m["p75"])
	st.Mean = anyFloat(m["mean"])
	st.Std = anyFloat(m["std"])
	return st
}

func anyInt(v any) int {
	switch x := v.(type) {
	case float64:
		return int(x)
	case string:
		f, _ := strconv.ParseFloat(x, 64)
		return int(f)
	default:
		return 0
	}
}

func anyFloat(v any) float64 {
	switch x := v.(type) {
	case float64:
		return x
	case string:
		f, _ := strconv.ParseFloat(x, 64)
		return f
	default:
		return 0
	}
}

// SortedIntKeys returns the integer keys of a histogram in ascending
// order. Helper shared by the renderer.
func SortedIntKeys(m map[int]int) []int {
	keys := make([]int, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	return keys
}
