package difffuzz

import (
	"encoding/base64"
	"encoding/json"
	"os"
	"path/filepath"
	"sort"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
)

// Config controls a fuzz run.
type Config struct {
	// Targets is the set of tool/subcommand invocations to fuzz.
	Targets []Target

	// Iterations is the number of fuzzed inputs per target.
	Iterations int

	// Seed seeds the per-target RNG (reproducible).
	Seed int64

	// Timeout bounds each individual binary invocation.
	Timeout time.Duration

	// CacheDir is where our tool binaries are built (reuses the pipeline cache).
	CacheDir string

	// WorkDir is a scratch directory for temp input files and coverage data.
	WorkDir string

	// Manifest provides valid seed-fixture bytes for the mutation strategy.
	Manifest *fixtures.Manifest

	// Coverage enables Go branch/statement coverage capture of our binaries.
	Coverage bool

	// MinimizeSteps bounds delta-debugging predicate evaluations per reproducer.
	MinimizeSteps int

	// MaxReprosPerTarget caps how many minimized reproducers are kept per target
	// (and per class) to keep the report bounded. 0 means a small default.
	MaxReprosPerTarget int

	// Logf, if set, receives progress logging.
	Logf func(format string, args ...any)
}

func (c Config) log(format string, args ...any) {
	if c.Logf != nil {
		c.Logf(format, args...)
	}
}

// Reproducer is a minimized divergence-triggering input recorded in the report.
type Reproducer struct {
	Class     DivergenceClass `json:"class"`
	Origin    Origin          `json:"origin"`
	Detail    string          `json:"detail"`
	InputLen  int             `json:"input_len"`    // minimized length
	OrigLen   int             `json:"original_len"` // pre-minimization length
	Input     string          `json:"input_b64,omitempty"`
	Printable string          `json:"input_text,omitempty"` // shown when printable
}

// TargetReport aggregates one target's results.
type TargetReport struct {
	Name        string                  `json:"name"`
	Tool        string                  `json:"tool"`
	Skipped     bool                    `json:"skipped"`
	SkipReason  string                  `json:"skip_reason,omitempty"`
	Inputs      int                     `json:"inputs"`
	ByClass     map[DivergenceClass]int `json:"by_class"`
	Reproducers []Reproducer            `json:"reproducers,omitempty"`
	Coverage    CoverageResult          `json:"coverage"`
}

// Report is the whole fuzz run's outcome.
type Report struct {
	Seed       int64          `json:"seed"`
	Iterations int            `json:"iterations"`
	Targets    []TargetReport `json:"targets"`
}

// Run executes the differential fuzz across all targets and returns the report.
func Run(cfg Config) (*Report, error) {
	if cfg.Iterations <= 0 {
		cfg.Iterations = 100
	}
	if cfg.Timeout <= 0 {
		cfg.Timeout = 10 * time.Second
	}
	if cfg.MinimizeSteps <= 0 {
		cfg.MinimizeSteps = 500
	}
	if cfg.MaxReprosPerTarget <= 0 {
		cfg.MaxReprosPerTarget = 5
	}
	if cfg.WorkDir == "" {
		cfg.WorkDir = filepath.Join(os.TempDir(), "difffuzz")
	}
	if err := os.MkdirAll(cfg.WorkDir, 0o755); err != nil {
		return nil, err
	}

	rep := &Report{Seed: cfg.Seed, Iterations: cfg.Iterations}
	for _, t := range cfg.Targets {
		tr := runTarget(cfg, t)
		rep.Targets = append(rep.Targets, tr)
	}
	return rep, nil
}

// runTarget fuzzes a single target.
func runTarget(cfg Config, t Target) TargetReport {
	tr := TargetReport{Name: t.Name, Tool: t.Tool, ByClass: map[DivergenceClass]int{}}

	rt, err := t.resolve(cfg.CacheDir)
	if err != nil {
		tr.Skipped = true
		tr.SkipReason = err.Error()
		cfg.log("target %s: SKIP (%v)", t.Name, err)
		return tr
	}

	// Coverage instrumentation of our binary (graceful no-op on failure).
	cb := newCoverageBuilder(cfg.Coverage, t.Tool, rt.ourBin, cfg.WorkDir)
	rt.ourBin = cb.bin
	ourEnv := cb.env()

	// Seed fixture bytes for the mutation strategy.
	var seedBytes []byte
	if t.SeedFixture != "" && cfg.Manifest != nil {
		if p := cfg.Manifest.Path(t.SeedFixture); p != "" {
			if b, rerr := os.ReadFile(p); rerr == nil {
				seedBytes = b
			}
		}
	}
	mut := NewMutator(t.Format, seedBytes, cfg.Seed+hashName(t.Name))

	// Per-target scratch dir for temp input files.
	dir, err := os.MkdirTemp(cfg.WorkDir, "t-"+sanitize(t.Name)+"-")
	if err != nil {
		tr.Skipped = true
		tr.SkipReason = "mkdtemp: " + err.Error()
		return tr
	}
	defer os.RemoveAll(dir)

	// Track which classes we have already kept a (minimized) reproducer for, so
	// the report holds a diverse, bounded set rather than N copies of one bug.
	keptByClass := map[DivergenceClass]int{}

	for i := 0; i < cfg.Iterations; i++ {
		in := mut.Next(i)
		_, _, class, detail := rt.execute(in.Data, dir, cfg.Timeout, ourEnv)
		if !IsDivergence(class) {
			continue
		}
		tr.ByClass[class]++

		if keptByClass[class] >= cfg.MaxReprosPerTarget {
			continue
		}
		// Minimize while the SAME class persists.
		pred := func(cand []byte) bool {
			_, _, c2, _ := rt.execute(cand, dir, cfg.Timeout, ourEnv)
			return c2 == class
		}
		minimized := Minimize(in.Data, pred, cfg.MinimizeSteps)
		tr.Reproducers = append(tr.Reproducers, makeReproducer(class, in.Origin, detail, minimized, len(in.Data)))
		keptByClass[class]++
		cfg.log("target %s: %s found (origin=%s), minimized %d->%d bytes",
			t.Name, class, in.Origin, len(in.Data), len(minimized))
	}

	tr.Inputs = cfg.Iterations
	tr.Coverage = cb.result()
	// Stable order for deterministic reports.
	sort.SliceStable(tr.Reproducers, func(a, b int) bool {
		return tr.Reproducers[a].Class < tr.Reproducers[b].Class
	})
	return tr
}

// makeReproducer builds a Reproducer, recording the input both base64 (always
// round-trippable) and as text when it is printable enough to read.
func makeReproducer(class DivergenceClass, origin Origin, detail string, minimized []byte, origLen int) Reproducer {
	r := Reproducer{
		Class:    class,
		Origin:   origin,
		Detail:   detail,
		InputLen: len(minimized),
		OrigLen:  origLen,
		Input:    base64Encode(minimized),
	}
	if isPrintable(minimized) {
		r.Printable = string(minimized)
	}
	return r
}

// WriteJSON writes the report as indented JSON to path.
func (r *Report) WriteJSON(path string) error {
	b, err := json.MarshalIndent(r, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0o644)
}

// HasParityDivergence reports whether any target found a STRICT parity
// divergence (i.e. a divergence other than both-crashed, which is a shared
// rejection of garbage rather than a drop-in mismatch).
func (r *Report) HasParityDivergence() bool {
	for _, t := range r.Targets {
		for c, n := range t.ByClass {
			if n > 0 && c != ClassNone && c != ClassBothCrashed {
				return true
			}
		}
	}
	return false
}

// hashName turns a target name into a small deterministic offset so each target
// gets a distinct (but reproducible) RNG stream.
func hashName(s string) int64 {
	var h int64
	for _, c := range s {
		h = h*131 + int64(c)
	}
	if h < 0 {
		h = -h
	}
	return h % 100000
}

func sanitize(s string) string {
	out := make([]byte, 0, len(s))
	for i := 0; i < len(s); i++ {
		c := s[i]
		if (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') {
			out = append(out, c)
		} else {
			out = append(out, '_')
		}
	}
	return string(out)
}

func isPrintable(b []byte) bool {
	if len(b) > 4096 {
		return false
	}
	for _, c := range b {
		if c == '\n' || c == '\t' || c == '\r' {
			continue
		}
		if c < 0x20 || c > 0x7e {
			return false
		}
	}
	return true
}

// base64Encode encodes minimized reproducer bytes for the JSON report.
func base64Encode(b []byte) string { return base64.StdEncoding.EncodeToString(b) }
