package realbench

import (
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// bedSubcommand maps each of our per-bed tool names to its upstream `bedtools
// <subcommand>` equivalent. A tool present here with a non-empty subcommand has
// a clean upstream pair; the empty-string entries (or absence) mean ours-only.
var bedSubcommand = map[string]string{
	"bedintersect":   "intersect",
	"bedmerge":       "merge",
	"bedsort":        "sort",
	"bedcoverage":    "coverage",
	"bedgenomecov":   "genomecov",
	"bedclosest":     "closest",
	"bedmap":         "map",
	"bedsubtract":    "subtract",
	"bedwindow":      "window",
	"bedslop":        "slop",
	"bedflank":       "flank",
	"bedshift":       "shift",
	"bedjaccard":     "jaccard",
	"bedfisher":      "fisher",
	"bedcomplement":  "complement",
	"bedmakewindows": "makewindows",
	"bedmulticov":    "multicov",
	"bedmultiinter":  "multiinter",
	"bednuc":         "nuc",
	"bedgetfasta":    "getfasta",
	"bed12tobed6":    "bed12tobed6",
	"bedbamtobed":    "bamtobed",
	"bedtobam":       "bedtobam",
	"bedexpand":      "expand",
	"bedgroupby":     "groupby",
	"bedannotate":    "annotate",
	"bedoverlap":     "overlap",
	"bedpairtobed":   "pairtobed",
	"bedpairtopair":  "pairtopair",
	"bedrandom":      "random",
	"bedreldist":     "reldist",
	"bedsample":      "sample",
	"bedshuffle":     "shuffle",
	"bedspacing":     "spacing",
	"bedsplit":       "split",
	"bedsummary":     "summary",
	"bedunionbedg":   "unionbedg",
	"bedcluster":     "cluster",
	"bedlinks":       "links",
	"bedigv":         "igv",
	"bedtag":         "tag",
}

// isBedTool reports whether tool is one of our per-bed binaries.
func isBedTool(tool string) bool {
	return strings.HasPrefix(tool, "bed") || tool == "bed12tobed6"
}

// binset resolves, per tool, the path to our binary and the upstream binary plus
// any leading upstream sub-argv (e.g. ["intersect"] when the upstream is
// `bedtools intersect`). An empty Path means the binary could not be resolved
// and cells using it SKIP.
type resolvedBin struct {
	Path     string   // absolute path to the binary, or "" if unresolved
	UpStub   []string // leading upstream argv (e.g. the bedtools subcommand)
	IsPerl   bool     // true when the binary is a Perl script (prinseq) run via `perl`
	NoteWhen string   // note recorded when Path is empty
}

// BinResolver resolves our + upstream binaries for every tool a cell may name,
// caching builds of our binaries. ourDir / upDir, when non-empty, are searched
// first; otherwise our binaries are built from source and upstream binaries are
// resolved from the vendored reference_code/ tree.
type BinResolver struct {
	ourDir   string
	upDir    string
	ourCache string
	notes    *[]string

	ourCacheMap map[string]string // logical tool -> our built/located path
	upCacheMap  map[string]resolvedBin
}

// NewBinResolver returns a resolver. notes accumulates human-readable resolution
// messages (missing binaries, build failures) for the report. ourDir/upDir, when
// non-empty, are directories searched first for the named tool binaries;
// otherwise our binaries are built from tools/*/cmd into ourCache and upstream
// binaries are resolved from the vendored reference_code/ tree.
func NewBinResolver(ourDir, upDir, ourCache string, notes *[]string) *BinResolver {
	return &BinResolver{
		ourDir:      ourDir,
		upDir:       upDir,
		ourCache:    ourCache,
		notes:       notes,
		ourCacheMap: map[string]string{},
		upCacheMap:  map[string]resolvedBin{},
	}
}

// note appends a resolution message (deduplicated by content).
func (r *BinResolver) note(msg string) {
	for _, n := range *r.notes {
		if n == msg {
			return
		}
	}
	*r.notes = append(*r.notes, msg)
}

// ourBinaryName is the name of our binary file for a logical tool. For bed*
// tools it is the tool name itself (tools/<tool>/cmd/<tool>).
func ourBinaryName(tool string) string { return tool }

// ourBinary resolves our binary path for tool, building it from
// tools/<tool>/cmd/<tool> when -our-bin is empty, else looking it up in the
// supplied dir. Returns "" (and records a note) when it cannot be resolved.
func (r *BinResolver) ourBinary(tool string) string {
	if p, ok := r.ourCacheMap[tool]; ok {
		return p
	}
	var path string
	if r.ourDir != "" {
		path = existingIn(r.ourDir, ourBinaryName(tool))
		if path == "" {
			r.note("our " + tool + " not found in " + r.ourDir)
		}
	} else {
		p, err := upstream.OurBinary(tool, r.ourCache)
		if err != nil {
			r.note("building our " + tool + ": " + err.Error())
		} else {
			path = p
		}
	}
	r.ourCacheMap[tool] = path
	return path
}

// upstreamBinary resolves the upstream binary (and any leading sub-argv) for a
// logical tool. bed* tools map to `bedtools <subcommand>`; samtools/bcftools/...
// map to their reference_code/ binaries; tools with no upstream pair return an
// empty Path (and the cell becomes ours-only).
func (r *BinResolver) upstreamBinary(tool string) resolvedBin {
	if rb, ok := r.upCacheMap[tool]; ok {
		return rb
	}
	rb := r.resolveUpstream(tool)
	r.upCacheMap[tool] = rb
	return rb
}

// resolveUpstream does the uncached upstream lookup for upstreamBinary.
func (r *BinResolver) resolveUpstream(tool string) resolvedBin {
	if isBedTool(tool) {
		sub, ok := bedSubcommand[tool]
		if !ok || sub == "" {
			return resolvedBin{NoteWhen: tool + " has no upstream bedtools pair (ours-only)"}
		}
		bin := r.locateUpstream("bedtools")
		if bin == "" {
			return resolvedBin{NoteWhen: "upstream bedtools not found"}
		}
		return resolvedBin{Path: bin, UpStub: []string{sub}}
	}
	bin := r.locateUpstream(tool)
	if bin == "" {
		return resolvedBin{NoteWhen: "upstream " + tool + " not found"}
	}
	return resolvedBin{Path: bin, IsPerl: strings.HasSuffix(bin, ".pl")}
}

// locateUpstream finds the upstream binary for an upstream key, preferring the
// -upstream-bin dir, then the vendored reference_code/ tree (via the upstream
// package). Returns "" when not found.
func (r *BinResolver) locateUpstream(key string) string {
	if r.upDir != "" {
		if p := existingIn(r.upDir, key); p != "" {
			return p
		}
	}
	if p, err := upstream.Binary(key); err == nil {
		return p
	} else if r.upDir == "" {
		r.note("upstream " + key + ": " + err.Error())
	}
	return ""
}

// existingIn returns DIR/name if it is an existing regular file, else "".
func existingIn(dir, name string) string {
	p := filepath.Join(dir, name)
	if fileExists(p) {
		return p
	}
	return ""
}
