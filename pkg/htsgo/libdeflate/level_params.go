package libdeflate

// Compression-level parameter tables ported from libdeflate's
// set_compression_level (deflate_compress.c:3920+). Slice 2 only
// exercises level 6, but we keep the full table so callers can pick
// other lazy/lazy2 levels once those are wired in.

// levelImpl identifies the parsing strategy a level selects.
type levelImpl int

const (
	implFastest levelImpl = iota // level 1
	implGreedy                   // levels 2-4
	implLazy                     // levels 5-7
	implLazy2                    // levels 8-9
)

// levelParams holds the per-level matchfinder/cost knobs.
type levelParams struct {
	impl            levelImpl
	maxSearchDepth  uint32 // unused for implFastest
	niceMatchLength uint32
}

// levelTable mirrors set_compression_level in deflate_compress.c. Index 0
// is the "no compression" entry (not used by Slice 2). Index 9 is the
// highest non-near-optimal level — levels 10-12 require near-optimal
// parsing which is out of scope for Slice 2 and will be addressed in a
// later slice.
var levelTable = [10]levelParams{
	0: {impl: implFastest, maxSearchDepth: 0, niceMatchLength: 0},
	1: {impl: implFastest, maxSearchDepth: 0, niceMatchLength: 32},
	2: {impl: implGreedy, maxSearchDepth: 6, niceMatchLength: 10},
	3: {impl: implGreedy, maxSearchDepth: 12, niceMatchLength: 14},
	4: {impl: implGreedy, maxSearchDepth: 16, niceMatchLength: 30},
	5: {impl: implLazy, maxSearchDepth: 16, niceMatchLength: 30},
	6: {impl: implLazy, maxSearchDepth: 35, niceMatchLength: 65},
	7: {impl: implLazy, maxSearchDepth: 100, niceMatchLength: 130},
	8: {impl: implLazy2, maxSearchDepth: 300, niceMatchLength: maxMatchLen},
	9: {impl: implLazy2, maxSearchDepth: 600, niceMatchLength: maxMatchLen},
}
