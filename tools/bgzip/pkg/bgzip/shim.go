// Package bgzip is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf`. See PR-C of
// the htsgo migration roadmap.
//
// The new package name is `bgzf` (matching htslib's terminology and
// the htsgo target tree); the shim keeps the legacy `bgzip` name so
// existing callers under `tools/bgzip/pkg/bgzip` continue to compile.
// New code should import `pkg/htsgo/bgzf` directly. The shim and the
// `bgzip` package name disappear together in PR-I.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf`.
package bgzip

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

type (
	Reader      = htsgo.Reader
	Writer      = htsgo.Writer
	Block       = htsgo.Block
	BlockOffset = htsgo.BlockOffset
)

const (
	DefaultCompression     = htsgo.DefaultCompression
	MaxBlockSize           = htsgo.MaxBlockSize
	MaxCompressedBlockSize = htsgo.MaxCompressedBlockSize
)

var (
	EOFBlock = htsgo.EOFBlock

	// Sentinel errors — kept as `var` so `errors.Is` callers continue
	// to compare against the same underlying *errorString values
	// regardless of which import path they use.
	ErrBadMagic     = htsgo.ErrBadMagic
	ErrNoExtra      = htsgo.ErrNoExtra
	ErrNoBCSubfield = htsgo.ErrNoBCSubfield
	ErrBadBSIZE     = htsgo.ErrBadBSIZE
	ErrTruncated    = htsgo.ErrTruncated
	ErrChecksum     = htsgo.ErrChecksum
	ErrISIZE        = htsgo.ErrISIZE

	NewReader            = htsgo.NewReader
	NewWriter            = htsgo.NewWriter
	NewWriterLevel       = htsgo.NewWriterLevel
	DecompressedSize     = htsgo.DecompressedSize
	UncompressedOffsetAt = htsgo.UncompressedOffsetAt
	Scan                 = htsgo.Scan
	ReadGZI              = htsgo.ReadGZI
	WriteGZI             = htsgo.WriteGZI
)
