// Package tabix is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix`. See PR-E of
// the htsgo migration roadmap.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix`.
package tabix

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

type (
	Bin      = htsgo.Bin
	Chunk    = htsgo.Chunk
	RefIndex = htsgo.RefIndex
	Config   = htsgo.Config
	Index    = htsgo.Index
	VOffset  = htsgo.VOffset
	CSI      = htsgo.CSI
	CSIBin   = htsgo.CSIBin
	CSIChunk = htsgo.CSIChunk
	CSIRef   = htsgo.CSIRef
)

const (
	MaxBin   = htsgo.MaxBin
	TileSize = htsgo.TileSize

	FormatGeneric = htsgo.FormatGeneric
	FormatSAM     = htsgo.FormatSAM
	FormatVCF     = htsgo.FormatVCF

	FlagZeroBased    = htsgo.FlagZeroBased
	FlagEndInclusive = htsgo.FlagEndInclusive

	PresetGFF = htsgo.PresetGFF
	PresetBED = htsgo.PresetBED
	PresetSAM = htsgo.PresetSAM
	PresetVCF = htsgo.PresetVCF
)

var (
	Magic    = htsgo.Magic
	CSIMagic = htsgo.CSIMagic

	Reg2bin              = htsgo.Reg2bin
	Reg2bins             = htsgo.Reg2bins
	LinearTile           = htsgo.LinearTile
	MakeVOffset          = htsgo.MakeVOffset
	NewIndex             = htsgo.NewIndex
	NewCSI               = htsgo.NewCSI
	PresetConfig         = htsgo.PresetConfig
	Build                = htsgo.Build
	BuildCSIFromDataFile = htsgo.BuildCSIFromDataFile
	Read                 = htsgo.Read
	ReadFile             = htsgo.ReadFile
	ReadCSI              = htsgo.ReadCSI
	ReadCSIFile          = htsgo.ReadCSIFile
	ChromIDInCSI         = htsgo.ChromIDInCSI
)
