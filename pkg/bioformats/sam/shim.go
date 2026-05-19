// Package sam is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam`. See PR-B of
// the htsgo migration roadmap.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam`.
package sam

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

type (
	Aux         = htsgo.Aux
	BAMReader   = htsgo.BAMReader
	BAMWriter   = htsgo.BAMWriter
	Cigar       = htsgo.Cigar
	CigarOp     = htsgo.CigarOp
	Header      = htsgo.Header
	HeaderField = htsgo.HeaderField
	HeaderLine  = htsgo.HeaderLine
	Program     = htsgo.Program
	ReadGroup   = htsgo.ReadGroup
	Reader      = htsgo.Reader
	Record      = htsgo.Record
	Reference   = htsgo.Reference
	SAMReader   = htsgo.SAMReader
	SAMWriter   = htsgo.SAMWriter
	Writer      = htsgo.Writer
)

const (
	FlagPaired        = htsgo.FlagPaired
	FlagProperPair    = htsgo.FlagProperPair
	FlagUnmapped      = htsgo.FlagUnmapped
	FlagMateUnmapped  = htsgo.FlagMateUnmapped
	FlagReverse       = htsgo.FlagReverse
	FlagMateReverse   = htsgo.FlagMateReverse
	FlagRead1         = htsgo.FlagRead1
	FlagRead2         = htsgo.FlagRead2
	FlagSecondary     = htsgo.FlagSecondary
	FlagQCFail        = htsgo.FlagQCFail
	FlagDuplicate     = htsgo.FlagDuplicate
	FlagSupplementary = htsgo.FlagSupplementary

	CigarMatch     = htsgo.CigarMatch
	CigarInsertion = htsgo.CigarInsertion
	CigarDeletion  = htsgo.CigarDeletion
	CigarSkipped   = htsgo.CigarSkipped
	CigarSoftClip  = htsgo.CigarSoftClip
	CigarHardClip  = htsgo.CigarHardClip
	CigarPadding   = htsgo.CigarPadding
	CigarEqual     = htsgo.CigarEqual
	CigarMismatch  = htsgo.CigarMismatch
	CigarBack      = htsgo.CigarBack
)

var (
	BAMMagic         = htsgo.BAMMagic
	ErrInvalidHeader = htsgo.ErrInvalidHeader
	ErrNotBAM        = htsgo.ErrNotBAM

	NewReader        = htsgo.NewReader
	NewSAMReader     = htsgo.NewSAMReader
	NewSAMWriter     = htsgo.NewSAMWriter
	NewBAMReader     = htsgo.NewBAMReader
	NewBAMBodyReader = htsgo.NewBAMBodyReader
	NewBAMWriter     = htsgo.NewBAMWriter
	ParseAux         = htsgo.ParseAux
	ParseCigar       = htsgo.ParseCigar
	ParseHeader      = htsgo.ParseHeader
	ParseHeaderText  = htsgo.ParseHeaderText
)
