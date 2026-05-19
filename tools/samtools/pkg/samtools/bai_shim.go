// Shim re-exporting the BAI format primitives from pkg/htsgo/bam.
// PR-D of the htsgo migration relocated bai.go + bai_test.go into the
// shared htsgo tree. Existing in-package code (and any external caller
// of the samtools package) keeps working through these aliases until
// PR-I deletes the shim.

package samtools

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
)

type (
	BAIChunk   = htsgo.BAIChunk
	BAIBin     = htsgo.BAIBin
	BAIRef     = htsgo.BAIRef
	BAIIndex   = htsgo.BAIIndex
	BAIBuilder = htsgo.BAIBuilder
)

const (
	MetaBin = htsgo.MetaBin
)

var (
	BAIMagic       = htsgo.BAIMagic
	ErrBadBAIMagic = htsgo.ErrBadBAIMagic
	NewBAIBuilder  = htsgo.NewBAIBuilder
	WriteBAI       = htsgo.WriteBAI
	ReadBAI        = htsgo.ReadBAI
)
