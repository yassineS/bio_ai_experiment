// Package bed is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed`. See PR-A of
// the htsgo migration roadmap.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed`.
package bed

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

type (
	Record       = htsgo.Record
	Reader       = htsgo.Reader
	Writer       = htsgo.Writer
	BEDPE        = htsgo.BEDPE
	BEDPEReader  = htsgo.BEDPEReader
	BEDPEWriter  = htsgo.BEDPEWriter
	IntervalNode = htsgo.IntervalNode
	IntervalTree = htsgo.IntervalTree
)

var (
	NewReader       = htsgo.NewReader
	NewWriter       = htsgo.NewWriter
	NewBEDPEReader  = htsgo.NewBEDPEReader
	NewBEDPEWriter  = htsgo.NewBEDPEWriter
	NewIntervalTree = htsgo.NewIntervalTree
)
