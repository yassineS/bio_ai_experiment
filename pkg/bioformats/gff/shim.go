// Package gff is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff`. See PR-A of
// the htsgo migration roadmap.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff`.
package gff

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
)

type (
	Strand  = htsgo.Strand
	Feature = htsgo.Feature
	Reader  = htsgo.Reader
)

const (
	StrandUnknown = htsgo.StrandUnknown
	StrandForward = htsgo.StrandForward
	StrandReverse = htsgo.StrandReverse
)

// NewReader is re-exported from pkg/htsgo/gff.
var NewReader = htsgo.NewReader
