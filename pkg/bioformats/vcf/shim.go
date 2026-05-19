// Package vcf is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf`. See PR-B of
// the htsgo migration roadmap.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf`.
package vcf

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

type (
	Header  = htsgo.Header
	Reader  = htsgo.Reader
	Sample  = htsgo.Sample
	Variant = htsgo.Variant
	Writer  = htsgo.Writer
)

var (
	NewReader = htsgo.NewReader
	NewWriter = htsgo.NewWriter
)
