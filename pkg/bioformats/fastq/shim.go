// Package fastq is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq`. See PR-A of
// the htsgo migration roadmap; new code should import the htsgo path
// directly.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq`.
package fastq

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// Type aliases preserve method sets and struct-literal compatibility,
// so callers using `fastq.Reader{}` etc. continue to work transparently.
type (
	QualityEncoding = htsgo.QualityEncoding
	Record          = htsgo.Record
	Reader          = htsgo.Reader
	Writer          = htsgo.Writer
)

// Enum constants must be re-declared at the matching type — Go const
// aliases (introduced 1.24) work, but explicit declaration keeps this
// shim buildable on Go 1.21 (the project's CI baseline).
const (
	Phred33 = htsgo.Phred33
	Phred64 = htsgo.Phred64
)

// NewReader is re-exported from pkg/htsgo/fastq.
var NewReader = htsgo.NewReader

// NewWriter is re-exported from pkg/htsgo/fastq.
var NewWriter = htsgo.NewWriter
