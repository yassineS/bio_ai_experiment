// Package fasta is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta`. See PR-A of
// the htsgo migration roadmap.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta`.
package fasta

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

type (
	Record       = htsgo.Record
	Reader       = htsgo.Reader
	Writer       = htsgo.Writer
	IndexEntry   = htsgo.IndexEntry
	Index        = htsgo.Index
	RandomAccess = htsgo.RandomAccess
)

var (
	NewReader                      = htsgo.NewReader
	NewWriter                      = htsgo.NewWriter
	LoadIndex                      = htsgo.LoadIndex
	ReadIndex                      = htsgo.ReadIndex
	BuildIndex                     = htsgo.BuildIndex
	BuildIndexFullHeader           = htsgo.BuildIndexFullHeader
	BuildIndexBytes                = htsgo.BuildIndexBytes
	BuildIndexFullHeaderBytes      = htsgo.BuildIndexFullHeaderBytes
	OpenRandomAccess               = htsgo.OpenRandomAccess
	OpenRandomAccessFullHeader     = htsgo.OpenRandomAccessFullHeader
	NewRandomAccess                = htsgo.NewRandomAccess
	OpenRandomAccessBGZF           = htsgo.OpenRandomAccessBGZF
	OpenRandomAccessBGZFFullHeader = htsgo.OpenRandomAccessBGZFFullHeader
)
