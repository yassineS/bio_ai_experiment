// Package iohelper is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper`. The real
// implementation moved as part of the htsgo migration (PR-A); this
// shim exists so existing callers under `pkg/bioformats/iohelper`
// continue to compile during the rolling migration. New code should
// import `pkg/htsgo/iohelper` directly.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper`.
package iohelper

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// OpenReader is re-exported from `pkg/htsgo/iohelper.OpenReader`.
var OpenReader = htsgo.OpenReader

// OpenWriter is re-exported from `pkg/htsgo/iohelper.OpenWriter`.
var OpenWriter = htsgo.OpenWriter
