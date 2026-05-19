// Shim re-exporting the region-parser surface from pkg/htsgo/region
// and the BAI-specific UnionChunks helper from pkg/htsgo/bam. PR-F of
// the htsgo migration relocated `Region`/`ParseRegion`/`ResolvedRegion`/
// `ResolveRegions` into the shared htsgo home; `UnionChunks` moved
// alongside BAIIndex/BAIChunk in `pkg/htsgo/bam` because it depends on
// them. Existing in-package code (coverage.go, consensus.go,
// mpileup.go, depth.go, view.go, parity_test.go, view_test.go) keeps
// working through these aliases until PR-I deletes the shim.

package samtools

import (
	htsgobam "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
)

type (
	Region         = region.Region
	ResolvedRegion = region.ResolvedRegion
)

var (
	ParseRegion    = region.ParseRegion
	ResolveRegions = region.ResolveRegions
	UnionChunks    = htsgobam.UnionChunks
)
