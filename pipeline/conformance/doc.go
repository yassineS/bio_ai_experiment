// Package conformance runs the originals' OWN test corpora through our
// re-implemented binaries and asserts we accept/reject and reproduce them
// consistently with upstream.
//
// Two corpora are wired in:
//
//   - reference_code/htslib/test/ — SAM/BAM/CRAM/VCF fixtures shipped by
//     htslib itself (empty files, CRLF input, padded/bad CIGAR, no-SEQ
//     records, unknown references, BGZF block-boundary cases, long
//     references). These are fed through our samtools (SAM→BAM→SAM round
//     trips and CRAM round trips) and compared against upstream samtools.
//
//   - reference_code/htscodecs/tests/ — the rANS / arithmetic / fqzcomp
//     compression compliance vectors. These are round-tripped through our
//     in-tree CRAM codec (pkg/htsgo/cram/codec) and asserted byte-identical
//     to the reference vectors.
//
// Every test SKIPs gracefully (never fails) when its corpus submodule is not
// initialised; docs/CONFORMANCE.md documents the `git submodule update`
// commands that populate them.
package conformance
