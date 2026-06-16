// Shared -W/--write-index support for the native plugins (contrast, isecGT,
// mendelian2, scatter, split). Upstream bcftools' plugins accept
// -W/--write-index[=FMT] and, after each output file is fully written, build a
// CSI (default) or TBI index next to it — exactly as `bcftools index` would.
// This file centralises the flag parsing and the post-write indexing so every
// plugin reproduces upstream's behaviour byte-for-byte, including the error it
// raises when the requested output is not BGZF-compressed (and therefore not
// indexable).
package bcftools

import (
	"fmt"
	"strings"
)

// writeIndexFmt selects the index flavour a plugin's -W option requests.
type writeIndexFmt int

const (
	// writeIndexOff means -W/--write-index was not given.
	writeIndexOff writeIndexFmt = iota
	// writeIndexCSI emits a `.csi` index (the upstream default for -W).
	writeIndexCSI
	// writeIndexTBI emits a `.tbi` index (-W=tbi); valid only for VCF.gz.
	writeIndexTBI
)

// parseWriteIndexArg interprets a `-W` / `--write-index` token (and its
// attached `-W=FMT` / `--write-index=FMT` form). It returns the requested
// index flavour. A bare flag selects CSI (upstream's default); an explicit
// `=csi` / `=tbi` selects that flavour. The handled bool is false when arg is
// not a write-index option at all, so callers can fall through to their own
// option handling.
//
// Upstream accepts the suffix case-insensitively and treats an empty suffix
// (`-W=`) as the CSI default.
func parseWriteIndexArg(arg string) (fmtSel writeIndexFmt, handled bool, err error) {
	var suffix string
	switch {
	case arg == "-W" || arg == "--write-index":
		return writeIndexCSI, true, nil
	case strings.HasPrefix(arg, "-W="):
		suffix = arg[len("-W="):]
	case strings.HasPrefix(arg, "--write-index="):
		suffix = arg[len("--write-index="):]
	default:
		return writeIndexOff, false, nil
	}
	switch strings.ToLower(suffix) {
	case "", "csi":
		return writeIndexCSI, true, nil
	case "tbi":
		return writeIndexTBI, true, nil
	default:
		return writeIndexOff, true, fmt.Errorf("unknown --write-index format %q (expect csi or tbi)", suffix)
	}
}

// outputIsIndexable reports whether an output container can be indexed. Only
// BGZF-framed outputs — VCF.gz (-Oz) and compressed BCF (-Ob) — qualify.
// Plain VCF (-Ov) and our uncompressed-BCF (-Ou) writer produce non-BGZF
// streams, mirroring upstream's "Indexing is only supported on BGZF-compressed
// files" rejection.
func outputIsIndexable(format OutputFormat) bool {
	switch format {
	case OutputVCFGz, OutputBCF:
		return true
	default:
		return false
	}
}

// writeIndexFor builds the index for a freshly written output file at path,
// matching upstream's per-file -W behaviour. format is the container the file
// was written in. When the format is not BGZF-indexable, it returns the same
// error upstream emits (verified against the real binary) so the plugin exits
// with a clean failure rather than producing a bogus or missing index.
//
// A TBI index (-W=tbi) is valid only for VCF.gz; upstream silently writes a CSI
// instead for BCF output (TBI cannot describe a BCF), so we mirror that by
// falling back to CSI for any non-VCF.gz container.
func writeIndexFor(path string, format OutputFormat, sel writeIndexFmt) error {
	if sel == writeIndexOff {
		return nil
	}
	if !outputIsIndexable(format) {
		// Reproduce htslib's bcf_idx_init message + the plugin-level wrapper,
		// matching `bcftools +<plugin> ... -W` on a plain/uncompressed output.
		return fmt.Errorf("Indexing is only supported on BGZF-compressed files: failed to initialise index for %s", path)
	}
	idxFmt := IndexCSI
	if sel == writeIndexTBI && format == OutputVCFGz {
		idxFmt = IndexTBI
	}
	if _, err := BuildIndex(path, IndexOptions{Format: idxFmt, Force: true}); err != nil {
		return fmt.Errorf("failed to initialise index for %s: %w", path, err)
	}
	return nil
}
