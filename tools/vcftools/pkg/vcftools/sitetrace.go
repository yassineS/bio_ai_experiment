// Site-trace outputs for --kept-sites and --removed-sites.
//
// Both flags emit a 2-column TSV (CHROM, POS) listing the sites that either
// survived (--kept-sites → `<prefix>.kept.sites`) or were filtered out
// (--removed-sites → `<prefix>.removed.sites`). They mirror upstream
// vcftools' output_kept_sites and output_removed_sites
// (variant_file_output.cpp:4285-4373) byte-for-byte: same header text,
// same column separator, same row order (input file order), same use of
// the 1-based POS.
//
// Upstream re-parses the entire VCF for each flag, applying apply_filters
// (entry_filters.cpp:23) to each entry. We instead piggy-back on the main
// filter pipeline in Run(): every site is either fed to stats.addVariant
// (kept) or hits a continue (removed), and we record both outcomes here.
// This is equivalent because our continue points cover the same filter
// gates as upstream's apply_filters; the row order also matches because
// we iterate the input in the same order.
package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

// siteTracker emits the .kept.sites and/or .removed.sites trace files for
// --kept-sites / --removed-sites. Either writer can be nil; recordKept and
// recordRemoved are no-ops when both are nil. close flushes and closes any
// underlying files (it is safe to call when there is nothing to close).
type siteTracker struct {
	keptW    *bufio.Writer
	keptF    io.Closer
	removedW *bufio.Writer
	removedF io.Closer
}

// newSiteTracker opens the .kept.sites and/or .removed.sites files
// requested by params and writes the upstream-matching header line
// ("CHROM\tPOS\n") into each. The header is emitted eagerly so that an
// invocation that filters every site (or no sites) still produces the
// header-only file upstream emits.
func newSiteTracker(prefix string, keep, remove bool) (*siteTracker, error) {
	t := &siteTracker{}
	if keep {
		f, err := os.Create(prefix + ".kept.sites")
		if err != nil {
			return nil, fmt.Errorf("opening --kept-sites output: %w", err)
		}
		t.keptF = f
		t.keptW = bufio.NewWriter(f)
		if _, err := t.keptW.WriteString("CHROM\tPOS\n"); err != nil {
			_ = f.Close()
			return nil, fmt.Errorf("writing --kept-sites header: %w", err)
		}
	}
	if remove {
		f, err := os.Create(prefix + ".removed.sites")
		if err != nil {
			// Roll back the kept-sites file so an error here doesn't
			// leak a half-initialised tracker.
			if t.keptF != nil {
				_ = t.keptF.Close()
			}
			return nil, fmt.Errorf("opening --removed-sites output: %w", err)
		}
		t.removedF = f
		t.removedW = bufio.NewWriter(f)
		if _, err := t.removedW.WriteString("CHROM\tPOS\n"); err != nil {
			_ = f.Close()
			if t.keptF != nil {
				_ = t.keptF.Close()
			}
			return nil, fmt.Errorf("writing --removed-sites header: %w", err)
		}
	}
	return t, nil
}

// recordKept logs a site that passed every filter. No-op if --kept-sites
// is disabled.
func (t *siteTracker) recordKept(chrom string, pos int) error {
	if t == nil || t.keptW == nil {
		return nil
	}
	if _, err := fmt.Fprintf(t.keptW, "%s\t%d\n", chrom, pos); err != nil {
		return fmt.Errorf("writing --kept-sites row: %w", err)
	}
	return nil
}

// recordRemoved logs a site that was dropped by some filter. No-op if
// --removed-sites is disabled.
func (t *siteTracker) recordRemoved(chrom string, pos int) error {
	if t == nil || t.removedW == nil {
		return nil
	}
	if _, err := fmt.Fprintf(t.removedW, "%s\t%d\n", chrom, pos); err != nil {
		return fmt.Errorf("writing --removed-sites row: %w", err)
	}
	return nil
}

// close flushes and closes any underlying files. Safe to call on a nil
// receiver or when no outputs were requested.
func (t *siteTracker) close() error {
	if t == nil {
		return nil
	}
	var firstErr error
	if t.keptW != nil {
		if err := t.keptW.Flush(); err != nil && firstErr == nil {
			firstErr = fmt.Errorf("flushing --kept-sites: %w", err)
		}
	}
	if t.keptF != nil {
		if err := t.keptF.Close(); err != nil && firstErr == nil {
			firstErr = fmt.Errorf("closing --kept-sites: %w", err)
		}
	}
	if t.removedW != nil {
		if err := t.removedW.Flush(); err != nil && firstErr == nil {
			firstErr = fmt.Errorf("flushing --removed-sites: %w", err)
		}
	}
	if t.removedF != nil {
		if err := t.removedF.Close(); err != nil && firstErr == nil {
			firstErr = fmt.Errorf("closing --removed-sites: %w", err)
		}
	}
	return firstErr
}
