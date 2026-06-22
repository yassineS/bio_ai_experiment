package main

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// resolveBins resolves the four binaries (our + upstream, samtools + bcftools).
//
//   - Our binaries come from -our-bin DIR when it names existing samtools/bcftools
//     files; otherwise they are built into ourCache via upstream.OurBinary.
//   - Upstream binaries come from -upstream-bin DIR when it names them; otherwise
//     they are resolved from the vendored reference_code/ tree via
//     upstream.Binary.
//
// A binary that cannot be resolved is left empty; cells needing it SKIP (with a
// clear reason) rather than failing the whole run. resolveBins only returns an
// error for a filesystem problem creating the build cache.
func resolveBins(ourBinDir, upBinDir, ourCache string) (binset, []string, error) {
	var b binset
	var notes []string

	// Ours.
	if ourBinDir != "" {
		b.oursSamtools = existingIn(ourBinDir, "samtools")
		b.oursBcftools = existingIn(ourBinDir, "bcftools")
		if b.oursSamtools == "" {
			notes = append(notes, fmt.Sprintf("our samtools not found in %s", ourBinDir))
		}
		if b.oursBcftools == "" {
			notes = append(notes, fmt.Sprintf("our bcftools not found in %s", ourBinDir))
		}
	} else {
		if p, err := upstream.OurBinary("samtools", ourCache); err == nil {
			b.oursSamtools = p
		} else {
			notes = append(notes, "building our samtools: "+err.Error())
		}
		if p, err := upstream.OurBinary("bcftools", ourCache); err == nil {
			b.oursBcftools = p
		} else {
			notes = append(notes, "building our bcftools: "+err.Error())
		}
	}

	// Upstream.
	if upBinDir != "" {
		b.upSamtools = existingIn(upBinDir, "samtools")
		b.upBcftools = existingIn(upBinDir, "bcftools")
		if b.upSamtools == "" {
			notes = append(notes, fmt.Sprintf("upstream samtools not found in %s", upBinDir))
		}
		if b.upBcftools == "" {
			notes = append(notes, fmt.Sprintf("upstream bcftools not found in %s", upBinDir))
		}
	} else {
		if p, err := upstream.Binary("samtools"); err == nil {
			b.upSamtools = p
		} else {
			notes = append(notes, "upstream samtools: "+err.Error())
		}
		if p, err := upstream.Binary("bcftools"); err == nil {
			b.upBcftools = p
		} else {
			notes = append(notes, "upstream bcftools: "+err.Error())
		}
	}
	return b, notes, nil
}

// existingIn returns DIR/name if it is an existing regular file, else "".
func existingIn(dir, name string) string {
	p := filepath.Join(dir, name)
	if st, err := os.Stat(p); err == nil && !st.IsDir() {
		return p
	}
	return ""
}
