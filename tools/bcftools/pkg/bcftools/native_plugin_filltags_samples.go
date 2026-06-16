// -S/--samples-file parsing for the fill-tags plugin, porting parse_samples
// from plugins/fill-tags.c.
//
// The file has one sample per line: the sample name in the first
// whitespace-separated column and a comma-separated list of population/group
// names in the second column, e.g.
//
//	NA12400 GRP1
//	NA18507 GRP1,GRP2
//
// Each distinct group becomes a population whose tags are suffixed with
// "_GROUP". Samples not present in the VCF are skipped with a warning;
// duplicate sample lines are skipped with a warning. The summary "ALL"
// population is added by the caller (Init), not here.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// parseFillTagsSamples reads the sample/group file and returns the populations
// (in first-seen group order). The active-sample mask of each population is
// sized to the number of samples in hdr. Warnings about missing or duplicate
// samples and a sample-count mismatch are written to stderr (if non-nil),
// matching upstream's fprintf(stderr, ...) diagnostics.
func parseFillTagsSamples(fname string, hdr *vcf.Header, stderr io.Writer) ([]population, error) {
	f, err := os.Open(fname)
	if err != nil {
		return nil, fmt.Errorf("fill-tags: could not read: %s", fname)
	}
	defer f.Close()

	nsmpl := len(hdr.Samples)
	smplIdx := make(map[string]int, nsmpl)
	for i, s := range hdr.Samples {
		smplIdx[s] = i
	}

	var pops []population
	popIdx := make(map[string]int)
	seen := make(map[string]bool)
	nseen := 0

	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), " \t\r\n")
		if line == "" {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("fill-tags: could not parse the file: %s", line)
		}
		smpl := fields[0]
		groupList := fields[len(fields)-1]

		idx, ok := smplIdx[smpl]
		if !ok {
			if stderr != nil {
				fmt.Fprintf(stderr, "Warning: The sample not present in the VCF: %s\n", smpl)
			}
			continue
		}
		if seen[smpl] {
			if stderr != nil {
				fmt.Fprintf(stderr, "Warning: The sample is listed twice in %s: %s\n", fname, smpl)
			}
			continue
		}
		seen[smpl] = true

		for _, g := range strings.Split(groupList, ",") {
			if g == "" {
				continue
			}
			pi, exists := popIdx[g]
			if !exists {
				pi = len(pops)
				popIdx[g] = pi
				pops = append(pops, population{
					name:   g,
					suffix: "_" + g,
					mask:   make([]bool, nsmpl),
				})
			}
			pops[pi].mask[idx] = true
		}
		nseen++
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("fill-tags: error reading %s: %w", fname, err)
	}

	if nseen != nsmpl && stderr != nil {
		fmt.Fprintf(stderr, "Warning: %d samples in the list, %d samples in the VCF.\n", nseen, nsmpl)
	}
	if len(pops) == 0 {
		return nil, fmt.Errorf("fill-tags: no populations given?")
	}
	return pops, nil
}
