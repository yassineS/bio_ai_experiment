// Package bcftools — `bcftools call -G` sample-group loader.
//
// Upstream `bcftools call -G FILE` partitions the input samples into
// sub-populations so the multiallelic EM scores each group's qsum
// separately and the per-group best-allele sets are unioned into the
// site's final allele set. The file is a two-column whitespace-
// separated `SAMPLE\tGROUP` table; a single "-" placeholder means
// "every sample is its own group". See mcall.c::init_sample_groups
// (~lines 253-355).
package bcftools

import (
	"bufio"
	"fmt"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// SampleGroups holds the parsed per-group sample partition. Indices
// are resolved against a sample list (typically the input VCF
// header's Samples field) when Resolve is called.
type SampleGroups struct {
	// AllOwnGroup is true when the file argument was the literal "-",
	// meaning every sample should be placed in its own group at
	// resolve time.
	AllOwnGroup bool
	// Tag is the FORMAT tag whose per-sample counts feed the per-
	// group qsum (mirrors -G's optional ",TAG" suffix and the
	// implicit QS/AD fallback in mcall.c:277-285).
	Tag string
	// SampleToGroup maps each sample name to its group name. Samples
	// absent from the file are an error at Resolve time.
	SampleToGroup map[string]string
	// GroupOrder records the first-seen group name per appearance
	// order so the resolved group indices are stable / reproducible.
	GroupOrder []string
}

// LoadSampleGroups parses the upstream `-G FILE[,TAG]` argument. The
// literal "-" returns an AllOwnGroup SampleGroups with empty maps.
func LoadSampleGroups(arg string) (*SampleGroups, error) {
	sg := &SampleGroups{}
	// `-G FILE,TAG` syntax: the optional ",TAG" suffix selects the
	// per-sample counts tag. Matches mcall.c parse-time handling of
	// the option string (call->sample_groups_tag).
	if i := strings.IndexByte(arg, ','); i >= 0 {
		sg.Tag = strings.TrimSpace(arg[i+1:])
		arg = arg[:i]
	}
	arg = strings.TrimSpace(arg)
	if arg == "-" {
		sg.AllOwnGroup = true
		return sg, nil
	}
	if arg == "" {
		return nil, fmt.Errorf("bcftools call -G: empty file argument")
	}
	r, err := iohelper.OpenReader(arg)
	if err != nil {
		return nil, fmt.Errorf("bcftools call -G: open %s: %w", arg, err)
	}
	defer r.Close()
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 4*1024*1024)
	seen := map[string]struct{}{}
	sg.SampleToGroup = map[string]string{}
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || line[0] == '#' {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("bcftools call -G: bad line in %s, expected `SAMPLE<tab>GROUP`: %q", arg, line)
		}
		smpl, grp := fields[0], fields[1]
		if _, exists := sg.SampleToGroup[smpl]; exists {
			return nil, fmt.Errorf("bcftools call -G: the sample %q is listed twice in %s", smpl, arg)
		}
		sg.SampleToGroup[smpl] = grp
		if _, ok := seen[grp]; !ok {
			seen[grp] = struct{}{}
			sg.GroupOrder = append(sg.GroupOrder, grp)
		}
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("bcftools call -G: %w", err)
	}
	if len(sg.SampleToGroup) == 0 {
		return nil, fmt.Errorf("bcftools call -G: no samples found in %s", arg)
	}
	return sg, nil
}

// Resolve materialises the per-group sample-index partition for the
// supplied sample-name list (typically vcf.Header.Samples). Returns
// one entry per group with the resolved sample indices. The order
// matches GroupOrder for FILE inputs, or input sample order when
// AllOwnGroup is true.
type ResolvedSampleGroup struct {
	Name    string
	Indices []int
}

// Resolve returns one ResolvedSampleGroup per output group. Samples
// absent from the file partition return an error matching upstream's
// "The sample %s is not listed" rejection.
func (sg *SampleGroups) Resolve(samples []string) ([]ResolvedSampleGroup, error) {
	if sg == nil {
		out := make([]ResolvedSampleGroup, 1)
		out[0].Name = ""
		out[0].Indices = make([]int, len(samples))
		for i := range samples {
			out[0].Indices[i] = i
		}
		return out, nil
	}
	if sg.AllOwnGroup {
		out := make([]ResolvedSampleGroup, len(samples))
		for i, s := range samples {
			out[i] = ResolvedSampleGroup{Name: s, Indices: []int{i}}
		}
		return out, nil
	}
	idxByName := make(map[string]int, len(samples))
	for i, s := range samples {
		idxByName[s] = i
	}
	groupIdx := map[string]int{}
	out := make([]ResolvedSampleGroup, 0, len(sg.GroupOrder))
	for _, g := range sg.GroupOrder {
		groupIdx[g] = len(out)
		out = append(out, ResolvedSampleGroup{Name: g})
	}
	// Walk the sample names in input order so each group's Indices
	// stays sorted ascending — mcall.c builds them the same way via
	// the per-sample smpl2grp[] linear pass.
	for i, s := range samples {
		g, ok := sg.SampleToGroup[s]
		if !ok {
			return nil, fmt.Errorf("bcftools call -G: the sample %q is not listed in the -G file", s)
		}
		gi := groupIdx[g]
		out[gi].Indices = append(out[gi].Indices, i)
	}
	return out, nil
}
