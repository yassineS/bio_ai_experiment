// Shared PED / PFM parsing helpers for the native trio plugins (trio-stats,
// trio-switch-rate). These mirror the tiny htslib utilities those plugins use:
// the per-sample index lookup (bcf_hdr_id2int(BCF_DT_SAMPLE,...)), the
// whitespace-delimited PED parse (ksplit_core + the column convention
// familyID,sampleID,paternalID,maternalID,...), the -P/--pfm P,F,M parse, and
// the trio sort by minimum sample index (cmp_trios).
package bcftools

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// sampleIndex builds a name->column map for the header samples, mirroring
// htslib's bcf_hdr_id2int(BCF_DT_SAMPLE,name) lookup.
func sampleIndex(hdr *vcf.Header) map[string]int {
	m := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		m[s] = i
	}
	return m
}

// parseTrioPFM parses the -P/--pfm "P,F,M" value (proband, father, mother) and
// resolves the three sample indices, mirroring the manual comma parse in
// trio-stats.c init_data(). All three samples must be present in the header.
func parseTrioPFM(pfm string, idx map[string]int) (trioStatsTrio, error) {
	parts := strings.Split(pfm, ",")
	if len(parts) != 3 {
		return trioStatsTrio{}, fmt.Errorf("could not parse -P %s", pfm)
	}
	child, ok := idx[parts[0]]
	if !ok {
		return trioStatsTrio{}, fmt.Errorf("no such sample: %q", parts[0])
	}
	father, ok := idx[parts[1]]
	if !ok {
		return trioStatsTrio{}, fmt.Errorf("no such sample: %q", parts[1])
	}
	mother, ok := idx[parts[2]]
	if !ok {
		return trioStatsTrio{}, fmt.Errorf("no such sample: %q", parts[2])
	}
	return trioStatsTrio{child: child, father: father, mother: mother}, nil
}

// parseTrioStatsPED reads a PED file and returns the trios whose father, mother
// AND child are all present in the header, sorted by the minimum sample index,
// mirroring parse_ped() + cmp_trios() in trio-stats.c. Columns are whitespace-
// delimited; column 1 is sampleID (child), 2 paternalID, 3 maternalID
// (0-based offsets off[1], off[2], off[3]). Rows whose parents/child are not in
// the VCF are skipped; a child listed in two trios is an error, and a duplicate
// (sample,father,mother) trio is skipped. At least one complete trio must
// resolve.
func parseTrioStatsPED(path string, idx map[string]int) ([]trioStatsTrio, error) {
	rows, err := parsePEDRows(path, idx)
	if err != nil {
		return nil, fmt.Errorf("trio-stats: %w", err)
	}
	hasChild := make(map[string]bool)
	hasTrio := make(map[string]bool)
	var trios []trioStatsTrio
	for _, r := range rows {
		key := r.childName + " " + r.fatherName + " " + r.motherName
		if hasTrio[key] {
			continue
		}
		if hasChild[r.childName] {
			return nil, fmt.Errorf("trio-stats: the child %q is listed in two trios", r.childName)
		}
		hasChild[r.childName] = true
		hasTrio[key] = true
		trios = append(trios, trioStatsTrio{child: r.child, father: r.father, mother: r.mother})
	}
	if len(trios) == 0 {
		return nil, fmt.Errorf("trio-stats: no complete trio identified")
	}
	sortTriosByMinIndex(trios)
	return trios, nil
}

// parseIndelStatsPED reads a PED file for indel-stats' -p de-novo mode and
// returns the trios whose father, mother AND child are all present in the
// header, sorted by the minimum sample index. It mirrors indel-stats.c's
// parse_ped() + cmp_trios(), which — unlike trio-stats.c — does NOT deduplicate
// trios and does NOT reject a child listed in two trios: every resolvable row is
// appended in file order before the stable sort. At least one complete trio must
// resolve.
func parseIndelStatsPED(path string, idx map[string]int) ([]trioStatsTrio, error) {
	rows, err := parsePEDRows(path, idx)
	if err != nil {
		return nil, fmt.Errorf("indel-stats: %w", err)
	}
	var trios []trioStatsTrio
	for _, r := range rows {
		trios = append(trios, trioStatsTrio{child: r.child, father: r.father, mother: r.mother})
	}
	if len(trios) == 0 {
		return nil, fmt.Errorf("indel-stats: no complete trio identified")
	}
	sortTriosByMinIndex(trios)
	return trios, nil
}

// pedRow is one resolved PED line (indices plus names) used while de-duplicating
// trios.
type pedRow struct {
	child, father, mother             int
	childName, fatherName, motherName string
	sex                               int    // 5th column when present (1=male,2=female), else 0
	popName                           string // 7th column when present (population grouping), else ""
}

// parsePEDRows reads the whitespace-delimited PED file and returns the rows
// whose child/father/mother all resolve against idx, preserving file order.
// Lines with fewer than 4 columns are an error (matching upstream's
// ncols<4 check). The first line is treated as data, exactly as upstream
// (which reads it via the same do/while loop, not as a header).
func parsePEDRows(path string, idx map[string]int) ([]pedRow, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	empty := true
	var rows []pedRow
	for sc.Scan() {
		line := sc.Text()
		empty = false
		fields := strings.Fields(line)
		if len(fields) < 4 {
			return nil, fmt.Errorf("could not parse the ped file: %s", line)
		}
		father, fok := idx[fields[2]]
		if !fok {
			continue
		}
		mother, mok := idx[fields[3]]
		if !mok {
			continue
		}
		child, cok := idx[fields[1]]
		if !cok {
			continue
		}
		sex := 0
		if len(fields) >= 5 {
			// upstream's trio-switch-rate / mendelian2 parse column 5 as the sex;
			// trio-stats ignores it. We keep it for callers that care.
			switch fields[4] {
			case "1":
				sex = 1
			case "2":
				sex = 2
			}
		}
		popName := ""
		if len(fields) > 6 {
			popName = fields[6]
		}
		rows = append(rows, pedRow{
			child: child, father: father, mother: mother,
			childName: fields[1], fatherName: fields[2], motherName: fields[3],
			sex: sex, popName: popName,
		})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if empty {
		return nil, fmt.Errorf("empty file: %s", path)
	}
	return rows, nil
}

// sortTriosByMinIndex sorts trios ascending by the minimum of their three
// sample indices, mirroring cmp_trios() in the trio plugins (a stable sort
// preserves the relative order of trios sharing a minimum, matching qsort on
// already-ordered input for the common single-trio / disjoint-trio cases).
func sortTriosByMinIndex(trios []trioStatsTrio) {
	minIdx := func(t trioStatsTrio) int {
		m := t.child
		if t.father < m {
			m = t.father
		}
		if t.mother < m {
			m = t.mother
		}
		return m
	}
	sort.SliceStable(trios, func(i, j int) bool {
		return minIdx(trios[i]) < minIdx(trios[j])
	})
}
