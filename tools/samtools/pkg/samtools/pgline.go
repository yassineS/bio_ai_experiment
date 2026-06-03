// Package samtools — shared @PG header-line injector.
//
// Upstream samtools subcommands inject one @PG line into the output
// header documenting the command that produced the output (program
// name + version + raw command-line, plus a PP "previous program"
// pointer chaining onto whatever @PG records the input already
// carried). The `--no-PG` flag suppresses the injection.
//
// InjectPG appends the new @PG record to the header's ordered Lines
// slice and its Programs view, chaining onto the last existing @PG.
// The ID is uniquified against the existing chain so a
// `samtools calmd | samtools calmd` chain produces ID=calmd then
// ID=calmd.1 in upstream's exact pattern (see bam_aux.c
// sam_hdr_add_pg).
package samtools

import (
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// InjectPG appends a fresh @PG header line documenting program pn
// (e.g. "samtools") subcommand id (e.g. "calmd") with version vn and
// the raw command-line cl. The new record's ID is unique within the
// header's existing @PG chain; PP points at the previous record.
// If hdr is nil or noPG is true the call is a no-op and the
// original header is returned. Otherwise the modified header is
// returned (mutated in place: hdr.Lines and hdr.Programs are
// appended).
func InjectPG(hdr *sam.Header, id, pn, vn, cl string, noPG bool) *sam.Header {
	if hdr == nil || noPG {
		return hdr
	}
	uniqID := uniquifyPGID(hdr, id)
	fields := []sam.HeaderField{
		{Tag: "ID", Value: uniqID},
		{Tag: "PN", Value: pn},
	}
	if vn != "" {
		fields = append(fields, sam.HeaderField{Tag: "VN", Value: vn})
	}
	// PP chains onto the last existing @PG line, matching upstream's
	// sam_hdr_add_pg behaviour. Multi-chain inputs always pick the
	// final entry as the parent.
	if len(hdr.Programs) > 0 {
		fields = append(fields, sam.HeaderField{Tag: "PP", Value: hdr.Programs[len(hdr.Programs)-1].ID})
	}
	if cl != "" {
		fields = append(fields, sam.HeaderField{Tag: "CL", Value: cl})
	}
	hl := sam.HeaderLine{Tag: "PG", Fields: fields}
	hdr.Lines = append(hdr.Lines, hl)
	hdr.Programs = append(hdr.Programs, sam.Program{ID: uniqID, Extra: fields})
	return hdr
}

// uniquifyPGID returns id when no existing @PG line carries it;
// otherwise it appends a numeric suffix ".N" matching upstream's
// pattern (id, id.1, id.2, ...).
func uniquifyPGID(hdr *sam.Header, id string) string {
	taken := map[string]struct{}{}
	for _, p := range hdr.Programs {
		taken[p.ID] = struct{}{}
	}
	if _, exists := taken[id]; !exists {
		return id
	}
	for n := 1; ; n++ {
		candidate := id + "." + strconv.Itoa(n)
		if _, exists := taken[candidate]; !exists {
			return candidate
		}
	}
}
