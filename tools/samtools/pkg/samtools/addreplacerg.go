package samtools

import (
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// AddReplaceRGMode controls how the per-record RG tag is set.
type AddReplaceRGMode int

const (
	// AddReplaceRGOrphanOnly only sets the RG tag on records that don't
	// already carry one (`-m orphan_only`). This is the default.
	AddReplaceRGOrphanOnly AddReplaceRGMode = iota
	// AddReplaceRGOverwriteAll replaces the RG tag on every record
	// (`-m overwrite_all`).
	AddReplaceRGOverwriteAll
)

// AddReplaceRGOptions configures AddReplaceRG. Either RGLine (a full
// "@RG\tID:...\t..." block) or RGID (an existing RG ID) must be set;
// when both are present, RGLine takes precedence and the line's ID
// becomes the per-record value.
type AddReplaceRGOptions struct {
	// RGLine is a full @RG line string with TAB-separated fields.
	// Mirrors `-r/--rg-line`.
	RGLine string
	// RGID is the RG ID to apply when adopting an existing line from
	// the header. Mirrors `-R/--rg-id`.
	RGID string
	// Mode selects orphan-only vs overwrite-all behaviour.
	Mode AddReplaceRGMode
	// NoPG suppresses the @PG line injection. By default, addreplacerg
	// appends a @PG line documenting the command via the shared
	// InjectPG helper.
	NoPG bool
	// PGCommand is the raw command-line stored under @PG:CL when NoPG
	// is false. The CLI populates this with os.Args.
	PGCommand string
	// OutputBAM forces BAM output. Upstream's `addreplacerg` defaults
	// to SAM text when writing to stdout; set this to true (via
	// -O bam or an .bam output suffix at the CLI layer) to emit BAM.
	OutputBAM bool
}

// AddReplaceRG copies records from in to out, ensuring an @RG record
// exists in the header (adding RGLine when supplied) and setting the
// per-record RG aux tag to the appropriate ID per Mode.
func AddReplaceRG(in io.Reader, out io.Writer, opts AddReplaceRGOptions) error {
	if opts.RGLine == "" && opts.RGID == "" {
		return errors.New("samtools addreplacerg: need -r RG-line or -R RG-id")
	}
	br, err := sam.NewReader(in)
	if err != nil {
		return err
	}
	hdr := br.Header()

	id := opts.RGID
	if opts.RGLine != "" {
		newRG, err := parseRGLine(opts.RGLine)
		if err != nil {
			return err
		}
		id = newRG.ID
		// Append to header if not already present.
		if findRG(hdr, id) < 0 {
			hl := sam.HeaderLine{Tag: "RG", Fields: append([]sam.HeaderField{{Tag: "ID", Value: id}}, newRG.Extra...)}
			hdr.Lines = append(hdr.Lines, hl)
			hdr.ReadGroups = append(hdr.ReadGroups, newRG)
		}
	} else {
		if findRG(hdr, id) < 0 {
			return fmt.Errorf("samtools addreplacerg: -R %q is not in @RG table", id)
		}
	}
	// In overwrite_all mode upstream's sam_hdr_remove_except("RG",
	// "ID", id) drops every @RG line except the chosen one — so the
	// header reflects the post-rewrite reality.
	if opts.Mode == AddReplaceRGOverwriteAll {
		removeOtherRGs(hdr, id)
	}

	hdr = InjectPG(hdr, "samtools", "samtools", "0.1.0", opts.PGCommand, opts.NoPG)
	var bw sam.Writer
	if opts.OutputBAM {
		bw = sam.NewBAMWriter(out)
	} else {
		bw = sam.NewSAMWriter(out)
	}
	if err := bw.WriteHeader(hdr); err != nil {
		return err
	}
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		setRecordRG(rec, id, opts.Mode)
		if err := bw.Write(rec); err != nil {
			return err
		}
	}
	return bw.Close()
}

// setRecordRG sets the RG aux tag on rec per mode. In orphan-only mode
// an existing RG tag is left untouched.
func setRecordRG(rec *sam.Record, id string, mode AddReplaceRGMode) {
	idx := -1
	for i, a := range rec.Aux {
		if a.Tag == "RG" {
			idx = i
			break
		}
	}
	if idx >= 0 {
		if mode == AddReplaceRGOverwriteAll {
			rec.Aux[idx].Value = id
			rec.Aux[idx].Type = 'Z'
		}
		return
	}
	rec.Aux = append(rec.Aux, sam.Aux{Tag: "RG", Type: 'Z', Value: id})
}

// parseRGLine parses a literal `@RG\tID:...\t...` line into a
// sam.ReadGroup. The leading `@RG\t` is optional — both
// `@RG\tID:foo\tSM:bar` and `ID:foo\tSM:bar` are accepted, matching
// upstream's permissiveness.
func parseRGLine(line string) (sam.ReadGroup, error) {
	s := line
	// Allow literal "\t" sequences in addition to real TABs (CLI
	// arguments often quote tabs as backslash-t).
	s = strings.ReplaceAll(s, `\t`, "\t")
	if strings.HasPrefix(s, "@RG\t") {
		s = strings.TrimPrefix(s, "@RG\t")
	} else {
		s = strings.TrimPrefix(s, "@RG")
		s = strings.TrimPrefix(s, "\t")
	}
	if s == "" {
		return sam.ReadGroup{}, errors.New("samtools addreplacerg: empty @RG line")
	}
	rg := sam.ReadGroup{}
	for _, part := range strings.Split(s, "\t") {
		if len(part) < 4 || part[2] != ':' {
			return sam.ReadGroup{}, fmt.Errorf("samtools addreplacerg: bad RG field %q", part)
		}
		tag := part[:2]
		val := part[3:]
		if tag == "ID" {
			rg.ID = val
		} else {
			rg.Extra = append(rg.Extra, sam.HeaderField{Tag: tag, Value: val})
		}
	}
	if rg.ID == "" {
		return sam.ReadGroup{}, errors.New("samtools addreplacerg: @RG line missing ID:")
	}
	return rg, nil
}

// removeOtherRGs drops every @RG header line except the one whose
// ID equals keepID. Mirrors upstream sam_hdr_remove_except("RG",
// "ID", keepID).
func removeOtherRGs(h *sam.Header, keepID string) {
	rgs := h.ReadGroups[:0]
	for _, rg := range h.ReadGroups {
		if rg.ID == keepID {
			rgs = append(rgs, rg)
		}
	}
	h.ReadGroups = rgs
	lines := h.Lines[:0]
	for _, l := range h.Lines {
		if l.Tag != "RG" {
			lines = append(lines, l)
			continue
		}
		keep := false
		for _, f := range l.Fields {
			if f.Tag == "ID" && f.Value == keepID {
				keep = true
				break
			}
		}
		if keep {
			lines = append(lines, l)
		}
	}
	h.Lines = lines
}

// findRG returns the index of the @RG entry with ID id, or -1.
func findRG(h *sam.Header, id string) int {
	for i, rg := range h.ReadGroups {
		if rg.ID == id {
			return i
		}
	}
	return -1
}
