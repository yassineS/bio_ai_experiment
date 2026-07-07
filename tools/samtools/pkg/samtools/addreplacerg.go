package samtools

import (
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
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
	// Mode selects orphan-only vs overwrite-all behaviour. Upstream's
	// default is overwrite_all (bam_addrprg.c sets retval->mode =
	// overwrite_all), so callers wanting parity should leave this at the
	// zero value only when they intend orphan_only.
	Mode AddReplaceRGMode
	// OverwriteHeaderRG mirrors upstream `-w`: when an @RG line with the
	// same ID already exists in the input header, replace it (remove the
	// old one and add the supplied -r line). Without it, an ID collision is
	// an error, matching upstream.
	OverwriteHeaderRG bool
	// NoPG is accepted; v1 never emits @PG lines so this is a no-op.
	NoPG bool
	// Threads is upstream's -@/--threads worker count. When > 1 it drives
	// block-parallel BGZF inflate on the input and block-parallel BGZF
	// deflate on the (always compressed) BAM output. Only the BGZF I/O is
	// parallelised — the RG-tagging pass is single-threaded, as in upstream
	// bam_addrprg.c — so the emitted records are byte-identical regardless of
	// the worker count. Parallel decode is opt-in (0/1 stays single-threaded)
	// because each worker adds block buffers to peak RSS.
	Threads int
}

// AddReplaceRG copies records from in to out, ensuring an @RG record
// exists in the header (adding RGLine when supplied) and setting the
// per-record RG aux tag to the appropriate ID per Mode.
func AddReplaceRG(in io.Reader, out io.Writer, opts AddReplaceRGOptions) error {
	if opts.RGLine == "" && opts.RGID == "" {
		return errors.New("samtools addreplacerg: need -r RG-line or -R RG-id")
	}
	br, err := alnio.NewReaderThreaded(in, "", ReadDecodeThreads(opts.Threads))
	if err != nil {
		return err
	}
	if rc, ok := br.(io.Closer); ok {
		defer rc.Close()
	}
	hdr := br.Header()

	id := opts.RGID
	if opts.RGLine != "" {
		newRG, err := parseRGLine(opts.RGLine)
		if err != nil {
			return err
		}
		id = newRG.ID
		// Upstream (bam_addrprg.c init_state): if an @RG with this ID is
		// already present, it is an error unless -w (OverwriteHeaderRG) was
		// given, in which case the existing line is removed first.
		if findRG(hdr, id) >= 0 {
			if !opts.OverwriteHeaderRG {
				return fmt.Errorf("samtools addreplacerg: @RG line with ID:%s already present in the header; use -w to overwrite", id)
			}
			removeRG(hdr, id)
		}
		// Add the new @RG line.
		hl := sam.HeaderLine{Tag: "RG", Fields: append([]sam.HeaderField{{Tag: "ID", Value: id}}, newRG.Extra...)}
		hdr.Lines = append(hdr.Lines, hl)
		hdr.ReadGroups = append(hdr.ReadGroups, newRG)
		// In overwrite_all mode upstream removes every other @RG line
		// (sam_hdr_remove_except keeps only the new ID).
		if opts.Mode == AddReplaceRGOverwriteAll {
			keepOnlyRG(hdr, id)
		}
	} else {
		if findRG(hdr, id) < 0 {
			return fmt.Errorf("samtools addreplacerg: -R %q is not in @RG table", id)
		}
	}

	// Output is always compressed BAM; -@ > 1 spreads its BGZF deflate across
	// the worker pool (byte-identical to the serial writer at the same level).
	var bw sam.Writer
	if opts.Threads > 1 {
		bw = sam.NewBAMWriterThreads(out, opts.Threads)
	} else {
		bw = sam.NewBAMWriter(out)
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

// AddReplaceRGFile is a path-based wrapper around AddReplaceRG. Under -@ >= 2
// it opens the input with iohelper.OpenRaw so the still-BGZF-framed bytes reach
// AddReplaceRG's block-parallel reader (NewReaderThreaded); the standard
// decompressing opener would inflate BGZF eagerly and the parallel input decode
// would never engage. This mirrors samtools calmd (CalmdFile) and sort. The
// decoded records are identical either way.
func AddReplaceRGFile(inPath string, out io.Writer, opts AddReplaceRGOptions) error {
	var in io.ReadCloser
	var err error
	if opts.Threads >= 2 {
		in, err = iohelper.OpenRaw(inPath)
	} else {
		in, err = iohelper.OpenReader(inPath)
	}
	if err != nil {
		return err
	}
	defer in.Close()
	return AddReplaceRG(in, out, opts)
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

// findRG returns the index of the @RG entry with ID id, or -1.
func findRG(h *sam.Header, id string) int {
	for i, rg := range h.ReadGroups {
		if rg.ID == id {
			return i
		}
	}
	return -1
}

// rgLineID returns the ID: field of an @RG header line, or "" if absent.
func rgLineID(hl sam.HeaderLine) string {
	for _, f := range hl.Fields {
		if f.Tag == "ID" {
			return f.Value
		}
	}
	return ""
}

// removeRG deletes the @RG line(s) with the given ID from both the ordered
// Lines slice (which drives header serialisation) and the typed ReadGroups
// slice, mirroring htslib's sam_hdr_remove_line_id.
func removeRG(h *sam.Header, id string) {
	lines := h.Lines[:0]
	for _, hl := range h.Lines {
		if hl.Tag == "RG" && rgLineID(hl) == id {
			continue
		}
		lines = append(lines, hl)
	}
	h.Lines = lines
	rgs := h.ReadGroups[:0]
	for _, rg := range h.ReadGroups {
		if rg.ID == id {
			continue
		}
		rgs = append(rgs, rg)
	}
	h.ReadGroups = rgs
}

// keepOnlyRG removes every @RG line whose ID is not id, mirroring htslib's
// sam_hdr_remove_except used by addreplacerg's overwrite_all mode.
func keepOnlyRG(h *sam.Header, id string) {
	lines := h.Lines[:0]
	for _, hl := range h.Lines {
		if hl.Tag == "RG" && rgLineID(hl) != id {
			continue
		}
		lines = append(lines, hl)
	}
	h.Lines = lines
	rgs := h.ReadGroups[:0]
	for _, rg := range h.ReadGroups {
		if rg.ID == id {
			rgs = append(rgs, rg)
		}
	}
	h.ReadGroups = rgs
}
