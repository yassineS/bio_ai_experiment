package samtools

import (
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// SplitOptions configures Split. The default output filename pattern is
// "%*_%!.bam" — i.e. `<basename>_<RG-id>.bam` next to the input.
type SplitOptions struct {
	// Pattern controls the output filename. Substitutions:
	//   %!  the read-group ID
	//   %*  the input file basename without extension
	//   %.  the input file extension (with leading '.', or "" when none)
	// Defaults to "%*_%!.bam" when empty.
	Pattern string
	// Unidentified is the path to write reads with no RG (or whose RG ID
	// is not in the @RG table). Empty means "discard them".
	Unidentified string
	// NoPG suppresses the @PG line injection. By default, split appends
	// a @PG line to every per-RG output documenting the command via the
	// shared InjectPG helper.
	NoPG bool
	// PGCommand is the raw command-line stored under @PG:CL when NoPG
	// is false. The CLI populates this with os.Args.
	PGCommand string
}

// splitOut wraps one open per-RG output file.
type splitOut struct {
	path string
	f    *os.File
	w    *sam.BAMWriter
}

// SplitFile splits the BAM at inPath into one output file per @RG ID
// recorded in its header. Records without an RG aux tag, or whose RG ID
// isn't in the @RG table, go to opts.Unidentified (or are dropped when
// that path is empty).
func SplitFile(inPath string, opts SplitOptions) error {
	in, err := os.Open(inPath)
	if err != nil {
		return err
	}
	defer in.Close()
	br, err := sam.NewBAMReader(in)
	if err != nil {
		return err
	}
	defer br.Close()
	hdr := br.Header()

	if opts.Pattern == "" {
		opts.Pattern = "%*_%!.bam"
	}
	base, ext := splitInputName(inPath)

	outs := make(map[string]*splitOut, len(hdr.ReadGroups))
	cleanup := func() {
		for _, of := range outs {
			_ = of.w.Close()
			_ = of.f.Close()
		}
	}
	for _, rg := range hdr.ReadGroups {
		path := expandSplitPattern(opts.Pattern, base, ext, rg.ID)
		f, ferr := os.Create(path)
		if ferr != nil {
			cleanup()
			return ferr
		}
		bw := sam.NewBAMWriter(f)
		// Upstream samtools split rewrites each per-RG output's header to
		// retain only the matching @RG line (sam_split.c). Match that so
		// downstream tools see a single-RG header per file.
		perRGHdr := headerKeepOnlyRG(hdr, rg.ID)
		perRGHdr = InjectPG(perRGHdr, "samtools", "samtools", "0.1.0", opts.PGCommand, opts.NoPG)
		if werr := bw.WriteHeader(perRGHdr); werr != nil {
			_ = f.Close()
			cleanup()
			return werr
		}
		outs[rg.ID] = &splitOut{path: path, f: f, w: bw}
	}

	var unident *splitOut
	if opts.Unidentified != "" {
		f, ferr := os.Create(opts.Unidentified)
		if ferr != nil {
			cleanup()
			return ferr
		}
		bw := sam.NewBAMWriter(f)
		unidentHdr := InjectPG(hdr, "samtools", "samtools", "0.1.0", opts.PGCommand, opts.NoPG)
		if werr := bw.WriteHeader(unidentHdr); werr != nil {
			_ = f.Close()
			cleanup()
			return werr
		}
		unident = &splitOut{path: opts.Unidentified, f: f, w: bw}
	}

	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			cleanup()
			if unident != nil {
				_ = unident.w.Close()
				_ = unident.f.Close()
			}
			return err
		}
		var rgID string
		if a, ok := rec.GetAux("RG"); ok {
			if s, sok := a.String(); sok {
				rgID = s
			}
		}
		if of, ok := outs[rgID]; ok {
			if werr := of.w.Write(rec); werr != nil {
				cleanup()
				if unident != nil {
					_ = unident.w.Close()
					_ = unident.f.Close()
				}
				return werr
			}
			continue
		}
		if unident != nil {
			if werr := unident.w.Write(rec); werr != nil {
				cleanup()
				_ = unident.w.Close()
				_ = unident.f.Close()
				return werr
			}
		}
	}
	for _, of := range outs {
		if cerr := of.w.Close(); cerr != nil {
			return cerr
		}
		if cerr := of.f.Close(); cerr != nil {
			return cerr
		}
	}
	if unident != nil {
		if cerr := unident.w.Close(); cerr != nil {
			return cerr
		}
		if cerr := unident.f.Close(); cerr != nil {
			return cerr
		}
	}
	return nil
}

// expandSplitPattern resolves the %!/%*/%. tokens in pat against the
// supplied basename, extension, and RG ID.
func expandSplitPattern(pat, base, ext, rgID string) string {
	var sb strings.Builder
	for i := 0; i < len(pat); i++ {
		if pat[i] == '%' && i+1 < len(pat) {
			switch pat[i+1] {
			case '!':
				sb.WriteString(rgID)
				i++
				continue
			case '*':
				sb.WriteString(base)
				i++
				continue
			case '.':
				sb.WriteString(ext)
				i++
				continue
			}
		}
		sb.WriteByte(pat[i])
	}
	return sb.String()
}

// headerKeepOnlyRG returns a shallow copy of h with every @RG line that
// does NOT match keepID dropped. The Lines slice is rebuilt so the BAM
// header text only carries the surviving @RG record; @SQ/@HD/@PG/@CO
// lines and ReadGroup order are otherwise preserved. The derived view
// slices (Refs, ReadGroups, Programs, Comments, HDFields) are updated to
// stay consistent with Lines.
func headerKeepOnlyRG(h *sam.Header, keepID string) *sam.Header {
	if h == nil {
		return nil
	}
	out := &sam.Header{
		Refs:     append([]sam.Reference(nil), h.Refs...),
		Programs: append([]sam.Program(nil), h.Programs...),
		Comments: append([]string(nil), h.Comments...),
		HDFields: append([]sam.HeaderField(nil), h.HDFields...),
	}
	for _, line := range h.Lines {
		if line.Tag == "RG" {
			id := ""
			for _, f := range line.Fields {
				if f.Tag == "ID" {
					id = f.Value
					break
				}
			}
			if id != keepID {
				continue
			}
		}
		out.Lines = append(out.Lines, line)
	}
	for _, rg := range h.ReadGroups {
		if rg.ID == keepID {
			out.ReadGroups = append(out.ReadGroups, rg)
		}
	}
	return out
}

// splitInputName returns (basename without extension, extension with
// leading dot or empty).
func splitInputName(p string) (string, string) {
	name := filepath.Base(p)
	dot := strings.LastIndex(name, ".")
	if dot < 0 {
		return name, ""
	}
	return name[:dot], name[dot:]
}
