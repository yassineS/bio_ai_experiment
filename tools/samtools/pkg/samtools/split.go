package samtools

import (
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
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
	// NoPG is accepted; v1 never injects @PG lines.
	NoPG bool
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
		if werr := bw.WriteHeader(hdr); werr != nil {
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
		if werr := bw.WriteHeader(hdr); werr != nil {
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
