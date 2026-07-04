package samtools

import (
	"bufio"
	"crypto/md5"
	"encoding/hex"
	"fmt"
	"io"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// DictOptions configures the dict subcommand. The defaults match upstream
// `samtools dict`: emit @HD followed by one @SQ per FASTA record, with
// SN/LN/M5 always populated and AS/SP/UR optional.
type DictOptions struct {
	// Assembly populates the @SQ AS: tag (`-a NAME`).
	Assembly string
	// Species populates the @SQ SP: tag (`-s SPECIES`).
	Species string
	// URI populates the @SQ UR: tag (`-u URI`). When empty and Path is a
	// real file path, upstream uses "file://<absolute-path>".
	URI string
	// AliasFromHeader, when true (-A), emits an AN: alias tag derived by
	// adding or removing a leading "chr" prefix on the sequence name,
	// matching upstream `samtools dict -A` (including the chrM/chrMT
	// mitochondrial synonyms).
	AliasFromHeader bool
	// NoHeader, when true (-H), suppresses the @HD line.
	NoHeader bool
}

// DictEntry is one @SQ row emitted by `samtools dict`.
type DictEntry struct {
	Name     string
	Length   int
	M5       string
	URI      string
	Assembly string
	Species  string
	Aliases  []string
}

// FormatSAM returns the textual @SQ line for this entry, terminated by a
// newline.
func (e DictEntry) FormatSAM() string {
	var sb strings.Builder
	sb.WriteString("@SQ\tSN:")
	sb.WriteString(e.Name)
	fmt.Fprintf(&sb, "\tLN:%d", e.Length)
	if e.M5 != "" {
		sb.WriteString("\tM5:")
		sb.WriteString(e.M5)
	}
	// Field order matches upstream `samtools dict` exactly:
	// SN, LN, M5, [AN], [UR], [AS], [SP]. (AH — for alt-locus
	// sequences — is not yet supported.)
	for _, al := range e.Aliases {
		sb.WriteString("\tAN:")
		sb.WriteString(al)
	}
	if e.URI != "" {
		sb.WriteString("\tUR:")
		sb.WriteString(e.URI)
	}
	if e.Assembly != "" {
		sb.WriteString("\tAS:")
		sb.WriteString(e.Assembly)
	}
	if e.Species != "" {
		sb.WriteString("\tSP:")
		sb.WriteString(e.Species)
	}
	sb.WriteByte('\n')
	return sb.String()
}

// aliasName builds the AN: alias value for name, matching upstream
// `samtools dict -A`: a leading "chr" prefix is stripped if present,
// otherwise one is added. Mitochondrial names additionally gain their
// well-known synonyms (M ↔ chrMT/MT and MT ↔ chrM/M) as comma-separated
// values within the single AN tag.
func aliasName(name string) string {
	var an string
	stripped := name
	if strings.HasPrefix(name, "chr") {
		stripped = name[3:]
		an = stripped
	} else {
		an = "chr" + name
	}
	switch stripped {
	case "M":
		an += ",chrMT,MT"
	case "MT":
		an += ",chrM,M"
	}
	return an
}

// Dict streams a FASTA file and emits the corresponding sequence
// dictionary (`@HD` + one `@SQ` per record) to w. The MD5 is computed
// over the uppercased, whitespace-stripped sequence — matching
// `samtools dict` (and Picard `CreateSequenceDictionary`).
func Dict(in io.Reader, w io.Writer, opts DictOptions) error {
	br := bufio.NewReader(in)
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	if !opts.NoHeader {
		if _, err := bw.WriteString("@HD\tVN:1.0\tSO:unsorted\n"); err != nil {
			return err
		}
	}

	var (
		entry  DictEntry
		hasher = md5.New()
		length int
		open   bool
	)

	flush := func() error {
		if !open {
			return nil
		}
		entry.Length = length
		entry.M5 = hex.EncodeToString(hasher.Sum(nil))
		_, err := bw.WriteString(entry.FormatSAM())
		return err
	}

	scan := bufio.NewScanner(br)
	scan.Buffer(make([]byte, 0, 64*1024), 32*1024*1024)
	for scan.Scan() {
		line := scan.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {
			if err := flush(); err != nil {
				return err
			}
			open = true
			length = 0
			hasher.Reset()
			head := strings.TrimSpace(string(line[1:]))
			fields := strings.Fields(head)
			entry = DictEntry{
				Assembly: opts.Assembly,
				Species:  opts.Species,
				URI:      opts.URI,
			}
			if len(fields) > 0 {
				entry.Name = fields[0]
				if opts.AliasFromHeader {
					entry.Aliases = []string{aliasName(entry.Name)}
				}
			}
			continue
		}
		// Sequence line — uppercase and feed to MD5.
		for _, c := range line {
			if c == ' ' || c == '\t' || c == '\r' || c == '\n' {
				continue
			}
			if c >= 'a' && c <= 'z' {
				c -= 32
			}
			hasher.Write([]byte{c})
			length++
		}
	}
	if err := scan.Err(); err != nil {
		return err
	}
	return flush()
}

// DictFile is the high-level CLI entry point: opens path (transparently
// gunzipping if needed), fills in the URI default when none was given,
// and writes the dictionary to w.
func DictFile(path string, w io.Writer, opts DictOptions) error {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return err
	}
	defer in.Close()
	if opts.URI == "" && path != "-" {
		if abs, err := filepath.Abs(path); err == nil {
			opts.URI = "file://" + abs
		}
	}
	return Dict(in, w, opts)
}
