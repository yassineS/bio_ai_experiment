// Package bedlinks implements `bedtools links`: it reads BED intervals and
// emits an HTML page of UCSC Genome Browser links, one row per record.
//
// Upstream reference:
//
//	reference_code/bedtools/src/linksBed/linksBed.cpp
//	reference_code/bedtools/src/linksBed/linksMain.cpp
//
// The output is a fixed HTML scaffold (head, intro paragraph, opening
// <table> tag, one <tr>/<td> block per record, closing </table>/</p>/
// </body>/</html>). The link target is built from -base/-org/-db; the
// link text uses 1-based start coordinates for display while the URL's
// `position=` query parameter uses 0-based start (matching upstream).
package bedlinks

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Defaults that mirror the upstream `links_main` defaults.
const (
	DefaultBase = "http://genome.ucsc.edu"
	DefaultOrg  = "human"
	DefaultDB   = "hg18"
)

// Options configures Run.
type Options struct {
	// Base is the UCSC mirror base URL (no trailing slash). Default:
	// "http://genome.ucsc.edu".
	Base string

	// Org is the UCSC organism token. Default: "human".
	Org string

	// DB is the UCSC build/db token. Default: "hg18".
	DB string

	// BedFile is the name shown in the HTML <title> tag. Upstream uses
	// the input filename verbatim (or "stdin" when reading from stdin);
	// the CLI wrapper passes that through.
	BedFile string
}

// Run reads BED records from r and writes the upstream-shaped HTML to w.
// Returns the number of <tr> data rows emitted.
func Run(r io.Reader, w io.Writer, opts Options) (int, error) {
	base := opts.Base
	if base == "" {
		base = DefaultBase
	}
	org := opts.Org
	if org == "" {
		org = DefaultOrg
	}
	db := opts.DB
	if db == "" {
		db = DefaultDB
	}
	bedFile := opts.BedFile
	if bedFile == "" {
		bedFile = "stdin"
	}

	// Mirror upstream's hgTracks URL prefix:
	//   <base>/cgi-bin/hgTracks?org=<org>&db=<db>&position=
	urlPrefix := base + "/cgi-bin/hgTracks?org=" + org + "&db=" + db + "&position="

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	// HTML header (verbatim from upstream).
	if _, err := fmt.Fprint(bw,
		"<html>\n",
		"\t<body>\n",
		"<title>", bedFile, "</title>\n",
		"<br>Firefox users: Press and hold the \"apple\" or \"alt\" key and click link to open in new tab.\n",
		"<p style=\"font-family:courier\">\n",
		"<table border=\"0\" align=\"justify\"\n",
		"<h3>BED Entries from: stdin </h3>\n",
	); err != nil {
		return 0, err
	}

	br := bed.NewReader(r)
	n := 0
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return n, err
		}
		if err := writeRow(bw, rec, urlPrefix); err != nil {
			return n, err
		}
		n++
	}

	// HTML footer (verbatim from upstream).
	if _, err := fmt.Fprint(bw,
		"</table>\n",
		"</p>\n",
		"\t</body>\n",
		"</html>\n",
	); err != nil {
		return n, err
	}
	return n, nil
}

// writeRow emits one <tr> block matching the upstream WriteURL switch
// over `_bed->bedType` (3/4/5/6/9/12). Upstream BED types are determined
// by the number of populated columns; we infer them from rec the same
// way bed.Reader does.
func writeRow(bw *bufio.Writer, rec *bed.Record, urlPrefix string) error {
	bedType := inferBedType(rec)

	// position=<chrom>:<start>-<end> in the URL (0-based start, matching
	// upstream's stringstream).
	position := fmt.Sprintf("%s:%d-%d", rec.Chrom, rec.ChromStart, rec.ChromEnd)

	if _, err := fmt.Fprintf(bw,
		"<tr>\n\t<td>\n\t\t<a href=%s%s>%s:%d-%d</a>\n\t</td>\n",
		urlPrefix, position, rec.Chrom, rec.ChromStart+1, rec.ChromEnd,
	); err != nil {
		return err
	}

	switch bedType {
	case 4:
		if _, err := fmt.Fprintf(bw, "\t<td>\n%s\n\t</td>\n", rec.Name); err != nil {
			return err
		}
	case 5:
		if _, err := fmt.Fprintf(bw, "\t<td>\n%s\n\t</td>\n\t<td>\n%d\n\t</td>\n",
			rec.Name, rec.Score); err != nil {
			return err
		}
	case 6, 9, 12:
		if _, err := fmt.Fprintf(bw,
			"\t<td>\n%s\n\t</td>\n\t<td>\n%d\n\t</td>\n\t<td>\n%s\n\t</td>\n",
			rec.Name, rec.Score, rec.Strand); err != nil {
			return err
		}
	}
	_, err := fmt.Fprint(bw, "</tr>\n")
	return err
}

// inferBedType returns the upstream `_bed->bedType` value: the number of
// populated columns in the record. Upstream's BedFile sets bedType from
// the FIRST data line; we mirror that per-record (bed.Reader populates
// optional fields in order, so a record with Strand populated implies
// columns 1..6 are all present).
func inferBedType(rec *bed.Record) int {
	switch {
	case len(rec.BlockStarts) > 0:
		return 12
	case rec.BlockCount != 0 || len(rec.BlockSizes) > 0:
		return 11
	case rec.ItemRGB != "":
		return 9
	case rec.ThickEnd != 0 || rec.ThickStart != 0:
		return 8
	case rec.Strand != "":
		return 6
	case rec.Score != 0:
		return 5
	case rec.Name != "":
		return 4
	default:
		return 3
	}
}
