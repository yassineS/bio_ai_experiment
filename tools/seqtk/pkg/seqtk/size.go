// size.go: implementation of `seqtk size`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_size (v1.5-r133).
// Behaviour: counts the number of records and the total sequence length
// across the input FASTA/FASTQ stream, then emits a single line:
//
//	<n_records>\t<total_bases>\n
//
// Upstream takes no flags. The only positional argument is the input
// file (or "-" for stdin); compression handling is done at the CLI
// layer via seqtk.OpenInput, exactly like the other ported subcommands.

package seqtk

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// Size reads every record in r (FASTA or FASTQ, auto-detected from the
// first non-whitespace byte) and writes a single "<n>\t<total_bases>\n"
// line to w, matching upstream "seqtk size" byte-for-byte.
func Size(r io.Reader, w io.Writer) error {
	br, isFastq := peekIsFastq(r)
	var n, total int64
	if isFastq {
		fr := fastq.NewReader(br, fastq.Phred33)
		for {
			rec, err := fr.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			n++
			total += int64(len(rec.Sequence))
		}
	} else {
		fr := fasta.NewReader(br)
		for {
			rec, err := fr.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			n++
			total += int64(len(rec.Sequence))
		}
	}
	_, err := fmt.Fprintf(w, "%d\t%d\n", n, total)
	return err
}
