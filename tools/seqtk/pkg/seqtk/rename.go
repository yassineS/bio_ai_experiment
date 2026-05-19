// rename.go: implementation of `seqtk rename`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_rename (v1.5-r133).
// Behaviour: rewrites each input record's name to "<prefix><N>" where
// N is a 1-based counter that increments on every new "fragment" (a
// singleton, or the second member of a paired-end pair). Two adjacent
// records are considered mates -- and therefore share the same N --
// when their names compare equal under upstream's stk_rename rule:
// identical byte-length, identical bytes except possibly the final two
// bytes when both records end with a '/<digit>' suffix.
//
// Algorithm (1:1 port of stk_rename):
//
//	last := empty
//	n := 1
//	for each record r in stream:
//	    if last is set:
//	        if same_pair(last, r):
//	            emit_renamed(last, n)
//	            emit_renamed(r,    n)       // r emitted directly, not via last
//	            n++
//	            last := empty
//	        else:
//	            emit_renamed(last, n)
//	            n++
//	            cpy_kseq(last, r)           // see "Sticky comment" below
//	    else:
//	        cpy_kseq(last, r)               // first record
//	if last is set:
//	    emit_renamed(last, n)
//
// Sticky-comment quirk (upstream BUG reproduced verbatim for parity):
// upstream's cpy_kseq -> cpy_kstr early-returns when the source kstring
// is empty (seqtk.c:1210, "if (src->l == 0) return;"), so the previous
// record's comment is NOT cleared when the new record has no comment.
// In stk_rename this means a singleton (or pair-first-member) whose
// description has a comment will "leak" that comment into every
// subsequent emission of `last` until a record with a non-empty comment
// overwrites it. We track this here via stickyComment so the Go port
// produces the same byte stream as upstream "seqtk rename". The
// comment from a record's own slot (not via the leaked `last`) is
// always its real comment: in the paired-second-member branch upstream
// passes `seq` (not `&last`) to stk_printseq_renamed, so that record's
// own comment is used.
//
// Output format mirrors the input (FASTA -> FASTA, FASTQ -> FASTQ) with
// sequences/qualities emitted on a single un-wrapped line, matching
// upstream's stk_printseq_renamed(..., line_len=0) call sites. The
// optional trailing comment is rendered as " <comment>" iff non-empty.

package seqtk

import (
	"bufio"
	"io"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// Rename reads every record from r and writes them back to w with names
// rewritten to "<prefix><N>" per upstream "seqtk rename". The prefix may
// be empty (the upstream default when no prefix argument is given), in
// which case names become bare integers.
func Rename(r io.Reader, w io.Writer, prefix string) error {
	br, isFastq := peekIsFastq(r)
	if isFastq {
		return renameFastq(br, w, prefix)
	}
	return renameFasta(br, w, prefix)
}

// commentOf returns the part of a description after the first run of
// whitespace separating the ID from any trailing comment, or "" if
// there is none. Matches kseq's split into (name, comment).
func commentOf(description, id string) string {
	if len(description) <= len(id) {
		return ""
	}
	rest := description[len(id):]
	// Upstream kseq splits the header on the first whitespace; the
	// comment is whatever follows that single separator. We preserve
	// every byte past the leading whitespace verbatim so the comment
	// round-trips byte-identically.
	i := 0
	for i < len(rest) && (rest[i] == ' ' || rest[i] == '\t') {
		i++
	}
	return rest[i:]
}

// writeRenamedHeader writes "<sigil><prefix><n>[ <comment>]\n" matching
// upstream stk_printseq_renamed's header layout.
func writeRenamedHeader(bw *bufio.Writer, sigil byte, prefix string, n int64, comment string) error {
	if err := bw.WriteByte(sigil); err != nil {
		return err
	}
	if prefix != "" {
		if _, err := bw.WriteString(prefix); err != nil {
			return err
		}
	}
	if _, err := bw.WriteString(strconv.FormatInt(n, 10)); err != nil {
		return err
	}
	if comment != "" {
		if err := bw.WriteByte(' '); err != nil {
			return err
		}
		if _, err := bw.WriteString(comment); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}

// updateSticky mirrors upstream cpy_kstr(&dst->comment, &src->comment):
// if the new comment is empty, the previous sticky value is left alone
// (upstream's `if (src->l == 0) return;` early-return).
func updateSticky(sticky, newComment string) string {
	if newComment == "" {
		return sticky
	}
	return newComment
}

func renameFasta(in io.Reader, w io.Writer, prefix string) error {
	fr := fasta.NewReader(in)
	bw := bufio.NewWriter(w)
	emitLast := func(rec *fasta.Record, n int64, sticky string) error {
		if err := writeRenamedHeader(bw, '>', prefix, n, sticky); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Sequence); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	emitSecond := func(rec *fasta.Record, n int64) error {
		// Paired-second-member: upstream passes `seq` (not `&last`) to
		// stk_printseq_renamed, so the comment is the record's own.
		if err := writeRenamedHeader(bw, '>', prefix, n, commentOf(rec.Description, rec.ID)); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Sequence); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	var last *fasta.Record
	var sticky string // upstream-style sticky `last.comment` (see top of file)
	var n int64 = 1
	for {
		rec, err := fr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if last != nil {
			if samePairName(last.ID, rec.ID) {
				if err := emitLast(last, n, sticky); err != nil {
					return err
				}
				if err := emitSecond(rec, n); err != nil {
					return err
				}
				n++
				// Upstream sets `last.name.l = 0` only, leaving the
				// comment buffer dirty. We model that by clearing
				// `last` (so the next iteration takes the cpy_kseq
				// branch via the `if last == nil` arm) without
				// touching sticky.
				last = nil
				continue
			}
			if err := emitLast(last, n, sticky); err != nil {
				return err
			}
			n++
			last = rec
			sticky = updateSticky(sticky, commentOf(rec.Description, rec.ID))
			continue
		}
		last = rec
		sticky = updateSticky(sticky, commentOf(rec.Description, rec.ID))
	}
	if last != nil {
		if err := emitLast(last, n, sticky); err != nil {
			return err
		}
	}
	return bw.Flush()
}

func renameFastq(in io.Reader, w io.Writer, prefix string) error {
	fr := fastq.NewReader(in, fastq.Phred33)
	bw := bufio.NewWriter(w)
	emitLast := func(rec *fastq.Record, n int64, sticky string) error {
		if err := writeRenamedHeader(bw, '@', prefix, n, sticky); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Sequence); err != nil {
			return err
		}
		if _, err := bw.WriteString("\n+\n"); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Quality); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	emitSecond := func(rec *fastq.Record, n int64) error {
		if err := writeRenamedHeader(bw, '@', prefix, n, commentOf(rec.Description, rec.ID)); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Sequence); err != nil {
			return err
		}
		if _, err := bw.WriteString("\n+\n"); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Quality); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	var last *fastq.Record
	var sticky string
	var n int64 = 1
	for {
		rec, err := fr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if last != nil {
			if samePairName(last.ID, rec.ID) {
				if err := emitLast(last, n, sticky); err != nil {
					return err
				}
				if err := emitSecond(rec, n); err != nil {
					return err
				}
				n++
				last = nil
				continue
			}
			if err := emitLast(last, n, sticky); err != nil {
				return err
			}
			n++
			last = rec
			sticky = updateSticky(sticky, commentOf(rec.Description, rec.ID))
			continue
		}
		last = rec
		sticky = updateSticky(sticky, commentOf(rec.Description, rec.ID))
	}
	if last != nil {
		if err := emitLast(last, n, sticky); err != nil {
			return err
		}
	}
	return bw.Flush()
}
