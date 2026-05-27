/*
 * btrace.h — instrumentation hooks for libdeflate's deflate encoder.
 *
 * Slice-0 oracle harness for bio_ai_experiment's libdeflate-parity port.
 * When the surrounding C source is compiled with -DBTRACE, every internal
 * decision in deflate_compress.c emits a structured event to the BTRACE_FP
 * stream.  When -DBTRACE is not set, the macros expand to nothing and
 * libdeflate behaves identically to upstream.
 *
 * Event format: whitespace-separated tokens, one event per line.
 *
 *     BTRACE <EVENT> key=value key=value ...
 *
 * The first token after BTRACE is the event tag (LIT, MATCH, BLOCK_HEADER,
 * HUFFMAN_LEN, ...).  Subsequent tokens are KEY=VALUE pairs.  Numeric
 * values are decimal unless prefixed with 0x.  Subsequent slices parse the
 * trace from Go and compare against the port's emissions.
 */

#ifndef BIO_AI_BTRACE_H
#define BIO_AI_BTRACE_H

#ifdef BTRACE

#include <stdio.h>

/*
 * BTRACE_FP is the destination stream.  By default we write to stderr; the
 * dtrace harness rebinds it to a dedicated trace file before invoking
 * libdeflate so the compressed output stays on stdout/the --out file.
 */
extern FILE *btrace_fp;
#ifndef BTRACE_FP
#define BTRACE_FP (btrace_fp ? btrace_fp : stderr)
#endif

#define BTRACE_EMIT(...) do { \
	fprintf(BTRACE_FP, __VA_ARGS__); \
	fputc('\n', BTRACE_FP); \
} while (0)

/* LZ77 literal/match emissions from the lazy parser. */
#define BTRACE_LIT(pos, byte) \
	BTRACE_EMIT("BTRACE LIT pos=%zu byte=0x%02x", (size_t)(pos), (unsigned)(byte) & 0xFFu)
#define BTRACE_MATCH(pos, len, off) \
	BTRACE_EMIT("BTRACE MATCH pos=%zu len=%u offset=%u", (size_t)(pos), (unsigned)(len), (unsigned)(off))

/* Hash-chain matchfinder operations. */
#define BTRACE_MF_LONGEST(pos, min_len, max_len, max_depth, best_len, best_off) \
	BTRACE_EMIT("BTRACE MF_LONGEST pos=%zu min_len=%u max_len=%u max_depth=%u best_len=%u best_off=%u", \
		(size_t)(pos), (unsigned)(min_len), (unsigned)(max_len), \
		(unsigned)(max_depth), (unsigned)(best_len), (unsigned)(best_off))
#define BTRACE_MF_SKIP(pos, count) \
	BTRACE_EMIT("BTRACE MF_SKIP pos=%zu count=%u", (size_t)(pos), (unsigned)(count))

/* Block lifecycle. */
#define BTRACE_BLOCK_BEGIN(pos) \
	BTRACE_EMIT("BTRACE BLOCK_BEGIN pos=%zu", (size_t)(pos))
#define BTRACE_BLOCK_HEADER(type, size, last) \
	BTRACE_EMIT("BTRACE BLOCK_HEADER type=%d size=%u last=%d", \
		(int)(type), (unsigned)(size), (int)(last))
#define BTRACE_BLOCK_COSTS(dynamic, static_, uncompressed) \
	BTRACE_EMIT("BTRACE BLOCK_COSTS dynamic=%u static=%u uncompressed=%u", \
		(unsigned)(dynamic), (unsigned)(static_), (unsigned)(uncompressed))

/* Huffman code construction (one event per nonzero codelen). */
#define BTRACE_HUFFMAN_LEN(table, sym, codelen) \
	BTRACE_EMIT("BTRACE HUFFMAN_LEN table=%s sym=%u codelen=%u", \
		(table), (unsigned)(sym), (unsigned)(codelen))

/* Block-split heuristic firing. */
#define BTRACE_SPLIT_CHECK(block_length, total_delta, cutoff, decision) \
	BTRACE_EMIT("BTRACE SPLIT_CHECK block_length=%u total_delta=%u cutoff=%u decision=%d", \
		(unsigned)(block_length), (unsigned)(total_delta), \
		(unsigned)(cutoff), (int)(decision))

/* min-match-len decisions. */
#define BTRACE_MIN_MATCH_LEN(stage, min_len) \
	BTRACE_EMIT("BTRACE MIN_MATCH_LEN stage=%s min_len=%u", (stage), (unsigned)(min_len))

/*
 * BTRACE_BITS is the noisiest event (one per bitstream write).  It is
 * gated on a separate -DBTRACE_BITS so day-to-day oracle generation
 * doesn't produce gigabytes of trace data.  Slices 1+3 enable it only
 * for small fixtures.
 */
#ifdef BTRACE_BITS
#define BTRACE_BIT_WRITE(value, len) \
	BTRACE_EMIT("BTRACE BITS value=0x%x len=%u", (unsigned)(value), (unsigned)(len))
#else
#define BTRACE_BIT_WRITE(value, len) do { (void)(value); (void)(len); } while (0)
#endif

#else /* !BTRACE */

#define BTRACE_LIT(pos, byte)                                  do { } while (0)
#define BTRACE_MATCH(pos, len, off)                            do { } while (0)
#define BTRACE_MF_LONGEST(pos, min_len, max_len, max_depth, best_len, best_off) do { } while (0)
#define BTRACE_MF_SKIP(pos, count)                             do { } while (0)
#define BTRACE_BLOCK_BEGIN(pos)                                do { } while (0)
#define BTRACE_BLOCK_HEADER(type, size, last)                  do { } while (0)
#define BTRACE_BLOCK_COSTS(dynamic, static_, uncompressed)     do { } while (0)
#define BTRACE_HUFFMAN_LEN(table, sym, codelen)                do { } while (0)
#define BTRACE_SPLIT_CHECK(block_length, total_delta, cutoff, decision) do { } while (0)
#define BTRACE_MIN_MATCH_LEN(stage, min_len)                   do { } while (0)
#define BTRACE_BIT_WRITE(value, len)                           do { } while (0)

#endif /* BTRACE */

#endif /* BIO_AI_BTRACE_H */
