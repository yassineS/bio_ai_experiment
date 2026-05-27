#!/usr/bin/env python3
"""
apply_btrace.py — produce instrumented copies of libdeflate sources.

Reads upstream libdeflate .c/.h files from --src-dir and writes copies to
--out-dir with BTRACE_*() macro invocations inserted at every internal
decision point in the deflate encoder.  Each transformation is anchored on
an exact substring match against the upstream source; if any anchor fails
to match, the script aborts with an error so we notice when libdeflate
upstream drifts.

Intended to be invoked from the Makefile:

    python3 patches/apply_btrace.py \\
        --src-dir reference_code/libdeflate \\
        --out-dir reference_code/libdeflate_trace/build/src

The output directory layout mirrors the input — `lib/deflate_compress.c`
becomes `<out>/lib/deflate_compress.c`, and so on.  Files that don't need
instrumentation are copied verbatim.

The goal is to keep all instrumentation in this script rather than carrying
a long patchfile.  When upstream rebases, this script's anchors are easier
to fix than a unified diff.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

# Files we copy verbatim (needed for compilation of the encoder + harness).
COPY_VERBATIM = [
    "common_defs.h",
    "libdeflate.h",
    "lib/adler32.c",
    "lib/crc32.c",
    "lib/crc32_multipliers.h",
    "lib/crc32_tables.h",
    "lib/cpu_features_common.h",
    "lib/decompress_template.h",
    "lib/deflate_compress.h",
    "lib/deflate_constants.h",
    "lib/deflate_decompress.c",
    "lib/gzip_compress.c",
    "lib/gzip_constants.h",
    "lib/gzip_decompress.c",
    "lib/lib_common.h",
    "lib/utils.c",
    "lib/zlib_compress.c",
    "lib/zlib_constants.h",
    "lib/zlib_decompress.c",
    "lib/bt_matchfinder.h",
    "lib/hc_matchfinder.h",
    "lib/ht_matchfinder.h",
    "lib/matchfinder_common.h",
    "lib/x86/adler32_impl.h",
    "lib/x86/adler32_template.h",
    "lib/x86/cpu_features.c",
    "lib/x86/cpu_features.h",
    "lib/x86/crc32_impl.h",
    "lib/x86/crc32_pclmul_template.h",
    "lib/x86/decompress_impl.h",
    "lib/x86/matchfinder_impl.h",
    "lib/arm/adler32_impl.h",
    "lib/arm/cpu_features.c",
    "lib/arm/cpu_features.h",
    "lib/arm/crc32_impl.h",
    "lib/arm/crc32_pmull_helpers.h",
    "lib/arm/crc32_pmull_wide.h",
    "lib/arm/matchfinder_impl.h",
    "lib/riscv/matchfinder_impl.h",
]


def _require_replace(text: str, old: str, new: str, label: str) -> str:
    """Replace the first occurrence of `old` with `new`; fail loudly if it
    isn't present (so upstream-rebase breakage is obvious)."""
    if old not in text:
        raise SystemExit(
            f"apply_btrace: anchor for {label!r} not found.  Has upstream "
            f"libdeflate drifted?\nLooked for:\n{old}"
        )
    if text.count(old) > 1:
        raise SystemExit(
            f"apply_btrace: anchor for {label!r} matched more than once; "
            f"refine the anchor.\nLooked for:\n{old}"
        )
    return text.replace(old, new, 1)


def instrument_deflate_compress(text: str) -> str:
    # 1) Inject btrace.h include after deflate_compress.h include.
    text = _require_replace(
        text,
        '#include "deflate_compress.h"\n#include "deflate_constants.h"\n',
        '#include "deflate_compress.h"\n#include "deflate_constants.h"\n'
        '#include "btrace.h"\n',
        label="include btrace.h",
    )

    # 2) Trace literal emissions in deflate_compress_lazy_generic.  There
    #    are three deflate_choose_literal call sites in the lazy parser; we
    #    instrument all of them.  Each anchor pins enough surrounding code
    #    to be unique.
    text = _require_replace(
        text,
        "\t\t\tif (cur_len < min_len ||\n"
        "\t\t\t    (cur_len == DEFLATE_MIN_MATCH_LEN &&\n"
        "\t\t\t     cur_offset > 8192)) {\n"
        "\t\t\t\t/* No match found.  Choose a literal. */\n"
        "\t\t\t\tdeflate_choose_literal(c, *in_next++, true,\n"
        "\t\t\t\t\t\t       seq);\n"
        "\t\t\t\tcontinue;\n"
        "\t\t\t}\n",
        "\t\t\tif (cur_len < min_len ||\n"
        "\t\t\t    (cur_len == DEFLATE_MIN_MATCH_LEN &&\n"
        "\t\t\t     cur_offset > 8192)) {\n"
        "\t\t\t\t/* No match found.  Choose a literal. */\n"
        "\t\t\t\tBTRACE_LIT((size_t)(in_next - in), *in_next);\n"
        "\t\t\t\tdeflate_choose_literal(c, *in_next++, true,\n"
        "\t\t\t\t\t\t       seq);\n"
        "\t\t\t\tcontinue;\n"
        "\t\t\t}\n",
        label="lazy: no-match literal",
    )

    text = _require_replace(
        text,
        "\t\t\t\t * becomes the current match.\n"
        "\t\t\t\t */\n"
        "\t\t\t\tdeflate_choose_literal(c, *(in_next - 2), true,\n"
        "\t\t\t\t\t\t       seq);\n",
        "\t\t\t\t * becomes the current match.\n"
        "\t\t\t\t */\n"
        "\t\t\t\tBTRACE_LIT((size_t)(in_next - 2 - in), *(in_next - 2));\n"
        "\t\t\t\tdeflate_choose_literal(c, *(in_next - 2), true,\n"
        "\t\t\t\t\t\t       seq);\n",
        label="lazy: better-next-match literal",
    )

    text = _require_replace(
        text,
        "\t\t\t\t\t * positions ahead, so use two literals.\n"
        "\t\t\t\t\t */\n"
        "\t\t\t\t\tdeflate_choose_literal(\n"
        "\t\t\t\t\t\tc, *(in_next - 3), true, seq);\n"
        "\t\t\t\t\tdeflate_choose_literal(\n"
        "\t\t\t\t\t\tc, *(in_next - 2), true, seq);\n",
        "\t\t\t\t\t * positions ahead, so use two literals.\n"
        "\t\t\t\t\t */\n"
        "\t\t\t\t\tBTRACE_LIT((size_t)(in_next - 3 - in), *(in_next - 3));\n"
        "\t\t\t\t\tdeflate_choose_literal(\n"
        "\t\t\t\t\t\tc, *(in_next - 3), true, seq);\n"
        "\t\t\t\t\tBTRACE_LIT((size_t)(in_next - 2 - in), *(in_next - 2));\n"
        "\t\t\t\t\tdeflate_choose_literal(\n"
        "\t\t\t\t\t\tc, *(in_next - 2), true, seq);\n",
        label="lazy2: two-literal lookahead",
    )

    # 3) Trace match emissions.  There are three deflate_choose_match call
    #    sites in lazy_generic (immediate "nice" match, lazy fallback,
    #    lazy2 fallback).  At each site, position is in_next - 1 - in (we
    #    advanced past the first byte already).
    text = _require_replace(
        text,
        "\t\t\tif (cur_len >= nice_len) {\n"
        "\t\t\t\tdeflate_choose_match(c, cur_len, cur_offset,\n"
        "\t\t\t\t\t\t     true, &seq);\n",
        "\t\t\tif (cur_len >= nice_len) {\n"
        "\t\t\t\tBTRACE_MATCH((size_t)(in_next - 1 - in), cur_len, cur_offset);\n"
        "\t\t\t\tdeflate_choose_match(c, cur_len, cur_offset,\n"
        "\t\t\t\t\t\t     true, &seq);\n",
        label="lazy: nice-len match",
    )

    text = _require_replace(
        text,
        "\t\t\t} else { /* !lazy2 */\n"
        "\t\t\t\t/*\n"
        "\t\t\t\t * No better match at the next position.  Output\n"
        "\t\t\t\t * the current match.\n"
        "\t\t\t\t */\n"
        "\t\t\t\tdeflate_choose_match(c, cur_len, cur_offset,\n"
        "\t\t\t\t\t\t     true, &seq);\n",
        "\t\t\t} else { /* !lazy2 */\n"
        "\t\t\t\t/*\n"
        "\t\t\t\t * No better match at the next position.  Output\n"
        "\t\t\t\t * the current match.\n"
        "\t\t\t\t */\n"
        "\t\t\t\tBTRACE_MATCH((size_t)(in_next - 2 - in), cur_len, cur_offset);\n"
        "\t\t\t\tdeflate_choose_match(c, cur_len, cur_offset,\n"
        "\t\t\t\t\t\t     true, &seq);\n",
        label="lazy: fallback match",
    )

    text = _require_replace(
        text,
        "\t\t\t\t * No better match at either of the next 2\n"
        "\t\t\t\t * positions.  Output the current match.\n"
        "\t\t\t\t */\n"
        "\t\t\t\tdeflate_choose_match(c, cur_len, cur_offset,\n"
        "\t\t\t\t\t\t     true, &seq);\n",
        "\t\t\t\t * No better match at either of the next 2\n"
        "\t\t\t\t * positions.  Output the current match.\n"
        "\t\t\t\t */\n"
        "\t\t\t\tBTRACE_MATCH((size_t)(in_next - 3 - in), cur_len, cur_offset);\n"
        "\t\t\t\tdeflate_choose_match(c, cur_len, cur_offset,\n"
        "\t\t\t\t\t\t     true, &seq);\n",
        label="lazy2: fallback match",
    )

    # 4) Trace block begin in lazy_generic.  The greedy path also has the
    #    same three-call pattern, so we anchor on the immediately preceding
    #    `next_recalc_min_len` declaration which is unique to the lazy
    #    parser.
    text = _require_replace(
        text,
        "\t\tconst u8 *next_recalc_min_len =\n"
        "\t\t\tin_next + MIN(in_end - in_next, 10000);\n"
        "\t\tstruct deflate_sequence *seq = c->p.g.sequences;\n"
        "\t\tu32 min_len;\n"
        "\n"
        "\t\tinit_block_split_stats(&c->split_stats);\n"
        "\t\tdeflate_begin_sequences(c, seq);\n"
        "\t\tmin_len = calculate_min_match_len(in_next,\n"
        "\t\t\t\t\t\t  in_max_block_end - in_next,\n"
        "\t\t\t\t\t\t  c->max_search_depth);\n",
        "\t\tconst u8 *next_recalc_min_len =\n"
        "\t\t\tin_next + MIN(in_end - in_next, 10000);\n"
        "\t\tstruct deflate_sequence *seq = c->p.g.sequences;\n"
        "\t\tu32 min_len;\n"
        "\n"
        "\t\tBTRACE_BLOCK_BEGIN((size_t)(in_next - in));\n"
        "\t\tinit_block_split_stats(&c->split_stats);\n"
        "\t\tdeflate_begin_sequences(c, seq);\n"
        "\t\tmin_len = calculate_min_match_len(in_next,\n"
        "\t\t\t\t\t\t  in_max_block_end - in_next,\n"
        "\t\t\t\t\t\t  c->max_search_depth);\n"
        "\t\tBTRACE_MIN_MATCH_LEN(\"initial\", min_len);\n",
        label="lazy: block begin + initial min_len",
    )

    text = _require_replace(
        text,
        "\t\t\tif (in_next >= next_recalc_min_len) {\n"
        "\t\t\t\tmin_len = recalculate_min_match_len(\n"
        "\t\t\t\t\t\t&c->freqs,\n"
        "\t\t\t\t\t\tc->max_search_depth);\n",
        "\t\t\tif (in_next >= next_recalc_min_len) {\n"
        "\t\t\t\tmin_len = recalculate_min_match_len(\n"
        "\t\t\t\t\t\t&c->freqs,\n"
        "\t\t\t\t\t\tc->max_search_depth);\n"
        "\t\t\t\tBTRACE_MIN_MATCH_LEN(\"recalc\", min_len);\n",
        label="lazy: recalc min_len",
    )

    # 5) Trace block header / cost decisions in deflate_flush_block.  We
    #    insert events before each block-type emission.
    text = _require_replace(
        text,
        "\tbest_cost = MIN(dynamic_cost, MIN(static_cost, uncompressed_cost));\n",
        "\tbest_cost = MIN(dynamic_cost, MIN(static_cost, uncompressed_cost));\n"
        "\tBTRACE_BLOCK_COSTS(dynamic_cost, static_cost, uncompressed_cost);\n",
        label="flush: block costs",
    )

    text = _require_replace(
        text,
        "\tif (best_cost == uncompressed_cost) {\n"
        "\t\t/*\n"
        "\t\t * Uncompressed block(s).  DEFLATE limits the length of\n",
        "\tif (best_cost == uncompressed_cost) {\n"
        "\t\tBTRACE_BLOCK_HEADER(0, block_length, is_final_block);\n"
        "\t\t/*\n"
        "\t\t * Uncompressed block(s).  DEFLATE limits the length of\n",
        label="flush: stored block header",
    )

    text = _require_replace(
        text,
        "\tif (best_cost == static_cost) {\n"
        "\t\t/* Static Huffman block */\n",
        "\tif (best_cost == static_cost) {\n"
        "\t\tBTRACE_BLOCK_HEADER(1, block_length, is_final_block);\n"
        "\t\t/* Static Huffman block */\n",
        label="flush: static block header",
    )

    text = _require_replace(
        text,
        "\t\t/* Dynamic Huffman block */\n"
        "\n"
        "\t\tcodes = &c->codes;\n",
        "\t\tBTRACE_BLOCK_HEADER(2, block_length, is_final_block);\n"
        "\t\t/* Dynamic Huffman block */\n"
        "\n"
        "\t\tcodes = &c->codes;\n",
        label="flush: dynamic block header",
    )

    # 6) Trace Huffman codeword lengths after deflate_make_huffman_codes.
    text = _require_replace(
        text,
        "static void\n"
        "deflate_make_huffman_codes(const struct deflate_freqs *freqs,\n"
        "\t\t\t   struct deflate_codes *codes)\n"
        "{\n"
        "\tdeflate_make_huffman_code(DEFLATE_NUM_LITLEN_SYMS,\n"
        "\t\t\t\t  MAX_LITLEN_CODEWORD_LEN,\n"
        "\t\t\t\t  freqs->litlen,\n"
        "\t\t\t\t  codes->lens.litlen,\n"
        "\t\t\t\t  codes->codewords.litlen);\n"
        "\n"
        "\tdeflate_make_huffman_code(DEFLATE_NUM_OFFSET_SYMS,\n"
        "\t\t\t\t  MAX_OFFSET_CODEWORD_LEN,\n"
        "\t\t\t\t  freqs->offset,\n"
        "\t\t\t\t  codes->lens.offset,\n"
        "\t\t\t\t  codes->codewords.offset);\n"
        "}\n",
        "static void\n"
        "deflate_make_huffman_codes(const struct deflate_freqs *freqs,\n"
        "\t\t\t   struct deflate_codes *codes)\n"
        "{\n"
        "\tunsigned _btrace_i;\n"
        "\tdeflate_make_huffman_code(DEFLATE_NUM_LITLEN_SYMS,\n"
        "\t\t\t\t  MAX_LITLEN_CODEWORD_LEN,\n"
        "\t\t\t\t  freqs->litlen,\n"
        "\t\t\t\t  codes->lens.litlen,\n"
        "\t\t\t\t  codes->codewords.litlen);\n"
        "\n"
        "\tdeflate_make_huffman_code(DEFLATE_NUM_OFFSET_SYMS,\n"
        "\t\t\t\t  MAX_OFFSET_CODEWORD_LEN,\n"
        "\t\t\t\t  freqs->offset,\n"
        "\t\t\t\t  codes->lens.offset,\n"
        "\t\t\t\t  codes->codewords.offset);\n"
        "\tfor (_btrace_i = 0; _btrace_i < DEFLATE_NUM_LITLEN_SYMS; _btrace_i++)\n"
        "\t\tif (codes->lens.litlen[_btrace_i])\n"
        "\t\t\tBTRACE_HUFFMAN_LEN(\"litlen\", _btrace_i, codes->lens.litlen[_btrace_i]);\n"
        "\tfor (_btrace_i = 0; _btrace_i < DEFLATE_NUM_OFFSET_SYMS; _btrace_i++)\n"
        "\t\tif (codes->lens.offset[_btrace_i])\n"
        "\t\t\tBTRACE_HUFFMAN_LEN(\"offset\", _btrace_i, codes->lens.offset[_btrace_i]);\n"
        "}\n",
        label="huffman codes: emit codelens",
    )

    # 7) Trace block-split heuristic firing.
    text = _require_replace(
        text,
        "\t\t/* Ready to end the block? */\n"
        "\t\tif (total_delta +\n"
        "\t\t    (block_length / 4096) * stats->num_observations >= cutoff)\n"
        "\t\t\treturn true;\n"
        "\t}\n"
        "\tmerge_new_observations(stats);\n"
        "\treturn false;\n"
        "}\n",
        "\t\t/* Ready to end the block? */\n"
        "\t\tif (total_delta +\n"
        "\t\t    (block_length / 4096) * stats->num_observations >= cutoff) {\n"
        "\t\t\tBTRACE_SPLIT_CHECK(block_length, total_delta, cutoff, 1);\n"
        "\t\t\treturn true;\n"
        "\t\t}\n"
        "\t\tBTRACE_SPLIT_CHECK(block_length, total_delta, cutoff, 0);\n"
        "\t}\n"
        "\tmerge_new_observations(stats);\n"
        "\treturn false;\n"
        "}\n",
        label="block split: check",
    )

    return text


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--src-dir", required=True, type=Path,
                   help="upstream libdeflate root")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="destination for instrumented sources")
    args = p.parse_args()

    src_root: Path = args.src_dir
    out_root: Path = args.out_dir
    out_root.mkdir(parents=True, exist_ok=True)

    # Copy verbatim files.
    for rel in COPY_VERBATIM:
        src = src_root / rel
        dst = out_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)

    # Instrument deflate_compress.c.
    src = src_root / "lib" / "deflate_compress.c"
    text = src.read_text()
    text = instrument_deflate_compress(text)
    dst = out_root / "lib" / "deflate_compress.c"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(text)

    print(f"apply_btrace: wrote instrumented sources to {out_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
