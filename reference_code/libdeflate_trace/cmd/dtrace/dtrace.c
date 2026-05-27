/*
 * dtrace.c — minimal CLI that exercises the instrumented libdeflate
 * encoder and dumps a structured trace alongside the compressed output.
 *
 *     ./dtrace --level 6 --in input.bin --out compressed.bin --trace trace.txt
 *
 * The compressor used is libdeflate's gzip variant (so the output is a
 * standalone .gz file that the htsgo decoder can round-trip).  When built
 * with -DBTRACE, every internal decision point in deflate_compress.c
 * writes a line to the file named by --trace.
 *
 * This is the slice-0 oracle: subsequent slices port the Go encoder
 * against the captured trace files.
 */

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libdeflate.h"

/*
 * BTRACE_FP target.  When the surrounding code is compiled without
 * -DBTRACE these symbols are still defined (the macros expand to no-ops
 * so they never read btrace_fp) but having the symbol present keeps the
 * link line uniform between trace and no-trace builds.
 */
FILE *btrace_fp = NULL;

struct args {
    int level;
    const char *in_path;
    const char *out_path;
    const char *trace_path;
    const char *format;  /* "gzip", "deflate", or "zlib"; default "gzip" */
};

static void usage(FILE *fp, const char *prog) {
    fprintf(fp,
        "usage: %s --in <input> --out <compressed> --trace <trace.txt> "
        "[--level N] [--format gzip|deflate|zlib]\n"
        "\n"
        "Slice-0 oracle harness for libdeflate-parity port.  Compresses\n"
        "<input> with libdeflate at the given level (default 6) and writes\n"
        "the compressed output to <compressed>.  When the harness is built\n"
        "with -DBTRACE the instrumented encoder emits one trace event per\n"
        "internal decision to <trace.txt>.\n", prog);
}

static int parse_args(int argc, char **argv, struct args *a) {
    a->level = 6;
    a->in_path = NULL;
    a->out_path = NULL;
    a->trace_path = NULL;
    a->format = "gzip";

    for (int i = 1; i < argc; i++) {
        const char *arg = argv[i];
        if (!strcmp(arg, "-h") || !strcmp(arg, "--help")) {
            usage(stdout, argv[0]);
            return 1;
        } else if (!strcmp(arg, "--level") && i + 1 < argc) {
            a->level = atoi(argv[++i]);
        } else if (!strcmp(arg, "--in") && i + 1 < argc) {
            a->in_path = argv[++i];
        } else if (!strcmp(arg, "--out") && i + 1 < argc) {
            a->out_path = argv[++i];
        } else if (!strcmp(arg, "--trace") && i + 1 < argc) {
            a->trace_path = argv[++i];
        } else if (!strcmp(arg, "--format") && i + 1 < argc) {
            a->format = argv[++i];
        } else {
            fprintf(stderr, "dtrace: unknown argument: %s\n", arg);
            usage(stderr, argv[0]);
            return -1;
        }
    }
    if (!a->in_path || !a->out_path) {
        fprintf(stderr, "dtrace: --in and --out are required\n");
        usage(stderr, argv[0]);
        return -1;
    }
    return 0;
}

static int slurp(const char *path, uint8_t **buf, size_t *len) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "dtrace: open %s: %s\n", path, strerror(errno));
        return -1;
    }
    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return -1; }
    long size = ftell(fp);
    if (size < 0) { fclose(fp); return -1; }
    rewind(fp);
    *buf = (uint8_t *)malloc((size_t)size + 1);
    if (!*buf) { fclose(fp); return -1; }
    size_t n = fread(*buf, 1, (size_t)size, fp);
    fclose(fp);
    if (n != (size_t)size) {
        fprintf(stderr, "dtrace: short read on %s\n", path);
        free(*buf);
        return -1;
    }
    *len = n;
    return 0;
}

int main(int argc, char **argv) {
    struct args a;
    int rc = parse_args(argc, argv, &a);
    if (rc != 0) return rc < 0 ? 1 : 0;

    if (a.trace_path) {
        btrace_fp = fopen(a.trace_path, "w");
        if (!btrace_fp) {
            fprintf(stderr, "dtrace: cannot open trace file %s: %s\n",
                    a.trace_path, strerror(errno));
            return 1;
        }
        /* Emit a header so the file is self-describing. */
        fprintf(btrace_fp,
                "BTRACE HEADER version=1 source=%s level=%d format=%s\n",
                a.in_path, a.level, a.format);
    }

    uint8_t *in_buf = NULL;
    size_t in_len = 0;
    if (slurp(a.in_path, &in_buf, &in_len) != 0) {
        if (btrace_fp) fclose(btrace_fp);
        return 1;
    }

    struct libdeflate_compressor *c = libdeflate_alloc_compressor(a.level);
    if (!c) {
        fprintf(stderr, "dtrace: libdeflate_alloc_compressor(%d) failed\n", a.level);
        free(in_buf);
        if (btrace_fp) fclose(btrace_fp);
        return 1;
    }

    size_t bound;
    if (!strcmp(a.format, "deflate")) {
        bound = libdeflate_deflate_compress_bound(c, in_len);
    } else if (!strcmp(a.format, "zlib")) {
        bound = libdeflate_zlib_compress_bound(c, in_len);
    } else {
        bound = libdeflate_gzip_compress_bound(c, in_len);
    }

    uint8_t *out_buf = (uint8_t *)malloc(bound == 0 ? 64 : bound);
    if (!out_buf) {
        fprintf(stderr, "dtrace: out-of-memory for %zu byte output buffer\n", bound);
        libdeflate_free_compressor(c);
        free(in_buf);
        if (btrace_fp) fclose(btrace_fp);
        return 1;
    }

    size_t out_len;
    if (!strcmp(a.format, "deflate")) {
        out_len = libdeflate_deflate_compress(c, in_buf, in_len, out_buf, bound);
    } else if (!strcmp(a.format, "zlib")) {
        out_len = libdeflate_zlib_compress(c, in_buf, in_len, out_buf, bound);
    } else {
        out_len = libdeflate_gzip_compress(c, in_buf, in_len, out_buf, bound);
    }
    if (out_len == 0) {
        fprintf(stderr, "dtrace: compression failed (output buffer too small)\n");
        free(out_buf);
        libdeflate_free_compressor(c);
        free(in_buf);
        if (btrace_fp) fclose(btrace_fp);
        return 1;
    }

    FILE *out_fp = fopen(a.out_path, "wb");
    if (!out_fp) {
        fprintf(stderr, "dtrace: open %s: %s\n", a.out_path, strerror(errno));
        free(out_buf);
        libdeflate_free_compressor(c);
        free(in_buf);
        if (btrace_fp) fclose(btrace_fp);
        return 1;
    }
    if (fwrite(out_buf, 1, out_len, out_fp) != out_len) {
        fprintf(stderr, "dtrace: short write on %s\n", a.out_path);
        fclose(out_fp);
        free(out_buf);
        libdeflate_free_compressor(c);
        free(in_buf);
        if (btrace_fp) fclose(btrace_fp);
        return 1;
    }
    fclose(out_fp);

    if (btrace_fp) {
        fprintf(btrace_fp,
                "BTRACE FOOTER in_bytes=%zu out_bytes=%zu\n",
                in_len, out_len);
        fclose(btrace_fp);
    }

    libdeflate_free_compressor(c);
    free(out_buf);
    free(in_buf);

    fprintf(stderr, "dtrace: %s (%zu bytes) -> %s (%zu bytes, level=%d, format=%s)\n",
            a.in_path, in_len, a.out_path, out_len, a.level, a.format);
    return 0;
}
