/*
 * dverify.c — round-trip sanity check for dtrace output.
 *
 *     ./dverify --in <original> --compressed <foo.gz> [--format gzip|deflate|zlib]
 *
 * Decompresses <foo.gz> with libdeflate's decoder and asserts the result
 * is byte-equal to <original>.  This is the sanity check that our BTRACE
 * instrumentation doesn't perturb the encoder's output.
 */

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libdeflate.h"

/*
 * Match dtrace.c so the BTRACE_FP macro in btrace.h links cleanly.
 * dverify never enables tracing — it just round-trip-checks the
 * compressed output — but the encoder objects are shared between both
 * binaries and reference this symbol.
 */
FILE *btrace_fp = NULL;

static int slurp(const char *path, uint8_t **buf, size_t *len) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "dverify: open %s: %s\n", path, strerror(errno));
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
    if (n != (size_t)size) { free(*buf); return -1; }
    *len = n;
    return 0;
}

int main(int argc, char **argv) {
    const char *in_path = NULL;
    const char *cmp_path = NULL;
    const char *format = "gzip";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--in") && i + 1 < argc) in_path = argv[++i];
        else if (!strcmp(argv[i], "--compressed") && i + 1 < argc) cmp_path = argv[++i];
        else if (!strcmp(argv[i], "--format") && i + 1 < argc) format = argv[++i];
        else {
            fprintf(stderr, "dverify: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }
    if (!in_path || !cmp_path) {
        fprintf(stderr, "dverify: --in and --compressed are required\n");
        return 1;
    }

    uint8_t *in_buf = NULL, *cmp_buf = NULL;
    size_t in_len = 0, cmp_len = 0;
    if (slurp(in_path, &in_buf, &in_len) != 0) return 1;
    if (slurp(cmp_path, &cmp_buf, &cmp_len) != 0) { free(in_buf); return 1; }

    struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
    if (!d) {
        fprintf(stderr, "dverify: libdeflate_alloc_decompressor failed\n");
        free(in_buf); free(cmp_buf); return 1;
    }

    /* Decode into a buffer the exact size of the original; the decoder
     * will report any mismatch via the actual_out_nbytes_ret pointer. */
    uint8_t *out_buf = (uint8_t *)malloc(in_len == 0 ? 64 : in_len);
    if (!out_buf) { libdeflate_free_decompressor(d); free(in_buf); free(cmp_buf); return 1; }

    enum libdeflate_result res;
    size_t actual_out = 0;
    if (!strcmp(format, "deflate")) {
        res = libdeflate_deflate_decompress(d, cmp_buf, cmp_len, out_buf, in_len, &actual_out);
    } else if (!strcmp(format, "zlib")) {
        res = libdeflate_zlib_decompress(d, cmp_buf, cmp_len, out_buf, in_len, &actual_out);
    } else {
        res = libdeflate_gzip_decompress(d, cmp_buf, cmp_len, out_buf, in_len, &actual_out);
    }
    if (res != LIBDEFLATE_SUCCESS) {
        fprintf(stderr, "dverify: decompression failed (code %d)\n", (int)res);
        libdeflate_free_decompressor(d); free(in_buf); free(cmp_buf); free(out_buf);
        return 1;
    }
    if (actual_out != in_len) {
        fprintf(stderr, "dverify: length mismatch: got %zu, want %zu\n",
                actual_out, in_len);
        libdeflate_free_decompressor(d); free(in_buf); free(cmp_buf); free(out_buf);
        return 1;
    }
    if (memcmp(out_buf, in_buf, in_len) != 0) {
        fprintf(stderr, "dverify: byte mismatch on %s\n", in_path);
        libdeflate_free_decompressor(d); free(in_buf); free(cmp_buf); free(out_buf);
        return 1;
    }

    libdeflate_free_decompressor(d);
    free(in_buf); free(cmp_buf); free(out_buf);
    fprintf(stderr, "dverify: round-trip OK for %s (%zu bytes -> %zu bytes compressed)\n",
            in_path, in_len, cmp_len);
    return 0;
}
