# reference_code/libdeflate_trace/Makefile
#
# Slice-0 oracle harness for the bio_ai_experiment libdeflate-parity port.
# Builds a small CLI ('dtrace') that wraps libdeflate's encoder and emits a
# structured trace of every internal decision when compiled with -DBTRACE.
#
# This is a test-time tool: it lives entirely under reference_code/ and
# does NOT affect the main Go module's dependency surface.  The Go side
# consumes only the artifacts produced by `make oracle`: the compressed
# binaries and the trace text files.
#
# Targets
# -------
#   build    : compile dtrace + dverify with BTRACE enabled
#   oracle   : (re)generate corpus inputs + oracle traces under
#              ../../pkg/htsgo/libdeflate/testdata/oracle/
#   test     : build + round-trip every corpus input through dverify
#   clean    : remove the build directory
#
# Toolchain: gcc or clang (set CC=).  No autotools, no cmake.

CC          ?= gcc
PYTHON      ?= python3

REPO_ROOT   := $(abspath ../..)
UPSTREAM    := $(REPO_ROOT)/reference_code/libdeflate
BUILD_DIR   := build
SRC_DIR     := $(BUILD_DIR)/src
OBJ_DIR     := $(BUILD_DIR)/obj
BIN_DIR     := $(BUILD_DIR)/bin
ORACLE_DIR  := $(REPO_ROOT)/pkg/htsgo/libdeflate/testdata/oracle
PATCHER     := patches/apply_btrace.py
BTRACE_HDR  := src/btrace.h

CFLAGS      ?= -O2 -Wall -Wno-unused-function -fPIC
CPPFLAGS    += -DBTRACE -I$(SRC_DIR) -I$(SRC_DIR)/lib -Isrc

# All .c files we compile.  Mirrors libdeflate's CMakeLists.txt minus the
# things we don't need (gzip program, tests, benchmarks).
LIB_SRCS := \
	lib/adler32.c \
	lib/crc32.c \
	lib/deflate_compress.c \
	lib/deflate_decompress.c \
	lib/gzip_compress.c \
	lib/gzip_decompress.c \
	lib/utils.c \
	lib/zlib_compress.c \
	lib/zlib_decompress.c \
	lib/x86/cpu_features.c \
	lib/arm/cpu_features.c

LIB_OBJS := $(patsubst lib/%.c,$(OBJ_DIR)/%.o,$(LIB_SRCS))

# Declare each generated .c as a prerequisite of .stamp, so make knows
# how to materialize them.
LIB_GENERATED_SRCS := $(patsubst %,$(SRC_DIR)/%,$(LIB_SRCS))

.PHONY: all build test oracle clean help

all: build

help:
	@echo "Targets: build | test | oracle | clean"
	@echo "  build  : compile dtrace + dverify (with -DBTRACE)"
	@echo "  oracle : regenerate corpus + traces under pkg/htsgo/libdeflate/testdata/oracle/"
	@echo "  test   : build, generate corpus, and round-trip-verify each input"

build: $(BIN_DIR)/dtrace $(BIN_DIR)/dverify

# Materialize the instrumented sources by copying upstream files into
# $(SRC_DIR) and patching deflate_compress.c.  The marker file lets make
# track whether we're up to date.
$(SRC_DIR)/.stamp: $(PATCHER) $(BTRACE_HDR) | $(UPSTREAM)/lib/deflate_compress.c
	@mkdir -p $(SRC_DIR)
	$(PYTHON) $(PATCHER) --src-dir $(UPSTREAM) --out-dir $(SRC_DIR)
	# Drop btrace.h alongside the patched source so #include "btrace.h"
	# resolves regardless of which directory is the C compiler's cwd.
	cp $(BTRACE_HDR) $(SRC_DIR)/btrace.h
	cp $(BTRACE_HDR) $(SRC_DIR)/lib/btrace.h
	touch $@

# The generated sources are produced by the stamp rule.
$(LIB_GENERATED_SRCS): $(SRC_DIR)/.stamp

$(OBJ_DIR)/%.o: $(SRC_DIR)/lib/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(BIN_DIR)/dtrace: cmd/dtrace/dtrace.c $(LIB_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LIB_OBJS)

$(BIN_DIR)/dverify: cmd/dverify/dverify.c $(LIB_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LIB_OBJS)

oracle: build
	@mkdir -p $(ORACLE_DIR)
	$(PYTHON) patches/make_corpus.py --out-dir $(ORACLE_DIR)
	@for name in empty single_byte repeated_a random_64k bgzf_payload; do \
	  echo "==> generating oracle trace for $$name"; \
	  $(BIN_DIR)/dtrace --level 6 --format gzip \
	    --in    $(ORACLE_DIR)/$$name.bin \
	    --out   $(ORACLE_DIR)/$$name.gz \
	    --trace $(ORACLE_DIR)/$$name.trace; \
	done
	@echo "oracle: corpus + traces written to $(ORACLE_DIR)"

test: oracle
	@for name in empty single_byte repeated_a random_64k bgzf_payload; do \
	  echo "==> round-trip-verifying $$name"; \
	  $(BIN_DIR)/dverify \
	    --in         $(ORACLE_DIR)/$$name.bin \
	    --compressed $(ORACLE_DIR)/$$name.gz \
	    --format     gzip; \
	done
	@echo "test: all corpus inputs round-trip cleanly"

clean:
	rm -rf $(BUILD_DIR)
