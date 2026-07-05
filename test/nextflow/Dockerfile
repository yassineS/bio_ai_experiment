# syntax=docker/dockerfile:1
#
# realbench container — ONE Linux image carrying BOTH sides of the benchmark on
# the same architecture:
#
#   /opt/ours      our Go re-implementations (tools/<tool>/cmd/<tool>)
#   /opt/upstream  the upstream oracle binaries built from reference_code/
#   /usr/local/bin/realbench   the Go matrix harness (pipeline/cmd/realbench)
#   awscli, jq, bgzip, tabix on PATH for staging + aggregation
#
# Build from the REPO ROOT so the build context includes tools/, pkg/,
# pipeline/, and the reference_code/ submodules:
#
#   git submodule update --init reference_code/htslib reference_code/samtools \
#       reference_code/bcftools reference_code/bedtools reference_code/seqtk \
#       reference_code/sickle reference_code/skewer reference_code/fastp \
#       reference_code/vcftools reference_code/prinseq reference_code/patches
#   docker build -f Dockerfile -t ghcr.io/yassines/bio-ai-realbench:latest .
#
# The image is intentionally single-arch: build it on the same arch you will run
# the benchmark on (linux/amd64 for AWS Batch) so ours and upstream are compared
# on identical hardware. Pass --platform=linux/amd64 if cross-building.

ARG GO_VERSION=1.24.9
ARG DEBIAN=bookworm

# =============================================================================
#  Stage 1: build OUR Go binaries + the realbench harness.
# =============================================================================
FROM golang:${GO_VERSION}-${DEBIAN} AS gobuild

WORKDIR /src
# Prime the module cache.
COPY go.mod go.sum ./
RUN --mount=type=cache,target=/root/.cache/go-build \
    --mount=type=cache,target=/go/pkg/mod \
    go mod download

# Copy the source we actually compile (tools + shared pkg + pipeline).
COPY tools/ tools/
COPY pkg/ pkg/
COPY pipeline/ pipeline/

# Build every tool binary into /opt/ours. tools/<tool>/cmd/<tool>/main.go is the
# entry point; iterate the cmd dirs and name the binary after the leaf dir.
RUN --mount=type=cache,target=/root/.cache/go-build \
    --mount=type=cache,target=/go/pkg/mod \
    mkdir -p /opt/ours && \
    for d in tools/*/cmd/*; do \
        name=$(basename "$d"); \
        echo "building our $name"; \
        go build -trimpath -o "/opt/ours/$name" "./$d"; \
    done

# Build the realbench matrix harness.
RUN --mount=type=cache,target=/root/.cache/go-build \
    --mount=type=cache,target=/go/pkg/mod \
    go build -trimpath -o /usr/local/bin/realbench ./pipeline/cmd/realbench

# =============================================================================
#  Stage 2: build the UPSTREAM oracle binaries from the submodules.
# =============================================================================
FROM debian:${DEBIAN}-slim AS upstream

# Build deps for htslib/samtools/bcftools/bedtools/vcftools/fastp/seqtk/sickle/skewer.
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        autoconf automake \
        zlib1g-dev libbz2-dev liblzma-dev \
        libdeflate-dev libcurl4-openssl-dev libssl-dev \
        libncurses-dev \
        libisal-dev \
        pkg-config \
        python3 \
        perl \
        ca-certificates \
        git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
# The reference_code/ submodules (already checked out on the host) are the
# upstream sources. Copy them in.
COPY reference_code/ reference_code/

RUN mkdir -p /opt/upstream

# --- htslib (provides bgzip / tabix / htsfile, and is a dep of sam/bcftools) --
RUN cd reference_code/htslib && \
    autoreconf -i 2>/dev/null || true && \
    ./configure && \
    make -j"$(nproc)" && \
    cp bgzip tabix htsfile /opt/upstream/

# --- samtools ----------------------------------------------------------------
RUN cd reference_code/samtools && \
    autoheader 2>/dev/null || true && \
    autoconf -Wno-syntax 2>/dev/null || true && \
    ./configure --with-htslib=../htslib && \
    make -j"$(nproc)" && \
    cp samtools /opt/upstream/

# --- bcftools ----------------------------------------------------------------
RUN cd reference_code/bcftools && \
    autoheader 2>/dev/null || true && \
    autoconf -Wno-syntax 2>/dev/null || true && \
    ./configure --with-htslib=../htslib && \
    make -j"$(nproc)" && \
    cp bcftools /opt/upstream/

# The "leaf" QC/format tools below are BEST-EFFORT: several are old C/C++ that
# may not compile with a modern (gcc-12) toolchain — e.g. skewer trips a
# std::map const-comparator static_assert. A failure here must NOT abort the
# image: the tool is simply absent from /opt/upstream and realbench SKIPs its
# cells. (htslib/samtools/bcftools above stay strict — they are the core.)

# --- bedtools ----------------------------------------------------------------
RUN cd reference_code/bedtools && make -j"$(nproc)" && cp bin/bedtools /opt/upstream/ \
    || echo "WARN: bedtools upstream build failed; realbench bed* cells will SKIP"

# --- seqtk -------------------------------------------------------------------
RUN cd reference_code/seqtk && make && cp seqtk /opt/upstream/ \
    || echo "WARN: seqtk upstream build failed; realbench seqtk cells will SKIP"

# --- sickle ------------------------------------------------------------------
RUN cd reference_code/sickle && make && cp sickle /opt/upstream/ \
    || echo "WARN: sickle upstream build failed; realbench sickle cells will SKIP"

# --- skewer (old C++; try a permissive std, else skip) -----------------------
RUN cd reference_code/skewer && \
    { make || make CXXFLAGS='-c -O2 -std=gnu++98 -fpermissive'; } && cp skewer /opt/upstream/ \
    || echo "WARN: skewer upstream build failed; realbench skewer cells will SKIP"

# --- fastp (needs libisal + libdeflate) --------------------------------------
RUN cd reference_code/fastp && make -j"$(nproc)" && cp fastp /opt/upstream/ \
    || echo "WARN: fastp upstream build failed; realbench fastp cells will SKIP"

# --- vcftools (autotools; apply the vendored tmpfile patch first) ------------
RUN cd reference_code/vcftools && \
    (git apply ../patches/vcftools-tmpfile-vla-off-by-one.patch 2>/dev/null || \
     patch -p1 < ../patches/vcftools-tmpfile-vla-off-by-one.patch 2>/dev/null || \
     echo "vcftools patch already applied or not needed") && \
    ./autogen.sh && ./configure && make -j"$(nproc)" && cp src/cpp/vcftools /opt/upstream/ \
    || echo "WARN: vcftools upstream build failed; realbench vcftools cells will SKIP"

# --- prinseq (Perl script, not a compiled binary) ----------------------------
# Installed as an executable wrapper so /opt/upstream/prinseq runs it via perl.
RUN cp reference_code/prinseq/prinseq-lite.pl /opt/upstream/prinseq-lite.pl && \
    printf '#!/bin/sh\nexec perl /opt/upstream/prinseq-lite.pl "$@"\n' > /opt/upstream/prinseq && \
    chmod +x /opt/upstream/prinseq

# NOTE: mosdepth is a Nim project distributed only as a linux/amd64 GitHub
# release asset (not built from source). The runtime stage below pulls the
# pinned release binary from the mosdepth repo into /opt/upstream/mosdepth so
# the mosdepth realbench cells run a real comparison. (The harness can still
# override it via the MOSDEPTH_BIN env var / release cache — see
# pipeline/internal/upstream.)

# =============================================================================
#  Stage 3: the runtime image.
# =============================================================================
FROM debian:${DEBIAN}-slim AS runtime

# Runtime shared libs the upstream binaries link against, plus the staging/
# aggregation tooling (awscli, jq, perl for prinseq, wget/curl for https/ftp).
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        zlib1g libbz2-1.0 liblzma5 \
        libdeflate0 libcurl4 libssl3 \
        libncurses6 \
        libisal2 \
        perl \
        awscli \
        jq \
        wget curl ca-certificates \
        gawk \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Our binaries, upstream binaries, and the harness.
COPY --from=gobuild   /opt/ours              /opt/ours
COPY --from=gobuild   /usr/local/bin/realbench /usr/local/bin/realbench
COPY --from=upstream  /opt/upstream          /opt/upstream

# Upstream mosdepth — a Nim project distributed only as a linux/amd64 GitHub
# release asset, so it is pulled from the mosdepth repo's release (not built
# from source). Baking it into /opt/upstream lets the mosdepth realbench cells
# run a real ours-vs-upstream comparison instead of SKIPping. Version pinned to
# match the harness's expected oracle (pipeline/internal/upstream MOSDEPTH cache
# name mosdepth_v0.3.14).
ARG MOSDEPTH_VERSION=v0.3.14
RUN curl -fsSL -o /opt/upstream/mosdepth \
        "https://github.com/brentp/mosdepth/releases/download/${MOSDEPTH_VERSION}/mosdepth" \
    && chmod +x /opt/upstream/mosdepth \
    && /opt/upstream/mosdepth --help >/dev/null 2>&1

# Put bgzip/tabix on PATH (Nextflow STAGE_* steps call them via $upstream_bin,
# but having them on PATH is convenient for interactive debugging).
RUN ln -sf /opt/upstream/bgzip /usr/local/bin/bgzip && \
    ln -sf /opt/upstream/tabix /usr/local/bin/tabix && \
    ln -sf /opt/upstream/samtools /usr/local/bin/samtools && \
    ln -sf /opt/upstream/bcftools /usr/local/bin/bcftools

# CRAM reference-free decode: REF_PATH=/dev/null stops samtools/htslib trying to
# fetch reference sequences over the network from the EBI ENA cache.
ENV REF_PATH=/dev/null \
    PATH="/opt/ours:/opt/upstream:${PATH}"

# Sanity smoke at build time. `realbench` on PATH and an executable
# /opt/upstream/mosdepth are HARD assertions: a missing binary must fail the
# build here rather than surface at runtime as an exit-127 "command not found"
# on AWS Batch. The -h/--version invocations are informational (tolerate their
# exit codes).
RUN set -e; \
    /opt/upstream/samtools --version | head -1; \
    /opt/ours/samtools --version 2>/dev/null | head -1 || true; \
    command -v realbench >/dev/null || { echo "FATAL: realbench not on PATH"; exit 1; }; \
    realbench -h 2>&1 | head -1 || true; \
    test -x /opt/upstream/mosdepth || { echo "FATAL: /opt/upstream/mosdepth missing or not executable"; exit 1; }; \
    /opt/upstream/mosdepth --version 2>/dev/null | head -1 || true

WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
