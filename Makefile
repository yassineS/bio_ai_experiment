# Makefile — build our Go re-implementations (and, best-effort, the vendored
# upstream oracle binaries) into bin/.
#
# bin/ is gitignored: nothing built here is ever committed.
#
#   bin/ours/      our Go tools + the realbench/realparity harnesses
#   bin/upstream/  the upstream reference binaries built from reference_code/
#
# Run `make help` for the target list. POSIX-make portable (no GNU-only
# constructs): the ours target shells out to a portable for-loop so it does not
# depend on `wildcard`/`foreach`.

# Go binary; override with `make GO=/path/to/go`.
GO ?= go

# Output layout.
BIN          := bin
OURS_DIR     := $(BIN)/ours
UPSTREAM_DIR := $(BIN)/upstream

.PHONY: all ours upstream build clean help

all: ours ## Default: build our Go binaries into bin/ours/

ours: ## Build every tools/*/cmd/* plus realbench/realparity into bin/ours/
	@mkdir -p $(OURS_DIR)
	@for d in tools/*/cmd/*; do \
		name=`basename "$$d"`; \
		echo "building ours/$$name"; \
		$(GO) build -o "$(OURS_DIR)/$$name" "./$$d" || exit 1; \
	done
	@echo "building ours/realbench"
	@$(GO) build -o "$(OURS_DIR)/realbench" ./pipeline/cmd/realbench
	@echo "building ours/realparity"
	@$(GO) build -o "$(OURS_DIR)/realparity" ./pipeline/cmd/realparity

upstream: ## Best-effort build the vendored upstream binaries into bin/upstream/
	@mkdir -p $(UPSTREAM_DIR)
	@echo "note: building upstream from reference_code/ submodules (best effort)."
	@echo "note: fastp/vcftools/prinseq/mosdepth need extra system deps and are"
	@echo "      NOT built here; the canonical full build is test/nextflow/Dockerfile."
	-cd reference_code/htslib && autoreconf -i 2>/dev/null; ./configure && $(MAKE) && cp bgzip tabix htsfile "$(CURDIR)/$(UPSTREAM_DIR)/"
	-cd reference_code/samtools && autoheader 2>/dev/null; autoconf -Wno-syntax 2>/dev/null; ./configure --with-htslib=../htslib && $(MAKE) && cp samtools "$(CURDIR)/$(UPSTREAM_DIR)/"
	-cd reference_code/bcftools && autoheader 2>/dev/null; autoconf -Wno-syntax 2>/dev/null; ./configure --with-htslib=../htslib && $(MAKE) && cp bcftools "$(CURDIR)/$(UPSTREAM_DIR)/"
	-cd reference_code/bedtools && $(MAKE) && cp bin/bedtools "$(CURDIR)/$(UPSTREAM_DIR)/"
	-cd reference_code/seqtk && $(MAKE) && cp seqtk "$(CURDIR)/$(UPSTREAM_DIR)/"
	-cd reference_code/sickle && $(MAKE) && cp sickle "$(CURDIR)/$(UPSTREAM_DIR)/"
	-cd reference_code/skewer && $(MAKE) && cp skewer "$(CURDIR)/$(UPSTREAM_DIR)/"

build: ours upstream ## Build both ours and upstream

clean: ## Remove the bin/ tree
	rm -rf $(BIN)

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-12s %s\n", $$1, $$2}'
