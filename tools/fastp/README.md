# Fastp - Go Implementation

An all-in-one FASTQ preprocessor that combines quality filtering, adapter trimming, and various other preprocessing steps, reimplemented in Go.

## Features

- **Adapter Trimming**: Removes adapter sequences from 3' and 5' ends
- **Quality Filtering**: Filters reads based on quality scores
- **Length Filtering**: Filters reads by minimum and maximum length
- **N Content Filtering**: Removes reads with excessive N bases
- **Poly-tail Trimming**: Removes poly-G and poly-X tails (common in NovaSeq)
- **Sliding-window Quality Trimming**: `--cut_front` / `--cut_tail` / `--cut_right` style trimming
- **Complexity Filtering**: Filters low-complexity sequences
- **Built-in Gzip Support**: Automatically handles .gz compressed files
- **Memory Efficient**: Streaming processing for large files
- **Detailed Statistics**: Comprehensive preprocessing statistics

## Installation

```bash
cd tools/fastp
go build ./cmd/fastp
```

## Usage

### Basic Usage (Single-End)

```bash
fastp -i input.fastq -o output.fastq
```

### Paired-End Usage

```bash
fastp -I read1.fastq -O out1.fastq --in2 read2.fastq --out2 out2.fastq
```

### With Adapter Trimming

```bash
fastp -i input.fastq -o output.fastq -x AGATCGGAAGAGC
```

### Comprehensive Preprocessing

```bash
fastp -i input.fastq -o output.fastq \
  -x AGATCGGAAGAGC \
  -q 20 \
  -l 30 \
  --trim-poly-g \
  --max-n-count 3
```

### Sliding-window Quality Trimming

```bash
# Drop low-quality bases from the 5' and 3' ends using a sliding window
fastp -i input.fastq -o output.fastq --cut-front --cut-tail -W 4 -M 20

# Cut the read at (and after) the first low-quality window scanning 5'->3'
fastp -i input.fastq -o output.fastq --cut-right --cut-window-size 4 --cut-mean-quality 20
```

### With Gzip Files

```bash
fastp -i input.fastq.gz -o output.fastq.gz -x AGATCGGAAGAGC
```

### Auto-detect Adapter

```bash
# Automatically detect and trim common adapters
fastp -i input.fastq -o output.fastq --detect-adapter
```

### UMI Extraction

```bash
# Extract an 8-base UMI from the front of each read (single-end / R1).
fastp -i input.fastq -o output.fastq --umi --umi_loc read1 --umi_len 8
```

See the [UMI processing](#umi-processing) section below for all supported
`--umi_loc` modes (`read1`, `read2`, `per_read`, `index1`, `index2`,
`per_index`).

### Duplication Evaluation and Dedup

```bash
# Track duplication rate, build the duplication-level histogram, and drop
# duplicate reads from the output FASTQ.
fastp -i input.fastq -o output.fastq \
      --dup_calc_accuracy 3 --dedup \
      --json report.json --html report.html
```

See the [Duplication evaluation](#duplication-evaluation) section below for
details on the accuracy/memory tradeoff and how the results show up in
the JSON and HTML reports.

### Base Correction

```bash
# Correct low-quality bases to N
fastp -i input.fastq -o output.fastq --base-correction --correction-threshold 20
```

### Overlap-based Base Correction (Paired-End)

```bash
# Correct mismatched bases in the PE overlap using the higher-quality mate.
# Verbatim port of upstream's --correction; --base-correction above is a
# separate (legacy) SE feature that masks low-quality bases to N.
fastp -I R1.fastq -O out1.fastq --in2 R2.fastq --out2 out2.fastq --correction \
  --overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20
```

### Overrepresented Sequence Analysis

```bash
# Sample 1 in every N reads (default 20) and report overrepresented sequences
# under read{1,2}_before_filtering in the JSON report.
fastp -i input.fastq -o output.fastq -p -P 20 --json report.json
```

### Splitting Output Into Multiple Files

```bash
# Split into 4 files (0001.out.fq .. 0004.out.fq).
fastp -i input.fastq -o out.fq -s 4

# Split so each file holds at most 4000 lines (1000 reads), 3-digit prefix.
fastp -i input.fastq -o out.fq -S 4000 -d 3
```

### Merge Overlapping Paired-End Reads

```bash
# Merge overlapping pairs into single reads (upstream -m/--merge).
# Merged reads go to --merged_out; merge auto-enables base correction.
fastp --in1 R1.fastq --in2 R2.fastq --merge --merged_out merged.fastq

# Also keep unmerged/unpaired survivors in the merge stream
fastp --in1 R1.fastq --in2 R2.fastq --merge --include_unmerged \
  --merged_out merged.fastq

# Legacy heuristic merge (project extension, not upstream-faithful)
fastp -I R1.fastq -O out1.fastq --in2 R2.fastq --out2 out2.fastq --merge-overlap
```

### Trim by a FASTA List of Adapters / Disable Adapter Trimming

```bash
# Trim read1 (and read2 if PE) by every sequence in adapters.fa (>=6 bp)
fastp -i input.fastq -o output.fastq --adapter_fasta adapters.fa

# Disable all adapter trimming
fastp -i input.fastq -o output.fastq -A   # --disable_adapter_trimming
```

### Poly-X Tail Trimming With a Dedicated Length

```bash
# Poly-X uses its own --poly_x_min_len (independent of poly-G's knob)
fastp -i input.fastq -o output.fastq --trim-poly-x --poly_x_min_len 12
```

### Multi-threaded Processing with HTML Report

```bash
# Use 4 threads and generate HTML report
fastp -i input.fastq -o output.fastq -w 4 --html report.html
```

## Reports

`fastp` can emit two complementary report files in addition to the cleaned FASTQ:

- `--html FILE` writes a **self-contained HTML report** with embedded CSS and
  inline SVG charts (per-base quality, per-base composition, length
  distribution, summary tables, filtering reasons, adapter trimming). The file
  contains no JavaScript and pulls in no external resources, so it can be
  emailed or archived as-is.
- `--json FILE` writes a **JSON report** whose schema is intentionally close to
  upstream fastp's `fastp.json` (top-level `summary`, `filtering_result`,
  `duplication`, `adapter_cutting`, and per-read `read{1,2}_before_filtering` /
  `read{1,2}_after_filtering` sections). This makes it directly consumable by
  tools such as MultiQC.

Both report formats work when reading from stdin / writing to stdout because
all statistics are collected in memory while the file streams.

```bash
# Paired-end run with both reports plus overlap-based PE adapter detection.
fastp -I R1.fq.gz -O clean_R1.fq.gz --in2 R2.fq.gz --out2 clean_R2.fq.gz \
      --detect_adapter_for_pe \
      --html fastp_report.html --json fastp_report.json
```

### Note on the `-h` short flag

Upstream fastp uses `-h` for help. To preserve that muscle memory the
`--html` flag is **long-only** in this implementation (no `-h` alias);
likewise `--json` is long-only. `-h`/`--help` prints usage and exits.

## Options

### Input/Output

- `-i, --input FILE` - Input FASTQ file (single-end)
- `-o, --output FILE` - Output FASTQ file (single-end)
- `-I, --in1 FILE` - Input FASTQ file read 1 (paired-end)
- `--in2 FILE` - Input FASTQ file read 2 (paired-end)
- `-O, --out1 FILE` - Output FASTQ file read 1 (paired-end)
- `--out2 FILE` - Output FASTQ file read 2 (paired-end)

### Adapter Trimming

- `-x, --adapter3 SEQ` - 3' adapter sequence
- `-y, --adapter5 SEQ` - 5' adapter sequence

### Quality Filtering

- `-q, --qual-threshold INT` - Quality threshold (default: 15)
- `--qual-percent INT` - Percent of bases meeting quality (default: 40)
- `-Q, --disable_quality_filtering` - Disable quality filtering. Upstream gates
  the entire quality-filter block on `qualfilter.enabled`, so this skips the
  low-quality-base percentage limit, the average-quality requirement, **and the
  N-base-count limit** (the N check lives inside the quality block in upstream
  `filter.cpp`).

### Length Filtering

- `-l, --min-length INT` - Minimum read length (default: 15)
- `--max-length INT` - Maximum read length (0 = no limit)
- `-L, --disable_length_filtering` - Disable length filtering. Skips both the
  too-short (`length_required`) and too-long (`length_limit`) discard checks.

### Content Filtering

- `--max-n-count INT` - Maximum N count (default: 5)
- `--max-n-percent FLOAT` - Maximum N percentage (default: 20.0)

Note: the N-base content check is part of upstream's **quality** filter, so
`--disable_quality_filtering` (`-Q`) also turns it off.

### Poly-tail Trimming

- `-g, --trim_poly_g` - Enable poly-G tail trimming
- `-G, --disable_trim_poly_g` - Disable poly-G tail trimming (overrides
  `--trim_poly_g`)
- `--trim-poly-x` - Enable poly-X tail trimming
- `--poly-g-min-len INT` - Minimum poly-G length (default: 10)

### Sliding-window Quality Trimming

These options mirror upstream fastp's `--cut_front` / `--cut_tail` / `--cut_right`.
A window of `--cut-window-size` bases is slid along the read and its mean Phred
quality is compared against `--cut-mean-quality`. If `--cut-window-size` is larger
than the read, the whole read is treated as a single (short) window. These apply to
single-end reads and to both reads of a pair.

**Ordering (load-bearing for parity):** the sliding-window cut runs **first** —
before poly-G, adapter trimming, and poly-X — exactly as upstream fastp does
(`seprocessor.cpp` runs `trimAndCut` before the adapter block). Running it after
adapter trimming makes the window see an already-shortened read and shifts the
window-boundary math, which diverges from upstream by ~1bp on any read where
adapter/poly trimming fired. The port preserves the upstream order so `--cut_tail`
is byte-exact even when combined with adapter trimming (parity Case 17).

- `-5, --cut-front` - Slide a window from the 5' end; drop the leading base while
  the window's mean quality is below the threshold (equivalently: trim everything
  before the first 5'->3' window whose mean quality meets the threshold)
- `-3, --cut-tail` - Slide a window from the 3' end; drop the trailing base while
  the window's mean quality is below the threshold (equivalently: trim everything
  after the first 3'->5' window whose mean quality meets the threshold)
- `-r, --cut-right` - Slide a window 5'->3'; the moment a window's mean quality
  drops below the threshold, cut the read there and discard everything to its right
- `-W, --cut-window-size INT` - Window size for the above (default: 4)
- `-M, --cut-mean-quality INT` - Mean Phred-quality threshold for the window (default: 20)

When more than one of `--cut-front` / `--cut-tail` / `--cut-right` is given, they
are applied in that order (matching upstream fastp).

### Complexity Filtering

- `--low-complexity` - Enable complexity filtering
- `--complexity-threshold FLOAT` - Complexity threshold (default: 0.3)

### UMI Processing

UMI (Unique Molecular Identifier) extraction supports the same set of
locations as upstream fastp. The UMI is appended to the read name as
`:UMI_<prefix><umi>` so the downstream aligner can preserve molecular
identity.

- `--umi` - Enable UMI processing
- `--umi_loc STRING` - UMI location:
  - `read1` - UMI is the prefix of read 1 (default for single-end)
  - `read2` - UMI is the prefix of read 2
  - `per_read` - UMI prefix on BOTH R1 and R2; the two UMIs are joined
    with `_` (default for paired-end)
  - `index1` - UMI is the i7 index parsed from the Illumina header
  - `index2` - UMI is the i5 index parsed from the Illumina header
  - `per_index` - UMI is `i7_i5` (both index fields combined)
- `--umi_len INT` - UMI length in bases (used by `read1`/`read2`/`per_read`)
- `--umi_prefix STRING` - Optional prefix prepended to the UMI in the read
  name (default: empty)
- `--umi_skip INT` - Bases to skip immediately after the UMI bases
  (default: 0)

```bash
# Per-read UMI extraction on paired-end Illumina reads.
fastp -I R1.fq.gz -O clean_R1.fq.gz --in2 R2.fq.gz --out2 clean_R2.fq.gz \
      --umi --umi_loc per_read --umi_len 6 --umi_prefix L_
```

In `read1`/`read2`/`per_read` modes the UMI bases plus `--umi_skip`
trailing bases are removed from the sequence and quality strings.
In `index1`/`index2`/`per_index` modes the sequence is left untouched
and the UMI is parsed from the description line of the FASTQ record.

The legacy `--umi-length`, `--umi-location`, and `--umi-skip` flags are
still accepted as aliases.

### Duplication Evaluation

Duplication rate is approximated with a fixed-size hash table, matching
upstream fastp's algorithm. Memory usage is roughly
`2 * 2^(17 + accuracy)` bytes — about 1 MB at `--dup_calc_accuracy 3`.

- `--dup_calc_accuracy INT` - Accuracy bucket in `[1, 6]`. Higher buckets
  use more memory but produce fewer spurious hash collisions on highly
  diverse libraries. `0` (the default) disables duplication tracking.
- `--dedup` - When set, the second and later occurrences of each
  duplicate read are dropped from the output FASTQ stream. `--dedup`
  implies duplication tracking; it will use accuracy `3` if you don't
  specify one explicitly.

```bash
# Duplication report only — no reads are dropped, but the JSON report's
# "duplication" section will contain the rate plus a per-occurrence-count
# histogram, and the HTML report will include a Duplication section.
fastp -i input.fastq -o output.fastq --dup_calc_accuracy 3 \
      --json report.json --html report.html
```

The hash key is the first 16 bytes of each read (or the full sequence
when it's shorter than 16 bp). In paired-end mode the R1 sequence is
hashed.

## Examples

### NovaSeq Data Preprocessing

```bash
# Remove poly-G tails common in NovaSeq
fastp -i input.fastq -o output.fastq --trim-poly-g -q 20
```

### Strict Quality Control

```bash
# High-quality reads only
fastp -i input.fastq -o output.fastq \
  -q 25 \
  -l 50 \
  --qual-percent 90 \
  --max-n-count 0
```

### Paired-End Preprocessing

```bash
# Comprehensive paired-end preprocessing
fastp -I R1.fastq.gz -O clean_R1.fastq.gz \
      --in2 R2.fastq.gz --out2 clean_R2.fastq.gz \
      -x AGATCGGAAGAGC \
      -q 20 -l 30 \
      --trim-poly-g \
      --max-n-count 2
```

### Complete Preprocessing Pipeline

```bash
# All-in-one preprocessing
fastp -i raw.fastq.gz -o clean.fastq.gz \
  -x AGATCGGAAGAGC \
  -q 20 \
  -l 30 \
  --trim-poly-g \
  --max-n-count 2 \
  --low-complexity
```

## Statistics Output

```
Fastp Processing Statistics:
  Total reads:           10000
  Total bases:           1500000
  Clean reads:           8543 (85.43%)
  Clean bases:           1234567 (82.30%)
  Adapter trimmed:       7543 (75.43%)
  Adapter bases removed: 234567
  Poly-G trimmed:        2345 (23.45%)
  Poly-G bases removed:  23456
  Too short filtered:    892 (8.92%)
  Too many N filtered:   345 (3.45%)
```

## Comparison with Original Fastp

This is a simplified Go implementation focusing on core preprocessing functionality.

### Implemented Features

- ✅ Adapter trimming (3' and 5')
- ✅ **Automatic adapter detection**
- ✅ Quality filtering
- ✅ Length filtering
- ✅ N content filtering
- ✅ Poly-G/X tail trimming
- ✅ Sliding-window quality trimming (`--cut_front`/`--cut_tail`/`--cut_right`)
- ✅ Complexity filtering
- ✅ Built-in gzip support
- ✅ Paired-end read support
- ✅ **HTML report generation**
- ✅ **UMI/barcode processing** (`read1`, `read2`, `per_read`, `index1`,
  `index2`, `per_index`)
- ✅ **Base correction**
- ✅ **Overlap-based PE base correction** (`-c`/`--correction`) and overlap
  knobs (`--overlap_len_require`/`--overlap_diff_limit`/`--overlap_diff_percent_limit`)
- ✅ **Overrepresentation analysis** (`-p`/`--overrepresentation_analysis`,
  `-P`/`--overrepresentation_sampling`)
- ✅ **Output splitting** (`-s`/`--split`, `-S`/`--split_by_lines`,
  `-d`/`--split_prefix_digits`) — the per-file boundaries reproduce upstream's
  multi-threaded pack/thread distribution (pack i -> thread i%`-w`, each thread
  owning a strided set of split files, rolling its file at a pack boundary once
  its accumulated count reaches the per-file size), so every split file is
  byte-for-byte identical to upstream for any `-w` (validated for `-w 1` and
  `-w 2/3/4`, SE and PE, `--split N` and `--split_by_lines`). Note: upstream's
  split is **not** thread-count-invariant — `-w 1` and `-w 4` send different
  reads to a given numbered file (and `--split_by_lines` even yields a
  different file *count* per `-w`); the contract is byte-parity with upstream
  **per thread count**, which the port meets exactly.
- ✅ **Multi-threading support**
- ✅ **Duplication evaluation** (`--dup_calc_accuracy`) and dedup
  (`--dedup`)

### Implemented since the original draft of this list

- ✅ Overlap-analysis-driven **merge writer** (`-m`/`--merge`,
  `--merged_out`, `--include_unmerged`) — a verbatim port of upstream's
  `OverlapAnalysis::merge`; merged FASTQ output is byte-for-byte identical.
  (The legacy `--merge-overlap` heuristic is retained as a separate flag.)
- ✅ **`--adapter_fasta`** (trim against a FASTA list of adapters).
- ✅ Separate **`--poly_x_min_len`** knob (independent of poly-G).
- ✅ The explicit **`--disable_adapter_trimming`** (`-A`),
  **`--disable_quality_filtering`** (`-Q`), **`--disable_length_filtering`**
  (`-L`), and **`--disable_trim_poly_g`** (`-G`) toggles — the complete set of
  upstream `--disable_*` flags. Each gates the same block upstream gates
  (quality filtering also covers the N-base limit; length filtering covers both
  too-short and too-long); validated byte-exact in parity Case 20 and the
  binary-free `TestUnitDisable*` unit tests.

### Parity validation of the trimming algorithms

The three trimming algorithms that are most sensitive to upstream's exact
logic are validated against the upstream binary with the comparison
appropriate to each:

- **Poly-G / poly-X trimming** (`-g`/`-x`): verbatim port of
  `polyx.cpp` (1 mismatch per 8 bases scanned, capped at 5, anchored on
  the last poly base) — validated **byte-exact** (parity Case 12,
  `TestUnitTrimPolyG`/`TestUnitTrimPolyX`).
- **Sliding-window quality trimming** (`--cut_front`/`--cut_tail`/
  `--cut_right`): verbatim port of `filter.cpp::trimAndCut` (keeps the
  leading high-Q bases of the offending window, applies the `t = t-w+1`
  rewind on the passing cut_tail window, skips boundary `N`s) and applied in
  upstream's order (the cut runs **before** poly-G/adapter/poly-X) —
  validated **byte-exact** on the small fixtures (parity Cases 13/14,
  `TestUnitSlidingWindowCut`) **and at scale** on 5000 varying-length reads
  with low-Q tails, both with and without adapter trimming (parity Cases
  17/18), plus the binary-free `TestUnitCutTailBoundary`. The at-scale cases
  guard the window-boundary edge the parity pipeline surfaced (~1% of reads
  off by ~1bp) that the small fixtures missed; cut_front/cut_right are
  confirmed unregressed at scale by Case 19.
- **Adapter auto-detection** (`--detect_adapter` / `--detect_adapter_for_pe`):
  verbatim port of upstream's kmer + nucleotide-tree evaluator. Because
  detection is a sampling-dependent heuristic, it is validated by a
  documented **similarity bound** against the upstream binary (detected
  adapter equals upstream's within a few bp; per-read trimmed-length
  agreement >= 99% with no read off by > 2bp; aggregate metrics within
  1%) rather than strict byte-equality (parity Case 16,
  `TestUnitDetectAdapterSE`). In practice the recovered adapter and the
  trimmed output are byte-identical to upstream.

### JSON report parity

The `--json` report reproduces upstream's `fastp.json` structure
field-for-field (`jsonreporter.cpp` + `Stats::reportJson`), validated against
the upstream binary in `parity_json_fields_test.go`:

- Per-read blocks (`read{1,2}_{before,after}_filtering`) emit, for **both** the
  before- and after-filtering streams: `q20_bases`/`q30_bases`/`q40_bases`,
  `total_cycles`, `quality_curves` (per-base `A`/`T`/`C`/`G` + `mean`),
  `content_curves` (`A`/`T`/`C`/`G`/`N`/`GC`), and the 1024-entry `kmer_count`
  5-mer histogram. Integer counts and `kmer_count` are **exact**; the curves
  are deterministic but compared within a `1e-4` tolerance because upstream
  emits them through C++ 6-significant-digit float formatting.
- `summary.sequencing` (the deterministic "single/paired end (N cycles[ + N
  cycles])" descriptor) — exact.
- Paired-end `insert_size` block (`peak`, `unknown`, full 512-bin `histogram`),
  reproduced from upstream's overlap-analysis insert-size binning — exact.
- The **`merged_and_filtered`** block rename in merge mode (`-m`/`--merge`) is
  emitted (the after-filtering block is renamed and `read2_after_filtering`
  dropped, matching upstream).

**Intentionally excluded from the parity comparison (non-reproducible):**
`summary.fastp_version` (upstream's version string vs our `tool.version`), the
top-level `command` string, and our `tool.time` wall-clock field — none are
deterministic across runs/builds.

**Go extension:** the per-read blocks also carry a `length_distribution`
object that upstream's JSON omits (upstream surfaces that data only in the HTML
report); it is `omitempty` and does not affect the upstream-field comparison.

### Documented residuals

None. Both formerly-open residuals — the multi-thread `--split` file-boundary
distribution and the unemitted JSON report sub-fields — are now closed and
validated against the upstream binary (see above).

## Testing

```bash
go test ./pkg/fastp -v
go test ./pkg/fastp -cover
```

Test coverage: **>85%**

## Performance

The Go implementation provides good performance for most use cases:

| Operation | Dataset | Time | Notes |
|-----------|---------|------|-------|
| Basic filtering | 1M reads | ~2.5s | Quality + length filtering |
| With adapter trim | 1M reads | ~3.2s | All filters enabled |
| Poly-G trimming | 1M reads | ~2.8s | NovaSeq data |

## Use Cases

### Preprocessing for Alignment

```bash
# Clean reads before alignment
fastp -i raw.fastq -o clean.fastq \
  -x AGATCGGAAGAGC -q 20 -l 50
```

### NovaSeq Data Cleanup

```bash
# Remove poly-G artifacts from NovaSeq
fastp -i novaseq.fastq -o clean.fastq \
  --trim-poly-g -q 15
```

### Quality Control Pipeline

```bash
# Comprehensive QC
fastp -i raw.fastq -o qc.fastq \
  -x AGATCGGAAGAGC \
  --trim-poly-g \
  -q 25 -l 40 \
  --max-n-count 2 \
  --low-complexity
```

## Development Roadmap

### Version 1.0.0 (Current)

- ✅ All-in-one preprocessing
- ✅ Adapter trimming
- ✅ Quality and length filtering
- ✅ N content filtering
- ✅ Poly-tail trimming
- ✅ Complexity filtering
- ✅ Built-in gzip support
- ✅ Comprehensive tests (>85% coverage)
- ✅ Paired-end read support

### Version 1.1.0 (Completed)

- ✅ Automatic adapter detection
- ✅ UMI/barcode processing
- ✅ Base correction
- ✅ Overlap analysis for paired-end
- ✅ Multi-threading support
- ✅ HTML report generation

### Version 1.2.0 (Completed)

- ✅ Duplication evaluation (`--dup_calc_accuracy`)
- ✅ Read deduplication (`--dedup`)
- ✅ Extended UMI locations (`per_read`, `index1`, `index2`, `per_index`)

### Version 1.3.0 (Future)

- [ ] Per-tile quality filtering
- [ ] Advanced quality profiling
- [ ] Support for additional sequencing platforms

## License

Apache License 2.0 - See [LICENSE](../../LICENSE) for details.

## References

- Original fastp: <https://github.com/OpenGene/fastp>
- Paper: Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics.

## Authors

- Original fastp by Shifu Chen
- Go implementation by Bio AI Experiment Team
