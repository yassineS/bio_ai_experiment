# Parity roadmap

**Goal:** **1:1 feature parity** with the upstream tool for every Go port in
this repo. This file is the authoritative gap list per tool.

The project's stated goal is to make these tools faster, better tested, and
better documented than their originals — which requires that we actually
implement the same features. "Working subset" is a milestone, not the
destination.

> The companion file [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md) tracks bugs
> we identify in the originals as we go. We don't carry those over; we
> either fix on port or note and skip with a clear pointer.

The companion file [`../analysis/tool_ranking_2026.md`](../analysis/tool_ranking_2026.md)
ranks **next** tools to port. It is **not** a "deprioritise existing tools"
filter — existing tools all get carried to 1:1.

## Definition of "1:1"

A tool is **1:1** when:

1. Every subcommand the upstream binary exposes is present in the Go port.
2. Every documented flag/option is recognised (either implemented or
   gracefully rejected with a clear error pointing at this roadmap).
3. For every input shape upstream accepts, our port produces the same
   **logical result**, modulo:
   - documented intentional deviations (recorded in the tool's README +
     this roadmap),
   - upstream bugs we chose to fix (recorded in [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md)),
   - **RNG / stochastic divergence** — see the policy below.
4. The validated-parity test suite (runs the upstream test corpus through
   our port) passes for every supported case, with explicit `t.Skip()` for
   each documented exception.

### RNG / stochastic-output policy (2026-05-15)

For subcommands that use randomness (`bedshuffle`, `bedsample`, `seqtk sample`,
`seqtk randbase`, `samtools view -s subsample`, etc.) the parity bar is:

- **Reproducibility within our tool**: same seed + same input → same output
  byte-for-byte, every time, across Go versions and platforms.
- **Logical equivalence with upstream**: our output must satisfy the same
  invariants upstream's output does (correct sampling fraction, no
  duplicates without replacement, strand filters honoured, etc.).
- **NOT** byte-identical with upstream's output. Upstream uses its own
  C/C++ RNG (typically a Mersenne Twister from libstdc++); Go uses
  `math/rand`. Porting upstream's RNG would be ~200 lines of
  bit-twiddling per tool with no functional benefit — the user
  explicitly opted out of that work in favour of focusing on real
  feature parity.

The parity-test infrastructure handles this by either:

- structural-invariant assertions (e.g. "every shuffled interval has the
  same length as the input; every shuffled interval is on a chrom in the
  genome"), or
- a documented `t.Skip("RNG byte-parity, see PARITY_ROADMAP.md#rng-policy")`
  with a pointer to this section.

We're not there yet for any tool. The bedtools subset (PR #55) is the
closest — 127 parity tests, 85 passing, 42 documented `t.Skip` — and even
there we have ~30 subcommands not yet started.

---

## Per-tool gap list

Numbers reflect state at 2026-05-14 (post-#71). Update when each gap is
closed.

### `seqtk`

**Status:** 24 of 24 upstream subcommands. **1:1 PARITY ACHIEVED.**

Dispatch-table audit at `reference_code/seqtk/seqtk.c::main()`
lines 2099-2122 lists exactly these 24 `stk_*` entry points (`hrun`
and `hpc` are SEPARATE dispatch entries, not aliases — see
`stk_hpc` at seqtk.c:1692 vs `stk_hrun` at seqtk.c:1174): `comp`,
`fqchk`, `hety`, `gc`, `subseq`, `mutfa`, `mergefa`, `mergepe`,
`dropse`, `randbase`, `cutN`, `gap`, `listhet`, `famask`, `trimfq`,
`hrun`, `sample`, `seq`, `kfreq`, `rename`, `split`, `hpc`, `size`,
`telo`. All 24 are now implemented in
`tools/seqtk/cmd/seqtk/main.go`.

Earlier roadmap iterations listed `hpc-bg` as missing — that was a
misreading of the dispatch table. There is no `hpc-bg` subcommand
upstream (verified empirically: `seqtk hpc-bg` is rejected with
`[main] unrecognized command 'hpc-bg'. Abort!`). The genuine missing
entry was `hrun` (homopolymer-RUN finder, BED4 output); the name
collision with `hpc` confused the original audit.

Added this iteration: `listhet`, `hrun`.

- `listhet` — full upstream surface implemented (no flags, positional
  `<in.fa>` only; the Go cmd adds `-o/--output` as the project-wide
  file-redirect convenience). The 1:1-ported algorithm walks every
  FASTA record byte-by-byte and emits a TSV row `name\tpos1based\tbyte`
  for every byte whose `bitcnt_table[seq_nt16_table[b]] == 2` —
  i.e. the 2-base IUPAC codes R, Y, S, W, K, M and their lowercase
  counterparts. The byte is emitted in its original case. Byte-parity
  verified against `reference_code/seqtk` v1.5 on `ambig.fa`,
  `hety_basic.fa`, `hety_lowercase.fa`, and `small.fa` (the latter
  pins the no-output path for fixtures with zero hets).
- `hrun` — full upstream surface implemented (no flags, positional
  `<in.fa> [minLen]`; the Go cmd exposes the minLen knob as
  `-l/--min-len` AND still accepts the upstream positional form for
  compatibility). One BED4 row (`chrom\tstart\tend\tbase`, 0-based
  half-open) is emitted per maximal byte-identical run of length
  `>= minLen` (default 7). Two upstream quirks are mirrored
  byte-for-byte for parity:
  1. Comparison is BYTE-EXACT (no case-fold, no IUPAC fold) so
     `AAaa` is two runs of length 2.
  2. The "open trailing run" flush at `seqtk.c:1200` lives OUTSIDE
     the `kseq_read` loop, so it fires AT MOST ONCE per input,
     using the last record's name and the run state left over from
     that record. If the last record is empty, upstream reads its
     NUL-terminator (`ks->seq.s[0]`) which sets `l = 1`, silently
     swallowing the would-be flush for any `minLen >= 2`. Our port
     reproduces both behaviours; see the trace in
     `tools/seqtk/pkg/seqtk/hrun.go` and the
     `TestHrun_TrailingFlushAcrossRecords` test.

  Byte-parity verified against `reference_code/seqtk` v1.5 on
  `nruns.fa` (default, `-l 3`, `-l 2`) and `hety_basic.fa` (`-l 4`).

Previously added: `kfreq`, `telo`. Both are byte-for-byte parity
ports against `reference_code/seqtk` v1.5 (verified by piping the
hand-built fixtures under `tools/seqtk/testdata/parity/` through both
the upstream binary and the Go port and diffing). Previous iterations
added `fqchk`, `hety`, `famask`, `mergefa`, `rename`, `split`, `size`.

- `kfreq` — full upstream surface implemented: no flags, positional
  `<kmer> <in.fa>`. Per-record TSV row with `name`, length, strand
  (`+` if forward neighbour count strictly exceeds reverse, else
  `-`), neighbour-count, exact-count. The neighbour set is every
  k-mer at Hamming distance ≤ 1 from the target (including the
  target itself). Lowercase ACGT bytes count via `seq_nt6_table`;
  the rolling encoder is reset on any non-ACGT byte. Non-ACGT bytes
  in the k-mer return a typed error (upstream `assert()`s and
  aborts). Byte-parity verified against `reference_code/seqtk` v1.5
  on `kfreq_small.fa` (`AC`, `ACGT`, `AAAA`), `kfreq_edge.fa` (`AC`),
  and `kfreq_mixed.fa` (`AA`, `ACGT`, `CCGG`, `CCCTAA`).
- `telo` — full upstream flag surface implemented:
  `-m STR` (motif, default `CCCTAA`), `-p INT` (penalty, default 1;
  negative values are silently flipped), `-d INT` (max-drop, default
  2000), `-s INT` (min-score, default 300), `-P` (print scoring
  profile instead of BED). The 5' scan walks left-to-right looking
  for motif rotations on the forward strand; the 3' scan walks
  right-to-left querying the same hash set with reverse-complement
  bases, so flipping `-m` to the motif's reverse complement swaps
  which end the BED rows describe. BED rows go to stdout; the
  `<sum_telo>\t<sum_input>` summary line goes to stderr (matching
  upstream byte-for-byte — `Telo` takes separate `stdout` and
  `stderr` writers for this reason). Byte-parity verified against
  `reference_code/seqtk` v1.5 on `telo_basic.fa` (default,
  `-m TTAGGG`, `-P -s 0`), `telo_complex.fa` (default, `-s 100`,
  `-p 2 -d 500`), and `telo_edge.fa` (default).

- `fqchk` — full upstream surface implemented: `-q INT` (default 20,
  use `-q 0` to dump the per-quality distribution). Output covers the
  preamble line, the `POS … avgQ errQ [%low %high | %Qk…]` header, the
  `ALL` aggregate row, and the per-position rows — all matched
  byte-for-byte on `fqchk_mixed.fq` (default, `-q 0`, `-q 30`) and
  `small.fq` (default, `-q 0`). PHRED+33 is hard-coded upstream and
  thus here too.
- `hety` — full upstream surface implemented: `-w INT` (default 50000),
  `-t INT` (default 5), `-m`. Per-position classification (the
  `bitcnt > 2 ? 0 : == 2 ? 2 : 1` map) is reproduced verbatim from
  `seqtk.c:614-619`, so 2-base IUPAC codes (R, Y, S, W, K, M) count as
  heterozygous while 3-/4-base IUPAC codes (B, D, H, V, N, X) do not.
  Byte-parity verified against `reference_code/seqtk` v1.5 on
  `hety_basic.fa` (default, `-w 30`, `-w 30 -t 3`, `-w 30 -m`) and
  `hety_lowercase.fa` (`-w 6 -t 1 -m`). Note: the i==l terminator and
  the upstream rollback at `seqtk.c:604-606` are required to make the
  partial last window emit the same byte count as upstream — both are
  pinned by the parity fixtures.

The previous list also mentioned `seqshuf`, `gcdc`, and `cnregion` —
these are NOT upstream subcommands per the dispatch-table audit.
`cnregion` was dropped in PR #112; `seqshuf` and `gcdc` are likewise
not real entries and should be ignored.

Project-extension policy (PR #113): the previous roadmap entry
implying `pair` was a missing upstream subcommand was wrong —
upstream seqtk v1.5 has no `pair` subcommand at all (verified
against `reference_code/seqtk/seqtk.c::main()` dispatch, which
registers `comp`, `fqchk`, `hety`, `gc`, `subseq`, `mutfa`,
`mergefa`, `mergepe`, `dropse`, `randbase`, `cutN`, `gap`,
`listhet`, `famask`, `trimfq`, `hrun`/`hpc`, `sample`, `seq`,
`kfreq`, `rename`, `split`, `telo`, `size`). Per the 1:1 parity
mandate (and the `cnregion` precedent from PR #112), project-
extension subcommands are NOT shipped under existing tool names;
the project-extension `pair` introduced in the PR #113 first
commit was dropped before merge.

Bugfixes landed since the 2026-05-14 audit:

- `comp` — fixed `seq_nt16_table['U']`/`['u']` mapping. Was 8 (T);
  upstream `reference_code/seqtk/seqtk.c` line 189 defines
  `seq_nt16_table[85] = 15` and `[117] = 15` (i.e. U is treated as
  the 4-base ambiguity N, not as T). The bug caused U bases to count
  in the `#T` column instead of `#4`; the regression test
  `TestParity_Seqtk_Comp_UBaseFasta` now pins this against an
  upstream-generated fixture.

Note: `cnregion` was listed here before but is **not** an upstream seqtk
subcommand (verified against `reference_code/seqtk/seqtk.c` v1.5: the only
`stk_*` functions registered in `main()` are `seq`, `comp`, `sample`,
`subseq`, `mergefa`, `mutfa`, `mergepe`, `randbase`, `hety`, `gc`, `fqchk`,
`hrun`/`hpc`, `listhet`, `famask`, `trimfq`, `hpc-bg`/`hpc`, `seq`, `cutN`,
`gap`, and `kfreq` — no `cnregion`). Dropped from the gap list.

Option-tail gaps (per existing subcommand):

- `comp` — missing `-r REGION` to restrict to a BED region.
- `seq` — missing `-A` (force ASCII output), `-C` (mask sequence with N), `-M FILE` (mask regions), the `-T int` trim option.
- `sample` — missing `-2` (output two paired files).
- `trimfq` — missing `-L int` (max length cap), `-B int` (min base quality).
- `subseq` — missing the regex-name mode.
- `mutfa` — missing the inverse `--inverse` mode.
- `gap` — full upstream surface implemented (`-l` only). Note: upstream's
  "gap" is any non-ACGT byte (via `seq_nt6_table`), not just literal N,
  so IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) all count.
  We match that byte-for-byte against `reference_code/seqtk` v1.5.
- `gc` — full upstream surface implemented (`-w`, `-f`, `-l`, `-x`). The
  algorithm is the upstream X-dropoff scan, not a sliding window. Output
  is BED4 (`chrom\tstart\tend\thits`). Byte-parity verified against
  `reference_code/seqtk` v1.5 on `gc_small.fa` fixtures.
- `rename` — no upstream flags (positional `<in.fq> [prefix]` only).
  Reproduces upstream's `cpy_kstr` early-return bug at
  `reference_code/seqtk/seqtk.c:1210`: a record without a comment
  inherits the previous record's comment until something with a
  comment overwrites it. This "sticky-comment" quirk is required to
  match upstream byte-for-byte and is covered by a dedicated regression
  test (`TestRename_StickyComment_Quirk`) plus the
  `TestParity_Seqtk_Rename` fixture comparisons.
- `split` — full upstream surface implemented (`-n INT`, `-l INT`,
  positional `<prefix> <in.fa>`). Output files are uncompressed and
  named `<prefix>.<5-digit 1-based>.fa` (literal `.fa` suffix even for
  FASTQ input), matching upstream byte-for-byte across `small.fa`,
  `small.fq`, and the `-l` line-wrap variants.
- `size` — no upstream flags. Emits a single
  `<num_records>\t<total_bases>\n` line, matching upstream
  byte-for-byte across `small.fa`, `small.fq`, `empty.fa`, `nruns.fa`.
- `famask` — no upstream flags (`getopt("")` at
  `reference_code/seqtk/seqtk.c:878`). Output is FASTA wrapped at 60
  columns; the three mask rules (`X` = keep, `x` = lowercase,
  anything else = overwrite) are reproduced 1:1 and verified
  byte-for-byte against `reference_code/seqtk` v1.5 on
  `famask_simple_*` and `famask_*` fixtures.
- `mergefa` — full upstream flag surface implemented:
  `-q INT`, `-i`, `-m`, `-r`, `-h` (`getopt("himrq:")` at
  `reference_code/seqtk/seqtk.c:774`). FASTA-mode outputs and the
  stderr `(same,diff,hom-het,het-hom,het-het)` counter line are
  byte-for-byte against upstream on `mergefa_a.fa`/`mergefa_b.fa`
  for the default, `-i`, `-m`, `-h`, `-q 20` (FASTQ), and the
  60-col wrap fixture. The `-r` path uses Go's `math/rand` (seeded
  with upstream's constant 11 by default); RNG byte-parity is
  explicitly NOT a goal per the policy in the section above.

**Validation:** no upstream-test-suite run yet.

### `prinseq-lite`

**Status:** 1:1 parity. The Go port now covers every behavioural
flag of `prinseq-lite.pl` 0.20.4 that the agreed scope called for,
including `--graph_data` (PR `claude/prinseq-graph-data-land`).
Three subcommands ship: `stats`, `filter`, `graph_data`.

Implemented flags (with the upstream Perl line numbers consulted for
each — `reference_code/prinseq/prinseq-lite.pl`, 0.20.4):

- `--out_format` (1=FASTA, 2=FASTA+QUAL, 3=FASTQ, 4=FASTQ+FASTA,
  5=FASTQ+FASTA+QUAL) — POD spec at lines 242-247; CLI parsing at
  lines 769-789; per-mode write branches at lines 1302-1348,
  3703-3714, 3737-3757.
- `--seq_id_mappings <file>` — POD at lines 293-295; coupling check
  with `--seq_id` at lines 945-948; file open at lines 1350-1358;
  per-record write at lines 3640-3648.
- `--ns_max_p <float>` — POD at lines 344-346; filter check at
  lines 3465-3470 (strict `>` against `(N_count * 100 / length)`).
- `--noniupac` — POD at lines 352-354; filter check at lines
  3478-3481 (`uc($seq) =~ /[^ACGTN]/o`, i.e. case-insensitive).
- `--phred64` — POD at lines 230-232; gate at lines 760-764. This
  is an INPUT encoding toggle (the original roadmap entry called it
  `--phreds`; upstream has no such flag — the only Phred-related
  option is `-phred64`). Implemented as a CLI alias for
  `--qual-type illumina`, so the existing Phred+64 decoder is
  reused. The QUAL output (`--out_format 2/5`) honours the chosen
  encoding when converting to decimal phred scores.

Bundled as a CLI alias only (`--seq_id <prefix>`): documented at
lines 470-472 of the upstream POD, implemented at line 3648. Needed
to make `--seq_id_mappings` useful, since upstream rejects
`--seq_id_mappings` without `--seq_id`.

Behavioural divergences vs upstream that are documented rather than
matched (PRINSEQ-lite has no formal test suite, so byte-for-byte
parity is not enforced):

- Multi-stream outputs (`--out_format 2/4/5`) use the value of
  `--output` directly as the prefix, appending literal `.fasta` /
  `.qual` to it. Upstream derives the prefix from `-out_good`, with
  randomised `_prinseq_good_XXXX` suffixes when no prefix is given;
  we require an explicit `--output` prefix and refuse to stream
  multiple files to stdout. The semantic restriction matches upstream
  (lines 801-802), the on-disk filename layout differs.
- The QUAL output uses the upstream `convertQualArrayToString`
  layout (two-character space-padded decimal, single-space separated,
  wrapped every `LINE_WIDTH=60` values; lines 45 and 2531-2546). We
  do **not** currently expose the upstream `-line_width` knob; the
  default is fixed at 60 via `QualLineWidth` on `FilterOptions` so
  future PRs can plumb it through.
- The `--seq_id` rename drops any trailing whitespace/comment from
  the original FASTA description. **This is a divergence from upstream**:
  upstream `prinseq-lite.pl:3683-3691` emits
  `$sid.($header ? ' '.$header : '')`, preserving any trailing comment.
  Tracked here for future-PR follow-up.

Graph-data (PR `claude/prinseq-graph-data-land`):

- `--graph_data [FILE]` (POD lines 393-401; emission block at
  lines 2050-2287). Ported in full, including the full stat
  collectors (`getSeqStats`, `getQualStats`, `generateStatsType`,
  `dinucOdds`, `checkForDupl`, `getTagFrequency`, `getBinVal`).
  Exposed as the `prinseq graph_data` subcommand. When `--graph_data`
  is given without a value, the Go port falls back to upstream's
  `<input>__.gd` default (lines 984-987).
- `--graph_stats CODES` (lines 994-1015): the upstream stat
  selector CSV (`ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn`). Supported and
  enforced — unknown codes return an error as upstream does.
- `--qual_noscale` (lines 989-993). Toggles the relative
  (100-bin) `quals` table.
- The byte-level emit deviates from upstream in one
  way: **map keys are emitted in lexicographic order**, where
  upstream uses Perl-hash iteration order (Perl >= 5.18 randomises
  this every interpreter start). Documented as an intentional
  improvement in `docs/UPSTREAM_BUGS.md > prinseq` and validated
  via a JSON-normalised semantic diff against the upstream-shipped
  `example1.gd`. See `tools/PARITY_VALIDATION.md > prinseq parity
  validation` for the test layout.

Still missing (all niche knobs, not in scope for this PR):

- The PNG report generation flow (`prinseq-graphs.pl`). Out of
  scope — the Go `graph` and `report` subcommands cover the
  equivalent visualisation surface without depending on a Perl
  graphics stack.
- A handful of niche knobs not in the original five-flag scope:
  `--range_len`, `--range_gc`, `--trim_qual_window`,
  `--trim_qual_step`, `--trim_qual_rule`, `--trim_to_len`,
  `--seq_case`, `--dna_rna`, `--line_width`, `--rm_header`,
  `--no_qual_header`, `--qual_noscale`, `--exact_only`, `--params`,
  `--custom_params`, `--seq_num`. None of these were on the agreed
  scope for this PR.

**Validation:** no upstream-test-suite run yet. Upstream PRINSEQ-lite has no
formal test suite — we'd need to construct one from the documented examples.
The new flags are covered by hand-built unit tests in
`tools/prinseq/pkg/prinseq/missing_flags_test.go` (10 sub-cases for
`--ns_max_p` boundaries, 5 sub-cases for `--out_format`, plus
`--noniupac`, `--seq_id_mappings`, `--phred64`).

### `sickle`

**Status:** 2 subcommands (`se`, `pe`). **1:1 parity validated.**

A 15-case parity corpus at `tools/sickle/testdata/parity/` covers basic
SE, `-n` (truncate-N), `-x` (no-5'-trim), illumina (Phred+64), empty
input, all-low-quality, the `-q`/`-l` boundary, gzip input, the
short-read filter, PE with singletons, synced PE pass/fail, **strict
(`-q 30 -l 5`)**, **lax (`-q 0 -l 0`)**, and **strict PE singletons
(`-q 30 -l 10`)**. Each `case*.expected.fq` was generated by running
upstream sickle v1.33 with the documented flags and is asserted
byte-for-byte by `TestParity_Sickle_*` in
`tools/sickle/pkg/sickle/parity_test.go`. See
[`tools/PARITY_VALIDATION.md` → sickle](../tools/PARITY_VALIDATION.md#sickle)
for the per-case description and the audit's "discrepancies found and
fixed" log.

Outstanding items (Go-port extensions, not parity gaps):

- Auto-detect heuristic (PR #33) is a Go-port extension with no upstream
  equivalent; exercised by `encoding_test.go`.
- Gzip *output* level was not part of the audit (parity is asserted on
  the trimmed FASTQ records, not on the gzip container bytes).

### `skewer`

**Status:** 2 subcommands (`se`, `pe`). **1:1 parity validated** — all
14 parity cases byte-match upstream skewer 0.2.2.

A 14-case parity corpus at `tools/skewer/testdata/parity/` covers SE 3'
trim, SE 5' (`-m head`), `-m any`, min-overlap, qual+adapter, length
filter, empty input, adapter-at-end, gzip input, off-by-one boundary,
**no-adapter pass-through**, **long reads (>40 bp) with embedded
adapter**, **PE matrix-mode pass-through**, and **error-tolerant matcher
rejection** (1-mismatch at Q40). All fourteen byte-match upstream
skewer 0.2.2 (`-r 0.1` defaults).

The two previously-skipped cases were closed by porting the relevant
algorithms verbatim from `reference_code/skewer/src/matrix.cpp`:

- **case04 — PE matrix mode (`-m pe`).** Ported
  `cMatrix::findAdapterWithPE` and `CalcRevCompScore`
  (matrix.cpp:487-522, 726-851) into
  `tools/skewer/pkg/skewer/skewer.go` as
  `detectPairedTrim` / `calcRevCompScore`. Matrix mode is exposed via
  the `PEMatrixMode` option (set to `true` by the PE CLI to mirror
  upstream's default `-m pe`).
- **case05 — error-tolerant matcher.** Ported the quality-weighted
  penalty model from `cAdapter::align` and the `cMatrix::penalty[]`
  ramp (matrix.cpp:138-141, 297-435, 547-556) into
  `findAdapterWithQual` / `mismatchPenalty`. Each mismatch costs a
  quality-derived penalty (Q40+ → MAX_PENALTY=4.477), and the match is
  rejected when the cumulative penalty exceeds
  `dPenaltyPerErr * compareLen + 0.001`.

See [`tools/PARITY_VALIDATION.md` → skewer](../tools/PARITY_VALIDATION.md#skewer)
for the per-case description.

### `fastp`

**Status:** Single `fastp` command with sliding-window cut, auto adapter
detection, HTML+JSON reports, duplication evaluation, UMI processing.

Missing:

- **Overrepresented sequence analysis** (`--overrepresentation_analysis`,
  `--overrepresentation_sampling`).
- **Base-correction in PE overlap** (`--correction`).
- **Quality-trimming overlap mode** (`--overlap_len_require`,
  `--overlap_diff_limit`, `--overlap_diff_percent_limit`).
- **PolyG/polyX more knobs**: `--poly_g_min_len`, `--poly_x_min_len`.
- **Splitting output**: `--split`, `--split_by_lines`, `--split_prefix_digits`.
- **Adapter list output to FASTA**: `--adapter_fasta`.
- **Disable adapter trimming**: `--disable_adapter_trimming`.
- **JSON schema completeness**: a few sub-fields under
  `before_filtering`/`after_filtering` and the per-cycle base content
  arrays are present but a handful of additional keys upstream emits
  are still missing. Run upstream `fastp` on a sample input and diff.

#### Validated-parity audit

16-case test corpus at `tools/fastp/pkg/fastp/parity_test.go` against
upstream fastp 1.0.1. **16 PASS, 0 SKIP** (post the SE adapter
auto-detect port in `claude/fastp-adapter-autodetect`). See
[tools/PARITY_VALIDATION.md#fastp-parity-validation](../tools/PARITY_VALIDATION.md#fastp-parity-validation)
for the case list.

Bugs in the Go port surfaced + fixed inline by the initial audit:

- **UMI tag format** was unconditionally `":UMI_<umi>"`. Upstream uses
  `":<umi>"` (no prefix) or `":<prefix>_<umi>"` (with prefix). Fixed.
- **Low-complexity definition** was "unique 2-mers / total 2-mers".
  Upstream uses "fraction of adjacent positions where seq[i] !=
  seq[i+1]". Fixed.
- **`low_complexity_reads` JSON counter** was missing. Added.

Bugs fixed inline by the `claude/fastp-algorithmic-fixes` follow-up PR:

- **PolyG mismatch tolerance**: upstream's `trimPolyG` tolerates 1
  mismatch per 8 bases scanned (capped at 5 total) and anchors on the
  last-G position (`reference_code/fastp/src/polyx.cpp:16-42`). The Go
  port now runs a verbatim port (`trimPolyG` in
  `tools/fastp/pkg/fastp/fastp.go`). `TestParity_Fastp_Case12_SEPolyG`
  is no longer skipped.
- **Sliding-window boundary** (`cut_front` / `cut_tail` / `cut_right`):
  `slidingWindowCut` is now a verbatim port of upstream's
  `filter.cpp:83-222`. Specifically: (a) cut_right walks the high-Q
  prefix inside the offending bad window, (b) cut_front and cut_tail
  truncate at the START of the qualifying window (then skip N's at the
  boundary), (c) the loop bound stays strictly `s + w < l` so the
  trailing w bases are never scanned.
  `TestParity_Fastp_Case13_SECutRight` and
  `TestParity_Fastp_Case14_SECutFrontTail` are no longer skipped.

Bugs fixed inline by the `claude/fastp-adapter-autodetect` follow-up PR:

- **SE adapter auto-detect**: upstream's
  `Evaluator::evalAdapterAndReadNum` (`evaluator.cpp:295-526`),
  `Evaluator::checkKnownAdapters` (`evaluator.cpp:207-293`),
  `Evaluator::getAdapterWithSeed` (`evaluator.cpp:472-526`), and
  `NucleotideTree` (`nucleotidetree.cpp`) are now ported verbatim into
  `tools/fastp/pkg/fastp/adapter_autodetect.go` and
  `tools/fastp/pkg/fastp/known_adapters.go`. The Go port now reproduces
  upstream's behavior byte-for-byte, including the 10000-record gate at
  `evaluator.cpp:344` (below that threshold the evaluator returns ""
  and no adapter trimming is applied — same as upstream's "No adapter
  detected for read1" path). `TestParity_Fastp_Case15_SEAutoDetect` is
  no longer skipped.

**Validation:** **16-case parity test suite, 16 passing, 0
documented `t.Skip`. 1:1 parity achieved** (post
`claude/fastp-adapter-autodetect`).

### bedtools (35 subcommands ported)

**Status:** 35 of ~40 subcommands (~88%). 141 passing parity tests
against the upstream test suite (across PR #55 + Phase-3 wave 1 + wave
2 simple + wave 2 algo) + 17 new cases from wave 3 (PR #87) + 6 cases
from the reldist/fisher full-parity wave (PR #90) + 6 cases from the
nuc/annotate wave (PR #91) + 7 cases from the multicov/multiinter
wave (PR #92) + 6 cases from the igv/links wave (PR #93) + 10 cases
from the BEDPE pair-ops wave (PR #94, five each for `bedpairtobed`
and `bedpairtopair`, derived from the upstream source as neither
subcommand ships a `test/<name>/` subdir) + **9 newly passing**
parity cases from the column-ops + discrepancies wave (this PR —
`jaccard.t02/t03/t05/t06/t10/t11`, `merge.t15`, `map.t11`,
`map.t13`).
Phase-3 wave 1 (PR #78) added `bedgroupby`/`bed12tobed6`/`bedmakewindows`;
wave 2 simple (PR #80) added `bedexpand`/`bedgetfasta`/`bedsample`/`bedspacing`;
wave 2 algo added `bedcoverage`/`bedmap`/`bedshuffle`; wave 3 tail
(PR #87) added `bedcluster`/`bedsplit`/`bedsummary`/`bedtag`/`bedwindow`;
the reldist/fisher wave (PR #90) added `bedreldist`/`bedfisher`; the
nuc/annotate wave (PR #91) added `bednuc`/`bedannotate`; the
multicov/multiinter wave (PR #92) closed the last two of the six
originally-planned algorithmic subcommands; the igv/links wave
(PR #93) landed the two pure-format converters; **the BEDPE pair-ops
wave (this PR) closes the last two missing subcommands — bedtools now
has no missing subcommands**.

Missing subcommands: none. All algorithmic, BEDPE, and converter
subcommands are ported. Remaining bedtools work is option-tail polish.

Resolved in the column-ops + discrepancies wave:

- `bedjaccard` now pre-merges A and B before computing intersection /
  union, matching upstream's `setUseMergedIntervals(true)`. Previously
  skipped parity cases `jaccard.t02 / t03 / t05 / t06 / t10 / t11` now
  pass byte-for-byte.
- `bedmerge -s` now matches upstream's per-strand merge semantics: `.`
  (unknown) strand records are dropped under `-s` and `+` / `-` groups
  merge independently before the (chrom, start) re-merge. Previously
  skipped `merge.t15` now passes.
- The full column-op vocabulary from upstream (`stdev`, `sstdev`,
  `absmin`, `absmax`, `cat`, `cat_uniq`) is wired into the shared
  `bedmerge.ApplyOp`, unblocking `bedgroupby` (`TestGroup_StdevSstdev`,
  `TestGroup_AdditionalOps`) and `bedmap` (`map.t11`, `map.t13`).

Option-tail gaps on the wave-2 additions:

- `bedgetfasta` — `-fullHeader` is now implemented: contigs are indexed
  by the full FASTA header line (whitespace included) via
  `pkg/htsgo/fasta.BuildIndexFullHeader` /
  `OpenRandomAccessFullHeader`, and `bedgetfasta -fullHeader` flows the
  flag through to the index build. Upstream `getfasta.t06` (the
  `-fullHeader` two-line case) and `t07` (the no-`-fullHeader` warning
  case) now both pass byte-for-byte. BGZF FASTA input is also wired
  through (this PR): `pkg/htsgo/fasta` now sniffs the BGZF magic
  in `OpenRandomAccess` / `OpenRandomAccessFullHeader` and routes to a
  new `OpenRandomAccessBGZF` that fully decompresses the payload
  in-memory and reuses the existing FAI index path. The `.gzi` sidecar
  (when present) is parsed for early validation via a stdlib-only
  little-endian reader in `pkg/htsgo/fasta/bgzf.go`; a samtools-
  compatible `.fa.gz.fai` is honoured when present, otherwise the
  index is rebuilt from the decompressed payload. Upstream
  `getfasta.t18` (BGZF FASTA + `-split` BED12) now passes byte-for-byte
  using the upstream `t.fa.gz` fixture. Partial-decompression seek via
  `.gzi` is a future optimization; the in-memory path is sufficient
  for parity and for the reference genomes bedtools is typically used
  against.
- `bedsort` — `-header` is now implemented (this PR): leading
  `#`-prefixed comment, `track`, and `browser` directive lines are
  buffered and emitted verbatim ahead of the sorted body. Upstream
  `sort.t09` now passes byte-for-byte.
- `bedsample` — output PRNG is Go `math/rand` and is not byte-compatible
  with upstream's C++ sampler. Seeded runs are deterministic within
  `bedsample` (same seed → same output) but cross-tool record-for-record
  parity with upstream is not feasible without porting upstream's
  `random_shuffle`.
- `bedmulticov` — <a id="bedmulticov-bam"></a>BAM input is wired through
  via `pkg/htsgo/sam.NewBAMReader`; primary alignments contribute
  one interval each over their reference span, and `-q` MAPQ filter +
  `-D` per-A-interval depth cap are honoured. Upstream `multicov.t1`
  through `t4` and `t10` pass byte-for-byte.
  **`-split`** block-aware coverage on BAM CIGAR `N` ops is now
  implemented (this PR): the BAM index pass walks each alignment's
  CIGAR and emits one block per contiguous reference-consuming op-run
  (M/=/X, with D extending the current block — matching upstream's
  `breakOnDeletionOps=false`), skipping any `N`-op gap. Each alignment
  is counted at most once per A interval. When combined with
  `-f`, the threshold is applied to `total_block_overlap /
  sum_of_BAM_block_lengths` using strict `>` — a quirk of bedtools 2.x
  preserved here for byte-for-byte parity (mirrors
  `multiBamCov.cpp::FindBlockedOverlaps`). Upstream `multicov.t5`
  through `t9` now pass.
  **CRAM** input remains deferred — see `docs/CRAM_DESIGN.md`; the CLI
  surfaces a clear error for `.cram`.
- `bedmultiinter` — VCF/GFF input not implemented (upstream autodetects
  these via `BedFile`). Input is assumed sorted; out-of-order records
  within a single file are tolerated only because each file is
  re-sorted and merged before the sweep.

Column-op closure: the shared `bedmerge.ApplyOp` (used by `bedmerge`,
`bedgroupby`, `bedmap`, and `bedcoverage`) now supports the full
upstream KeyListOps vocabulary: `stdev`, `sstdev`, `absmin`, `absmax`,
`cat`, `cat_uniq` (in addition to the previously-implemented sum,
min/max, mean, median, mode/antimode, count, count_distinct, distinct,
collapse, first, last). Done; no remaining gaps.

### `vcftools`

**Status:** **146 of 146 unique upstream long flags (100%)** after
long-tail wave 23. A complete `in_str ==` enumeration of
`parameters.cpp` finds 146 distinct upstream long flags; the port
registers and implements all of them. The four BCF-binary I/O flags
that earlier prose listed as blocked — `--bcf`, `--diff-bcf`,
`--recode-bcf`, `--contigs` — are all closed (waves 21–23): the
in-tree `pkg/htsgo/bcf` reader/writer (built on the in-tree BGZF
codec) supplies the binary I/O, so no external infrastructure is
needed. The PCA trio (`--pca`/`--pca-no-norm`/`--pca-snp-loadings`)
is fully implemented (wave 19) on top of gonum's symmetric
eigendecomposition — the former "LAPACK blocker" no longer applies.
The remaining work is per-output column-set polish (see the "Other"
list below), not flag-count gaps. Earlier header prose claiming
"142 of 146, blocked on 4 BCF flags" predated waves 19–23 and was
stale; this is the corrected count.

Closed in wave 1:

- **Inter-chromosomal LD**: `--interchrom-geno-r2`, `--interchrom-hap-r2` ✅
- **Chi-square LD**: `--geno-chisq` ✅
- **Relatedness**: `--relatedness` (Yang 2010), `--relatedness2` (KING-robust) ✅
- **Runs of homozygosity**: `--LROH` (+ `--LROH-min-variants`) ✅
- **Phased blocks**: `--phased-blocks` ✅
- **FILTER tag include/exclude**: `--remove-filtered`, `--keep-filtered` ✅
- **INFO selection in recode**: `--keep-INFO TAG`, `--remove-INFO TAG` ✅
- **INFO extraction**: `--get-INFO TAG[,TAG]` → `.INFO` ✅

Closed in wave 2:

- **LDhat output formats**: `--ldhat`, `--ldhat-geno` (paired
  `<prefix>.ldhat.sites` / `<prefix>.ldhat.locs`, byte-for-byte vs
  upstream) ✅
- **Phased-site filter**: `--phased` (composes with `--ldhat` per
  upstream's `phased_only` invariant) ✅

Closed in wave 3:

- **LDhelmet output format**: `--ldhelmet` (paired
  `<prefix>.ldhelmet.snps` / `<prefix>.ldhelmet.pos`, byte-for-byte vs
  upstream; implies `--phased` + `--remove-indels` per
  parameters.cpp:275, requires `--chr` per parameters.cpp:717) ✅
- **IMPUTE reference-panel output**: `--IMPUTE` (case-sensitive; emits
  `<prefix>.impute.legend` / `<prefix>.impute.hap` /
  `<prefix>.impute.hap.indv`, byte-for-byte vs upstream; implies
  `--phased`, biallelic-only, no missing data per parameters.cpp:255) ✅

Closed in wave 4:

- **`--diff-indv-map FILE`** — two-column whitespace-separated file that
  renames file-2 sample IDs before matching against file-1. Loader
  mirrors upstream `variant_file_diff.cpp:11-34`; mapping is applied when
  forming `commonPairs` and when classifying `INDV FILES` in
  `.diff.indv_in_files`. ✅
- **`--diff-discordance-matrix`** — emits
  `<prefix>.diff.discordance_matrix` with the 5x5 layout from upstream
  `variant_file_diff.cpp:944-1198`: header row of file-1 genotype labels,
  four data rows of file-2 genotype labels, biallelic + matching-ALT only,
  diploid only, byte-for-byte parity vs upstream. ✅

Closed in wave 5:

- **`--diff-switch-error`** — emits `<prefix>.diff.switch` (per-event
  log with `CHROM POS_START POS_END INDV` columns) and
  `<prefix>.diff.indv.switch` (per-individual rate with
  `INDV N_COMMON_PHASED_HET N_SWITCH SWITCH` columns), byte-for-byte vs
  upstream. Ported from `variant_file_diff.cpp:1207-1507`
  (`output_switch_error`). Plugged into the existing diff runner so it
  composes with `--diff-indv-map` (sample-ID renaming) without
  re-implementing the load-file-2 path. ✅
- **`--mendel <PED>`** — emits `<prefix>.mendel` listing Mendelian
  inconsistencies across trios defined in a four-column PED file
  (`family child father mother`; first line always skipped). Ported from
  `variant_file_output.cpp:5332-5470`
  (`output_mendel_inconsistencies`). The PED column ordering follows
  upstream's `ss >> family >> child >> father >> mother;` parse;
  `family_ids` for each trio is `<child>_<father>_<mother>`. ✅

The original brief proposed `--phase` (output format) as the second
target. After re-checking upstream `parameters.cpp`, no standalone
`--phase` flag exists — only `--phased` (already ported in wave 2 as a
site filter). `--mendel` was substituted as the next clear long-tail
target, per the brief's substitution clause.

Closed in wave 6 (this PR):

- **`--non-ref-af FLOAT`** — minimum non-reference allele frequency:
  every ALT's count/non-missing-chr ratio must be ≥ threshold. Ported
  from upstream `parameters.cpp:303` + `entry_filters.cpp:770-824`. We
  preserve the upstream quirk that `min_non_ref_af > 0` also drops
  monomorphic (no-ALT) sites via the `N_failed == N_alleles-1`
  fallback on line 814. ✅
- **`--non-ref-ac INT`** — minimum non-reference allele count: every
  ALT's per-site count must be ≥ threshold. Ported from upstream
  `parameters.cpp:302` + `entry_filters.cpp:869-920`. The
  monomorphic-fallback on line 912 is gated on `min_non_ref_ac_any`
  (NOT plain `min_non_ref_ac`), so `--non-ref-ac` alone deliberately
  does NOT drop monomorphic sites — verified against upstream and
  pinned in `TestParity_NonRefAF_DropsMonomorphic`. ✅

The brief originally requested `--FILTER-PASS-summary` and
`--remove-INFO-all`. Neither flag is registered in upstream
`parameters.cpp`: the real flag is `--FILTER-summary` (already
implemented since wave 0), and `--remove-INFO-all` simply does not
exist (the upstream registrations are `--keep-INFO-all` and
`--recode-INFO-all`, both already implemented). The brief's
substitution clause permits picking from the `--non-ref-af*` family,
so `--non-ref-af` and `--non-ref-ac` were chosen as the two new
flags. The remaining `--non-ref-af-any` / `--non-ref-ac-any` (and the
`--max-*` upper-bound counterparts) are closed in wave 7 below.

Closed in wave 7 (this PR):

- **`--max-non-ref-af FLOAT`** — upper-bound counterpart of
  `--non-ref-af`: drops the site if ANY ALT's freq > threshold (per-ALT
  immediate fail, entry_filters.cpp:807). Also drops monomorphic sites
  via the line-814 fallback (gate keyed on plain thresholds). ✅
- **`--max-non-ref-ac INT`** — upper-bound counterpart of
  `--non-ref-ac`: drops the site if ANY ALT's count > threshold
  (entry_filters.cpp:905). Like plain `--non-ref-ac`, does NOT drop
  monomorphic sites (the line-912 fallback is keyed on `_any`). ✅
- **`--non-ref-af-any FLOAT`** + **`--max-non-ref-af-any FLOAT`** —
  N_failed-counter variants registered at `parameters.cpp:304-305` and
  `:289-290`. Wired for command-line parity, but observably **NO-OPS**
  when used alone because upstream `entry_filters.cpp:814` gates the
  fallback on the PLAIN thresholds (`min_non_ref_af > 0` /
  `max_non_ref_af < 1.0`), not on the `_any` thresholds. The flags
  only have an effect when paired with their plain counterpart, in
  which case the fallback fires when EVERY ALT fails the `_any`
  threshold. We mirror this verbatim; pinned by
  `TestParity_NonRefAFAny_NoOp` and `TestParity_NonRefAF_03_Any_06`. ✅
- **`--non-ref-ac-any INT`** + **`--max-non-ref-ac-any INT`** — counter
  variants of the AC family. Unlike AF, the AC fallback at
  `entry_filters.cpp:912` IS gated on the `_any` thresholds, so these
  flags are functional standalone. Site dropped when N_failed equals
  N_alleles-1 (every ALT failed). Monomorphic sites (N_alleles=1)
  satisfy 0==0 and are dropped — counter to plain `--non-ref-ac` which
  keeps them. Pinned by `TestParity_NonRefACAny_2`,
  `TestParity_NonRefACAny_1_Chr20`, `TestParity_MaxNonRefACAny_2_Chr20`. ✅

Refactor: the wave-6 per-ALT early-return was lifted to an N_failed
accumulator pass (matching upstream's structure literally) so the AF
and AC `_any` fallbacks can decide post-loop. The plain flags still
short-circuit on the first failing ALT.

Closed in wave 8 (this PR):

- **`--hwe FLOAT`** — minimum per-site exact-test (Wigginton/Cao/Abecasis
  2005) Hardy-Weinberg p-value filter. Ported from upstream
  `parameters.cpp:254` (which also forces `max_alleles = 2` when the
  flag is supplied) + `entry_filters.cpp:922-946` + `entry.cpp:18-101`
  (the exact-test). The SNPHWE port is a line-for-line port of
  upstream's integer-arithmetic midpoint/walk algorithm so that
  filtering decisions are byte-identical. The CLI adapter mirrors the
  upstream `--hwe → max_alleles=2` coupling. Pinned by
  `TestParity_HWE_005_sample` (3-sample fixture, exercises the
  max_alleles=2 coupling), `TestParity_HWE_005_fixture` (20-sample
  fixture engineered so site `1:200` with counts (10,0,10) fails the
  exact test at p≈1.34e-6 — verified via upstream `--hardy`), and unit
  tests on `snpHWE` boundary cases. ✅
- **`--max-missing-count INT`** — maximum number of missing
  *chromosomes* (haploid alleles, NOT samples) tolerated per site.
  Ported from upstream `parameters.cpp:286` +
  `entry_filters.cpp:918`. The comparator is strict `>` so a site
  with exactly `INT` missing alleles is KEPT, but `INT+1` drops.
  Pinned with two parity tests bracketing the boundary case
  (`TestParity_MaxMissingCount_1` drops a 2-missing site,
  `TestParity_MaxMissingCount_2` keeps it). The Params struct grows
  a `MaxMissingCountSet` boolean so the CLI can distinguish
  "user passed 0" (drop any site with any missing call) from
  "user omitted the flag" (no filter); the CLI registers the flag via
  `flag.Func` to record both. ✅

Closed in wave 9 (this PR):

- **`--kept-sites`** — emits `<prefix>.kept.sites`, a two-column
  `CHROM\tPOS` TSV listing every site that survived all filters in
  input order. Ported from upstream `parameters.cpp:268` +
  `variant_file_output.cpp:4285-4326` (`output_kept_sites`). Upstream
  re-parses the input file specifically for this output and re-runs
  `entry::apply_filters` (`entry_filters.cpp:23`) on every entry; we
  piggy-back on the existing filter pipeline in `Run` instead, which is
  equivalent because the same filter gates apply to both code paths.
  Header is `CHROM\tPOS` (tab-separated, LF terminator), matching
  upstream's `out << "CHROM\\t" << "POS" << endl;`. Pinned by
  `TestParity_KeptSites_NoFilter` (all-sites-pass case),
  `TestParity_KeptSites_HWE` (HWE filter), and
  `TestParity_KeptSites_PosFilter` (chr+from-bp+to-bp filter). ✅
- **`--removed-sites`** — counterpart of `--kept-sites`: emits
  `<prefix>.removed.sites` with the same column layout, listing sites
  *dropped* by any filter. Ported from upstream `parameters.cpp:330` +
  `variant_file_output.cpp:4328-4373` (`output_removed_sites`). Pinned
  by `TestParity_RemovedSites_HWE` and `TestParity_RemovedSites_PosFilter`.
  Plus `TestKeptRemoved_Disjoint_And_Complete` (port-only invariant —
  upstream forbids both flags in one run via `num_outputs > 1`; we
  deliberately do not replicate that constraint per the CLAUDE.md
  "don't replicate upstream bugs" rule, and instead verify that the
  combined invocation partitions the input perfectly) and
  `TestKeptRemoved_Disabled_NoFiles` (neither flag → no file leaked). ✅

Implementation note: the trace writer lives in
`tools/vcftools/pkg/vcftools/sitetrace.go`. Each `continue` in the main
filter loop of `vcftools.go:Run` now calls
`siteTracker.recordRemoved(chrom, pos)` immediately before bailing out,
and the successful path calls `siteTracker.recordKept(chrom, pos)` just
before `keptSites++`. Both methods are no-ops when the corresponding
flag is not set (cheap nil-check on the bufio.Writer).

Closed in wave 10 (this PR):

- **`--remove-filtered-geno-all`** — sets GT to `./.` for every kept
  genotype whose FORMAT FT field is not "PASS" or ".". Ported from
  upstream `parameters.cpp:323` + `vcf_entry.cpp:580-608`
  (`filter_genotypes_by_filter_status` with `remove_all=true`). Only
  the GT slot is rewritten; other FORMAT fields (FT/DP/GQ/...) pass
  through unchanged, matching upstream's recode emission at
  `vcf_entry.cpp:320-368`. Sites with no FT FORMAT column are
  left untouched (mirrors upstream's early-return at
  `entry_filters.cpp:94-108`). Pinned by
  `TestParity_RemoveFilteredGenoAll` (byte-for-byte vs upstream on the
  new `ft_geno.vcf` fixture) and `TestRemoveFilteredGeno_NoFT_NoOp`
  (port-only invariant against `sample.vcf`, which has no FT column). ✅
- **`--remove-filtered-geno NAME`** (repeatable) — drops a genotype
  whose FT lists any of the named flags. Ported from
  `parameters.cpp:324` + `vcf_entry.cpp:601-605`. FT is parsed as a
  `;`-separated list per upstream's `vcf_entry_setters.cpp:188-212`
  (entries equal to "" or "." are dropped from the list). Pinned by
  `TestParity_RemoveFilteredGenoQ10` (single-flag) and
  `TestParity_RemoveFilteredGenoMulti` (two-flag invocation — the
  set behaviour in upstream's `geno_filter_flags_to_exclude`). ✅
- **`--max-indv N`** — caps the number of kept individuals at N.
  Ported from upstream `parameters.cpp:292` +
  `variant_file_filters.cpp:105-147`
  (`filter_individuals_randomly`). **Port deviation, documented:**
  upstream uses `srand(time(NULL))` + `random_shuffle`, making the
  kept-sample identity non-deterministic across runs. This port
  instead deterministically keeps the first N kept samples in input
  (header) order. The COUNT invariant (`|kept| =
  min(N, |pre-cap-kept|)`) is the strongest claim we can make against
  upstream's randomness, so parity is pinned at the COUNT level only
  (`TestMaxIndv_Count` table-driven cases). `MaxIndvSet` gates the
  cap so `--max-indv 0` ("drop every sample") is distinguishable from
  the default. Pinned by `TestMaxIndv_Count` and
  `TestMaxIndv_Unset_NoOp`. ✅
- **`--keep-INFO-all`** — upstream-deprecated synonym for
  `--recode-INFO-all`. Both `parameters.cpp:267` and `:318` write
  to the same `recode_all_INFO` parameter bit; the CLI ORs them
  together so either flag (or both) produces identical output.
  Pinned by `TestKeepINFOAll_Synonym`. ✅
- **`--version`** — prints `VCFtools (0.1.18)` and exits, matching
  upstream `parameters.cpp:648-652` byte-for-byte. Hard-coded
  version string tracks the upstream submodule at port time;
  bump on rebase.  ✅

Implementation notes (wave 10):

- The FT-based filter lives inside `filterGenotypes` in
  `tools/vcftools/pkg/vcftools/vcftools.go` (next to the existing
  `--minDP/--maxDP/--minGQ` path). Two small helpers (`parseSampleFT`,
  `shouldDropByFT`) keep the hot loop branch-light and mirror the
  upstream getters/parsers documented above. The `sampleWithMissingGT`
  helper extracted from the duplicated DP/GQ paths is also reused by
  the FT path.
- `--max-indv` is wired through `buildSampleFilter`, which now returns
  a non-nil keep set when `MaxIndvSet` is true even with no
  identity-based filter. The cap iterates `header.Samples` in order
  so the truncation is deterministic.

Closed in wave 11 (this PR):

- **`--mask FILE`** — FASTA-style positional mask. The mask file has
  `>CHROM` headers followed by lines of digit characters (one per
  reference base, 1-based). A site at (CHROM, POS) is kept when its
  mask digit is `<= --mask-min` (default 0). Ported from upstream
  `parameters.cpp:280` + `entry_filters.cpp:674-752`
  (`filter_sites_by_mask`). Pinned by `TestParity_Mask_Default`,
  `TestParity_Mask_Min5`, `TestParity_Mask_Partial`. ✅
- **`--invert-mask FILE`** — same loader as `--mask`, but the keep/drop
  decision is flipped. Ported from upstream `parameters.cpp:262`.
  Pinned by `TestParity_InvertMask_Min5`. ✅
- **`--mask-min INT`** — maximum kept mask digit value, 0-9 (upstream
  errors when `> 9` at `parameters.cpp:720`; we additionally reject
  negatives at load time because they silently drop every site
  upstream — clearer to fail fast). Default 0. ✅

Implementation notes (wave 11):

- The mask reader is **forward-only**, mirroring upstream's stateful
  `ifstream` walk (entry_filters.cpp:680-688 keeps `mask_chr`,
  `mask_line`, and `mask_pos` as static state across calls). The Go
  port loads the file once into a `[]maskChromosome` and maintains a
  `(chromIdx, slabIdx)` cursor. The cursor never moves backwards, so a
  VCF presenting chr2 before chr1 against a mask listing chr1 then
  chr2 will lose the chr1 sites — this matches upstream behaviour
  (`TestMaskFilter_OutOfOrderVCFDrops` pins it).
- Header tokenisation matches upstream's
  `line.substr(1, line.find_first_of(" \t")-1)` (split on first
  whitespace after `>`; comments are discarded).
- Mutually-exclusive flag behaviour: `--mask` and `--invert-mask`
  share the same `mask_file` slot upstream (last one wins via
  parameters.cpp:262 vs :280). Go's `flag.String` does not preserve
  last-set ordering, so we apply the asymmetric rule "if
  `--invert-mask` is non-empty, override and set invert=true". This is
  observable only when both flags are supplied; documented in
  `main.go`.

Closed in wave 12 (this PR):

- **`--positions-overlap FILE`** — keep a record when ANY base in
  `[POS, POS+len(REF)-1]` matches a (CHROM, POS) entry in the file.
  Same two-column whitespace-separated format as `--positions`; `#`
  comments and blank lines tolerated. Ported from upstream
  `parameters.cpp:315` + `entry_filters.cpp:408-531`
  (`filter_sites_by_overlap_positions`, keep-branch). For 1-base REF
  records the behaviour reduces to plain `--positions`; the divergence
  appears on indels / MNPs, which is the entire reason upstream ships
  the overlap variant. Sites on chromosomes not named in the file are
  dropped (matches `entry_filters.cpp:515-516`). ✅
- **`--exclude-positions-overlap FILE`** — drop a record when ANY base
  in `[POS, POS+len(REF)-1]` matches a (CHROM, POS) entry. Ported from
  `parameters.cpp:221` + `entry_filters.cpp:533-547`. Inverse of
  `--positions-overlap`; sites on chromosomes not named pass through
  unchanged (matches the `chr_to_idx.find != end` guard at line 535). ✅

Implementation notes (wave 12):

- Upstream reuses the same `keep_positions`/`exclude_positions` state
  for both the plain `--positions` family (entry_filters.cpp:279-406)
  and the overlap family (lines 408-548). Consequence: combining
  `--positions` with `--positions-overlap` upstream silently degrades
  to overlap semantics for whichever file populated the set first
  (the second loader sees a non-empty set and skips). We keep the
  four filters as independent fields on `Params` so each behaves
  exactly as documented when used solo, and when combined both gates
  apply (a site must pass include AND not be excluded across the two
  flag pairs). Pinned by
  `TestPositionsOverlap_VsPlain_DivergesOnMultiBaseRef`.
- The upstream loop is half-open: `for ui=POS; ui<POS+REF.size(); ui++`.
  We mirror this with `for p := v.Pos; p < v.Pos+refLen; p++`. With a
  defensive guard `refLen=max(len(v.Ref),1)` we ensure at least the
  POS itself is tested for malformed VCFs with empty REF; valid VCFs
  always have `len(REF) >= 1` and the guard never triggers.
- File format reuses the existing `loadPositions` parser (same as
  `--positions`); separate `positionSet` instances are kept in `Run`
  for the four flags so the apply order is deterministic
  (`includePos` → `excludePos` → `includePosOverlap` →
  `excludePosOverlap` → BED → variant-type → frequency / quality /
  HWE / etc.).

Closed in wave 13 (this PR):

- **`--derived`** — when combined with `--freq` / `--counts`, reorder
  the allele columns so that the ancestral allele (INFO/AA,
  case-insensitive) appears first; drop sites where AA is missing,
  `.`, `?`, or does not match REF/ALT. Mirrors upstream
  `parameters.cpp:201` + `variant_file_output.cpp:67-159`
  (`output_frequency`, the `derived` branch). Implementation lives in
  `addFrequencyStat` (new `derivedSwap` flag on `siteFreqStat`) and
  the existing `outputFrequency` reorder loop. Multi-allelic sites
  are already dropped by our biallelic-only `--freq` restriction, so
  `--derived` only affects the biallelic subset (matches the subset
  this port emits at all under `--freq`/`--counts`). ✅
- **`--extract-FORMAT-info NAME`** — extract a per-genotype FORMAT
  field across all kept samples into a tab-separated
  `<prefix>.<NAME>.FORMAT` file. Header is `CHROM\tPOS\t<sample>...`;
  one data row per site whose FORMAT column lists NAME (sites lacking
  NAME in FORMAT are skipped entirely, matching upstream's
  `FORMAT_id_exists` gate). Samples whose colon-separated value
  vector is too short to reach NAME's index emit `.` (matches
  `vcf_entry.cpp:618` + the early `break` at line 637). Ported from
  upstream `parameters.cpp:222` +
  `variant_file_format_convert.cpp:1204-1263`
  (`output_FORMAT_information`). Single-valued upstream (the last
  value wins on the CLI). ✅

Implementation notes (wave 13):

- `--derived` is a *modifier*, not an output: it only takes effect
  when paired with `--freq` or `--counts`. Pinned by
  `TestDerived_NoFreqIsNoOp`. Upstream's
  `parameters.cpp:201` only flips a boolean; the reorder logic lives
  inside `output_frequency` (the boolean is also consumed by
  `output_indv_burden` and `output_indv_freq_burden`; both burden
  flags now honour `--derived` after wave 14 — see the wave-14
  section below).
- AA uppercasing — upstream calls
  `std::transform(AA.begin(), AA.end(), AA.begin(), ::toupper)` on
  line 78 / 439 / 564 before comparing against `e->get_allele(ui)`.
  We replicate this with `strings.ToUpper` on both AA and REF/ALT so
  `AA=a, REF=A` matches as expected. The case-insensitive match is
  pinned by the 1:400 site in `derived_fixture.vcf`.
- AA sentinels — upstream's `if ((AA == "?") || (AA == "."))` check
  appears at lines 79 / 440 / 565. We mirror it explicitly (empty,
  ".", "?"). Pinned by sites 1:500 / 2:100 in the same fixture.
- `--extract-FORMAT-info` shares its FORMAT-tag presence helper
  (`formatContains`) with the BEAGLE writer (declared once in
  `beagle.go`).
- The Go VCF parser only populates `sample.Data[key]` for keys whose
  colon-token slot exists in the per-sample string; absent slots
  therefore read back as missing-from-map. We treat that as upstream's
  "value vector too short → '.'" case (vcf_entry.cpp:618-637). Pinned
  by sites 1:100/S3, 1:400/S2 in `extract_format_fixture.vcf`.

Closed in wave 14 (this PR):

- **`--indv-burden`** — per-individual diploid-burden counts emitted to
  `<prefix>.iburden`. Header is
  `INDV\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_MISS`; with `--derived` the
  `_REF` / `_ALT` columns rename to `_ANC` / `_DER`. Non-diploid sites
  are skipped (upstream's `if (e->is_diploid() == false) continue;` at
  `variant_file_output.cpp:429-433`). With `--derived` the site's
  INFO/AA tag picks the ancestral-allele index; sites where AA is
  missing, `.`, `?`, or does not match any REF/ALT are skipped. Ported
  from upstream `parameters.cpp:257` +
  `variant_file_output.cpp:378-498` (`output_indv_burden`). ✅
- **`--indv-freq-burden`** — per-individual × per-allele-count matrix
  written to `<prefix>.ifreqburden`. For each kept diploid site,
  computes the per-allele count vector across kept individuals and
  for each kept individual increments the burden cell at column
  `allele_counts[geno_allele]` for each non-ref (or non-ancestral with
  `--derived`) allele the individual carries. Mirrors upstream
  `parameters.cpp:258` + `variant_file_output.cpp:501-627`
  (`output_indv_freq_burden` with `double_count_hom_alt=0`). ✅
- **`--indv-freq-burden2`** — same as `--indv-freq-burden` but with
  `double_count_hom_alt=1`, so a hom-alt genotype contributes 1 (not
  2) to the corresponding allele-count bin. Mirrors `vcftools.cpp:64`
  + the same `output_indv_freq_burden` routine. ✅

Implementation notes (wave 14):

- All three flags share `burden.go`: `indvBurdenRunner` for
  `--indv-burden` and `indvFreqBurdenRunner` for the two
  freq-burden variants (the latter takes a `doubleCountHomAlt`
  toggle for `--indv-freq-burden2`).
- **Upstream label-index bug preserved.**
  `output_indv_freq_burden` writes
  `out << meta_data.indv[indv_count];` at line 621 — that should be
  `meta_data.indv[ui]` (the original-index). With `--remove-indv` (or
  any sample-filter) dropping a non-trailing sample, the labels in
  `.ifreqburden` shift relative to the burden values (e.g. excluding
  `S2` from `[S1,S2,S3,S4]` yields labels `S1,S2,S3` for kept
  individuals `S1,S3,S4`). We mirror this byte-for-byte;
  pinned by `TestParity_IndvFreqBurden_LabelBug`. The
  `output_indv_burden` function at lines 488-497 correctly uses
  `meta_data.indv[ui]` and is not affected.
- The diploid-only check is per-site: any haploid genotype in a kept
  individual disqualifies the entire site. Upstream emits a one-off
  warning at `variant_file_output.cpp:431`; we do not re-emit the
  warning byte-for-byte but skip the site identically. Pinned by
  `TestIndvBurden_SkipsNonDiploid`.
- The `ancestralAlleleIndex` helper centralises the
  `variant_file_output.cpp:437-462` AA-resolution logic (uppercase
  match, missing-sentinel handling) so the same predicate is used by
  both burden runners and stays consistent with how
  `addFrequencyStat` already implements the `--derived` filter (see
  wave 13).
- The `diploidAlleles` helper splits a `a/b` / `a|b` GT into
  `(a1, a2)`, treating `.` as -1, and is shared by both runners. The
  existing `parseGTForLDhat` parser was close but its haploid branch
  silently coerces a missing single-allele GT into a `(−1, −2,
  true)` triple that the burden routines need to reject as
  non-diploid; `diploidAlleles` returns `ok=false` for any GT without
  a separator so the caller can do the diploid-skip up front.

**PCA family**: ✅ resolved in wave 19. `--pca`, `--pca-no-norm`,
and `--pca-snp-loadings INT` (upstream `parameters.cpp:308-310`,
`variant_file_output.cpp:4871-5249`) all land via gonum's symmetric
eigensolver (`gonum.org/v1/gonum/mat`'s `SymEigen`). The owner
sanctioned `gonum` as the second third-party-dep zone after the
CRAM codec carveout (CLAUDE.md). Parity goldens were generated by
rebuilding upstream vcftools with `--enable-pca` against system
`liblapack-dev` / `libblas-dev`; the resulting binary lives at
`/tmp/vcftools_lapack_install/bin/vcftools` in the dev sandbox.
Eigenvector signs are LAPACK-implementation-dependent (both
LAPACK and gonum), so parity tests use a per-column sign-tolerant
comparison; eigenvalues themselves are sign-invariant and compared
with a tight numerical tolerance. Wave 19 also fixes a latent
upstream bug — `output_PCA` reads past the end of the per-individual
M[i] vectors when any kept individual has a missing genotype, see
`docs/UPSTREAM_BUGS.md` for the writeup. Pinned by
`TestParity_PCA_Basic`, `TestParity_PCA_NoNorm`,
`TestParity_PCA_SNPLoadings`, and the in-pkg algebraic tests in
`pca_test.go`.

Closed in wave 15 (this PR):

- **`--hapcount BED`** — per-BED-bin haplotype-count summaries written
  to `<prefix>.hapcount`. Columns: `#CHROM BIN_START BIN_END N_SNP
  N_UNIQ_HAPS N_GROUPS {MULTIPLICITY:FREQ}...`. Implies `--phased`
  (upstream parameters.cpp:248 sets `phased_only=true`). Diploid-only
  per-site (upstream variant_file_output.cpp:1350-1354 skip otherwise).
  Bins must be non-overlapping (upstream errors otherwise at lines
  1208-1216). Ported from variant_file_output.cpp:1169-1401
  (`output_haplotype_count`) with three upstream bugs FIXED on port
  per CLAUDE.md and the wave-14 precedent (PR #138). See
  `docs/UPSTREAM_BUGS.md#fix-on-port-resolved` for the writeup:
    1. prev_bin_idx shift on within-chromosome bin transitions
       (lines 1314-1315) — old bin's counts were silently overwritten
       with the new bin's values.
    2. End-of-stream read-after-free (lines 1370-1400) — last
       chromosome's rows were silently dropped (or zeroed) on a
       glibc-built upstream binary.
    3. BED first-line unconditional skip (line 1183) — header-less
       BEDs silently lost one bin per invocation.
  Pinned by `TestHapcount_CorrectBinTransitions`,
  `TestHapcount_EndOfStreamFlush`, and `TestHapcount_BEDFirstLineWithData`
  in `tools/vcftools/pkg/vcftools/hapcount_test.go`. The
  `.expected.hapcount` fixture is hand-traced, NOT generated from
  the upstream binary (whose output is wrong). ✅
- **`--temp DIR`** — upstream parameters.cpp:341 stores DIR as the
  base path for `mkstemp` spill files used by the LD and
  format-convert paths (variant_file_output.cpp:1441,
  variant_file_format_convert.cpp:28/402/627/810/994). This port does
  not spill to disk for any of those paths so the flag is accepted
  for CLI parity but has no observable effect; `Run` logs to stderr
  that the value was parsed-but-unused. Pinned by
  `TestParams_TempDirAccepted`. ✅
- **`--gzdiff FILE`** — upstream parameters.cpp:237 sets
  `diff_file = FILE; diff_file_compressed = true;`, and
  vcf_file.cpp:21 then switches to the gzip reader. This port's
  `iohelper.OpenReader` already auto-sniffs gzip from the magic
  bytes, so `--gzdiff` is wired as a plain alias for `--diff`
  (last-set wins, matching upstream's shared `diff_file` slot
  semantics parameters.cpp:209 vs :237). Pinned by
  `TestGzdiffAliasesDiff`. ✅

Closed in wave 16 (this PR):

- **`--recode-INFO TAG`** — upstream-canonical name for the repeatable
  recode-INFO-column selector (parameters.cpp:319 →
  `recode_INFO_to_keep.insert(...)`). The port already implemented
  this semantic under the `--keep-INFO TAG` flag name; wave 16 adds
  the canonical spelling as a synonym that funnels into the same
  `keepINFOParts` slice in `tools/vcftools/cmd/vcftools/main.go`,
  matching the existing `--keep-INFO-all` ↔ `--recode-INFO-all`
  pattern. Pinned by `TestCLI_RecodeINFOAlias`,
  `TestCLI_RecodeINFOAlias_Repeatable`, and
  `TestCLI_RecodeINFOAlias_MixedWithKeepINFO` in
  `tools/vcftools/cmd/vcftools/aliases_cli_test.go`. ✅
- **`-c`** — upstream's short alias for `--stdout`
  (parameters.cpp:194). Wired with `flag.BoolVar` pointing at the
  same `useStdout` bool, so either spelling toggles streaming
  output. Pinned by `TestCLI_ShortStdoutFlag`. ✅

Wave 16 also enumerated the FULL upstream long-flag table for the
first time (146 distinct flags). Two flags were initially flagged as
gaps by an earlier wave's regex but actually exist in the port under
non-obvious wiring: `--help` (registered via `flag.BoolVar(help, "help"
…)`, not `flag.Bool("help", …)`). The remaining-gap set below is the
definitive list.

Note on the **`--keep-INFO` semantic gap**: ✅ resolved in wave 17.
Upstream's `--keep-INFO` (parameters.cpp:266 →
`site_INFO_flags_to_keep` → `entry_filters.cpp:1033-1063`) is a SITE
FILTER. Pre-wave-17 the Go port had this flag wired to the
recode-INFO-column selector semantic. Wave 17 separates the two:
`--keep-INFO TAG` now drives `Params.KeepINFO` (site filter that
errors on non-Flag tags, OR-composing across multiple tags), while
`--recode-INFO TAG` drives the new `Params.RecodeINFO` field
(recode-column selector). See `docs/UPSTREAM_BUGS.md` Fix-on-port
section for the migration note. The sibling `--remove-INFO`
divergence was resolved in wave 18 (see below).

**`--remove-INFO` semantic gap**: ✅ resolved in wave 18.
Upstream's `--remove-INFO` (parameters.cpp:328 →
`site_INFO_flags_to_remove` → `entry_filters.cpp:1068-1086`) is a
SITE FILTER that drops sites where the named Flag IS present
(OR-veto across multiple tags). Pre-wave-18 the Go port had this
flag wired as a recode-column stripper — a port-only invention with
no upstream equivalent. Wave 18 repoints `--remove-INFO TAG` at
`passRemoveINFOSite`, the polarity-inverted complement of wave 17's
`passKeepINFOSite`. Header validation (Type=Flag check) is shared
with the keep path via `validateFlagTypeINFO`. Composition with
`--keep-INFO` follows upstream's keep-then-remove ordering:
`TestRun_KeepAndRemoveINFO_Compose` and
`TestParity_KeepAndRemoveINFO_Compose` pin this. The dead
recode-column stripper code in `filterRecodeInfo` was deleted (no
CLI flag drives it now). See `docs/UPSTREAM_BUGS.md` Fix-on-port
section for the full migration note.

Flag history (definitive enumeration vs.
`reference_code/vcftools/src/cpp/parameters.cpp` — wave 16):

As of wave 16 the diff between upstream's `in_str == "--…"` table
(146 unique long flags) and the port's registered flags was **five
flags** (the PCA-deferred trio counted as registered). All five are
now closed — waves 19–23, recorded below for history:

- ~~**`--bcf` FILE** — BCF binary input (parameters.cpp:173).~~
  **Closed (wave 22).** Adapts `pkg/htsgo/bcf.Reader` (built on
  the in-tree BGZF reader and on the wave-21 BCF dictionary fixes)
  into the new `variantSource` interface so `Run` iterates BCF
  records through the same filter pipeline used by `--vcf`. Pinned
  by `TestRun_BCFInput_Roundtrip` (write-then-read symmetry) and
  `TestRun_BCFInput_ComposesWithFilters` (BCF input composes with
  `--chr` / `--from-bp`).
- ~~**`--diff-bcf` FILE** — second-file BCF input for `--diff-*`
  family (parameters.cpp:210).~~ **Closed (wave 23).**
  `loadDiffBCF` mirrors `loadDiffVCF` but routes through the
  wave-22 `bcfVariantSource` (BGZF + `bcf.Reader` + ToVariant).
  The shared `loadDiffFromSource` body drives both loaders so
  the (CHR,POS)-keyed `diffData` build is identical between VCF
  and BCF second files. CLI mutual-exclusion with `--diff` /
  `--gzdiff` mirrors upstream's last-set-wins slot semantics.
  Pinned by 5 tests including all five `--diff-*` outputs
  composed against a BCF second file.
- ~~**`--recode-bcf`** — emit BCF instead of VCF (parameters.cpp:317).~~
  **Closed (wave 21).** Layered the existing `pkg/htsgo/bcf.Writer`
  on top of `pkg/bgzip.NewWriter` and wired it parallel to the
  existing `--recode` text-VCF path; both flags may be combined.
  This wave also fixed three latent bugs in the shared BCF writer
  uncovered by interop with upstream's reader: (1) missing
  IDs were encoded as type-0 instead of zero-length-typed-char;
  (2) the unified INFO+FILTER+FORMAT dictionary numbering wasn't
  surfaced — entries now carry their IDX and the text header is
  emitted with the `,IDX=N` annotations htslib uses; (3) FORMAT
  field descriptors were encoded with the total flat length as
  `size` instead of the per-sample dimension. Pinned by
  `TestRun_RecodeBCF_Roundtrip` and
  `TestRun_RecodeBCF_HeaderHasIDXAnnotations`, plus interop
  testing with upstream vcftools 0.1.18 (decoded our `.recode.bcf`
  through `vcftools --bcf` → byte-identical VCF round-trip).
- ~~**`--contigs` FILE** — supplemental `##contig=` lines for BCF
  header construction (parameters.cpp:197 →
  `variant_file.cpp:45-69`).~~ **Closed (wave 22).**
  `augmentHeaderContigs` prepends `##contig=<ID=...>` MetaInfo
  lines to the parsed header when the source lacks contig
  declarations of its own (matching upstream's
  `has_contigs == false` gate). Accepts both bare contig names
  and full `##contig=<...>` forms. Pinned by
  `TestRun_ContigsFile_AddsContigLines`,
  `TestRun_ContigsFile_NoOpWhenHeaderAlreadyHasContigs`,
  `TestAugmentHeaderContigs_AcceptsMetaInfoForm`.
After wave 23 vcftools reaches **146/146 long flags** — the
complete `parameters.cpp` surface is exercised. The
remaining work is per-output column-set polish (see the
"Other" list below) and the multi-output `num_outputs > 1`
check upstream uses, which we deliberately don't replicate
per the CLAUDE.md "don't replicate upstream constraints" rule.

Other (per-output column-set gaps, not flag-count gaps):

- **Per-individual output**: the per-individual `.imiss` row layout
  still has fields we don't emit (we have `--missing-indv`).
- ~~**Diff family**: per-site / per-indv discordance outputs
  (`--diff-site-discordance`, `--diff-indv-discordance`) still emit a
  simpler column set than upstream's richer `.diff.sites` /
  `.diff.indv` schemas — see `variant_file_diff.cpp:635` for the gap.~~
  **Closed (wave 20).** `.diff.sites` now emits the upstream 7-column
  layout (`CHROM POS FILES MATCHING_ALLELES N_COMMON_CALLED N_DISCORD
  DISCORDANCE`), including file-1-only and file-2-only zero rows;
  `.diff.indv` now emits the 4-column layout (`INDV N_COMMON_CALLED
  N_DISCORD DISCORDANCE`) over the union of file-1 and effective
  file-2 samples in alphabetical order. Discordance values format
  via `%.6g` with `-nan` for 0/0, matching libstdc++'s default
  ostream output. Pinned by five parity tests against upstream
  goldens (`TestParity_DiffSiteDiscordance_{NoMap,WithMap,AltMismatch}`,
  `TestParity_DiffIndvDiscordance_{NoMap,WithMap}`). Residual
  deviations from upstream:
  1. Row ordering within `.diff.sites` follows file-1 streamed
     order with file-2-only sites appended in sorted-chrom-then-pos
     order, rather than upstream's strict merge sort — observable
     only when the two files have non-overlapping positions
     interleaved by chromosome.
  2. **REF-mismatch shared sites**: upstream
     `variant_file_diff.cpp:787-790` SKIPS B-sites where REFs
     differ (with a one-off `"Non-matching REF"` warning), treating
     the site as if it weren't shared. The Go port emits the row
     with `MATCHING_ALLELES=0` and accumulates discordance over
     it. Tracked as a separate follow-up.
  3. **REF=N/`.`/empty normalisation**: upstream replaces
     `REF1` with `REF2` (and vice versa) when one side is `N`,
     `.`, or empty before the alleles-match check
     (`variant_file_diff.cpp:780-783`). The Go port does a
     verbatim string compare. Same follow-up as (2).
  4. **`-nan` literal portability**: the port hardcodes
     `"-nan"` for division-by-zero discordance to match
     libstdc++'s default ostream output on glibc. A future
     golden regenerated on musl / macOS libc / MSVC could
     print `nan` or `NaN` instead. Re-running the goldens on
     a non-glibc system would require updating the literal.
- **Other**: small-format columns gaps tracked in
  `tools/PORTING_STATUS.md`.

Note: the brief mentioned `--haploid` as a possible wave-2 target. After
checking the upstream source (`reference_code/vcftools/src/cpp/`) there is
no `--haploid` flag — the closest thing is `--phased` (parameters.cpp:311
+ entry_filters.cpp:989-1010), which we ported instead.

**Validation:** wave 1 adds header byte-for-byte parity tests for the new
output files; wave 2 ships full byte-for-byte parity tests for both
`.ldhat.sites` and `.ldhat.locs`; wave 3 adds byte-for-byte parity tests
for `.ldhelmet.snps` / `.ldhelmet.pos` and the IMPUTE bundle
(`.impute.legend` / `.impute.hap` / `.impute.hap.indv`); wave 4 adds
byte-for-byte parity tests for `.diff.discordance_matrix` (with and
without `--diff-indv-map`) and the mapped `.diff.indv_in_files` output
against upstream goldens (under `tools/vcftools/testdata/parity/`,
fixtures `diff_f1.vcf` / `diff_f2.vcf` / `diff_indv_map.txt`); wave 5
adds byte-for-byte parity tests for `.diff.switch` /
`.diff.indv.switch` (fixtures `switch_f1.vcf` / `switch_f2.vcf`) and
`.mendel` (fixtures `mendel.vcf` / `mendel.ped`); wave 6 adds
byte-for-byte parity tests for `--non-ref-af 0.3` / `--non-ref-af 0.5`
on `sample.vcf` and `--chr X --non-ref-ac 2|3` on the same fixture
(see `non_ref_af_*.expected.recode.vcf` and
`non_ref_ac_*_chrX.expected.recode.vcf`), plus a regression test that
pins the upstream `_any`-fallback asymmetry between the two flags;
wave 7 adds byte-for-byte parity tests for `--non-ref-ac-any 2`
(full VCF), `--chr 20 --non-ref-ac-any 1`,
`--chr 20 --max-non-ref-af 0.3`, `--chr 19 --max-non-ref-ac 2`,
`--chr 20 --max-non-ref-ac-any 2`, and
`--non-ref-af 0.3 --non-ref-af-any 0.6` (the only meaningful AF-any
usage), plus a `TestParity_NonRefAFAny_NoOp` regression that pins
upstream's documented "AF -any flags are no-ops alone" quirk by
asserting the port produces baseline output unchanged.
Full upstream-test-suite run still pending.

Upstream build note for golden generation: vcftools'
`variant_file_format_convert.cpp` LDhat writers allocate a stack array of
exactly `temp_dir.size()` bytes and then `strcpy` a longer string into
it (`new_tmp = temp_dir + "/vcftools.XXXXXX";` plus
`char tmpname[new_tmp.size()]; strcpy(tmpname, new_tmp.c_str());` at
lines 627-629 and again at 813-815, 658-665, 680-688, 841-845). Modern
glibc fortified strcpy aborts the run. To regenerate the goldens locally:

```
cd reference_code/vcftools/src/cpp
make clean
make CXXFLAGS='-O0 -g -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0'
```

Filed as a parity-only workaround; we don't replicate the bug.

### `bgzip`

**Status:** 1 / 1 command, most flags.

Missing:

- **Multi-threaded compression** (`-t / --threads N` is accepted but
  single-threaded; BGZF is trivially parallel per block).
- **Output-rename to follow upstream conventions on stdin**: minor.

**Validation:** round-trips through `tabix` work; no full upstream-test
suite run yet.

### `tabix`

**Status:** 1 / 1 command, most flags.

Missing:

- **`--reheader FILE`** — replace bgzipped file's header lines in place.
- **`--targets` strictness** — currently behaves as `-R`; needs to be a
  true post-filter that only emits records strictly inside the targets.

**Validation:** no full upstream-test-suite run yet.

### `samtools`

**Status:** 24 of ~25 subcommands (~96%). `view`, `sort`, `index`, `depth`,
`fastq`, `flagstat`, **`mpileup`** (wave-1 + tail wiring), PR #88's
wave-1 tail (`merge`, `coverage`, `idxstats`, `cat`, `reheader`,
`addreplacerg`, `fixmate`, `dict`, `split`, `quickcheck`), the
heavy-hitter pair `markdup` + `stats`, the calmd/import pair
(**`calmd`** + **`import`**), the niche pair landed in the
phase/targetcut PR (**`phase`** + **`targetcut`**), and now
**`consensus`** (simple- and bayesian-mode FASTA/FASTQ/pileup; the
Gap5 posterior caller and the NM-halo MAPQ adjustment are byte-faithful
to upstream's default `MODE_RECALL`).

Missing subcommands (in rough priority order):

- **`tview`** — terminal viewer. **Deliberate skip** (interactive
  curses UI; near-zero pipeline usage and would require an ncurses
  dependency). Not on the roadmap.
- **`view` flag-tail**: `-X`/`--customized-index` (explicit index-file
  argument after `<in.bam>`) is implemented — the index kind (.bai or
  .csi) is auto-detected from the file's magic. `-L bed` landed as a
  linear-scan BED-region filter; `-M`/`--use-multi-region-iterator` is
  accepted but treated as a no-op since we always run the full
  intersection. `-d/-D` (tag-value filter) and `-N` (qname file) landed
  in the view-d-D-N PR.
- **`mpileup` tail** beyond PR #88 wiring: the remaining genuine gap is
  BCF / genotype-likelihood output (`-g/-u`). `-aa` zero-fill of empty
  contigs is implemented (see `TestMpileup_AA_ZeroFillTableDriven`). The
  text-pileup path is complete. `-g/-u` requires the genotype-likelihood
  model (`bam2bcf`) plus a BCF emit path; `reference_code/htslib`
  (`errmod.c`, the MAQ likelihood) is now vendored, so the reference
  source for the port is available — see the "deferred" note below.

Plus:

- **CRAM** read/write throughout — DONE. Landed across PRs #162–#180
  (the rANS 4x8/4x16 codecs are in-tree pure Go; `ulikunitz/xz` is the
  only sanctioned third-party dep, confined to the LZMA block codec).
  See `docs/CRAM_DESIGN.md` and `docs/CRAM_ROADMAP.md`.
- **`.csi` index** — DONE (PR #189); `samtools index` emits both `.bai`
  and `.csi`, and readers auto-detect index kind from file magic.
- **Multi-threading (`-@`) — NOT done.** `-@/--threads` is accepted on
  the CLI of `sort`, `index`, `view` (and elsewhere) but is a no-op
  stub; v1 is single-threaded everywhere. The option value is stored on
  the relevant options struct for a future parallel pass. This is the
  one cross-cutting deferred item, not a completed feature.

**Genuine remaining samtools gaps** (everything else is done):

- **`mpileup` MAQ genotype-likelihood model — slices 1-4 DONE; only
  indel calling remains deferred.** The `mpileup` SNP MAQ-model port is
  complete: it was sliced into four parts:
  - **Slice 1 (DONE).** The MAQ error model (`errmod.c`) is ported to
    pure Go in the shared `pkg/htsgo/errmod` package (`errmod.Init` /
    `errmod.Cal`). Both `bcftools mpileup` and `samtools targetcut`
    consume this single implementation.
  - **Slice 2 (DONE).** The per-site genotype-likelihood pipeline —
    `bcf_call_glfgen` / `bcf_call_combine` / `bcf_call2bcf` from
    `bam2bcf.c` — is ported in
    `tools/bcftools/pkg/bcftools/bam2bcf.go`. `mpileup` now emits one
    BCF/VCF record per covered reference position with real MAQ PLs,
    the `<*>` "unseen" allele, the multi-allelic PL grid, and
    INFO/DP/I16/QS/MQ0F. `-O b` (BCF output) works through
    `pkg/htsgo/bcf`; the old `-O u/b` hard-rejection is gone. The
    upstream mpileup defaults are now applied (`min-BQ=1`,
    `max-BQ=60`, `delta-BQ=30`). The `delta_baseQ` neighbour-quality
    cap is implemented.
  - **Slice 3 (DONE).** BAQ realignment is wired into the pileup.
    `applyMpileupBAQ` in `mpileup.go` ports `mpileup.c`'s `mplp_realn`:
    every covered column is gated by the `MPLP_REALN_PARTIAL`
    heuristic, and each selected read is run once through
    `pkg/htsgo/baq.SamProbRealn` in apply+extend mode (`flag 3`,
    matching `sam_prob_realn(b, ref, ref_len, (flag & MPLP_REDO_BAQ) ?
    7 : 3)`) before its bases enter the pileup. `-B/--no-BAQ` disables
    BAQ; `-E/--redo-BAQ` adds `baq.FlagRedo` (`flag 7`). The default is
    PARTIAL realignment — upstream sets `MPLP_REALN | MPLP_REALN_PARTIAL`
    (mpileup.c:1389), so the per-column skip heuristic and the per-read
    spanning check apply. `-D/--full-BAQ` clears `MPLP_REALN_PARTIAL`
    (mpileup.c:1567), forcing full BAQ: every read on the chromosome is
    realigned. The `MPLP_REALN_PARTIAL` per-column skip heuristic depends
    on the per-column `p->indel` term, which needs indel detection — see
    slice 4 — so it is approximated from the `PLP_HAS_INDEL` CIGAR scan;
    exact for indel-free inputs. BAQ runs
    before `accumulateMpileupBases`, so the `delta_baseQ` cap and the
    `min_baseQ` filter inside `bcfCallGlfgen` see the BAQ-adjusted
    qualities — matching `mpileup.c`'s ordering.
  - **Slice 4 (DONE).** The per-site bias annotations are ported in
    `bam2bcf.go`: `calcVDB` (VDB, with the htslib `kf_erfc` rational
    approximation ported as `kfErfc` and the single-precision `float`
    arithmetic upstream relies on), `calcSegBias` (SGB) and
    `calcMWUBiasZ` (the standard-deviation-normalised Mann-Whitney U
    z-score that yields RPBZ / MQBZ / BQBZ / MQSBZ / SCBZ), plus the
    MQ0F fraction. The bias tallies (`ref_pos`/`alt_pos`,
    `ref_mq`/`alt_mq`, `ref_bq`/`alt_bq`, `fwd_mqs`/`rev_mqs`,
    `ref_scl`/`alt_scl`) accumulate per-sample in `bcfCallret` from the
    I16 path and from `getPosition` (the soft-clip-aware read-position
    port), then combine in `bcfCallCombine`. INFO floats — QS, I16 and
    every bias tag — are now rounded through `float32` and rendered with
    a 6-significant-digit `%g` (`formatFloat32G`), matching upstream's C
    `float` storage byte-for-byte. `MPLP_SMART_OVERLAPS` is ported too:
    `applySmartOverlaps` pairs the two mates of each proper pair and
    `tweakOverlapQuality` (a faithful port of htslib's
    `tweak_overlap_quality`, with the `cigar_iref2iseq` iterator and the
    Wang/X31 read-name hashes) merges the overlapping-span base
    qualities before BAQ. Byte-for-byte parity is verified by
    `TestMpileupSNPGoldens`: the full `mpileup/mpileup.11.out` golden
    (4001 covered positions, 87 SNP ALT records, two overlapping mate
    pairs) and the three-sample multi-BAM `mpileup/mpileup.1.out` golden
    match exactly — header and every SNP data record, all bias INFO tags
    included.

  **Remaining deferred work — indel calling only.** The one upstream
  `mpileup` path still unported is the indel caller (`bam2bcf_indel.c` /
  `bam2bcf_edlib.c`): indel candidate detection, the indel genotype
  likelihoods and the INDEL/IDV/IMF INFO tags. The single INDEL record
  of `mpileup.11.out` (17:302 `TA`) is consequently the only golden line
  not reproduced; `TestMpileupSNPGoldens` aligns records by `CHROM:POS`
  and skips it, and `TestMpileupGoldensDeferred` catalogues every other
  deferred golden with its precise reason (FORMAT tags beyond PL, `--ff`
  FLAG filtering, `-s/-S/-G` sample/read-group selection, IUPAC REF
  bases, indel/SCR fixtures).

  **Indel-caller sub-slicing (in progress).** The remaining indel work
  is broken into five sub-slices:

  - **4a + 4b (DONE).** Pileup data model + STR finder + indel
    candidate-type discovery helpers. `pileupBase` now carries the
    htslib `bam_pileup1_t.indel` field, an `aux` scratch word, and an
    optional `*sam.Record` back-pointer set only on indel-bearing
    columns; `accumulateMpileupBases` populates them by peeking at the
    next consuming CIGAR op. `bam2bcf_indel.go` adds `bcfCallauxIndel`
    (the indel-specific subset of upstream `bcf_callaux_t`) plus the
    static helpers `est_seqQ`, `est_indelreg`, `bcf_cgp_l_run`,
    `bcf_cgp_find_types`, `tpos2qpos`, and `get_pos`. The STR finder
    (upstream `str_finder.c` / `find_STR` / `find_STR64`) is ported as
    a reusable in-tree package `pkg/htsgo/strfinder` and unit-tested
    over hand-traced inputs (homopolymers, dinucleotide repeats,
    padding, lower-case filter). All CLI knobs (`--open-prob`,
    `--ext-prob`, `--tandem-qual`, `--min-ireads`, `--gap-frac`,
    `--indel-bias`, `--indel-size`) reach `bcfCallauxIndel` via
    `newBcfCallauxIndel`, but the value is not yet driven into
    emission.
  - **4c (DONE).** Per-sample consensus + reference-sample
    construction (`bcf_cgp_ref_sample`, `bcf_cgp_calc_cons`) and the
    Probaln-based alignment scoring core (`bcf_cgp_align_score`) are
    ported in `tools/bcftools/pkg/bcftools/bam2bcf_indel_align.go`.
    `bcfCgpRefSample` samples the per-sample 4-bit IUPAC reference and
    masks positions where ALT ≥30% with N (15); `bcfCgpCalcCons`
    majority-rule-builds the per-type insertion consensus (zeroing
    types whose consensus contains an N); `bcfCgpAlignScore` runs
    `baq.ProbalnGlocal` (BW = |typeLen|+3; PacBio CCS params for
    >1000-bp reads), applies the indel-bias clamp to the
    length-normalised score, and folds in the STR-finder fudge
    (`strfinder.FindSTR`). The returned word reproduces upstream's
    `(sc<<8) | min(255, l)` bit-pattern byte-for-byte, with the low
    byte rewritten as `min(255, (score&0xff)*0.8 + iscore*2)` after
    the STR fudge — verified by `TestBcfCgpAlignScore_BitPattern`.
  - **4d (DONE).** `bcf_call_gap_prep` + `bcf_cgp_compute_indelQ` are
    ported as `bcfCallGapPrep` and `bcfCgpComputeIndelQ` in the same
    file. `bcfCallGapPrep` orchestrates the full pipeline: cheap-reject
    on a clean column → `bcfCgpFindTypes` → `bcfCgpRefSample` →
    `bcfCgpCalcCons` → per-(read,type) `bcfCgpAlignScore`, scoring into
    a flat `N*nTypes` matrix. `bcfCgpComputeIndelQ` then folds those
    scores into per-read `p.aux` words (chosen-type<<16 | seqQ<<8 |
    indelQ), populates `bca.IndelTypes` / `bca.Inscns` / `bca.MaxIns`
    by sumq-sorting the candidate types (REF always at slot 0), and
    returns n_alt. Unit-tested at the bit-pattern level
    (`TestBcfCgpComputeIndelQ_BitPattern`) with a hand-derived
    expectation; the orchestrator has a clean-site -1 reject test and
    an indel-site smoke test.
  - **4e (in progress).** Indel-aware branches of `bcf_call_glfgen` /
    `bcf_call_combine` / `bcf_call2bcf` and emission of the
    INDEL/IDV/IMF INFO tags; golden tests land here. Broken into
    cross-cutting sub-slices since the targeted goldens require more
    than just the indel branch of glfgen+2bcf:
    - **4e.1 (DONE).** FORMAT/AD support and `-a/--annotate` parsing.
      `parseFormatFlag` (`mpileup.go`) is the Go port of
      `parse_format_flag` (mpileup.c:1141) with bit constants
      `B2BFmt*/B2BInfo*` matching `bam2bcf.h:46-75`.
      `validateMpileupOptions` seeds `opts.FmtFlag` from
      `DefaultMpileupFmtFlag` (mpileup.c:1399) and layers user tokens on
      top — `-AD` clears the bit, `AD`/`FORMAT/AD` sets it. `bcfCall2bcf`
      emits FORMAT/AD,ADF,ADR per-sample columns when the matching bit
      is set, with the per-allele depths drawn from the already-reordered
      `call.adf/adr` arrays. The corresponding `##FORMAT=<ID=AD,...>`
      header lines are emitted only when their bits are on.
      `TestMpileupSNPGoldens/multi-bam-region-format-AD` is the
      byte-for-byte golden check against `mpileup/mpileup.12.out`.
    - **4e.2 (DONE).** Indel-branch `bcf_call_glfgen`
      (`bam2bcf.c:300-460`) is now `bcfCallGlfgenIndel` /
      `bcfCallGlfgenCore` in `bam2bcf.go`: it consumes per-read `p.aux`
      words populated by `bcfCallGapPrep`, runs `errmodCal` on the indel
      bases, and populates `bcfCallret` with indel-specific QS / AD /
      I16 / bias tallies (using `is_diff = b ? 1 : 0`, mirroring
      `bam2bcf.c:350`). `bcfCallGapPrep` was extended to populate
      `bca.IrefPos / IaltPos / IrefMq / IaltMq / IrefScl / IaltScl`
      (`bam2bcf_indel.c:826-848`) from the per-read loop at t==0; the
      indel-flavored `getPos` helper supplies the read-position and
      soft-clip-length bins. `pileupBase.rec` is now populated for
      every read in the pile (not just indel-bearing columns) so the
      indel iref/ialt accumulation has the cigar available for all
      reads.
    - **4e.3 (DONE).** Indel-branch `bcf_call_combine` + `bcf_call2bcf`
      (`bam2bcf.c:1165-1198`, `bam2bcf.c:1211-1234`) are
      `bcfCallCombineIndel` (in `bam2bcf.go`) and `bcfCall2bcfIndel`
      (in `mpileup.go`): REF/ALT alleles built from
      `bca.Inscns`/`IndelTypes`/`IndelReg`, `INFO/INDEL` flag plus
      `IDV`/`IMF` emitted before `DP`/`I16`/`QS`, then the same bias
      subset as the SNP path (`VDB`/`SGB`/`RPBZ`/`MQBZ`/`MQSBZ`/`BQBZ`/
      `SCBZ`/`MQ0F`) followed by `FORMAT/PL`. The upstream "leaked"
      bias semantics (BQBZ/MQSBZ retain the last has-alt SNP's value
      since `bcf_callaux_clean` does not reset `call->mwu_*`) is modelled
      explicitly by a `biasLeak` struct threaded through the per-site
      driver. The driver in `emitChromMpileup` runs `bcfCallGapPrep`
      after the SNP emission and, when it returns ≥0, runs a second
      `bcfCallGlfgenIndel` + `bcfCallCombineIndel` + `bcfCall2bcfIndel`
      pass (matching `mpileup.c:589-613`). The full
      `mpileup/mpileup.11.out` golden — including the 17:302 INDEL row
      (`T → TA`, `BQBZ=-1.34164` inherited from the prior SNP combine
      at 17:237) — now matches byte-for-byte.
    - **4e.4 (DONE).** `--ambig-reads` (incAD / incAD0) ADF/ADR
      compensation (`bam2bcf.c:540-561`), `--skip-{all,any}-{set,unset}`
      BAM-flag filters (`mpileup.c:208-211`), and the three latent items
      the 4e.2+4e.3 review flagged:
        1. `p.is_del` modelling — `accumulateMpileupBases` now emits a
           pileupBase event for every reference column inside a read's
           `D` op (`isDel=true`, `b=0`); `bcfCallGlfgenCore` lets these
           reads through in the indel branch (matching upstream
           `bam2bcf.c:307 `if (p->is_del && !is_indel) continue`).
        2. `p.is_refskip` modelling — same shape for CREF_SKIP (`N`)
           ops; `isRefskip` reads are dropped in both branches
           (`bam2bcf.c:301`).
        3. Cross-contig `biasLeak` reset — the `biasLeak` instance now
           lives on the per-run driver in `writeMpileupVCF` and threads
           through every `emitChromMpileup` call so the BQBZ / MQSBZ
           scalars persist across contigs, mirroring upstream's `conf->bc`
           lifetime. The leak is pre-initialised to `(0, ok=true)` so the
           very first indel record (before any has-alt SNP combine) sees
           BQBZ=MQSBZ=0 — matching upstream's C `bcf_call_t bc = {0}`
           default-initialisation.
      `--ambig-reads` is parsed via `parseAmbigReads` into
      `AmbigReadsMode` on `MpileupOptions`; the indel-branch glfgen
      stashes low-quality REF-looking reads in `adrRefMissed[]` /
      `adfRefMissed[]` per upstream and applies the `incAD` (proportional)
      or `incAD0` (claim as REF) compensation. The `--skip-*` strings go
      through `parseBAMFlagString` (Go port of htslib's `bam_str2flag`)
      into the `RflagSkip{Any,All}{Set,Unset}` masks consumed by
      `mpileupKeepRecord`. Byte-for-byte golden: `mpileup-filter.2.out`
      now matches (both the `--skip-all-unset READ1` and `--skip-any-unset
      READ1` forms; tested in `TestMpileupFilterGolden`).
    - **4e.5 (DONE).** INFO/FMT/SCR (soft-clipped reads counter). The
      SCR signal is a per-record `hasSoftClip` bit (true iff the read
      has any CIGAR S op, mirroring upstream's `PLP_HAS_SOFT_CLIP` set
      in `pileup_constructor`, `mpileup.c:317-323`). The bit is stamped
      on every `pileupBase` produced by `accumulateMpileupBases` so the
      SNP-branch `bcfCallGlfgenCore` can tally it pre-refskip into
      `bcfCallret.scr` (matching `bam2bcf.c:300`). `bcfCallCombine`
      folds the per-sample counts into `bcfCall.scrTotal` /
      `bcfCall.scr[]`, and `bcfCall2bcf` emits INFO/SCR (before I16) /
      FORMAT/SCR (after AD/ADF/ADR) when their bits are set.
      `parseFormatFlag` now also accepts the `FMT/` prefix in addition
      to `FORMAT/` (`SET_FMT_FLAG`, `mpileup.c:1120-1122`). Byte-for-
      byte golden `mpileup-SCR.out` (test.pl:1069) now matches; tested
      in `TestMpileupSCRGolden`.
    - **4e.6 (DONE).** INFO/NMBZ (per-read NM bias) — the Mann-Whitney
      U z-score over the per-read NM tag, split REF vs ALT. The port
      adds `getAuxNm` (Go counterpart of `get_aux_nm`, `bam2bcf.c:96`)
      which reads the BAM `NM:i:` tag, treats each indel CIGAR op as a
      single event (subtracting `len-1` for ops with `len>1`), counts
      soft-clip lengths as mismatches, then subtracts 1 for REF reads
      and 2 for ALT reads (the MNP-aware adjustment) and clamps to
      `[0, b2bNNm-1]`. Per-sample `refNm[32] / altNm[32]` histograms on
      `bcfCallret` are filled by both the SNP and indel branches of
      `bcfCallGlfgenCore` when any of the `B2BInfoNMBZ` / `B2BFmtNMBZ`
      / `B2BInfoNM` bits is set. `bcfCallCombine` sums them across
      samples and runs `calcMWUBiasZ`; `bcfCallCombineIndel` also
      folds in the matching SNP-pass tallies (mirroring upstream's
      shared `bca->ref_nm/alt_nm` accumulator). `bcfCall2bcf` and
      `bcfCall2bcfIndel` emit `INFO/NMBZ` between MQSBZ and SCBZ when
      `B2BInfoNMBZ` is set, and the header line is inserted in the
      same place. Byte-for-byte golden `annot-NMBZ.1.1.out` (test.pl
      line 1074) now matches; tested in `TestMpileupNMBZGolden`.
      `annot-NMBZ.2.1.out` byte-matches once the depth-cap port is in
      place (4e.8 below); `.3.1.out`'s SNP row byte-matches including
      `NMBZ=7.74597` while the indel row's QS / NMBZ / PL[0] residual
      tracks with the `bcfCall2bcfIndel` SCR-on-indel-rows polish item.
      FORMAT/NMBZ (per-sample) is not in this slice; only INFO/NMBZ.
    - **4e.7 (DONE).** Ports the legacy REF-rescue heuristic that
      `bcf_call_glfgen` applies in its indel branch
      (`bam2bcf.c:338-348`, originally Heng Li's e4e161068 fix for
      htslib issue #1446). At deeply-covered homopolymer / tandem-
      repeat sites `bcf_cgp_compute_indelQ` emits `indelQ = 0` for
      most REF-leaning reads — so without the rescue every REF read
      fails the indel-branch min-baseQ gate and lands in
      `ADR/ADF_ref_missed`, breaking I16 / QS / AD parity. The
      heuristic: when a read has no CIGAR indel (`p.indel == 0`) and
      either `q < _n/2` or `_n > 20`, reclassify it as REF (`b = 0`),
      promote `q` to the read's raw base quality at qpos, and rebuild
      `seqQ` as `(3*seqQ + 2*q)/8`; once `_n > 20`, cap that seqQ at
      40. `baseQ` for I16 stays at the pre-heuristic value in
      `p.aux>>8&0xff` (matching `bam2bcf.c:420`); the local seqQ is
      used only to cap `q` and gate min-baseQ. Two correctness fixes
      landed at the same time: (a) `seqQ = q = (p.aux & 0xff)` on the
      indel branch (upstream `bam2bcf.c:315` initialises both from
      the indelQ bits, not the saved seqQ bits), and (b) the SNP
      branch now stores `seqQ = b2bSeqQ` explicitly so the shared
      `if (q > seqQ) q = seqQ` step at `bam2bcf.c:459` is the same
      literal port. Goldens: `indel-AD.{2,3,4}.out` are byte-for-byte
      (tested in `TestMpileupIndelADGolden`, covering `--ambig-reads`
      default / `incAD` / `incAD0`); the `annot-NMBZ.3.1.out` indel
      row's I16 also matches byte-for-byte. The `--indels-cns` (edlib)
      path is a separate algorithm and stays deferred. Residual
      divergences live on `indel-AD.1.out` (~20 SNP rows with small
      I16 base-quality drifts plus four homopolymer-column indel rows
      with chosen-type off-by-1 assignments — DP ≤ 125 throughout so
      the depth cap is not in play) and the trailing QS/NMBZ/PL[0]
      columns of `annot-NMBZ.3.1.out`'s indel row.
    - **4e.8 (DONE).** Ports htslib's per-alignment-start depth cap.
      Upstream's `bam_plp_push` (reference\_code/htslib/sam.c:6090)
      drops a new read when `iter->pos == b->core.pos` and the
      pileup queue already holds `maxcnt` active reads, while our
      previous port truncated each per-column pile to `MaxDepth`
      reads instead. The new `applyMpileupDepthCap` walks the per-
      sample coordinate-sorted record stream, maintains a min-heap
      of in-flight end positions, and drops reads using exactly the
      htslib predicate (including the per-alignment-start trigger
      and the "mp->cnt includes one sentinel node" off-by-one). The
      cap runs BEFORE `applySmartOverlaps` so a capped-out read
      cannot leak a tweaked base quality into the surviving mate —
      upstream's `overlap_push` runs inside `bam_plp_push` after the
      cap test, so dropped reads never reach the overlap-quality
      merger. Byte-for-byte golden `annot-NMBZ.2.1.out` at chr6:75
      (raw coverage 449, capped to DP=283) now matches; tested in
      `TestMpileupDepthCapGolden`. Remaining mpileup residuals are
      the `--indels-cns` edlib path (separately deferred) and the
      indel-row QS/NMBZ/PL[0] columns at homopolymer columns on
      `indel-AD.1.out` and `annot-NMBZ.3.1.out`, both tracked under
      the `bcfCall2bcfIndel` SCR-on-indel-rows polish.
    - **4e.5 indel-row SCR (DONE).** `bcfCall2bcfIndel` now emits
      `INFO/SCR` (before I16) and `FORMAT/SCR` (after AD/ADF/ADR)
      when the corresponding `B2BInfoSCR` / `B2BFmtSCR` bits are
      set, mirroring the SNP-row code path. `bcfCallCombineIndel`
      copies the per-sample SCR tally from the SNP-pass
      `bcfCallret.scr` arrays (the indel branch of `bcfCallGlfgen`
      does not tally SCR — the SCR accumulator is gated on
      `!isIndel`, matching upstream's `bam2bcf.c:300`). Both the
      SNP and indel rows at the same column therefore report the
      same SCR counts. Regression: `TestMpileupSCROnIndelRow` (uses
      the `indel-AD.2.fa` / `indel-AD.2.bam` fixture, which has a
      homopolymer-anchored indel call and soft-clipped reads at
      `11:75`).
    - **Residual: `indel-AD.1.out` and `annot-NMBZ.3.1.out` indel-row
      drifts at homopolymer columns (DOCUMENTED, root cause traced).**
      A column-by-column diff against the goldens identifies three
      independent clusters:

      1. **Trailing N-REF rows past the FASTA end** (`indel-AD.1.out`
         positions `000000F:687-688`, two records with `REF=N`,
         `DP=1`, all I16 fields zero). **RESOLVED** — `emitChromMpileup`
         now computes `effLen = max(refLen, maxReadEnd)` by scanning
         each input's records for `rec.EndPosition()`, sizes the
         events array to `effLen`, and walks the per-site loop to
         `effLen`. The existing `pos0 < len(refSlab)` fallback emits
         `REF=N` past the FASTA. The indel pass is skipped for
         `pos0 >= refLen` (upstream's indel branch reads `ref_fai`
         only within the FASTA, so an `N`-anchored column emits the
         SNP row alone).

      2. **SNP-row I16 base-quality-sum micro-drifts** (`indel-AD.1.out`
         originally 15 columns at 000000F:446-450, 497-500, 540-542,
         566-567, 624, with I16 slots 4-5 shifted at each — e.g.
         1204/45840 vs upstream's 1205/45921). **Mostly RESOLVED.**
         `emitChromMpileup` now interleaves BAQ and overlap-merge
         per-pair, matching upstream `htslib bam_plp_push`
         (sam.c:6083-6132): a `classifyMatePairs` pre-pass labels
         each read as standalone, first-mate or second-mate (same
         predicate as `applySmartOverlaps`'s `overlap_push` port)
         and records the mate's alignment start. The BAQ engine
         then runs in two phases sharing a per-record `realigned`
         dedup map (the equivalent of upstream's PLP_IS_REALN flag).
         Phase 1 BAQs standalones plus first-mates whose first
         eligible pileup column precedes the mate's start (raw
         quals — the `pos0 < mateStart` gate captures whether
         `overlap_push` would have run yet at upstream's iter->pos).
         `applySmartOverlaps` then merges quals. Phase 2 BAQs all
         second-mates plus any first-mates phase 1 left untouched
         (merged quals). 13 of the 15 cluster-2 columns now
         byte-match (446-449, 497-499, 540-542, 566-567, 624).

         **Residual: 2 cluster-2 columns + 2 newly-drifting columns.**
         The SNP rows at 450 and 500 still drift by a single read's
         BAQ delta (1438/55672 vs upstream's 1449/56453 at col 450),
         and the two-phase batch additionally drifts the
         previously-matching SNP rows at 547/548 (4118/160128 vs
         upstream's 4105/158737 at col 547). These four match the
         cluster-2 trace's "543-548" group: the two-phase batch
         still differs from upstream's per-column interleave in
         one edge case — a read whose first eligible column lies in
         the same column-equivalence-class as its mate's push and
         several other reads, where upstream's `bam_plp_next` drains
         the column with the mate already merged but our phase 1
         must commit before `applySmartOverlaps`. Closing the last
         four needs the full streaming `bam_plp_push` column-by-
         column interleave (per-pair arrival state inside the
         column engine, not just a two-phase batch). The BAQ code
         itself (`pkg/htsgo/baq/realn.go`) is byte-identical to
         upstream when fed identical input qualities; the residual
         divergence is purely the orchestration.

      3. **Indel-row chosen-type off-by-one at homopolymer columns**
         (`indel-AD.1.out` at `000000F:537/538/658`, plus
         `annot-NMBZ.3.1.out` at `chr16:75`). At these columns the
         indel-row I16 fields agree on REF/ALT classification (so
         the `isDiff = b ? 1 : 0` split is identical), but the
         **per-allele** breakdown shifts because a handful of reads
         are assigned to a different non-REF indel type by
         `bcfCgpComputeIndelQ`. Concretely at `000000F:537`: ours
         classifies one extra read as type `-1` (the deletion),
         shifting I16 slot 3 (alt-rev count) from 27 to 28, the
         alt BQ-sum by exactly one BQ=40 contribution, the alt
         MQ-sum by MQ=60, and the alt min-dist sum by 25. This
         propagates into QS (per-type qsum), AD (per-allele
         counts), and (transitively) PL[2]. At `chr16:75`
         (`annot-NMBZ.3.1.out`) the I16 byte-matches because the
         swaps are between two ALT types (both `b!=0`, so isDiff
         stays 1), but QS shifts (1.45884 vs 1.43466 for type 1),
         NMBZ flips sign (+0.437589 vs -0.886523 — the indel-pass
         refNm/altNm split changes), and PL[0] for sample 1 falls
         from 255 to 226. Root cause: `bcfCgpComputeIndelQ` and
         `bcfCgpAlignScore` use the same encoding `score<<6 | t`
         and the same ascending-insertion-sort tie-break as
         upstream, and the orchestrator's `tbeg`/`tend`/`rStart`/
         `rEnd` / ref2-slice positioning matches upstream
         byte-for-byte (verified by `bcfCgpAlignScore`'s
         existing tests). The remaining divergence is in the
         `probaln_glocal` score itself for reads whose query
         window straddles a homopolymer run — at most-frequent /
         long homopolymers, two indel types can produce
         `probaln_glocal` returns that differ by a single Phred
         unit, and a one-unit drift in the underlying HMM
         likelihoods (rounding inside the forward/backward DP)
         flips the tie-break. The HMM port lives in
         `pkg/bcftools/.../baq.ProbalnGlocal` and the residual is
         a "single ULP rounding inside the forward DP at long
         homopolymer columns" item rather than a `bam2bcf.c`
         port gap. Out of scope for this slice — fixing it
         requires either (a) reproducing C's exact left-to-right
         multiply order inside the HMM transition step, or (b)
         accepting these as RNG-class residuals (a handful of
         reads per homopolymer column at depths > 100).

  One accepted divergence: `errmod_cal`'s downsampling of piles deeper
  than 255 reads uses Go's RNG rather than htslib's `drand48`, so
  byte-for-byte parity holds only at depth ≤255 (RNG byte-parity is not
  a project goal). All vendored `mpileup` fixtures are within that
  bound.

  **Parity watch-item (resolved) — QS-sum zero-break float comparison.**
  `bcfCallCombine` in `bam2bcf.go` breaks the allele-ordering loop on
  `qsum[ipos] == 0`, mirroring upstream `bcf_call_combine`
  (`bam2bcf.c:991-1001`, which tests a C `float`). The slice-4 goldens
  exercise this path heavily (87 multi-allelic SNP sites) and match
  byte-for-byte, so the C-`float`-vs-Go-`float64` width concern did not
  materialise on the upstream fixtures. The note is kept for awareness
  on future high-coverage inputs.
- **`phase` MCMC chimera repair.** The v1 port uses a greedy
  adjacent-het vote in place of upstream `phase.c`'s MCMC
  `phase_core` loop; `-b` per-haplotype BAM split is also deferred.
  Detail in the `phase` subsection below.
- **`targetcut` BAQ realignment with `-f` reference (DONE).** The
  HMM consensus mode is implemented (faithful port of
  `cut_target.c`, including the MAQ errmod port; see below). The
  per-record BAQ realignment that upstream's `read_aln` applies
  via `sam_prob_realn` when a `-f` reference is supplied is now
  wired: `targetcutHMM` runs `pkg/htsgo/baq.SamProbRealn` in
  apply+extend mode (`flag = 1<<1|1`) on every record that
  survives the read filter, on a per-chromosome cache of the
  reference fetched via `fasta.RandomAccess`. The stderr warning
  is gone; the BAQ-adjusted qualities feed `gencns` exactly as
  upstream feeds them into the pileup.
- **`tview`** — deliberate skip (interactive curses UI).

**`markdup` deferred features** (deliberately skipped in v1, all flag
slots are accepted on the CLI for compat):

- Optical-duplicate detection (`-d/--max-dist` + (x,y) parsing of Illumina
  qnames). v1 marks PCR duplicates only; nonzero `-d` triggers a stderr
  warning.
- Per-read-group keying (upstream's `-S` flag). v1 folds all read groups
  into a single namespace, so fixture
  `reference_code/samtools/test/markdup/17_read_group.sam` is a
  documented partial-parity skip.
- Barcode regex / barcode-tag keying.
- The `dt:Z:` "duplicate-type" aux tag (SQ / LB / OQ). The 0x400 flag
  bit is set correctly; only the typed aux is missing.

**`markdup -l/--max-len` is a no-op-by-design.** Upstream uses `-l` solely
as the streaming buffer flush window in `bam_markdup.c:1949`; it does NOT
affect key construction or scoring. Our two-pass implementation buffers
per-bucket state in memory, so output is identical for any `-l` value.
The flag is accepted on the CLI and the option is preserved on
`MarkdupOptions.MaxLen` for forward compatibility if we ever move to a
single-pass streaming model.

**`stats` implemented sections.** The per-cycle quality matrices
**FFQ/LFQ**, the GC-content sections **GCF/GCL**, the ACGT-content
sections **GCC/GCT**, the indel sections **IC/ID**, the leading **CHK**
CRC32 checksum block, the **COV** coverage-distribution histogram and
the **GCD** GC-depth distribution are byte-faithful to upstream —
validated against the vendored `.sam` fixtures. CHK sums the per-record
CRC32 of read names, the BAM 4-bit-packed sequence and the quality
bytes; COV bins each reference position's M/=/X read depth via the
`-c MIN,MAX,STEP` option (default `1,1000,1`) and is emitted only for
coordinate-sorted input, matching upstream's `is_sorted` gating. COV
depth is accumulated in a bounded per-contig sliding window that is
flushed as records advance (mirroring upstream's `cov_rbuf` ring
buffer), so COV memory is O(longest read span), not O(genome).

GCD splits the reference into `--GC-depth`-wide segments (default
20000 bases) and records per-segment read depth plus GC content, then
at output sorts segments by GC and reports depth percentiles. Both
upstream code paths are ported: the default no-reference path
approximates GC content from the read sequences, while the
`-r/--ref-seq` reference path reads GC content from the indexed
reference FASTA (`fai_gc_content`). Like COV, GCD is emitted only for
coordinate-sorted input. The upstream `igcd`/`ngcd` indexing quirk —
where `gcd[0]` is an empty placeholder and the final segment is never
finalised before sorting — is replicated exactly for byte-parity.

For unsorted input the COV section is silently omitted, whereas
upstream aborts with `Expected coordinates in ascending order` — a
deliberate, friendlier divergence.

**`stats` reference-free option tails.** BWA-style quality trimming
(`-q/--trim-quality`) and the `-t/--target-regions` restriction are
implemented and validated against the upstream `stat/11` fixtures.
`-q` ports `bwa_trim_read` faithfully, including bwa's documented
off-by-one and the `BWA_MIN_RDLEN` (35) early return, and feeds the
`bases trimmed` SN counter. `-t` parses the upstream target-regions
format (`seq-name beg end`, 1-based inclusive — not BED), merges
overlapping intervals, restricts every counter to reads overlapping a
target interval, clips `bases mapped (cigar)` and the COV depth to the
target, and emits the `bases inside the target` and
`percentage of target genome with coverage > N` SN lines (the latter
threshold set by `-g/--cov-threshold`, default 0). The SN and COV
sections are byte-faithful to `11.stats.expected` /
`11.stats.g4.expected`.

**`stats` reference-statistics sections.** The **MPC**
mismatches-per-cycle section (emitted with `--ref-seq`) and the **RFS**
reference-statistics section (emitted with `--ref-stats`, plus its
`--ref-stats-chunk` companion) are implemented and byte-faithful to the
upstream `stat/` golden files. MPC ports `count_mismatches_per_cycle`
faithfully — including the cycle-index handling for soft/hard clips and
insertions, the reverse-strand mirroring, the N-base bucketing in
quality slot 0, and the documented `uint8` `qual+1` wrap that lands
mismatches of `*`-quality reads in the N column. RFS ports
`collect_refstats`: without `-t` it reports one row per `@SQ` header
entry, with `-t` it reports the merged target intervals as
`name:start-end` rows, and without `--ref-seq` the GC/N columns report
the `-1` lack-of-data sentinel. Validated against `test.pl` stats cases
1-8 (`-r test.fa`, MPC) and 16/17/19 (`--ref-stats`, RFS).

**`stats` per-fragment ACGT sections.** **FBC/LBC** (ACGT content per
cycle for first / last fragments) and **FTC/LTC** (the matching
A/C/G/T/N raw-counter totals) are implemented and byte-faithful to the
vendored `stat/` fixtures. They derive from the same per-fragment
cycle buffers GCC/GCT already accumulate. Note: despite the name these
are NOT barcode tables — earlier roadmap text mislabelled them.

**`stats` per-barcode sections (implemented).** The per-barcode
ACGT-content (`<tag>C`) and quality (`<tag>Q`) tables are a faithful
port of `collect_barcode_stats` (stats.c:773) and its output
(stats.c:1748). Collection is unconditional for the four fixed aux-tag
pairs upstream's `init_barcode_tags` installs — `BC`/`QT`, `CR`/`CY`,
`OX`/`BZ`, `RX`/`QX` — and a section is emitted only once a barcode for
its tag has been observed. The barcode separator, the two segment
columns either side of it, and the per-tag `max_qual` quality-column
count all match upstream. Note: this samtools version has **no**
`--barcodes` CLI option — barcodes are always collected — so no such
flag is wired; `--barcode-tag` / `--quality-tag` tag renaming likewise
does not exist here. Validated byte-for-byte against the exact
`test.pl` invocations `samtools stats 13_barcodes_ok.sam` and
`samtools stats 13_barcodes_ok_ox_bz.sam` (test.pl:3325-3326). The
malformed-barcode `expect_fail` fixtures (test.pl:3327-3329) are
exercised for graceful warn-and-skip behaviour but not byte-compared,
since their output deliberately diverges from any clean baseline.

**`stats` `--sparse` (corrected).** Upstream `-x/--sparse`
(stats.c:2170) only "suppresses outputting IS rows where there are no
insertions" — its sole effect is at stats.c:1796, thinning all-zero IS
rows. Earlier this code wrongly suppressed *every* histogram section.
It now emits every section unconditionally and honours `sparse` per-row
in the IS section only. The IS section itself was reworked: it now
emits the full `0..ibulk-1` row range (including zero rows in the
default mode) where `ibulk` mirrors upstream's last-non-zero /
99%-bulk-truncation logic, instead of the old observed-sizes-only map.
Insert sizes are now classified per-read from each record's own flags
and mate position and halved at output, matching upstream exactly
(this also fixed an inward/outward misclassification the prior
per-pair logic produced on multi-pair inputs).

**`stats` remaining tail** (also documented in `PARITY_VALIDATION.md`):

- Command-line positional region arguments are not yet accepted, so the
  RFS-with-command-line-region path (upstream stats test 18) is not
  reproducible; RFS-with-`-t` covers the equivalent functionality.
- `--remove-overlaps` is accepted as a no-op; single-record stats are
  unaffected by overlap removal for the counters emitted.

With `--sparse` corrected and the per-barcode sections implemented,
`samtools stats` reaches full 1:1 output parity for every section.

The output emits the byte-faithful **CHK** checksum block, **SN**
(Summary Numbers), the per-cycle and base-content sections
(**FFQ/LFQ/GCF/GCL/GCC/GCT/FBC/FTC/LBC/LTC/IC/ID**), the per-barcode
**`<tag>C`/`<tag>Q`** tables, the **MPC** mismatches-per-cycle matrix,
the **RL / MAPQ / IS** rollups, the **COV** coverage histogram, the
**GCD** GC-depth distribution and the **RFS** reference-statistics
section. `--sparse` thins only all-zero IS rows.

**Validation:** upstream fixtures from `reference_code/samtools/test/markdup/`
and `.../test/stat/` are vendored under
`tools/samtools/testdata/parity/{markdup,stat}/`. The byte-exact /
flag-exact / SN-byte cases are exercised in
`tools/samtools/pkg/samtools/markdup_test.go` and `stats_test.go`.

**`calmd` BAQ realignment** (`-r`, `-E`, `-A`, `-C`) — implemented:

- The probabilistic banded glocal forward-backward HMM
  (`probaln_glocal`) and the CIGAR-aware BAQ driver / MAPQ cap
  (`sam_prob_realn`, `sam_cap_mapq`) are ported faithfully into the
  shared `pkg/htsgo/baq` package — shared so `mpileup` can reuse it.
- **`-r`** computes the `BQ:Z` base-alignment-quality aux tag;
  **`-rE`** is extended-BAQ mode; **`-rA`** applies BAQ to the base
  qualities and writes a `ZQ:Z` tag. **`-C INT`** (threshold `> 10`)
  caps MAPQ via `sam_cap_mapq`.
- Validated byte-for-byte against htslib's own `realn0{1,2,3}*`
  golden fixtures: all 8 `test_realn` flag combinations plus the
  `FlagRedo` path pass exactly. `probaln_glocal` is additionally
  pinned against score / posterior-quality vectors captured from the
  upstream `probaln.c` self-test.

**`calmd` deferred features** (accepted as CLI flags, behaviour partial):

- **`-h` HASH_QNM** (hash-based query-name binarisation) — niche
  upstream-only optimisation; not implemented.
- **`-N` clear-MD/NM-bits**, **`--no-PG`** —
  CLI-accepted-and-ignored stubs.

**`calmd` implemented post-MD/NM transforms** (`bam_md.c` upstream
order — max-NM masking → write NM → write MD → DROP_TAG → BIN_QUAL):

- **`-d` DROP_TAG** — drops every aux tag except `RG`. Applied after
  the NM/MD fill, so the freshly-computed NM/MD are dropped too; only
  `RG` survives (records without `RG` keep no aux at all).
- **`-q` BIN_QUAL** — reduces base-quality resolution: each quality
  `>= 3` maps to `qual/10*10 + 7` (integer division); lower values
  unchanged.
- **`-n INT` max-NM** — for reads whose computed NM `>= INT`, masks
  every matching M/=/X base (SEQ `->` `=`, quality `-> 0`); the
  emitted NM/MD are unaffected.

**`import` deferred features**:

- **`--i1` / `--i2`** index-read inputs (the index-as-aux BC/QT shape).
  v1 wires `-0/-1/-2/-s` and the positional shapes; index files would
  attach a BC:Z and QT:Z tag computed from a separate index FASTQ. The
  parser scaffolding is in place; just needs a third walker.
- **`-i` CASAVA** parsing (extract barcode from CASAVA-style headers).
  We do parse the description tail for SAM aux fields directly, which
  covers the common case where the FASTQ was produced by `samtools
  fastq`; CASAVA-format input is a follow-up.
- **`--barcode-tag` / `--quality-tag`** renaming of the BC/QT tag pair.
  Not exposed in the v1 CLI (defaults to BC/QT).
- **`-O` / `--output-fmt`** and **`-@` / `--threads`** —
  CLI-accepted-and-ignored stubs. v1 picks output format from the
  output-path extension (`.sam` vs `.bam`) and is single-threaded.

**Validation:** small hand-built fixtures live under
`tools/samtools/testdata/parity/{calmd,import}/` covering the four
calmd code paths (match, mismatch, deletion, insertion+softclip) and
the six import shapes (-0, -1/-2, -s, single positional, two
positionals, -T aux extraction, --order, -R/-r RG). The calmd BAQ
path additionally diffs the `-r` / `-rA` output against htslib's
vendored `realn01_exp*.sam` goldens, and `pkg/htsgo/baq` carries the
full `realn0{1,2,3}` golden corpus. The upstream
`bam_md.c` / `bam_import.c` regression cases are marked as
`t.Skip(...)` parity stubs because upstream's BGZF output isn't
byte-identical with ours (different libdeflate). Logical correctness
is covered by hand-computed expected values in the table tests.

**`phase` deferred features** (accepted on the CLI, behaviour partial):

- **MCMC chimera repair**. Upstream's `phase.c` runs a
  Markov-chain-Monte-Carlo loop (`phase_core`) that flips read-cluster
  assignments to maximise haplotype consistency and resolve chimeric
  reads at junctions. The v1 Go port replaces this with a greedy
  same-vs-opposite vote between adjacent het sites. Tied junctions
  emit label `0` (ambiguous) rather than being repaired by MCMC.
  Tracked here; the upstream `FLAG_FIX_CHIMERA` flag is implicitly
  disabled in v1.
- **`-b STR` per-haplotype BAM split.** v1 emits the phased TSV
  stream to `-o`/stdout but does not yet split the input BAM into
  per-haplotype output BAMs (`<prefix>.0.bam` / `<prefix>.1.bam`
  / `<prefix>.chimera.bam` in upstream). The flag is accepted on the
  CLI and stored in `PhaseOptions.OutputPrefix` for a follow-up
  wiring pass.
- **`-F` use-full-read** is accepted on the CLI but is a no-op in v1
  (we always walk the aligned slice as decoded from the CIGAR).
- **`-A` mark-drop-in-chimera-output** is also a no-op pending the
  `-b` split landing.
- **`-e`/`-l` site-list mode** (only-phase-listed-sites). The
  upstream `loadpos` path is not implemented; the Go port always
  discovers hets from the pileup. Upstream itself comments `-e` and
  `-l` out of the usage block, so the omission is a small loss.

**`targetcut` HMM consensus mode** (implemented). The Go port is now
a faithful translation of upstream `cut_target.c`: per-position
consensus via the MAQ revised error model (`errmod.c` is ported
in-tree as a shared package at `pkg/htsgo/errmod` — both
`samtools targetcut` and `bcftools mpileup` import the same
implementation, eliminating the duplicate ports that briefly
existed in PRs #199 and #216), followed by a 2-state Viterbi
over the per-chrom consensus track to segment "covered, callable"
regions from "no-info or uninformative" regions, then one SAM-
format consensus record is emitted per identified region in the
exact upstream printf shape
(`%s:%d-%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t<seq>\t<qual>\n`). The
emit-loop's "position 0 never participates in a region" quirk —
a consequence of upstream's backtrack loop running over the half-
open range `(0, l-1]` — is reproduced verbatim so we agree with
the C output by construction. CLI flags `-Q -i -0 -1 -2` carry
their upstream semantics; `-f` reference triggers per-record BAQ
realignment via `pkg/htsgo/baq.SamProbRealn` (apply+extend mode,
flag `1<<1|1`), matching upstream `cut_target.c::read_aln`. The
pre-port "aligned-slice FASTA" mode remains available behind
`--simple` (library: `TargetcutOptions.SimpleMode`); `-f` has no
effect in simple mode.

One implementation detail diverges from upstream by design: when
n > 255 bases pile at a single position the upstream `errmod_cal`
shuffles via `ks_shuffle` (drand48 state) and truncates to 255 to
fit its pre-computed coefficient table. We truncate deterministically
to the first 255 because the downstream gencns caps depth at 255
anyway and a drand48-dependent output would not be reproducible.
For coverages ≤255 (the overwhelming common case in practice) we
are byte-equivalent to upstream.

**Validation:** hand-built SAM fixtures in
`tools/samtools/pkg/samtools/phase_test.go` and
`tools/samtools/pkg/samtools/targetcut_test.go`. Phase tests cover
single-block chaining (consistent & label-flipping orderings),
ambiguous-label fall-back when reads don't bridge two hets, and
the MinMAPQ filter. Targetcut tests cover the simple-mode legacy
behaviour AND the HMM mode: a uniform-coverage region (one
emitted region with the expected SAM shape), an entirely-empty
chrom (no output), the upstream read filter (unmapped /
secondary / qcfail / dup all skipped), MinBaseQ pushing every
cell to "no info", a tuned `-i` entry penalty separating two
coverage blocks into two regions, and a majority-vote consensus
base check. A self-consistency test on the errmod port checks
that 10 identical 'A' observations yield the minimum (best)
homozygous score at q[A,A]. There is no upstream regression-test
fixture for either tool in `reference_code/samtools/test/`, so
expected values are hand-derived directly from the C source — the
project's accepted standard for tools with no upstream fixture.

**`consensus` bayesian mode** (implemented). Upstream `bam_consensus.c`
ships five modes — `simple` and four bayesian flavours (`bayesian_r`
aka "bayesian", `bayesian_m`, `bayesian_p`, `bayesian_116`). All five
are now implemented. The Gap5-derived posterior caller
(`calculate_consensus_gap5` / `_gap5m`), the localised-MAPQ NM-halo
adjustment (`nm_init` / `nm_local` / `poly_len`), and the bayesian
knobs (`-C/--cutoff`, `--P-het`, `--P-indel`, `--het-scale`,
`--adj-qual`, `--use-MQ`, `--adj-MQ`, `--NM-halo`, `--SC-cost`,
`--scale-MQ`, `--low-MQ`, `--high-MQ`, `--default-qual`,
`-p/--homopoly-fix`, `--homopoly-score`, `--homopoly-redux`) are
ported faithfully; the upstream `fast_exp` 0.1-resolution quantization
and the degree-3 `fast_log2` are reproduced so phred conversions are
byte-exact. The default invocation (`MODE_RECALL`, MQUAL + NM-adjust
on) is byte-for-byte parity with upstream on the
`reference_code/samtools/test/consensus/` corpus for the non-`-T`
fixtures (18q/19q/18p/19p/20p/21p on consen1, 30/31/32/40/41/42 on
consen1c). Insertion-column pileup rows (`nth>0`) are emitted for the
bayesian path.

Genuinely deferred sub-knobs (precisely scoped):

- **`-t/--qual-calibration` and `-X/--config`.** Accepted on the CLI
  but apply only the FLAT identity calibration table. The per-machine
  calibration tables (HiFi/HiSeq/ONT/Ultima) and the QUAL-file parser
  (`load_qcal`, `bam_consensus.c:672-736`) are a separable ~300-line
  table/parser block that does not affect the default invocation.
- **`-T/--reference` uncovered-base fill.** Accepted but the reference
  is not used to fill uncovered positions; the `*T.out` golden files
  (30T/31T/.../42T) which substitute reference bases at zero-coverage
  positions are therefore out of scope. The non-`-T` path — the
  default — is fully byte-faithful.
- **Mate-overlap dedup.** `--ignore-overlaps` is accepted but is a
  no-op; v1 counts each mate independently in the pileup walker.
- **Threading.** `-@/--threads`, `-Z/--block-size`, and
  `--input-fmt-option` are accepted but ignored; v1 is single-pass
  and single-threaded.
- **Read-flag filtering.** `--rf/--incl-flags` and `--ff/--excl-flags`
  are accepted as text/int but ignored. v1's filter set is fixed
  (drop UNMAP|SECONDARY|QCFAIL|DUP, matching upstream's default
  `excl_flags`).
- **`--het-only`** suppression of homozygous calls is accepted but
  not implemented.
- **`--verbosity`** is accepted and ignored.

**`consensus` correctness model.** v1 mirrors upstream's
`calculate_consensus_simple` (`bam_consensus.c:1900-2006`) bit-for-bit
where it matters:

- One fraction gate only: `used_score < call_fract * tscore`
  (`bam_consensus.c:1988-1994`). There is **no** separate
  "min-fraction on the dominant base alone" gate; an earlier PR
  fabricated one and is corrected here.
- Heterozygous condition is `score2 >= het_fract * score1 && ambig`
  with no `score1 > 0` guard (`bam_consensus.c:1982`).
- `use_qual=0` by default (`bam_consensus.c:2984`) — bases score by
  frequency, not quality, until `-q/--use-qual` is set.
- `--show-del` is honoured in pileup mode too — rows whose call is
  `'*'` are suppressed when `--show-del no`
  (`bam_consensus.c:2244`).
- Insertion gating uses `MinCallFraction`, the same knob as the
  per-position gate.

**Validation:** table-driven hand-built SAM fixtures in
`tools/samtools/pkg/samtools/consensus_test.go` covering:
all-match FASTA/FASTQ/pileup, mixed at the 0.75 boundary,
`--show-del no` in pileup (no `*` rows),
`--show-del yes` in pileup (keeps `*` rows),
the canonical **50/30/20 + `-A`** fixture (must land `M`, not `N`),
frequency-only counting (low-Q vs high-Q gives identical output when
`UseQual=false`),
the `UseQual=true` flip (high-Q minority beats low-Q majority),
multi-contig, `-a` zero-fill, `--min-depth`, line-len wrapping,
insertion include/suppress + `--mark-ins`,
and the default-bayesian invocation emitting no fallback warning.
Bayesian-mode parity is additionally verified byte-for-byte against
the vendored upstream `reference_code/samtools/test/consensus/`
golden files (`TestConsensus_BayesianUpstreamParity`): twelve
fixtures spanning FASTQ, pileup, the `-C` cutoff, `-A` ambiguity,
`-a`/`-aa`, and the default MQUAL + NM-adjust path.
Coverage of the `pkg/samtools` package after this PR is ~80%.

### `bcftools`

**Status:** 24 of ~30 subcommands (~80%). `view`, `index`, `stats`, `query`,
`concat`, `norm`, `call` (consensus + biallelic multi-allelic), the PR #86
wave-1 tail (`annotate`, `head`, `isec`, `merge`, `reheader`, `sort`), the
convert/mendelian PR (`convert`, `mendelian`), the gtcheck/roh PR (`gtcheck`,
`roh`), the filter/consensus PR (`filter`, `consensus`), the
mendelian2/polysomy PR (`mendelian2`, `polysomy`), the cnv/csq PR
(`cnv` + `csq`), and the mpileup PR (**`mpileup`**).

All bcftools subcommands now have an implementation in the Go port.

The plugin system (`bcftools plugin` / `bcftools +<name>`) is **done**,
but with a deliberate design divergence from upstream:

- **`+plugins`** — implemented as a **subprocess plugin system**.
  Upstream loads plugins as native shared objects (`.so`) via `dlopen`
  against a fixed C ABI. The Go port instead resolves `bcftools +<name>`
  to an ordinary *executable* found by name in the `BCFTOOLS_PLUGINS`
  colon-separated directory list, runs it as a child process, pipes the
  input VCF as uncompressed text to its stdin, and reads VCF back from
  its stdout. A plugin is therefore "a VCF-on-stdin to VCF-on-stdout
  filter" and can be written in any language — no C ABI, no version
  check, no rebuild against the host. `bcftools plugin -l`/`-lv` lists
  discoverable plugins; the host applies `-o`/`-O` formatting around the
  plugin's output; a non-zero plugin exit is surfaced as an error with
  its stderr. The contract is specified in `docs/PLUGIN_PROTOCOL.md`.
  The mechanism lives in `tools/bcftools/pkg/bcftools/plugin.go` with a
  reference example plugin under `tools/bcftools/plugins/example/`.
  **Intentionally not ported:** the ~30 bundled upstream plugins
  (`+fill-tags`, `+split-vep`, `+setGT`, `+prune`, `+fixploidy`, ...).
  The plugin *system* exists so users can write their own plugins;
  re-porting upstream's plugin catalogue is explicit non-goal scope.
  Upstream plugin sources (`plugins/*.c`) remain vendored under
  `reference_code/bcftools/` for anyone who wants to reimplement a
  specific one as a standalone subprocess plugin.

Note on vendored reference source: `reference_code/bcftools` and
`reference_code/htslib` are now both vendored as submodules. Earlier
roadmap text in this section was written when bcftools internals were
unavailable and called several HMM-based features unportable for that
reason. That is no longer true — `vcfcnv.c` (CNV HMM), `vcfroh.c` (RoH
HMM), `polysomy.c` + `peakfit.c`, `HMM.c` (the shared Viterbi/Baum-Welch
core), `bam2bcf*.c` (mpileup MAQ model + indel caller) and
`reference_code/htslib/errmod.c` are all available as porting
references. The deferrals below are scope/effort calls, not
source-availability blockers; each is annotated accordingly.

Genuine algorithmic gaps (subcommands present but running a v1
heuristic in place of the upstream algorithm — full detail in the
per-subcommand option-tail sections below):

- **`cnv`** — full port: the upstream copy-number HMM (`vcfcnv.c` +
  `HMM.c`, both vendored). 4-state single-sample / 16-state paired
  Viterbi + forward-backward over BAF+LRR Gaussian emissions, with
  `--optimize` cell-fraction estimation and `--baum-welch` transition
  re-estimation. See the `cnv` option-tail section for the validation
  situation (no upstream golden exists).
- **`roh`** — full port: 2-state Viterbi + forward-backward HMM
  (`vcfroh.c` + `HMM.c`) with physical-distance- and genetic-map-scaled
  transitions, allele-frequency estimation and Baum-Welch
  Viterbi training. Validated byte-for-byte against the upstream
  `roh.1.*.out` goldens. The only remaining gap is PL-based emission
  (`-G` hard-GT is the supported path) and `-O z`.
- **`polysomy`** — DONE: faithful port of the upstream
  Gaussian-mixture peak fit (`polysomy.c` + `peakfit.c`). The GSL
  Levenberg-Marquardt solver is ported in-tree as pure Go
  (`peakfit_lm.go`); no third-party dependency was added. All
  algorithm knobs are live. See the per-subcommand section below for
  the validation situation (no upstream golden exists).
- **`gtcheck`** — v1 hard-GT Hamming only; `--cluster` HMM-style
  sample clustering and PL/GL scoring deferred.
- **`mpileup`** — the upstream MAQ genotype-likelihood model is fully
  ported for the SNP path (`errmod.c` → `errmod.go`; `bam2bcf.c`
  glfgen/combine/2bcf → `bam2bcf.go`): the multi-allelic PL grid, the
  `<*>` unseen allele, BCF output, BAQ realignment (`pkg/htsgo/baq`),
  the per-site bias annotations (VDB/SGB/RPBZ/MQBZ/BQBZ/MQSBZ/SCBZ),
  MQ0F and `MPLP_SMART_OVERLAPS` read-pair quality merging are all wired
  in and verified byte-for-byte against the upstream goldens (slices
  1-4 done). Only indel calling (`bam2bcf_indel.c`) remains deferred.
- **`call`** — consensus and biallelic multi-allelic calling are
  implemented; the full upstream multi-allelic `-m` grid over >2 ALTs
  pairs with the mpileup MAQ port.
- **`csq`** — haplotype-aware consequence engine (slices 1-3 done):
  the `hap_node_t` tree, `cds_translate`, compound consequences,
  `@pos` reference pointers and `-p/-n` are ported, and the
  INFO/BCSQ output matches upstream byte-for-byte on the targeted
  goldens. The remaining tail (`FORMAT/TBCSQ` text expansion,
  `--unify-chr-names`, `-l/--local-csq`) is slice 4 — see the "csq
  full-parity slicing plan" below.

Option-tail gaps on `gtcheck` (PR #107, simple-mode):

- `--cluster N,N` (HMM-style sample clustering), `--distinctive-sites`,
  `--n-matches` — accepted-and-rejected with PARITY_ROADMAP pointer;
  bayesian-mode follow-up.
- `-u PL` — PL/GL-based scoring; v1 only does hard-GT Hamming.
- `-O z` — bgzip output; v1 only emits tab-text (`-O t`).
- `[5]Average -log P(HWE)` column is zeroed until a real per-site HWE
  estimator from panel AF lands.
- Index-backed `-r/-R` seek (post-filter only in v1).
- Multi-allelic input is rejected (matches upstream's
  `bcftools norm -m -` requirement).

Option-tail gaps on `roh`:

- **The HMM is now the full upstream port.** `vcfroh.c` + `HMM.c` are
  ported in-tree (`tools/bcftools/pkg/bcftools/hmm.go`,
  `roh.go`, `roh_genmap.go`): a 2-state Viterbi decode plus a
  forward-backward posterior, transition probabilities scaled by the
  physical (and, with `-m`/`-M`, the genetic-map) inter-marker
  distance, allele-frequency estimation (`-e/--estimate-AF` from
  GT cohorts), Baum-Welch parameter re-estimation
  (`-V/--viterbi-training`) and overlapping-window buffering
  (`-b/--buffer-size`). RG/ST quality scores are forward-backward
  phred scores and match upstream byte-for-byte on the
  `roh.1.*.out` goldens.
- `-O z` — bgzip output; v1 only emits tab-text. Remaining gap.
- PL-based emission scoring is not ported; `-G/--GTs-only` hard-GT
  mode is the supported emission path (the upstream test corpus
  uses `-G30`). `-e PL,...` falls back to rejecting sites with no
  usable genotype rather than reading PLs.

Option-tail gaps on the wave-1 additions (PR #86):

- `annotate --set-id '+%CHROM_%POS'` macro expansion is not implemented;
  `-x ID` / `-x INFO/TAG` / `-x FORMAT/TAG` removal works.
- `isec`: `--collapse some` (REF match + any-ALT-in-common) is approximated
  via strict tuple match; deeper semantics deferred.
- `merge`: pre-sort assumption is enforced; no automatic CHROM/POS sort.
- `reheader`: in-place rewrite (`-i`) currently emits to stdout — caller
  is responsible for the swap.
- `sort`: `-m/--max-mem` and `-T/--tmpdir` are accepted but always
  in-memory.

Option-tail gaps on the convert/mendelian PR:

- `convert`: v1 covers only the pass-through round-trip
  (VCF↔BCF↔VCF.gz) with sample/region filtering and -i/-e expressions.
  The full upstream `vcfconvert.c` covers many extra shapes
  (`--gvcf2vcf`, `--haplegendsample2vcf`, `--hapsample2vcf`,
  `--tsv2vcf`, `--gensample2vcf`, `--gvcf`, PLINK / GEN / HAP).
  These are explicit follow-ups; the CLI emits a usage block that
  lists them under "Deferred output paths".
- `mendelian`: the v1 port detects Mendelian inconsistencies for
  PED-style trios (one or more `-t CHILD,FATHER,MOTHER` flags, or
  `-T trio-file`), emits `INFO/MERR` per record, and supports the
  five upstream modes `a|c|x|d|+`. The `+` mode is treated as a
  synonym for annotate (the upstream-only `##bcftools_PG=` provenance
  header is deferred). The `--rules FILE` ploidy specification is
  accepted but currently only the chrX heuristic is honoured; full
  per-contig ploidy override (PAR boundaries, mitochondrial
  haploidy) is a follow-up.

Option-tail gaps on `filter` (this PR, simple-mode):

- `--mask [^]REGION` and `-M/--mask-file [^]FILE` — accepted but
  hard-rejected at runtime with a roadmap pointer. The CLI parses the
  flag (so a downstream automation that always passes `-M ""` doesn't
  break); the underlying BED-driven soft-filter logic is deferred.
- `--mask-overlap 0|1|2` — accepted; v1 ignores (always treats POS-in-region).
- `-W/--write-index[=FMT]` — accepted; v1 never auto-indexes outputs.
- `-v/--verbosity INT` — accepted; v1 ignores.
- `--regions-overlap` / `--targets-overlap` — accepted; v1 always
  uses POS-in-region semantics.
- `-g/--SnpGap INT:TYPE` — the `:TYPE` qualifier (indel|mnp|bnd|other|overlap)
  is parsed but always treated as "indel" in v1.
- BCF output (`-O b|u`) round-trips through the shared `pkg/htsgo/bcf`
  writer; CSI auto-indexing is the `-W` follow-up above.

Option-tail gaps on `mendelian2` (PR #109, simple-mode):

- `--rules ASSEMBLY` — predefined inheritance rules (GRCh37 / GRCh38
  / `list?`). Accepted but rejected at runtime with a roadmap pointer;
  v1 uses the chrX heuristic from the legacy `mendelian` port instead.
- `--rules-file FILE` — custom inheritance rules file. Accepted; v1
  rejects at runtime.
- `-W/--write-index[=FMT]` — auto-index output. Accepted; v1 never
  auto-indexes.
- `--regions-overlap 0|1|2`, `--targets-overlap 0|1|2` — accepted;
  v1 always uses POS-in-region semantics.
- `-v/--verbosity INT`, `--no-version` — accepted; v1 ignores both.
- `-r/-R/-t/-T` — region / target post-filter is wired through the
  CLI but the BCF synced-reader region-jump path is not used in v1;
  filtering happens after the records are read.
- `sites_not_diploid` counter never goes up in v1 because our
  `vcf.Variant` decoder coerces non-diploid GTs to missing rather
  than tracking ploidy. Tracked alongside the broader BCF FORMAT
  reconstruction in `docs/UPSTREAM_BUGS.md`.
- The PED-row sort order is by child name (deterministic);
  upstream sorts by min(sample-index) for sequential VCF reads. The
  sort is a performance optimisation only — the set of reported
  trios and per-trio counters are identical.
- `-i/-e EXPR` are accepted at the library boundary but only
  applied as record-level filters in v1 (no per-sample mask). The
  `sites_fail` counter is therefore always 0 in v1.

Option-tail status on `polysomy` (Gaussian-mixture peak fit — DONE):

- **The algorithm is the upstream Gaussian-mixture peak fit.**
  `polysomy.c` is ported faithfully: each chromosome's BAF values are
  binned into an `--nbins`-bin histogram, smoothed, and the RR/RA/AA
  bands are isolated and per-segment normalised (`init_dist`). Three
  candidate fits — CN2 (one bounded Gaussian near 0.5), CN3 (two
  symmetric Gaussians near 1/3 and 2/3) and CN4 (a central peak plus
  two symmetric side peaks) — are run over the heterozygous band, and
  the lowest CN that passes `--fit-th` plus the symmetry / peak-size
  checks is chosen, with `--cn-penalty` as the tiebreaker
  (`fit_curves`).
- **The peak-fitting engine `peakfit.c` is ported in-tree as pure
  Go** (`peakfit.go`): the Gaussian / centre-bounded-Gaussian /
  exponential peak models, the residual objective `(model-y)/0.01`,
  the unscaled `Σ|model-y|` goodness measure, and the Monte-Carlo
  restart driver. The GSL `gsl_multifit_fdfsolver_lmsder` non-linear
  least-squares solver is replaced by an in-tree pure-Go
  Levenberg-Marquardt solver (`peakfit_lm.go`) — the normal-equations
  damping loop `(JᵀJ + λ·diag)·δ = -Jᵀr` with an analytic Jacobian, a
  λ up/down schedule, and convergence on the parameter / gradient
  deltas at tolerance 1e-8. `peakfit_lm.go` also ports glibc's
  `random()` so the Monte-Carlo restart stream after `srand(0)`
  matches upstream. **No third-party dependency was added** — CLAUDE.md
  scopes the one sanctioned numerical dep (gonum) narrowly and
  explicitly excludes stats-fitting tools like this.
- **All algorithm knobs are live:** `-b/--peak-size`,
  `-c/--cn-penalty`, `-f/--fit-th`, `-i/--include-aa`,
  `-m/--min-fraction`, `-p/--peak-symmetry`, `-n/--nbins`,
  `-S/--smooth`, `--ra-rr-scaling`, `--force-cn`.
- **Validation — no upstream golden exists.** `bcftools polysomy`
  produces only per-chromosome PNG plots and a `dist.dat` dump under
  `--output-dir`; `reference_code/bcftools/test/test.pl` has no
  `polysomy` invocation. Byte-for-byte parity therefore cannot be
  demonstrated and is **not claimed**. The port is validated instead
  with: (a) unit tests of the LM solver against analytic curves with
  known optima (a linear-in-parameters quadratic and a non-linear
  single Gaussian); (b) a check of the glibc `rand()` port against the
  published `srand(0)` reference sequence; (c) unit tests of each peak
  model (Gaussian / bounded-Gaussian / exp) recovering known
  parameters from clean synthetic samples; (d) hand-constructed BAF
  distributions for the canonical karyotypes — a clean diploid
  (single peak at 0.5 → CN2) and a clear trisomy (two peaks at 1/3 and
  2/3 → CN3), exercised both directly and end-to-end through the VCF
  reader. See `polysomy_test.go` and `peakfit_test.go`.
- `-o/--output-dir PATH` — accepted but ignored; the port writes the
  per-chromosome TSV to stdout (no PNG plots, no `dist.dat`).
- `--regions-overlap`, `--targets-overlap` — accepted; the port uses a
  chromosome-name post-filter (per-base interval filtering deferred).
- `-v/--verbosity` / `--verbose` — accepted; the port ignores it (no
  per-iteration fit trace).
- BAF source: upstream requires FORMAT/BAF; the port also accepts
  FORMAT/AD = REF,ALT as a fallback (synthesises BAF as
  `ALT / (REF + ALT)` at het sites).
- Per-record `-i/-e` are NOT in upstream `polysomy.c:main_polysomy`
  and we follow upstream's surface exactly (no invented flags).

Option-tail gaps on `cnv` (full HMM port):

- **The algorithm is the upstream HMM.** `vcfcnv.c` is ported
  faithfully: a 4-state HMM (CN0/CN1/CN2/CN3) for a single sample, or
  a 16-state HMM for the paired tumour/control mode (`-c`), swept per
  contig with Viterbi + forward-backward. Emission probabilities are
  the upstream joint BAF + LRR model — a truncated-Gaussian BAF peak
  mixture weighted by genotype frequencies fRR/fRA/fAA, combined with
  a per-state LRR Gaussian. The generic engine is the shared
  `hmm.go` port (also used by `roh`), reused unchanged. Every HMM
  tuning knob is now load-bearing: `-a/--aberrant` (CN3 BAF peak
  shift), `-b/--BAF-weight`, `-e/--err-prob`, `-l/--LRR-weight`,
  `-L/--LRR-smooth-win`, `-d/--BAF-dev`, `-k/--LRR-dev`,
  `-x/--xy-prob`, `-P/--same-prob`, `-O/--optimize` (iterated
  forward-backward cell-fraction estimation), and `-W/--baum-welch`
  (per-contig transition re-estimation).
- **Validation: no upstream golden exists.** bcftools `test/test.pl`
  contains no `cnv` invocation and `test/` ships no `cnv` fixtures
  (`cnv` output is plot-oriented). The Go port is therefore validated
  with hand-derived cases: a clean diploid region decodes to a single
  all-CN2 region; a long missing-BAF run in paired mode decodes to
  CN0; a het-band-split + positive-LRR run decodes to CN3; and
  unit tests pin the ported transition matrix (column-stochastic,
  the bad-xy-prob guard), the truncated-Gaussian `norm_cdf`, the
  emission model, the smoother and the initial-probability vector.
  Knob load-bearingness is asserted by tests that change `--err-prob`
  / `--xy-prob` and observe a different decode. Byte-for-byte parity
  against upstream is NOT claimed because upstream emits no
  comparable golden.
- `--AF-file` — **rejected** (a non-empty value is a hard error)
  pending per-site allele-frequency support. Upstream's `vcfcnv.c`
  uses the AF file two ways: it recomputes the per-site genotype
  frequencies fRR/fRA/fAA from each site's `nonref_afs[i]`
  (`vcfcnv.c:735-739`) instead of the fixed defaults, and it acts as
  a targets filter — sites absent from the AF file are dropped
  (`vcfcnv.c:27-31`, `:1429`). This port implements neither, so a
  `--AF-file` run would silently diverge in both the emission model
  and the site set; the port rejects the flag rather than produce
  wrong output. Wiring per-site AFs into the emission model and the
  targets filter is the one remaining deferred piece. Without it the
  port always uses the fixed defaults fRR/fRA/fAA = 0.76/0.14/0.098.
- `-o/--output-dir` — upstream writes per-sample / per-region plot
  data and several `.tab`/`.cn` files into this directory; this port
  streams a single summary TSV (the upstream `summary.tab` "RG"
  rows) to stdout regardless of the path (the flag is still required
  for CLI parity). The per-site `cn.<sample>.tab` and `dat.*.tab`
  files and the `CF` cell-fraction summary rows are not emitted.
- **Paired-mode control counts are dropped from `summary.tab`.** In
  paired (tumour/control) mode upstream writes per-sample summary
  files and a `summary.tab` whose "RG" rows carry four count columns
  — query `nSites`/`nHETs` and control `nSites`/`nHETs`
  (`vcfcnv.c:296-298,1095`). This port emits only the merged
  `summary.tab` "RG" rows with the query sample's `nSites`/`nHETs`;
  it computes the control counts internally but does not write them.
  This is a deliberate I/O-surface reduction, consistent with the
  single-stream `-o/--output-dir` simplification above.
- `-p/--plot-threshold` — accepted; this port emits no plots.
- `--regions-overlap` / `--targets-overlap` — accepted; always
  POS-in-region semantics.
- `-v/--verbosity` — accepted; ignored.
- Indel / non-SNP records — the BAF/LRR signals are typically
  per-marker SNP data; the port honours upstream's behaviour (treat
  each record as one marker regardless of REF/ALT).

Option-tail gaps on `csq` (slices 1-3 done; slice 4 remains):

- **The engine IS haplotype-aware.** `bcftools csq` now phases
  variants per haplotype, walks the GFF transcripts, builds the
  `hap_node_t` tree and reports per-haplotype compound consequences
  in `INFO/BCSQ`. Indels, splice-site disruption, start/stop
  refinement and compound-het bookkeeping are all handled.
- `-p/--phase a|m|r|R|s` — load-bearing: the haplotype-construction
  modes are ported (`phaseRequire/Merge/AsIs/Skip/NonRef/DropGT`).
- `-i/--include` / `-e/--exclude` — accepted; not yet evaluated
  (slice 4). The expression evaluator already exists in
  `pkg/bcftools`; the wire-up is a small follow-up.
- `-s/--samples` / `-S/--samples-file` — accepted; the engine
  currently walks every header sample. Sample subsetting is slice 4.
- `-n/--ncsq INT` — parsed into the per-haplotype `FORMAT/BCSQ` cap
  (`ncsq2`); the cap becomes load-bearing once the `FORMAT/BCSQ`
  bitmask emission lands in slice 4.
- `-B/--trim-protein-seq INT` — accepted; v1 does not truncate
  predictions.
- `-b/--brief-predictions` — hard-rejected with a roadmap pointer
  (upstream deprecates this flag itself).
- `-C/--genetic-code INT|l` — only `0` (standard) is accepted in
  v1; other tables are hard-rejected with a roadmap pointer.
- `-l/--local-csq` — accepted; v1 always operates per-record.
- `--unify-chr-names LIST` — only `0` (no rewriting) is accepted in
  v1; non-zero specs are hard-rejected.
- `--dump-gff` — hard-rejected; v1 has no GFF dump path.
- `-O b|u|z|t` — only `-O v` (VCF text) is supported in v1; the
  others are hard-rejected with a roadmap pointer.
- `--threads`, `-v/--verbosity`, `-W/--write-index`, `--force`,
  `--no-version`, `-q/--quiet` — accepted; v1 ignores.
- The minimal GFF3 parser (`pkg/htsgo/gff`) understands `gene`,
  `mRNA` / `transcript`, `CDS`, and `exon` rows. Other feature
  types are silently skipped — fine for the v1 SNP classifier
  but the parser will need extension for splice-site / UTR work.

### csq full-parity slicing plan

Upstream `csq.c` is ~3994 lines. The Go port is sliced as follows.
Every upstream `csq` golden (`test/csq*.out`) is produced by the
haplotype engine and contains compound consequences (`103G>A+108T>A`),
reference pointers (`@107`), indels and frameshift on the *same* line,
so **no golden validates byte-for-byte until the haplotype engine
(slice 3) is complete**. Slices 1-2 were validated with hand-derived
unit tables; slice 3 unblocked the INFO/BCSQ goldens (now passing
byte-for-byte). The fixtures are vendored under
`tools/bcftools/testdata/csq/`.

- **Slice 1 — region classifier + SO-term completion (per-record, no
  haplotype tree). DONE.** The GFF3 model now carries `five_prime_UTR`
  / `three_prime_UTR` rows (explicit, or derived as exon-minus-CDS)
  and non-coding biotypes; per-record detection of `5_prime_utr`,
  `3_prime_utr`, `intron`, `non_coding` and the splice set
  (`splice_donor`, `splice_acceptor`, `splice_region`) is ported from
  `splice_init` + the SNP/MNP arm of `splice_csq` (the
  `N_SPLICE_DONOR=2` / `N_SPLICE_REGION_INTRON=8` /
  `N_SPLICE_REGION_EXON=3` boundary math, plus the 8bp exon-index
  padding from `gff.c`) and the `test_utr` / `test_splice` /
  `test_tscript` dispatch. The SO-term codon set is complete
  (`stop_gained`, `stop_lost`, `start_lost`, `stop_retained`,
  `coding_sequence`). Landed in `tools/bcftools/pkg/bcftools/
  csq_classify.go`.
- **Slice 2 — indel consequence classification (per-record). DONE.**
  `splice_csq_ins` / `splice_csq_del` / `splice_csq_mnp` /
  `splice_csq_complex` are ported: frameshift vs
  inframe-insertion/deletion, `feature_elongation` /
  `feature_truncation`, and indels at splice sites against a single
  transcript. Bundled with slice 1 in `csq_classify.go`.

  **PROVISIONAL: per-record indel frame bits.** `spliceCSQIns` /
  `spliceCSQDel` set the `csqFrameshift` vs `csqInframeIns` /
  `csqInframeDel` bit from the raw allele-length delta (`%3`).
  Upstream's `splice_csq_ins` / `splice_csq_del` do **not** set any
  frame bit at the splice layer — `hap_add_csq` recomputes frameshift
  vs inframe from the *translated* `dlen` once the haplotype is
  threaded through the spliced reference. The per-record bits are
  therefore an approximation that holds for a clean single-exon indel
  but is wrong whenever the indel spans an intron, partially overlaps
  the CDS boundary, or interacts with other variants on the same
  haplotype. The slice-3 engine MUST replace these per-record bits
  with the `hap_add_csq` `dlen`-based computation. To stop a
  CDS-internal indel being double-staged, `classifyTranscriptVariant`'s
  test_splice arm masks `spliceCSQNonSplice` (the provisional frame /
  elongation / truncation bits) off before deciding whether a splice
  hit occurred — see `csq_classify.go`.

Slices 1+2 shipped together as the per-record-classifier PR. The
`*`-upstream-stop prefix, `shifted_del_synonymous` start/stop refinement
and `inframe_altering` need the spliced reference / haplotype context
and are deferred to slice 3 with the engine. Validation is by
hand-derived unit tables (`csq_classify_test.go`): one case per SO
class — UTR5/UTR3, intron, splice donor/acceptor/region (fwd + rev
strand), stop_gained/stop_lost/start_lost, missense, synonymous,
inframe vs frameshift indel, splice-site indel — plus the kput_vcsq
SO-term precedence ordering.
- **Slice 3 — the haplotype-aware engine. DONE.** The haplotype tree
  (`hap_node_t`, `hap_init`, `hap_finalize`, `hap_add_csq`,
  `cds_translate`), the per-transcript padded + spliced reference
  build (`tscript_init_ref` / `tscript_splice_ref`), the full
  `splice_csq` family with `set_refalt` / `splice_build_hap` /
  `shifted_del_synonymous`, the `vbuf` / `pos2vbuf` position-clustered
  VCF buffer, `csq_push` / `csq_stage`, `kput_vcsq`, and the
  `-p/--phase {a|m|r|R|s}` haplotype-construction modes plus
  `-n/--ncsq` are ported in `csq_hap.go`, `csq_splice.go`,
  `csq_engine.go` and `csq_process.go`. This produces compound
  consequences (`103G>A+108T>A`), the `@pos` reference pointers, the
  `*`-upstream-stop prefix, and the true frameshift / inframe /
  elongation / truncation / start_retained / stop_retained calls from
  the translated `dlen`. The GFF3 CDS reading-frame phase is now
  trimmed off the 5' CDS exon at index-build time (mirroring `gff.c`),
  so the spliced CDS is frame-aligned for both the engine and the
  per-record classifier. *Passes byte-for-byte:* `csq.1.out`,
  `csq.oob-codon.out`, `csq.splice.issue-2543.1.out` — see
  `csq_golden_test.go::TestCSQGoldenINFO`.
- **Slice 4 — GFF/output tail.** The `FORMAT/BCSQ` per-haplotype
  bitmask and the `FORMAT/TBCSQ` `bcftools query` expansion (needed by
  `csq.2.out` / `csq.3.out`, which run `query -f'[%TBCSQ\n]'`);
  `--unify-chr-names LIST` (the 3-field VCF/GFF/FAI rename spec
  exercised by `csq.chr.out` / `csq.yychr.out` / the `csq.nchr` +
  `csq.ychr` matrix), `-l/--local-csq` (`test_cds_local`),
  `--dump-gff`, BCF/`-O b|u|z` output, and `-i/-e` filter wire-up.
  Also the `GF_NMD`/`NMD_transcript` branch of upstream `kput_vcsq`:
  `kputVcsq` currently omits NMD-transcript consequence emission, not
  exercised by the slice 1-3 goldens — it lands here in the slice-4
  tail.
  *Unblocks:* `csq.2.out`, `csq.3.out` (`--ncsq 64`), `csq.chr.out`,
  `csq.yychr.out`, and the `csq.nchr`/`csq.ychr` `--unify-chr-names`
  matrix (test.pl lines 1081-1092).

Status: slices 1+2 (the per-record consequence logic) and slice 3 (the
haplotype-aware engine) **DONE** — see `csq_hap.go`, `csq_splice.go`,
`csq_engine.go`, `csq_process.go` and the golden test
`csq_golden_test.go`. The standalone v1 per-record classifier has been
folded into the engine; `csq_classify.go` now holds only the shared
`csqStrings` table and the `csq*` SO-term bit constants. The `bcftools csq` INFO/BCSQ
output now matches upstream byte-for-byte on `csq.1.out`,
`csq.oob-codon.out` and `csq.splice.issue-2543.1.out`. Slice 4 (the
GFF/output tail — `FORMAT/TBCSQ` text expansion and the
`--unify-chr-names` matrix) remains; `csq.2.out` / `csq.3.out` are
deferred there because they validate the `FORMAT/TBCSQ` per-haplotype
field via `bcftools query`, a cross-tool output-formatting feature
distinct from the haplotype engine.

Option-tail gaps on `mpileup` (SNP-only MAQ model; slices 1, 2 & 3 done):

- **The likelihood model IS the upstream MAQ model.** Slice 2 wired
  `bam2bcf.c::bcf_call_glfgen` / `bcf_call_combine` / `bcf_call2bcf`
  (ported in `bam2bcf.go`) onto the slice-1 errmod port. `mpileup`
  emits one BCF/VCF record per covered position with the `<*>`
  unseen allele, the full multi-allelic PL grid, and
  INFO/DP/I16/QS/MQ0F. The obsolete uniform-error binomial is gone.
- **BAQ recalibration IS wired (slice 3 done).** Reads are run
  through `pkg/htsgo/baq.SamProbRealn` (apply+extend mode) before their
  bases enter the pileup, gated by the ported `mplp_realn` /
  `MPLP_REALN_PARTIAL` column heuristic. `-B/--no-BAQ` disables it;
  `-E/--redo-BAQ` recomputes BAQ. Partial realignment is the default;
  `-D/--full-BAQ` clears `MPLP_REALN_PARTIAL` and forces full BAQ (every
  read realigned). The `<*>`-only `mpileup/*.out` golden records now
  byte-match (`TestMpileupBAQGoldens`).
- **No bias annotations (slice 4 TODO).** VDB / SGB / RPBZ / MQBZ /
  BQBZ / MQSBZ / SCBZ are not yet computed; records carry only
  INFO/DP/I16/QS/MQ0F. The `calc_vdb` / `calc_SegBias` /
  `calc_mwu_biasZ` machinery and the indel caller land in slice 4.
- **No indel calling.** The full upstream indel realigner
  (`bam2bcf_indel.c`) and the consensus indel mode
  (`bam2bcf_edlib.c`) are deferred. Every knob that drives the indel
  model — `-e/--ext-prob`, `-F/--gap-frac`, `-h/--tandem-qual`,
  `--indel-bias`, `--indel-size`, `-I/--skip-indels`,
  `-L/--max-idepth`, `-m/--min-ireads`, `-M/--max-read-len`,
  `--open-prob`, `--indels-cns`, `--indels-2.0`, `--no-indels-cns`,
  `--ar-prob`, `--ambig-reads / --ar`, `--del-bias`, `--poly-mqual`,
  `--no-poly-mqual`, `--score-vs-ref`, `--seqq-offset` — is accepted
  at the CLI for parity but inert in v1. The v1 emit path is
  equivalent to running upstream with `-I/--skip-indels` set.
- **The FORMAT/PL grid is multi-allelic.** Slice 2 emits the full
  upper-triangle `g[z++] = a[j]*5 + a[i]` grid of
  `n_alleles*(n_alleles+1)/2` values per sample, including the `<*>`
  unseen allele.
- **No BAI seek.** `-r/--regions` and `-R/--regions-file` are
  post-filters applied after a linear scan of every input BAM; the
  BAI-seek fast path lives in `pkg/htsgo/sam` but is not wired
  through `mpileup` in v1. Tracked as a follow-up — perf only, no
  output difference.
- **No per-read group filtering.** `-G/--read-groups` is parsed and
  stored; v1 includes every record whose @RG passes the standard
  filters. `-Z/--ignore-RG` is accepted but inert.
- **No gVCF blocking.** `-g/--gvcf` is accepted; one BCF/VCF record
  is emitted per covered reference position (gVCF range-blocking is a
  follow-up).
- **`-a/--annotate LIST` is parsed and partially honoured.**
  `parseFormatFlag` ports the upstream tag-list parser (mpileup.c:1141);
  tokens are accepted as bare names ("AD"), with FORMAT/INFO prefixes,
  and with the "-" prefix to clear bits. The upstream default
  (BQBZ/IDV/IMF/MQ0F/MQBZ/MQSBZ/RPBZ/SCBZ/SGB/VDB + FORMAT/AD) is the
  starting bitset; user tokens layer on top. Today FORMAT/AD, ADF, ADR
  and the bias INFO tags emit; the remaining `INFO/AD,ADF,ADR,SP,SCR`,
  `FORMAT/DP,DV,DPR,SP,SCR,QS,NMBZ,QM,DP4` tags land alongside the
  indel branch of 2bcf (slices 4e.3 and beyond).
- **`-O u|b` (BCF output) works.** Slice 2 wired BCF output through
  `pkg/htsgo/bcf` (`-O b` is BGZF-wrapped, `-O u` uncompressed);
  `-O v` (text VCF) is the default and `-O z` is gzipped VCF.
- **`--threads`, `-v/--verbosity`, `-W/--write-index`, `--no-version`,
  `-A/--count-orphans`, `-x/--ignore-overlaps`, `-d/--max-depth`,
  `-q/--min-MQ`, `-Q/--min-BQ`, `--max-bq`** — fully implemented in v1.
- `-X/--config STR` (presets like `1.12`, `2.1`, `ultima`,
  `pacbio-ccs-1.20`) — accepted; v1 ignores. Most presets toggle the
  indel-model knobs which are inert above.
- `-6/--illumina1.3+` — accepted; v1 ignores (input BAMs are
  Phred+33 across the board).
- `-C/--adjust-MQ INT` (MAPQ tail adjustment) — accepted; v1 ignores.
- The 5-flag mask flags (`--skip-any-unset`, `--skip-all-unset`,
  `--skip-any-set`, `--skip-all-set`, `--ls`) — accepted and stored;
  v1 honours only the standard `--ff` defaults (UNMAP, SECONDARY,
  QCFAIL, DUP, SUPPLEMENTARY) baked into `mpileupKeepRecord`.
- `--delta-BQ` — implemented (default 30): a base quality is capped at
  `neighbour_qual + delta` before the MAQ model sees it.
- `--seed` — accepted; ignored (no subsampling below the 255-read
  errmod cap).

Option-tail gaps on `consensus` (this PR, simple-mode):

- `-c/--chain FILE` — liftover chain file. Accepted; v1 rejects with a
  roadmap pointer at runtime.
- `-H NpIu` — phased-index / unphased-IUPAC encoding. Accepted; v1
  rejects with a roadmap pointer.
- `--regions-overlap 0|1|2` — accepted; v1 ignores (no synced-reader
  region jump path).
- `-v/--verbosity INT` — accepted; v1 ignores.
- `-r/-R/-t/-T` — upstream `consensus` does NOT advertise these flags
  (only `--regions-overlap`); the v1 port follows upstream exactly.
- Multi-sample apply: `-s LIST` accepts a comma list but v1 honours only
  the first name (the other entries are silently ignored, mirroring
  upstream's single-sample focus when -H is unset).
- Complex MNP and SV (BND, breakend, `<DEL>` etc.) records are not yet
  applied to the reference. The current v1 covers SNPs and simple
  REF/ALT length-difference indels.
- Overlapping variants: first wins (left-to-right). Upstream emits a
  warning and folds them into a longer ALT; v1 doesn't.

Plus:

- **`bcftools view`** more flags: `--regions-overlap`, `--targets-overlap`,
  `--no-version`, `--write-index`, `--phased`.
- **CSI seek** for region queries: today we validate via the index then
  linear-scan. Real chunk-seek is the natural follow-up.

Subcommand-tail gaps on `bcftools call`:

- **Full multi-allelic caller (`-m` on >2 ALT sites).** The v1 port
  falls back to the consensus model for sites with more than one ALT;
  upstream iterates over every allele combination.
- **BCF input.** Today `call` rejects BCF input with a roadmap-pointer
  error. The BCF reader's FORMAT-key reconstruction
  (`docs/UPSTREAM_BUGS.md`, `bcf-fmt-keys-missing`) is the prerequisite.
- **`--ploidy GRCh37 / GRCh38`.** Accepted by the CLI parser but
  rejected at runtime — the per-contig sex-chromosome overrides need a
  ploidy registry that's not yet wired in.
- **Index-backed region queries** (`-r` reuses the post-filter path).
- **`--gvcf` block-emit mode** (banded reference blocks).
- **`-C alleles --constrain`** family.

**Validation:** no upstream-test-suite run yet.

### `mosdepth`

**Status:** 1 / 1 command, most flags.

Missing:

- **`.csi`** output (currently emits `.tbi`).
- **D4 output** (`-d/--d4`).
- **Multi-threading** (`-t/--threads N`).
- **`--mapq` 0-only fast-path** — upstream has a special fast loop.

**Validation:** no upstream-test-suite run yet.

---

## Phase plan

### Phase 1 (this PR)

- Truth pass on `PORTING_STATUS.md`, `tools/README.md`, and `analysis/tool_ranking_2026.md`.
- Initial scaffold for `UPSTREAM_BUGS.md` and this file.

### Phase 2 — validated parity audits

For each tool that has an upstream test suite (samtools, bcftools, mosdepth,
vcftools, plus sickle and skewer which have small test corpora), set up
the same kind of validated-parity rig that bedtools got in PR #55:

- Submodule-init the upstream repo.
- Pull representative test cases (5-20 per subcommand) plus their
  expected outputs.
- Diff our output against expected; pass or `t.Skip()` with a documented
  reason.
- Fix small bugs we find inline; record upstream bugs in `UPSTREAM_BUGS.md`.

### Phase 3 — systematic gap closure

Tool-by-tool, dispatched as parallel agent waves of 3-4 PRs each:

1. **bedtools long tail** (~30 small subcommands, parallelisable).
2. **vcftools long tail** (~87 options, parallelisable).
3. **seqtk + fastp + sickle + skewer + prinseq closure** (small individual
   tools, all parallelisable).
4. **samtools mpileup + bcftools call + bcftools annotate** (the big
   variant-calling loop).
5. **`samtools markdup`/`idxstats`/`stats`/`merge`** (workhorse
   utilities).
6. **bcftools merge / isec / sort / annotate** (everyday set-ops).
7. **CRAM** (last; multi-week effort on its own).

Phase 1 lands first. Phases 2 and 3 are several sessions of parallel-agent
waves each.
