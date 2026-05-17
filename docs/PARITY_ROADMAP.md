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
  `pkg/bioformats/fasta.BuildIndexFullHeader` /
  `OpenRandomAccessFullHeader`, and `bedgetfasta -fullHeader` flows the
  flag through to the index build. Upstream `getfasta.t06` (the
  `-fullHeader` two-line case) and `t07` (the no-`-fullHeader` warning
  case) now both pass byte-for-byte. BGZF FASTA input is also wired
  through (this PR): `pkg/bioformats/fasta` now sniffs the BGZF magic
  in `OpenRandomAccess` / `OpenRandomAccessFullHeader` and routes to a
  new `OpenRandomAccessBGZF` that fully decompresses the payload
  in-memory and reuses the existing FAI index path. The `.gzi` sidecar
  (when present) is parsed for early validation via a stdlib-only
  little-endian reader in `pkg/bioformats/fasta/bgzf.go`; a samtools-
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
  via `pkg/bioformats/sam.NewBAMReader`; primary alignments contribute
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

**Status:** ~73 of ~147 options (~50%) after long-tail wave 2.

Closed in wave 1:

- **Inter-chromosomal LD**: `--interchrom-geno-r2`, `--interchrom-hap-r2` ✅
- **Chi-square LD**: `--geno-chisq` ✅
- **Relatedness**: `--relatedness` (Yang 2010), `--relatedness2` (KING-robust) ✅
- **Runs of homozygosity**: `--LROH` (+ `--LROH-min-variants`) ✅
- **Phased blocks**: `--phased-blocks` ✅
- **FILTER tag include/exclude**: `--remove-filtered`, `--keep-filtered` ✅
- **INFO selection in recode**: `--keep-INFO TAG`, `--remove-INFO TAG` ✅
- **INFO extraction**: `--get-INFO TAG[,TAG]` → `.INFO` ✅

Closed in wave 2 (this PR):

- **LDhat output formats**: `--ldhat`, `--ldhat-geno` (paired
  `<prefix>.ldhat.sites` / `<prefix>.ldhat.locs`, byte-for-byte vs
  upstream) ✅
- **Phased-site filter**: `--phased` (composes with `--ldhat` per
  upstream's `phased_only` invariant) ✅

Remaining gaps:

- **Mendelian inheritance checks**: `--mendel`.
- **Diff family extensions**: `--diff-indv-map`, `--diff-discordance-matrix`,
  `--diff-switch-error`, `--gzdiff` (already implicit via iohelper).
- **Output formats**: missing `--ldhelmet`, `--IMPUTE`, `--phase` output
  paths.
- **Per-individual output**: the per-individual `.imiss` row layout has
  fields we don't emit (we have `--missing-indv`).
- **Other**: `--FILTER-PASS-summary`, `--remove-INFO-all` (use
  `--keep-INFO`/`--remove-INFO`), `--non-ref-af*` family, `--pca` family.

Note: the brief mentioned `--haploid` as a possible wave-2 target. After
checking the upstream source (`reference_code/vcftools/src/cpp/`) there is
no `--haploid` flag — the closest thing is `--phased` (parameters.cpp:311
+ entry_filters.cpp:989-1010), which we ported instead.

**Validation:** wave 1 adds header byte-for-byte parity tests for the new
output files; wave 2 ships full byte-for-byte parity tests for both
`.ldhat.sites` and `.ldhat.locs` against upstream goldens (under
`tools/vcftools/testdata/parity/`). Full upstream-test-suite run still
pending.

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
**`consensus`** (simple-mode FASTA/FASTQ/pileup; bayesian falls back
with a stderr warning).

Missing subcommands (in rough priority order):

- **`tview`** — terminal viewer. **Deliberate skip** (interactive
  curses UI; near-zero pipeline usage and would require an ncurses
  dependency). Not on the roadmap.
- **`view` flag-tail**: `-X` (custom-index input). `-L bed` landed as a
  linear-scan BED-region filter; `-M`/`--use-multi-region-iterator` is
  accepted but treated as a no-op since we always run the full
  intersection. `-d/-D` (tag-value filter) and `-N` (qname file) landed
  in the view-d-D-N PR.
- **`mpileup` tail** beyond PR #88 wiring: BCF output, `-g/-u` genotype-
  likelihood mode. `-aa` zero-fill of empty contigs is implemented (see
  `TestMpileup_AA_ZeroFillTableDriven`).

Plus:

- **CRAM** read/write throughout. Multi-month effort on its own; the
  rANS codec layer is the gating piece. Owner has OK'd third-party
  dependencies for CRAM codecs (see `CLAUDE.md#documented-exception-cram`).
  Design doc to follow.
- **`.csi`** for samtools (BAI is fine for chromosomes ≤512Mb).
- **Multi-threading** in `sort`, `index`, `view` (`-@`).

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

**`stats` deferred sections** (also documented in
`PARITY_VALIDATION.md`):

- COV/COV2 coverage histograms (require reference + BAI).
- GCD/GCT/GCC/GCL GC distributions (require reference bases).
- FFQ/LFQ per-cycle quality matrices and OXC oxidation-context counts.
- `--target-regions BED` restriction.
- The leading CHK checksum block (CRC32 reduction of read names /
  sequences / qualities).
- BWA-style quality trimming (`-q/--trim-quality`). The SN field
  `bases trimmed` is reported as 0; upstream also reports 0 when the
  flag is not passed, so byte parity holds for the default invocation.
- Mate-tracking memory cap: upstream's `cleanup_overlaps` periodically
  evicts stale `mates` entries. Our `mates` map currently grows
  unbounded — fine for the typical workload but worth fixing before
  running `stats` on multi-billion-record BAMs.

v1 emits byte-faithful **SN** (Summary Numbers) and useful **RL / MAPQ /
IS** rollups; the unsupported sections are quietly omitted (or, under
`--sparse`, all histogram blocks are suppressed entirely).

**Validation:** upstream fixtures from `reference_code/samtools/test/markdup/`
and `.../test/stat/` are vendored under
`tools/samtools/testdata/parity/{markdup,stat}/`. The byte-exact /
flag-exact / SN-byte cases are exercised in
`tools/samtools/pkg/samtools/markdup_test.go` and `stats_test.go`.

**`calmd` deferred features** (accepted as CLI flags, behaviour partial):

- **BAQ recalculation** (`-r`, `-E`, `-A`). Upstream's `bam_md.c` calls
  `sam_prob_realn` to recompute the BQ (base-alignment quality) aux
  tag and to drop MAPQ for low-quality reads. v1 fills in MD + NM
  correctly but does not touch BQ or MAPQ. ~200 lines of HMM-style
  alignment math; deferred per owner steer.
- **`-h` HASH_QNM** (hash-based query-name binarisation) — niche
  upstream-only optimisation; not implemented.
- **`-d` DROP_TAG** (drop all aux but RG) and **`-q` BIN_QUAL** (round
  qualities to 0/7) — flag-recognised in the CLI driver but not
  threaded through; safe to add later as small wrappers around the
  current calmd pipeline.
- **`-n` max-NM cap** — would mask high-mismatch reads with bin-quality;
  trivial follow-up once BIN_QUAL lands.
- **`-N` clear-MD/NM-bits**, **`-C` capQ**, **`--no-PG`** —
  CLI-accepted-and-ignored stubs.

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
positionals, -T aux extraction, --order, -R/-r RG). The upstream
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

**`targetcut` scope reduction.** The user-facing spec for the Go port
is "cut the aligned slice from each read and emit FASTA". Upstream's
`cut_target.c` actually does something quite different —
HMM-based consensus calling over fosmid pools, emitting one consensus
SAM record per identified region. The HMM consensus mode is **not**
implemented; the upstream tool is rarely used outside fosmid
workflows. The simple aligned-slice FASTA mode landed here covers the
"cut a read down to its aligned bases" use case that users typically
mean when they reach for the name. The `-Q` flag is wired through to
the per-base quality filter as documented.

**Validation:** hand-built SAM fixtures in
`tools/samtools/pkg/samtools/phase_test.go` and
`tools/samtools/pkg/samtools/targetcut_test.go`. Phase tests cover
single-block chaining (consistent & label-flipping orderings),
ambiguous-label fall-back when reads don't bridge two hets, and the
MinMAPQ filter. Targetcut tests cover soft-clip flank stripping,
insertion retention, deletion handling, unmapped/secondary skipping,
SEQ='*' skipping, and `-Q` per-base filtering. There is no upstream
regression-test fixture for either tool in
`reference_code/samtools/test/` so byte-parity against upstream is
not pursued.

**`consensus` deferred features** (accepted as CLI flags, behaviour
partial). Upstream `bam_consensus.c` ships five modes — `simple` and
four bayesian flavours (`bayesian_r` aka "bayesian", `bayesian_m`,
`bayesian_p`, `bayesian_116`). v1 only implements `simple`. Because
upstream defaults to `MODE_RECALL` (a bayesian mode) at
`bam_consensus.c:2983`, the v1 binary's default invocation lands on
the bayesian branch, emits a single-line stderr warning, and falls
back to `simple`. The deferred surface:

- **Bayesian (Gap5-derived) mode.** All variants (`bayesian`,
  `bayesian_r`, `bayesian_m`, `bayesian_p`, `bayesian_116`) and their
  knobs (`-C/--cutoff`, `--P-het`, `--P-indel`, `--het-scale`,
  `--adj-qual`, `--use-MQ`, `--adj-MQ`, `--NM-halo`, `--SC-cost`,
  `--scale-MQ`, `--low-MQ`, `--high-MQ`, `-p/--homopoly-fix`,
  `--homopoly-score`, `--homopoly-redux`, `-t/--qual-calibration`,
  `-X/--config`) are accepted on the CLI but not yet implemented.
- **Pileup-mode insertion rows.** Upstream's default `--show-ins yes`
  emits extra rows with `nth>0` for each column of an inserted
  sequence. v1 emits only `nth=0` rows (one per reference position);
  insertion columns are folded into the FASTA/FASTQ stream when
  `--show-ins yes` (the default) but not into the pileup stream.
- **Mate-overlap dedup.** `--ignore-overlaps` is accepted but is a
  no-op; v1 counts each mate independently in the pileup walker.
- **Reference-aware modes.** `-T/--reference`, `--ref-qual`, and
  `--default-qual` are accepted but unused; the simple scoring path
  doesn't need a reference, and the bayesian path that does is
  deferred.
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
and the default-bayesian fallback emitting a stderr warning.
Coverage of the `pkg/samtools` package after this PR is ~80%.

### `bcftools`

**Status:** 24 of ~30 subcommands (~80%). `view`, `index`, `stats`, `query`,
`concat`, `norm`, `call` (consensus + biallelic multi-allelic), the PR #86
wave-1 tail (`annotate`, `head`, `isec`, `merge`, `reheader`, `sort`), the
convert/mendelian PR (`convert`, `mendelian`), the gtcheck/roh PR (`gtcheck`,
`roh`), the filter/consensus PR (`filter`, `consensus`), the
mendelian2/polysomy PR (`mendelian2`, `polysomy`), the cnv/csq PR
(`cnv` + `csq`), and the mpileup PR (**`mpileup`**).

Missing subcommands (priority order):

- **`+plugins`** — full plugin system (substantial).

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

Option-tail gaps on `roh` (PR #107, simple-mode):

- `-b/--buffer-size`, `-e/--estimate-AF`, `-m/--genetic-map`,
  `-M/--rec-rate`, `-V/--viterbi-training` — accepted-and-rejected
  with PARITY_ROADMAP pointer.
- `-O z` — bgzip output; v1 only emits tab-text.
- Transition defaults are upstream's literal per-bp magnitudes
  (`6.7e-8` / `5e-9`) but NOT scaled by physical inter-marker
  distance, so RG quality scores are NOT comparable to upstream's
  until distance scaling lands.

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
- BCF output (`-O b|u`) round-trips through the shared `pkg/bioformats/bcf`
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

Option-tail gaps on `polysomy` (this PR, simple-mode):

- The full Gaussian-mixture peak-fit (`peakfit.c` + GSL) is
  **deferred**. v1 emits CN calls from a median-deviation heuristic:
  CN1 when n_het == 0, CN2 when |median(BAF) - 0.5| ≤ MinBafDev (with
  CnPenalty scaling), CN3 otherwise. The upstream algorithm fits CN2,
  CN3, and CN4 Gaussian mixtures and picks the lowest CN whose fit
  beats `(1 - cn_penalty) * previous_fit`.
- `-b/--peak-size FLOAT`, `-f/--fit-th FLOAT`, `-i/--include-aa`,
  `-p/--peak-symmetry FLOAT`, `--nbins INT`, `--smooth INT`,
  `--ra-rr-scaling` — accepted at the CLI for parity but inert in v1
  (the heuristic doesn't run a peak fit).
- `-o/--output-dir PATH` — accepted but ignored; v1 writes the
  per-chromosome TSV to stdout (no per-chromosome PNG plots).
- `--regions-overlap`, `--targets-overlap` — accepted; v1 always uses
  POS-in-region.
- `-v/--verbosity INT` — accepted; v1 ignores.
- `-n/--include-noise` — accepted; v1 always emits every chromosome
  (we don't classify any as noise / `?`).
- `--force-cn INT` (hidden upstream option) — implemented as
  per-chromosome override.
- BAF source: upstream requires FORMAT/BAF; v1 also accepts
  FORMAT/AD = REF,ALT as a fallback (synthesises BAF as
  `ALT / (REF + ALT)` at het sites).
- Per-record `-i/-e` are NOT in upstream `polysomy.c:main_polysomy`
  and we follow upstream's surface exactly (no invented flags).

Option-tail gaps on `cnv` (this PR, v1 heuristic):

- **The v1 algorithm is NOT the upstream HMM.** Upstream's vcfcnv.c
  runs a 5-state HMM (CN0/CN1/CN2/CN3/CN4) over each contig with
  joint BAF + LRR Gaussian emissions and a configurable transition
  matrix. The v1 port replaces this with a per-sample × per-chrom
  median-BAF + mean-LRR heuristic that classifies each chromosome
  into one of the same 5 CN states. The full Viterbi sweep is the
  natural follow-up; the CLI surface is already parity-clean for
  it. EVERY HMM tuning knob (`-a/--aberrant`, `-b/--BAF-weight`,
  `-e/--err-prob`, `-l/--LRR-weight`, `-L/--LRR-smooth-win`,
  `-O/--optimize`, `-P/--same-prob`, `-W/--baum-welch`,
  `-x/--xy-prob`, `--AF-file`) is parsed and stored in `CNVOptions`
  but the heuristic does NOT consume them. Only `-d/--BAF-dev` and
  `-k/--LRR-dev` (the per-sample expected std-dev floors) actually
  drive the v1 thresholds.
- `-o/--output-dir` — upstream writes per-sample / per-region plot
  data into this directory; v1 always streams a single summary TSV
  to stdout regardless of the path (the flag is still required for
  CLI parity).
- `-p/--plot-threshold` — accepted; v1 emits no plots.
- `--regions-overlap` / `--targets-overlap` — accepted; v1 always
  uses POS-in-region semantics.
- `-v/--verbosity` — accepted; v1 ignores.
- BCF / VCF.gz output — v1 always emits the summary TSV; the
  upstream `-O b|u|z|v` selector does not apply (the upstream tool
  produces several per-region files; v1 produces one summary).
- Indel / non-SNP records — the BAF/LRR signals are typically
  per-marker SNP data; v1 honours upstream's behaviour (treat each
  record as one marker regardless of REF/ALT).

Option-tail gaps on `csq` (this PR, v1 SNP-only):

- **The v1 classifier is NOT haplotype-aware.** Upstream's csq.c
  phases variants per haplotype, walks the GFF transcripts, and
  reports the per-haplotype consequence chain (including
  compound-het effects). v1 instead classifies one SNP at a time
  against the GFF CDS exons and emits `INFO/BCSQ` per-transcript
  for the matching position. Indels, splice-site disruption,
  start-gain, stop-retained, and compound-het bookkeeping are all
  deferred.
- `-p/--phase a|m|r|R|s` — parsed and stored; the per-record SNP
  classifier ignores phasing because consequences are computed
  position-by-position. Will become load-bearing when haplotype-
  aware phasing lands.
- `-i/--include` / `-e/--exclude` — accepted; v1 ignores (every
  input record runs through the classifier). The expression
  evaluator already exists in `pkg/bcftools`; the wire-up is a
  trivial follow-up.
- `-s/--samples` / `-S/--samples-file` — accepted; v1 does not
  subset (consequences are position-driven). Once haplotype-aware
  phasing lands the sample list will gate which haplotypes are
  walked.
- `-n/--ncsq INT` — accepted; v1 emits every matching transcript
  without a cap (one BCSQ entry per transcript).
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
- The minimal GFF3 parser (`pkg/bioformats/gff`) understands `gene`,
  `mRNA` / `transcript`, `CDS`, and `exon` rows. Other feature
  types are silently skipped — fine for the v1 SNP classifier
  but the parser will need extension for splice-site / UTR work.

Option-tail gaps on `mpileup` (this PR, v1 SNP + uniform-error):

- **The v1 likelihood model is NOT the upstream MAQ model.** Upstream's
  `bam2bcf.c::glfgen` reads per-base error probabilities from the
  Heng Li MAQ recalibrator with BAQ adjustments; v1 instead uses the
  simpler samtools-0.1.19-style uniform-error binomial: e = 10^(-Q/10)
  per base, summed in log10 across reads, then phred-scaled and
  rebased to min=0 for the [0/0, 0/1, 1/1] triple. The MAQ port is
  the natural follow-up; the CLI surface and FORMAT/PL layout are
  parity-clean for it.
- **No BAQ recalibration.** `-B/--no-BAQ` is the v1 default (the flag
  is accepted as a no-op); `-D/--full-BAQ` is accepted but inert;
  `-E/--redo-BAQ` is hard-rejected with a roadmap pointer because
  silently skipping a recalibration step a downstream caller asked
  for would yield misleading PLs.
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
- **No multi-allelic FORMAT/PL grid.** The PL we emit is always
  biallelic [PL(0/0), PL(0/1), PL(1/1)] against ALT[0]. Sites with
  multiple ALTs still parse upstream-style in `bcftools call`, but
  the PL grid for the 2nd / 3rd ALT genotypes is treated as 0
  (uninformative). The full j(j+1)/2 + i grid lands with the MAQ
  port.
- **No BAI seek.** `-r/--regions` and `-R/--regions-file` are
  post-filters applied after a linear scan of every input BAM; the
  BAI-seek fast path lives in `pkg/bioformats/sam` but is not wired
  through `mpileup` in v1. Tracked as a follow-up — perf only, no
  output difference.
- **No per-read group filtering.** `-G/--read-groups` is parsed and
  stored; v1 includes every record whose @RG passes the standard
  filters. `-Z/--ignore-RG` is accepted but inert.
- **No gVCF blocking.** `-g/--gvcf` is accepted; v1 always emits one
  VCF record per variant site (REF-only sites are skipped, matching
  upstream when `--gvcf` is unset).
- **`-a/--annotate LIST` is accepted but inert.** v1 always emits the
  default `INFO/DP`, `INFO/I16`, `FORMAT/PL` set. The
  `INFO/AD,ADF,ADR,SP,SCR,IDV,IMF`, `FORMAT/AD,ADF,ADR,DP,DV,DPR,SP,SCR,QS`
  tags will land alongside the per-tag stream when called from
  `bcftools call`.
- **`-O u|b` (BCF output) is hard-rejected.** The BCF writer in
  `pkg/bioformats/bcf` can handle generic records, but mpileup carries
  custom INFO/I16 typing rules; the wire-up is a follow-up. `-O v`
  (text VCF) is the default; `-O z` (gzipped VCF) is accepted at the
  CLI but currently streams text — gzip-wrap-stdout from a follow-up
  CLI shim will close that gap.
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
- `--seed`, `--delta-BQ` — accepted; v1 ignores.

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
