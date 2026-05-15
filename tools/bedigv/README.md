# bedigv

Pure-Go reimplementation of `bedtools igv`. Reads BED intervals and emits an
IGV batch-mode script that drives one snapshot per interval.

## Usage

```bash
# Defaults: png snapshots into the current dir
bedigv -i regions.bed > snapshots.batch

# Custom output dir + an existing IGV session to load first
bedigv -i regions.bed -path /data/snaps -sess my.session.xml > snapshots.batch

# Slop the locus by 50 bp, sort by base, collapse reads, write SVG
bedigv -i regions.bed -slop 50 -sort base -clps -img svg > snapshots.batch
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-i` | `--input` | BED input (default stdin; `-` = stdin). Transparent gzip. |
| `-o` | `--output` | Output file (default stdout; `-` = stdout). |
| | `-path` | Snapshot directory. Emitted as `snapshotDirectory <path>`. Default `./`. |
| | `-sess` | IGV session file. When set, a `load <path>` line is emitted. |
| | `-sort` | BAM sort directive: `base`, `position`, `strand`, `quality`, `sample`, `readGroup`. Default: no sort. |
| | `-clps` | Emit a `collapse` line per record. |
| | `-name` | Append the BED record's name column to the snapshot filename. Errors if the name is empty. |
| | `-slop` | Flank each interval by N bp in the `goto` locus (filename keeps the original coords). Default 0. |
| | `-img` | Snapshot extension: `png`, `eps`, `svg`, `jpg`. Default `png`. |
| `-h` | `--help` | |
| `-v` | `--version` | |

## Output shape

For each BED record the port emits (matching upstream `bedToIgv.cpp`):

```text
snapshotDirectory <path>
[load <session>]
goto <chrom>:<start-slop>-<end+slop>
[sort <type>]
[collapse]
snapshot <chrom>_<start>_<end>[_<name>][_slop<N>].<img>
```

## Behaviour

- Coordinates in the `goto` locus are `start-slop` / `end+slop`. Upstream
  emits the raw arithmetic result and so do we; if `start-slop` underflows
  past zero the locus contains a negative number (matches upstream).
- The snapshot filename uses the *original* (non-slopped) coordinates,
  optionally suffixed with `_<name>` (when `-name`) and/or `_slop<N>`
  (when `-slop > 0`).
- `-name` with an empty BED-name column is an error (upstream exits 1
  with the same intent).

## Deviations from the task brief

- The task brief proposed a `-name {NUM|POS}` value flag. Upstream's
  `-name` is a boolean: when set it appends the BED name column to the
  snapshot filename. We follow upstream so the output is byte-for-byte
  compatible.
- The task brief described an alternative output header
  (`new` / `genome ...` / `load ...` / `maxPanelHeight ...`). Upstream
  emits only `snapshotDirectory <path>` (and `load <sess>` when `-sess`
  is given). We match upstream.
- The task brief described `-clps` as "collapse single-record output
  into a single goto". Upstream emits one `collapse` line per record
  immediately after `goto` (and after `sort` when both are set). We
  match upstream.

## Validated parity

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md). Upstream ships
no `igv/` test subdir, so the three parity cases derive expected output
from `reference_code/bedtools/src/bedToIgv/bedToIgv.cpp` directly.
