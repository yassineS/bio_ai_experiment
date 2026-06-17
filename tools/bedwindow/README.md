# bedwindow - Windowed BED Overlap

A Go re-implementation of `bedtools window` (aka `windowBed`). For each feature
in A it examines a "window" — A's interval expanded by `-w` (or asymmetrically
by `-l`/`-r`) base pairs on each side — and reports every feature in B that
overlaps that window. For each overlap the entire A and B records are reported.

## Build

```bash
go build ./tools/bedwindow/cmd/bedwindow
```

## Usage

```bash
bedwindow -a A.bed -b B.bed [options]
```

## Options

| Option | Description |
|--------|-------------|
| `-a, --input-a FILE` | Input BED file A (required) |
| `-b, --input-b FILE` | Input BED file B (required) |
| `-o, --output FILE` | Output file (default: stdout) |
| `-w, --window INT` | Bp added upstream and downstream of each A entry (**default 1000**) |
| `-l, --left INT` | Bp added upstream (left of) each A entry (default 1000) |
| `-r, --right INT` | Bp added downstream (right of) each A entry (default 1000) |
| `-sw` | Define `-l`/`-r` relative to A's strand (swap for `-` strand) |
| `-sm` | Only report B hits on the **same** strand as A |
| `-Sm` | Only report B hits on the **opposite** strand to A |
| `-u, --unique` | Write each A entry once if it has any B overlap |
| `-c, --count` | Append the B-hit count to each A entry (0 included) |
| `-v, --invert` | Report only A entries with **no** B overlap |
| `-wa, --write-a` | (extension) Write the original A entry only |
| `-wb, --write-b` | (extension) Write the original B entry only |
| `-h, --help` | Show help |
| `--version` | Show version |

`-w` cannot be combined with `-l`/`-r`, and `-l`/`-r` must be given together —
matching upstream. `-sm` and `-Sm` are mutually exclusive.

`-wa`/`-wb` are Go-only conveniences; upstream `bedtools window` has no such
flags (its default already prints the full A and B records).

## Parity notes

This port matches `bedtools window` (v2.31.1) byte-for-byte, including the
behaviours that earlier diverged:

- **Window is added to A, not B.** Upstream expands the *A* feature's
  coordinates and then queries the B database. This is observable for the
  asymmetric `-l`/`-r` (and strand `-sw`) cases: `-l` grows the window upstream
  (lower coordinates) of A, `-r` downstream. A B record downstream of A is only
  reached by `-r` (or `-w`), never by `-l` alone.
- **Default window is 1000 bp.** With no `-w`/`-l`/`-r`, each A interval is
  searched with a 1000 bp window on each side. (The previous port defaulted to
  0.)
- **Per-A B-hit order follows upstream's UCSC bin traversal.** B is indexed in
  the UCSC binning tree by its *original* coordinates; for each A the hits are
  emitted finest-bin-level first, then by bin number ascending, then in B-file
  order — not plain file order and not B-start order. See `binorder.go`.
- **Full columns round-trip verbatim.** Records are kept as their raw input
  text, so BED12 B records (and any extra columns) are echoed unchanged. (The
  previous port re-rendered B from a typed record and truncated BED12 block
  columns to 6 fields.)
- The `-c` (count) and `-v` (invert) A-only outputs are unchanged and remain
  byte-exact.

## Testing

```bash
go test ./tools/bedwindow/...
```

The parity tests assert byte-for-byte equality against the live upstream
`bedtools` binary (built from the `reference_code/bedtools` submodule) for the
bin hit-order, default-window, BED12-B, and strand cases; they `t.Fatalf`
(never skip) when the binary is unavailable. Binary-free `TestUnit*` cases pin
the bin comparator.
```
