# bedlinks

Pure-Go reimplementation of `bedtools links`. Reads BED intervals and emits
an HTML page of UCSC Genome Browser links, one row per record.

## Usage

```bash
# Defaults: human / hg18 against the main UCSC mirror
bedlinks -i regions.bed > links.html

# Point at a local mouse mm9 mirror
bedlinks -i regions.bed \
  -base http://mymirror.example.edu \
  -org mouse \
  -db mm9 > links.html
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-i` | `--input` | BED input (default stdin; `-` = stdin). Transparent gzip. |
| `-o` | `--output` | Output file (default stdout; `-` = stdout). |
| | `-base` | UCSC mirror base URL. Default `http://genome.ucsc.edu`. |
| | `-org` | UCSC organism. Default `human`. |
| | `-db` | UCSC build / db. Default `hg18` (matches upstream; the task brief asked for `hg19` but upstream's hard-coded default is `hg18`). |
| `-h` | `--help` | |
| `-v` | `--version` | |

## Output shape

A fixed HTML scaffold (head + intro + opening table tag) followed by one
`<tr>` per BED record, then a fixed footer:

```html
<html>
        <body>
<title>stdin</title>
<br>Firefox users: Press and hold the "apple" or "alt" key and click link to open in new tab.
<p style="font-family:courier">
<table border="0" align="justify"
<h3>BED Entries from: stdin </h3>
<tr>
        <td>
                <a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg18&position=chr1:100-200>chr1:101-200</a>
        </td>
        ... (name/score/strand <td>s for BED4/5/6+)
</tr>
...
</table>
</p>
        </body>
</html>
```

The URL's `position=` query uses 0-based half-open coordinates; the link
text uses 1-based start to match how UCSC displays the locus. Both match
upstream byte-for-byte.

## Per-record column emission

Upstream switches on `bedType` (the number of columns the input has):

| bedType | Extra `<td>` blocks |
| ------- | ------------------- |
| 3       | (link only)         |
| 4       | `name`              |
| 5       | `name`, `score`     |
| 6/9/12  | `name`, `score`, `strand` |

We follow the same switch.

## Deviations from the task brief

- The task brief proposed a single `<A HREF=...>CHR=...; START=...; END=...;
  NAME=...; SCORE=...; STRAND=...</A>` line per record. Upstream emits a
  full HTML page with a `<table>` of `<tr>`/`<td>` rows; we match upstream
  for byte-for-byte parity.
- The task brief proposed default `db=hg19`. Upstream's hard-coded default
  is `hg18`; we follow upstream.

## Validated parity

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md). Upstream ships no
`links/` test subdir, so the three parity cases derive expected output from
`reference_code/bedtools/src/linksBed/linksBed.cpp` directly.
