# mosdepth oracle binary (not committed)

The mosdepth live-parity tests (`tools/mosdepth/pkg/mosdepth/live_oracle_test.go`)
compare our port against the genuine upstream binary. mosdepth is written in
Nim and isn't buildable in this environment, so the binary is fetched on demand
rather than committed (it's ~19 MB).

To enable the live-oracle tests locally:

    curl -sL https://github.com/brentp/mosdepth/releases/download/v0.3.14/mosdepth \
      -o reference_code/mosdepth/mosdepth
    chmod +x reference_code/mosdepth/mosdepth

The tests `t.Skip` cleanly when the binary is absent (e.g. in CI).
