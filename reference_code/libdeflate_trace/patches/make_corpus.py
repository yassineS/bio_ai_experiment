#!/usr/bin/env python3
"""
make_corpus.py — vendor a small reproducible test corpus for the
libdeflate-parity slice-0 oracle.

Emits five binary inputs under <out-dir>:

  empty.bin         0 bytes
  single_byte.bin   one byte 0xAB
  repeated_a.bin    100 bytes of 'A'
  random_64k.bin    65536 bytes from a deterministic LCG (seed 0xC0FFEE)
  bgzf_payload.bin  65280 bytes (BGZF_MAX_BLOCK_SIZE - HEADER_OVERHEAD)
                    of pseudo-typical SAM-like ASCII content, also
                    deterministically generated so reruns are byte-exact.

The PRNG is rolled in-line (LCG) rather than using random.* so the corpus
is portable across Python versions — the bytes are 100% reproducible.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def lcg_bytes(n: int, seed: int = 0xC0FFEE) -> bytes:
    """Knuth MMIX-style LCG byte stream; chosen for reproducibility, not
    statistical quality."""
    state = seed & ((1 << 64) - 1)
    out = bytearray(n)
    a = 6364136223846793005
    c = 1442695040888963407
    mask = (1 << 64) - 1
    for i in range(n):
        state = (state * a + c) & mask
        out[i] = (state >> 33) & 0xFF
    return bytes(out)


def bgzf_payload_bytes(n: int, seed: int = 0xBA5E) -> bytes:
    """Produce SAM-flavoured ASCII roughly the size of a BGZF block.

    The content mimics the leading bytes of a typical SAM record stream
    (read names + DNA sequence + base qualities) so the deflate encoder
    exercises realistic literal/match distributions.  Still fully
    deterministic for reproducibility.
    """
    bases = b"ACGT"
    qual = b"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLM"  # Phred+33 0..44
    out = bytearray()
    state = seed & ((1 << 64) - 1)
    a = 2862933555777941757
    c = 3037000493
    mask = (1 << 64) - 1
    rec = 0
    while len(out) < n:
        # Pseudo-record: name<TAB>flag<TAB>chr1<TAB>pos<TAB>60<TAB>76M<TAB>*<TAB>0<TAB>0<TAB>SEQ<TAB>QUAL\n
        state = (state * a + c) & mask
        pos = (state >> 16) % 250_000_000
        rec += 1
        name = f"r{rec:08d}".encode("ascii")
        seq = bytearray(76)
        qb = bytearray(76)
        for i in range(76):
            state = (state * a + c) & mask
            seq[i] = bases[(state >> 23) & 3]
            qb[i] = qual[(state >> 11) % len(qual)]
        line = (
            name + b"\t0\tchr1\t" + str(pos).encode("ascii") +
            b"\t60\t76M\t*\t0\t0\t" + bytes(seq) + b"\t" + bytes(qb) + b"\n"
        )
        out.extend(line)
    return bytes(out[:n])


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    out: Path = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    (out / "empty.bin").write_bytes(b"")
    (out / "single_byte.bin").write_bytes(bytes([0xAB]))
    (out / "repeated_a.bin").write_bytes(b"A" * 100)
    (out / "random_64k.bin").write_bytes(lcg_bytes(64 * 1024))
    # 65280 = BGZF_MAX_BLOCK_SIZE (65536) - typical header/footer overhead.
    (out / "bgzf_payload.bin").write_bytes(bgzf_payload_bytes(65280))

    for name in ("empty", "single_byte", "repeated_a", "random_64k", "bgzf_payload"):
        path = out / f"{name}.bin"
        print(f"make_corpus: {path} ({path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
