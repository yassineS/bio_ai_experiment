#!/usr/bin/env python3
"""Regenerate the FASTQ input fixtures for the fastp parity test.

The fixtures are checked in (so the parity test is reproducible without
needing Python at test time), but this script makes the contents
auditable. Run it from this directory:

    python3 generate.py

Every read is constructed deterministically (random.seed(42)).
"""
import random

random.seed(42)

ADAPTER = "AGATCGGAAGAGC"


def fmt(records):
    out = []
    for n, seq, qual in records:
        out.append(f"@{n}\n{seq}\n+\n{qual}\n")
    return "".join(out)


def write(path, records):
    with open(path, "w") as f:
        f.write(fmt(records))


def hi_q(n):
    return "".join(chr(33 + random.randint(30, 40)) for _ in range(n))


def med_q(n):
    return "".join(chr(33 + random.randint(20, 35)) for _ in range(n))


def lo_q(n):
    return "".join(chr(33 + random.randint(2, 10)) for _ in range(n))


def main():
    # 1) clean SE reads: 30 records, 100bp, high quality, random ACGT.
    clean = []
    for i in range(30):
        seq = "".join(random.choice("ACGT") for _ in range(100))
        clean.append((f"clean_{i}", seq, hi_q(100)))
    write("se_clean.fq", clean)

    # 2) SE reads with the canonical Illumina TruSeq adapter embedded mid-read.
    adapter = []
    for i in range(20):
        insert_len = random.randint(30, 60)
        insert = "".join(random.choice("ACGT") for _ in range(insert_len))
        rest = "".join(
            random.choice("ACGTN")
            for _ in range(100 - insert_len - len(ADAPTER))
        )
        seq = (insert + ADAPTER + rest)[:100]
        adapter.append((f"adp_{i}", seq, hi_q(100)))
    write("se_adapter.fq", adapter)

    # 3) Reads with N's scattered through; tests --n_base_limit boundary.
    ns = []
    for i in range(10):
        seq = "".join(
            random.choice("ACGTN") if random.random() < 0.3 else random.choice("ACGT")
            for _ in range(80)
        )
        ns.append((f"ns_{i}", seq, hi_q(80)))
    write("se_ns.fq", ns)

    # 4) Low-complexity mix: alternating AT-only and full-complexity reads.
    lowcomp = []
    for i in range(10):
        if i % 2 == 0:
            seq = ("AT" * 50)[:80]  # 100% alternation: complexity=1.0.
        else:
            seq = "".join(random.choice("ACGT") for _ in range(80))
        lowcomp.append((f"lc_{i}", seq, hi_q(80)))
    write("se_lowcomp.fq", lowcomp)

    # 5) Low-complexity homopolymer mix: pure homopolymers (complexity=0) plus
    # mixed reads, so the filter actually rejects something.
    homop = []
    for i in range(10):
        if i % 2 == 0:
            seq = ("A" * 80)  # 0% adjacent differences → complexity=0.
        else:
            seq = "".join(random.choice("ACGT") for _ in range(80))
        homop.append((f"hp_{i}", seq, hi_q(80)))
    write("se_homopolymer.fq", homop)

    # 6) Low-quality 3' tail reads (60 hi-Q + 40 lo-Q). For sliding-window
    # cut_right/cut_tail/cut_front parity.
    lq = []
    for i in range(15):
        seq = "".join(random.choice("ACGT") for _ in range(100))
        lq.append((f"lq_{i}", seq, hi_q(60) + lo_q(40)))
    write("se_lqtail.fq", lq)

    # 7) UMI reads: 8bp leading UMI + 72bp insert.
    umi = []
    for i in range(15):
        umi_seq = "".join(random.choice("ACGT") for _ in range(8))
        rest = "".join(random.choice("ACGT") for _ in range(72))
        umi.append((f"umi_{i}", umi_seq + rest, hi_q(80)))
    write("se_umi.fq", umi)

    # 8) Reads for Q20/Q30 filtering: a mix of fully-hi-Q and partially-lo-Q.
    qmix = []
    for i in range(20):
        seq = "".join(random.choice("ACGT") for _ in range(80))
        if i < 10:
            qual = hi_q(80)
        else:
            qual = med_q(40) + lo_q(40)
        qmix.append((f"qm_{i}", seq, qual))
    write("se_qmix.fq", qmix)

    # 9) PE: R1 and R2 are the reverse complement of each other where they
    # overlap (synthetic short-insert pairs).
    r1 = []
    r2 = []
    for i in range(40):
        seq = "".join(random.choice("ACGT") for _ in range(80))
        r1.append((f"pe_{i}", seq, hi_q(80)))
        comp = seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
        r2.append((f"pe_{i}", comp, hi_q(80)))
    write("pe_r1.fq", r1)
    write("pe_r2.fq", r2)

    # 10) Reads explicitly suitable for duplication evaluation: 50 reads
    # where 20 are duplicates of one another.
    dup = []
    template = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"[:60]
    for i in range(50):
        if i < 30:
            seq = "".join(random.choice("ACGT") for _ in range(60))
        else:
            seq = template
        dup.append((f"d_{i}", seq, hi_q(60)))
    write("se_dup.fq", dup)


if __name__ == "__main__":
    main()
