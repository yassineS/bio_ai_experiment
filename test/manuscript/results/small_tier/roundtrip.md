# Round-trip validation — scale small

| check | format | status | detail |
|---|---|---|---|
| bgzf-compress-decompress | BGZF | PASS | byte-identical compress→decompress; cross-decodes upstream bgzip |
| bam-reencode | BAM | PASS | records identical after BAM re-encode/decode |
| bam-via-cram | CRAM | PASS | BAM→CRAM→BAM agrees with upstream |
| vcf-via-bcf | BCF | PASS | data rows identical after VCF→BCF→VCF |
| fastq-idempotent | FASTQ | PASS | seqtk seq idempotent |
| bgzf-interop-bidirectional | BGZF | PASS | BGZF interop both ways (ours↔upstream) byte-identical |
| bam-interop-bidirectional | BAM | PASS | BAM interop both ways (ours↔upstream) records identical |
| cram-interop-bidirectional | CRAM | PASS | CRAM interop both ways (ours↔upstream) records identical per fixed reader |
| vcfgz-interop-bidirectional | VCF.gz | PASS | VCF.gz interop both ways (ours↔upstream) records identical |
| bcf-interop-bidirectional | BCF | PASS | BCF interop both ways (ours↔upstream) records identical |
| fastq-interop-bidirectional | FASTQ | PASS | FASTQ-over-BGZF interop both ways (ours↔upstream) byte-identical |
| bai-index-interop | BAI | PASS | BAI index interop both ways: region chr2 queried identically |
| csi-index-interop | CSI | PASS | CSI index interop both ways: region chr2 queried identically |
| tbi-index-interop | TBI | PASS | TBI index interop both ways: region chr2 queried identically |
