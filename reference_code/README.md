# Reference Code Directory

This directory contains the original source code for bioinformatics tools identified for improvement in our analysis. These tools have been added as git submodules to allow easy access to their source code for analysis, comparison, and reference during the rewrite process.

## Overview

Based on the comprehensive analysis in [`analysis/top_50_packages_for_improvement.md`](../analysis/top_50_packages_for_improvement.md), we identified 50 bioinformatics packages that would benefit most from code rewrites and improved documentation. Of these 50 packages, 42 are available as open-source projects on GitHub and have been added as submodules. Additionally, 8 more tools were added (3 per user request and 5 from the extended ranked list) to reach a total of 50 submodules.

## Successfully Added Tools (50 total)

The following 42 tools have been successfully added as git submodules:

### Annotation Tools (7)
1. **Augustus** - Gene prediction (C++)
   - Repository: [Gaius-Augustus/Augustus](https://github.com/Gaius-Augustus/Augustus)
   - Path: `reference_code/augustus/`

2. **BRAKER** - Gene prediction pipeline (Perl)
   - Repository: [Gaius-Augustus/BRAKER](https://github.com/Gaius-Augustus/BRAKER)
   - Path: `reference_code/braker/`

3. **EvidenceModeler** - Gene structure combination (Perl)
   - Repository: [EVidenceModeler/EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler)
   - Path: `reference_code/evidencemodeler/`

4. **MAKER** - Genome annotation pipeline (Perl)
   - Repository: [Yandell-Lab/maker](https://github.com/Yandell-Lab/maker)
   - Path: `reference_code/maker/`

5. **PASA** - Transcript assembly (Perl/C++)
   - Repository: [PASApipeline/PASApipeline](https://github.com/PASApipeline/PASApipeline)
   - Path: `reference_code/pasa/`

6. **Prokka** - Prokaryotic annotation (Perl)
   - Repository: [tseemann/prokka](https://github.com/tseemann/prokka)
   - Path: `reference_code/prokka/`

7. **SNAP** - Gene prediction (C)
   - Repository: [KorfLab/SNAP](https://github.com/KorfLab/SNAP)
   - Path: `reference_code/snap/`

### Alignment Tools (5)
8. **Bowtie2** - Fast read alignment (C++)
   - Repository: [BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)
   - Path: `reference_code/bowtie2/`

9. **BWA** - DNA sequence alignment (C)
   - Repository: [lh3/bwa](https://github.com/lh3/bwa)
   - Path: `reference_code/bwa/`

10. **DIAMOND** - Fast protein alignment (C++)
    - Repository: [bbuchfink/diamond](https://github.com/bbuchfink/diamond)
    - Path: `reference_code/diamond/`

11. **minimap2** - Long-read alignment (C)
    - Repository: [lh3/minimap2](https://github.com/lh3/minimap2)
    - Path: `reference_code/minimap2/`

12. **Subread** - Read alignment and counting (C)
    - Repository: [ShiLab-Bioinformatics/subread](https://github.com/ShiLab-Bioinformatics/subread)
    - Path: `reference_code/subread/`

### Assembly Tools (4)
13. **Canu** - Long-read assembly (Perl/C++)
    - Repository: [marbl/canu](https://github.com/marbl/canu)
    - Path: `reference_code/canu/`

14. **IDBA** - Multiple k-mer assembly (C++)
    - Repository: [loneknightpy/idba](https://github.com/loneknightpy/idba)
    - Path: `reference_code/idba/`

15. **Ray** - Parallel assembly (C++)
    - Repository: [sebhtml/ray](https://github.com/sebhtml/ray)
    - Path: `reference_code/ray/`

16. **wtdbg2** - Long-read assembly (C)
    - Repository: [ruanjue/wtdbg2](https://github.com/ruanjue/wtdbg2)
    - Path: `reference_code/wtdbg2/`

### Epigenomics Tools (3)
17. **Bismark** - Bisulfite-seq alignment (Perl)
    - Repository: [FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark)
    - Path: `reference_code/bismark/`

18. **methylKit** - DNA methylation (R)
    - Repository: [al2na/methylKit](https://github.com/al2na/methylKit)
    - Path: `reference_code/methylkit/`

19. **Segway** - Genome segmentation (Python)
    - Repository: [hoffmangroup/segway](https://github.com/hoffmangroup/segway)
    - Path: `reference_code/segway/`

### Metagenomics Tools (2)
20. **MEGAN** - Metagenome analysis (Java)
    - Repository: [husonlab/megan-ce](https://github.com/husonlab/megan-ce)
    - Path: `reference_code/megan/`

21. **mothur** - Microbial ecology (C++)
    - Repository: [mothur/mothur](https://github.com/mothur/mothur)
    - Path: `reference_code/mothur/`

### Phylogenetics Tools (4)
22. **IQ-TREE** - Phylogenetic inference (C++)
    - Repository: [iqtree/iqtree2](https://github.com/iqtree/iqtree2)
    - Path: `reference_code/iq_tree/`

23. **MrBayes** - Bayesian phylogenetics (C)
    - Repository: [NBISweden/MrBayes](https://github.com/NBISweden/MrBayes)
    - Path: `reference_code/mrbayes/`

24. **PhyML** - Maximum likelihood (C)
    - Repository: [stephaneguindon/phyml](https://github.com/stephaneguindon/phyml)
    - Path: `reference_code/phyml/`

25. **RAxML** - Maximum likelihood trees (C)
    - Repository: [stamatak/standard-RAxML](https://github.com/stamatak/standard-RAxML)
    - Path: `reference_code/raxml/`

### Population Genetics Tools (3)
26. **EIGENSOFT** - Population stratification (C)
    - Repository: [DReichLab/EIG](https://github.com/DReichLab/EIG)
    - Path: `reference_code/eigensoft/`

27. **PHASE** - Haplotype reconstruction (C)
    - Repository: [stephens999/phase](https://github.com/stephens999/phase)
    - Path: `reference_code/phase/`

28. **VCFtools** - VCF manipulation (C++/Perl)
    - Repository: [vcftools/vcftools](https://github.com/vcftools/vcftools)
    - Path: `reference_code/vcftools/`

### QC Tools (6)
29. **fastx_toolkit** - FASTA/Q processing (C++)
    - Repository: [agordon/fastx_toolkit](https://github.com/agordon/fastx_toolkit)
    - Path: `reference_code/fastx_toolkit/`

30. **PRINSEQ** - Sequence quality control (Perl)
    - Repository: [uwb-linux/prinseq](https://github.com/uwb-linux/prinseq)
    - Path: `reference_code/prinseq/`

31. **seqtk** - FASTA/Q processing (C)
    - Repository: [lh3/seqtk](https://github.com/lh3/seqtk)
    - Path: `reference_code/seqtk/`

32. **Sickle** - Quality trimming (C)
    - Repository: [najoshi/sickle](https://github.com/najoshi/sickle)
    - Path: `reference_code/sickle/`

33. **Skewer** - Adapter trimming (C++)
    - Repository: [relipmoc/skewer](https://github.com/relipmoc/skewer)
    - Path: `reference_code/skewer/`

34. **Trim Galore** - Quality and adapter trimming (Perl)
    - Repository: [FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore)
    - Path: `reference_code/trim_galore/`

### RNA-seq Tools (2)
35. **Cuffdiff** - Differential expression (C++)
    - Repository: [cole-trapnell-lab/cufflinks](https://github.com/cole-trapnell-lab/cufflinks)
    - Path: `reference_code/cuffdiff/`

36. **Cufflinks** - Transcript assembly (C++)
    - Repository: [cole-trapnell-lab/cufflinks](https://github.com/cole-trapnell-lab/cufflinks)
    - Path: `reference_code/cufflinks/`

37. **limma** - Differential expression (R)
    - Repository: [cran/limma](https://github.com/cran/limma)
    - Path: `reference_code/limma/`

### Variant Calling Tools (4)
38. **Control-FREEC** - CNV detection (C++)
    - Repository: [BoevaLab/FREEC](https://github.com/BoevaLab/FREEC)
    - Path: `reference_code/control_freec/`

39. **DeepVariant** - Deep learning variant caller (Python/C++)
    - Repository: [google/deepvariant](https://github.com/google/deepvariant)
    - Path: `reference_code/deepvariant/`

40. **LoFreq** - Low-frequency variant caller (C)
    - Repository: [CSB5/lofreq](https://github.com/CSB5/lofreq)
    - Path: `reference_code/lofreq/`

41. **LUMPY** - Structural variant caller (C++)
    - Repository: [arq5x/lumpy-sv](https://github.com/arq5x/lumpy-sv)
    - Path: `reference_code/lumpy/`

42. **Strelka** - Small variant caller (C++)
    - Repository: [Illumina/strelka](https://github.com/Illumina/strelka)
    - Path: `reference_code/strelka/`

### Additional Tools (3)

43. **DISCOVAR** - Variant detection assembly (C++)
    - Repository: [bayolau/discovardenovo](https://github.com/bayolau/discovardenovo) (unofficial mirror)
    - Path: `reference_code/discovar/`
    - Note: Unofficial Git tracking of DISCOVAR-denovo from Broad Institute

44. **fineSTRUCTURE** - Population structure analysis (C)
    - Repository: [danjlawson/finestructure4](https://github.com/danjlawson/finestructure4)
    - Path: `reference_code/finestructure/`
    - Note: Includes Chromopainter for chromosome painting analysis

45. **AdmixTools** - Admixture testing and population genetics (C)
    - Repository: [DReichLab/AdmixTools](https://github.com/DReichLab/AdmixTools)
    - Path: `reference_code/admixtools/`
    - Note: Tools to test whether admixture occurred and related population genetics analyses

### Extended List Tools (5 tools to reach 50)

46. **StringTie** - RNA-seq transcript assembly and quantification (C++)
    - Repository: [gpertea/stringtie](https://github.com/gpertea/stringtie)
    - Path: `reference_code/stringtie/`
    - Note: Transcript assembly and quantification for RNA-Seq data

47. **Minia** - Short-read de Bruijn graph assembler (C++)
    - Repository: [GATB/minia](https://github.com/GATB/minia)
    - Path: `reference_code/minia/`
    - Note: Memory-efficient genome assembler based on a de Bruijn graph

48. **PLINK** - Whole genome association analysis toolset (C)
    - Repository: [chrchang/plink-ng](https://github.com/chrchang/plink-ng)
    - Path: `reference_code/plink/`
    - Note: Free, open-source whole genome association analysis toolset

49. **STAR** - Ultrafast RNA-seq aligner (C)
    - Repository: [alexdobin/STAR](https://github.com/alexdobin/STAR)
    - Path: `reference_code/star/`
    - Note: Spliced Transcripts Alignment to a Reference

50. **Manta** - Structural variant and indel caller (C++)
    - Repository: [Illumina/manta](https://github.com/Illumina/manta)
    - Path: `reference_code/manta/`
    - Note: Structural variant and indel caller for mapped sequencing data

## Unavailable Tools (8/50+)

The following 8 tools could not be added as submodules for the reasons listed:

### Proprietary/Commercial (3)
1. **NovoAlign** - Commercial short read aligner, no public repository
2. **GeneMark** - Commercial gene prediction software, license required
3. **FGENESH** - Commercial gene prediction software, license required

### Not Available on GitHub (5)
4. **MaxBin** - Original version on SourceForge, GitHub repos appear to be removed or made private
5. **STRUCTURE** - Proprietary population genetics software from Stanford (web.stanford.edu/group/pritchardlab/structure.html); no public repository available
6. **Homer** - Official distribution at homer.ucsd.edu, no maintained GitHub mirror, distributed as tarball
7. **RAST** - Web service-based annotation system, standalone repository has access restrictions
8. **sourcefind** - Not a bioinformatics tool (appears to be radio astronomy software); request may have been in error

**Note**: vcftools was already included in the original 42 tools. DISCOVAR was successfully added as an unofficial mirror.

## Working with Submodules

### Initializing and Updating Submodules

When you first clone this repository, the submodule directories will be empty. To fetch the actual code:

```bash
# Initialize and fetch all submodules
git submodule init
git submodule update

# Or do both in one command
git submodule update --init --recursive
```

### Updating a Specific Submodule

To update a specific tool to its latest version:

```bash
cd reference_code/<tool-name>
git pull origin master  # or main, depending on the repository
cd ../..
git add reference_code/<tool-name>
git commit -m "Update <tool-name> to latest version"
```

### Updating All Submodules

To update all submodules to their latest versions:

```bash
git submodule update --remote --merge
```

## Purpose and Usage

These submodules serve several purposes in the Bio AI Experiment project:

1. **Reference Implementation**: Provide access to the original implementations for comparison
2. **Code Analysis**: Enable detailed analysis of code quality, structure, and patterns
3. **Testing Data**: Extract test cases and example data from the original implementations
4. **Feature Completeness**: Ensure our rewrites implement all features of the original tools
5. **Bug Reference**: Identify known bugs and issues to fix in the rewrite

## Next Steps

With the reference code in place, the project can proceed to:

1. **Detailed Analysis**: Deep dive into each tool's implementation
2. **Feature Extraction**: Document all features and behaviors
3. **Test Case Development**: Create comprehensive test suites based on original tools
4. **Go Rewrites**: Begin implementing improved versions in Go
5. **Performance Benchmarking**: Compare new implementations against originals

## License Considerations

Each submodule retains its original license. Before using code from any of these tools, please review the license in each respective repository. Most are open source (GPL, MIT, Apache, etc.), but terms vary.

## Contributing

When working with these submodules:

- Do not modify the code within submodule directories
- Document any findings in the `analysis/` directory
- Create issue references in our main repository for bugs found in original tools
- Respect the licenses of each tool

## Repository Structure

```
reference_code/
├── README.md                    # This file
├── annotation/                  # Annotation tools (via submodules)
├── alignment/                   # Alignment tools (via submodules)
├── assembly/                    # Assembly tools (via submodules)
├── epigenomics/                 # Epigenomics tools (via submodules)
├── metagenomics/                # Metagenomics tools (via submodules)
├── phylogenetics/               # Phylogenetics tools (via submodules)
├── population_genetics/         # Population genetics tools (via submodules)
├── qc/                          # QC tools (via submodules)
├── rnaseq/                      # RNA-seq tools (via submodules)
└── variant_calling/             # Variant calling tools (via submodules)
```

---

**Last Updated**: 2025-10-20 (Updated with 5 more tools: StringTie, Minia, PLINK, STAR, Manta)  
**Total Submodules**: 50 ✓ Target Reached  
**Status**: Active Development
