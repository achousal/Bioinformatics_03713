# Cross-Species Regulatory Element Analysis Pipeline

**Team:** Andres Chousal, Zahin Peerzade, Siddharth Sabata, Yinuo Yang

**Course:** Bioinformatics Data Practicum

**Tissues:** Liver, Pancreas (Human \& Mouse)

---

## Table of Contents

- [Introduction](#introduction)
- [Workflow](#workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
- [Usage](#usage)
- [Output Files](#outputs)
- [Pitfalls and Limitations](#pitfalls-and-limitations)

---

## Introduction

This pipeline performs a comparative analysis of transcriptional regulatory elements (promoters and enhancers) between liver and pancreas tissues in human and mouse. It maps open chromatin regions across species, identifies conserved and tissue-specific elements, classifies peaks as promoters or enhancers, and performs motif discovery and functional enrichment.

We aim to answer three primary questions with these results:

**Question 1**: Is transcriptional regulatory element activity more conserved across tissues or species?

**Question 2**: To what extent does the transcriptional regulatory code differ between enhancers and promoters?

**Question 3**: To what extent are the biological processes upregulated in tissue conserved across species?

(*Questions directly quoted from project description document*)

---
## Workflow

![](https://github.com/achousal/Bioinformatics_03713/blob/main/pipeline.png)
1. **Data Quality Control**
ATAC-seq quality reports (provided) from human and mouse liver, pancreas, and ovary were used to determine which datasets to use for further analysis. The ovarian dataset consisted of short read lengths for both human and mouse, and was thrown out. This step is essential for ensuring informative results. 
2. **Cross-Species Mapping**
Open chromatin regions identified from ATAC-seq are mapped between human and mouse genomes using HALPER (HAL Liftover Post-processing for Epigenomic Regions) with Cactus whole-genome alignments. This step is crucial for identifying orthologous regulatory elements between species, allowing direct comparison of regulatory activity at corresponding genomic locations. The mapping process accounts for genomic rearrangements and evolutionary changes that have occurred since the divergence of humans and mice. 
3. **Conservation and Specificity Analysis (Question 1)**
This step identifies and categorizes regulatory elements based on their conservation patterns: elements conserved between species for the same tissue, elements shared across tissues within a species, and elements specific to a tissue or species. By quantifying these different categories, we can directly address whether regulatory element activity is more conserved across tissues or species. This analysis provides insights into the evolutionary constraints on gene regulation.
4. **Functional Annotation (Question 3)** 
Regulatory elements identified in previous steps are annotated and analyzed for GO term enrichment with ChIPseeker. This analysis connects regulatory elements to their potential target genes and biological functions, revealing which biological processes are regulated by conserved versus species-specific elements. The functional analysis helps determine whether similar biological processes are regulated by conserved elements across species, despite potential differences in the specific regulatory elements.
5. **Enhancer/Promoter Classification (Question 2)** 
Open chromatin regions are classified as promoters (within 2kb upstream and 200bp downstream of TSS) or enhancers (all other regions) using BEDTools and genome annotations. This classification is essential for understanding how conservation patterns differ between these two types of regulatory elements. Promoters, which are proximal to genes, may be under different evolutionary constraints than enhancers, which can act at a distance and may evolve more rapidly.
6. **Motif Discovery (Question 2)** 
Sequences from classified regulatory elements are extracted and analyzed using MEME-ChIP for de novo motif discovery. This step identifies enriched transcription factor binding motifs in different categories of regulatory elements, revealing how the transcriptional regulatory code differs between tissues and species. By comparing motifs between enhancers and promoters across species and tissues, we can understand how the language of transcription factors has evolved.

---

## Requirements

- [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing)
    - Follow [these](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md) directions exactly
- [BEDTools](https://bedtools.readthedocs.io/)
- [MEME Suite](https://meme-suite.org/)
- [ChIPseeker (R package)](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
- [Anaconda3](https://www.anaconda.com/)
- HPC with SLURM

---

## Installation

Install dependencies using conda. Closely follow installation instructions for installing each required tool. 

---

## Input Files

- **ATAC-seq peak files:**
    - 4 total: human liver, human pancreas, mouse liver, mouse pancreas
    - `idr.optimal_peak.narrowPeak.gz`
    - Pipeline automatically decrompresses 
- **Cactus (multi-species) alignment:**
    - `alignment.hal`
- **Genome annotations:**
    - 2 total: human, mouse
    - `annotation.gff3.gz`
    - Pipeline automatically decrompresses 
- **Reference genomes:**
    - 2 total: human, mouse
    - `ref_genome.fa`
- **Motif database:**
    - 2 total: human, mouse
    - `motif_db.meme`

---

## Usage

First, makes ure to edit the `CONFIGURATION` section of `run_full_pipeline.sh` based on 
your file paths. 

Ensure you install HALPER exactly based on [this](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md) guide. If the guide is not 
followed exactly, this pipeline will not work. 

For the `HALPER CONFIG` step, make sure to alter the following accordingly in the 
configuration section of `run_full_pipeline.sh`:

```
export PATH=[repos dir]/hal/bin:${PATH}
export PYTHONPATH=[repos dir]/halLiftover-postprocessing:${PYTHONPATH}
```

Next, add all of the full paths for the GTF files, peak files, genome files, motif databases, and multi-species alignment. 

Finally, ensure the base directory and halper directory are modified accordingly. 

**All intermediate and final outputs will be organized in the `results/` directory.**

---

## Outputs

- `results/`: Main output directory (set as `$OUTPUT_DIR`)
    - `mapped_peaks/`: Contains HALPER-mapped peak files between species
    - `cross_species/`: Results of cross-species conservation analysis
        - `processed_files/`: BED files processed from HALPER outputs for downstream intersection; CHANGE
    - `tissue_comparison/`: Results of within-species tissue comparisons
    - `sequences/`: FASTA files for motif analysis
    - `meme_results`: Output from MEME-ChIP motif discovery
- `promoters_enhancers/`: Contains classified promoter/enhancer BEDs for each species and tissue (set as `$CLASSIFIED_DIR`)
- `logs/`: SLURM and pipeline run logs

---

## Pitfalls and Limitations

- **Memory limits:** Ensure your SLURM job requests do not exceed 2GB/core on RM-shared.
- **File paths:** Update all input/output paths as needed for your environment.
- **Genome builds:** All files must use the same genome build (e.g., hg38 for human).
- **HALPER and MEME Suite:** Ensure correct Python and software versions for compatibility.

---
## Cite Us

If you use this pipeline in your research, please cite:

Andres Chousal, Zahin Peerzade, Siddharth Sabata, Yinuo Yang. *RegulatoryElementAnalysisPipeline: A pipeline for Cross-Species Regulatory Element Analysis of Human and Mouse ATAC data*. GitHub repository. Available at: https://github.com/achousal/Bioinformatics_03713 (Accessed 29 Apr 2025).


---

## References

- HALPER: Xiaoyu Zhang, Irene Kaplow, Morgan Wirthlin, Tyler Park, Andreas Pfenning. HALPER facilitates the identification of regulatory element orthologs across species. *Bioinformatics*, Volume 36, Issue 15, 1 August 2020, Pages 4339-4340. 
- MEME Suite: Bailey TL, Johnson J, Grant CE, Noble WS. The MEME Suite. *Nucleic Acids Res*. 2015;43(W1):W39-W49. doi:10.1093/nar/gkv416 
- ChIPseeker: Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. *Bioinformatics*. 2015;31(14):2382-2383. doi:10.1093/bioinformatics/btv145 
- BEDTools: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26(6):841-842. doi:10.1093/bioinformatics/btq033 

---

*GenAI assistance used to create README*
