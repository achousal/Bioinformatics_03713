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
- [Output Files](#output-files)
- [Parameters](#parameters)
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
Open chromatin regions identified from ATAC-seq are mapped between human and mouse genomes using HALPER (HAL Liftover Post-processing for Epigenomic Regions) with Cactus whole-genome alignments. This step is crucial for identifying orthologous regulatory elements between species, allowing direct comparison of regulatory activity at corresponding genomic locations. The mapping process accounts for genomic rearrangements and evolutionary changes that have occurred since the divergence of humans and mice. ***add something about IDR optimal narrow peaks?***
3. **Conservation and Specificity Analysis** (Question 1)
This step identifies and categorizes regulatory elements based on their conservation patterns: elements conserved between species for the same tissue, elements shared across tissues within a species, and elements specific to a tissue or species. By quantifying these different categories, we can directly address whether regulatory element activity is more conserved across tissues or species. This analysis provides insights into the evolutionary constraints on gene regulation.
4. **Functional Annotation** (Question 3)
Regulatory elements identified in previous steps are annotated using ChIPseeker and analyzed for GO term enrichment with GREAT. This analysis connects regulatory elements to their potential target genes and biological functions, revealing which biological processes are regulated by conserved versus species-specific elements. The functional analysis helps determine whether similar biological processes are regulated by conserved elements across species, despite potential differences in the specific regulatory elements.
5. **Enhancer/Promoter Classification** (Question 2)
Open chromatin regions are classified as promoters (within 2kb upstream and 200bp downstream of TSS) or enhancers (all other regions) using BEDTools and genome annotations. This classification is essential for understanding how conservation patterns differ between these two types of regulatory elements. Promoters, which are proximal to genes, may be under different evolutionary constraints than enhancers, which can act at a distance and may evolve more rapidly.
6. **Motif Discovery** (Question 2)
Sequences from classified regulatory elements are extracted and analyzed using MEME-ChIP for de novo motif discovery. This step identifies enriched transcription factor binding motifs in different categories of regulatory elements, revealing how the transcriptional regulatory code differs between tissues and species. By comparing motifs between enhancers and promoters across species and tissues, we can understand how the language of transcription factors has evolved.

---

## Requirements

- [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing)
- [BEDTools](https://bedtools.readthedocs.io/)
- [MEME Suite](https://meme-suite.org/)
- [GREAT](http://great.stanford.edu/)
- [ChIPseeker (R package)](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
- [Anaconda3](https://www.anaconda.com/)
- HPC with SLURM

---

## Installation

Install dependencies using conda. Closely follow installation instructions for installing each required tool. 

---

## Input Files

- **ATAC-seq peak files:**
    - Human: `/ocean/projects/bio230007p/ikaplow/HumanAtac/[tissue]/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`
    - Mouse: `/ocean/projects/bio230007p/ikaplow/MouseAtac/[tissue]/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`
- **Cactus alignment:**
    - `/ocean/projects/bio230007p/ikaplow/Alignments/12-mammals.hal`
- **Genome annotations:**
    - Human: `/ocean/projects/bio230007p/ssabata/HumanGenomeInfo/gencode.v47.annotation.gff3.gz`
    - Mouse: `/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM25.annotation.gtf.gz`
- **Reference genomes:**
    - Human: `/ocean/projects/bio230007p/ssabata/HumanGenomeInfo/hg38.fa`
    - Mouse: `/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/GRCm38.primary_assembly.genome.fa`
- **Motif database:**
    - `/ocean/projects/bio230007p/ikaplow/CIS-BP_2.00`

---

## Usage

Run the entire pipeline with:

```bash
bash run_full_pipeline.sh
```

This script will:

- Perform cross-species mapping with HALPER
- Identify conserved and tissue-specific peaks
- Classify peaks as promoters/enhancers
- Extract sequences and run motif discovery
- Generate summary statistics and output files

**All intermediate and final outputs will be organized in the `results/` directory.**

---

## Output Files

- `mapped_peaks/` — Cross-species mapped peaks (BED)
- `conserved_peaks/`, `tissue_specific/`, `shared_peaks/` — Conservation analysis (BED)
- `classified/` — Promoter/enhancer classification (BED)
- `sequences/` — FASTA files for motif analysis
- `meme_results/` — Motif discovery results (MEME-ChIP)
- `functional_analysis/` — Functional annotation and enrichment results
- `statistics.txt` — Summary of peak counts and classifications

---

## Parameters

Key parameters can be set at the top of `run_full_pipeline.sh`, including:

- Input file locations
- Output directory
- Number of CPUs and memory
- Promoter region definition (default: 2kb upstream, 200bp downstream of TSS)

---

## Pitfalls and Limitations

- **Memory limits:** Ensure your SLURM job requests do not exceed 2GB/core on RM-shared.
- **File paths:** Update all input/output paths as needed for your environment.
- **Genome builds:** All files must use the same genome build (e.g., hg38 for human).
- **HALPER and MEME Suite:** Ensure correct Python and software versions for compatibility.

---

## References

- [HALPER: https://github.com/pfenninglab/halLiftover-postprocessing](https://github.com/pfenninglab/halLiftover-postprocessing)
- [MEME Suite: https://meme-suite.org/](https://meme-suite.org/)
- [GREAT: http://great.stanford.edu/](http://great.stanford.edu/)
- [ChIPseeker: https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
- [BEDTools: https://bedtools.readthedocs.io/](https://bedtools.readthedocs.io/)

---

*GenAI assistance used to create README*
