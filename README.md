# Shotgun Metagenomic Analysis of Gut Microbiomes in Vegans and Omnivores

## Introduction 
In recent years, the human gut microbiome has been a growing point of conversation due to the high genomic diversity [1]. The current analysis will be based on research by De Filippis et al., who investigated whether specific _Prevotella copri_ strains were associated with specific diets [1]. Their results found that those with high-fiber vegan diets had specific strains present that play a role in carbohydrate catabolism [1]. Whereas subjects on an omnivorous diet showed greater expression of the leuB gene, which contributes to the pathway of branched-chain amino acid biosynthesis [1]. Overall, these initial findings suggest a correlation between dietary patterns and the composition of the gut microbiome, underscoring the importance of this field of study. Specifically, this comparison is important because of the potential implications when evaluating Western and non-Western dietary patterns. These populations typically differ in macronutrient intake: Western diets are generally higher in protein and fat, while non-Western diets are often higher in dietary fiber. To further evaluate these differences, shotgun metagenomic analysis will be used to examine microbial diversity within the gut microbiomes of three vegans and three omnivores. This analysis will include measures of alpha diversity, beta diversity, and differential abundance to compare microbial communities between the two dietary groups.

### Design Rationale 
For this analysis, gut microbiome data will be taken from the study by De Filippis et al., using datasets from three vegans and three omnivores. Each sample contains DNA from millions of microorganisms, and this analysis will focus on identifying the presence of _Prevotella copri_ strains across the different dietary groups [1]. The sequencing data will be obtained from the NCBI Sequence Read Archive (SRA), and the `SRA Toolkit` will be used in the Linux terminal to download and process the raw sequencing files.

Before analysis, the sequencing reads will first be checked for quality using `FastQC` [2]. This program provides information about sequencing quality, GC content, and possible adapter contamination. The `FastQC` results will then be aggregated and summarized using `MultiQC` [3], allowing for efficient comparison across all samples. These results help determine whether trimming is needed and guide how the trimming should be performed. If necessary, reads will be trimmed using `fastp` [4] to remove low-quality bases and adapter sequences. Host DNA contamination has already been performed on the SRR NCBI datasets for all samples.

Taxonomic classification will then be used to identify which microorganisms are present in the metagenomic samples. Several types of tools can be used for this step, including DNA-based classifiers such as `Kraken2` and `Centrifuge`, protein-based classifiers such as `Kaiju`, and marker-based methods such as `MetaPhlAn2`. In the original study used for this analysis, De Filippis et al. performed taxonomic profiling using the marker-based tool `MetaPhlAn2` [1]. However, benchmarking studies have shown that DNA-based classifiers generally perform well for whole-genome metagenomic datasets [5]. Because tool performance can vary depending on factors such as database choice, computational resources, and dataset type, careful tool selection is important when designing a metagenomic analysis pipeline [6].

For this analysis, `Kraken2` was selected because it provides fast classification and performs well with large shotgun sequencing datasets. One limitation of `Kraken2` is that some reads cannot be assigned unambiguously at the species level because closely related organisms may share nearly identical DNA sequences. In these cases, reads are assigned to higher taxonomic levels such as genus or family. To address this, `Bracken` will be used after `Kraken2` to improve abundance estimation. `Bracken` redistributes reads assigned to higher taxonomic levels to more specific taxa, allowing more accurate species- and genus-level abundance estimates from metagenomic datasets [7].

Following taxonomic classification and abundance estimation, microbial community diversity will be analyzed to compare the gut microbiomes of the vegan and omnivore groups. These analyses will be performed in `R` using the `phyloseq` package, which provides a framework for organizing and analyzing microbiome data and integrates with commonly used packages such as `vegan` and `ggplot2` for ecological analysis and visualization [8]. Alpha diversity metrics will be calculated to measure diversity within each sample, while beta diversity will be used to compare microbial community composition between groups. Differences between samples will be evaluated using Bray–Curtis dissimilarity, which incorporates taxonomic abundance information and has been shown to be a sensitive metric for detecting differences between microbial communities [9].

The last step is the differential abundance analysis, which will identify taxa that differ significantly between the vegan and omnivore samples. Differential abundance testing is a central component of microbiome analysis. However, in microbiome datasets, only relative abundances of taxa are observed, which makes it challenging to identify truly differential taxa [10]. Due to this, `ANCOM-BC` (Analysis of Compositions of Microbiomes with Bias Correction) will be used for differential abundance, which was specifically developed to account for compositional bias in microbiome sequencing data. Methods originally developed for RNA-seq analysis, such as `DESeq2` and `edgeR`, have been shown to perform poorly on microbiome taxonomic profiles, except for limma, due to differences in the statistical properties of the data [10]. In contrast, methods such as `ANCOM` and `ANCOM-BC` better control false discovery rates when identifying differential taxa in microbiome datasets [11]. Therefore, `ANCOM-BC` was selected to identify taxa whose abundance differs significantly between the vegan and omnivore groups in this study.

## Methods
The metagenomics analysis was conducted using the _Digital Research Alliance of Canada Nibi Cluster_ and `RStudio`. All software and modules on the cluster were run through `Docker` container images executed with `Apptainer`, ensuring a consistent and reproducible environment. The workflow is outlined below, with all associated scripts available in the [`scripts`](scripts/) directory. All .sh scripts were executed on the Nibi scratch environment, while downstream analysis in `RStudio` was performed outside of the cluster.

### 1.0 - Data Acquisition & Software Setup 
#### 1.1 - Nibi Cluster Containers
To ensure reproducibility and maintain consistent software versions, all command-line tools were executed within containerized environments. `Apptainer` was used to retrieve pre-built `Docker` images and convert them into `.sif` container files, which were organized within a dedicated `containers/` directory. The process of building these containers was automated using the [`00_buildcontainers.sh`](scripts/00_buildcontainers.sh) script, which includes version numbers.

#### 1.2 - R Environment Setup 
Data analysis was performed in `R` version 4.5.1. All necessary `CRAN` and `Bioconductor` dependencies, including their respective versions, can be installed using the [`00_packages.R`](scripts/00_packages.R) script included in this repository.

#### 1.3 - Data Acquisition 
The data used in this analysis were obtained from the NCBI Sequence Read Archive (SRA) and consisted of three vegan and three omnivore samples from the study by De Filippis et al. [1]. 

* **Vegan Samples** - SRR8146944, SRR8146951, and SRR8146954 
* **Omnivore Samples** - SRR8146935, SRR8146936, and SRR8146938 
  
Raw SRA files were obtained using the `prefetch` command from the `SRA Toolkit` (v3.2.1) container, as outlined in [`01_data.sh`](scripts/01_data.sh). The downloaded samples were then validated to ensure complete retrieval in [`02_data.sh`](scripts/02_data.sh). Conversion to `FASTQ` format was subsequently performed using `fasterq-dump` in [`03_data.sh`](scripts/03_data.sh). The resulting FASTQ files were then compressed using `pigz`, enabling parallel compression to reduce computational time. 

### 2.0 - Quality Control & Trimming 
#### 2.1 - Quality Control with FastQC

#### 2.2 - Quality Control with MultiQC

#### 2.3 - Trimming with fastp

### 3.0 - Taxonomic Classification with Kraken2 
#### 3.1 - Standard 16GB Database, 0.05 * 0.10
#### 3.2 - Standard Full Database, 0.15

### 4.0 - Bracken for Normalization 
#### 4.1 - Bracken BIOM table 


## References
[1] F. De Filippis et al., “Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets,” Cell Host & Microbe, vol. 25, no. 3, pp. 444-453.e3, Mar. 2019, doi: https://doi.org/10.1016/j.chom.2019.01.004. <br/>
[2] s-andrews, “s-andrews/FastQC,” GitHub, Nov. 20, 2018. https://github.com/s-andrews/FastQC <br/>
[3]“MultiQC,” Introduction to RNA-seq using high-performance computing, Nov. 19, 2021. https://hbctraining.github.io/Intro-to-rnaseq-fasrc-salmon-flipped/lessons/11_multiQC.html (accessed Mar. 17, 2026). <br/>
[4] S. Chen, “fastp 1.0: An ultra‐fast all‐round tool for FASTQ data quality control and preprocessing,” iMeta, vol. 4, no. 5, Sep. 2025, doi: https://doi.org/10.1002/imt2.70078.  <br/>
[5] L. C. Terrón-Camero, F. Gordillo-González, E. Salas-Espejo, and E. Andrés-León, “Comparison of Metagenomics and Metatranscriptomics Tools: A Guide to Making the Right Choice,” Genes, vol. 13, no. 12, p. 2280, Dec. 2022, doi: https://doi.org/10.3390/genes13122280. <br/>
[6] I. B. Martins, J. M. Silva, and J. R. Almeida, “A systematic review and benchmarking of modern metagenomic tools for taxonomic classification,” Computers in Biology and Medicine, vol. 206, p. 111600, Apr. 2026, doi: https://doi.org/10.1016/j.compbiomed.2026.111600. <br/>
[7] J. Lu, F. P. Breitwieser, P. Thielen, and S. L. Salzberg, “Bracken: estimating species abundance in metagenomics data,” PeerJ Computer Science, vol. 3, p. e104, Jan. 2017, doi: https://doi.org/10.7717/peerj-cs.104. <br/>
[8] T. Wen, G. Niu, T. Chen, Q. Shen, J. Yuan, and Y. Liu, “The best practice for microbiome analysis using R,” Protein & Cell, vol. 14, no. 10, pp. 713–725, May 2023, doi: https://doi.org/10.1093/procel/pwad024. <br/>
[9] J. G. Kers and E. Saccenti, “The Power of Microbiome Studies: Some Considerations on Which Alpha and Beta Metrics to Use and How to Report Results,” Frontiers in Microbiology, vol. 12, p. 796025, Mar. 2022, doi: https://doi.org/10.3389/fmicb.2021.796025.  <br/>
[10] H. Zhou, K. He, J. Chen, and X. Zhang, “LinDA: linear models for differential abundance analysis of microbiome compositional data,” Genome Biology, vol. 23, no. 1, Apr. 2022, doi: https://doi.org/10.1186/s13059-022-02655-5. <br/>
[11] J. Wirbel, M. Essex, S. K. Forslund, and G. Zeller, “A realistic benchmark for differential abundance testing and confounder adjustment in human microbiome studies,” Genome Biology, vol. 25, no. 1, Sep. 2024, doi: https://doi.org/10.1186/s13059-024-03390-9. <br/>

