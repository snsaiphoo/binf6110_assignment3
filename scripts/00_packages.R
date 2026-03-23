# Install and load required packages
# Do this before running 10_metanalysis.R
# Versions are
# [1] dplyr_1.2.0       microbiome_1.32.0 ANCOMBC_2.12.1    ggplot2_4.0.2     phyloseq_1.54.2   vegan_2.7-3      
# [7] permute_0.9-10    biomformat_1.38.3


# Install CRAN packages
install.packages(c("ggplot2", "dplyr", "vegan"))

# Install Bioconductor manager
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("phyloseq", "ANCOMBC", "microbiome", "biomformat"))

# Load libraries
library(biomformat)
library(vegan)
library(phyloseq)
library(ggplot2)
library(ANCOMBC)
library(microbiome)
library(dplyr)