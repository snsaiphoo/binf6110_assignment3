# Install and load required packages
# Do this before running 10_metanalysis.R

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