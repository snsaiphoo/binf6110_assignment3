library(biomformat)

# Read BIOM
biom <- read_biom("bracken_0.15_full.biom")

# Extract OTU table
otu <- as.matrix(biom_data(biom))

# Transpose so samples = rows (for vegan)
otu <- t(otu)

library(vegan)

# Filter low-abundance taxa (THIS is why yours was slow)
otu <- otu[, colSums(otu) > 10]

# Rarefaction (faster + stable)
rarecurve(
  otu,
  step = 10000,
  sample = min(rowSums(otu)),
  label = TRUE
)


# Build phyloseq properly
library(phyloseq)

# Create metadata (you were missing this)
sample_names <- rownames(otu)

metadata <- data.frame(
  Sample = sample_names,
  Diet = ifelse(grepl("omnivore", sample_names), "Omnivore", "Vegan")
)

rownames(metadata) <- metadata$Sample
metadata$Sample <- NULL

# IMPORTANT: phyloseq expects taxa as rows
otu_ps <- otu_table(t(otu), taxa_are_rows = TRUE)

sample_data_ps <- sample_data(metadata)

physeq <- phyloseq(otu_ps, sample_data_ps)

physeq
