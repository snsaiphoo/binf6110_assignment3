library(dplyr)

files <- list.files("bracken_05", pattern = "*.bracken", full.names = TRUE)
s
data_list <- lapply(files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  df <- df[, c("name", "new_est_reads")]
  
  # Clean sample name
  sample_name <- gsub(".bracken", "", basename(f))
  colnames(df)[2] <- sample_name
  
  return(df)
})

merged <- Reduce(function(x, y) full_join(x, y, by = "name"), data_list)

merged[is.na(merged)] <- 0

rownames(merged) <- merged$name
merged$name <- NULL

dim(merged)
colnames(merged)


sample_names <- colnames(merged)

metadata <- data.frame(
  Sample = sample_names,
  Diet = ifelse(grepl("omnivore", sample_names), "Omnivore", "Vegan")
)

rownames(metadata) <- metadata$Sample
metadata$Sample <- NULL


library(phyloseq)

otu <- otu_table(as.matrix(merged), taxa_are_rows = TRUE)
sample_data <- sample_data(metadata)

physeq <- phyloseq(otu, sample_data)

physeq

# Rarefaction
library(vegan)
otu_table <- as.data.frame(t(otu_table(physeq)))
rare_curve <- rarecurve(otu_table, step = 1000)


plot_richness(physeq, x = "Diet", measures = c("Shannon", "Simpson"))