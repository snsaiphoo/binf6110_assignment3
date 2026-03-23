library(biomformat)
library(vegan)
library(phyloseq)
library(ggplot2)
library(ANCOMBC)
library(microbiome)

# Load BIOM
biom <- read_biom("../table_full_0.15.biom")

# Extract OTU table
otu <- as.matrix(biom_data(biom))

# Filter taxa (rows = taxa)
otu <- otu[rowSums(otu) > 10, ]

# DO NOT transpose here
# biom already has taxa as rows and samples as columns
groups <- ifelse(grepl("vegan", colnames(otu)), "vegan", "omnivore")
colors <- ifelse(groups == "vegan", "forestgreen", "orange")

# Rarefaction 
rarecurve(
  t(otu),
  step = 100000,
  col = colors,
  label = FALSE
)

legend(
  "bottomright",
  legend = c("Vegan", "Omnivore"),
  col = c("forestgreen", "orange"),
  lty = 1,
  cex = 0.8
)

# Create phyloseq OTU table
otu_ps <- otu_table(otu, taxa_are_rows = TRUE)

tax <- observation_metadata(biom)
tax <- as.data.frame(tax)

# Clean prefixes
tax <- apply(tax, 2, function(x) gsub("^[a-z]__", "", x))
tax <- as.data.frame(tax)

# Align with OTU taxa
tax <- tax[rownames(otu), ]

tax_ps <- tax_table(as.matrix(tax))

# Metadata
sample_names_vec <- colnames(otu)   # ← correct source

metadata <- data.frame(
  Sample = sample_names_vec,
  Diet = ifelse(grepl("omnivore", sample_names_vec, ignore.case = TRUE),
                "Omnivore", "Vegan")
)

rownames(metadata) <- metadata$Sample
metadata$Sample <- NULL

sample_ps <- sample_data(metadata)

# Build phyloseq object

physeq <- phyloseq(otu_ps, tax_ps, sample_ps)
colnames(tax_table(physeq)) <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

# convert to relative abundance
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
physeq_phy <- tax_glom(physeq_rel, taxrank = "Phylum")
df <- psmelt(physeq_phy)

# Get top 10 phyla by total abundance
top_phyla <- names(sort(tapply(df$Abundance, df$Phylum, sum), decreasing = TRUE))[1:10]

# Replace others with "Other"
df$Phylum <- ifelse(df$Phylum %in% top_phyla, df$Phylum, "Other")

# plotting relative abundance
p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(
    title = "Taxonomic Composition by Diet Group",
    x = "Sample",
    y = "Relative Abundance"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave(
#   "../figures/taxonomic_abundance.png",
#   plot = p,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

# plot alpha diversity 
alpha_plot <- plot_richness(
  physeq,
  x = "Diet",
  color = "Diet",
  measures = c("Observed", "Shannon")
) +
  geom_boxplot(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Alpha Diversity by Diet",
    x = "Diet",
    y = "Diversity"
  ) 

# ggsave(
#   filename = "../figures/alpha_diversity.png",
#   plot = alpha_plot,
#   width = 8,
#   height = 5,
#   dpi = 300
# )

# beta diversity

# bray-curtis - uses abundance
bray_dist <- phyloseq::distance(physeq, method = "bray")

ord <- ordinate(physeq, method = "PCoA", distance = bray_dist)

beta_df <- plot_ordination(physeq, ord, justDF = TRUE)
beta_df$Sample <- rownames(beta_df)
meta_df <- as(sample_data(physeq), "data.frame")
meta_df$Sample <- rownames(meta_df)
beta_df <- merge(beta_df, meta_df, by = "Sample")
head(beta_df)

beta_plot <- ggplot(beta_df, aes(x = Axis.1, y = Axis.2, color = Diet)) +
  geom_point(size = 4) +
  theme_minimal() + 
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCoA of Bray-Curtis Dissimilarity",
    x = "PCoA 1",
    y = "PCoA 2"
  )

beta_plot

# ggsave(
#   "../figures/beta_diversity_pcoa.png",
#   plot = beta_plot,
#   width = 7,
#   height = 5,
#   dpi = 300
# )

# jaccard - cares about presence and absence

jaccard_dist <- phyloseq::distance(physeq, method = "jaccard")
ord_jaccard <- ordinate(physeq, method = "PCoA", distance = jaccard_dist)
jaccard_df <- plot_ordination(physeq, ord_jaccard, justDF = TRUE)

jaccard_df$Sample <- rownames(jaccard_df)
meta_df <- as(sample_data(physeq), "data.frame")
meta_df$Sample <- rownames(meta_df)
jaccard_df <- merge(jaccard_df, meta_df, by = "Sample")

jaccard_plot <- ggplot(jaccard_df, aes(x = Axis.1, y = Axis.2, color = Diet)) +
  geom_point(size = 4) +
  theme_minimal() +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCoA of Jaccard Distance",
    x = "PCoA 1",
    y = "PCoA 2"
  )

jaccard_plot
# 
# ggsave("../figures/jaccard_pcoa.png", plot = jaccard_plot, width = 7, height = 5, dpi = 300)

# nmds

ord_nmds <- ordinate(physeq, method = "NMDS", distance = "bray")

nmds_df <- plot_ordination(physeq, ord_nmds, justDF = TRUE)

nmds_df$Sample <- rownames(nmds_df)

nmds_df <- merge(nmds_df, meta_df, by = "Sample")

nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Diet)) +
  geom_point(size = 4) +
  theme_minimal() + 
  stat_ellipse(level = 0.95) +
  labs(
    title = "NMDS (Bray-Curtis)",
    x = "NMDS1",
    y = "NMDS2"
  )

nmds_plot

# ggsave("../figures/nmds_bray.png", plot = nmds_plot, width = 7, height = 5, dpi = 300)

# ANCOM-BC2 - Differential Abundance
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ANCOMBC)

# 1. Clean taxonomy
tax <- as.data.frame(tax_table(physeq))
tax[is.na(tax)] <- "Unknown"
tax[tax == ""] <- "Unknown"
tax_table(physeq) <- tax_table(as.matrix(tax))

# 2. Aggregate to Family level 
physeq_family <- tax_glom(physeq, "Family")

physeq_family <- prune_taxa(
  taxa_sums(physeq_family) > 10,   # stricter than before
  physeq_family
)

# 3. Run ANCOM-BC2
ancombc.out <- ancombc2(
  data = physeq_family,
  tax_level = "Family",
  fix_formula = "Diet",
  rand_formula = NULL,
  group = "Diet",
  p_adj_method = "BH",
  pseudo_sens = TRUE,   
  prv_cut = 0.20,
  lib_cut = 1000,
  s0_perc = 0.05,
  struc_zero = TRUE,
  neg_lb = TRUE
)

# 4. Extract results
res <- ancombc.out$res

# Significant taxa
sig <- subset(res, q_DietVegan < 0.05)

# View results
head(res)
sig

top_res <- res %>%
  arrange(q_DietVegan) %>%
  head(20)

ggplot(top_res, aes(x = lfc_DietVegan, y = reorder(taxon, lfc_DietVegan))) +
  geom_point(aes(color = p_DietVegan < 0.05), size = 3) +
  geom_vline(xintercept = 0, color = "red") +
  labs(
    title = "Top Differentially Abundant Families",
    x = "Log Fold Change (Vegan vs Omnivore)",
    y = "Family"
  ) +
  theme_minimal()







