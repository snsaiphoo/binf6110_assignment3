library(biomformat)
library(vegan)
library(phyloseq)
library(ggplot2)
library(ANCOMBC)
library(microbiome)
library(dplyr)

# Load the data 
# 10 vegans and 10 omnivores
biom <- read_biom("../table_full_0.15.biom")
otu <- as.matrix(biom_data(biom))

# Filter low abundance taxa
otu <- otu[rowSums(otu) > 10, ]

# Groups for plotting
groups <- ifelse(grepl("vegan", colnames(otu)), "vegan", "omnivore")
colors <- ifelse(groups == "vegan", "forestgreen", "orange")

# Rarefraction curve for base visualization 

rarecurve(
  t(otu),
  step = 100000,
  col = colors,
  label = FALSE
)

title("Vegan vs Omnivore Samples")

legend(
  "bottomright",
  legend = c("Vegan", "Omnivore"),
  col = c("forestgreen", "orange"),
  lty = 1,
  cex = 0.8
)

# Creating the phyloseq object 

otu_ps <- otu_table(otu, taxa_are_rows = TRUE)

tax <- observation_metadata(biom)
tax <- as.data.frame(tax)

# Clean taxonomy prefixes
tax <- apply(tax, 2, function(x) gsub("^[a-z]__", "", x))
tax <- as.data.frame(tax)

tax <- tax[rownames(otu), ]
tax_ps <- tax_table(as.matrix(tax))

# Creating the metadata
sample_names_vec <- colnames(otu)

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

# changing names from taxonomy1 to the correct labels
colnames(tax_table(physeq)) <- c(
  "Kingdom","Phylum","Class","Order","Family","Genus","Species"
)

# relative abundance plot

physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
physeq_phy <- tax_glom(physeq_rel, taxrank = "Phylum")

df <- psmelt(physeq_phy)

# identify the top phyla
top_phyla <- names(sort(tapply(df$Abundance, df$Phylum, sum), decreasing = TRUE))[1:10]
df$Phylum <- ifelse(df$Phylum %in% top_phyla, df$Phylum, "Other")
df$Sample <- sub("^(SRR[0-9]+).*", "\\1", df$Sample)

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Taxonomic Composition by Diet Group",
       x = "Sample", y = "Relative Abundance") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p

# relative abundance with species
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# aggregate at species level
physeq_species <- tax_glom(physeq_rel, taxrank = "Species")

df <- psmelt(physeq_species)

# identify top species
top_species <- names(sort(tapply(df$Abundance, df$Species, sum), decreasing = TRUE))[1:10]

df$Species <- ifelse(df$Species %in% top_species, df$Species, "Other")

df$Sample <- sub("^(SRR[0-9]+).*", "\\1", df$Sample)

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Species-Level Composition by Diet Group",
       x = "Sample", y = "Relative Abundance") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p

# ggsave(
#   "../figures/taxonomic_abundance_species.png",
#   plot = p,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

# alpha diversity

alpha_plot <- plot_richness(
  physeq,
  x = "Diet",
  color = "Diet",
  measures = c("Observed", "Shannon")
) +
  geom_boxplot(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Alpha Diversity by Diet",
       x = "Diet", y = "Diversity")

alpha_plot


# ggsave(
#   filename = "../figures/alpha_diversity.png",
#   plot = alpha_plot,
#   width = 8,
#   height = 5,
#   dpi = 300
# )

# Perform Wilcoxon rank-sum tests to further assess differences in alpha diversity between diet groups

alpha_df <- estimate_richness(physeq, measures = c("Observed", "Shannon"))
alpha_df$Diet <- sample_data(physeq)$Diet

wilcox.test(Observed ~ Diet, data = alpha_df)
wilcox.test(Shannon ~ Diet, data = alpha_df)


# beta diversity and validation 

bray_dist <- phyloseq::distance(physeq, method = "bray")

# PERMANOVA
set.seed(123)
permanova <- adonis2(bray_dist ~ Diet, data = metadata)
print(permanova)

# PCoA
ord <- ordinate(physeq, method = "PCoA", distance = bray_dist)

beta_df <- plot_ordination(physeq, ord, justDF = TRUE)
beta_df$Sample <- rownames(beta_df)

meta_df <- as(sample_data(physeq), "data.frame")
meta_df$Sample <- rownames(meta_df)

beta_df <- merge(beta_df, meta_df, by = "Sample")

beta_plot <- ggplot(beta_df, aes(x = Axis.1, y = Axis.2, color = Diet)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCoA of Bray-Curtis",
       x = "PCoA 1", y = "PCoA 2")

beta_plot

# ggsave(
#   "../figures/beta_diversity_pcoa.png",
#   plot = beta_plot,
#   width = 7,
#   height = 5,
#   dpi = 300
# )

# Clean taxonomy

tax <- as.data.frame(tax_table(physeq))
tax[is.na(tax)] <- "Unknown"
tax[tax == ""] <- "Unknown"
tax_table(physeq) <- tax_table(as.matrix(tax))

# Differential abundance with ANCOMBC2 with original phyloseq object

physeq_family <- tax_glom(physeq, "Family")

physeq_family <- prune_taxa(
  taxa_sums(physeq_family) > 50,
  physeq_family
)

physeq_family <- subset_taxa(physeq_family, Kingdom == "Bacteria")

ancombc.family <- ancombc2(
  data = physeq_family,
  tax_level = "Family",
  fix_formula = "Diet",
  group = "Diet",
  p_adj_method = "BH",
  prv_cut = 0.20,
  lib_cut = 1000,
  pseudo_sens = TRUE,
  struc_zero = TRUE,
  neg_lb = TRUE
)

res_family <- ancombc.family$res

sig_family <- subset(res_family, q_DietVegan < 0.05)

library(dplyr)
library(ggplot2)

# Take top 15 strongest effects
sig_family_plot <- sig_family %>%
  arrange(desc(lfc_DietVegan)) %>%
  head(15)

# Plot
f <- ggplot(sig_family_plot,
            aes(x = lfc_DietVegan,
                y = reorder(taxon, lfc_DietVegan))) +
  geom_point(color = "steelblue", size = 3) +
  geom_vline(xintercept = 0, color = "red") +
  labs(title = "Differentially Abundant Families (ANCOM-BC2, BH)",
       x = "Log Fold Change (Vegan vs Omnivore)",
       y = "Family") +
  theme_minimal()

f

ggsave(
  "../figures/sda_families.png",
  plot = f,
  width = 7,
  height = 5,
  dpi = 300
)
# Differential Abundance at the Genus Level

physeq_genus <- tax_glom(physeq, "Genus")

physeq_genus <- prune_taxa(
  taxa_sums(physeq_genus) > 20,
  physeq_genus
)

physeq_genus <- subset_taxa(physeq_genus, Kingdom == "Bacteria")

ancombc.genus <- ancombc2(
  data = physeq_genus,
  tax_level = "Genus",
  fix_formula = "Diet",
  group = "Diet",
  p_adj_method = "BH",
  prv_cut = 0.20,
  lib_cut = 1000,
  pseudo_sens = TRUE,
  struc_zero = TRUE,
  neg_lb = TRUE
)

res_genus <- ancombc.genus$res


sig_genus <- subset(res_genus, q_DietVegan < 0.05)
nrow(sig_genus)

head(res_genus)
top_genus <- res_genus %>%
  mutate(abs_lfc = abs(lfc_DietVegan)) %>%
  arrange(desc(abs_lfc)) %>%
  head(20)

g <- ggplot(top_genus, aes(x = reorder(taxon, lfc_DietVegan), y = lfc_DietVegan)) +
  geom_point(aes(color = lfc_DietVegan > 0), size = 3) +
  geom_hline(yintercept = 0, color = "black") +
  coord_flip() +
  scale_color_manual(values = c("red", "blue"),
                     labels = c("Omnivore-enriched", "Vegan-enriched"),
                     name = "Direction") +
  labs(
    title = "Top Genera by Effect Size (ANCOM-BC2)",
    subtitle = "No statistically significant Genera",
    x = "Genus",
    y = "Log Fold Change"
  ) +
  theme_minimal()

g

ggsave(
  "../figures/top_genera_effect_size.png",
  plot = g,
  width = 7,
  height = 5,
  dpi = 300
)
