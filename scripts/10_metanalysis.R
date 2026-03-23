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

# Assign colors for plotting
colors <- ifelse(groups == "vegan", "forestgreen", "orange")

# Rarefaction curve to check sequencing depth across samples
# Helps see if samples are comparable
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

# Create phyloseq object for downstream analysis
# OTU table
otu_ps <- otu_table(otu, taxa_are_rows = TRUE)

# Extract taxonomy
tax <- observation_metadata(biom)
tax <- as.data.frame(tax)

# Remove prefixes like k__, p__, etc. for cleaner labels
tax <- apply(tax, 2, function(x) gsub("^[a-z]__", "", x))
tax <- as.data.frame(tax)

# Match taxonomy rows to OTU table
tax <- tax[rownames(otu), ]
tax_ps <- tax_table(as.matrix(tax))

# Creating the metadata
sample_names_vec <- colnames(otu)

metadata <- data.frame(
  Sample = sample_names_vec,
  Diet = ifelse(grepl("omnivore", sample_names_vec, ignore.case = TRUE),
                "Omnivore", "Vegan")
)

# Set sample names properly
rownames(metadata) <- metadata$Sample
metadata$Sample <- NULL

sample_ps <- sample_data(metadata)

# Combine all into one phyloseq object
physeq <- phyloseq(otu_ps, tax_ps, sample_ps)

# Changing names from taxonomy1 to the correct labels
colnames(tax_table(physeq)) <- c(
  "Kingdom","Phylum","Class","Order","Family","Genus","Species"
)

# Relative abundance plot at family level

# Convert counts to relative abundance (proportions)
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Aggregate taxa at family level
physeq_fam <- tax_glom(physeq_rel, taxrank = "Family")

# Convert to dataframe for plotting
df <- psmelt(physeq_fam)

# Keep top 10 most abundant family, group rest as "Other"
top_fam <- names(sort(tapply(df$Abundance, df$Family, sum), decreasing = TRUE))[1:10]
df$Family <- ifelse(df$Family %in% top_fam, df$Family, "Other")

# Clean sample names for readability
df$Sample <- sub("^(SRR[0-9]+).*", "\\1", df$Sample)

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Family Level Composition by Diet Group",
       x = "Sample", y = "Relative Abundance") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p

ggsave(
  "../figures/taxonomic_abundance_family.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

# Relative abundance plot at genus level

# Convert counts to relative abundance (proportions)
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Aggregate taxa at family level
physeq_gen <- tax_glom(physeq_rel, taxrank = "Genus")

# Convert to dataframe for plotting
df <- psmelt(physeq_gen)

# Keep top 10 most abundant family, group rest as "Other"
top_gen <- names(sort(tapply(df$Abundance, df$Genus, sum), decreasing = TRUE))[1:10]
df$Genus <- ifelse(df$Genus %in% top_gen, df$Genus, "Other")

# Clean sample names for readability
df$Sample <- sub("^(SRR[0-9]+).*", "\\1", df$Sample)

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Genus Level Composition by Diet Group",
       x = "Sample", y = "Relative Abundance") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p

ggsave(
  "../figures/taxonomic_abundance_genus.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

# Relative abundance at species level (more detailed view)
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# aggregate at species level
physeq_species <- tax_glom(physeq_rel, taxrank = "Species")

df <- psmelt(physeq_species)

# Keep top 10 species, group others
top_species <- names(sort(tapply(df$Abundance, df$Species, sum), decreasing = TRUE))[1:10]

df$Species <- ifelse(df$Species %in% top_species, df$Species, "Other")

df$Sample <- sub("^(SRR[0-9]+).*", "\\1", df$Sample)

# Plot species-level composition
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

# Alpha diversity analysis
# Plot Observed richness and Shannon diversity by diet
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

# Statistical test to compare diversity between diets
# Wilcoxon test used since data may not be normally distributed

alpha_df <- estimate_richness(physeq, measures = c("Observed", "Shannon"))
alpha_df$Diet <- sample_data(physeq)$Diet

wilcox.test(Observed ~ Diet, data = alpha_df)
wilcox.test(Shannon ~ Diet, data = alpha_df)

# Beta diversity and validation 

# Calculate Bray-Curtis distance (measures differences between samples)
bray_dist <- phyloseq::distance(physeq, method = "bray")

# PERMANOVA test to check if microbiome composition differs by diet
# Tests whether vegan vs omnivore groups are significantly different
set.seed(123)
permanova <- adonis2(bray_dist ~ Diet, data = metadata)
print(permanova)

# PCoA ordination to visualize differences between samples
# Reduces complex distance data into 2 dimensions

ord <- ordinate(physeq, method = "PCoA", distance = bray_dist)

# Convert ordination results to dataframe for plotting
beta_df <- plot_ordination(physeq, ord, justDF = TRUE)
beta_df$Sample <- rownames(beta_df)

# Extract metadata and match with ordination results
meta_df <- as(sample_data(physeq), "data.frame")
meta_df$Sample <- rownames(meta_df)

# Merge ordination data with metadata (to get Diet labels)
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
# Clean taxonomy to avoid NA or empty values
# Replace missing taxonomy with "Unknown" so analysis does not break

tax <- as.data.frame(tax_table(physeq))
tax[is.na(tax)] <- "Unknown"
tax[tax == ""] <- "Unknown"
tax_table(physeq) <- tax_table(as.matrix(tax))

# Differential abundance at the Family level using ANCOM-BC2
# Aggregate taxa at family level
physeq_family <- tax_glom(physeq, "Family")

# Remove very low abundance taxa to reduce noise
physeq_family <- prune_taxa(
  taxa_sums(physeq_family) > 50,
  physeq_family
)

# Keep only bacteria 
physeq_family <- subset_taxa(physeq_family, Kingdom == "Bacteria")

# Run ANCOM-BC2 using BH correction
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

# Select top 15 families with strongest effect
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
# Differential abundance at the Genus level

physeq_genus <- tax_glom(physeq, "Genus")

physeq_genus <- prune_taxa(
  taxa_sums(physeq_genus) > 20,
  physeq_genus
)

physeq_genus <- subset_taxa(physeq_genus, Kingdom == "Bacteria")

# Run ANCOM-BC2 on Genus (BH correction)
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

# Get significant genera, empty
sig_genus <- subset(res_genus, q_DietVegan < 0.05)

# Check if any genera are significant
nrow(sig_genus)

head(res_genus)

# Since no genera are significant, show top ones by effect size instead
top_genus <- res_genus %>%
  mutate(abs_lfc = abs(lfc_DietVegan)) %>%
  arrange(desc(abs_lfc)) %>%
  head(20)

# Plot top genera by effect size (not statistically significant)
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
