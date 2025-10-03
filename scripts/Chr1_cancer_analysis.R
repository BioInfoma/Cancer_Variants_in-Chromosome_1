# 📦 1. Load Libraries
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("vcd")) install.packages("vcd")

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(vcd)

# 📁 2. Load PCA Data
pca <- fread("pca_chr1.eigenvec")
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

# 🧬 3. Load Metadata (Download separately if not available)
meta <- fread("integrated_call_samples_v3.20130502.ALL.panel", skip = 1,
              col.names = c("IID", "Population", "Superpopulation", "Gender"))

# 🧠 4. Merge PCA + Metadata
merged <- inner_join(pca, meta, by = "IID")

# 📈 5. PCA Plot
ggplot(merged, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA of Chromosome 1 Cancer Variants", x = "PC1", y = "PC2")

# 🍩 6. Donut Chart
pop_counts <- meta %>% count(Superpopulation)
ggplot(pop_counts, aes(x = "", y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Sample Distribution by Superpopulation")

# 📜 7. Load Genotypes
geno <- fread("genotypes.tsv", header = FALSE)
sample_ids <- fread("sample_ids.txt", header = FALSE)[[1]]
colnames(geno)[1:4] <- c("CHR", "POS", "REF", "ALT")
colnames(geno)[5:(4 + length(sample_ids))] <- sample_ids
geno$ID <- paste0(geno$CHR, "_", geno$POS, "_", geno$REF, "_", geno$ALT)

# 🔁 8. Reshape Genotype Table
geno_long <- melt(geno, id.vars = c("CHR", "POS", "REF", "ALT", "ID"),
                  variable.name = "Sample", value.name = "GT")

# 🧬 9. Join with Metadata
geno_meta <- inner_join(geno_long, meta, by = c("Sample" = "IID"))

# 🧪 10. Classify Genotypes
geno_meta <- geno_meta %>%
  mutate(Genotype = case_when(
    GT %in% c("0/0", "0|0") ~ "REF",
    GT %in% c("0/1", "1/0", "0|1", "1|0") ~ "HET",
    GT %in% c("1/1", "1|1") ~ "ALT",
    TRUE ~ NA_character_
  ))

# 📊 11. Barplot of Genotype Distribution
counts <- geno_meta %>%
  filter(!is.na(Genotype)) %>%
  count(ID, Superpopulation, Genotype) %>%
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0)

counts_long <- counts %>%
  pivot_longer(cols = c("REF", "HET", "ALT"),
               names_to = "Genotype",
               values_to = "Count")

ggplot(counts_long, aes(x = Superpopulation, y = Count, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Genotype Distribution by Superpopulation",
       x = "Superpopulation", y = "Count") +
  scale_fill_manual(values = c("REF" = "#1b9e77", "HET" = "#d95f02", "ALT" = "#7570b3")) +
  theme_minimal()

# 📐 12. Chi-square Test
geno_matrix <- counts %>%
  group_by(Superpopulation) %>%
  summarise(across(c(REF, HET, ALT), sum)) %>%
  column_to_rownames("Superpopulation") %>%
  as.matrix()


chi_result <- chisq.test(geno_matrix)
print(chi_result)

# 📉 13. Mosaic Plot of Chi-square Result
mosaic(geno_matrix, shade = TRUE, legend = TRUE,
       main = "Genotype Distribution by Superpopulation")

# 🔥 14. Allele Frequency Heatmap
freq_table <- geno_meta %>%
  filter(!is.na(Genotype)) %>%
  group_by(ID, Superpopulation) %>%
  summarise(ALT_freq = sum(Genotype %in% c("ALT", "HET")) / n(),
            .groups = "drop")

freq_matrix <- freq_table %>%
  pivot_wider(names_from = Superpopulation, values_from = ALT_freq) %>%
  column_to_rownames("ID") %>%
  as.matrix()

# Remove variants with all NA or zero variance
freq_matrix_filtered <- freq_matrix[rowSums(is.na(freq_matrix)) == 0, ]  # remove rows with NAs
freq_matrix_filtered <- freq_matrix_filtered[apply(freq_matrix_filtered, 1, var) > 0, ]  # remove no-variance rows

# Select top 20 most variable variants (you can change to top 50, etc.)
top_variants <- order(apply(freq_matrix_filtered, 1, var), decreasing = TRUE)[1:20]
freq_matrix_top <- freq_matrix_filtered[top_variants, ]

# Now plot a clean heatmap
pheatmap(freq_matrix_top, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Top 20 Variant Allele Frequencies by Superpopulation")

# From your full filtered genotype dataset (geno_meta)
all_variant_ids <- unique(geno_meta$ID)

variant_info <- strsplit(all_variant_ids, "_")
variant_df <- do.call(rbind, variant_info) %>% as.data.frame()
colnames(variant_df) <- c("CHR", "POS", "REF", "ALT")
variant_df$POS <- as.integer(as.character(variant_df$POS))

# Optional: Save for annotation
write.table(variant_df, "all_variants.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#Annotate All Variants Using biomaRt
variant_df <- unique(geno[, .(CHR, POS, REF, ALT)])
variant_df[, `:=`(POS = as.integer(POS), start = POS, end = POS)]

vep_input_chr_pos <- paste(variant_df$CHR, variant_df$POS, sep = "\t")
writeLines(vep_input_chr_pos, "vep_input_chr_pos.txt")

# Save raw annotation
write.csv(annot, "all_variants_annotated.csv", row.names = FALSE)
