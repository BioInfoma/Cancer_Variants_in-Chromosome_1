# Install these once if not installed
install.packages("pheatmap")


library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)


# Load PCA
pca <- fread("pca_chr1.eigenvec")
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

# Load metadata

meta <- fread("integrated_call_samples_v3.20130502.ALL.panel", skip = 1,
              col.names = c("IID", "Population", "Superpopulation", "Gender"))


# Merge PCA + metadata
merged <- inner_join(pca, meta, by = "IID")

#Plot PCA
ggplot(merged, aes(x = PC1, y = PC2, color = Superpopulation)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA of Chromosome 1 Cancer Variants", x = "PC1", y = "PC2")


# Plot Donut Chart (Sample Distribution)
pop_counts <- meta %>% count(Superpopulation)

ggplot(pop_counts, aes(x = "", y = n, fill = Superpopulation)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Sample Distribution by Superpopulation")

#Prepare Genotype Matrix

geno <- fread("genotypes.tsv")
colnames(geno)[1:5] <- c("CHR", "POS", "ID", "REF", "ALT")

# Reshape geno to long format 
geno_long <- geno %>%
  pivot_longer(
    cols = starts_with("V"),  # all sample columns like V6, V7, ...
    names_to = "Sample",
    values_to = "GT"
  )

# Step 2: Fix sample names (replace "V6", "V7", etc. with real sample IDs)
col_names <- colnames(geno)[6:ncol(geno)]  # should be 2504 columns
sample_ids <- meta$IID                    # should be 2503 sample IDs

# Drop last column name to match sample_ids length
col_names_corrected <- col_names[-length(col_names)]  # now 2503
name_map <- setNames(sample_ids, col_names_corrected)

# Apply mapping
geno_long$Sample <- name_map[as.character(geno_long$Sample)]

# Step 3: Join with sample metadata
geno_meta <- inner_join(geno_long, meta, by = c("Sample" = "IID"))

# Step 4: Categorize genotype into REF, HET, ALT
geno_meta <- geno_meta %>%
  mutate(Genotype = case_when(
    GT %in% c("0/0", "0|0") ~ "REF",
    GT %in% c("0/1", "1/0", "0|1", "1|0") ~ "HET",
    GT %in% c("1/1", "1|1") ~ "ALT",
    TRUE ~ NA_character_
  ))

# counts <- summarized table with columns: Superpopulation, ALT, HET, REF

# First, pivot the data to long format for ggplot
counts_long <- counts %>%
  pivot_longer(cols = c("REF", "HET", "ALT"), 
               names_to = "Genotype", 
               values_to = "Count")

# Plot
ggplot(counts_long, aes(x = Superpopulation, y = Count, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Genotype Distribution by Superpopulation",
       x = "Superpopulation",
       y = "Number of Genotypes") +
  scale_fill_manual(values = c("REF" = "#1b9e77", "HET" = "#d95f02", "ALT" = "#7570b3")) +  # Color-blind friendly
  theme_minimal(base_size = 14)




#chi square test:

# Count genotype per population
counts <- geno_meta %>%
  filter(!is.na(Genotype)) %>%
  count(ID, Superpopulation, Genotype) %>%
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0)

# Chi-square test for each variant
# Count genotype per population


# Reformat into a matrix of genotype counts by Superpopulation
geno_matrix <- counts %>%
  select(Superpopulation, REF, HET, ALT) %>%
  column_to_rownames("Superpopulation") %>%
  as.matrix()

# Run chi-square test
chi_result <- chisq.test(geno_matrix)

# See result
chi_result

chi_result$FDR <- p.adjust(chi_result$p_value, method = "BH")
chi_result$Significant <-chi_result$FDR < 0.05

#Significance Plot
# Load the required package
library(vcd)  # for mosaic plot


# Create mosaic plot
mosaic(
  geno_matrix,
  shade = TRUE,
  legend = TRUE,
  main = "Distribution of Genotype Counts by Superpopulation"
)

#Allele Frequency Heatmap

colnames(geno)[1:5] <- c("CHR", "POS", "ID", "REF", "ALT")
geno$ID <- paste0(geno$CHR, "_", geno$POS, "_", geno$REF, "_", geno$ALT)

# Recompute freq_table
geno_long <- melt(geno, id.vars = 1:5, variable.name = "Sample", value.name = "GT")
geno_meta <- inner_join(geno_long, meta, by = c("Sample" = "IID"))

geno_meta <- geno_meta %>%
  mutate(Genotype = case_when(
    GT %in% c("0/0", "0|0") ~ "REF",
    GT %in% c("0/1", "1/0", "0|1", "1|0") ~ "HET",
    GT %in% c("1/1", "1|1") ~ "ALT",
    TRUE ~ NA_character_
  ))

freq_table <- geno_meta %>%
  filter(!is.na(Genotype)) %>%
  group_by(ID, Superpopulation) %>%
  summarise(
    ALT_freq = sum(Genotype == "ALT" | Genotype == "HET") / n(),
    .groups = "drop"
  )

# Now pivot and plot
freq_matrix <- freq_table %>%
  pivot_wider(names_from = Superpopulation, values_from = ALT_freq) %>%
  column_to_rownames("ID") %>%
  as.matrix()

pheatmap(freq_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Allele Frequencies of Cancer Variants")
