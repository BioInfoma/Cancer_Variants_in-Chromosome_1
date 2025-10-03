# 📦 1. Load Required Libraries
if (!require("httr")) install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("dplyr")) install.packages("dplyr")
if (!require("data.table")) install.packages("data.table")

library(httr)
library(jsonlite)
library(dplyr)
library(data.table)

#✅ Load variant_df from earlier step
geno <- fread("genotypes.tsv", header = FALSE)
sample_ids <- fread("sample_ids.txt", header = FALSE)[[1]]
colnames(geno)[1:4] <- c("CHR", "POS", "REF", "ALT")
colnames(geno)[5:(4 + length(sample_ids))] <- sample_ids
geno$ID <- paste0(geno$CHR, "_", geno$POS, "_", geno$REF, "_", geno$ALT)

variant_df <- geno[, .(CHR, POS)] %>% unique()
variant_df[, POS := as.integer(POS)]
variant_df[, CHR_POS := paste0(CHR, ":", POS)]

# ✅ Load previous results if they exist
if (file.exists("partial_annotation_results.tsv")) {
  rs_df <- fread("partial_annotation_results.tsv")
  done_positions <- rs_df$CHR_POS
} else {
  rs_df <- data.table()
  done_positions <- character(0)
}

# ✅ Define Query Function
query_ensembl <- function(pos) {
  url <- paste0("https://rest.ensembl.org/overlap/region/human/", pos, "?feature=variation")
  res <- GET(url, content_type("application/json"))
  if (status_code(res) == 200) {
    dat <- fromJSON(content(res, "text", encoding = "UTF-8"))
    if (length(dat) > 0) {
      data.frame(CHR_POS = pos,
                 rsID = sapply(dat, function(x) x$id),
                 consequence = sapply(dat, function(x) paste(x$consequence_terms, collapse = ",")),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(CHR_POS = pos, rsID = NA, consequence = NA)
    }
  } else {
    data.frame(CHR_POS = pos, rsID = NA, consequence = NA)
  }
}

# ✅ Run Annotation in Batches with Resume Support
positions <- setdiff(variant_df$CHR_POS, done_positions)
batch_size <- 1000  # ✅ Set higher batch size for speed
all_results <- list()

for (i in seq(1, length(positions), by = batch_size)) {
  cat("🔄 Processing:", i, "of", length(positions), "\n")
  batch <- positions[i:min(i + batch_size - 1, length(positions))]
  result_batch <- lapply(batch, query_ensembl)
  result_df <- do.call(rbind, result_batch)
  
  # ✅ Append to cumulative results
  rs_df <- rbind(rs_df, result_df)
  
  # ✅ Save progress to disk after each batch
  fwrite(rs_df, "partial_annotation_results.tsv", sep = "\t")
  
  Sys.sleep(1)  # API rate respect
}

# ✅ Combine Final Results
fwrite(rs_df, "annotated_variants_full.tsv", sep = "\t")

# ✅ Merge Annotation Back
annotated_variants <- left_join(variant_df, rs_df, by = "CHR_POS")

# ✅ Optional: Filter for Functional Variants
functional_hits <- annotated_variants %>%
  filter(!is.na(consequence)) %>%
  filter(grepl("missense|stop_gained|frameshift|splice", consequence))

# ✅ Optional: Export rsIDs to Lookup in ClinVar
write.table(functional_hits$rsID,
            "clinvar_rsids.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ✅ Summary Stats
summary_stats <- table(unlist(strsplit(paste(annotated_variants$consequence, collapse = ","), ",")))
print(summary_stats)