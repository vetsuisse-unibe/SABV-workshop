#' ---
#' title: "Dataset Retrieval and Preparation: GTEx Liver"
#' subtitle: "EU-SABV Workshop"
#' author: "Vidhya Jagannathan, University of Bern, Switzerland"
#' format:
#'   html:
#'     toc: true
#'     toc-depth: 3
#'     toc-location: left
#'     number-sections: true
#'     code-fold: false
#'     code-tools: true
#'     theme: cosmo
#'     css: styles.css
#' execute:
#'   echo: true
#'   warning: true
#'   message: true
#'   error: true
#'   eval: true
#' ---
#' 
#' # Overview
#' 
#' **Script 01: Download and Prepare GTEx Liver Data**
#' 
#' This script downloads GTEx V8 liver tissue RNA-seq data and prepares a workshop subset for the hands-on exercises.
#' 
#' The full GTEx v8 gene read-count matrix contains **56,200 genes × 17,382 RNA-seq samples** across many tissues. For liver, the metadata contains **251 samples (178 Male, 73 Female)**; those present in the count matrix are used to create a balanced **50 Male / 50 Female** workshop subset. Liver is a strongly sexually dimorphic tissue.
#' 
#' The prepared file, `workshop_gtex_liver.RData`, is available in the folder `hands_on_datasets`.
#' 
#' ## 0. Install/Load Packages
#' 
## -----------------------------------------------------------------------------
required_pkgs <- c("data.table", "dplyr", "readr", "curl", "BiocManager")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# DESeq2 for later use (install now to avoid delays during workshop)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", ask = FALSE)
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", ask = FALSE)
}

library(data.table)
library(dplyr)
library(readr)

#' 
#' ## 1. Set Up Directories
#' 
## -----------------------------------------------------------------------------
# Raw GTEx files live in ../gtex_data/ — one level above quarto_tutorials/
data_dir <- file.path("..", "gtex_data")
dir.create(data_dir, showWarnings = FALSE)

#' 
#' ## 2. Download GTEx V8 Data
#' 
#' Files come from the [GTEx Portal](https://gtexportal.org/home/datasets) (V8 release). Each download is skipped if the file already exists. The gene-read-count matrix is ~875 MB, so the first run can take several minutes.
#' 
#' ## 2a. Gene read counts (all tissues)
#' 
## -----------------------------------------------------------------------------
counts_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
counts_file <- file.path(data_dir, "GTEx_gene_reads.gct.gz")

if (!file.exists(counts_file)) {
  curl::curl_download(counts_url, counts_file, quiet = FALSE)
}

#' 
#' ## 2b. Sample attributes (includes tissue type, sex, etc.)
#' 
## -----------------------------------------------------------------------------
sample_url <- "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
sample_file <- file.path(data_dir, "GTEx_SampleAttributes.txt")

if (!file.exists(sample_file)) {
  download.file(sample_url, sample_file, mode = "wb")
}

#' 
#' ## 2c. Subject phenotypes (includes sex, age, death classification)
#' 
## -----------------------------------------------------------------------------
subject_url <- "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
subject_file <- file.path(data_dir, "GTEx_SubjectPhenotypes.txt")

if (!file.exists(subject_file)) {
  download.file(subject_url, subject_file, mode = "wb")
}

#' 
#' ## 3. Load and Parse Metadata
#' 
## -----------------------------------------------------------------------------
# Subject phenotypes (SEX: 1 = Male, 2 = Female)
subjects <- read_tsv(subject_file, show_col_types = FALSE)
subjects <- subjects %>%
  mutate(
    SEX_LABEL = ifelse(SEX == 1, "Male", "Female"),
    AGE_DECADE = AGE   # AGE is stored as a decade range, e.g. "60-69"
  )
print(sprintf("Total subjects: %d (%d M, %d F)",
              nrow(subjects), sum(subjects$SEX == 1), sum(subjects$SEX == 2)))

# Sample attributes -> liver samples only
samples <- read_tsv(sample_file, show_col_types = FALSE)
liver_samples <- samples %>%
  filter(SMTSD == "Liver") %>%
  mutate(SUBJID = gsub("^(GTEX-[^-]+).*", "\\1", SAMPID))
print(sprintf("Total liver samples: %d", nrow(liver_samples)))

# Merge with subject phenotypes
liver_meta <- liver_samples %>% left_join(subjects, by = "SUBJID")
print(sprintf("Liver samples with phenotype data: %d (%d M, %d F)",
              nrow(liver_meta),
              sum(liver_meta$SEX == 1, na.rm = TRUE),
              sum(liver_meta$SEX == 2, na.rm = TRUE)))

# Restrict to samples actually present in the RNA-seq count matrix
# (read only the GCT header first).
counts_header <- names(fread(counts_file, skip = 2, nrows = 0, header = TRUE))
count_sample_ids <- counts_header[-c(1, 2)]   # drop Name and Description columns
liver_meta <- liver_meta %>% filter(SAMPID %in% count_sample_ids)
print(sprintf("Liver samples matched in count matrix: %d (%d M, %d F)",
              nrow(liver_meta),
              sum(liver_meta$SEX == 1, na.rm = TRUE),
              sum(liver_meta$SEX == 2, na.rm = TRUE)))

#' 
#' ## 4. Create a GTEx Liver Subset
#' 
#' Create a GTEx liver cohort for the workshop by subsetting to 50 male and 50 female samples that are roughly age-balanced.
#' 
## -----------------------------------------------------------------------------
set.seed(42)  # for reproducibility

male_samples <- liver_meta %>% filter(SEX == 1)
female_samples <- liver_meta %>% filter(SEX == 2)

# Select 50 of each, balanced across age decades where possible
select_balanced <- function(df, target_n = 50) {
  if (nrow(df) <= target_n) {
    return(df)
  }

  # Try to balance by age decade
  age_decades <- unique(df$AGE_DECADE)
  n_decades <- length(age_decades)
  target_per_decade <- ceiling(target_n / n_decades)

  # Sample from each decade explicitly (avoids n() in min() issue)
  selected_list <- list()
  for (decade in age_decades) {
    decade_data <- df %>% filter(AGE_DECADE == decade)
    n_to_sample <- min(nrow(decade_data), target_per_decade)
    if (n_to_sample > 0) {
      sampled <- decade_data %>% slice_sample(n = n_to_sample)
      selected_list[[decade]] <- sampled
    }
  }

  selected <- bind_rows(selected_list)

  # If we got more than target_n, randomly trim
  if (nrow(selected) > target_n) {
    selected <- selected %>% slice_sample(n = target_n)
  }

  # If we got fewer than target_n, add more randomly
  if (nrow(selected) < target_n) {
    remaining <- df %>% filter(!SAMPID %in% selected$SAMPID)
    n_more <- target_n - nrow(selected)
    if (nrow(remaining) > 0) {
      extra <- remaining %>% slice_sample(n = min(nrow(remaining), n_more))
      selected <- bind_rows(selected, extra)
    }
  }

  return(selected)
}

selected_males <- select_balanced(male_samples, 50)
selected_females <- select_balanced(female_samples, 50)
selected_samples <- bind_rows(selected_males, selected_females)

print(sprintf("Selected %d males and %d females",
            nrow(selected_males), nrow(selected_females)))

cat("Age distribution in selected samples:\n")
print(table(selected_samples$SEX_LABEL, selected_samples$AGE_DECADE))

#' 
#' ## 5. Load and Subset Gene Expression Counts
#' 
## -----------------------------------------------------------------------------
# GCT format: skip 2 header lines, then tab-delimited with Name, Description, samples...
counts_raw <- fread(counts_file, skip = 2, header = TRUE)

print(sprintf("Full dataset: %d genes × %d samples",
            nrow(counts_raw), ncol(counts_raw) - 2))

# Extract gene info
gene_info <- counts_raw[, .(Name, Description)]
colnames(gene_info) <- c("gene_id", "gene_name")

# Subset to our selected samples. These were pre-matched to the count matrix,
# so all selected IDs should be available here.
sample_ids <- selected_samples$SAMPID
available_ids <- intersect(sample_ids, colnames(counts_raw))

print(sprintf("Matched %d of %d selected samples in counts matrix",
            length(available_ids), length(sample_ids)))

if (length(available_ids) != length(sample_ids)) {
  warning("Some selected samples were not found in the count matrix; dropping unmatched samples.")
  selected_samples <- selected_samples %>% filter(SAMPID %in% available_ids)
}

# Subset counts
counts_subset <- counts_raw[, c("Name", available_ids), with = FALSE]
rownames_vec <- counts_subset$Name
counts_matrix <- as.matrix(counts_subset[, -1])
rownames(counts_matrix) <- rownames_vec

#' 
#' ## 6. Filter to Expressed Genes
#' 
#' We filter out low-expressed genes, keeping those with a mean count > 10.
#' 
## -----------------------------------------------------------------------------
# Keep only genes with a baseline expression level (mean count > 10).
expressed_idx <- which(rowMeans(counts_matrix) > 10)
counts_final <- counts_matrix[expressed_idx, ]
genes_final  <- gene_info[expressed_idx, ]
print(sprintf("Genes passing the expression filter: %d", nrow(counts_final)))

# Check if following X/Y-linked genes and other interesting sex-biased immune
# and metabolic genes survived the filter.
check_genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY",
                 "KDM6A", "DDX3X", "KDM5C", "EIF1AX", "ZFX", "USP9X", "SMC1A",
                 "CYP3A4", "CYP1A2", "ACE2", "TLR7", "FOXP3",
                 "IL6", "TNF", "ESR1", "AR", "CYP19A1")
survived <- check_genes %in% genes_final$gene_name
print(sprintf("Genes of interest surviving the filter: %d of %d",
              sum(survived), length(check_genes)))
if (any(!survived)) print(sprintf("  Dropped: %s", paste(check_genes[!survived], collapse = ", ")))

#' 
#' ## 7. Create Clean Workshop Metadata
#' 
## -----------------------------------------------------------------------------
workshop_meta <- selected_samples %>%
  filter(SAMPID %in% colnames(counts_final)) %>%
  dplyr::select(
    sample_id = SAMPID,
    subject_id = SUBJID,
    sex = SEX_LABEL,
    sex_numeric = SEX,
    age_decade = AGE_DECADE,
    ischemic_time = SMTSISCH,
    rna_integrity = SMRIN,
    tissue = SMTSD
  ) %>%
  mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    # Create a simulated "treatment" for the factorial design exercise
    # This allows the sex:treatment interaction demo
    # We'll assign based on ischemic time median split as a proxy
    ischemic_group = ifelse(ischemic_time > median(ischemic_time, na.rm = TRUE),
                            "High_Ischemic", "Low_Ischemic")
  ) %>%
  as.data.frame()

rownames(workshop_meta) <- workshop_meta$sample_id

# Ensure column order matches counts
counts_final <- counts_final[, workshop_meta$sample_id]

print(sprintf("Final dataset: %d genes × %d samples (%d M, %d F)",
              nrow(counts_final), ncol(counts_final),
              sum(workshop_meta$sex == "Male"),
              sum(workshop_meta$sex == "Female")))

#' 
#' ## 8. Save Workshop Data
#' 
## -----------------------------------------------------------------------------
# Save the prepared dataset where the exercises load it from.
save(counts_final, workshop_meta, genes_final,
     file = file.path("..", "hands_on_datasets", "workshop_gtex_liver.RData"))
cat("Saved: hands_on_datasets/workshop_gtex_liver.RData\n")

#' 
#' ## 9. Preliminary QC
#' 
## -----------------------------------------------------------------------------
# Check XIST (should be female-biased)
xist_idx <- which(genes_final$gene_name == "XIST")
if (length(xist_idx) > 0) {
  xist_male <- mean(counts_final[xist_idx, workshop_meta$sex == "Male"])
  xist_female <- mean(counts_final[xist_idx, workshop_meta$sex == "Female"])
  print(sprintf("XIST mean counts - Male: %.1f, Female: %.1f  ✓",
              xist_male, xist_female))
}

# Check a Y-linked gene
y_genes <- c("RPS4Y1", "DDX3Y", "KDM5D")
for (yg in y_genes) {
  yg_idx <- which(genes_final$gene_name == yg)
  if (length(yg_idx) > 0) {
    yg_male <- mean(counts_final[yg_idx, workshop_meta$sex == "Male"])
    yg_female <- mean(counts_final[yg_idx, workshop_meta$sex == "Female"])
    print(sprintf("%s mean counts - Male: %.1f, Female: %.1f  ✓",
                yg, yg_male, yg_female))
    break
  }
}

#' 
#' The prepared dataset now lives in `hands_on_datasets/`; tutorial 2b (`02b_workshop_gtex_liver_exploration`) loads it for the hands-on exercises.
