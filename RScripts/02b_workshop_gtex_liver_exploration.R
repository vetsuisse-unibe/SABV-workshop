#' ---
#' title: "Hands-on Exercise: GTEx Liver SABV Exploration — Age"
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
#' > Hands-on tutorial for the GTEx liver SABV exercises: the variability myth, PCA, sex-biased expression, factorial sex × age modeling, and the common analysis mistakes (pooling, disaggregation, subgroup p-value comparison).
#' 
#' # Overview
#' 
#' Script 02b: GTEx Liver Exploration — Sex as Biological Variable in Genomics
#' 
#' Uses the GTEx liver subset prepared by tutorial 01 (`01_download_prepare_GTEx_liver`). The design-phase power simulations live in the companion tutorial `02a_workshop_simulation_exercises`.
#' 
#' ## Common analysis pitfalls and the factorial-model solution
#' 
#' This tutorial demonstrates why pooling, disaggregation, and comparing p-values across subgroups are not valid substitutes for a factorial model (the **SIRF** framework: Karp et al. 2025, *Nature Communications*; the misconception percentages quoted throughout come from the barriers survey of Gaskill et al. 2025, *eLife*).
#' 
#' ::: {.callout-note}
#' ## Primary literature for this module
#' 
#' - **Reynolds (2024)**, *A Guide to Sample Size for Animal-based Studies* (Wiley). Source for: the experimental vs biological unit, pseudo-replication, effect sizes, the lower power of interaction terms relative to main effects, factorial/screening efficiency, and "right-sizing not significance" (3Rs).
#' - **Feuvrier, Bohacek & Germain (2026, draft)**, *Addressing pitfalls of differential -omics response analysis by investigating sex-specificity of stress-dependent transcriptional changes in the hippocampus.* Source for: the Venn-overlap fallacy, baseline/intercept coefficient misinterpretation, and the no-intercept contrast model.
#' - **Karp et al. (2025)**, *Nature Communications* — the **SIRF** framework for analysing sex as a biological variable. **Gaskill et al. (2025)**, *eLife* — the barriers/misconceptions survey behind the percentage claims (81% pooling, 80% difference-in-significance, 93% variability). **Garcia-Sifuentes & Maney (2021)**; **Nieuwenhuis et al. (2011)**, *Nature Neuroscience*.
#' :::
#' 
#' ## Load Packages
#' 
## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggrepel)
library(ggVennDiagram)   # Venn diagrams
library(patchwork)       # combine plots with + and / (see install instructions)

# Set a consistent theme
theme_set(theme_bw(base_size = 14))

#' 
#' ## Load Data
#' 
## -----------------------------------------------------------------------------
data_file <- file.path("..", "hands_on_datasets", "workshop_gtex_liver.RData")
load(data_file)

print(sprintf("Counts: %d genes × %d samples (%d M, %d F)",
              nrow(counts_final), ncol(counts_final),
              sum(workshop_meta$sex == "Male"),
              sum(workshop_meta$sex == "Female")))

#' 
#' # Testing misconceptions on real data: are females more variable?
#' 
#' We now leave the simulations behind (companion tutorial `02a`) and put a common claim to the test on the real GTEx liver data: *is female gene expression inherently more variable than male?* This is the reasoning most often used to justify studying only one sex — so it is worth checking directly.
#' 
## -----------------------------------------------------------------------------
# KEY MISCONCEPTION ADDRESSED: "Female data is inherently more variable"
# (Gaskill et al. 2025: 93% incorrectly believed inclusion increases variability).
# We test this misconception with real GTEx liver expression data.

male_idx   <- workshop_meta$sex == "Male"
female_idx <- workshop_meta$sex == "Female"

# Variance-stabilize the counts (VST) so a gene's spread is comparable across the
# whole expression range. Raw counts are heteroscedastic — highly expressed genes
# vary more in absolute terms — which would otherwise dominate the comparison.
dds_var <- DESeqDataSetFromMatrix(counts_final, workshop_meta, design = ~ 1)
vsd_mat <- assay(vst(dds_var, blind = TRUE))

# On the VST (log-like) scale we compare each gene's standard deviation within
# each sex. We use SD, not the coefficient of variation: CV (sd/mean) is only
# meaningful on a ratio scale (raw counts), not on a log scale.
sd_male   <- apply(vsd_mat[, male_idx],   1, sd)
sd_female <- apply(vsd_mat[, female_idx], 1, sd)

var_comparison <- data.frame(sd_male = sd_male, sd_female = sd_female)
var_comparison$more_variable <- ifelse(var_comparison$sd_female > var_comparison$sd_male,
                                       "Female more variable",
                                       "Male more variable")

sprintf("Females more variable: %d genes (%.1f%%) | Males more variable: %d genes (%.1f%%).",
        sum(var_comparison$more_variable == "Female more variable"),
        mean(var_comparison$more_variable == "Female more variable") * 100,
        sum(var_comparison$more_variable == "Male more variable"),
        mean(var_comparison$more_variable == "Male more variable") * 100)

# Paired test: the same genes, male vs female spread.
wtest <- wilcox.test(var_comparison$sd_male, var_comparison$sd_female, paired = TRUE)
sprintf("Paired Wilcoxon test p-value: %.4f", wtest$p.value)

# Each point = one gene; points on the red line have equal male/female spread.
p_variability <- ggplot(var_comparison, aes(x = sd_male, y = sd_female)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Gene Expression Variability: Males vs Females",
    subtitle = "GTEx Liver (VST scale) | Each point = one gene | Red line = equal variability",
    x = "Standard Deviation (Males)",
    y = "Standard Deviation (Females)"
  ) +
  coord_fixed(clip = "off") +
  annotate("text", x = max(var_comparison$sd_male, na.rm = TRUE) * 0.95,
           y = max(var_comparison$sd_female, na.rm = TRUE) * 0.15,
           label = sprintf("Female > Male: %.0f%%\nMale > Female: %.0f%%",
                           mean(var_comparison$more_variable == "Female more variable") * 100,
                           mean(var_comparison$more_variable == "Male more variable") * 100),
           size = 4, hjust = 1) +
  theme(plot.margin = margin(10, 35, 10, 10))
print(p_variability)

#' 
#' Each point is a gene; the dashed line marks equal male/female spread. The cloud sits **below** the line for most genes — in liver, variability is actually **higher in males** (more variable in ~67% of genes, vs 33% for females). The per-gene differences are small, but the direction is clear: if either sex is "noisier" here it is **males**, the opposite of the common assumption that female data are more variable. The paired Wilcoxon test is significant (tens of thousands of genes give enormous power), but that reflects a small, consistent lean rather than a large effect. **Bottom line: female expression is not the more variable sex.**
#' 
#' This matches the rodent literature that tackled the "oestrous cycle makes females noisier" argument head-on: across traits, females are **not** more variable than males, and unstaged females are no more variable than males (Prendergast et al. 2014; Beery & Zucker 2011) — so the cycle is not a valid reason to exclude them.
#' 
#' # PART 3: HANDS-ON ANALYSIS
#' 
#' ## EXERCISE 3.1: Data Exploration
#' 
## -----------------------------------------------------------------------------
# Create DESeqDataSet for normalization
dds <- DESeqDataSetFromMatrix(
  countData = counts_final,
  colData = workshop_meta,
  design = ~ sex
)

# Variance stabilizing transformation for PCA
vsd <- vst(dds, blind = TRUE)

#' 
#' ## PCA
#' 
#' The quickest look is DESeq2's built-in `plotPCA()` — one line, and enough to answer the basic question: **do samples separate by sex?**
#' 
## -----------------------------------------------------------------------------
plotPCA(vsd, intgroup = "sex")

#' 
#' Samples do **not** split cleanly by sex — as expected, since sex-linked genes are only a small fraction of the transcriptome. One thing still worth checking is whether **age is confounded with sex** in this subset. The optional block below redraws the PCA with custom styling and colours the points by age, for anyone who wants a publication-style figure.
#' 
#' ::: {.callout-note collapse="true"}
#' ## Optional: a nicer-looking PCA (and the age-confounding check)
#' 
## -----------------------------------------------------------------------------
pca_data <- plotPCA(vsd, intgroup = c("sex"), returnData = TRUE)
pca_var <- attr(pca_data, "percentVar")

# Enhanced PCA plot, coloured and shaped by sex
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sex, shape = sex)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Male" = "#2166AC", "Female" = "#B2182B")) +
  labs(
    title = "PCA of GTEx Liver Samples",
    subtitle = "Do samples separate by sex?",
    x = sprintf("PC1 (%.1f%% variance)", pca_var[1] * 100),
    y = sprintf("PC2 (%.1f%% variance)", pca_var[2] * 100),
    color = "Sex", shape = "Sex"
  )
print(p1)

# Coloured by age decade, shaped by sex: is age confounded with sex here?
pca_data$age_decade <- workshop_meta[pca_data$name, "age_decade"]
p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = age_decade, shape = sex)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "YlOrRd") +
  labs(
    title = "PCA colored by Age Decade, shaped by Sex",
    subtitle = "Check: Is age confounded with sex in our subset?",
    x = sprintf("PC1 (%.1f%% variance)", pca_var[1] * 100),
    y = sprintf("PC2 (%.1f%% variance)", pca_var[2] * 100)
  )
print(p2)

#' :::
#' 
#' ## EXERCISE 3.2: Identify sex chromosome genes
#' 
## -----------------------------------------------------------------------------
# XIST (X-inactivation marker — expressed in females)
xist_idx <- which(genes_final$gene_name == "XIST")
# Y-linked genes
y_gene_names <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")
y_gene_idx <- which(genes_final$gene_name %in% y_gene_names)

# X-inactivation escapee genes
escape_genes <- c("KDM6A", "DDX3X", "KDM5C")
escape_idx <- which(genes_final$gene_name %in% escape_genes)

# Build a data frame for plotting
sex_chr_genes <- c(xist_idx, y_gene_idx, escape_idx)
sex_chr_names <- genes_final$gene_name[sex_chr_genes]

# Get normalized counts
# Note: size factors must be estimated before normalized counts are available.
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

plot_data <- data.frame()
for (i in seq_along(sex_chr_genes)) {
  idx <- sex_chr_genes[i]
  gene_name <- genes_final$gene_name[idx]

  df <- data.frame(
    gene = gene_name,
    expression = norm_counts[idx, ],
    sex = workshop_meta$sex,
    category = case_when(
      gene_name == "XIST" ~ "X-inactivation (XIST)",
      gene_name %in% y_gene_names ~ "Y-linked",
      gene_name %in% escape_genes ~ "X-escape"
    )
  )
  plot_data <- rbind(plot_data, df)
}

p_sex_chr <- ggplot(plot_data, aes(x = sex, y = log2(expression + 1), fill = sex)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Male" = "#2166AC", "Female" = "#B2182B")) +
  labs(
    title = "Sex Chromosome Gene Expression in GTEx Liver",
    subtitle = "XIST: female-specific | Y genes: male-specific | X-escape: expressed from both X's",
    y = "log2(Normalized Counts + 1)",
    x = ""
  ) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
print(p_sex_chr)

#' 
#' These nine genes show the three textbook classes of sex-chromosome expression: **XIST** is female-specific (the X-inactivation switch, off in males); the **Y-linked genes** (DDX3Y, EIF1AY, KDM5D, RPS4Y1, UTY) are male-specific; and the **X-escape genes** (DDX3X, KDM5C, KDM6A) are modestly higher in females, which express them from both X chromosomes. These are dosage-driven differences seen in essentially every tissue, so they serve as **positive controls** — XIST up in females and Y-genes up in males confirm the sex labels and pipeline are correct. (The Y-genes aren't quite zero in females: a low signal from reads cross-mapping between the Y genes and their near-identical X paralogs — a common real-data artifact.)
#' 
#' ### Quantifying X-inactivation escape
#' 
#' The X-escape genes are expressed from *both* X chromosomes in females, so they run higher in females (typically 1.5–2×). We can put a number on it — and screen a wider panel of known escapees — with the female/male expression ratio. These are genuine biological sex differences, not technical artifacts.
#' 
## -----------------------------------------------------------------------------
escape_gene_list <- c("KDM6A", "DDX3X", "KDM5C", "XIST",
                       "EIF1AX", "ZFX", "USP9X", "SMC1A")

escape_results <- data.frame()
for (gene in escape_gene_list) {
  idx <- which(genes_final$gene_name == gene)
  if (length(idx) > 0) {
    male_expr <- mean(norm_counts[idx, workshop_meta$sex == "Male"])
    female_expr <- mean(norm_counts[idx, workshop_meta$sex == "Female"])
    ratio <- female_expr / male_expr

    escape_results <- rbind(escape_results, data.frame(
      gene = gene,
      male_mean = round(male_expr, 1),
      female_mean = round(female_expr, 1),
      F_M_ratio = round(ratio, 2),
      likely_escape = ifelse(ratio > 1.2 & gene != "XIST", "Yes",
                              ifelse(gene == "XIST", "X-inactivation marker", "No"))
    ))
  }
}

print("F/M ratio > 1.2 flags likely X-escape genes:")
print(escape_results, row.names = FALSE)

#' 
#' ## EXERCISE 3.3: Differential Expression with DESeq2 — sex as factor
#' 
## -----------------------------------------------------------------------------
#| fig-width: 12
#| fig-height: 6
# Simple model: ~ sex
dds_sex <- DESeq(dds)

res_sex <- results(dds_sex, contrast = c("sex", "Female", "Male"),
                    alpha = 0.05)

cat("DESeq2 results (Female vs Male):\n")
summary(res_sex)

# Annotate with gene names
res_sex_df <- as.data.frame(res_sex) %>%
  mutate(gene_id = rownames(res_sex)) %>%
  left_join(genes_final, by = c("gene_id" = "gene_id")) %>%
  arrange(padj)

# Top 20 sex-biased genes:
print(head(res_sex_df[, c("gene_name", "log2FoldChange", "padj")], 20))

# Volcano plot — colour points by the DIRECTION of the significant sex bias,
# so male-biased genes (higher in males, negative LFC) show up in blue.
res_sex_df$direction <- case_when(
  is.na(res_sex_df$padj) | res_sex_df$padj >= 0.05 ~ "Not significant",
  res_sex_df$log2FoldChange > 0                    ~ "Female-biased",
  TRUE                                             ~ "Male-biased"
)

# Label the key sex-linked genes plus the few most significant genes. Keeping the
# set small (and spread across both sides) avoids the crowded, smudged labels.
key_genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "KDM6A", "DDX3X", "UTY", "EIF1AY")
top_genes <- head(res_sex_df$gene_name[order(res_sex_df$padj)], 8)
label_genes <- unique(c(key_genes, top_genes))
res_sex_df$label <- ifelse(res_sex_df$gene_name %in% label_genes,
                           res_sex_df$gene_name, "")

# The most sex-linked genes (XIST, Y-genes) have astronomically small p-values that
# stretch the y-axis and squash everything else. Cap -log10(p) for display and flag the
# capped points with a triangle, so nothing is hidden — just compressed at the top.
y_cap <- 50
res_sex_df$neglog10p <- -log10(res_sex_df$pvalue)
res_sex_df$y_display <- pmin(res_sex_df$neglog10p, y_cap)
res_sex_df$capped    <- res_sex_df$neglog10p > y_cap

volcano_colors <- c("Male-biased" = "#2166AC", "Female-biased" = "#B2182B",
                    "Not significant" = "grey80")
point_sizes <- c("Male-biased" = 2.4, "Female-biased" = 2.4, "Not significant" = 1.1)

# --- Volcano plot (y-axis capped) ---
p_volcano <- ggplot(res_sex_df, aes(x = log2FoldChange, y = y_display)) +
  geom_point(aes(color = direction, size = direction, shape = capped), alpha = 0.6) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3, fontface = "bold", color = "black",
    segment.color = "grey60", segment.size = 0.3,
    min.segment.length = 0, box.padding = 0.8, point.padding = 0.3,
    max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = volcano_colors) +
  scale_size_manual(values = point_sizes, guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "shown", "TRUE" = sprintf("capped at %d", y_cap)),
                     name = NULL) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5, linewidth = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, linewidth = 0.7) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Volcano plot",
    subtitle = sprintf("Female vs Male | y capped at %d (triangles) | blue = male-biased, red = female-biased", y_cap),
    x = "log2 Fold Change (Female / Male)",
    y = "-log10(p-value), capped",
    color = ""
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 11),
    legend.position = "top"
  )

# --- MA plot: mean expression vs fold change (the standard companion to the volcano) ---
p_ma <- ggplot(res_sex_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = direction, size = direction), alpha = 0.6) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3, fontface = "bold", color = "black",
    segment.color = "grey60", segment.size = 0.3,
    min.segment.length = 0, box.padding = 0.6, max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = volcano_colors, guide = "none") +
  scale_size_manual(values = point_sizes, guide = "none") +
  scale_x_log10() +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.5, linewidth = 0.6) +
  labs(
    title = "MA plot",
    subtitle = "Mean expression vs fold change | same colours",
    x = "Mean of normalized counts (log scale)",
    y = "log2 Fold Change (Female / Male)"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 11)
  )

# Volcano (capped) and MA plot side by side.
print(p_volcano + p_ma + plot_layout(widths = c(1, 1)))

#' 
#' The two panels show the same result two ways: the **volcano** emphasises statistical significance (capped, so XIST and the Y-genes don't flatten everything else), while the **MA plot** shows that the strongest sex biases sit at genes with high mean expression. Both use the same colour scheme.
#' 
#' ## EXERCISE 3.4: Factorial Model — sex × age interaction
#' 
#' The factorial model `~ sex + condition + sex:condition` tests both main effects *and* their interaction in one fit. Here **condition = age group** (under 50 vs 50+). The model asks the question: does expression change with age, and does that change differ by sex?
#' 
## -----------------------------------------------------------------------------
# condition = age group (under 50 vs 50+); model: ~ sex + condition + sex:condition.
workshop_meta$condition <- factor(
  ifelse(as.numeric(sub("-.*", "", workshop_meta$age_decade)) >= 50, "Older", "Younger"),
  levels = c("Younger", "Older")
)

# Create DESeq dataset with factorial design
dds_factorial <- DESeqDataSetFromMatrix(
  countData = counts_final,
  colData = workshop_meta,
  design = ~ sex + condition + sex:condition  # THE KEY MODEL
)

dds_factorial <- DESeq(dds_factorial)

print(resultsNames(dds_factorial))

# Main effect of sex, adjusted for condition.
res_sex_adj <- results(dds_factorial, name = "sex_Female_vs_Male", alpha = 0.05)
summary(res_sex_adj)

# Main effect of condition.
res_condition <- results(dds_factorial, name = "condition_Older_vs_Younger",
                          alpha = 0.05)
summary(res_condition)

# Interaction: genes where the age effect differs between sexes.
res_interaction <- results(dds_factorial, name = "sexFemale.conditionOlder",
                            alpha = 0.05)
summary(res_interaction)

# Interaction volcano
res_int_df <- as.data.frame(res_interaction) %>%
  mutate(gene_id = rownames(res_interaction)) %>%
  left_join(genes_final, by = c("gene_id" = "gene_id")) %>%
  arrange(padj)

# Top interaction genes (age effect differs by sex):
print(head(res_int_df[, c("gene_name", "log2FoldChange", "padj")], 10))

# --- Real sex × age interaction plots: three GTEx liver genes ---
# Genes taken from the top of res_int_df (ranked by the interaction term). These are
# illustrative, SELECTION-DEPENDENT examples (see the caveats below the chunk), not
# independent confirmation. All three are FDR-significant for sexFemale.conditionOlder.
panel_genes <- c("IGKJ5", "SKI", "DACT2")
panel_ids   <- genes_final$gene_id[match(panel_genes, genes_final$gene_name)]

# Per-sample log2 normalized counts for the three genes.
nc <- counts(dds_factorial, normalized = TRUE)
panel_long <- do.call(rbind, lapply(seq_along(panel_genes), function(i) {
  data.frame(
    gene      = factor(panel_genes[i], levels = panel_genes),  # keep panel order
    sex       = workshop_meta$sex,
    condition = factor(workshop_meta$condition, levels = c("Younger", "Older")),
    log2count = log2(nc[panel_ids[i], ] + 1)
  )
}))

# Cell means +/- SEM on the log2 scale (the scale the model tested on, so
# "non-parallel lines = interaction" reads correctly).
panel_summary <- panel_long %>%
  group_by(gene, sex, condition) %>%
  summarise(mean = mean(log2count), se = sd(log2count) / sqrt(n()), .groups = "drop")

# SABV-standard interaction plot: group means by sex as separate lines, SEM error bars,
# raw points jittered behind. Non-parallel (here, crossing) lines = interaction.
p_interaction <- ggplot(panel_long,
        aes(x = condition, y = log2count, color = sex, group = sex)) +
  geom_point(position = position_jitter(width = 0.06), alpha = 0.25, size = 1.2) +
  geom_line(data = panel_summary, aes(y = mean), linewidth = 1.1) +
  geom_point(data = panel_summary, aes(y = mean), size = 3) +
  geom_errorbar(data = panel_summary,
                aes(y = mean, ymin = mean - se, ymax = mean + se),
                width = 0.08, linewidth = 0.7) +
  scale_color_manual(values = c("Male" = "#2166AC", "Female" = "#B2182B")) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(
    title = "Real sex × age interactions in GTEx liver",
    subtitle = paste0("Non-parallel (here, crossing) lines = the age effect differs by sex.\n",
                      "Three FDR-significant genes from res_interaction (selection-dependent illustration)."),
    y = "log2 normalized count", x = "Age group", color = "Sex"
  )
print(p_interaction)

#' 
#' These are **real** GTEx liver genes whose age effect differs by sex: the per-sex lines cross (IGKJ5, SKI, DACT2), the visual signature of a sex × age interaction. In this subset **35 genes pass FDR** (`padj < 0.05`) for the `sexFemale.conditionOlder` term, so interactions *are* detectable here — though far fewer than the abundant age *main* effects, consistent with the **interaction-power penalty** quantified in Exercise 2.1b (an interaction needs roughly 4× the N of a main effect of the same size).
#' 
#' Two cautions come with these panels. **(1) Selection bias:** these genes were chosen *because* they ranked top on the very interaction statistic being illustrated, so they demonstrate the pattern but are **not** independent confirmation — to *quantify* a sex × age effect you would validate the selected genes in an independent set of donors, not re-test the same samples used to rank them. **(2) Imbalance:** the female cells are unequal (15 Younger vs 35 Older, against 25/25 in males), and an imbalanced design hands more apparent power to the larger group (Feuvrier et al.) — a reminder that the remedy is balance and modelling, not gene-picking. Treat the panels as "here is what an interaction looks like in these data," then test and validate it properly.
#' 
#' ## EXERCISE 3.5: Model parameterization — baselines, intercepts, and contrasts
#' 
#' A factorial model isn't automatically interpreted correctly: *how* it is parameterized changes which coefficient answers which question (Feuvrier et al. 2026). One gene (KL, Klotho) makes the point.
#' 
#' ::: {.callout-warning}
#' ## The trap
#' 
#' A student fits `~ sex * condition`, reads the `conditionOlder` coefficient as "the age effect," then re-runs with males (not females) as the reference sex — and the coefficient **changes**. Is the model broken? What does `conditionOlder` actually estimate?
#' :::
#' 
## -----------------------------------------------------------------------------
# Self-contained example gene (KL, Klotho) as a continuous outcome.
kl_idx <- which(genes_final$gene_name == "KL")
age_grp <- factor(
  ifelse(as.numeric(sub("-.*", "", workshop_meta$age_decade)) >= 50, "Older", "Younger"),
  levels = c("Younger", "Older"))
param_df <- data.frame(
  expression = log2(counts_final[kl_idx, ] + 1),
  sex        = factor(workshop_meta$sex, levels = c("Female", "Male")),
  condition  = age_grp
)

# (1) Intercept model, FEMALE baseline: 'conditionOlder' = age effect IN FEMALES only.
fit_Fbase <- lm(expression ~ sex * condition, data = param_df)
co_F <- summary(fit_Fbase)$coefficients["conditionOlder", c(1, 4)]

# (2) Intercept model, MALE baseline: same data, releveled. 'conditionOlder' = age effect IN MALES.
param_df_M <- param_df; param_df_M$sex <- relevel(param_df_M$sex, ref = "Male")
fit_Mbase <- lm(expression ~ sex * condition, data = param_df_M)
co_M <- summary(fit_Mbase)$coefficients["conditionOlder", c(1, 4)]

sprintf("conditionOlder | Female-baseline: beta=%.3f p=%.4f  ||  Male-baseline: beta=%.3f p=%.4f",
        co_F[1], co_F[2], co_M[1], co_M[2])

# (3) No-intercept ("means / cell-means") model + explicit contrasts (Feuvrier et al.'s recommendation).
param_df$group <- factor(paste(param_df$sex, param_df$condition, sep = "_"))
fit_cell <- lm(expression ~ 0 + group, data = param_df)
b <- coef(fit_cell)
F_age   <- b["groupFemale_Older"] - b["groupFemale_Younger"]   # age effect in females
M_age   <- b["groupMale_Older"]   - b["groupMale_Younger"]     # age effect in males
avg_age <- ((b["groupFemale_Older"] + b["groupMale_Older"]) -
            (b["groupFemale_Younger"] + b["groupMale_Younger"])) / 2  # true average age effect
interaction <- F_age - M_age                                   # sex x age interaction
sprintf("Contrasts | age-in-F=%.3f  age-in-M=%.3f  AVERAGE age=%.3f  interaction(F-M)=%.3f",
        F_age, M_age, avg_age, interaction)

#' 
#' **The point.** The `conditionOlder` coefficient is **not** the overall age effect — it is the age effect *in the baseline sex*, which is exactly why swapping the reference changes it. Nothing is broken; the coefficient just doesn't mean what the student thought. The no-intercept model `~ 0 + group` plus explicit contrasts lets you ask for *exactly* the comparison you want — age in females, in males, the true **average**, or the **interaction** — with no baseline ambiguity.
#' 
## -----------------------------------------------------------------------------
# Visualize the four cell means as an interaction plot with SEM error bars
cell_summary <- param_df %>%
  group_by(sex, condition) %>%
  summarise(mean = mean(expression), se = sd(expression) / sqrt(n()), .groups = "drop")

p_param <- ggplot(cell_summary, aes(x = condition, y = mean, color = sex, group = sex)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.08, linewidth = 0.8) +
  scale_color_manual(values = c("Male" = "#2166AC", "Female" = "#B2182B")) +
  labs(title = "KL: cell means as an interaction plot",
       subtitle = "The slope difference between lines IS the sex x age interaction",
       x = "Age group", y = "log2(counts + 1)", color = "Sex")
print(p_param)

#' 
#' # Common analysis mistakes and the factorial-model solution
#' 
## -----------------------------------------------------------------------------
# Reference: Gaskill et al. 2025 — 81% thought pooling was OK.
# Reference: Garcia-Sifuentes & Maney 2021 — analysis errors are endemic.

# We first illustrate the three mistakes with two example genes, then repeat each
# one genome-wide. Both are autosomal and expressed in both sexes:
#   KL    — Klotho, an aging / longevity-associated gene
#   WFDC2 — WAP four-disulfide core domain 2 (HE4), a secreted epithelial protein
age_group <- factor(
  ifelse(as.numeric(sub("-.*", "", workshop_meta$age_decade)) >= 50, "Older", "Younger"),
  levels = c("Younger", "Older")
)
make_example <- function(gene_name) {
  idx <- which(genes_final$gene_name == gene_name)
  data.frame(expression = log2(counts_final[idx, ] + 1),
             sex        = workshop_meta$sex,
             condition  = age_group)
}
example_data  <- make_example("KL")
example_data2 <- make_example("WFDC2")

#' 
#' ## MISTAKE 1: POOLING — Ignoring sex entirely
#' 
#' Pooling fits `expression ~ condition` (here condition = age group) and ignores sex (81% of researchers in Gaskill et al. thought this was fine). KL (Klotho) differs markedly between the sexes in liver. Ignore sex and that variation is absorbed into the residual noise, inflating it enough to mask the age effect. Adding sex to the regression model recovers it.
#' 
## -----------------------------------------------------------------------------
# Pooled model (ignores sex) vs. model that accounts for sex.
fit_pooled <- lm(expression ~ condition, data = example_data)
residual_se_pooled <- summary(fit_pooled)$sigma
p_pooled <- summary(fit_pooled)$coefficients[2, 4]
sprintf("Model 1 (~ condition):       residual SE = %.4f, condition p = %.4f",
        residual_se_pooled, p_pooled)

fit_correct <- lm(expression ~ sex + condition, data = example_data)
residual_se_correct <- summary(fit_correct)$sigma
p_correct <- summary(fit_correct)$coefficients["conditionOlder", 4]
sprintf("Model 2 (~ sex + condition): residual SE = %.4f, condition p = %.4f",
        residual_se_correct, p_correct)

#' 
#' Ignoring sex left KL's age effect **non-significant** (p ≈ 0.12), but with sex in the model the same effect is clearly significant (p ≈ 0.02). Leave out *any* real source of variation and it does the same — inflating the noise.
#' 
#' **And it's not just one gene.** Test every gene for an age effect, first without sex in the model and then with it, and compare the resulting ranked gene lists.
#' 
## -----------------------------------------------------------------------------
# Genome-wide: does ignoring sex change which genes look "age-associated"?
res_pool <- results(DESeq(DESeqDataSetFromMatrix(counts_final, workshop_meta, ~ condition)),
                    name = "condition_Older_vs_Younger")
res_adj  <- results(DESeq(DESeqDataSetFromMatrix(counts_final, workshop_meta, ~ sex + condition)),
                    name = "condition_Older_vs_Younger")
ord_pool <- rownames(res_pool)[order(res_pool$padj)]
ord_adj  <- rownames(res_adj)[order(res_adj$padj)]
for (N in c(10, 50, 100)) {
  shared <- length(intersect(head(ord_pool, N), head(ord_adj, N)))
  cat(sprintf("Top %3d age genes: %d shared, %d change when sex is added\n", N, shared, N - shared))
}
cat(sprintf("Age genes at FDR<0.05: %d pooled vs %d sex-adjusted; %d 'age' hits become non-significant once sex is added to the model.\n",
            sum(res_pool$padj < 0.05, na.rm = TRUE), sum(res_adj$padj < 0.05, na.rm = TRUE),
            sum(res_pool$padj < 0.05 & res_adj$padj >= 0.05, na.rm = TRUE)))

#' 
#' Because the female samples skew older, sex and age are **confounded** in this subset. Ignore sex and some sex differences get counted as age effects: 7 of the top-10 "age genes" change, and the FDR<0.05 age list shrinks from **375 genes (pooled) to 240 (sex-adjusted)** — 162 "age" hits become non-significant once sex is added. So pooling distorts the age list **both ways** — inventing false age genes (most here) and hiding real ones like KL.
#' 
#' ## MISTAKE 2: DISAGGREGATION — Separate analyses by sex
#' 
#' Analysing males and females separately splits your cohort in two — halving the power in each group — and it never actually tests whether the sexes differ. Worse, it tempts the classic error: "significant in one sex but not the other, so the effect is sex-specific." But a **difference in significance is not a significant difference** (80% of researchers in Gaskill et al. fell for this; Nieuwenhuis et al. 2011, *Nature Neuroscience*).
#' 
## -----------------------------------------------------------------------------
male_data <- example_data %>% filter(sex == "Male")
female_data <- example_data %>% filter(sex == "Female")

fit_male <- lm(expression ~ condition, data = male_data)
fit_female <- lm(expression ~ condition, data = female_data)

p_male <- summary(fit_male)$coefficients[2, 4]
p_female <- summary(fit_female)$coefficients[2, 4]

sprintf("Males-only condition p = %.4f (n = %d) | Females-only condition p = %.4f (n = %d)",
        p_male, nrow(male_data), p_female, nrow(female_data))

#' 
#' It isn't limited to one gene. Run the **whole** age analysis separately in each sex:
#' 
## -----------------------------------------------------------------------------
# Genome-wide disaggregation: age genes found in a males-only vs females-only analysis.
res_age_M <- results(DESeq(DESeqDataSetFromMatrix(counts_final[, workshop_meta$sex == "Male"],
                          workshop_meta[workshop_meta$sex == "Male", ], ~ condition)),
                     name = "condition_Older_vs_Younger")
res_age_F <- results(DESeq(DESeqDataSetFromMatrix(counts_final[, workshop_meta$sex == "Female"],
                          workshop_meta[workshop_meta$sex == "Female", ], ~ condition)),
                     name = "condition_Older_vs_Younger")
cat(sprintf("Age genes (FDR<0.05): %d in a males-only analysis, %d in a females-only analysis.\n",
            sum(res_age_M$padj < 0.05, na.rm = TRUE), sum(res_age_F$padj < 0.05, na.rm = TRUE)))

#' 
#' At first glance this looks like a big sex difference — "age affects the female liver, not the male." It's an illusion: you can't conclude a sex difference by comparing two separate analyses, and the counts are unstable because each run uses only half the data.
#' 
#' **The tempting Venn diagram and convincing heatmap.** Venn diagrams are commonly used to compare gene lists from separate sub-group analysis. Heatmaps are also common in papers because they make expression patterns easy to see. Both are useful visualizations, but both are descriptive: neither one tests whether the age effect is statistically different between males and females.
#' 
#' Below, the Venn diagram compares the male-only and female-only age-hit lists. Next to it, the heatmap shows genes that appear significant in the female-only age analysis but not in the male-only analysis — exactly the kind of list that is often labelled "female-specific" in papers. The figure may look persuasive, but at this point we have still only compared separate subgroup analyses.
#' 
## -----------------------------------------------------------------------------
#| fig-width: 13
#| fig-height: 6
# Hits from separate age analyses in males and females.
male_age_hits <- rownames(res_age_M)[which(res_age_M$padj < 0.05)]
female_age_hits <- rownames(res_age_F)[which(res_age_F$padj < 0.05)]

venn_sets <- list(
  "Male-only" = male_age_hits,
  "Female-only" = female_age_hits
)

p_venn_pitfall <- ggVennDiagram(venn_sets, label_alpha = 0, set_size = 3.2) +
  scale_fill_gradient(low = "white", high = "#4575B4") +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_y_continuous(expand = expansion(mult = 0.10)) +
  labs(subtitle = "Subgroup age-hit lists") +
  theme(legend.position = "none")

# Genes that might be incorrectly called "female-specific age genes" from subgroup tests.
female_only_age_hits <- setdiff(female_age_hits, male_age_hits)

# Rank them by the female-only adjusted p-value and keep a readable number for the heatmap.
female_only_ranked <- female_only_age_hits[order(res_age_F[female_only_age_hits, "padj"])]
heatmap_gene_ids <- head(female_only_ranked, 30)

# If there are too few female-only hits, fall back to the top female age-associated genes
# so the plotting code still demonstrates the visual pitfall.
if (length(heatmap_gene_ids) < 5) {
  heatmap_gene_ids <- head(rownames(res_age_F)[order(res_age_F$padj)], 30)
}

heatmap_mat <- vsd_mat[heatmap_gene_ids, , drop = FALSE]
heatmap_mat <- heatmap_mat - rowMeans(heatmap_mat)

# Use gene symbols where available; make them unique for display.
heatmap_labels <- genes_final$gene_name[match(rownames(heatmap_mat), genes_final$gene_id)]
heatmap_labels[is.na(heatmap_labels) | heatmap_labels == ""] <- rownames(heatmap_mat)[is.na(heatmap_labels) | heatmap_labels == ""]
rownames(heatmap_mat) <- make.unique(heatmap_labels)

annotation_col <- data.frame(
  Sex = workshop_meta$sex,
  Age_group = age_group
)
rownames(annotation_col) <- colnames(heatmap_mat)

# Order samples so the eye sees the four groups clearly.
sample_order <- order(annotation_col$Sex, annotation_col$Age_group)

p_heatmap <- pheatmap(
  heatmap_mat[, sample_order, drop = FALSE],
  annotation_col = annotation_col[sample_order, , drop = FALSE],
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  # Symmetric breaks so white sits exactly at 0 (equal expression) and the red/blue
  # scale is centered — high and low deviations are shown on the same footing.
  breaks = seq(-max(abs(heatmap_mat)), max(abs(heatmap_mat)), length.out = 101),
  main = "Subgroup-selected age genes",
  silent = TRUE
)

# Draw the Venn diagram and heatmap side by side with patchwork. wrap_elements() lets
# us drop the pheatmap (a grid object, not a ggplot) straight into the layout.
print(p_venn_pitfall + wrap_elements(full = p_heatmap$gtable) +
      plot_layout(widths = c(1, 1.4)))

cat(sprintf(
  "Male-only age hits: %d | Female-only age hits: %d | Shared hits: %d\n",
  length(male_age_hits), length(female_age_hits), length(intersect(male_age_hits, female_age_hits))
))

#' 
#' The Venn diagram compares *lists*, not effect sizes. The heatmap is useful for visualization, but it does not prove a female-specific age effect. A gene outside the Venn overlap, or a gene that looks patterned in the heatmap, may simply have crossed an arbitrary FDR threshold in one smaller analysis but not the other. That selection rule is not an interaction test.
#' 
#' ::: {.callout-important}
#' ## Feuvrier et al. (2026) — three reasons the Venn overlap misleads, and a published example
#' 
#' The hippocampal-stress study dissects this exact figure and adds nuance every SABV analyst should internalize:
#' 
#' 1. **Absence of evidence is not evidence of absence.** Genes significant in only one sex usually show a *similar* response in the other sex that simply did not clear the threshold — visible directly in the per-sex log-fold-change heatmaps. Non-significance in one subgroup is not a sex-specific effect.
#' 2. **Imbalance fakes a sex difference.** When the dataset is imbalanced, the over-represented sex has more power and yields more DEGs, so it looks "more responsive" purely as an artifact. Feuvrier et al. show the asymmetry largely disappears when the dataset is re-balanced — yet genes still fall outside the overlap, because of point (1).
#' 3. **Even shared (overlap) genes can differ.** A gene significant in *both* sexes may still respond to very different *magnitudes* — which the Venn diagram cannot show and which only an interaction/contrast test quantifies.
#' 
#' **Published cautionary tale.** Marrocco et al. (2017, *Nat. Commun.*) reported sexually dimorphic stress transcription in CA3 pyramidal neurons using exactly this pattern — non-overlapping male/female DEG lists, with a very small pooled sample. A re-analysis (Ziegler et al. 2020) could not reproduce the sex specificity. The combination of a Venn-overlap interpretation and a small sample made it into a high-impact journal — which is why the interaction test, not the gene-list comparison, must be the standard.
#' :::
#' 
#' The next section introduces the missing statistical step: a formal **sex × age interaction** test.
#' 
#' ## MISTAKE 3: Comparing p-values between subgroups
#' 
#' Mistake 2 was **splitting** the data into two separate analyses. This mistake is the faulty *inference* people then draw from that split: "significant in males, not in females, therefore the effect differs by sex." The split gave us two p-values (`p_male`, `p_female` from above) — and comparing them is *not* a test. **A difference in significance is not a significant difference.** Only a formal interaction term, fit on both sexes in one model, actually asks whether the effect differs.
#' 
## -----------------------------------------------------------------------------
#| fig-width: 9
#| fig-height: 7
sprintf("Males p = %.4f (%s) | Females p = %.4f (%s)",
        p_male, ifelse(p_male < 0.05, "significant", "not significant"),
        p_female, ifelse(p_female < 0.05, "significant", "not significant"))

fit_interaction <- lm(expression ~ sex * condition, data = example_data)
interaction_p <- summary(fit_interaction)$coefficients["sexFemale:conditionOlder", 4]
sprintf("Interaction p-value: %.4f -> %s", interaction_p,
        ifelse(interaction_p < 0.05,
               "the age effect differs by sex",
               "no evidence the age effect differs by sex"))

# Visual comparison: the three wrong approaches vs. the factorial model.
# Build four small ggplots and arrange them 2x2 with patchwork.
yr <- range(example_data$expression)                 # shared y-axis across panels
panel_theme <- theme(legend.position = "none",
                     plot.title = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9))

g_pooled <- ggplot(example_data, aes(condition, expression)) +
  geom_boxplot(fill = "grey80") + coord_cartesian(ylim = yr) +
  labs(title = "WRONG: Pooled (sex ignored)",
       subtitle = sprintf("p = %.4f", summary(fit_pooled)$coefficients[2, 4]),
       x = "Age group", y = "Expression") + panel_theme

g_male <- ggplot(male_data, aes(condition, expression)) +
  geom_boxplot(fill = "#2166AC") + coord_cartesian(ylim = yr) +
  labs(title = "WRONG: Males only",
       subtitle = sprintf("p = %.4f (n=%d)", p_male, nrow(male_data)),
       x = "Age group", y = "Expression") + panel_theme

g_female <- ggplot(female_data, aes(condition, expression)) +
  geom_boxplot(fill = "#B2182B") + coord_cartesian(ylim = yr) +
  labs(title = "WRONG: Females only",
       subtitle = sprintf("p = %.4f (n=%d)", p_female, nrow(female_data)),
       x = "Age group", y = "Expression") + panel_theme

g_factorial <- ggplot(example_data, aes(condition, expression, fill = sex)) +
  geom_boxplot() + coord_cartesian(ylim = yr) +
  scale_fill_manual(values = c("Male" = "#2166AC", "Female" = "#B2182B")) +
  labs(title = "CORRECT: Factorial model",
       subtitle = sprintf("Interaction p = %.4f", interaction_p),
       x = "Age group", y = "Expression", fill = "Sex") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9))

(g_pooled | g_male) / (g_female | g_factorial)

#' 
#' For KL the per-sex p-values look very different — significant in females (≈0.03), not in males (≈0.37) — exactly the pattern that tempts a "female-specific" conclusion. But the interaction test is **non-significant** (≈0.15): the data give **no evidence** that KL's age effect differs between the sexes, so we cannot call it sex-specific. (Note the direction of the claim: non-significance does not *prove* the effects are equal — it just means splitting the data produced two different-looking p-values without the interaction test finding a real difference.) The four-panel figure shows why each "wrong" view misleads, while the factorial model asks the question directly.
#' 
#' **In short:** pooling inflates residual variance, disaggregation halves N and can't test interaction, and comparing subgroup p-values is invalid. The factorial model tests the main effects *and* the sex × age interaction in a single fit — the correct approach.
#' 
#' The same trap appears with **WFDC2**, only in the *opposite* direction — disaggregation makes the age effect look **male**-specific, yet the interaction test again finds no real sex difference:
#' 
## -----------------------------------------------------------------------------
zm  <- example_data2
p_m <- summary(lm(expression ~ condition, data = subset(zm, sex == "Male")))$coefficients[2, 4]
p_f <- summary(lm(expression ~ condition, data = subset(zm, sex == "Female")))$coefficients[2, 4]
p_i <- summary(lm(expression ~ sex * condition, data = zm))$coefficients["sexFemale:conditionOlder", 4]
sprintf("WFDC2 — males-only age p = %.4f (%s) | females-only age p = %.4f (%s) | interaction p = %.4f",
        p_m, ifelse(p_m < 0.05, "sig", "ns"),
        p_f, ifelse(p_f < 0.05, "sig", "ns"), p_i)

#' 
#' So whether a gene looks "female-specific" (KL) or "male-specific" (WFDC2), disaggregation manufactures the appearance of a sex difference that the interaction test does not support. **Only the interaction term answers the question.**
#' 
#' **The interaction test settles it — genome-wide.** Of all genes, how many *actually* have an age effect that differs by sex?
#' 
## -----------------------------------------------------------------------------
# The proper question, asked once on the full data (the factorial model from Exercise 3.4).
interaction_hits <- rownames(res_interaction)[which(res_interaction$padj < 0.05)]
n_int    <- length(interaction_hits)
n_tested <- sum(!is.na(res_interaction$padj))
cat(sprintf("Genes with a real sex x age interaction (FDR<0.05): %d of %d tested.\n", n_int, n_tested))
cat(sprintf("Compare this with the subgroup Venn diagram above: %d male-only age hits and %d female-only age hits came from separate analyses, but only %d genes pass the formal interaction test.\n",
            length(male_age_hits), length(female_age_hits), n_int))

#' 
#' So the hundreds of genes that looked "significant in one sex only" shrink to a **handful** once you test properly — the rest were an artifact of splitting the data in half, not real biology. The rule: to ask whether an effect differs between the sexes, **put both sexes in one model and test the interaction** — don't run and compare separate analyses.
#' 
#' ## Final Summary
#' 
#' Key messages, mapped to their sources:
#' 
#' - **Sample size is shared, not doubled**, when both sexes are included (Gaskill et al. 2025), and N is counted in **experimental units, not measurements** — technical replicates do not add N (Reynolds 2024, Ch. 2; avoid pseudo-replication).
#' - **Interactions cost power.** A design powered for a main effect is typically underpowered for the sex × treatment interaction (Reynolds 2024, Ch. 19; ≈4× rule of thumb). Screen broadly, then size a confirmatory study for the interaction ("right-sizing, not significance"; 3Rs).
#' - **Females are not substantially more variable**; apparent differences in "responsiveness" often come from **imbalance**, which gives the larger group more power (Feuvrier et al. 2026).
#' - **Pooling and disaggregation are not valid substitutes for modelling sex**, and a **difference in significance is not a significant difference** (Nieuwenhuis et al. 2011). The Venn-overlap of separate analyses is descriptive, not a test.
#' - **Parameterization matters:** the `condition` coefficient in an intercept model is the *baseline-sex* effect, not the overall effect; a **no-intercept model with explicit contrasts** (`~ 0 + group`) asks exactly the intended question (Feuvrier et al. 2026).
#' - **Factorial / interaction models** such as `~ sex + treatment + sex:treatment` (or the contrast equivalent) are the correct way to evaluate main effects *and* interactions in a single coherent fit.
