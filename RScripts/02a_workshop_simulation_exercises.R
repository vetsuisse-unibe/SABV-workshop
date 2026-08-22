#' ---
#' title: "Hands-on Exercise: Experimental Design & Power Simulation"
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
#' > Design-phase tutorial: simulate cohorts to see how including both sexes, cohort size, and the main-effect-vs-interaction distinction change your statistical power — *before* any real data are collected.
#' 
#' # Overview
#' 
#' Script 02a: Experimental Design — Sex as Biological Variable in Genomics.
#' 
#' This tutorial is purely a **simulation** exercise. No sequencing data are needed: we generate synthetic cohorts, run simple linear models, and repeat the experiment thousands of times to measure statistical **power** directly. The next tutorial (`02b_workshop_gtex_liver_exploration`) then applies these ideas to real GTEx liver data.
#' 
#' ## Learning objectives
#' 
#' By the end of this tutorial, you will be able to:
#' 
#' 1. Define statistical power as the probability of detecting a genuine effect.
#' 2. Differentiate total sample size from sample size per group.
#' 3. Compare statistical power when detecting a main effect versus an interaction.
#' 4. Understand why introducing both sexes does not necessarily require doubling your animal count.
#' 
#' Key terms used in this module:
#' 
#' - **Effect size:** how large the biological signal is relative to background noise.
#' - **Power:** how consistently a study detects a true effect across repeated trials.
#' - **Interaction:** when a treatment outcome depends on another variable — such as a drug affecting males and females differently.
#' 
#' ::: {.callout-note}
#' ## Primary literature for this module
#' 
#' - **Reynolds (2024)**, *A Guide to Sample Size for Animal-based Studies* (Wiley). Source for: the experimental vs biological unit, pseudo-replication, effect sizes, the lower power of interaction terms relative to main effects, factorial/screening efficiency, and "right-sizing not significance" (3Rs).
#' - **Gaskill et al. (2025)**, *eLife* — the survey of researcher *barriers and misconceptions* about studying sex as a variable; source of the percentage claims (e.g. "including both sexes doubles the sample size") the simulations are designed to correct.
#' - **Karp et al. (2025)**, *Nature Communications* — the **SIRF** framework for analysing sex as a biological variable (the framework name belongs to this paper, not to the Gaskill survey).
#' :::
#' 
#' ## Load Packages
#' 
## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(ggplot2)
library(dplyr)
library(tidyr)

# Set a consistent theme
theme_set(theme_bw(base_size = 14))

#' 
#' # EXPERIMENTAL DESIGN - Power Simulation Exercise
#' 
#' ## EXERCISE 2.1: Compare designs — n=50M/50F vs n=100M
#' 
## -----------------------------------------------------------------------------
#| fig-alt:
#|   - "Boxplot of simulated liver function in control and drug-treated mice. The group distributions overlap, while a red X marks each group mean."
#|   - "Grouped bar chart comparing statistical power for drug and sex effects in a male-only design, a same-size two-sex design, and a larger two-sex design. A dashed line marks 80 percent power."
# KEY MISCONCEPTION ADDRESSED: "Including both sexes doubles sample size"
# Reality: N is SHARED between sexes in exploratory inclusion
# (Gaskill et al. 2025: 80% of researchers held this misconception)
# STEP 1: Look at ONE experiment first, to understand what's happening.

# --- Phase 1: Experimental Design Setup ---
# Define cohort size: N = 100
n_mice <- 100
# Create the treatment vector:
# 0 = Control group (n = 50 mice)
# 1 = Experimental drug-treatment group (n = 50 mice)
treatment <- c(rep(0, 50), rep(1, 50))

# --- Phase 2: Define the biological signal and the noise ---
# Define the biological effect of the drug. On average, treated mice experience
# a 0.5-unit improvement in liver function.
true_drug_effect <- 0.5
# Natural phenotypic variation: inter-individual variability between mice due to
# genetics, environment, or handling.
mouse_sd <- 1.0
# Standardized effect size (Cohen's d) = the signal-to-noise ratio of the treatment effect.
# A value of 0.5 ("medium") means the drug effect is half the size of natural individual
# variation — it works at the group level, but treated and untreated mice still look
# alike one-to-one.
sprintf("Effect size (Cohen's d) = %.2f", true_drug_effect / mouse_sd)

# --- Phase 3: Simulate the data ---
set.seed(42)  # reproducible noise
# Outcome = drug effect (in treated mice) + individual biological variability
liver_function <- true_drug_effect * treatment + rnorm(n_mice, mean = 0, sd = mouse_sd)

# --- Phase 4: Data Wrangling ---
# Create a dataframe and assign factor labels for analysis and plotting.
experiment_data <- data.frame(
  mouse_id = 1:n_mice,
  treatment = factor(treatment, levels = c(0, 1), labels = c("Control", "Drug")),
  liver_function = liver_function
)

# --- Phase 5: Descriptive Statistics ---
# Check group means against the true effect size.
means <- experiment_data %>%
  group_by(treatment) %>%
  summarise(mean_lf = mean(liver_function), .groups = "drop")
sprintf("Control: %.2f | Drug: %.2f | Diff: %.2f",
        means$mean_lf[1], means$mean_lf[2], diff(means$mean_lf))

# Boxplot of the single experiment, with group means marked
p_single <- ggplot(experiment_data, aes(x = treatment, y = liver_function, fill = treatment)) +
  geom_boxplot(alpha = 0.6, outlier.size = 1) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 2) +
  # Automatically calculates and plots the mean as a red 'X'
  stat_summary(fun = mean, geom = "point", color = "red", size = 4, shape = 4) +
  scale_fill_manual(values = c("Control" = "lightblue", "Drug" = "lightcoral")) +
  labs(
    title = "Drug Effect on Liver Function",
    subtitle = "Red X = group mean",
    y = "Liver Function (units)", x = ""
  ) +
  theme_minimal() +               # cleans up the background grid
  theme(legend.position = "none")

print(p_single)

# Run a linear regression to test for differences between treatment groups.
fit_single <- lm(liver_function ~ treatment, data = experiment_data)
# p-value for the treatment effect
p_value_single <- summary(fit_single)$coefficients[2, 4]
sprintf("Experimental p-value: %.4f", p_value_single)

# STEP 2: Simulate the same experiment 1,000 times to determine power. Power is the
# probability of detecting a true biological effect as statistically significant (p < 0.05).
run_one_experiment <- function(n_control = 50, n_drug = 50, true_effect = 0.5) {
  treatment <- c(rep(0, n_control), rep(1, n_drug))
  outcome   <- true_effect * treatment + rnorm(n_control + n_drug, 0, 1)
  fit       <- lm(outcome ~ treatment)
  p_value   <- summary(fit)$coefficients[2, 4]
  detected  <- p_value < 0.05
  return(detected)
}

# Run the experiment 1000 times
results_A <- replicate(1000, run_one_experiment(n_control = 50, n_drug = 50, true_effect = 0.5))
power_A <- mean(results_A)
sprintf("Detected drug effect in %d / 1000 experiments | POWER = %.1f%%",
        sum(results_A), power_A * 100)

# STEP 3: Compare three study designs.
# The function runs many simulated experiments using cohorts of both sexes and
# counts how often it picks up two signals:
#   - the drug effect — does the treatment change liver function?
#   - the sex difference — do males and females have different baseline liver function?
# The three designs:
#   A: 100 males only          (N = 100) — can test the drug, but not sex
#   B: 50 males + 50 females   (N = 100) — same total N, split across sexes
#   C: 100 males + 100 females (N = 200) — both sexes, double the sample size
# Comparing them shows how cohort size and including both sexes change our
# ability to detect each signal.
simulate_power_design <- function(n_male, n_female, true_treatment_effect,
                                   true_sex_effect = 0, n_sims = 1000) {

  # Track, per simulation, whether each main effect was significant (p < 0.05).
  detection <- data.frame(treatment = logical(n_sims), sex = logical(n_sims))

  for (i in 1:n_sims) {
    # Split by sex and randomize treatment within each sex (balanced arms).
    sex <- c(rep(0, n_male), rep(1, n_female))
    treatment <- c(sample(rep(0:1, c(n_male / 2, n_male / 2))),
                   sample(rep(0:1, c(n_female / 2, n_female / 2))))

    # Simulate data: treatment effect + sex difference + noise (no true interaction).
    outcome <- true_treatment_effect * treatment +
               true_sex_effect * sex +
               rnorm(n_male + n_female, 0, 1)

    # Factorial model; each main effect is estimated adjusted for the other.
    pvals <- summary(lm(outcome ~ treatment + sex))$coefficients[, 4]
    detection$treatment[i] <- pvals["treatment"] < 0.05
    detection$sex[i]       <- pvals["sex"] < 0.05
  }

  colMeans(detection)   # power to detect each main effect
}

# Design A — 100 males only: 50 control + 50 drug.
set.seed(42)
power_A_simple <- mean(replicate(1000, {
  y <- 0.5 * c(rep(0, 50), rep(1, 50)) + rnorm(100)
  summary(lm(y ~ c(rep(0, 50), rep(1, 50))))$coefficients[2, 4] < 0.05
}))
sprintf("Design A — drug-effect power: %.1f%%", power_A_simple * 100)

# Designs B and C — both sexes, analysed with the factorial model lm(outcome ~ treatment + sex).
# Reuse simulate_power_design() (true sex effect = 0.3 SD, no true interaction):
#   B: 50 M + 50 F   (N = 100, same total as A)
#   C: 100 M + 100 F (N = 200, confirmatory)
power_B <- simulate_power_design(n_male = 50,  n_female = 50,
                                 true_treatment_effect = 0.5, true_sex_effect = 0.3)
power_C <- simulate_power_design(n_male = 100, n_female = 100,
                                 true_treatment_effect = 0.5, true_sex_effect = 0.3)

# Pull out the scalars used by the plot and summaries below.
power_B_treatment <- power_B["treatment"]; power_B_sex <- power_B["sex"]
power_C_treatment <- power_C["treatment"]; power_C_sex <- power_C["sex"]

sprintf("Design B (N=100) — drug %.1f%%, sex %.1f%%", power_B_treatment * 100, power_B_sex * 100)
sprintf("Design C (N=200) — drug %.1f%%, sex %.1f%%", power_C_treatment * 100, power_C_sex * 100)

# Visualization comparing designs
comparison_data <- data.frame(
  Design = c("A: Males only\n(n=100)",
             "B: Both sexes\n(n=100 shared)",
             "C: Both sexes\n(n=200)"),
  Drug_Effect = c(power_A_simple * 100, power_B_treatment * 100, power_C_treatment * 100),
  Sex_Difference = c(NA, power_B_sex * 100, power_C_sex * 100)
)

plot_data <- comparison_data %>%
  pivot_longer(cols = -Design, names_to = "Outcome", values_to = "Power") %>%
  filter(!is.na(Power))

p_power <- ggplot(plot_data, aes(x = Design, y = Power, fill = Outcome)) +
  geom_col(position = "dodge", alpha = 0.7) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("Drug_Effect" = "#2166AC", "Sex_Difference" = "#B2182B")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Comparing Study Designs: Power to Detect Effects",
    subtitle = "Red dashed line = 80% power (standard in biology)",
    y = "Statistical Power (%)",
    x = "",
    fill = ""
  ) +
  theme(axis.text = element_text(size = 11),
        legend.position = "bottom")

print(p_power)

sprintf(paste0(
  "Design B (N=100 shared): treatment power %.1f%%, sex-difference power %.1f%%. ",
  "Design C (N=200): treatment power %.1f%%, sex-difference power %.1f%%."),
  power_B_treatment * 100, power_B_sex * 100, power_C_treatment * 100, power_C_sex * 100)

#' 
#' ::: {.callout-note collapse="true"}
#' ## Interpretation — try to reason it out yourself first, then expand
#' 
#' Before opening this: look at the three bars. Which design reaches 80% power for the drug effect? Should Designs A and B differ on the *drug* effect at all? Why is the sex-difference bar so low?
#' 
#' - **Design A** (males only, n = 100): reaches only 70% power — below the 80% threshold — and by design cannot detect sex differences or flag whether the drug works differently in females.
#' - **Design B** (both sexes, n = 100 shared): the drug-effect power (~72%) is **essentially identical to Design A** — and it *should* be, since the two designs have the same total N and the same treatment effect, and putting sex in the model removes its variance rather than adding noise. The small A-vs-B gap is Monte-Carlo (sampling) error from running finitely many simulations, not a real benefit of adding a second sex; with a fixed seed and enough replicates it would disappear. The point is the *opposite* of a cost: including both sexes does **not** dilute the drug-effect analysis, and it hands you the sex-difference test for free. The low power for that sex effect (32%) simply reflects that the baseline male-to-female difference here is small.
#' - **Design C** (both sexes, n = 200): a larger cohort shows that doubling the sample size predictably amplifies both signals (pushing drug-effect power to 94% and sex-difference power to 54%).
#' :::
#' 
#' ## BONUS: What if the sex effect was larger?
#' 
#' What if sex differences were as large as the drug effect — say 0.5 SD of baseline biology instead of 0.3? Re-run designs B and C and compare the sex-difference power.
#' 
## -----------------------------------------------------------------------------
# Re-run designs B and C with a larger true sex effect (0.5 SD instead of 0.3),
# using the same factorial simulator.
set.seed(42)
power_B_alt <- simulate_power_design(n_male = 50,  n_female = 50,
                                     true_treatment_effect = 0.5, true_sex_effect = 0.5)
power_C_alt <- simulate_power_design(n_male = 100, n_female = 100,
                                     true_treatment_effect = 0.5, true_sex_effect = 0.5)

power_B_alt_treatment <- power_B_alt["treatment"]; power_B_alt_sex <- power_B_alt["sex"]
power_C_alt_treatment <- power_C_alt["treatment"]; power_C_alt_sex <- power_C_alt["sex"]

sprintf(paste0(
  "Sex effect 0.3 SD — sex-difference power: Design B %.1f%%, Design C %.1f%%. ",
  "Sex effect 0.5 SD — sex-difference power: Design B %.1f%%, Design C %.1f%%."),
  power_B_sex * 100, power_C_sex * 100,
  power_B_alt_sex * 100, power_C_alt_sex * 100)

#' 
#' **Takeaway.**
#' 
#' 1. **Subtle sex differences (0.3 SD)** — comparable to minor baseline weight shifts. Power stays low (~54%) even at n = 200. This is an inherent trait of small biological effect sizes, not a failure of the dual-sex design.
#' 2. **Strong sex differences (0.5 SD)** — equal to the drug effect itself. Power reaches ~69% at n = 100 and ~94% at n = 200. Prominent dimorphism is easily captured using standard, well-powered study layouts.
#' 
#' **Practical implication:** include both sexes in the initial experiment to screen for sex differences at low cost. If a sex-specific signal emerges, estimate its effect size and use that to power a focused follow-up study — rather than doubling cohort sizes upfront, which wastes animals and adds no statistical benefit. This also keeps the study aligned with the **3Rs** principles, and with Reynolds's "**right-sizing, not significance**" and **sequential/screening** logic (Reynolds 2024, §3.10–3.11): screen broadly, then run a definitive, well-powered confirmatory study only on the signals that survive.
#' 
#' ## EXERCISE 2.1b: The interaction costs more power than a main effect
#' 
#' In Exercise 2.1, Designs A–C powered the *main effects* — the overall treatment effect and baseline sex differences. However, the SABV question that usually matters most is the **interaction**: does the drug work *differently* in males and females? As highlighted by Reynolds (2024, Ch. 19, citing Jones & Nachtsheim 2011), statistical power in factorial designs is highest for main effects and substantially lower for interactions. As a rule of thumb, detecting an interaction of a given magnitude requires roughly **four times** the total sample size needed to detect a main effect of that same magnitude.
#' 
#' A common point of confusion is that these two questions measure fundamentally different quantities. Suppose a drug's true effect size is **+0.25 SD in males** and **+0.75 SD in females**:
#' 
#' - **Main effect (the average):** (0.25 + 0.75) / 2 = **0.50**. Answers: *"Does the treatment work on average across both sexes?"*
#' - **Interaction (the gap):** 0.75 − 0.25 = **0.50**. Answers: *"Does the treatment response differ by sex?"*
#' 
#' Both equal 0.50 in this simulation on purpose, so we can compare their statistical power at the exact same effect-size magnitude (0.5 SD). We evaluate each question using standard **0/1 (dummy) coding**:
#' 
#' 1. **Main effect ($Q_1$):** fit the additive model `lm(y ~ treatment + sex)` and test `treatment`. Because there is no interaction term in this model, R estimates a single treatment effect averaged across sexes.
#' 2. **Interaction ($Q_2$):** fit the full model `lm(y ~ treatment * sex)` and test the `treatment:sex` coefficient, which directly isolates the male–female response gap.
#' 
## -----------------------------------------------------------------------------
#| fig-alt: "Line chart of statistical power against total sample size. The treatment main-effect line remains above the sex-by-treatment interaction line at every sample size."
# The two sex-specific drug effects that everything below is built from:
effect_male   <- 0.25   # drug effect in males
effect_female <- 0.75   # drug effect in females

# Main effect = the AVERAGE of the two ; Interaction = the DIFFERENCE (gap) between them.
cat(sprintf("Main effect (average) = %.2f  |  Interaction (difference) = %.2f\n",
            (effect_male + effect_female) / 2, effect_female - effect_male))

# We compare power for the main effect and the interaction at the SAME size (0.5 SD),
# each asked with the model that matches the question -- all in 0/1 dummy coding.
# We never read a "main effect" out of a model that also contains an interaction term,
# so there is no reference-coding pitfall (contrast Exercise 3.5 in the next tutorial).

# Q1: does the treatment work ON AVERAGE?  -> additive model, test `treatment`.
power_main_effect <- function(n_per_cell, n_sims = 1000) {
  avg_effect <- (effect_male + effect_female) / 2      # 0.50 -- the average effect
  mean(replicate(n_sims, {
    sex       <- rep(c(0, 1), each = 2 * n_per_cell)     # 0 = Male, 1 = Female
    treatment <- rep(rep(c(0, 1), each = n_per_cell), 2) # 0 = Control, 1 = Treated
    outcome   <- treatment * avg_effect + rnorm(4 * n_per_cell)  # SAME effect in both sexes
    summary(lm(outcome ~ treatment + sex))$coefficients["treatment", "Pr(>|t|)"] < 0.05
  }))
}

# Q2: does the treatment effect DIFFER by sex?  -> interaction model, test `treatment:sex`.
power_interaction <- function(n_per_cell, n_sims = 1000) {
  mean(replicate(n_sims, {
    sex       <- rep(c(0, 1), each = 2 * n_per_cell)     # 0 = Male, 1 = Female
    treatment <- rep(rep(c(0, 1), each = n_per_cell), 2) # 0 = Control, 1 = Treated
    # If treated & male -> 0.25 | if treated & female -> 0.75
    trt_effect <- ifelse(sex == 0, effect_male, effect_female)
    outcome    <- treatment * trt_effect + rnorm(4 * n_per_cell)
    summary(lm(outcome ~ treatment * sex))$coefficients["treatment:sex", "Pr(>|t|)"] < 0.05
  }))
}

set.seed(42)
int_tab <- data.frame(n_per_cell = c(25, 50, 100, 200))
int_tab$total_N     <- int_tab$n_per_cell * 4
int_tab$treatment   <- sapply(int_tab$n_per_cell, power_main_effect) * 100
int_tab$interaction <- sapply(int_tab$n_per_cell, power_interaction) * 100
print(int_tab)

# Two power curves across total N.
int_long <- int_tab %>%
  pivot_longer(c(treatment, interaction), names_to = "Effect", values_to = "Power")

p_int_power <- ggplot(int_long, aes(x = total_N, y = Power, color = Effect, group = Effect)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("treatment" = "#2166AC", "interaction" = "#B2182B"),
                     labels = c("treatment" = "Treatment main effect",
                                "interaction" = "Sex x treatment interaction")) +
  labs(title = "Power Hierarchy: Main Effect vs Interaction",
       subtitle = "Same true magnitude (0.5 SD) for both; the interaction needs far more animals",
       x = "Total N (balanced 2x2)", y = "Statistical Power (%)", color = "") +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position = "bottom")
print(p_int_power)

#' 
#' As shown in the plot, the interaction curve sits well below the main-effect curve at every sample size. A cohort size that achieves >90% power to detect an overall drug effect (`total_N = 200`) yields only ~40% power to confirm whether that effect varies by sex.
#' 
#' **Practical implication:** a single exploratory experiment can *screen* for sex differences, but *confirming* an interaction requires a study explicitly powered for the interaction term — not just the main effect.
#' 
#' **A note on context.** The 4× sample-size rule represents an equal-magnitude comparison (0.5 SD vs 0.5 SD). In biological translation, the interactions that matter most are often larger — such as a drug working strongly in one sex but showing zero effect (or an opposite effect) in the other. Qualitative interactions like these are easier to detect, so treat the 4× benchmark as an *upper bound* on sample-size requirements rather than a fixed rule.
#' 
#' ::: {.callout-tip}
#' ## Where to go next
#' 
#' The next tutorial **`02b_workshop_gtex_liver_exploration`** applies these design principles to real GTEx liver RNA-seq data: testing the "females are more variable" myth, running PCA and differential expression, fitting the factorial sex × age model, and demonstrating the common analysis mistakes (pooling, disaggregation, comparing subgroup p-values).
#' :::
