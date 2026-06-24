# =============================================================================
# Two-sample Mendelian randomization
#   Exposure : brain cell type-specific cis-caQTLs (PsychENCODE2;
#              NeuN+ neuronal / NeuN- non-neuronal open chromatin regions)
#   Outcome  : Migraine
#              - Discovery   : FinnGen Release 12 (G6_MIGRAINE)
#              - Replication : UK Biobank WGS migraine GWAS (GCST90473326)
#              - Replication : Million Veteran Program migraine GWAS (GCST90475837)
#              - Subtype     : MVP migraine with aura (GCST90477543) [exploratory]
#
# Design note:
#   Each open chromatin region (OCR) is instrumented by ONE or TWO cis-caQTL
#   SNPs. The causal estimate therefore uses the Wald ratio (single-SNP OCRs)
#   or inverse-variance weighted (IVW) method (two-SNP OCRs). With only 1-2
#   instruments per exposure, multi-instrument pleiotropy / sensitivity
#   estimators (MR-Egger, MR-PRESSO, weighted median, weighted mode,
#   leave-one-out) are NOT applicable and are deliberately not used; pleiotropy
#   is instead assessed at the instrument level (PheWAS) and through Bayesian
#   colocalization in separate scripts. Cochran's Q heterogeneity is reported
#   only for the two-SNP IVW OCRs.
#
# Pipeline:
#   1. Load caQTL instruments; select one OCR (cell type + peak)
#   2. Format exposure and outcome data
#   3. Harmonise
#   4. MR (Wald ratio / IVW) + OR
#   5. Cochran's Q heterogeneity (IVW / >=2 SNPs only)
#   6. Steiger directionality test
#   7. Export results
#
# To run across all OCRs and all outcome cohorts, wrap section
# "PER-OCR / PER-COHORT ANALYSIS" in a loop over your OCR list and the four
# outcome files (see the loop skeleton at the bottom).
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(TwoSampleMR)     # MR estimators, harmonisation, Steiger test
library(tidyverse)
library(data.table)
library(openxlsx)

# ---- User settings ----------------------------------------------------------
# >>> CHECK THESE PATHS AND COLUMN NAMES AGAINST YOUR OWN FILES <<<

# --- input files ---
path_caqtl  <- "data/psychencode2_brain_caQTL.txt"   # caQTL instruments (all OCRs)
out_root    <- "output/caQTL_migraine_result"        # results root directory

# Outcome summary-statistics files (one per cohort)
path_outcome <- list(
  FinnGen = "data/finngen_R12_G6_MIGRAINE.txt",       # discovery
  UKB     = "data/GCST90473326_UKB_migraine.txt",     # replication 1
  MVP     = "data/GCST90475837_MVP_migraine.txt",     # replication 2
  MVP_MA  = "data/GCST90477543_MVP_migraine_aura.txt" # exploratory subtype (with aura)
)
# Total sample size (cases + controls) per cohort, for Steiger test
outcome_samplesize <- list(
  FinnGen = 26894 + 374605,
  UKB     = 25393 + 433047,
  MVP     = 31836 + 405831,
  MVP_MA  =  4445 + 441492
)
# Number of cases per cohort (for Steiger binary-trait prevalence; optional)
outcome_ncase <- list(FinnGen = 26894, UKB = 25393, MVP = 31836, MVP_MA = 4445)

# --- target instrument (single OCR; change cell type + peak to analyse) ---
target_celltype <- "non-neuron"   # one of: neuron / non-neuron
target_peak     <- "Peak_7234"    # OCR (peak) identifier to analyse
which_outcome   <- "FinnGen"      # which cohort in path_outcome to run here
set.seed(5201314)

# --- caQTL (exposure) column names in path_caqtl ---
# The caQTL file is expected to contain, for every OCR, its cis instrument SNP(s).
caqtl_celltype_col <- "cell_type"      # "neuron" / "non-neuron"
caqtl_peak_col     <- "peak"           # OCR / peak identifier (e.g. Peak_7234)
caqtl_snp_col      <- "SNP"            # rsID of the cis-caQTL instrument
caqtl_beta_col     <- "beta"           # caQTL effect size (per unit normalized accessibility)
caqtl_se_col       <- "se"
caqtl_ea_col       <- "effect_allele"
caqtl_oa_col       <- "other_allele"
caqtl_eaf_col      <- "eaf"            # exposure allele frequency (for Steiger); set NULL if absent
caqtl_pval_col     <- "pval"
caqtl_samplesize   <- 1932            # PsychENCODE2 caQTL sample size (ATAC-seq libraries)

# --- outcome (migraine GWAS) column names ---
# Defaults follow a typical GWAS summary-statistics layout; edit to match each file.
out_snp_col   <- "SNP"
out_beta_col  <- "beta"
out_se_col    <- "se"
out_ea_col    <- "effect_allele"
out_oa_col    <- "other_allele"
out_eaf_col   <- "eaf"
out_pval_col  <- "pval"
out_chr_col   <- "chr"
out_pos_col   <- "pos"

# =============================================================================
# PER-OCR / PER-COHORT ANALYSIS
# =============================================================================
t1 <- Sys.time()

# ---- Load caQTL instruments and select the target OCR -----------------------
caqtl <- as.data.frame(fread(path_caqtl))
exp_sel <- caqtl[caqtl[[caqtl_celltype_col]] == target_celltype &
                 caqtl[[caqtl_peak_col]]     == target_peak, ]
if (nrow(exp_sel) == 0) stop("No instruments found for ", target_celltype, " ", target_peak)

# Exposure phenotype label, e.g. "non-neuron_Peak_7234"
exp_sel$exposure_label <- paste0(target_celltype, "_", target_peak)

# ---- Format exposure --------------------------------------------------------
exp_dat <- format_data(
  dat               = exp_sel,
  type              = "exposure",
  phenotype_col     = "exposure_label",
  snp_col           = caqtl_snp_col,
  beta_col          = caqtl_beta_col,
  se_col            = caqtl_se_col,
  effect_allele_col = caqtl_ea_col,
  other_allele_col  = caqtl_oa_col,
  eaf_col           = caqtl_eaf_col,
  pval_col          = caqtl_pval_col
)
exp_dat$samplesize.exposure <- caqtl_samplesize

# ---- Outcome ----------------------------------------------------------------
oc_file <- path_outcome[[which_outcome]]
outcome <- as.data.frame(fread(oc_file))
outcome$phenotype  <- paste0("MIGRAINE_", which_outcome)
outcome$samplesize <- outcome_samplesize[[which_outcome]]

out_dat <- format_data(
  dat               = outcome,
  type              = "outcome",
  snps              = exp_dat$SNP,
  phenotype_col     = "phenotype",
  snp_col           = out_snp_col,
  beta_col          = out_beta_col,
  se_col            = out_se_col,
  effect_allele_col = out_ea_col,
  other_allele_col  = out_oa_col,
  pval_col          = out_pval_col,
  eaf_col           = out_eaf_col,
  chr_col           = out_chr_col,
  pos_col           = out_pos_col,
  samplesize_col    = "samplesize"
)
out_dat <- out_dat %>% subset(., !duplicated(SNP))

# ---- Harmonise --------------------------------------------------------------
mydata <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
mydata <- mydata[which(mydata$mr_keep == TRUE), ]
n_iv   <- nrow(mydata)
if (n_iv == 0) stop("No instruments remain after harmonisation for ",
                    target_celltype, " ", target_peak, " / ", which_outcome)

# ---- MR: Wald ratio (1 SNP) or IVW (>=2 SNPs) -------------------------------
# Method is chosen automatically by the number of instruments retained.
if (n_iv == 1) {
  method_list <- "mr_wald_ratio"
} else {
  method_list <- "mr_ivw"   # inverse-variance weighted for two-SNP OCRs
}
res    <- mr(mydata, method_list = method_list)
res_or <- generate_odds_ratios(res)   # adds OR and 95% CI

# ---- Cochran's Q heterogeneity (only meaningful with >=2 instruments) -------
if (n_iv >= 2) {
  het <- mr_heterogeneity(mydata, method_list = "mr_ivw")
} else {
  het <- data.frame(method = "Wald ratio",
                    Q = NA, Q_df = NA, Q_pval = NA,
                    note = "single instrument: heterogeneity not estimable")
}

# ---- Steiger directionality test --------------------------------------------
# Confirms the exposure explains more variance than the outcome (correct
# exposure-to-outcome direction). Requires allele frequencies and sample sizes.
mydata$samplesize.exposure <- caqtl_samplesize
mydata$samplesize.outcome  <- outcome_samplesize[[which_outcome]]
# If outcome is binary and you have ncase/ncontrol + prevalence, you may instead
# use get_r_from_lor() to derive r.outcome; here we use the default directionality test.
steiger <- tryCatch(
  directionality_test(mydata),
  error = function(e) data.frame(snp_r2.exposure = NA, snp_r2.outcome = NA,
                                 correct_causal_direction = NA, steiger_pval = NA)
)

# ---- Export -----------------------------------------------------------------
out_dir <- file.path(out_root, paste0(target_celltype, "_", target_peak, "_", which_outcome))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.table(mydata,  file.path(out_dir, "harmonised_instruments.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)   # instruments used
write.table(res,     file.path(out_dir, "mr_result.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)   # MR estimate (b, se, pval)
write.table(res_or,  file.path(out_dir, "mr_result_OR.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)   # MR estimate as OR + 95% CI
write.table(het,     file.path(out_dir, "heterogeneity.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)   # Cochran's Q (IVW / >=2 SNP only)
write.table(steiger, file.path(out_dir, "steiger.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)   # Steiger directionality

t2 <- Sys.time()
message("MR for ", target_celltype, " ", target_peak, " -> MIGRAINE_", which_outcome,
        " (", n_iv, " IV", ifelse(n_iv > 1, "s", ""), ", ",
        ifelse(n_iv == 1, "Wald ratio", "IVW"), ") done in ",
        round(difftime(t2, t1, units = "mins"), 3), " minutes.")

# =============================================================================
# OPTIONAL: loop over all OCRs x all cohorts
# -----------------------------------------------------------------------------
# Uncomment and adapt to reproduce the full discovery + replication + subtype
# analyses. `ocr_list` should enumerate every (cell_type, peak) pair tested.
#
# ocr_list <- unique(caqtl[, c(caqtl_celltype_col, caqtl_peak_col)])
# all_res <- list()
# for (ci in seq_len(nrow(ocr_list))) {
#   ct  <- ocr_list[[caqtl_celltype_col]][ci]
#   pk  <- ocr_list[[caqtl_peak_col]][ci]
#   for (co in names(path_outcome)) {
#     ## ... set target_celltype <- ct; target_peak <- pk; which_outcome <- co
#     ## ... rerun the PER-OCR / PER-COHORT ANALYSIS block above
#     ## ... collect res_or / steiger into all_res
#   }
# }
# combined <- dplyr::bind_rows(all_res)
# write.xlsx(combined, file.path(out_root, "all_OCR_all_cohorts_MR.xlsx"))
# =============================================================================
