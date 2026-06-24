# =============================================================================
# Bayesian colocalization: brain cis-caQTL vs migraine GWAS
#   Method   : coloc.abf (R package coloc), +/- 100 kb window, default priors
#   Exposure : brain cell type-specific cis-caQTLs (PsychENCODE2, non-neuronal)
#   Outcome  : Migraine, FinnGen Release 12 discovery GWAS (G6_MIGRAINE)
#   Tested   : the 7 OCRs Bonferroni-replicated in BOTH UK Biobank + MVP
#   Cutoff   : PP.H4 > 0.80 = strong evidence of a shared causal variant
#
# Adapted from the lab colocalization pipeline. Key changes vs the previous
# (sc-eQTL / NGF) version of this script:
#   - OUTCOME is now the FinnGen R12 migraine GWAS (case-control), not sc-eQTL.
#   - dataset1 = migraine GWAS (type "cc", with case fraction s);
#     dataset2 = caQTL (type "quant"). Matches caQTL = exposure, migraine = outcome.
#   - Fixed the coloc.abf call so dataset2's MAF is INSIDE the dataset2 list().
#   - peak_list now enumerates the 7 Bonferroni-replicated OCRs.
# =============================================================================

library(coloc)
library(locuscomparer)
library(data.table)
library(openxlsx)
library(tidyverse)
library(TwoSampleMR)
library(plinkbinr)

# ---- User settings ----------------------------------------------------------
# >>> CHECK PATHS / COLUMN NAMES AGAINST YOUR OWN FILES <<<

# caQTL full summary stats (ALL cis-SNPs per OCR) -- non-neuronal here
path_caqtl    <- "H:/caQTL_Brain/gila_caQTL.txt.gz"
# Migraine GWAS (FinnGen R12 discovery) -- REPLACES the old sc-eQTL outcome file
path_migraine <- "H:/Migraine_GWAS/finngen_R12_G6_MIGRAINE.txt.gz"
# LD reference panel for clumping
path_bfile    <- "D:/R_ref_data/MR_ref/EUR"

# The 7 OCRs Bonferroni-replicated in both UKB + MVP (edit IDs to match your file)
peak_list <- c("Peak_7234", "Peak_8343", "Peak_124742",
               "Peak_163823", "Peak_202875", "Peak_65776", "Peak_261762")

# Migraine (outcome) is case-control: need total N and case fraction s
migraine_N     <- 26894 + 374605   # FinnGen R12 total (cases + controls)
migraine_ncase <- 26894
migraine_s     <- migraine_ncase / migraine_N

# caQTL (exposure) sample size
caqtl_N <- 1932

set.seed(5201314)

# ---- Read caQTL full data ---------------------------------------------------
data2 <- fread(path_caqtl)
cat("caQTL data loaded, rows:", nrow(data2), "\n")

# ---- Read & process migraine GWAS outcome (shared by all peaks) -------------
mig <- fread(path_migraine)
head(mig)

# >>> EDIT these column names to match your FinnGen file <<<
# Need: SNP, chrom, pos, effect/other allele, MAF, beta, se, P
mig <- mig %>%
  dplyr::select(SNP, chrom, pos, effect_allele, other_allele, maf, beta, se, pval)
colnames(mig) <- c('SNP', 'chrom', 'pos', 'effect_allele',
                   'other_allele', 'MAF', 'beta', 'se', 'P')

# MAF lookup (used to supply MAF to the caQTL dataset if it lacks one)
mafdata <- mig[, c("SNP", "MAF")]

mig$samplesize <- migraine_N
mig$varbeta    <- mig$se^2
mig$z          <- mig$beta / mig$se

# drop P == 0 to avoid coloc errors
GWASdata0 <- mig %>% na.omit() %>% filter(P > 0)
cat("Migraine outcome ready, SNPs:", nrow(GWASdata0), "\n\n")

### ============================================
### Batch loop: each OCR runs the full pipeline
### ============================================
for (peak in peak_list) {

  cat("========================================\n")
  cat("Processing:", peak, "\n")
  cat("========================================\n")

  out_dir <- paste0("./", peak)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # ---- select this OCR from the full caQTL data ----
  # If the OCR column "OCR" holds values like "Peak_7234", this works as is.
  # If values are "gila_Peak_7234", use: data2[OCR == paste0("gila_", peak)]
  data_peak <- data2[OCR == peak]
  if (nrow(data_peak) == 0) {
    cat("WARNING:", peak, "not found in caQTL data; skipped.\n\n"); next
  }
  cat(peak, "rows after selection:", nrow(data_peak), "\n")

  # ---- significant SNPs + LD clump to define lead SNP(s) ----
  gwas1 <- data_peak[data_peak$P < 5e-8, ]
  if (nrow(gwas1) == 0) {
    cat("WARNING:", peak, "has no P < 5e-8 SNP; skipped.\n\n"); next
  }
  gwas1 <- as.data.frame(gwas1)

  # >>> EDIT caQTL column names (CHR/BP/A1/A2/BETA/SE/P) to match your file <<<
  gwas1$N <- caqtl_N
  exp_dat <- format_data(
    dat = gwas1, type = "exposure",
    chr_col = "CHR", pos_col = "BP", snp_col = "SNP",
    effect_allele_col = "A1", other_allele_col = "A2",
    beta_col = "BETA", se_col = "SE", pval_col = "P",
    samplesize_col = "N", phenotype_col = "phenotype"
  )

  exp_dat_clump <- ld_clump(
    clump_kb = 10000, clump_r2 = 0.01, pop = "EUR",
    dplyr::tibble(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure, id = exp_dat$id.exposure),
    plink_bin = get_plink_exe(), bfile = path_bfile
  )
  exp_dat_clump <- subset(exp_dat, SNP %in% exp_dat_clump$rsid)
  gwas2 <- exp_dat_clump %>% arrange(pval.exposure)

  # ---- write the +/- 100 kb region(s) ----
  kb100 <- gwas2[, c("chr.exposure", "SNP", "pos.exposure")]
  kb100$start <- kb100$pos.exposure - 100000
  kb100$end   <- kb100$pos.exposure + 100000
  write.csv(kb100, paste0(out_dir, "/", peak, "_kb100.csv"), row.names = FALSE)
  cat(peak, "lead SNP count:", nrow(kb100), "\n")

  # ---- prepare this OCR's cis-caQTL data (ALL region SNPs, NOT clumped) ----
  data_qtl <- data_peak
  data_qtl$N <- caqtl_N
  # supply MAF from the GWAS lookup if the caQTL file lacks a usable MAF column
  data_qtl <- merge.data.table(data_qtl, mafdata, by = "SNP")
  # >>> EDIT caQTL column names below to match your file <<<
  data_qtl <- data_qtl %>%
    dplyr::select(SNP, CHR, BP, A1, A2, MAF, BETA, SE, P, N)
  colnames(data_qtl) <- c('SNP', 'chrom', 'pos', 'effect_allele',
                          'other_allele', 'MAF', 'beta', 'se', 'P', 'samplesize')
  data_qtl$varbeta <- data_qtl$se^2
  data_qtl$z       <- data_qtl$beta / data_qtl$se

  # ---- coloc per lead SNP ----
  kb100res <- list(); final_snp_res1 <- c()
  pph0_4_all <- c(); pph0_4_snp_all <- c(); pph3_4 <- c()

  for (i in 1:nrow(kb100)) {

    leadchr <- as.numeric(kb100$chr.exposure[i])

    QTLdata <- data_qtl[data_qtl$chrom == leadchr, ]
    QTLdata <- QTLdata[QTLdata$pos > kb100$start[i] & QTLdata$pos < kb100$end[i], ]
    QTLdata <- subset(QTLdata, !duplicated(SNP)) %>% na.omit()

    sameSNP <- intersect(QTLdata$SNP, GWASdata0$SNP)
    if (length(sameSNP) == 0) {
      cat("  lead SNP", kb100$SNP[i], "no shared SNP; skipped\n"); next
    }

    QTLdata  <- QTLdata[QTLdata$SNP %in% sameSNP, ] %>% arrange(SNP) %>% na.omit()
    QTLdata[QTLdata$P == 0, ]$P <- 1e-300
    GWASdata <- GWASdata0[GWASdata0$SNP %in% sameSNP, ] %>% arrange(SNP) %>% na.omit()

    # ---- coloc.abf : dataset1 = migraine (cc), dataset2 = caQTL (quant) ----
    # default priors (p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    res1 <- coloc.abf(
      dataset1 = list(pvalues = GWASdata$P, snp = GWASdata$SNP,
                      type = "cc", s = migraine_s,
                      N = GWASdata$samplesize[1], MAF = GWASdata$MAF),
      dataset2 = list(pvalues = QTLdata$P, snp = QTLdata$SNP,
                      type = "quant",
                      N = QTLdata$samplesize[1], MAF = QTLdata$MAF)
    )

    # ---- locuscompare plot ----
    gwas_fn <- GWASdata[, c('SNP', 'P')] %>% rename(rsid = SNP, pval = P)
    pqtl_fn <- QTLdata[,  c('SNP', 'P')] %>% rename(rsid = SNP, pval = P)
    pdf(paste0(out_dir, "/", kb100$SNP[i], "_100kb_coloc.pdf"))
    print(locuscompare(in_fn1 = gwas_fn, in_fn2 = pqtl_fn,
                       title1 = 'Migraine GWAS', title2 = peak))
    dev.off()

    # ---- tidy results ----
    a <- as.data.frame(res1[1]) %>%
      tibble::rownames_to_column("PP0_4") %>%
      pivot_wider(names_from = "PP0_4", values_from = "summary")
    a$chr <- kb100$chr.exposure[i]; a$pos <- kb100$pos.exposure[i]; a$lead_snp <- kb100$SNP[i]

    b <- as.data.frame(res1[3]) %>%
      tibble::rownames_to_column("pp") %>%
      pivot_wider(names_from = "pp", values_from = "priors")
    a <- cbind(a, b)

    cc <- as.data.frame(res1[2])
    cc$chr <- kb100$chr.exposure[i]; cc$pos <- kb100$pos.exposure[i]; cc$lead_snp <- kb100$SNP[i]

    pph0_4_all     <- rbind(pph0_4_all, a)
    pph0_4_snp_all <- rbind(pph0_4_snp_all, cc)
    kb100res[[i]]  <- res1

    # PP.H4 > 0.8 -> keep candidate causal SNPs
    if (kb100res[[i]][["summary"]][6] > 0.8) {
      need_result <- kb100res[[i]]$results %>% filter(SNP.PP.H4 > 0.0001)
      need_result$lead_snp <- kb100$SNP[i]
      final_snp_res1 <- rbind(final_snp_res1, need_result)
    }

    # extra reference ratios
    n1 <- a$PP.H4.abf + a$PP.H3.abf
    n2 <- a$PP.H4.abf / a$PP.H3.abf
    n3 <- a$PP.H4.abf / (a$PP.H3.abf + a$PP.H4.abf)
    n  <- data.frame(n1, n2, n3)
    colnames(n) <- c("PPH3_PPH4", "PPH4/PPH3", "PPH4/PPH3_PPH4")
    n$chr <- kb100$chr.exposure[i]; n$pos <- kb100$pos.exposure[i]; n$lead_snp <- kb100$SNP[i]
    pph3_4 <- rbind(pph3_4, n)

    cat("  done lead SNP", i, "/", nrow(kb100), ":", kb100$SNP[i], "\n")
  }

  # ---- save this OCR's results ----
  write.csv(pph0_4_all,     paste0(out_dir, "/", peak, "_100KB_leadsnp.csv"), row.names = FALSE)
  write.csv(final_snp_res1, paste0(out_dir, "/", peak, "_100KB_snp(PPH4_0.8).csv"), row.names = FALSE)
  write.csv(pph0_4_snp_all, paste0(out_dir, "/", peak, "_100KB_ALL_leadsnp.csv"), row.names = FALSE)
  write.csv(pph3_4,         paste0(out_dir, "/", peak, "_100KB_pph3_4_ref.csv"), row.names = FALSE)

  cat(peak, "done. Results in:", out_dir, "\n\n")
}
