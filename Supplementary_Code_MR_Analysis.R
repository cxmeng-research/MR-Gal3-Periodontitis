# ==============================================================================
# Supplementary Code for: 
# "The Causal Relationship between Circulating Galectin-3, Adiposity, and 
# Periodontitis: A Mendelian Randomization Study"
# 
# Description: This script reproduces the Univariable MR (UVMR), 
# Multivariable MR (MVMR), and Reverse MR analyses described in the manuscript.
# R version requirement: R 4.x+
# Required packages: TwoSampleMR, data.table, dplyr
# ==============================================================================

library(data.table)
library(dplyr)
library(TwoSampleMR)

# Note to reviewers/users: Please insert your OpenGWAS API token here
options(ieugwasr_token = "YOUR_API_TOKEN_HERE")

# Define relative file paths for summary statistics
# Users must update these paths to their local directories
path_gal3 <- "path/to/Galectin-3_ebi-a-GCST90012009.vcf.gz"
path_bmi  <- "path/to/BMI_ukb-b-19953.vcf.gz"
path_t2dm <- "path/to/T2DM_ebi-a-GCST006867.vcf.gz"
path_out  <- "path/to/finngen_R11_K11_GINGIVITIS_PERIODONTAL.gz"

# ------------------------------------------------------------------------------
# Function: Safe VCF parsing and quality control
# ------------------------------------------------------------------------------
process_vcf_safe <- function(path, trait_name) {
  d <- fread(path, skip = "#CHROM")
  setnames(d, 10, "stats")
  d[, c("b", "se", "lp", "af") := tstrsplit(stats, ":", fixed=TRUE)[1:4]]
  df <- data.frame(
    SNP = d$ID, beta = as.numeric(d$b), se = as.numeric(d$se),
    pval = 10^(-as.numeric(d$lp)), eaf = as.numeric(d$af),
    effect_allele = d$ALT, other_allele = d$REF, exposure = trait_name,
    stringsAsFactors = FALSE
  )
  # Exclude missing values and duplicated variants
  df <- df[complete.cases(df[, c("beta", "se", "pval")]), ]
  return(df[!duplicated(df$SNP), ])
}

# Read full summary statistics
exp_gal3_full <- process_vcf_safe(path_gal3, "Galectin-3")
exp_bmi_full  <- process_vcf_safe(path_bmi, "BMI")
exp_t2dm_full <- process_vcf_safe(path_t2dm, "T2DM")
out_raw <- fread(path_out)

# ==============================================================================
# Part 1: Univariable MR (Galectin-3 -> Periodontitis)
# ==============================================================================
# Extract significant SNPs (P < 5e-8) and clump (kb=10000, r2=0.001)
ivs_gal3 <- clump_data(
  format_data(exp_gal3_full[exp_gal3_full$pval < 5e-8, ], type="exposure"), 
  clump_kb = 10000, clump_r2 = 0.001
)

# Extract outcome data matching the IVs
out_uvmr <- format_data(
  as.data.frame(out_raw[rsids %in% ivs_gal3$SNP]), type = "outcome",
  snp_col = "rsids", beta_col = "beta", se_col = "sebeta", 
  eaf_col = "af_alt", effect_allele_col = "alt", 
  other_allele_col = "ref", pval_col = "pval"
)

# Harmonize and compute causal estimates
dat_uvmr <- harmonise_data(ivs_gal3, out_uvmr)
res_uvmr <- mr(dat_uvmr)
res_uvmr_or <- generate_odds_ratios(res_uvmr)

# ==============================================================================
# Part 2: Multivariable MR (Galectin-3 + BMI + T2DM -> Periodontitis)
# ==============================================================================
# Extract union of significant SNPs across all exposures
union_snps <- unique(c(
  exp_gal3_full$SNP[exp_gal3_full$pval < 5e-8], 
  exp_bmi_full$SNP[exp_bmi_full$pval < 5e-8], 
  exp_t2dm_full$SNP[exp_t2dm_full$pval < 5e-8]
))

# Format datasets explicitly to avoid automated parsing errors
fmt_gal3 <- format_data(exp_gal3_full[exp_gal3_full$SNP %in% union_snps, ], type="exposure") %>% mutate(id.exposure="Gal3")
fmt_bmi  <- format_data(exp_bmi_full[exp_bmi_full$SNP %in% union_snps, ], type="exposure") %>% mutate(id.exposure="BMI")
fmt_t2dm <- format_data(exp_t2dm_full[exp_t2dm_full$SNP %in% union_snps, ], type="exposure") %>% mutate(id.exposure="T2DM")

# Perform global clumping on the combined data
combined_mv <- as.data.frame(rbind(fmt_gal3, fmt_bmi, fmt_t2dm))
mv_clumped <- clump_data(combined_mv, clump_kb = 10000, clump_r2 = 0.001)

# Format outcome for MVMR
out_mvmr <- format_data(
  as.data.frame(out_raw[rsids %in% mv_clumped$SNP]), type="outcome",
  snp_col="rsids", beta_col="beta", se_col="sebeta", pval_col="pval"
)

# Harmonize MVMR data
mv_dat <- mv_harmonise_data(list(
  Gal3 = as.data.frame(mv_clumped[mv_clumped$id.exposure=="Gal3", ]),
  BMI  = as.data.frame(mv_clumped[mv_clumped$id.exposure=="BMI", ]),
  T2DM = as.data.frame(mv_clumped[mv_clumped$id.exposure=="T2DM", ])
), out_mvmr)

# Execute MVMR calculation
res_mvmr <- mv_multiple(mv_dat)

# ==============================================================================
# Part 3: Reverse MR (BMI -> Galectin-3)
# ==============================================================================
# Extract BMI IVs explicitly mapped
exp_bmi_fmt_rev <- format_data(
  as.data.frame(exp_bmi_full[exp_bmi_full$pval < 5e-8, ]), type = "exposure",
  snp_col = "SNP", beta_col = "beta", se_col = "se", pval_col = "pval", 
  eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele"
)
ivs_bmi_rev <- clump_data(exp_bmi_fmt_rev, clump_kb = 10000, clump_r2 = 0.001)

# Extract Gal-3 outcome explicitly mapped
out_gal3_sub <- exp_gal3_full[exp_gal3_full$SNP %in% ivs_bmi_rev$SNP, ]
out_gal3_fmt_rev <- format_data(
  as.data.frame(out_gal3_sub), type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se", pval_col = "pval", 
  eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele"
)

# Strict harmonization discarding intermediate palindromic SNPs (action = 2)
dat_rev_mr <- harmonise_data(ivs_bmi_rev, out_gal3_fmt_rev, action = 2)
res_rev_mr <- mr(dat_rev_mr)