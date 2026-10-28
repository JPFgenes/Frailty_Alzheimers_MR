# Frailty and Alzheimer's disease: MR and colocalisation code

1. **`01_univariable_mr.R`** — forward GF/F1–F6 to AD (Kunkle and Wightman), reverse AD to frailty; IVW and sensitivity estimators, MR-PRESSO, heterogeneity, leave-one-out, and BH-FDR summaries.
2. **`02_multivariable_mr.R`** — separate frailty + education + income models, then frailty + education + income + longevity, each against both AD outcomes; union of significant instruments, joint clumping and MVMR conditional instrument strength.
3. **`03_colocalisation.R`** — F2/F5/F6 versus AD within ±1 Mb windows, coloc.abf and H0–H4 posterior probabilities.

## Usage

All GWAS files must be converted to these columns:

`SNP, CHR, BP, effect_allele, other_allele, beta, se, pval`

Use an optional `eaf` column for UVMR; colocalisation requires genuine frequency and sample size metadata (see script). Ensure every file uses the *same genome build* (hg19/GRCh37 for the hard-coded APOE coordinates). Binary AD effect sizes should be **log odds ratios**, not ORs.

Provide a local PLINK executable and 1000 Genomes EUR PLINK reference files, then change `plink_bin` and `ld_reference_prefix` in the R configuration. Install the declared R packages, including `MRPRESSO` for its optional sensitivity analysis. Run each script using `Rscript <filename>`.


The uploaded notebook contains repeated factor-specific UVMR and MVMR blocks, `MRPRESSO` will be skipped if it is not installed; review warnings.
