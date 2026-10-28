packages <- c('data.table','dplyr','TwoSampleMR','MVMR')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install required packages: ', paste(missing, collapse=', '))

cfg <- list(
  input_dir = 'data', output_dir = 'results/mvmr',
  frailty = setNames(paste0('frailty_',c('GF',paste0('F',1:6)),'.tsv'),c('GF',paste0('F',1:6))),
  education = 'education.tsv', income = 'income.tsv', longevity = 'longevity.tsv',
  outcomes = c(Kunkle='ad_kunkle.tsv', Wightman='ad_wightman.tsv'),
  plink_bin = 'plink', ld_reference_prefix = 'reference/EUR',
  p_threshold = 5e-8, clump_kb = 10000, clump_r2 = .001,
  apoe_chr = 19L, apoe_start = 44909039L, apoe_end = 45912650L,
  # Select population-appropriate genome build; coordinates above are hg19/GRCh37.
  # All input data must be mapped to the same genome build and variant identifiers.
  min_snps = 5L
)

# ALL inputs must already use these standard columns:
# SNP, CHR, BP, effect_allele, other_allele, beta, se, pval
# Additional columns are ignored. SNP must use a common identifier across files.
# beta must be aligned to effect_allele (log-OR for binary outcome).
read_sumstats <- function(filename) {
  path <- file.path(cfg$input_dir,filename)
  if(!file.exists(path)) stop('Missing file: ',path)
  x <- data.table::fread(path, data.table=FALSE)
  required <- c('SNP','CHR','BP','effect_allele','other_allele','beta','se','pval')
  if(!all(required %in% names(x))) stop('Missing columns in ',path,': ',paste(setdiff(required,names(x)),collapse=', '))
  x <- x[,required]
  x$SNP <- as.character(x$SNP)
  x$CHR <- as.integer(x$CHR); x$BP <- as.integer(x$BP)
  x$effect_allele <- toupper(as.character(x$effect_allele))
  x$other_allele <- toupper(as.character(x$other_allele))
  x$beta <- as.numeric(x$beta); x$se <- as.numeric(x$se); x$pval <- as.numeric(x$pval)
  x <- x[!is.na(x$SNP) & !is.na(x$CHR) & !is.na(x$BP) &
         is.finite(x$beta) & is.finite(x$se) & x$se>0 &
         is.finite(x$pval) & x$pval>=0 & x$pval<=1 &
         x$effect_allele %in% c('A','C','G','T') & x$other_allele %in% c('A','C','G','T') &
         x$effect_allele != x$other_allele,]
  x <- x[!(x$CHR==cfg$apoe_chr & x$BP>=cfg$apoe_start & x$BP<=cfg$apoe_end),]
  if(anyDuplicated(x$SNP)) stop('Duplicate SNP IDs in ',path,'; resolve before running')
  x
}

# Align each exposure and outcome to the same SNP effect allele.
# Conservatively discard palindromic A/T and C/G SNPs in absence of EAF data.
orient <- function(target, reference) {
  z <- merge(reference[,c('SNP','effect_allele','other_allele')], target,
             by='SNP',suffixes=c('.ref',''),sort=FALSE)
  pal <- (z$effect_allele.ref=='A' & z$other_allele.ref=='T') |
         (z$effect_allele.ref=='T' & z$other_allele.ref=='A') |
         (z$effect_allele.ref=='C' & z$other_allele.ref=='G') |
         (z$effect_allele.ref=='G' & z$other_allele.ref=='C')
  same <- z$effect_allele==z$effect_allele.ref & z$other_allele==z$other_allele.ref
  swapped <- z$effect_allele==z$other_allele.ref & z$other_allele==z$effect_allele.ref
  # Complement matches support differing strand representations for non-palindromic SNPs.
  comp <- c(A='T',T='A',C='G',G='C')
  comp_same <- comp[z$effect_allele]==z$effect_allele.ref & comp[z$other_allele]==z$other_allele.ref
  comp_swap <- comp[z$effect_allele]==z$other_allele.ref & comp[z$other_allele]==z$effect_allele.ref
  ok <- !pal & (same|swapped|comp_same|comp_swap)
  z <- z[ok,]
  flip <- (swapped|comp_swap)[ok]
  z$beta[flip] <- -z$beta[flip]
  z[,c('SNP','beta','se')]
}

# Joint clumping of the UNION of significant instruments for all included exposures.
# The shared id.exposure prevents independent per-trait clumping.
clump_union <- function(exposures) {
  hits <- dplyr::bind_rows(lapply(names(exposures), function(nm) {
    d <- exposures[[nm]]
    d <- d[d$pval < cfg$p_threshold,]
    d[,c('SNP','CHR','BP','pval')]
  }))
  if(!nrow(hits)) stop('No genome-wide significant SNPs for this model')
  hits <- hits[order(hits$pval),]
  hits <- hits[!duplicated(hits$SNP),]
  hits$id.exposure <- 'joint'; hits$exposure <- 'joint'
  hits$chr_name <- hits$CHR; hits$chrom_start <- hits$BP
  hits$pval.exposure <- hits$pval
  hits$beta.exposure <- 0; hits$se.exposure <- 1
  hits$effect_allele.exposure <- 'A'; hits$other_allele.exposure <- 'G'
  out <- TwoSampleMR::clump_data(hits,
    clump_kb=cfg$clump_kb, clump_r2=cfg$clump_r2,
    clump_p1=cfg$p_threshold,
    bfile=cfg$ld_reference_prefix, plink_bin=cfg$plink_bin)
  unique(out$SNP)
}

run_model <- function(frailty_name, exposure_set, outcome_name, outcome) {
  model_name <- paste(names(exposure_set),collapse='_')
  selected <- clump_union(exposure_set)
  # Reference orientation from the first exposure containing each SNP, not
  # necessarily the frailty GWAS (union SNPs can originate from covariates).
  ref <- dplyr::bind_rows(exposure_set)
  ref <- ref[ref$SNP %in% selected,]
  ref <- ref[!duplicated(ref$SNP),c('SNP','effect_allele','other_allele')]
  aligned <- lapply(exposure_set, orient, reference=ref)
  aligned_outcome <- orient(outcome, ref)
  shared <- Reduce(intersect,c(lapply(aligned,`[[`,'SNP'),list(aligned_outcome$SNP)))
  if(length(shared)<cfg$min_snps) {
    warning('Too few complete harmonised instruments: ',frailty_name,' / ',model_name,' / ',outcome_name)
    return(NULL)
  }
  get_values <- function(d,col) d[[col]][match(shared,d$SNP)]
  BXGs <- do.call(cbind,lapply(aligned,get_values,'beta'))
  seBXGs <- do.call(cbind,lapply(aligned,get_values,'se'))
  colnames(BXGs) <- colnames(seBXGs) <- names(exposure_set)
  BYG <- get_values(aligned_outcome,'beta')
  seBYG <- get_values(aligned_outcome,'se')
  fm <- MVMR::format_mvmr(BXGs=BXGs,BYG=BYG,seBXGs=seBXGs,seBYG=seBYG,RSID=shared)
  fit <- MVMR::ivw_mvmr(fm)
  strength <- tryCatch(MVMR::strength_mvmr(fm),error=function(e) {
    warning('Conditional F-statistic unavailable: ',conditionMessage(e));NULL
  })
  directory <- file.path(cfg$output_dir,frailty_name,model_name,outcome_name)
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(as.data.frame(fit),file.path(directory,'ivw_mvmr.csv'),row.names=TRUE)
  if(!is.null(strength)) utils::write.csv(as.data.frame(strength),file.path(directory,'conditional_strength.csv'),row.names=TRUE)
  utils::write.csv(data.frame(SNP=shared),file.path(directory,'instrument_snps.csv'),row.names=FALSE)
  data.frame(frailty=frailty_name,model=model_name,outcome=outcome_name,n_snps=length(shared))
}

frailty <- lapply(cfg$frailty,read_sumstats)
ea <- read_sumstats(cfg$education)
income <- read_sumstats(cfg$income)
longevity <- read_sumstats(cfg$longevity)
outcomes <- lapply(cfg$outcomes,read_sumstats)

run_log <- list()
for(f in names(frailty)) {
  sets <- list(primary=list(frailty=frailty[[f]],education=ea,income=income),
               longevity_adjusted=list(frailty=frailty[[f]],education=ea,income=income,longevity=longevity))
  for(mode in names(sets)) for(outcome_name in names(outcomes)) {
    message('Running ',f,' / ',mode,' / ',outcome_name)
    result <- tryCatch(run_model(f,sets[[mode]],outcome_name,outcomes[[outcome_name]]),
                       error=function(e) {warning('Model failed: ',conditionMessage(e));NULL})
    if(!is.null(result)) {result$adjustment <- mode;run_log[[length(run_log)+1L]] <- result}
  }
}
if(length(run_log)) {
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(dplyr::bind_rows(run_log),file.path(cfg$output_dir,'run_log.csv'),row.names=FALSE)
}
