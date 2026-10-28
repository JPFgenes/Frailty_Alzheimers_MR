packages <- c('data.table','dplyr','TwoSampleMR')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install missing packages: ', paste(missing,collapse=', '))

cfg <- list(
  input_dir='data', output_dir='results/uvmr',
  frailty=setNames(paste0('frailty_',c('GF',paste0('F',1:6)),'.tsv'),c('GF',paste0('F',1:6))),
  ad=c(Kunkle='ad_kunkle.tsv',Wightman='ad_wightman.tsv'),
  plink_bin='plink', ld_reference_prefix='reference/EUR',
  p_threshold=5e-8, clump_kb=10000, clump_r2=.001,
  apoe_start=44909039L, apoe_end=45912650L,
  min_F=10, seed=123, presso_bootstraps=1000
)

read_gwas <- function(file) {
  path <- file.path(cfg$input_dir,file)
  if(!file.exists(path)) stop('File not found: ',path)
  d <- data.table::fread(path,data.table=FALSE)
  req <- c('SNP','CHR','BP','effect_allele','other_allele','beta','se','pval')
  if(!all(req %in% names(d))) stop('Missing columns in ',path,': ',paste(setdiff(req,names(d)),collapse=', '))
  if(!'eaf' %in% names(d)) d$eaf <- NA_real_
  d <- d[,c(req,'eaf'),drop=FALSE]
  d$SNP <- as.character(d$SNP)
  d$CHR <- as.integer(d$CHR); d$BP <- as.integer(d$BP)
  d$effect_allele <- toupper(as.character(d$effect_allele))
  d$other_allele <- toupper(as.character(d$other_allele))
  for(x in c('beta','se','pval','eaf')) d[[x]] <- as.numeric(d[[x]])
  d <- d[!is.na(d$SNP) & !is.na(d$CHR) & !is.na(d$BP) &
         d$CHR %in% 1:22 & is.finite(d$beta) & is.finite(d$se) & d$se>0 &
         is.finite(d$pval) & d$pval>=0 & d$pval<=1 &
         d$effect_allele %in% c('A','C','G','T') & d$other_allele %in% c('A','C','G','T') &
         d$effect_allele!=d$other_allele,]
  d <- d[!(d$CHR==19L & d$BP>=cfg$apoe_start & d$BP<=cfg$apoe_end),]
  if(anyDuplicated(d$SNP)) stop('Duplicate SNP identifiers in ',path)
  d
}

# TwoSampleMR expects exposure/outcome-specific column names. The F check
# is the simple per-SNP approximation beta^2/se^2 (not conditional MVMR F).
prepare_exposure <- function(d,trait) {
  d <- d[d$pval<cfg$p_threshold & (d$beta/d$se)^2>cfg$min_F,]
  if(!nrow(d)) return(NULL)
  x <- TwoSampleMR::format_data(d,type='exposure',snp_col='SNP',beta_col='beta',
    se_col='se',effect_allele_col='effect_allele',other_allele_col='other_allele',
    eaf_col='eaf',pval_col='pval',chr_col='CHR',pos_col='BP')
  x$exposure <- trait; x$id.exposure <- trait
  # Local PLINK binary and 1000 Genomes EUR reference prefix are user-supplied.
  TwoSampleMR::clump_data(x,clump_kb=cfg$clump_kb,clump_r2=cfg$clump_r2,
    clump_p1=cfg$p_threshold,bfile=cfg$ld_reference_prefix,plink_bin=cfg$plink_bin)
}

prepare_outcome <- function(d,trait,selected_snps) {
  d <- d[d$SNP %in% selected_snps,]
  if(!nrow(d)) return(NULL)
  x <- TwoSampleMR::format_data(d,type='outcome',snp_col='SNP',beta_col='beta',
    se_col='se',effect_allele_col='effect_allele',other_allele_col='other_allele',
    eaf_col='eaf',pval_col='pval',chr_col='CHR',pos_col='BP')
  x$outcome <- trait; x$id.outcome <- trait
  x
}

write_table <- function(x,path) if(!is.null(x) && NROW(x)>0)
  utils::write.csv(as.data.frame(x),path,row.names=FALSE)

methods <- c('mr_ivw','mr_egger_regression','mr_weighted_median',
             'mr_penalised_weighted_median','mr_weighted_mode')

run_pair <- function(exp_data,exp_name,out_data,out_name) {
  dir <- file.path(cfg$output_dir,exp_name,paste0('to_',out_name))
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  instrument <- prepare_exposure(exp_data,exp_name)
  if(is.null(instrument)||!nrow(instrument)) {warning('No instruments: ',exp_name);return(NULL)}
  outcome <- prepare_outcome(out_data,out_name,instrument$SNP)
  if(is.null(outcome)) {warning('No outcome matches: ',exp_name,' -> ',out_name);return(NULL)}
  h <- TwoSampleMR::harmonise_data(instrument,outcome,action=2)
  h <- h[h$mr_keep %in% TRUE,]
  if(!nrow(h)) {warning('No harmonised instruments: ',exp_name,' -> ',out_name);return(NULL)}
  write_table(h,file.path(dir,'harmonised.csv'))
  avail <- methods[c(1,2,3,4,5)]
  if(nrow(h)<3) avail <- 'mr_ivw' # for a single SNP, use Wald ratio instead
  if(nrow(h)==1) avail <- 'mr_wald_ratio'
  mr <- tryCatch(TwoSampleMR::mr(h,method_list=avail),error=function(e){warning(conditionMessage(e));NULL})
  write_table(mr,file.path(dir,'mr_estimates.csv'))
  if(nrow(h)>=3) {
    het <- tryCatch(TwoSampleMR::mr_heterogeneity(h),error=function(e)NULL)
    if(!is.null(het) && nrow(het)) het$I2 <- pmax(0,(het$Q-het$Q_df)/het$Q)
    write_table(het,file.path(dir,'heterogeneity.csv'))
    write_table(tryCatch(TwoSampleMR::mr_pleiotropy_test(h),error=function(e)NULL),
      file.path(dir,'egger_intercept.csv'))
    loo <- tryCatch(TwoSampleMR::mr_leaveoneout(h),error=function(e)NULL)
    write_table(loo,file.path(dir,'leave_one_out.csv'))
  }
  if(nrow(h)>=4 && requireNamespace('MRPRESSO',quietly=TRUE)) {
    set.seed(cfg$seed)
    presso <- tryCatch(MRPRESSO::mr_presso(BetaOutcome='beta.outcome',
      BetaExposure='beta.exposure',SdOutcome='se.outcome',SdExposure='se.exposure',
      OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=h,
      NbDistribution=max(cfg$presso_bootstraps,nrow(h)+1L),SignifThreshold=.05),
      error=function(e){warning(conditionMessage(e));NULL})
    if(!is.null(presso)) capture.output(print(presso),file=file.path(dir,'mr_presso.txt'))
  }
  if(!requireNamespace('MRPRESSO',quietly=TRUE))
    message('MRPRESSO not installed; MR-PRESSO skipped: ',exp_name,' -> ',out_name)
  if(is.null(mr)) return(NULL)
  mr$direction <- if(exp_name %in% names(frailty)) 'frailty_to_AD' else 'AD_to_frailty'
  mr
}

frailty <- lapply(cfg$frailty,read_gwas)
ad <- lapply(cfg$ad,read_gwas)
results <- list()
for(f in names(frailty)) for(a in names(ad)) {
  for(z in list(list(frailty[[f]],f,ad[[a]],a),list(ad[[a]],a,frailty[[f]],f))) {
    message('UVMR: ',z[[2]],' -> ',z[[4]])
    ans <- tryCatch(run_pair(z[[1]],z[[2]],z[[3]],z[[4]]),
      error=function(e){warning('Failed: ',conditionMessage(e));NULL})
    if(!is.null(ans)) results[[length(results)+1L]] <- ans
  }
}
if(length(results)) {
  all <- dplyr::bind_rows(results)
  # FDR for the primary IVW method separately by direction and outcome.
  all$p_fdr_ivw <- NA_real_
  for(direction in unique(all$direction)) for(out in unique(all$outcome)) {
    ii <- which(all$direction==direction & all$outcome==out & all$method=='Inverse variance weighted')
    if(length(ii)) all$p_fdr_ivw[ii] <- stats::p.adjust(all$pval[ii],method='BH')
  }
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  write_table(all,file.path(cfg$output_dir,'all_mr_results.csv'))
}
message('Completed UVMR; inspect warnings and QC outputs before interpreting results.')
