packages <- c('data.table','coloc')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install missing packages: ',paste(missing,collapse=', '))

cfg <- list(
  input_dir='data', output_dir='results/coloc',
  frailty=c(F2='frailty_F2.tsv',F5='frailty_F5.tsv',F6='frailty_F6.tsv'),
  outcomes=c(Kunkle='ad_kunkle.tsv',Wightman='ad_wightman.tsv'),
  # Case counts must match the EXACT outcome data files used.
  cases=c(Kunkle=21982,Wightman=39918),
  controls=c(Kunkle=41944,Wightman=358140),
  window_bp=1000000L, p_threshold=5e-8,
  apoe_chr=19L,apoe_start=44909039L,apoe_end=45912650L,
  mhc_chr=6L,mhc_start=25000000L,mhc_end=34000000L,
  min_snps=50L
)
# Input standard columns: SNP, CHR, BP, effect_allele, other_allele,
# beta, se, pval, MAF, N. MAF must be empirical and between 0 and 0.5.
# Binary outcomes additionally require accurate N and case fraction above.
# Frailty traits are treated as quantitative; provide sdY if NOT on SD scale.
# All positions must be on the same build (APOE coordinates above: GRCh37).
read_stats <- function(filename) {
  path <- file.path(cfg$input_dir,filename)
  if(!file.exists(path)) stop('Missing: ',path)
  x <- data.table::fread(path,data.table=FALSE)
  required <- c('SNP','CHR','BP','effect_allele','other_allele','beta','se','pval','MAF','N')
  if(!all(required %in% names(x))) stop('Missing columns in ',path,': ',paste(setdiff(required,names(x)),collapse=', '))
  x <- x[,unique(c(required,intersect('sdY',names(x))))]
  x$CHR <- as.integer(x$CHR);x$BP <- as.integer(x$BP)
  x$SNP <- as.character(x$SNP)
  x$effect_allele <- toupper(x$effect_allele);x$other_allele <- toupper(x$other_allele)
  numeric_fields <- c('beta','se','pval','MAF','N')
  for(nm in numeric_fields) x[[nm]] <- as.numeric(x[[nm]])
  x <- x[!is.na(x$CHR)&!is.na(x$BP)&!is.na(x$SNP)&
         x$CHR %in% 1:22 & x$effect_allele %in% c('A','C','G','T') &
         x$other_allele %in% c('A','C','G','T') & x$effect_allele != x$other_allele & is.finite(x$beta)&is.finite(x$se)&x$se>0&
         is.finite(x$pval)&x$pval>=0&x$pval<=1&
         is.finite(x$MAF)&x$MAF>0&x$MAF<=.5&
         is.finite(x$N)&x$N>0,]
  x <- x[!(x$CHR==cfg$apoe_chr & x$BP>=cfg$apoe_start & x$BP<=cfg$apoe_end) &
         !(x$CHR==cfg$mhc_chr & x$BP>=cfg$mhc_start & x$BP<=cfg$mhc_end),]
  if(anyDuplicated(x$SNP)) stop('Resolve duplicate SNPs in ',path)
  x
}

# Non-palindromic allele orientation: flip outcome beta when alleles reversed.
harmonise_region <- function(exposure,outcome) {
  z <- merge(exposure,outcome,by='SNP',suffixes=c('.x','.y'))
  same <- z$effect_allele.x==z$effect_allele.y & z$other_allele.x==z$other_allele.y
  swap <- z$effect_allele.x==z$other_allele.y & z$other_allele.x==z$effect_allele.y
  comp <- c(A='T',T='A',C='G',G='C')
  csame <- comp[z$effect_allele.y]==z$effect_allele.x & comp[z$other_allele.y]==z$other_allele.x
  cswap <- comp[z$effect_allele.y]==z$other_allele.x & comp[z$other_allele.y]==z$effect_allele.x
  pal <- (z$effect_allele.x=='A' & z$other_allele.x=='T') |
         (z$effect_allele.x=='T' & z$other_allele.x=='A') |
         (z$effect_allele.x=='C' & z$other_allele.x=='G') |
         (z$effect_allele.x=='G' & z$other_allele.x=='C')
  keep <- !pal&(same|swap|csame|cswap)&z$CHR.x==z$CHR.y & z$BP.x==z$BP.y
  z <- z[keep,]
  flip <- (swap|cswap)[keep]
  z$beta.y[flip] <- -z$beta.y[flip]
  z
}

make_dataset <- function(z,which,trait,outcome_name=NULL) {
  suffix <- paste0('.',which)
  get <- function(name) z[[paste0(name,suffix)]]
  dat <- list(beta=get('beta'),varbeta=get('se')^2,
              snp=z$SNP,MAF=get('MAF'),type=trait)
  ns <- get('N')
  if(length(unique(ns))==1L) dat$N <- ns[1] else dat$N <- ns
  if(trait=='cc') {
    dat$s <- cfg$cases[[outcome_name]]/(cfg$cases[[outcome_name]]+cfg$controls[[outcome_name]])
    if(any(ns != cfg$cases[[outcome_name]]+cfg$controls[[outcome_name]]))
      warning('Outcome N differs from configured case/control counts: ',outcome_name)
  } else if(paste0('sdY',suffix) %in% names(z)) {
    sd_vals <- get('sdY')
    if(length(unique(sd_vals))==1L) dat$sdY <- sd_vals[1]
  }
  dat
}

frailty <- lapply(cfg$frailty,read_stats)
outcomes <- lapply(cfg$outcomes,read_stats)
dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
result_rows <- list()

for(f in names(frailty)) {
  expdat <- frailty[[f]]
  hits <- expdat[expdat$pval < cfg$p_threshold,]
  if(!nrow(hits)) {warning('No significant lead SNP candidates for ',f);next}
  hits <- hits[order(hits$pval),]
  # Avoid overlapping lead-centred windows without requiring extra LD data.
  leads <- hits[0,]
  for(i in seq_len(nrow(hits))) {
    candidate <- hits[i,,drop=FALSE]
    if(nrow(leads) && any(leads$CHR==candidate$CHR &
                          abs(leads$BP-candidate$BP)<=2*cfg$window_bp)) next
    leads <- rbind(leads,candidate)
  }
  for(outcome_name in names(outcomes)) for(i in seq_len(nrow(leads))) {
    lead <- leads[i,,drop=FALSE]
    lo <- max(1L, lead$BP-cfg$window_bp);hi <- lead$BP+cfg$window_bp
    a <- expdat[expdat$CHR==lead$CHR & expdat$BP>=lo & expdat$BP<=hi,]
    b <- outcomes[[outcome_name]]
    b <- b[b$CHR==lead$CHR & b$BP>=lo & b$BP<=hi,]
    z <- harmonise_region(a,b)
    if(nrow(z)<cfg$min_snps) {warning('Insufficient shared SNPs: ',f,' ',outcome_name,' ',lead$SNP);next}
    d1 <- make_dataset(z,'x','quant')
    d2 <- make_dataset(z,'y','cc',outcome_name)
    ans <- tryCatch(coloc::coloc.abf(dataset1=d1,dataset2=d2),
      error=function(e) {warning('Coloc failed: ',conditionMessage(e));NULL})
    if(is.null(ans)) next
    s <- ans$summary
    result_rows[[length(result_rows)+1L]] <- data.frame(
      frailty=f,outcome=outcome_name,lead_snp=lead$SNP,
      chr=lead$CHR,window_start=lo,window_end=hi,nsnps=length(z$SNP),
      PP_H0=unname(s['PP.H0.abf']),PP_H1=unname(s['PP.H1.abf']),
      PP_H2=unname(s['PP.H2.abf']),PP_H3=unname(s['PP.H3.abf']),
      PP_H4=unname(s['PP.H4.abf']))
  }
}
if(length(result_rows)) {
  results <- do.call(rbind,result_rows)
  write.csv(results,file.path(cfg$output_dir,'coloc_abf_summary.csv'),row.names=FALSE)
} else warning('No evaluable colocalisation regions')
