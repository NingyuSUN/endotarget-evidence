# Shared functions for the secretory-phase program-level disease contrast
# (docs/SECRETORY_PROGRAM_PROTOCOL.md). Builds on scripts/mde_lib.R for hashing, I/O, freeze
# verification, the positive-control design, the v1 estimator copy and the MDE solver.
# Sourcing this file runs nothing.
source("scripts/mde_lib.R")

# edgeR expression filter with the phase group and TMM factors, exactly the steps that open
# mde_fit_subset (and the frozen positive-control fits); returns the DGEList for the log-CPM.
sp_filter_norm <- function(counts,group,cfg) {
  f <- cfg$filtering
  keep <- edgeR::filterByExpr(counts,group=group,lib.size=colSums(counts),min.count=f$min_count,
                              min.total.count=f$min_total_count,large.n=f$large_n,min.prop=f$min_prop)
  if(sum(keep)<2L) stop("too few expression-pass genes")
  y <- edgeR::DGEList(counts=counts[keep,,drop=FALSE],group=group)
  n <- cfg$normalization
  y <- edgeR::calcNormFactors(y,method=n$method,refColumn=NULL,logratioTrim=n$logratio_trim,sumTrim=n$sum_trim,
                              doWeighting=n$do_weighting,Acutoff=n$A_cutoff,p=n$p)
  list(keep=keep,y=y,raw_library_size=colSums(counts),filtered_library_size=y$samples$lib.size,norm_factor=y$samples$norm.factors)
}
# voom log-CPM: log2((count + 0.5) / (lib.size * norm.factor + 1) * 1e6)
sp_logcpm <- function(y) {
  lib <- y$samples$lib.size*y$samples$norm.factors
  t(log2(t(y$counts+0.5)/(lib+1)*1e6))
}
# The voom object with the configured settings (same call as mde_fit_subset) for the set test.
sp_voom <- function(y,design,cfg) {
  v <- cfg$voom
  va <- list(counts=y,design=design,lib.size=NULL,normalize.method=v$normalize_method,block=NULL,correlation=NULL,
             weights=NULL,span=v$span,plot=v$plot,save.plot=v$save_plot)
  if("adaptive.span" %in% names(formals(limma::voom))) va$adaptive.span <- v$adaptive_span else if(as.character(packageVersion("limma"))!="3.62.1") stop("unreviewed voom implementation without adaptive.span")
  do.call(limma::voom,va)
}
# Signature from the other cohort's frozen phase table: estimable genes at BH q < threshold,
# with their log2 fold changes as weights; flags membership in the test cohort's family.
sp_signature <- function(other,own_family,q_threshold) {
  need <- c("gene_id","logFC","adj.P.Val","fit_estimable")
  if(!all(need %in% names(other))) stop("signature table lacks the required columns")
  est <- as.logical(other$fit_estimable) & is.finite(other$logFC) & is.finite(other$adj.P.Val)
  sel <- which(est & other$adj.P.Val<q_threshold & other$logFC!=0)
  g <- as.character(other$gene_id[sel])
  if(anyDuplicated(g)) stop("duplicate gene IDs in the signature")
  data.frame(gene_id=g,beta=other$logFC[sel],sign=sign(other$logFC[sel]),in_family=g %in% own_family,stringsAsFactors=FALSE)
}
# Donor-level scores over the signature genes present in the family: the projection onto the
# signature direction (1 = the signature cohort's full phase shift) and the mean signed log2 shift.
sp_scores <- function(logcpm,sig,pro_ix) {
  s <- sig[sig$in_family,,drop=FALSE]
  if(nrow(s)<2L) stop("fewer than two signature genes in the family")
  x <- logcpm[s$gene_id,,drop=FALSE]
  m <- rowMeans(x[,pro_ix,drop=FALSE]); d <- x-m
  data.frame(donor_id=colnames(logcpm),projection=as.numeric(crossprod(s$beta,d))/sum(s$beta^2),
             mean_signed_shift=as.numeric(crossprod(s$sign,d))/nrow(s),stringsAsFactors=FALSE)
}
# Ordinary least squares on a numeric design; one coefficient's estimate, SE, t, df, two-sided p
# and confidence interval. Rank deficiency or too few residual df gives NOT_ESTIMABLE.
sp_lm <- function(y,design,coefficient,mcfg) {
  if(length(y)!=nrow(design) || any(!is.finite(y))) stop("response/design mismatch")
  co <- match(coefficient,colnames(design)); if(is.na(co)) stop(paste("coefficient missing:",coefficient))
  rank <- qr(design,tol=mcfg$rank_tolerance)$rank; df <- nrow(design)-rank
  base <- list(n=nrow(design),columns=colnames(design),rank=rank,residual_df=df,coefficient=coefficient)
  if(rank!=ncol(design)) return(c(base,list(status="NOT_ESTIMABLE",reason="rank_deficient")))
  if(df<mcfg$minimum_residual_df) return(c(base,list(status="NOT_ESTIMABLE",reason="residual_df_below_minimum")))
  fit <- lm.fit(design,y); s2 <- sum(fit$residuals^2)/df
  V <- s2*solve(crossprod(design)); est <- unname(fit$coefficients[co]); se <- sqrt(V[co,co])
  tval <- est/se; p <- 2*pt(-abs(tval),df); crit <- qt(1-(1-mcfg$confidence_level)/2,df)
  c(base,list(status="ESTIMATED",estimate=est,se=se,t=tval,p=p,ci_low=est-crit*se,ci_high=est+crit*se,
              residual_sd=sqrt(s2),coefficients=as.list(setNames(unname(fit$coefficients),colnames(design)))))
}
# Area under the ROC curve of a score for cases versus controls (Mann-Whitney form, ties 0.5).
sp_auc <- function(score,is_case) {
  a <- score[is_case]; b <- score[!is_case]
  if(!length(a) || !length(b)) return(NA_real_)
  mean(outer(a,b,">")+0.5*outer(a,b,"=="))
}
# Add the phase-by-disease product column to a positive-control design.
sp_add_interaction <- function(design) {
  if(!all(c("SEC_minus_PRO","ENDO_minus_non_ENDO") %in% colnames(design))) stop("design lacks the phase or disease column")
  cbind(design,SEC_x_ENDO=design[,"SEC_minus_PRO"]*design[,"ENDO_minus_non_ENDO"])
}
# Treatment contrasts for the secretory sub-phase over the configured order; reference = first level present.
sp_subphase_contrasts <- function(sub,levels_cfg) {
  sub <- as.character(sub)
  if(any(is.na(sub)) || any(!sub %in% levels_cfg)) stop("sub-phase label outside the configured levels")
  present <- levels_cfg[levels_cfg %in% unique(sub)]
  if(length(present)<2L) return(list(contrasts=NULL,reference=present,levels_present=present))
  m <- sapply(present[-1],function(l) as.numeric(sub==l))
  m <- matrix(m,nrow=length(sub)); colnames(m) <- paste0("sub_phase:",gsub(" ","_",present[-1]),"_vs_",gsub(" ","_",present[1]))
  list(contrasts=m,reference=present[1],levels_present=present)
}
# Secretory-only design: intercept, disease, age centered on the donors entering the model.
sp_secretory_design <- function(clinical_sec,cfg,extra=NULL) {
  age <- as.numeric(clinical_sec[[cfg$design$age_column]]); disease <- as.numeric(clinical_sec$disease=="ENDO")
  if(any(!is.finite(age)) || length(unique(disease))!=2) stop("secretory design needs both disease groups and ages")
  d <- cbind(intercept=1,ENDO_minus_non_ENDO=disease,age_centered_years=age-mean(age))
  if(!is.null(extra)) d <- cbind(d,extra)
  rownames(d) <- clinical_sec$donor_id; d
}
# Drop constant nuisance columns from a design subset (never the coefficient of interest).
sp_drop_constant <- function(design,keep_cols) {
  const <- apply(design,2,function(v) length(unique(v))==1L); const["intercept"] <- FALSE
  const[keep_cols] <- FALSE
  list(design=design[,!const,drop=FALSE],dropped=colnames(design)[const])
}
# Rotation gene-set test on one coefficient with directional gene weights.
sp_fry <- function(vo,design,coefficient,index,weights_all) {
  co <- match(coefficient,colnames(design)); if(is.na(co)) stop("fry coefficient missing")
  if(length(weights_all)!=nrow(vo)) stop("gene weights must cover every gene of the voom object")
  r <- limma::fry(vo,index=list(signature=index),design=design,contrast=co,gene.weights=weights_all)
  list(n_genes=unname(r[1,"NGenes"]),direction=as.character(r[1,"Direction"]),p=unname(r[1,"PValue"]),
       p_mixed=unname(r[1,"PValue.Mixed"]))
}
