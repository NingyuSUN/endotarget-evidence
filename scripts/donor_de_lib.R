# Donor-level count model. All numerical choices are read from the pre-effect config.
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(statmod); library(jsonlite)
})

validate_donor_inputs <- function(counts, design, clinical, cohort, cfg) {
  if (!cohort %in% names(cfg$cohorts)) stop("unknown cohort")
  if (!is.matrix(counts) || !is.numeric(counts) || nrow(counts)<2L) stop("invalid counts matrix")
  if (any(!is.finite(counts)) || any(counts<0)) stop("counts must be finite and nonnegative")
  if (isTRUE(cfg$cohorts[[cohort]]$integer_counts) && any(counts!=floor(counts))) stop("integer count contract failed")
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts)) ||
      any(!grepl("^ENSG[0-9]{11}$",rownames(counts)))) stop("invalid unique stable genes")
  if (is.null(colnames(counts)) || anyDuplicated(colnames(counts))) stop("invalid count donor IDs")
  if (!all(c("donor_id","disease") %in% names(clinical)) ||
      anyDuplicated(clinical$donor_id)) stop("invalid clinical donor IDs")
  if (!identical(colnames(counts),as.character(clinical$donor_id)) ||
      !identical(rownames(design),as.character(clinical$donor_id))) stop("exact donor order required")
  if (any(!clinical$disease %in% c("non-ENDO","ENDO"))) stop("unknown disease")
  ng <- table(factor(clinical$disease,levels=c("non-ENDO","ENDO")))
  if (any(ng<cfg$cohorts[[cohort]]$minimum_per_group)) stop("disease group donor minimum failed")
  cell_column <- cfg$cohorts[[cohort]]$cell_column
  if (!cell_column %in% names(clinical) ||
      any(!is.finite(clinical[[cell_column]])) ||
      any(clinical[[cell_column]]<cfg$design$min_cells)) stop("donor cell threshold failed")
  if (any(colSums(counts)<=0)) stop("zero count library")
  if (!is.matrix(design) || !is.numeric(design) || any(!is.finite(design))) stop("invalid numeric design")
  coefficient <- cfg$disease_coefficient
  if (sum(colnames(design)==coefficient)!=1L ||
      !identical(as.numeric(design[,coefficient]),as.numeric(clinical$disease=="ENDO"))) stop("disease coefficient mismatch")
  rank <- qr(design,tol=cfg$fit$rank_tolerance)$rank
  if (rank!=ncol(design)) stop("design is not full rank")
  df <- nrow(design)-rank
  if (df<cfg$minimum_residual_df) stop("insufficient residual degrees of freedom")
  list(n_donors=ncol(counts),input_genes=nrow(counts),rank=rank,residual_df=df,
       non_ENDO=unname(ng[1]),ENDO=unname(ng[2]))
}

.fit_donor_model <- function(counts, design, clinical, cohort, cfg) {
  audit <- validate_donor_inputs(counts,design,clinical,cohort,cfg)
  f <- cfg$filtering
  raw_lib <- colSums(counts)
  support_n <- min(table(factor(clinical$disease,levels=c("non-ENDO","ENDO"))))
  if (support_n>f$large_n) support_n <- f$large_n+(support_n-f$large_n)*f$min_prop
  cpm_cutoff <- f$min_count/median(raw_lib)*1e6
  filter_audit <- list(cpm_cutoff=cpm_cutoff,minimum_support_donors=support_n,
                      min_total_count=f$min_total_count,median_raw_library_size=median(raw_lib))
  group <- factor(clinical$disease,levels=c("non-ENDO","ENDO"))
  keep <- edgeR::filterByExpr(counts,group=group,lib.size=colSums(counts),
             min.count=f$min_count,min.total.count=f$min_total_count,
             large.n=f$large_n,min.prop=f$min_prop)
  if (sum(keep)<2L) stop("too few expression-pass genes")
  y <- edgeR::DGEList(counts=counts[keep,,drop=FALSE],group=group)
  if (any(y$samples$lib.size<=0)) stop("zero filtered library")
  n <- cfg$normalization
  y <- edgeR::calcNormFactors(y,method=n$method,refColumn=NULL,
             logratioTrim=n$logratio_trim,sumTrim=n$sum_trim,
             doWeighting=n$do_weighting,Acutoff=n$A_cutoff,p=n$p)
  v <- cfg$voom
  va <- list(counts=y,design=design,lib.size=NULL,normalize.method=v$normalize_method,
             block=NULL,correlation=NULL,weights=NULL,span=v$span,plot=v$plot,save.plot=v$save_plot)
  adaptive_supported <- "adaptive.span" %in% names(formals(limma::voom))
  if (adaptive_supported) va$adaptive.span <- v$adaptive_span
  # limma 3.62.1 has the historical fixed-span implementation.
  if (!adaptive_supported && as.character(packageVersion("limma"))!="3.62.1")
    stop("unreviewed voom implementation without adaptive.span")
  vo <- do.call(limma::voom,va)
  validate_voom_weights(vo$weights)
  fit <- limma::lmFit(vo,design=design,method=cfg$fit$method)
  e <- cfg$ebayes
  ea <- list(fit=fit,proportion=e$proportion,stdev.coef.lim=e$stdev_coef_lim,
             trend=e$trend,robust=e$robust,winsor.tail.p=e$winsor_tail_p)
  legacy_supported <- "legacy" %in% names(formals(limma::eBayes))
  if (legacy_supported) ea$legacy <- e$legacy
  if (!legacy_supported && as.character(packageVersion("limma"))!="3.62.1")
    stop("unreviewed eBayes implementation without legacy argument")
  fit <- do.call(limma::eBayes,ea)
  co <- match(cfg$disease_coefficient,colnames(design))
  logfc <- fit$coefficients[,co]
  statistic <- fit$t[,co]
  p <- fit$p.value[,co]
  se <- sqrt(fit$s2.post)*fit$stdev.unscaled[,co]
  degrees <- rep_len(fit$df.total,length(logfc))
  finite <- is.finite(logfc)&is.finite(statistic)&is.finite(p)&
            is.finite(se)&se>0&is.finite(degrees)&degrees>0&
            is.finite(fit$df.residual)&fit$df.residual>=cfg$minimum_residual_df
  q <- rep(NA_real_,length(p))
  q[finite] <- p.adjust(p[finite],method=cfg$testing$adjust_method,n=length(p))
  critical <- qt(1-(1-cfg$testing$confidence_level)/2,df=degrees)
  tab <- data.frame(gene_id=rownames(y),logFC=logfc,AveExpr=fit$Amean,t=statistic,
            P.Value=p,adj.P.Val=q,moderated_SE=se,df_residual=fit$df.residual,df_total=degrees,
            CI_low=logfc-critical*se,CI_high=logfc+critical*se,
            direction=sign(logfc),fit_estimable=finite,
            null_matching_logCPM=rowMeans(edgeR::cpm(y,log=TRUE,prior.count=.5,
                                                    normalized.lib.sizes=TRUE)),
            row.names=NULL,check.names=FALSE)
  norm <- data.frame(donor_id=colnames(counts),raw_library_size=colSums(counts),
             filtered_library_size=y$samples$lib.size,norm_factor=y$samples$norm.factors,
             effective_library_size=y$samples$lib.size*y$samples$norm.factors,
             row.names=NULL)
  list(table=tab,keep=keep,filter_audit=filter_audit,normalization=norm,voom=vo,fit=fit,audit=audit,
       settings=list(adaptive_span_argument_supported=adaptive_supported,
                     effective_adaptive_span=FALSE,legacy_argument_supported=legacy_supported,
                     effective_legacy_ebayes=TRUE))
}

fit_donor_model <- function(...) {
  warnings <- character()
  result <- withCallingHandlers(.fit_donor_model(...), warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  result$warnings <- warnings
  result
}

validate_voom_weights <- function(weights) {
  if (!is.matrix(weights) || any(!is.finite(weights)) || any(weights<=0))
    stop("voom precision weights must be finite and positive")
  invisible(TRUE)
}
