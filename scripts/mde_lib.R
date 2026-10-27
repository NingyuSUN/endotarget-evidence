# Shared functions for the minimum detectable effect analysis (docs/MDE_ANALYSIS_PROTOCOL.md).
# The estimator steps in mde_fit_subset are copied from scripts/donor_de_lib.R
# (.fit_donor_model) with cycle phase as the filter group, exactly as the frozen
# scripts/phase_positive_control_fit_cohort.R did; every numerical setting is read from
# config/statistical_analysis.json. Sourcing this file runs nothing.
suppressPackageStartupMessages({library(edgeR); library(limma); library(statmod); library(jsonlite)})

mde_sha <- function(path) {
  ans <- system2("sha256sum",shQuote(path),stdout=TRUE,stderr=TRUE)
  if(!is.null(attr(ans,"status")) || length(ans)!=1L || !grepl("^[0-9a-f]{64} ",ans)) stop(paste("cannot hash",path))
  substr(ans,1,64)
}
mde_read_tsv <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,na.strings="NA",quote="",comment.char="")
mde_write_tsv <- function(x,path) write.table(x,path,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,na="NA")
mde_write_json <- function(x,path) {
  tmp <- paste0(path,".partial.",Sys.getpid())
  jsonlite::write_json(x,tmp,pretty=TRUE,auto_unbox=TRUE,null="null",na="null",digits=NA)
  if(!file.rename(tmp,path)) stop("atomic JSON rename failed")
}
mde_parse_opts <- function(args) {
  if(length(args)%%2!=0) stop("arguments must be --key value pairs")
  opts <- list()
  if(length(args)) for(i in seq(1,length(args),by=2)) {
    k <- substring(args[i],3)
    if(!startsWith(args[i],"--") || !is.null(opts[[k]])) stop("bad/duplicate argument")
    opts[[k]] <- args[i+1]
  }
  opts
}
mde_software <- function() list(R=as.character(getRversion()),edgeR=as.character(packageVersion("edgeR")),
                                limma=as.character(packageVersion("limma")),statmod=as.character(packageVersion("statmod")),
                                jsonlite=as.character(packageVersion("jsonlite")))
# Freeze verification. Data, configuration and code artifacts must match; documentation
# artifacts are recorded rather than enforced because the analysis host keeps pre-revision
# bytes of published documents (same convention as the positive-control scripts).
mde_verify_freeze <- function(path,phase) {
  fr <- jsonlite::fromJSON(path,simplifyVector=FALSE)
  if(!identical(fr$phase,phase) || length(fr$artifacts)==0L) stop(paste("wrong freeze phase/schema:",path))
  exempt <- unlist(fr$host_verification_exempt)
  check <- list()
  for(p in names(fr$artifacts)) {
    ok <- file.exists(p) && identical(mde_sha(p),fr$artifacts[[p]])
    check[[p]] <- ok
    if(!ok && !startsWith(p,"docs/") && !(p %in% exempt)) stop(paste("frozen artifact changed or missing:",p))
  }
  list(freeze=fr,check=check)
}
# Hash chain of a saved fit: either a pinned object hash, or a pinned provenance hash whose
# provenance records the object hash. Returns the verified chain for the provenance record.
mde_verify_fit_chain <- function(entry) {
  if(!is.null(entry$model_sha256)) {
    expected <- entry$model_sha256; via <- entry$model_sha256_source
  } else {
    if(is.null(entry$provenance) || is.null(entry$provenance_sha256)) stop("fit entry lacks a hash chain")
    actual_prov <- mde_sha(entry$provenance)
    if(!identical(actual_prov,entry$provenance_sha256)) stop(paste("provenance hash differs from the pinned value:",entry$provenance))
    pv <- jsonlite::fromJSON(entry$provenance,simplifyVector=FALSE)
    expected <- pv$output_sha256[[entry$model]]
    if(is.null(expected)) stop("provenance does not record the model object hash")
    via <- paste0(entry$provenance_sha256_source," -> ",entry$provenance,":output_sha256")
  }
  actual <- mde_sha(entry$model)
  if(!identical(actual,expected)) stop(paste("model object hash differs from the pinned chain:",entry$model))
  list(model=entry$model,sha256=actual,pinned_via=via)
}
# Two-sided rejection probability of a t statistic with realized SE and df under true effect delta.
mde_power <- function(delta,se,df,tcrit) {
  ncp <- delta/se
  pt(-tcrit,df,ncp=ncp)+1-pt(tcrit,df,ncp=ncp)
}
# Exact MDE: the delta whose rejection probability equals `power` at level alpha, per gene.
# NA where SE or df is not finite and positive. Also returns the closed-form approximation.
mde_solve <- function(se,df,alpha,power=0.8,tol=1e-8) {
  n <- length(se); out <- rep(NA_real_,n); approx <- rep(NA_real_,n)
  if(!is.finite(alpha) || alpha<=0 || alpha>=1) stop("alpha must be in (0,1)")
  ok <- is.finite(se)&se>0&is.finite(df)&df>0
  tcrit <- qt(1-alpha/2,df)
  approx[ok] <- (tcrit[ok]+qt(power,df[ok]))*se[ok]
  warn <- 0L; failed <- 0L
  for(i in which(ok)) {
    f <- function(d) mde_power(d,se[i],df[i],tcrit[i])-power
    upper <- 2*approx[i]; k <- 0L
    repeat { v <- suppressWarnings(f(upper)); if(is.finite(v) && v>0) break; upper <- upper*2; k <- k+1L; if(k>20L) break }
    if(k>20L) { failed <- failed+1L; next }
    r <- tryCatch(withCallingHandlers(uniroot(f,c(0,upper),tol=tol,maxiter=1000)$root,
                   warning=function(w){warn<<-warn+1L;invokeRestart("muffleWarning")}),error=function(e)NA_real_)
    if(is.na(r)) failed <- failed+1L
    out[i] <- r
  }
  list(mde=out,approx=approx,precision_warnings=warn,failed=failed)
}
# Design for a (sub)cohort of the phase positive control. `primary` mirrors the frozen
# positive-control design (phase, disease, centered age, batch treatment contrasts over the
# sorted levels present, adenomyosis); `reduced` omits the batch contrasts. Constant nuisance
# columns are dropped with a receipt; rank deficiency or residual df below the minimum is an error.
mde_rebuild_design <- function(clinical,cohort,cfg,pcfg,variant=c("primary","reduced")) {
  variant <- match.arg(variant)
  levels_phase <- c(pcfg$phase_factor$reference,pcfg$phase_factor$case)
  batch_col <- cfg$cohorts[[cohort]]$batch_column
  required <- c("donor_id","phase","disease",cfg$design$age_column,cfg$design$adenomyosis_column,batch_col)
  if(!all(required %in% names(clinical))) stop("missing clinical design field")
  if(anyDuplicated(clinical$donor_id) || anyNA(clinical[,required])) stop("duplicate donor IDs or missing covariate")
  if(any(!clinical$phase %in% levels_phase) || length(unique(clinical$phase))!=2) stop("both phase groups required")
  age <- as.numeric(clinical[[cfg$design$age_column]]); adeno <- as.numeric(clinical[[cfg$design$adenomyosis_column]])
  if(any(!is.finite(age)) || any(!adeno %in% c(0,1))) stop("invalid age/adenomyosis")
  batch <- as.character(clinical[[batch_col]]); if(any(!nzchar(batch))) stop("missing batch")
  levels <- sort(unique(batch))
  cross <- table(factor(batch,levels=levels),factor(clinical$phase,levels=levels_phase))
  design <- cbind(intercept=rep(1,nrow(clinical)),SEC_minus_PRO=as.numeric(clinical$phase==pcfg$phase_factor$case))
  dropped <- character()
  disease <- as.numeric(clinical$disease=="ENDO")
  if(length(unique(disease))==1) dropped <- c(dropped,"disease") else design <- cbind(design,ENDO_minus_non_ENDO=disease)
  if(length(unique(age))==1) dropped <- c(dropped,"age") else design <- cbind(design,age_centered_years=age-mean(age))
  if(variant=="primary") {
    if(!any(rowSums(cross>0)==2)) stop("no mixed phase batch support")
    if(length(levels)==1) dropped <- c(dropped,batch_col) else {
      contrasts <- model.matrix(~factor(batch,levels=levels))[,-1,drop=FALSE]
      colnames(contrasts) <- paste0(batch_col,":",levels[-1],"_vs_",levels[1])
      design <- cbind(design,contrasts)
    }
  }
  if(length(unique(adeno))==1) dropped <- c(dropped,"adenomyosis") else design <- cbind(design,adenomyosis=adeno)
  rownames(design) <- clinical$donor_id
  rank <- qr(design,tol=cfg$fit$rank_tolerance)$rank
  if(rank!=ncol(design)) stop("design not full rank")
  df <- nrow(design)-rank
  if(df<cfg$minimum_residual_df) stop("residual df below minimum")
  list(design=design,receipt=list(variant=variant,rows=nrow(design),columns=ncol(design),rank=rank,residual_df=df,
       age_center=mean(age),batch_levels=levels,n_batch_levels=length(levels),dropped_constant_nuisance=dropped,
       phase_batch_table=unclass(cross)))
}
# The v1 estimator on a (sub)cohort: expression filter with the phase group, TMM, voom, lmFit,
# robust eBayes, BH over the subset family. Returns the table with the columns the frozen
# positive-control tables carry, the filter mask and the limma fit object.
mde_fit_subset <- function(counts,design,group,cfg,coefficient) {
  if(!is.matrix(counts) || !is.numeric(counts) || nrow(counts)<2L || any(!is.finite(counts)) || any(counts<0) || any(colSums(counts)<=0)) stop("invalid counts/library")
  if(!identical(colnames(counts),rownames(design)) || anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts))) stop("counts/design donor order invalid")
  warnings <- character()
  res <- withCallingHandlers({
    f <- cfg$filtering
    raw_lib <- colSums(counts)
    support_n <- min(table(group))
    if(support_n>f$large_n) support_n <- f$large_n+(support_n-f$large_n)*f$min_prop
    keep <- edgeR::filterByExpr(counts,group=group,lib.size=colSums(counts),min.count=f$min_count,
                                min.total.count=f$min_total_count,large.n=f$large_n,min.prop=f$min_prop)
    if(sum(keep)<2L) stop("too few expression-pass genes")
    y <- edgeR::DGEList(counts=counts[keep,,drop=FALSE],group=group)
    if(any(y$samples$lib.size<=0)) stop("zero filtered library")
    n <- cfg$normalization
    y <- edgeR::calcNormFactors(y,method=n$method,refColumn=NULL,logratioTrim=n$logratio_trim,sumTrim=n$sum_trim,
                                doWeighting=n$do_weighting,Acutoff=n$A_cutoff,p=n$p)
    v <- cfg$voom
    va <- list(counts=y,design=design,lib.size=NULL,normalize.method=v$normalize_method,block=NULL,correlation=NULL,
               weights=NULL,span=v$span,plot=v$plot,save.plot=v$save_plot)
    adaptive <- "adaptive.span" %in% names(formals(limma::voom))
    if(adaptive) va$adaptive.span <- v$adaptive_span else if(as.character(packageVersion("limma"))!="3.62.1") stop("unreviewed voom implementation without adaptive.span")
    vo <- do.call(limma::voom,va)
    if(!is.matrix(vo$weights) || any(!is.finite(vo$weights)) || any(vo$weights<=0)) stop("voom precision weights must be finite and positive")
    fit <- limma::lmFit(vo,design=design,method=cfg$fit$method)
    e <- cfg$ebayes
    ea <- list(fit=fit,proportion=e$proportion,stdev.coef.lim=e$stdev_coef_lim,trend=e$trend,robust=e$robust,winsor.tail.p=e$winsor_tail_p)
    legacy <- "legacy" %in% names(formals(limma::eBayes))
    if(legacy) ea$legacy <- e$legacy else if(as.character(packageVersion("limma"))!="3.62.1") stop("unreviewed eBayes implementation without legacy argument")
    fit <- do.call(limma::eBayes,ea)
    co <- match(coefficient,colnames(design)); if(is.na(co)) stop("coefficient missing from design")
    logfc <- fit$coefficients[,co]; statistic <- fit$t[,co]; p <- fit$p.value[,co]
    se <- sqrt(fit$s2.post)*fit$stdev.unscaled[,co]; degrees <- rep_len(fit$df.total,length(logfc))
    finite <- is.finite(logfc)&is.finite(statistic)&is.finite(p)&is.finite(se)&se>0&is.finite(degrees)&degrees>0&
              is.finite(fit$df.residual)&fit$df.residual>=cfg$minimum_residual_df
    q <- rep(NA_real_,length(p)); q[finite] <- p.adjust(p[finite],method=cfg$testing$adjust_method,n=length(p))
    critical <- qt(1-(1-cfg$testing$confidence_level)/2,df=degrees)
    list(table=data.frame(gene_id=rownames(y),logFC=logfc,AveExpr=fit$Amean,t=statistic,P.Value=p,adj.P.Val=q,
                          moderated_SE=se,df_residual=fit$df.residual,df_total=degrees,CI_low=logfc-critical*se,CI_high=logfc+critical*se,
                          direction=sign(logfc),fit_estimable=finite,row.names=NULL,check.names=FALSE),
         keep=keep,fit=fit,family_size=sum(keep),
         filter_audit=list(cpm_cutoff=f$min_count/median(raw_lib)*1e6,minimum_support_donors=support_n,
                           min_total_count=f$min_total_count,median_raw_library_size=median(raw_lib)))
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  res$warnings <- warnings
  res
}
mde_read_counts <- function(path) {
  ct <- mde_read_tsv(path)
  if(names(ct)[1]!="gene_id") stop("counts first column must be gene_id")
  counts <- as.matrix(ct[,-1,drop=FALSE]); storage.mode(counts) <- "double"; rownames(counts) <- ct$gene_id
  counts
}
