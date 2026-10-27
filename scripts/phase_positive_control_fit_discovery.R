#!/usr/bin/env Rscript
# Phase positive control: discovery SEC_minus_PRO donor model with the v1 estimator.
# Every numerical setting is read from config/statistical_analysis.json unchanged. The only
# differences from scripts/donor_de_lib.R are those the frozen protocol specifies: the
# coefficient (SEC_minus_PRO) and the grouping variable used by filterByExpr and the group
# minimum (cycle phase instead of disease). Run from the project root in the statistical runtime.
suppressPackageStartupMessages({library(edgeR); library(limma); library(statmod); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)%%2!=0) stop("arguments must be --key value pairs")
opts <- list()
if(length(args)) for(i in seq(1,length(args),by=2)) {
  if(!startsWith(args[i],"--") || !is.null(opts[[substring(args[i],3)]])) stop("bad/duplicate argument")
  opts[[substring(args[i],3)]] <- args[i+1]
}
required <- c("counts","design","clinical","config","phase-config","freeze","output-dir")
if(!all(required %in% names(opts))) stop(paste("required:",paste(required,collapse=" ")))
sha <- function(path) {
  ans <- system2("sha256sum",shQuote(path),stdout=TRUE,stderr=TRUE)
  if(!is.null(attr(ans,"status")) || length(ans)!=1L || !grepl("^[0-9a-f]{64} ",ans)) stop(paste("cannot hash",path))
  substr(ans,1,64)
}
# F0 freeze: data and configuration artifacts must be unchanged. Documentation artifacts are
# recorded rather than enforced here, because the analysis host keeps the pre-revision bytes of
# published documents; they are verified on the published checkout and recorded in the F1 receipt.
fr <- jsonlite::fromJSON(opts$freeze,simplifyVector=FALSE)
if(!identical(fr$phase,"PHASE_POSITIVE_CONTROL_PROTOCOL_FROZEN_PRE_DATA") || length(fr$artifacts)==0L) stop("wrong freeze phase/schema")
freeze_check <- list()
for(p in names(fr$artifacts)) {
  ok <- file.exists(p) && identical(sha(p),fr$artifacts[[p]])
  freeze_check[[p]] <- ok
  if(!ok && !startsWith(p,"docs/")) stop(paste("frozen artifact changed:",p))
}
cfg <- jsonlite::fromJSON(opts$config,simplifyVector=TRUE)
pcfg <- jsonlite::fromJSON(opts[["phase-config"]],simplifyVector=TRUE)
coefficient <- pcfg$coefficient
if(!identical(coefficient,"SEC_minus_PRO")) stop("unexpected coefficient in phase config")
levels_phase <- c(pcfg$phase_factor$reference,pcfg$phase_factor$case)
min_group <- pcfg$cohorts$discovery$minimum_per_phase_group
readtab <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("NA"),quote="",comment.char="")
ct <- readtab(opts$counts)
if(names(ct)[1]!="gene_id") stop("counts first column must be gene_id")
counts <- as.matrix(ct[,-1,drop=FALSE]); storage.mode(counts)<-"double"; rownames(counts) <- ct$gene_id
dt <- readtab(opts$design)
if(names(dt)[1]!="donor_id") stop("design first column must be donor_id")
design <- as.matrix(dt[,-1,drop=FALSE]); storage.mode(design)<-"double"; rownames(design)<-dt$donor_id
clinical <- readtab(opts$clinical)
out <- opts[["output-dir"]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("output directory is not empty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
started <- format(Sys.time(),tz="UTC",usetz=TRUE)

# ---- validation: mirrors validate_donor_inputs with the phase factor as the group ----
if (!is.matrix(counts) || !is.numeric(counts) || nrow(counts)<2L) stop("invalid counts matrix")
if (any(!is.finite(counts)) || any(counts<0)) stop("counts must be finite and nonnegative")
if (isTRUE(pcfg$cohorts$discovery$integer_counts) && any(counts!=floor(counts))) stop("integer count contract failed")
if (is.null(rownames(counts)) || anyDuplicated(rownames(counts)) || any(!grepl("^ENSG[0-9]{11}$",rownames(counts)))) stop("invalid unique stable genes")
if (is.null(colnames(counts)) || anyDuplicated(colnames(counts))) stop("invalid count donor IDs")
if (!all(c("donor_id","phase","disease","primary_lineage_cells") %in% names(clinical)) || anyDuplicated(clinical$donor_id)) stop("invalid clinical ledger")
if (!identical(colnames(counts),as.character(clinical$donor_id)) || !identical(rownames(design),as.character(clinical$donor_id))) stop("exact donor order required")
if (any(!clinical$phase %in% levels_phase)) stop("unknown phase label")
group <- factor(clinical$phase,levels=levels_phase)
ng <- table(group)
if (any(ng<min_group)) stop("phase group donor minimum failed")
if (any(!is.finite(clinical$primary_lineage_cells)) || any(clinical$primary_lineage_cells<cfg$design$min_cells)) stop("donor cell threshold failed")
if (any(colSums(counts)<=0)) stop("zero count library")
if (!is.matrix(design) || !is.numeric(design) || any(!is.finite(design))) stop("invalid numeric design")
if (sum(colnames(design)==coefficient)!=1L || !identical(as.numeric(design[,coefficient]),as.numeric(group==pcfg$phase_factor$case))) stop("phase coefficient mismatch")
if (sum(colnames(design)=="ENDO_minus_non_ENDO")!=1L || !identical(as.numeric(design[,"ENDO_minus_non_ENDO"]),as.numeric(clinical$disease=="ENDO"))) stop("disease covariate mismatch")
rank <- qr(design,tol=cfg$fit$rank_tolerance)$rank
if (rank!=ncol(design)) stop("design is not full rank")
df <- nrow(design)-rank
if (df<cfg$minimum_residual_df) stop("insufficient residual degrees of freedom")
audit <- list(n_donors=ncol(counts),input_genes=nrow(counts),rank=rank,residual_df=df,
              proliferative=unname(ng[[1]]),secretory=unname(ng[[2]]),group_column="phase",coefficient=coefficient,
              design_columns=colnames(design))

# ---- estimator: identical to .fit_donor_model in scripts/donor_de_lib.R ----
warnings <- character()
model <- withCallingHandlers({
  f <- cfg$filtering
  raw_lib <- colSums(counts)
  support_n <- min(ng)
  if (support_n>f$large_n) support_n <- f$large_n+(support_n-f$large_n)*f$min_prop
  cpm_cutoff <- f$min_count/median(raw_lib)*1e6
  filter_audit <- list(cpm_cutoff=cpm_cutoff,minimum_support_donors=support_n,
                       min_total_count=f$min_total_count,median_raw_library_size=median(raw_lib))
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
  if (!adaptive_supported && as.character(packageVersion("limma"))!="3.62.1")
    stop("unreviewed voom implementation without adaptive.span")
  vo <- do.call(limma::voom,va)
  if (!is.matrix(vo$weights) || any(!is.finite(vo$weights)) || any(vo$weights<=0))
    stop("voom precision weights must be finite and positive")
  fit <- limma::lmFit(vo,design=design,method=cfg$fit$method)
  e <- cfg$ebayes
  ea <- list(fit=fit,proportion=e$proportion,stdev.coef.lim=e$stdev_coef_lim,
             trend=e$trend,robust=e$robust,winsor.tail.p=e$winsor_tail_p)
  legacy_supported <- "legacy" %in% names(formals(limma::eBayes))
  if (legacy_supported) ea$legacy <- e$legacy
  if (!legacy_supported && as.character(packageVersion("limma"))!="3.62.1")
    stop("unreviewed eBayes implementation without legacy argument")
  fit <- do.call(limma::eBayes,ea)
  co <- match(coefficient,colnames(design))
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
            null_matching_logCPM=rowMeans(edgeR::cpm(y,log=TRUE,prior.count=.5,normalized.lib.sizes=TRUE)),
            row.names=NULL,check.names=FALSE)
  norm <- data.frame(donor_id=colnames(counts),raw_library_size=colSums(counts),
             filtered_library_size=y$samples$lib.size,norm_factor=y$samples$norm.factors,
             effective_library_size=y$samples$lib.size*y$samples$norm.factors,row.names=NULL)
  list(table=tab,keep=keep,filter_audit=filter_audit,normalization=norm,voom=vo,fit=fit,
       settings=list(adaptive_span_argument_supported=adaptive_supported,effective_adaptive_span=FALSE,
                     legacy_argument_supported=legacy_supported,effective_legacy_ebayes=TRUE))
}, warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })

write_tsv <- function(x,name) write.table(x,file=file.path(out,name),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,na="NA")
write_tsv(model$table,"de_all.tsv")
write_tsv(model$table[model$table$fit_estimable,,drop=FALSE],"de_fitted.tsv")
write_tsv(data.frame(gene_id=rownames(counts),expression_pass=model$keep),"expression_filter.tsv")
write_tsv(model$normalization,"normalization.tsv")
write_tsv(dt,"design.tsv")
saveRDS(model,file=file.path(out,"model.rds"),version=3,compress="gzip")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
input_paths <- c(opts$counts,opts$design,opts$clinical,opts$config,opts[["phase-config"]],opts$freeze,
                 "scripts/phase_positive_control_fit_discovery.R","scripts/donor_de_lib.R")
outputs <- list.files(out,full.names=TRUE)
sig <- model$table$adj.P.Val<.05
prov <- list(started_utc=started,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
             cohort="discovery",contrast=coefficient,group_column="phase",data_kind="public_cohort",
             scientific_effect_run=TRUE,protocol="docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md",
             freeze_phase=fr$phase,freeze_artifact_check=freeze_check,audit=audit,
             tested_family_size=nrow(model$table),estimable_genes=sum(model$table$fit_estimable),
             significant_BH05=sum(sig,na.rm=TRUE),
             significant_BH05_up=sum(sig & model$table$logFC>0,na.rm=TRUE),
             significant_BH05_down=sum(sig & model$table$logFC<0,na.rm=TRUE),
             effective_compatibility=model$settings,filter_audit=model$filter_audit,
             captured_warnings=warnings,
             software=list(R=R.version.string,edgeR=as.character(packageVersion("edgeR")),
                           limma=as.character(packageVersion("limma")),statmod=as.character(packageVersion("statmod"))),
             reference_implementation="scripts/donor_de_lib.R (hash recorded in input_sha256; not sourced)",
             input_sha256=setNames(lapply(input_paths,sha),input_paths),
             output_sha256=setNames(lapply(outputs,sha),outputs),
             heca_phase_effects_computed=FALSE,
             scientific_scope="observational donor-level phase association in a positive-control design; not a disease effect, not causal")
jsonlite::write_json(prov,file.path(out,"provenance.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
cat("PHASE_DE_COMPLETE tested",nrow(model$table),"estimable",sum(model$table$fit_estimable),"BH05",prov$significant_BH05,"\n")
