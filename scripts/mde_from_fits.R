#!/usr/bin/env Rscript
# Minimum detectable effect analysis, Part A (docs/MDE_ANALYSIS_PROTOCOL.md section 3):
# per-gene MDE at fixed significance scenarios from a saved limma fit, CI-containment counts
# and expression-stratified summaries. Reads the fit object only; fits nothing; writes nothing
# outside --output-dir. Real runs require the frozen rules and verify the fit's hash chain.
source("scripts/mde_lib.R")
opts <- mde_parse_opts(commandArgs(trailingOnly=TRUE))
synthetic <- identical(opts[["synthetic-test"]],"true")
required <- c("fit-name","config","mde-config","output-dir")
if(!synthetic) required <- c(required,"freeze")
if(!all(required %in% names(opts))) stop(paste("required:",paste(required,collapse=" ")))
cfg <- jsonlite::fromJSON(opts$config,simplifyVector=TRUE)
mcfg <- jsonlite::fromJSON(opts[["mde-config"]],simplifyVector=TRUE)
fit_name <- opts[["fit-name"]]
started <- format(Sys.time(),tz="UTC",usetz=TRUE)
freeze_check <- NULL; chain <- NULL
if(synthetic) {
  for(k in c("model","coefficient")) if(is.null(opts[[k]])) stop("synthetic mode requires --model and --coefficient")
  if(!startsWith(opts$model,"work/synthetic_mde/")) stop("synthetic mode requires a work/synthetic_mde/ fixture")
  model_path <- opts$model; coefficient <- opts$coefficient; cohort <- "synthetic"
} else {
  fz <- mde_verify_freeze(opts$freeze,"MDE_ANALYSIS_RULES_FROZEN_PRE_COMPUTATION"); freeze_check <- fz$check
  if(!identical(fz$freeze$artifacts[["config/mde_analysis.json"]],mde_sha(opts[["mde-config"]]))) stop("mde config is not the frozen one")
  entry <- mcfg$fits[[fit_name]]
  if(is.null(entry)) stop(paste("unknown fit entry:",fit_name))
  chain <- mde_verify_fit_chain(entry)
  model_path <- entry$model; coefficient <- entry$coefficient; cohort <- entry$cohort
}
out <- opts[["output-dir"]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("output directory is not empty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

m <- readRDS(model_path)
fit <- m$fit; tab <- m$table
if(is.null(fit) || is.null(tab) || !inherits(fit,"MArrayLM")) stop("model object lacks the limma fit or the table")
co <- match(coefficient,colnames(fit$design))
if(is.na(co)) stop(paste("coefficient not in the fit design:",coefficient))
n <- nrow(tab)
if(!identical(as.character(tab$gene_id),rownames(fit$coefficients)) || nrow(fit$coefficients)!=n) stop("table and fit gene order differ")
se <- sqrt(fit$s2.post)*fit$stdev.unscaled[,co]; df <- rep_len(fit$df.total,n)
if(max(abs(se-tab$moderated_SE),na.rm=TRUE)>1e-10 || max(abs(df-tab$df_total),na.rm=TRUE)>1e-10 ||
   max(abs(fit$coefficients[,co]-tab$logFC),na.rm=TRUE)>1e-10) stop("fit object and table disagree on SE, df or logFC")
est <- as.logical(tab$fit_estimable)
m_family <- n
power <- mcfg$power
alphas <- c(mcfg$nominal_alpha,mcfg$bh_q*mcfg$bh_scenarios_k/m_family)
names(alphas) <- c("alpha_0.05",paste0("bh_k",mcfg$bh_scenarios_k))
if(!identical(names(alphas),mcfg$scenario_names)) stop("scenario names differ from the frozen configuration")
if(any(alphas>=1)) stop("a scenario alpha is not below 1 for this family size")

per <- data.frame(gene_id=tab$gene_id,AveExpr=tab$AveExpr,logFC=tab$logFC,moderated_SE=se,df_total=df,abs_t=abs(tab$t),
                  P.Value=tab$P.Value,adj.P.Val=tab$adj.P.Val,CI_low=tab$CI_low,CI_high=tab$CI_high,fit_estimable=est,
                  stringsAsFactors=FALSE,check.names=FALSE)
solver_notes <- list()
for(s in names(alphas)) {
  sol <- mde_solve(ifelse(est,se,NA_real_),df,alphas[[s]],power=power,tol=mcfg$mde_solver$tolerance)
  per[[paste0("mde_",s)]] <- sol$mde
  if(s=="alpha_0.05") per[["mde_approx_alpha_0.05"]] <- sol$approx
  solver_notes[[s]] <- list(alpha=alphas[[s]],precision_warnings=sol$precision_warnings,failed=sol$failed,
                            max_abs_exact_minus_approx=max(abs(sol$mde-sol$approx),na.rm=TRUE))
}
# family-wide Bonferroni interval
crit_b <- qt(1-(1-mcfg$ci_levels$family_bonferroni)/(2*m_family),df=df)
per$CI_low_bonferroni <- tab$logFC-crit_b*se; per$CI_high_bonferroni <- tab$logFC+crit_b*se

qs <- mcfg$mde_summary_quantiles
quant <- function(x) { v <- x[est & is.finite(x)]; if(!length(v)) return(NULL); setNames(as.list(quantile(v,qs,type=7,names=FALSE)),paste0("q",qs)) }
scenario_summary <- list()
for(s in names(alphas)) {
  x <- per[[paste0("mde_",s)]]
  qq <- quant(x)
  scenario_summary[[s]] <- list(alpha=alphas[[s]],quantiles_log2=qq,quantiles_fold=lapply(qq,function(v)2^v),
    share_at_most=setNames(lapply(mcfg$mde_at_most_log2,function(d)mean(x[est]<=d,na.rm=TRUE)),paste0("le_",mcfg$mde_at_most_log2)),
    n_estimable=sum(est),n_solved=sum(est & is.finite(x)))
}
contain <- function(lo,hi,d) sum(est & is.finite(lo) & is.finite(hi) & lo>-d & hi<d)
ci_contain <- list(per_gene=setNames(lapply(mcfg$ci_containment_log2_thresholds,function(d)contain(per$CI_low,per$CI_high,d)),paste0("within_",mcfg$ci_containment_log2_thresholds)),
                   family_bonferroni=setNames(lapply(mcfg$ci_containment_log2_thresholds,function(d)contain(per$CI_low_bonferroni,per$CI_high_bonferroni,d)),paste0("within_",mcfg$ci_containment_log2_thresholds)),
                   n_estimable=sum(est),per_gene_level=mcfg$ci_levels$per_gene,family_level=mcfg$ci_levels$family_bonferroni,
                   bonferroni_per_gene_alpha=(1-mcfg$ci_levels$family_bonferroni)/m_family)
# expression bins
breaks <- quantile(per$AveExpr[est],probs=seq(0,1,length.out=mcfg$expression_bins+1),type=7,names=FALSE)
bin <- rep(NA_integer_,n); bin[est] <- as.integer(cut(per$AveExpr[est],breaks=breaks,include.lowest=TRUE,labels=FALSE))
per$expression_bin <- bin
by_bin <- do.call(rbind,lapply(seq_len(mcfg$expression_bins),function(b) {
  ix <- which(bin==b)
  row <- data.frame(expression_bin=b,AveExpr_min=min(per$AveExpr[ix]),AveExpr_max=max(per$AveExpr[ix]),n_genes=length(ix),
                    median_SE=median(per$moderated_SE[ix]),median_df_total=median(per$df_total[ix]),check.names=FALSE)
  for(s in names(alphas)) row[[paste0("median_mde_",s)]] <- median(per[[paste0("mde_",s)]][ix],na.rm=TRUE)
  row
}))
sig <- est & !is.na(tab$adj.P.Val) & tab$adj.P.Val<0.05
observed <- list(genes_BH_q_lt_0_05=sum(sig),max_abs_t=max(abs(tab$t[est]),na.rm=TRUE),
                 median_abs_logFC_all=median(abs(tab$logFC[est])),
                 abs_logFC_quantiles_among_significant=if(sum(sig)) setNames(as.list(quantile(abs(tab$logFC[sig]),qs,type=7,names=FALSE)),paste0("q",qs)) else NULL,
                 SE_quantiles=setNames(as.list(quantile(se[est],qs,type=7,names=FALSE)),paste0("q",qs)),
                 df_total_range=range(df[est]),df_residual=unique(fit$df.residual),
                 s2_prior=fit$s2.prior,df_prior_range=range(fit$df.prior))
summary <- list(fit_name=fit_name,cohort=cohort,coefficient=coefficient,synthetic=synthetic,family_size=m_family,n_estimable=sum(est),
                power=power,two_sided=TRUE,scenarios=scenario_summary,ci_containment=ci_contain,observed=observed,
                solver=solver_notes,expression_bins=mcfg$expression_bins)
mde_write_tsv(per,file.path(out,"mde_per_gene.tsv"))
mde_write_tsv(by_bin,file.path(out,"mde_by_expression_bin.tsv"))
mde_write_json(summary,file.path(out,"summary.json"))
figure <- tryCatch({
  png(file.path(out,"mde_ecdf.png"),width=1400,height=900,res=160)
  cols <- c("#1b4f72","#2e86c1","#48c9b0","#f39c12","#c0392b")
  plot(NA,xlim=c(0,6),ylim=c(0,1),xlab="minimum detectable effect (log2 fold change, 80% power)",ylab="share of tested genes",
       main=paste0(fit_name,": ",coefficient," (",m_family," genes)"))
  for(i in seq_along(alphas)) { x <- per[[paste0("mde_",names(alphas)[i])]][est]; x <- x[is.finite(x)]; lines(ecdf(x),col=cols[i],lwd=2,do.points=FALSE,verticals=TRUE) }
  abline(v=c(0.585,1,2),lty=3,col="grey50")
  legend("bottomright",legend=sprintf("%s (alpha=%.2e)",names(alphas),alphas),col=cols,lwd=2,bty="n",cex=0.8)
  dev.off(); "mde_ecdf.png"
},error=function(e){try(dev.off(),silent=TRUE);paste("figure_not_rendered:",conditionMessage(e))})
outputs <- setdiff(list.files(out,full.names=TRUE),file.path(out,"provenance.json"))
prov <- list(started_utc=started,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),analysis="MDE_PART_A_FROM_SAVED_FIT",
             protocol=mcfg$protocol,fit_name=fit_name,cohort=cohort,coefficient=coefficient,synthetic=synthetic,
             data_kind=if(synthetic)"synthetic_test_only" else "saved_public_cohort_fit",
             model=list(path=model_path,sha256=mde_sha(model_path),chain=chain),
             freeze=if(synthetic) NULL else list(path=opts$freeze,sha256=mde_sha(opts$freeze),phase="MDE_ANALYSIS_RULES_FROZEN_PRE_COMPUTATION",artifact_check=freeze_check),
             input_sha256=setNames(lapply(c(opts$config,opts[["mde-config"]],"scripts/mde_lib.R","scripts/mde_from_fits.R","scripts/statistical_runtime_env.sh"),mde_sha),
                                   c(opts$config,opts[["mde-config"]],"scripts/mde_lib.R","scripts/mde_from_fits.R","scripts/statistical_runtime_env.sh")),
             family_size=m_family,scenario_alphas=as.list(alphas),power=power,figure=figure,software=mde_software(),
             output_sha256=setNames(lapply(outputs,mde_sha),outputs),
             no_model_fitted=TRUE,v1_artifacts_modified=FALSE,
             scientific_scope=if(synthetic)"invented fixture, no biological content" else "design sensitivity under realized variance estimates; not evidence about any gene")
mde_write_json(prov,file.path(out,"provenance.json"))
cat("MDE_FROM_FITS_COMPLETE",fit_name,"family",m_family,"estimable",sum(est),
    "median_mde_alpha_0.05",signif(scenario_summary$alpha_0.05$quantiles_log2$q0.5,4),
    "median_mde_bh_k100",signif(scenario_summary$bh_k100$quantiles_log2$q0.5,4),"\n")
