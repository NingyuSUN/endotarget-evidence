#!/usr/bin/env Rscript
# Synthetic fixture for the MDE analysis code: invented negative-binomial counts for a
# discovery-shaped phase cohort (ledger columns of the positive-control discovery ledger),
# a planted phase program, and a full-cohort "frozen" model object produced with the same
# estimator. Everything lives under work/synthetic_mde/ and carries no biological content.
source("scripts/mde_lib.R")
args <- mde_parse_opts(commandArgs(trailingOnly=TRUE))
cfg <- jsonlite::fromJSON(if(is.null(args$config)) "config/statistical_analysis.json" else args$config,simplifyVector=TRUE)
pcfg <- jsonlite::fromJSON(if(is.null(args[["phase-config"]])) "config/phase_positive_control.json" else args[["phase-config"]],simplifyVector=TRUE)
out <- "work/synthetic_mde"; dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261027)
G <- 800L; n_sec <- 10L; n_pro <- 12L; N <- n_sec+n_pro
genes <- sprintf("ENSG%011d",seq_len(G))
donors <- sprintf("syn%02d",seq_len(N))
phase <- c(rep(pcfg$phase_factor$case,n_sec),rep(pcfg$phase_factor$reference,n_pro))
disease <- sample(c("ENDO","non-ENDO"),N,replace=TRUE); disease[1:2] <- c("ENDO","non-ENDO"); disease[n_sec+1:2] <- c("ENDO","non-ENDO")
age <- sample(25:42,N,replace=TRUE)
batch <- sample(c("batchA","batchB","batchC"),N,replace=TRUE); batch[c(1,n_sec+1)] <- "batchA"; batch[c(2,n_sec+2)] <- "batchB"
adeno <- rbinom(N,1,0.25); if(length(unique(adeno))==1) adeno[1] <- 1-adeno[1]
mu <- exp(rnorm(G,mean=4,sd=1.8)); lib <- exp(rnorm(N,mean=0,sd=0.25))
effect <- numeric(G); planted <- sample(G,100); effect[planted] <- sample(c(-1,1),100,replace=TRUE)*runif(100,0.8,3)
counts <- matrix(0,G,N,dimnames=list(genes,donors))
for(j in seq_len(N)) {
  m <- mu*lib[j]*2^(effect*(phase[j]==pcfg$phase_factor$case))
  counts[,j] <- rnbinom(G,mu=m,size=1/0.15)
}
counts <- counts[rowSums(counts)>0,,drop=FALSE]
clinical <- data.frame(donor_id=donors,phase=phase,sub_phase=ifelse(phase==pcfg$phase_factor$case,"mid","proliferative"),disease=disease,age=age,
                       capture_batch=batch,sequencing_batch="first",adenomyosis=adeno,primary_lineage_cells=sample(200:900,N,replace=TRUE),stringsAsFactors=FALSE)
ct <- data.frame(gene_id=rownames(counts),counts,check.names=FALSE)
gz <- gzfile(file.path(out,"counts.tsv.gz"),"w"); write.table(ct,gz,sep="\t",row.names=FALSE,quote=FALSE); close(gz)
mde_write_tsv(clinical,file.path(out,"clinical.tsv"))
built <- mde_rebuild_design(clinical,"discovery",cfg,pcfg,"primary")
fit <- mde_fit_subset(counts,built$design,factor(clinical$phase,levels=c(pcfg$phase_factor$reference,pcfg$phase_factor$case)),cfg,pcfg$coefficient)
model <- list(table=fit$table,keep=fit$keep,filter_audit=fit$filter_audit,fit=fit$fit,audit=built$receipt,synthetic=TRUE)
saveRDS(model,file.path(out,"model.rds"),version=3,compress="gzip")
mde_write_json(list(synthetic=TRUE,genes=nrow(counts),donors=N,planted_effects=100,family_size=fit$family_size,
                    n_sig=sum(fit$table$adj.P.Val<0.05,na.rm=TRUE),design_columns=colnames(built$design),
                    note="invented data for code testing only; never a scientific result"),file.path(out,"fixture.json"))
cat("SYNTHETIC_FIXTURE_WRITTEN genes",nrow(counts),"family",fit$family_size,"sig",sum(fit$table$adj.P.Val<0.05,na.rm=TRUE),"\n")
