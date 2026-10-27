#!/usr/bin/env Rscript
# Synthetic fixture for the secretory program analysis code: two invented cohorts shaped like
# the positive-control ledgers (discovery-like capture batches, HECA-like library lanes), a
# shared planted phase program, a planted attenuation of that program in the ENDO secretory
# donors of the first cohort only, and "frozen" phase-fit tables produced with the same
# estimator. Everything lives under work/synthetic_secretory/ and carries no biological content.
source("scripts/secretory_program_lib.R")
args <- mde_parse_opts(commandArgs(trailingOnly=TRUE))
cfg <- jsonlite::fromJSON(if(is.null(args$config)) "config/statistical_analysis.json" else args$config,simplifyVector=TRUE)
pcfg <- jsonlite::fromJSON(if(is.null(args[["phase-config"]])) "config/phase_positive_control.json" else args[["phase-config"]],simplifyVector=TRUE)
out <- "work/synthetic_secretory"; dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261027)
PRO <- pcfg$phase_factor$reference; SEC <- pcfg$phase_factor$case
G <- 900L; genes <- sprintf("ENSG%011d",seq_len(G))
mu <- exp(rnorm(G,mean=4,sd=1.8))
planted <- sample(G,150); effect <- numeric(G); effect[planted] <- sample(c(-1,1),150,replace=TRUE)*runif(150,0.8,3)
make_cohort <- function(key,n_pro_endo,n_pro_non,n_sec_endo,n_sec_non,batches,batch_col,sub_col,sub_levels,pro_label,depth_col,scale,attenuation) {
  N <- n_pro_endo+n_pro_non+n_sec_endo+n_sec_non
  donors <- sprintf("%s%03d",toupper(substr(key,1,3)),seq_len(N))
  phase <- c(rep(PRO,n_pro_endo+n_pro_non),rep(SEC,n_sec_endo+n_sec_non))
  disease <- c(rep("ENDO",n_pro_endo),rep("non-ENDO",n_pro_non),rep("ENDO",n_sec_endo),rep("non-ENDO",n_sec_non))
  age <- sample(22:44,N,replace=TRUE); adeno <- rbinom(N,1,0.2); if(length(unique(adeno))==1) adeno[1] <- 1-adeno[1]
  batch <- rep_len(batches,N)[sample(N)]; batch[1] <- batches[1]; batch[n_pro_endo+n_pro_non+1] <- batches[1]
  sub <- ifelse(phase==SEC,sample(sub_levels,N,replace=TRUE),pro_label)
  depth <- sample(c(120,300,600,900,2000,5000),N,replace=TRUE); depth[c(1,2,n_pro_endo+n_pro_non+1,N)] <- c(900,2000,5000,600)
  lib <- exp(rnorm(N,mean=0,sd=0.25))
  counts <- matrix(0,G,N,dimnames=list(genes,donors))
  for(j in seq_len(N)) {
    shift <- effect*scale*(phase[j]==SEC)*(1-attenuation*(disease[j]=="ENDO"))
    counts[,j] <- rnbinom(G,mu=mu*lib[j]*2^shift,size=1/0.15)
  }
  counts <- counts[rowSums(counts)>0,,drop=FALSE]
  clinical <- data.frame(donor_id=donors,phase=phase,stringsAsFactors=FALSE)
  clinical[[sub_col]] <- sub; clinical$disease <- disease; clinical$age <- age; clinical[[batch_col]] <- batch
  if(key=="discovery") clinical$sequencing_batch <- "first"
  clinical$adenomyosis <- adeno; clinical[[depth_col]] <- depth
  d <- file.path(out,key); dir.create(d,showWarnings=FALSE)
  ct <- data.frame(gene_id=rownames(counts),counts,check.names=FALSE)
  gz <- gzfile(file.path(d,"counts.tsv.gz"),"w"); write.table(ct,gz,sep="\t",row.names=FALSE,quote=FALSE); close(gz)
  mde_write_tsv(clinical,file.path(d,"clinical.tsv"))
  group <- factor(phase,levels=c(PRO,SEC))
  built <- mde_rebuild_design(clinical,key,cfg,pcfg,"primary")
  fit <- mde_fit_subset(counts,built$design,group,cfg,pcfg$coefficient)
  fn <- sp_filter_norm(counts,group,cfg)
  if(!identical(unname(fn$keep),unname(fit$keep))) stop("fixture filter mismatch")
  mde_write_tsv(fit$table,file.path(d,"de_all.tsv"))
  mde_write_tsv(data.frame(gene_id=rownames(counts),expression_pass=unname(fit$keep)),file.path(d,"expression_filter.tsv"))
  mde_write_tsv(data.frame(donor_id=donors,raw_library_size=fn$raw_library_size,filtered_library_size=fn$filtered_library_size,
                           norm_factor=fn$norm_factor,effective_library_size=fn$filtered_library_size*fn$norm_factor),file.path(d,"normalization.tsv"))
  mde_write_tsv(data.frame(donor_id=donors,built$design,check.names=FALSE),file.path(d,"design.tsv"))
  list(dir=d,n=N,family=fit$family_size,n_sig=sum(fit$table$adj.P.Val<0.05,na.rm=TRUE),design_columns=colnames(built$design),
       sub_col=sub_col,sub_levels=sub_levels,depth_col=depth_col)
}
A <- make_cohort("discovery",12,10,7,6,paste0("batch",11:16),"capture_batch","sub_phase",c("secretory_early","secretory_mid","secretory_late"),"proliferative","primary_lineage_cells",1.0,0.4)
B <- make_cohort("heca",15,7,16,4,c("20","42058","42196","43103","43281","43937"),"library_group_Lane","stage_fine",c("Secretory Early","Secretory Mid","Secretory"),"Proliferative","primary_nuclei",0.8,0.0)
for(pair in list(list("discovery",A,B),list("heca",B,A))) {
  key <- pair[[1]]; own <- pair[[2]]; other <- pair[[3]]
  mde_write_json(list(synthetic=TRUE,cohort=key,counts=file.path(own$dir,"counts.tsv.gz"),clinical=file.path(own$dir,"clinical.tsv"),
                      phase_fit_table=file.path(own$dir,"de_all.tsv"),expression_filter=file.path(own$dir,"expression_filter.tsv"),
                      normalization=file.path(own$dir,"normalization.tsv"),design=file.path(own$dir,"design.tsv"),
                      signature_table=file.path(other$dir,"de_all.tsv"),sub_phase_column=own$sub_col,sub_phase_levels=own$sub_levels,depth_column=own$depth_col),
                 file.path(out,paste0("inputs_",key,".json")))
}
mde_write_json(list(synthetic=TRUE,genes=G,planted_program_genes=150,planted_attenuation=list(discovery=0.4,heca=0.0),scale=list(discovery=1.0,heca=0.8),
                    cohorts=list(discovery=A[c("n","family","n_sig","design_columns")],heca=B[c("n","family","n_sig","design_columns")]),
                    expectation="theta_discovery about -0.4 x (scale ratio) with an interval excluding zero; theta_heca near zero",
                    note="invented data for code testing only; never a scientific result"),file.path(out,"fixture.json"))
cat("SYNTHETIC_FIXTURE_WRITTEN discovery family",A$family,"sig",A$n_sig,"| heca family",B$family,"sig",B$n_sig,"\n")
