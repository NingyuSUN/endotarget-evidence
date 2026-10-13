#!/usr/bin/env Rscript
# Run from the project root; every real run requires a matching pre-effect freeze.
source("scripts/donor_de_lib.R")
args <- commandArgs(trailingOnly=TRUE)
if (length(args)%%2!=0) stop("arguments must be --key value pairs")
opts <- list()
if(length(args)) for(i in seq(1,length(args),by=2)) {
  if(!startsWith(args[i],"--") || !is.null(opts[[substring(args[i],3)]])) stop("bad/duplicate argument")
  opts[[substring(args[i],3)]] <- args[i+1]
}
required <- c("counts","design","clinical","cohort","config","output-dir","pre-effect-freeze")
if(!all(required %in% names(opts))) stop(paste("required:",paste(required,collapse=" ")))
sha <- function(path) {
  ans <- system2("sha256sum",shQuote(path),stdout=TRUE,stderr=TRUE)
  if(!is.null(attr(ans,"status")) || length(ans)!=1L ||
     !grepl("^[0-9a-f]{64} ",ans)) stop(paste("cannot hash",path))
  substr(ans,1,64)
}
verify_freeze <- function(path,expected_phase) {
  fr <- jsonlite::fromJSON(path,simplifyVector=FALSE)
  if(!identical(fr$phase,expected_phase) || length(fr$artifacts)==0L) stop("wrong freeze phase/schema")
  for(p in names(fr$artifacts))
    if(!identical(sha(p),fr$artifacts[[p]])) stop(paste("frozen artifact changed:",p))
  fr
}
synthetic <- identical(opts[["synthetic-test"]],"true")
expected_phase <- if(synthetic) "SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE" else "STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT"
pre <- verify_freeze(opts[["pre-effect-freeze"]],expected_phase)
for(p in c(opts$counts,opts$design,opts$clinical,opts$config,
           "scripts/fit_donor_de.R","scripts/donor_de_lib.R"))
  if(!p %in% names(pre$artifacts)) stop(paste("input/code omitted from pre-effect freeze",p))
if(opts$cohort=="heca" && !synthetic) {
  if(is.null(opts[["discovery-freeze"]])) stop("HECA effects require prior discovery rank freeze")
  dfreeze <- verify_freeze(opts[["discovery-freeze"]],"DISCOVERY_FROZEN_PRE_EXTERNAL")
  roles <- c("discovery_de_all","discovery_de_fitted","universe","induced_edges",
             "all_U_directions","null_bins","null_status",
             "null_memberships","statistical_config","runtime_lock","network",
             "discovery_counts","discovery_design","discovery_clinical",
             unlist(lapply(c("B0","B1","G1"),function(x)paste0(x,c("_full","_top20","_top50")))))
  if(!all(roles %in% names(dfreeze$products)) ||
     !all(unlist(dfreeze$products[roles]) %in% names(dfreeze$artifacts)))
    stop("discovery freeze lacks mandatory ranked/null/input products")
}
cfg <- jsonlite::fromJSON(opts$config,simplifyVector=TRUE)
if (!opts$cohort %in% names(cfg$cohorts)) stop("unknown cohort")
if (!synthetic && !identical(opts$design,cfg$cohorts[[opts$cohort]]$primary_design))
  stop("primary run requires the locked cohort-specific design path")
if (!synthetic && (!identical(opts$counts,paste0("data/processed/",opts$cohort,"/counts_stable.tsv.gz")) ||
    !identical(opts$clinical,paste0("manifests/final_",opts$cohort,"_clinical_selected.tsv"))))
  stop("primary run requires the locked cohort-specific count and clinical paths")
if (synthetic && (!all(startsWith(c(opts$counts,opts$design,opts$clinical),"work/synthetic/")) ||
    !all(startsWith(normalizePath(c(opts$counts,opts$design,opts$clinical),mustWork=TRUE),
                   paste0(normalizePath("work/synthetic",mustWork=TRUE),"/")))))
  stop("synthetic mode requires work/synthetic/ input fixtures")
readtab <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,
                                   na.strings=c("NA"),quote="",comment.char="")
ct <- readtab(opts$counts)
if(names(ct)[1]!="gene_id") stop("counts first column must be gene_id")
counts <- as.matrix(ct[,-1,drop=FALSE]); storage.mode(counts)<-"double"
rownames(counts) <- ct$gene_id
dt <- readtab(opts$design)
if(names(dt)[1]!="donor_id") stop("design first column must be donor_id")
design <- as.matrix(dt[,-1,drop=FALSE]); storage.mode(design)<-"double"
rownames(design)<-dt$donor_id
clinical <- readtab(opts$clinical)
out <- opts[["output-dir"]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("output directory is not empty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
started <- format(Sys.time(),tz="UTC",usetz=TRUE)
model <- fit_donor_model(counts,design,clinical,opts$cohort,cfg)
write_tsv <- function(x,name) write.table(x,file=file.path(out,name),sep="\t",
                       row.names=FALSE,col.names=TRUE,quote=FALSE,na="NA")
write_tsv(model$table,"de_all.tsv")
write_tsv(model$table[model$table$fit_estimable,,drop=FALSE],"de_fitted.tsv")
write_tsv(data.frame(gene_id=rownames(counts),expression_pass=model$keep),"expression_filter.tsv")
write_tsv(model$normalization,"normalization.tsv")
write_tsv(dt,"design.tsv")
saveRDS(model,file=file.path(out,"model.rds"),version=3,compress="gzip")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
input_paths <- c(opts$counts,opts$design,opts$clinical,opts$config,
                 opts[["pre-effect-freeze"]],"scripts/fit_donor_de.R","scripts/donor_de_lib.R")
if(!is.null(opts[["discovery-freeze"]]))input_paths<-c(input_paths,opts[["discovery-freeze"]])
outputs <- list.files(out,full.names=TRUE)
prov <- list(started_utc=started,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
             cohort=opts$cohort,data_kind=if(synthetic)"synthetic_test_only" else "public_cohort",
             scientific_effect_run=!synthetic,audit=model$audit,
             tested_family_size=nrow(model$table),estimable_genes=sum(model$table$fit_estimable),
             significant_BH05=sum(model$table$adj.P.Val<.05,na.rm=TRUE),
             effective_compatibility=model$settings,filter_audit=model$filter_audit,
             captured_warnings=model$warnings,
             input_sha256=setNames(lapply(input_paths,sha),input_paths),
             output_sha256=setNames(lapply(outputs,sha),outputs),
             scientific_scope=if(synthetic)"invented fixture, no biological evidence" else "observational donor-level association, not causal or therapeutic validation")
jsonlite::write_json(prov,file.path(out,"provenance.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
cat("DONOR_DE_COMPLETE",opts$cohort,"tested",nrow(model$table),
    "estimable",sum(model$table$fit_estimable),"\n")
