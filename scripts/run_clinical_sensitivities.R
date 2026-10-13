#!/usr/bin/env Rscript
# Additive post-primary implementation of previously adopted clinical subsets.
source("scripts/donor_de_lib.R")
clinical_sensitivity_scenarios <- function() list(
  list(cohort="discovery",scenario="strict_histology_conflict_exclusion",expected_estimable=TRUE),
  list(cohort="discovery",scenario="adenomyosis_free",expected_estimable=TRUE),
  list(cohort="heca",scenario="adenomyosis_free",expected_estimable=TRUE),
  list(cohort="heca",scenario="strict_histology_conflict_exclusion",expected_estimable=FALSE),
  list(cohort="discovery",scenario="strict_and_adenomyosis_free",expected_estimable=FALSE))
sens_read <- function(path)read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,quote="",comment.char="",na.strings="NA")
sens_bool <- function(x){v<-tolower(as.character(x));if(anyNA(v)||any(!v %in% c("true","false")))stop("missing/invalid eligibility flag");v=="true"}
sens_sha <- function(path){x<-system2("sha256sum",shQuote(path),stdout=TRUE);if(length(x)!=1L||!is.null(attr(x,"status"))||!grepl("^[0-9a-f]{64} ",x))stop("hash failure");substr(x,1,64)}
sens_json <- function(x,path){tmp<-paste0(path,".partial.",Sys.getpid());jsonlite::write_json(x,tmp,pretty=TRUE,auto_unbox=TRUE,null="null",na="null",digits=16);if(!file.rename(tmp,path))stop("atomic JSON failed")}
sens_tsv <- function(x,path)write.table(x,path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
sens_verify_freeze <- function(path,phase){fr<-jsonlite::fromJSON(path,simplifyVector=FALSE);if(!identical(fr$phase,phase)||!length(fr$artifacts))stop("wrong freeze");for(p in names(fr$artifacts))if(!identical(sens_sha(p),fr$artifacts[[p]]))stop(paste("freeze changed",p));fr}
sensitivity_membership <- function(clinical,cohort,scenario){
  if(cohort=="discovery"){
    strict<-sens_bool(clinical$strict_histology_sensitivity)
    free<-sens_bool(clinical$adenomyosis_free_sensitivity)
    if(any(free!=(clinical$adenomyosis==0)))stop("adenomyosis-free flag conflicts with clinical covariate")
  } else if(cohort=="heca"){
    strict<-!sens_bool(clinical$histology_conflict);free<-clinical$adenomyosis==0
  } else stop("unknown cohort")
  switch(scenario,strict_histology_conflict_exclusion=strict,adenomyosis_free=free,strict_and_adenomyosis_free=strict&free,stop("unknown scenario"))
}
select_sensitivity_inputs <- function(counts,clinical,design,cohort,scenario,cfg){
  if(!identical(colnames(counts),as.character(clinical$donor_id)) || anyDuplicated(clinical$donor_id))stop("original exact count/clinical donor order required")
  if(!is.matrix(counts)||!is.numeric(counts)||any(!is.finite(counts))||any(counts<0)||any(colSums(counts)<=0)||anyDuplicated(rownames(counts))||any(!grepl("^ENSG[0-9]{11}$",rownames(counts))))stop("invalid source counts")
  if(cfg$cohorts[[cohort]]$integer_counts&&any(counts!=floor(counts)))stop("integer count contract failed")
  expected<-clinical$donor_id[sensitivity_membership(clinical,cohort,scenario)]
  if(is.null(rownames(design))||anyDuplicated(rownames(design))||!setequal(expected,rownames(design)))stop("design donors differ from exact frozen eligibility flags")
  ix<-match(rownames(design),clinical$donor_id)
  list(counts=counts[,ix,drop=FALSE],clinical=clinical[ix,,drop=FALSE],design=design,cohort=cohort,scenario=scenario)
}
sensitivity_eligibility <- function(x,cfg){
  cl<-x$clinical;design<-x$design;cohort<-x$cohort
  groups<-table(factor(cl$disease,levels=c("ENDO","non-ENDO")))
  rank<-qr(design,tol=cfg$fit$rank_tolerance)$rank;df<-nrow(design)-rank
  reason<-character()
  if(any(groups<cfg$cohorts[[cohort]]$minimum_per_group))reason<-c(reason,"disease group donor minimum failed")
  cells<-cl[[cfg$cohorts[[cohort]]$cell_column]]
  if(any(!is.finite(cells))||any(cells<cfg$design$min_cells))reason<-c(reason,"primary cell/nucleus minimum failed")
  if(any(!is.finite(design))||rank!=ncol(design))reason<-c(reason,"full rank design failed")
  if(df<cfg$minimum_residual_df)reason<-c(reason,"residual df minimum failed")
  if(!cfg$disease_coefficient %in% colnames(design)||!identical(as.numeric(design[,cfg$disease_coefficient]),as.numeric(cl$disease=="ENDO")))reason<-c(reason,"disease contrast coding mismatch")
  batch<-as.character(cl[[cfg$cohorts[[cohort]]$batch_column]])
  cross<-table(batch,factor(cl$disease,levels=c("non-ENDO","ENDO")))
  if(!any(rowSums(cross>0)==2))reason<-c(reason,"mixed disease batch support failed")
  list(estimable=!length(reason),reason=paste(reason,collapse="; "),n_donors=nrow(cl),ENDO=unname(groups[1]),non_ENDO=unname(groups[2]),rank=rank,residual_df=df,
       minimum_per_group=cfg$cohorts[[cohort]]$minimum_per_group,minimum_cells=min(cells),mixed_batch_levels=rownames(cross)[rowSums(cross>0)==2],disease_batch_table=unclass(cross))
}
clinical_sensitivity_main <- function(args){
  if(length(args)%%2)stop("use --key value pairs")
  opt<-list();if(length(args))for(i in seq(1,length(args),by=2)){key<-substring(args[i],3);if(!startsWith(args[i],"--")||!is.null(opt[[key]]))stop("bad option");opt[[key]]<-args[i+1]}
  defaults<-list(`counts-root`="data/processed",`clinical-root`="manifests",`design-root`="reports",config="config/statistical_analysis.json",`output-dir`="work/statistical/clinical_sensitivities")
  for(key in names(defaults))if(is.null(opt[[key]]))opt[[key]]<-defaults[[key]]
  synthetic<-identical(opt[["synthetic-test"]],"true")
  if(synthetic){for(key in c("counts-root","clinical-root","design-root","output-dir"))if(!startsWith(opt[[key]],"work/synthetic/"))stop("synthetic fixtures must stay in work/synthetic")} else {
    for(key in names(defaults))if(!identical(opt[[key]],defaults[[key]]))stop("production input/output paths must be fixed")
    pre<-sens_verify_freeze("receipts/statistical_pre_effect_freeze.json","STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT")
    sens_verify_freeze("receipts/discovery_rank_freeze.json","DISCOVERY_FROZEN_PRE_EXTERNAL")
  }
  if(Sys.getenv("OPENBLAS_NUM_THREADS")!="1"||Sys.getenv("OMP_NUM_THREADS")!="1")stop("source frozen statistical runtime first: BLAS/OMP must be1")
  cfg<-jsonlite::fromJSON(opt$config)
  scenarios<-clinical_sensitivity_scenarios();sources<-list();paths<-c(opt$config,"scripts/run_clinical_sensitivities.R","scripts/donor_de_lib.R","scripts/statistical_runtime_env.sh","requirements/statistical-runtime.explicit.txt")
  for(cohort in c("discovery","heca")){
    cp<-file.path(opt[["counts-root"]],cohort,"counts_stable.tsv.gz");lp<-file.path(opt[["clinical-root"]],paste0("final_",cohort,"_clinical_selected.tsv"))
    tab<-sens_read(cp);if(names(tab)[1]!="gene_id")stop("count schema")
    ct<-as.matrix(tab[,-1,drop=FALSE]);storage.mode(ct)<-"double";rownames(ct)<-tab$gene_id
    sources[[cohort]]<-list(counts=ct,clinical=sens_read(lp));paths<-c(paths,cp,lp)
  }
  prepared<-list()
  for(spec in scenarios){
    cohort<-spec$cohort;scenario<-spec$scenario;id<-paste0(cohort,"_",scenario)
    dp<-file.path(opt[["design-root"]],paste0("final_clinical_design_",id,".tsv"));tab<-sens_read(dp)
    if(names(tab)[1]!="donor_id")stop("design schema")
    design<-as.matrix(tab[,-1,drop=FALSE]);storage.mode(design)<-"double";rownames(design)<-tab$donor_id
    x<-select_sensitivity_inputs(sources[[cohort]]$counts,sources[[cohort]]$clinical,design,cohort,scenario,cfg)
    x$audit<-sensitivity_eligibility(x,cfg)
    if(!identical(x$audit$estimable,spec$expected_estimable))stop(paste("pre-existing viability ruling differs",id))
    if("age_centered_years" %in% colnames(design)&&any(abs(design[,"age_centered_years"]-(x$clinical$age-mean(x$clinical$age)))>1e-12))stop("frozen design age centering mismatch")
    if("adenomyosis" %in% colnames(design)&&!identical(as.numeric(design[,"adenomyosis"]),as.numeric(x$clinical$adenomyosis)))stop("frozen adenomyosis design mismatch")
    x$design_path<-dp;prepared[[id]]<-x;paths<-c(paths,dp)
  }
  if(!synthetic)paths<-c(paths,"receipts/statistical_pre_effect_freeze.json","receipts/discovery_rank_freeze.json","docs/PRE_EFFECT_METHOD_DECISIONS.md","docs/STATISTICAL_IMPLEMENTATION_DECISIONS.md","docs/REMAINING_SCIENTIFIC_WORK.md")
  hashes<-setNames(lapply(paths,sens_sha),paths)
  out<-opt[["output-dir"]]
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("existing sensitivity outputs preserved; do not repeat successful fits")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  design_history<-if(synthetic)NULL else setNames(lapply(vapply(prepared,function(x)x$design_path,""),function(p)system2("git",c("log","-1","--format=%H","--",shQuote(p)),stdout=TRUE)),names(prepared))
  snapshot<-list(status="OLD_CLINICAL_DESIGNS_AND_ELIGIBILITY_SNAPSHOTTED_BEFORE_SENSITIVITY_FITS",snapshot_utc=format(Sys.time(),tz="UTC",usetz=TRUE),input_sha256=hashes,
       design_last_commit=design_history,scenario_donors=lapply(prepared,function(x)x$clinical$donor_id),scenario_audits=lapply(prepared,function(x)x$audit),
       primary_external_effects_previously_observed=!synthetic,implementation_timing="new wrapper after primary effects; clinical subset definitions and designs pre-exist; not preregistered",synthetic=synthetic)
  sens_json(snapshot,file.path(out,"input_snapshot.json"))
  process<-list(pid=Sys.getpid(),host=Sys.info()[["nodename"]],started_utc=format(Sys.time(),tz="UTC",usetz=TRUE),workers=1,BLAS_threads=1,status="RUNNING",exit_status=NULL)
  sens_json(process,file.path(out,"process.json"));completed<-FALSE
  on.exit({process$status<-if(completed)"COMPLETED" else "FAILED_OR_INTERRUPTED";process$exit_status<-if(completed)0 else 1;process$finished_utc<-format(Sys.time(),tz="UTC",usetz=TRUE);sens_json(process,file.path(out,"process.json"))},add=TRUE)
  statuses<-list()
  for(id in names(prepared)){
    x<-prepared[[id]];target<-file.path(out,id);dir.create(target)
    sens_tsv(x$clinical,file.path(target,"clinical_selected.tsv"));sens_tsv(data.frame(donor_id=rownames(x$design),x$design,check.names=FALSE),file.path(target,"design.tsv"))
    if(!x$audit$estimable){status<-list(status="NOT_ESTIMABLE",fit_executed=FALSE,cohort=x$cohort,scenario=x$scenario,audit=x$audit);sens_json(status,file.path(target,"status.json"));statuses[[id]]<-status;cat("SENSITIVITY_NOT_ESTIMABLE",id,x$audit$reason,"\n");next}
    started<-proc.time()[["elapsed"]]
    model<-fit_donor_model(x$counts,x$design,x$clinical,x$cohort,cfg)
    sens_tsv(model$table,file.path(target,"de_all.tsv"));sens_tsv(model$table[model$table$fit_estimable,],file.path(target,"de_fitted.tsv"))
    sens_tsv(data.frame(gene_id=rownames(x$counts),expression_pass=model$keep),file.path(target,"expression_filter.tsv"));sens_tsv(model$normalization,file.path(target,"normalization.tsv"))
    saveRDS(model,file.path(target,"model.rds"),version=3,compress="gzip");writeLines(capture.output(sessionInfo()),file.path(target,"sessionInfo.txt"))
    status<-list(status="COMPLETE_SENSITIVITY_ASSOCIATION",fit_executed=TRUE,cohort=x$cohort,scenario=x$scenario,audit=model$audit,clinical_audit=x$audit,elapsed_seconds=proc.time()[["elapsed"]]-started,
       expression_pass_family=nrow(model$table),estimable_genes=sum(model$table$fit_estimable),BH05=sum(model$table$adj.P.Val<.05,na.rm=TRUE),minimum_BH_q=min(model$table$adj.P.Val,na.rm=TRUE),captured_warnings=model$warnings,effective_settings=model$settings)
    sens_json(status,file.path(target,"status.json"));statuses[[id]]<-status
    sens_json(list(scientific_effect_run=!synthetic,synthetic=synthetic,cohort=x$cohort,scenario=x$scenario,input_snapshot_sha256=sens_sha(file.path(out,"input_snapshot.json")),input_sha256=hashes,
         output_sha256=setNames(lapply(list.files(target,full.names=TRUE),sens_sha),basename(list.files(target,full.names=TRUE))),implementation_timing=snapshot$implementation_timing,
         source_primary_ledgers_preserved=TRUE,scientific_scope="observational sensitivity association, not causal or therapeutic validation"),file.path(target,"provenance.json"))
    cat("SENSITIVITY_COMPLETE",id,"genes",nrow(model$table),"BH05",status$BH05,"\n")
  }
  for(p in paths)if(!identical(sens_sha(p),hashes[[p]]))stop(paste("sensitivity input/code changed",p))
  sens_json(list(status="COMPLETED",scenarios=statuses,models_fit=sum(vapply(statuses,function(x)x$fit_executed,TRUE)),scientific_validation=FALSE,
       original_lists_replaced=FALSE,synthetic=synthetic,input_snapshot_sha256=sens_sha(file.path(out,"input_snapshot.json"))),file.path(out,"summary.json"))
  completed<-TRUE
}
if(sys.nframe()==0L)clinical_sensitivity_main(commandArgs(trailingOnly=TRUE))
