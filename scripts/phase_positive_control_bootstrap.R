#!/usr/bin/env Rscript
# Phase positive control, endpoint S2: paired external-donor bootstrap of the G1-minus-B0
# top-20 concordance, with the v1 resampling semantics. The v1 module is sourced unchanged
# and its fixed-family fit, mask validation, paired evaluation, checkpoint and summary
# functions are reused. Two things are phase-specific and defined here: draws are stratified
# by cycle phase, and the design is rebuilt with SEC_minus_PRO as the contrast and disease as
# a nuisance covariate. The reused fit extracts the coefficient named in cfg$disease_coefficient,
# so a copy of the configuration with that name set to SEC_minus_PRO is passed to it.
source("scripts/resample_donor_models.R")  # guarded entry point; sourcing runs nothing
rebuild_phase_design <- function(clinical,cfg,pcfg) {
  levels_phase <- c(pcfg$phase_factor$reference,pcfg$phase_factor$case)
  batch_col <- cfg$cohorts$heca$batch_column
  required <- c("donor_id","phase","disease",cfg$design$age_column,cfg$design$adenomyosis_column,batch_col)
  if(!all(required %in% names(clinical))) stop("missing clinical design field")
  if(anyDuplicated(clinical$donor_id) || anyNA(clinical[,required])) stop("duplicate draw IDs or missing covariate")
  if(any(!clinical$phase %in% levels_phase) || length(unique(clinical$phase))!=2) stop("both phase groups required")
  age <- as.numeric(clinical[[cfg$design$age_column]]); adeno <- as.numeric(clinical[[cfg$design$adenomyosis_column]])
  if(any(!is.finite(age)) || any(!adeno %in% c(0,1))) stop("invalid age/adenomyosis")
  batch <- as.character(clinical[[batch_col]]); if(any(!nzchar(batch))) stop("missing batch")
  levels <- sort(unique(batch))
  cross <- table(factor(batch,levels=levels),factor(clinical$phase,levels=levels_phase))
  if(!any(rowSums(cross>0)==2)) stop("no mixed phase batch support")
  design <- cbind(Intercept=rep(1,nrow(clinical)),SEC_minus_PRO=as.numeric(clinical$phase==pcfg$phase_factor$case))
  dropped <- character()
  disease <- as.numeric(clinical$disease=="ENDO")
  if(length(unique(disease))==1) dropped<-c(dropped,"disease") else design<-cbind(design,ENDO_minus_non_ENDO=disease)
  if(length(unique(age))==1) dropped<-c(dropped,"age") else design<-cbind(design,age_centered=age-mean(age))
  if(length(levels)==1) dropped<-c(dropped,batch_col) else {
    contrasts <- model.matrix(~factor(batch,levels=levels))
    design<-cbind(design,contrasts[,-1,drop=FALSE])
  }
  if(length(unique(adeno))==1) dropped<-c(dropped,"adenomyosis") else design<-cbind(design,adenomyosis=adeno)
  rownames(design)<-clinical$donor_id
  rank<-qr(design,tol=cfg$fit$rank_tolerance)$rank
  if(rank!=ncol(design)) stop("design not full rank; nonconstant nuisance cannot be dropped")
  df <- nrow(design)-rank
  if(df<cfg$minimum_residual_df) stop("nominal residual df below 5")
  list(design=design,receipt=list(nominal_rows=nrow(design),rank=rank,nominal_residual_df=df,age_center=mean(age),
       observed_sorted_batch_levels=levels,dropped_constant_nuisance=dropped,phase_batch_table=unclass(cross),
       unique_original_donors=length(unique(clinical$donor_id)),duplicated_rows_are_not_new_people=TRUE))
}
generate_phase_draws <- function(clinical,cfg,pcfg) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(cfg$bootstrap$seed)
  groups<-lapply(c(pcfg$phase_factor$case,pcfg$phase_factor$reference),function(g)which(clinical$phase==g))
  do.call(rbind,lapply(seq_len(cfg$bootstrap$n),function(r){
    ix<-unlist(lapply(groups,function(v)sample(v,length(v),replace=TRUE)),use.names=FALSE)
    data.frame(replicate=r,draw_position=seq_along(ix),source_index=ix,source_donor_id=clinical$donor_id[ix],phase=clinical$phase[ix],stringsAsFactors=FALSE)
  }))
}
verify_phase_freeze <- function(path,phase) {
  fr <- jsonlite::fromJSON(path,simplifyVector=FALSE)
  if(!identical(fr$phase,phase) || !length(fr$artifacts)) stop("freeze schema/phase mismatch")
  for(p in names(fr$artifacts)) {
    ok <- file.exists(p) && identical(resample_sha(p),fr$artifacts[[p]])
    if(!ok && !startsWith(p,"docs/")) stop(paste("freeze changed",p))
  }
  fr
}
phase_bootstrap_main <- function(args) {
  if(length(args)%%2) stop("use --key value pairs")
  opt<-list(); for(i in seq(1,length(args),by=2)){k<-substring(args[i],3);if(!startsWith(args[i],"--")||!is.null(opt[[k]]))stop("bad option");opt[[k]]<-args[i+1]}
  required<-c("counts","clinical","config","phase-config","primary-model","eval-mask","f1-freeze","output-dir")
  if(!all(required %in% names(opt))) stop(paste("required:",paste(required,collapse=" ")))
  f1<-verify_phase_freeze(opt[["f1-freeze"]],"PHASE_POSITIVE_CONTROL_F1_DISCOVERY_FROZEN_PRE_HECA")
  cfg<-jsonlite::fromJSON(opt$config);pcfg<-jsonlite::fromJSON(opt[["phase-config"]])
  if(!identical(as.integer(cfg$bootstrap$n),500L) || cfg$bootstrap$seed!=42 || !isTRUE(cfg$bootstrap$fixed_external_gene_family) || !isFALSE(cfg$bootstrap$refilter) || !isFALSE(cfg$bootstrap$redraw_invalid))stop("bootstrap contract differs")
  cfg_phase<-cfg; cfg_phase$disease_coefficient<-pcfg$coefficient  # the reused fit reads this name only to select the coefficient column
  clinical<-read_resample_tsv(opt$clinical); ct<-read_resample_tsv(opt$counts)
  if(names(ct)[1]!="gene_id")stop("counts schema")
  counts<-as.matrix(ct[,-1,drop=FALSE]);storage.mode(counts)<-"double";rownames(counts)<-ct$gene_id
  if(!identical(colnames(counts),as.character(clinical$donor_id)) || anyDuplicated(clinical$donor_id) || anyDuplicated(rownames(counts)) || any(!grepl("^ENSG[0-9]{11}$",rownames(counts))) || any(!is.finite(counts)) || any(counts<0) || any(colSums(counts)<=0))stop("source input contract failed")
  if(isTRUE(pcfg$cohorts$heca$integer_counts) && any(counts!=floor(counts)))stop("source count integer contract failed")
  cells<-clinical[[cfg$cohorts$heca$cell_column]]
  if(any(!is.finite(cells)) || any(cells<cfg$design$min_cells))stop("original donor cells failed")
  phases<-table(factor(clinical$phase,levels=c(pcfg$phase_factor$reference,pcfg$phase_factor$case)))
  if(any(phases<pcfg$cohorts$heca$minimum_per_phase_group))stop("phase group minimum failed")
  paths<-c(opt$counts,opt$clinical,opt$config,opt[["phase-config"]],opt[["f1-freeze"]],"scripts/phase_positive_control_bootstrap.R","scripts/resample_donor_models.R","scripts/donor_de_lib.R","scripts/statistical_runtime_env.sh","requirements/statistical-runtime.explicit.txt")
  out<-opt[["output-dir"]];dir.create(out,recursive=TRUE,showWarnings=FALSE)
  primary<-readRDS(opt[["primary-model"]]);family<-primary$table$gene_id
  if(!identical(family,rownames(counts)[primary$keep]) || !identical(primary$normalization$donor_id,colnames(counts)))stop("primary original family/order does not match source")
  mask<-read_resample_tsv(opt[["eval-mask"]]);validate_mask(mask)
  if(any(!mask$gene_id[mask$original_evaluable] %in% family))stop("originally evaluable mask outside original family")
  original_ix<-match(mask$gene_id,primary$table$gene_id)
  expected_eval<-!is.na(original_ix) & mask$discovery_sign!=0
  present<-which(!is.na(original_ix)); t0<-primary$table[original_ix[present],]
  expected_eval[present]<-expected_eval[present] & t0$fit_estimable & is.finite(t0$logFC) & is.finite(t0$t) & is.finite(t0$P.Value) & t0$logFC!=0
  if(!identical(as.logical(mask$original_evaluable),as.logical(expected_eval)))stop("mask differs from original external evaluability")
  for(m in c("B0","B1","G1"))for(k in c(20,50)){
    frozen<-read_resample_tsv(f1$products[[paste0(m,"_top",k)]])
    x<-mask[mask$method==m & mask$k==k,]
    if(!identical(x$gene_id,frozen$gene_id) || !identical(as.numeric(x$discovery_sign),as.numeric(frozen$discovery_sign)))stop("mask lists/signs differ from the F1 freeze")
  }
  paths<-c(paths,opt[["primary-model"]],opt[["eval-mask"]])
  primary_prov<-file.path(dirname(opt[["primary-model"]]),"provenance.json")
  pp<-jsonlite::fromJSON(primary_prov,simplifyVector=FALSE)
  if(!isTRUE(pp$scientific_effect_run) || pp$cohort!="heca" || pp$contrast!="SEC_minus_PRO")stop("HECA phase primary provenance absent")
  for(p in c(opt$counts,opt$clinical,opt$config))if(!identical(pp$input_sha256[[p]],resample_sha(p)))stop("primary input hash mismatch")
  for(p in names(pp$output_sha256))if(!identical(pp$output_sha256[[p]],resample_sha(p)))stop("HECA primary output hash mismatch")
  paths<-c(paths,primary_prov)
  binding<-list(mode="phase_bootstrap",inputs=setNames(lapply(paths,resample_sha),paths),family=family,mask_rows=nrow(mask),R_version=as.character(getRversion()),edgeR=as.character(packageVersion("edgeR")),limma=as.character(packageVersion("limma")),statmod=as.character(packageVersion("statmod")),stratification="phase",coefficient=pcfg$coefficient)
  manifest<-file.path(out,"binding.rds");ledger_path<-file.path(out,"draw_ledger.tsv")
  if(file.exists(manifest)){
    saved<-readRDS(manifest)
    if(!identical(binding,saved$binding) || !file.exists(ledger_path) || resample_sha(ledger_path)!=saved$ledger_sha256)stop("unsafe resume binding/ledger changed")
    ledger<-read_resample_tsv(ledger_path)
    if(!identical(ledger,generate_phase_draws(clinical,cfg,pcfg)))stop("draw ledger not exact frozen RNG draws")
  } else {
    if(length(list.files(out,all.files=TRUE,no..=TRUE)))stop("nonempty unbound bootstrap output")
    ledger<-generate_phase_draws(clinical,cfg,pcfg)
    temp<-paste0(ledger_path,".partial");write_resample_tsv(ledger,temp);if(!file.rename(temp,ledger_path))stop("ledger rename failed")
    saveRDS(list(binding=binding,ledger_sha256=resample_sha(ledger_path)),paste0(manifest,".partial"),version=3);if(!file.rename(paste0(manifest,".partial"),manifest))stop("binding rename failed")
    write_atomic_json(list(binding=binding,ledger_sha256=resample_sha(ledger_path),ledger_saved_before_any_fit=TRUE,RNG="Mersenne-Twister/Inversion/Rejection",seed=42,stratification=as.list(phases)),file.path(out,"draw_freeze.json"))
  }
  fingerprint<-resample_sha(manifest)
  process<-list(pid=Sys.getpid(),host="SRV01",started_utc=format(Sys.time(),tz="UTC",usetz=TRUE),binding_sha256=fingerprint,status="RUNNING",exit_status=NULL)
  write_atomic_json(process,file.path(out,"process.json"))
  completed<-FALSE
  on.exit({process$status<-if(completed)"COMPLETED" else "INTERRUPTED_OR_FAILED";process$exit_status<-if(completed)0 else 1;process$finished_utc<-format(Sys.time(),tz="UTC",usetz=TRUE);write_atomic_json(process,file.path(out,"process.json"))},add=TRUE)
  counts<-counts[match(family,rownames(counts)),,drop=FALSE]
  differences<-rep(NA_real_,cfg$bootstrap$n);reasons<-character(cfg$bootstrap$n)
  workers<-if(is.null(opt$workers))1L else as.integer(opt$workers)
  if(is.na(workers) || workers<1 || workers>4)stop("--workers must be 1..4")
  process$workers<-workers;process$BLAS_threads<-1;write_atomic_json(process,file.path(out,"process.json"))
  Sys.setenv(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
  fit_one<-function(r) {
    path<-file.path(out,sprintf("replicate_%03d.json",r))
    if(file.exists(path)){
      result<-read_bootstrap_checkpoint(path,fingerprint,r,mask,length(family))
    } else {
      rows<-ledger[ledger$replicate==r,];ix<-rows$source_index
      clin<-clinical[ix,,drop=FALSE];clin$donor_id<-paste0("draw",seq_along(ix));cts<-counts[,ix,drop=FALSE];colnames(cts)<-clin$donor_id
      start<-proc.time()[["elapsed"]];design_receipt<-NULL
      result<-tryCatch({built<-rebuild_phase_design(clin,cfg,pcfg);design_receipt<-built$receipt;fit<-fit_fixed_family(cts,built$design,cfg_phase);ans<-evaluate_fixed_masks(fit$table,mask);ans$warnings<-fit$warnings;ans$family_size<-fit$family_size;ans},error=function(e)list(valid=FALSE,reason=conditionMessage(e)))
      result$worker_pid<-Sys.getpid();result$paired_replicate_valid<-isTRUE(result$valid)
      result$replicate<-r;result$binding_sha256<-fingerprint;result$elapsed_seconds<-proc.time()[["elapsed"]]-start
      result$unique_actual_people<-length(unique(rows$source_donor_id));result$unique_actual_secretory<-length(unique(rows$source_donor_id[rows$phase==pcfg$phase_factor$case]));result$unique_actual_proliferative<-length(unique(rows$source_donor_id[rows$phase==pcfg$phase_factor$reference]))
      if(!is.null(design_receipt))design_receipt$unique_original_donors<-result$unique_actual_people
      result$design_receipt<-design_receipt
      write_bootstrap_checkpoint(result,path)
    }
    cat("BOOTSTRAP_CHECKPOINT",r,"valid",isTRUE(result$valid),"pid",Sys.getpid(),"\n")
    result
  }
  results<-if(workers==1L)lapply(seq_len(cfg$bootstrap$n),fit_one) else parallel::mclapply(seq_len(cfg$bootstrap$n),fit_one,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  for(r in seq_along(results)){result<-results[[r]];if(inherits(result,"try-error") || is.null(result$valid))stop("worker/checkpoint failed; preserve completed draw files");if(isTRUE(result$valid))differences[r]<-result$difference else reasons[r]<-result$reason}
  for(p in paths)if(!identical(resample_sha(p),binding$inputs[[p]]))stop(paste("input/code changed during bootstrap",p))
  if(!identical(resample_sha(manifest),fingerprint))stop("binding changed during bootstrap")
  summary<-bootstrap_summary(differences,cfg$bootstrap$n);summary$binding_sha256<-fingerprint;summary$failure_reasons<-as.list(table(reasons[nzchar(reasons)]));summary$stratification<-"phase";summary$contrast<-pcfg$coefficient
  summary$point_estimate_fixed_sample<-NULL
  write_resample_tsv(data.frame(replicate=seq_along(differences),valid=is.finite(differences),G1_minus_B0_top20=differences,reason=reasons),file.path(out,"paired_differences.tsv"))
  write_atomic_json(summary,file.path(out,"summary.json"));completed<-TRUE
  cat("BOOTSTRAP_COMPLETE",summary$valid,"valid",summary$invalid,"invalid\n")
}
if(sys.nframe()==0L)phase_bootstrap_main(commandArgs(trailingOnly=TRUE))
