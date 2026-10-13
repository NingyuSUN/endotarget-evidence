#!/usr/bin/env Rscript
# New resampling module. Frozen primary library is sourced unchanged.
source("scripts/donor_de_lib.R")
read_resample_tsv <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,quote="",comment.char="",na.strings="NA")
resample_sha <- function(path) {
  ans <- system2("sha256sum",shQuote(path),stdout=TRUE)
  if (!is.null(attr(ans,"status")) || length(ans)!=1 || !grepl("^[0-9a-f]{64} ",ans)) stop(paste("hash failed",path))
  substr(ans,1,64)
}
write_atomic_json <- function(x,path) {
  temporary <- paste0(path,".partial.",Sys.getpid())
  jsonlite::write_json(x,temporary,pretty=TRUE,auto_unbox=TRUE,null="null",na="null",digits=16)
  if(!file.rename(temporary,path)) stop("atomic JSON rename failed")
}
write_resample_tsv <- function(x,path) write.table(x,path,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
verify_resample_freeze <- function(path,phase) {
  fr <- jsonlite::fromJSON(path,simplifyVector=FALSE)
  if(!identical(fr$phase,phase) || !length(fr$artifacts)) stop("freeze schema/phase mismatch")
  for(p in names(fr$artifacts)) if(!identical(resample_sha(p),fr$artifacts[[p]])) stop(paste("freeze changed",p))
  fr
}
rebuild_resample_design <- function(clinical,cohort,cfg) {
  required <- c("donor_id","disease",cfg$design$age_column,cfg$design$adenomyosis_column,cfg$cohorts[[cohort]]$batch_column)
  if(!all(required %in% names(clinical))) stop("missing clinical design field")
  if(anyDuplicated(clinical$donor_id) || anyNA(clinical[,required])) stop("duplicate draw IDs or missing covariate")
  if(any(!clinical$disease %in% c("ENDO","non-ENDO")) || length(unique(clinical$disease))!=2) stop("both disease groups required")
  age <- as.numeric(clinical[[cfg$design$age_column]])
  adeno <- as.numeric(clinical[[cfg$design$adenomyosis_column]])
  if(any(!is.finite(age)) || any(!adeno %in% c(0,1))) stop("invalid age/adenomyosis")
  batch <- as.character(clinical[[cfg$cohorts[[cohort]]$batch_column]])
  if(any(!nzchar(batch))) stop("missing batch")
  levels <- sort(unique(batch))
  cross <- table(factor(batch,levels=levels),factor(clinical$disease,levels=c("non-ENDO","ENDO")))
  if(!any(rowSums(cross>0)==2)) stop("no mixed disease batch support")
  design <- cbind(Intercept=rep(1,nrow(clinical)),ENDO_minus_non_ENDO=as.numeric(clinical$disease=="ENDO"))
  dropped <- character()
  if(length(unique(age))==1) dropped<-c(dropped,"age") else design<-cbind(design,age_centered=age-mean(age))
  if(length(levels)==1) dropped<-c(dropped,cfg$cohorts[[cohort]]$batch_column) else {
    contrasts <- model.matrix(~factor(batch,levels=levels))
    design<-cbind(design,contrasts[,-1,drop=FALSE])
  }
  if(length(unique(adeno))==1) dropped<-c(dropped,"adenomyosis") else design<-cbind(design,adenomyosis=adeno)
  rownames(design)<-clinical$donor_id
  rank<-qr(design,tol=cfg$fit$rank_tolerance)$rank
  if(rank!=ncol(design)) stop("design not full rank; nonconstant nuisance cannot be dropped")
  df <- nrow(design)-rank
  if(df<cfg$minimum_residual_df) stop("nominal residual df below 5")
  list(design=design,receipt=list(nominal_rows=nrow(design),rank=rank,nominal_residual_df=df,
       age_center=mean(age),observed_sorted_batch_levels=levels,dropped_constant_nuisance=dropped,
       disease_batch_table=unclass(cross),unique_original_donors=length(unique(clinical$donor_id)),
       duplicated_rows_are_not_new_people=TRUE))
}
fit_fixed_family <- function(counts,design,cfg) {
  if(!is.matrix(counts) || !is.numeric(counts) || nrow(counts)<2 || any(!is.finite(counts)) || any(counts<0) || any(colSums(counts)<=0)) stop("invalid fixed-family counts/library")
  if(!identical(colnames(counts),rownames(design)) || anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts))) stop("fixed family names/order invalid")
  if(any(!is.finite(design)) || qr(design,tol=cfg$fit$rank_tolerance)$rank!=ncol(design) || nrow(design)-ncol(design)<cfg$minimum_residual_df) stop("fixed-family design invalid")
  warnings <- character()
  result <- withCallingHandlers({
    y<-edgeR::DGEList(counts=counts)
    n<-cfg$normalization
    y<-edgeR::calcNormFactors(y,method=n$method,refColumn=NULL,logratioTrim=n$logratio_trim,sumTrim=n$sum_trim,doWeighting=n$do_weighting,Acutoff=n$A_cutoff,p=n$p)
    v<-cfg$voom
    va<-list(counts=y,design=design,lib.size=NULL,normalize.method=v$normalize_method,block=NULL,correlation=NULL,weights=NULL,span=v$span,plot=v$plot,save.plot=v$save_plot)
    adaptive<-"adaptive.span" %in% names(formals(limma::voom))
    if(adaptive) va$adaptive.span<-v$adaptive_span else if(as.character(packageVersion("limma"))!="3.62.1") stop("unreviewed voom compatibility")
    vo<-do.call(limma::voom,va); validate_voom_weights(vo$weights)
    fit<-limma::lmFit(vo,design=design,method=cfg$fit$method)
    e<-cfg$ebayes
    ea<-list(fit=fit,proportion=e$proportion,stdev.coef.lim=e$stdev_coef_lim,trend=e$trend,robust=e$robust,winsor.tail.p=e$winsor_tail_p)
    legacy<-"legacy" %in% names(formals(limma::eBayes))
    if(legacy) ea$legacy<-e$legacy else if(as.character(packageVersion("limma"))!="3.62.1") stop("unreviewed eBayes compatibility")
    fit<-do.call(limma::eBayes,ea)
    co<-match(cfg$disease_coefficient,colnames(design)); if(is.na(co)) stop("disease coefficient missing")
    logfc<-fit$coefficients[,co]; statistic<-fit$t[,co]; p<-fit$p.value[,co]
    se<-sqrt(fit$s2.post)*fit$stdev.unscaled[,co]; degrees<-rep_len(fit$df.total,length(logfc))
    finite<-is.finite(logfc)&is.finite(statistic)&is.finite(p)&is.finite(se)&se>0&is.finite(degrees)&degrees>0&is.finite(fit$df.residual)&fit$df.residual>=cfg$minimum_residual_df
    q<-rep(NA_real_,length(p)); q[finite]<-p.adjust(p[finite],method=cfg$testing$adjust_method,n=length(p))
    list(table=data.frame(gene_id=rownames(counts),logFC=logfc,t=statistic,P.Value=p,adj.P.Val=q,fit_estimable=finite),
         family_size=nrow(counts),normalization=data.frame(donor_id=colnames(counts),filtered_library_size=y$samples$lib.size,norm_factor=y$samples$norm.factors),
         effective_settings=list(adaptive_span_supported=adaptive,adaptive_span=FALSE,legacy_supported=legacy,legacy=TRUE))
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  result$warnings<-warnings; result
}
generate_draws <- function(clinical,cfg) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(cfg$bootstrap$seed)
  groups<-lapply(c("ENDO","non-ENDO"),function(g)which(clinical$disease==g))
  do.call(rbind,lapply(seq_len(cfg$bootstrap$n),function(r){
    ix<-unlist(lapply(groups,function(v)sample(v,length(v),replace=TRUE)),use.names=FALSE)
    data.frame(replicate=r,draw_position=seq_along(ix),source_index=ix,source_donor_id=clinical$donor_id[ix],disease=clinical$disease[ix],stringsAsFactors=FALSE)
  }))
}
validate_mask <- function(mask) {
  fields<-c("method","k","rank","gene_id","discovery_sign","original_evaluable")
  if(!all(fields %in% names(mask))) stop("bootstrap mask schema missing")
  if(!is.logical(mask$original_evaluable) || anyNA(mask$original_evaluable)) stop("mask must use TRUE/FALSE")
  if(any(!mask$discovery_sign %in% c(-1,0,1))) stop("invalid frozen discovery sign")
  for(m in c("B0","B1","G1")) for(k in c(20,50)) {
    x<-mask[mask$method==m & mask$k==k,,drop=FALSE]
    if(nrow(x)!=k || anyDuplicated(x$gene_id) || !identical(as.integer(x$rank),seq_len(k))) stop("fixed mask must contain exact complete top20/top50")
    if(any(x$original_evaluable & x$discovery_sign==0) || !any(x$original_evaluable)) stop("invalid/empty originally evaluable list")
  }
  if(nrow(mask)!=210) stop("unexpected mask rows")
  for(g in unique(mask$gene_id)) {
    x<-mask[mask$gene_id==g,]
    if(length(unique(x$original_evaluable))!=1 || length(unique(x$discovery_sign))!=1) stop("inconsistent gene evaluability/sign across lists")
  }
  invisible(TRUE)
}
evaluate_fixed_masks <- function(tab,mask) {
  union<-unique(mask$gene_id[mask$original_evaluable]); ix<-match(union,tab$gene_id)
  bad<-is.na(ix)
  present<-which(!bad)
  bad[present]<-!tab$fit_estimable[ix[present]] | !is.finite(tab$logFC[ix[present]]) | !is.finite(tab$t[ix[present]]) | !is.finite(tab$P.Value[ix[present]]) | tab$logFC[ix[present]]==0
  if(any(bad)) return(list(valid=FALSE,reason="invalid_or_zero_originally_evaluable_union_gene",failed_genes=union[bad]))
  rates<-list()
  for(m in c("B0","B1","G1")) for(k in c(20,50)) {
    x<-mask[mask$method==m & mask$k==k & mask$original_evaluable,]
    agrees<-sum(sign(tab$logFC[match(x$gene_id,tab$gene_id)])==x$discovery_sign)
    rates[[paste0(m,"_top",k)]]<-list(agree=agrees,evaluable=nrow(x),selected=k,rate=agrees/nrow(x))
  }
  list(valid=TRUE,reason="",rates=rates,difference=rates$G1_top20$rate-rates$B0_top20$rate)
}
bootstrap_summary <- function(differences,n=500) {
  invalid<-sum(!is.finite(differences));complete<-length(differences)==n
  suppress<-!complete || invalid>n*.1 || !any(is.finite(differences))
  interval<-if(suppress)NULL else as.numeric(quantile(differences[is.finite(differences)],c(.025,.975),type=7,names=FALSE))
  list(total_required=n,completed=length(differences),valid=sum(is.finite(differences)),invalid=invalid,
       interval_suppressed=suppress,interval=interval,quantile_type=7,
       valid_only_conditioning=invalid>0,scientific_validation=FALSE)
}
write_bootstrap_checkpoint <- function(result,path) {
  write_atomic_json(result,path)
  checksum<-paste0(path,".sha256")
  temporary<-paste0(checksum,".partial.",Sys.getpid())
  writeLines(resample_sha(path),temporary)
  if(!file.rename(temporary,checksum))stop("checkpoint checksum atomic rename failed")
}
read_bootstrap_checkpoint <- function(path,fingerprint,replicate,mask,family_size) {
  checksum<-paste0(path,".sha256")
  if(!file.exists(checksum))stop("uncommitted checkpoint: checksum absent; preserve and inspect")
  saved_sha<-readLines(checksum,warn=FALSE)
  if(length(saved_sha)!=1L || !grepl("^[0-9a-f]{64}$",saved_sha) || !identical(saved_sha,resample_sha(path)))stop("bootstrap checkpoint checksum mismatch")
  result<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  if(!identical(result$binding_sha256,fingerprint) || !identical(as.integer(result$replicate),as.integer(replicate)) || !is.logical(result$valid) || length(result$valid)!=1L || is.na(result$valid))stop("unsafe replicate checkpoint")
  if(isTRUE(result$valid)) {
    if(!identical(as.integer(result$family_size),as.integer(family_size)))stop("checkpoint fixed-family size changed")
    keys<-unlist(lapply(c("B0","B1","G1"),function(m)paste0(m,"_top",c(20,50))))
    if(!setequal(names(result$rates),keys))stop("checkpoint paired-rate keys changed")
    for(m in c("B0","B1","G1"))for(k in c(20,50)){
      rate<-result$rates[[paste0(m,"_top",k)]]
      expected_n<-sum(mask$method==m & mask$k==k & mask$original_evaluable)
      values<-unlist(rate[c("agree","evaluable","selected","rate")],use.names=FALSE)
      if(length(values)!=4L || !is.numeric(values) || any(!is.finite(values)) || rate$evaluable!=expected_n || rate$selected!=k || rate$agree<0 || rate$agree>expected_n || rate$agree!=floor(rate$agree) || abs(rate$rate-rate$agree/expected_n)>1e-14)stop("checkpoint fixed-denominator/count/rate mismatch")
    }
    if(!is.numeric(result$difference) || length(result$difference)!=1L || !is.finite(result$difference) || abs(result$difference)>1 || abs(result$difference-(result$rates$G1_top20$rate-result$rates$B0_top20$rate))>1e-14)stop("checkpoint paired difference inconsistent")
  } else if(!is.character(result$reason) || length(result$reason)!=1L || !nzchar(result$reason))stop("invalid checkpoint needs explicit failure reason")
  result
}
resampling_main <- function(args) {
  if(length(args)%%2) stop("use --key value pairs")
  opt<-list(); for(i in seq(1,length(args),by=2)){k<-substring(args[i],3);if(!startsWith(args[i],"--")||!is.null(opt[[k]]))stop("bad option");opt[[k]]<-args[i+1]}
  mode<-opt$mode; if(!mode %in% c("bootstrap","lodo-fold"))stop("--mode bootstrap or lodo-fold required")
  cohort<-if(mode=="bootstrap")"heca" else "discovery"
  defaults<-list(counts=paste0("data/processed/",cohort,"/counts_stable.tsv.gz"),clinical=paste0("manifests/final_",cohort,"_clinical_selected.tsv"),config="config/statistical_analysis.json",`primary-model`="work/statistical/heca_primary/de/model.rds",`eval-mask`="work/statistical/heca_primary/evaluation/bootstrap_interface.tsv",`pre-effect-freeze`="receipts/statistical_pre_effect_freeze.json",`discovery-freeze`="receipts/discovery_rank_freeze.json")
  for(k in names(defaults))if(is.null(opt[[k]]))opt[[k]]<-defaults[[k]]
  if(is.null(opt[["output-dir"]]))stop("--output-dir required")
  synthetic<-identical(opt[["synthetic-test"]],"true")
  if(synthetic){
    for(p in c(opt$counts,opt$clinical,opt[["output-dir"]]))if(!startsWith(p,"work/synthetic/"))stop("synthetic paths must be work/synthetic/")
  } else {
    pre<-verify_resample_freeze(opt[["pre-effect-freeze"]],"STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT")
    discovery<-verify_resample_freeze(opt[["discovery-freeze"]],"DISCOVERY_FROZEN_PRE_EXTERNAL")
    for(p in c(opt$counts,opt$clinical,opt$config,"scripts/donor_de_lib.R"))if(!p %in% names(pre$artifacts))stop("input not frozen")
    if(!identical(opt$counts,defaults$counts) || !identical(opt$clinical,defaults$clinical) || !identical(opt$config,defaults$config))stop("production requires frozen input paths")
  }
  cfg<-jsonlite::fromJSON(opt$config)
  clinical<-read_resample_tsv(opt$clinical); ct<-read_resample_tsv(opt$counts)
  if(names(ct)[1]!="gene_id")stop("counts schema")
  counts<-as.matrix(ct[,-1,drop=FALSE]);storage.mode(counts)<-"double";rownames(counts)<-ct$gene_id
  if(!identical(colnames(counts),as.character(clinical$donor_id)) || anyDuplicated(clinical$donor_id) || anyDuplicated(rownames(counts)) || any(!grepl("^ENSG[0-9]{11}$",rownames(counts))) || any(!is.finite(counts)) || any(counts<0) || any(colSums(counts)<=0))stop("source input contract failed")
  if(cfg$cohorts[[cohort]]$integer_counts && any(counts!=floor(counts)))stop("source count integer contract failed")
  cells<-clinical[[cfg$cohorts[[cohort]]$cell_column]]
  if(any(!is.finite(cells)) || any(cells<cfg$design$min_cells))stop("original donor cells failed")
  paths<-c(opt$counts,opt$clinical,opt$config,"scripts/resample_donor_models.R","scripts/donor_de_lib.R","scripts/statistical_runtime_env.sh","requirements/statistical-runtime.explicit.txt")
  if(!synthetic)paths<-c(paths,opt[["pre-effect-freeze"]],opt[["discovery-freeze"]])
  out<-opt[["output-dir"]];dir.create(out,recursive=TRUE,showWarnings=FALSE)
  if(mode=="lodo-fold") {
    if(is.null(opt[["drop-donor"]]) || !opt[["drop-donor"]] %in% clinical$donor_id)stop("unknown drop donor")
    if(length(list.files(out,all.files=TRUE,no..=TRUE)))stop("preserve existing fold output")
    keep<-clinical$donor_id!=opt[["drop-donor"]];clinical<-clinical[keep,,drop=FALSE];counts<-counts[,keep,drop=FALSE]
    built<-rebuild_resample_design(clinical,cohort,cfg)
    # Frozen primary validator reapplies clinical minima and refilters; no repaired model.
    model<-fit_donor_model(counts,built$design,clinical,cohort,cfg)
    write_resample_tsv(model$table,file.path(out,"de_all.tsv"));write_resample_tsv(model$table[model$table$fit_estimable,],file.path(out,"de_fitted.tsv"))
    write_resample_tsv(data.frame(donor_id=rownames(built$design),built$design,check.names=FALSE),file.path(out,"design.tsv"))
    write_resample_tsv(data.frame(gene_id=rownames(counts),expression_pass=model$keep),file.path(out,"expression_filter.tsv"))
    saveRDS(model,file.path(out,"model.rds"),version=3)
    write_atomic_json(list(dropped_donor=opt[["drop-donor"]],audit=model$audit,design_receipt=built$receipt,input_sha256=setNames(lapply(paths,resample_sha),paths),warnings=model$warnings,synthetic=synthetic),file.path(out,"provenance.json"))
    return(invisible(NULL))
  }
  if(!identical(as.integer(table(factor(clinical$disease,levels=c("ENDO","non-ENDO")))),c(15L,7L)) && !synthetic)stop("bootstrap original cohort must be 15/7")
  if(!identical(as.integer(cfg$bootstrap$n),500L) || cfg$bootstrap$seed!=42 || !isTRUE(cfg$bootstrap$fixed_external_gene_family) || !isFALSE(cfg$bootstrap$refilter) || !isFALSE(cfg$bootstrap$redraw_invalid))stop("bootstrap contract differs")
  primary<-readRDS(opt[["primary-model"]]);family<-primary$table$gene_id
  if(!identical(family,rownames(counts)[primary$keep]) || !identical(primary$normalization$donor_id,colnames(counts)))stop("primary original family/order does not match source")
  mask<-read_resample_tsv(opt[["eval-mask"]]);validate_mask(mask)
  if(any(!mask$gene_id[mask$original_evaluable] %in% family))stop("originally evaluable mask outside original family")
  original_ix<-match(mask$gene_id,primary$table$gene_id)
  expected_eval<-!is.na(original_ix) & mask$discovery_sign!=0
  present<-which(!is.na(original_ix)); t0<-primary$table[original_ix[present],]
  expected_eval[present]<-expected_eval[present] & t0$fit_estimable & is.finite(t0$logFC) & is.finite(t0$t) & is.finite(t0$P.Value) & t0$logFC!=0
  if(!identical(as.logical(mask$original_evaluable),as.logical(expected_eval)))stop("mask differs from original external evaluability")
  if(!synthetic)for(m in c("B0","B1","G1"))for(k in c(20,50)){
    frozen<-read_resample_tsv(discovery$products[[paste0(m,"_top",k)]])
    x<-mask[mask$method==m & mask$k==k,]
    if(!identical(x$gene_id,frozen$gene_id) || !identical(as.numeric(x$discovery_sign),as.numeric(frozen$discovery_sign)))stop("mask lists/signs differ from frozen discovery")
  }
  paths<-c(paths,opt[["primary-model"]],opt[["eval-mask"]])
  primary_prov<-file.path(dirname(opt[["primary-model"]]),"provenance.json")
  if(!synthetic){
    pp<-jsonlite::fromJSON(primary_prov,simplifyVector=FALSE)
    if(!isTRUE(pp$scientific_effect_run) || pp$cohort!="heca")stop("HECA primary provenance absent")
    for(p in c(opt$counts,opt$clinical,opt$config))if(!identical(pp$input_sha256[[p]],resample_sha(p)))stop("primary input hash mismatch")
    for(p in names(pp$output_sha256))if(!identical(pp$output_sha256[[p]],resample_sha(p)))stop("HECA primary output hash mismatch")
    paths<-c(paths,primary_prov)
  }
  binding<-list(mode=mode,synthetic=synthetic,inputs=setNames(lapply(paths,resample_sha),paths),family=family,mask_rows=nrow(mask),R_version=as.character(getRversion()),edgeR=as.character(packageVersion("edgeR")),limma=as.character(packageVersion("limma")),statmod=as.character(packageVersion("statmod")))
  manifest<-file.path(out,"binding.rds");ledger_path<-file.path(out,"draw_ledger.tsv")
  if(file.exists(manifest)){
    saved<-readRDS(manifest)
    if(!identical(binding,saved$binding) || !file.exists(ledger_path) || resample_sha(ledger_path)!=saved$ledger_sha256)stop("unsafe resume binding/ledger changed")
    ledger<-read_resample_tsv(ledger_path)
    if(!identical(ledger,generate_draws(clinical,cfg)))stop("draw ledger not exact frozen RNG draws")
  } else {
    if(length(list.files(out,all.files=TRUE,no..=TRUE)))stop("nonempty unbound bootstrap output")
    ledger<-generate_draws(clinical,cfg)
    temp<-paste0(ledger_path,".partial");write_resample_tsv(ledger,temp);if(!file.rename(temp,ledger_path))stop("ledger rename failed")
    saveRDS(list(binding=binding,ledger_sha256=resample_sha(ledger_path)),paste0(manifest,".partial"),version=3);if(!file.rename(paste0(manifest,".partial"),manifest))stop("binding rename failed")
    write_atomic_json(list(binding=binding,ledger_sha256=resample_sha(ledger_path),ledger_saved_before_any_fit=TRUE,RNG="Mersenne-Twister/Inversion/Rejection",seed=42),file.path(out,"draw_freeze.json"))
  }
  fingerprint<-resample_sha(manifest)
  process<-list(pid=Sys.getpid(),host=Sys.info()[["nodename"]],started_utc=format(Sys.time(),tz="UTC",usetz=TRUE),binding_sha256=fingerprint,status="RUNNING",exit_status=NULL)
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
      result<-tryCatch({built<-rebuild_resample_design(clin,cohort,cfg);design_receipt<-built$receipt;fit<-fit_fixed_family(cts,built$design,cfg);ans<-evaluate_fixed_masks(fit$table,mask);ans$warnings<-fit$warnings;ans$family_size<-fit$family_size;ans},error=function(e)list(valid=FALSE,reason=conditionMessage(e)))
      result$worker_pid<-Sys.getpid();result$paired_replicate_valid<-isTRUE(result$valid)
      result$replicate<-r;result$binding_sha256<-fingerprint;result$elapsed_seconds<-proc.time()[["elapsed"]]-start
      result$unique_actual_people<-length(unique(rows$source_donor_id));result$unique_actual_ENDO<-length(unique(rows$source_donor_id[rows$disease=="ENDO"]));result$unique_actual_non_ENDO<-length(unique(rows$source_donor_id[rows$disease=="non-ENDO"]))
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
  summary<-bootstrap_summary(differences,cfg$bootstrap$n);summary$binding_sha256<-fingerprint;summary$failure_reasons<-as.list(table(reasons[nzchar(reasons)]));summary$synthetic<-synthetic
  write_resample_tsv(data.frame(replicate=seq_along(differences),valid=is.finite(differences),G1_minus_B0_top20=differences,reason=reasons),file.path(out,"paired_differences.tsv"))
  write_atomic_json(summary,file.path(out,"summary.json"));completed<-TRUE
  cat("BOOTSTRAP_COMPLETE",summary$valid,"valid",summary$invalid,"invalid\n")
}
if(sys.nframe()==0L)resampling_main(commandArgs(trailingOnly=TRUE))
