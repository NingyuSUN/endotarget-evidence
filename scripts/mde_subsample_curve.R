#!/usr/bin/env Rscript
# Minimum detectable effect analysis, Part B (docs/MDE_ANALYSIS_PROTOCOL.md section 4):
# donor-subsampling sensitivity curve on the phase positive-control signal. Each draw is a
# smaller study analyzed from scratch with the v1 estimator (refilter, TMM, voom, robust
# eBayes, BH over the subset family) under the primary and the reduced design. An identity
# check on all donors must reproduce the frozen full-cohort fit before any subsample is fitted.
# Draws are generated and written before any fit, never redrawn; invalid draws are retained.
source("scripts/mde_lib.R")
curve_main <- function(args) {
  opts <- mde_parse_opts(args)
  synthetic <- identical(opts[["synthetic-test"]],"true")
  required <- c("cohort","config","phase-config","mde-config","output-dir")
  if(!synthetic) required <- c(required,"freeze")
  if(!all(required %in% names(opts))) stop(paste("required:",paste(required,collapse=" ")))
  cohort <- opts$cohort
  if(!cohort %in% c("discovery","heca")) stop("cohort must be discovery or heca")
  cfg <- jsonlite::fromJSON(opts$config,simplifyVector=TRUE)
  pcfg <- jsonlite::fromJSON(opts[["phase-config"]],simplifyVector=TRUE)
  mcfg <- jsonlite::fromJSON(opts[["mde-config"]],simplifyVector=TRUE)
  sub <- mcfg$subsampling
  coefficient <- sub$coefficient
  if(!identical(coefficient,pcfg$coefficient)) stop("coefficient differs between configurations")
  started <- format(Sys.time(),tz="UTC",usetz=TRUE)
  freeze_check <- NULL; chain <- NULL
  if(synthetic) {
    for(k in c("counts","clinical","full-model")) if(is.null(opts[[k]]) || !startsWith(opts[[k]],"work/synthetic_mde/")) stop("synthetic mode requires work/synthetic_mde/ --counts, --clinical and --full-model")
    counts_path <- opts$counts; clinical_path <- opts$clinical; model_path <- opts[["full-model"]]
    levels <- if(is.null(opts$levels)) sub$levels_per_arm[[cohort]] else as.integer(strsplit(opts$levels,",")[[1]])
    draws <- if(is.null(opts$draws)) sub$draws_per_level else as.integer(opts$draws)
  } else {
    fz <- mde_verify_freeze(opts$freeze,"MDE_ANALYSIS_RULES_FROZEN_PRE_COMPUTATION"); freeze_check <- fz$check
    if(!identical(fz$freeze$artifacts[["config/mde_analysis.json"]],mde_sha(opts[["mde-config"]]))) stop("mde config is not the frozen one")
    entry <- mcfg$fits[[sub$reference_fit[[cohort]]]]
    chain <- mde_verify_fit_chain(entry)
    counts_path <- sub$counts[[cohort]]; clinical_path <- sub$clinical[[cohort]]; model_path <- entry$model
    for(p in c(counts_path,clinical_path)) if(!identical(fz$freeze$artifacts[[p]],mde_sha(p))) stop(paste("input differs from the frozen hash:",p))
    levels <- sub$levels_per_arm[[cohort]]; draws <- sub$draws_per_level
    if(!is.null(opts$levels) || !is.null(opts$draws)) stop("levels and draws are frozen for real runs")
  }
  out <- opts[["output-dir"]]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
  ckdir <- file.path(out,"checkpoints"); dir.create(ckdir,showWarnings=FALSE)
  workers <- if(is.null(opts$workers)) 1L else as.integer(opts$workers)
  if(is.na(workers) || workers<1 || workers>sub$max_workers) stop(paste("--workers must be 1..",sub$max_workers))

  # ---- inputs and contract ----
  counts <- mde_read_counts(counts_path); clinical <- mde_read_tsv(clinical_path)
  levels_phase <- c(pcfg$phase_factor$reference,pcfg$phase_factor$case)
  cell_col <- cfg$cohorts[[cohort]]$cell_column
  if(!all(c("donor_id","phase","disease",cell_col) %in% names(clinical)) || anyDuplicated(clinical$donor_id)) stop("invalid clinical ledger")
  if(!identical(colnames(counts),as.character(clinical$donor_id))) stop("exact donor order required")
  if(anyDuplicated(rownames(counts)) || any(!grepl("^ENSG[0-9]{11}$",rownames(counts)))) stop("invalid stable gene IDs")
  if(any(!is.finite(counts)) || any(counts<0) || any(colSums(counts)<=0)) stop("counts contract failed")
  if(isTRUE(pcfg$cohorts[[cohort]]$integer_counts) && any(counts!=floor(counts))) stop("integer count contract failed")
  if(any(!clinical$phase %in% levels_phase)) stop("unknown phase label")
  cells <- as.numeric(clinical[[cell_col]]); if(any(!is.finite(cells)) || any(cells<cfg$design$min_cells)) stop("donor cell threshold failed")
  sec_ix <- which(clinical$phase==pcfg$phase_factor$case); pro_ix <- which(clinical$phase==pcfg$phase_factor$reference)
  if(any(levels>min(length(sec_ix),length(pro_ix))) || any(levels<cfg$minimum_residual_df)) stop("a level exceeds an arm or is below the df minimum")
  group_all <- factor(clinical$phase,levels=levels_phase)

  # ---- reference sets from the frozen full-cohort fit ----
  full <- readRDS(model_path); ft <- full$table
  if(!identical(rownames(counts)[full$keep],as.character(ft$gene_id))) stop("full fit family does not match the source matrix")
  f_est <- as.logical(ft$fit_estimable)
  S <- as.character(ft$gene_id[f_est & !is.na(ft$adj.P.Val) & ft$adj.P.Val<0.05])
  S_sign <- sign(ft$logFC[match(S,ft$gene_id)])
  full_sign <- setNames(sign(ft$logFC[f_est]),ft$gene_id[f_est])
  top_full <- head(as.character(ft$gene_id[f_est][order(-abs(ft$t[f_est]))]),20)

  # ---- identity check: the subsampling estimator on all donors must reproduce the frozen fit ----
  ic <- sub$identity_check
  built_all <- mde_rebuild_design(clinical,cohort,cfg,pcfg,"primary")
  fit_all <- mde_fit_subset(counts,built_all$design,group_all,cfg,coefficient)
  ta <- fit_all$table
  same_family <- identical(as.character(ta$gene_id),as.character(ft$gene_id))
  d_logfc <- if(same_family) max(abs(ta$logFC-ft$logFC)) else Inf
  d_t <- if(same_family) max(abs(ta$t-ft$t)) else Inf
  d_p <- if(same_family) max(abs(log10(ta$P.Value)-log10(ft$P.Value))) else Inf
  identity <- list(same_family=same_family,max_abs_diff_logFC=d_logfc,max_abs_diff_t=d_t,max_abs_diff_log10_p=d_p,
                   tolerances=list(logFC=ic$tolerance_logFC,t=ic$tolerance_t,log10_p=ic$tolerance_log10_p),
                   design_columns=colnames(built_all$design),full_n_sig=length(S),full_top20=top_full,
                   passed=same_family && d_logfc<=ic$tolerance_logFC && d_t<=ic$tolerance_t && d_p<=ic$tolerance_log10_p)
  mde_write_json(identity,file.path(out,"identity_check.json"))
  if(!isTRUE(identity$passed)) stop("identity check failed; the subsampling estimator does not reproduce the frozen fit")
  rm(fit_all)

  # ---- draw ledger, written before any subsample is fitted ----
  generate <- function() {
    RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(sub$seed)
    do.call(rbind,lapply(levels,function(lv) do.call(rbind,lapply(seq_len(draws),function(r) {
      s <- sample(sec_ix,lv,replace=FALSE); p <- sample(pro_ix,lv,replace=FALSE)
      ix <- c(s,p)
      data.frame(level_n_per_arm=lv,replicate=r,draw_position=seq_along(ix),arm=c(rep(pcfg$phase_factor$case,lv),rep(pcfg$phase_factor$reference,lv)),
                 source_index=ix,source_donor_id=clinical$donor_id[ix],stringsAsFactors=FALSE)
    }))))
  }
  input_paths <- c(counts_path,clinical_path,model_path,opts$config,opts[["phase-config"]],opts[["mde-config"]],"scripts/mde_lib.R","scripts/mde_subsample_curve.R","scripts/statistical_runtime_env.sh")
  if(!synthetic) input_paths <- c(input_paths,opts$freeze)
  binding <- list(cohort=cohort,coefficient=coefficient,levels=levels,draws=draws,seed=sub$seed,replace=FALSE,
                  inputs=setNames(lapply(input_paths,mde_sha),input_paths),software=mde_software(),synthetic=synthetic)
  ledger_path <- file.path(out,"draw_ledger.tsv"); binding_path <- file.path(out,"draw_freeze.json")
  ledger <- generate()
  if(file.exists(binding_path)) {
    saved <- jsonlite::fromJSON(binding_path,simplifyVector=FALSE)
    if(!identical(saved$ledger_sha256,mde_sha(ledger_path))) stop("ledger changed since it was frozen")
    tmp <- tempfile(); mde_write_tsv(ledger,tmp)
    if(!identical(mde_sha(tmp),saved$ledger_sha256)) stop("regenerated draws differ from the frozen ledger")
    for(p in names(saved$binding$inputs)) if(!identical(saved$binding$inputs[[p]],binding$inputs[[p]])) stop(paste("input changed since the ledger freeze:",p))
    unlink(tmp)
  } else {
    if(length(list.files(ckdir))) stop("checkpoints exist without a draw freeze")
    mde_write_tsv(ledger,ledger_path)
    mde_write_json(list(binding=binding,ledger_sha256=mde_sha(ledger_path),ledger_saved_before_any_fit=TRUE,
                        rng="Mersenne-Twister/Inversion/Rejection",seed=sub$seed,arm_sizes=list(secretory=length(sec_ix),proliferative=length(pro_ix))),binding_path)
  }
  fingerprint <- mde_sha(binding_path)

  # ---- per-draw analysis ----
  metrics_of <- function(tab,receipt) {
    est <- as.logical(tab$fit_estimable); q <- tab$adj.P.Val
    sig <- est & !is.na(q) & q<0.05
    idx <- match(S,tab$gene_id); inS <- !is.na(idx)
    hit <- rep(FALSE,length(S)); hit_any <- rep(FALSE,length(S))
    if(any(inS)) { hit_any[inS] <- sig[idx[inS]]; hit[inS] <- sig[idx[inS]] & sign(tab$logFC[idx[inS]])==S_sign[inS] }
    dsig <- as.character(tab$gene_id[sig]); fs <- full_sign[dsig]; keep <- !is.na(fs) & fs!=0
    disagree <- if(any(keep)) mean(sign(tab$logFC[sig][keep])!=fs[keep]) else NA_real_
    ord <- order(-abs(tab$t[est])); top <- head(as.character(tab$gene_id[est][ord]),20)
    med_se <- median(tab$moderated_SE[est]); med_df <- median(tab$df_total[est]); family <- nrow(tab)
    a100 <- mcfg$bh_q*100/family
    list(family_size=family,estimable_genes=sum(est),n_sig=sum(sig),n_sig_up=sum(sig & tab$logFC>0),n_sig_down=sum(sig & tab$logFC<0),
         recall_full=if(length(S)) mean(hit) else NA_real_,recall_full_any_sign=if(length(S)) mean(hit_any) else NA_real_,
         sign_disagreement=disagree,n_sig_with_full_sign=sum(keep),top20_overlap=length(intersect(top,top_full)),
         median_se=med_se,median_df_total=med_df,
         typical_mde_alpha_0.05=mde_solve(med_se,med_df,mcfg$nominal_alpha,power=mcfg$power)$mde,
         typical_mde_bh_k100=if(a100<1) mde_solve(med_se,med_df,a100,power=mcfg$power)$mde else NA_real_,
         residual_df=receipt$residual_df,rank=receipt$rank,n_batch_levels=receipt$n_batch_levels,
         dropped_constant_nuisance=paste(receipt$dropped_constant_nuisance,collapse=";"))
  }
  tasks <- unique(ledger[,c("level_n_per_arm","replicate")])
  fit_one <- function(i) {
    lv <- tasks$level_n_per_arm[i]; r <- tasks$replicate[i]
    path <- file.path(ckdir,sprintf("level_%02d_rep_%03d.json",lv,r))
    if(file.exists(path)) {
      res <- jsonlite::fromJSON(path,simplifyVector=FALSE)
      if(!identical(res$binding_sha256,fingerprint) || res$level_n_per_arm!=lv || res$replicate!=r) stop(paste("checkpoint does not belong to this run:",path))
      cat("CURVE_CHECKPOINT reused",lv,r,"\n"); return(res)
    }
    rows <- ledger[ledger$level_n_per_arm==lv & ledger$replicate==r,]
    ix <- rows$source_index
    clin <- clinical[ix,,drop=FALSE]; cts <- counts[,ix,drop=FALSE]
    group <- factor(clin$phase,levels=levels_phase)
    res <- list(level_n_per_arm=lv,replicate=r,binding_sha256=fingerprint,donors=as.character(clin$donor_id),designs=list())
    for(variant in c("primary","reduced")) {
      t0 <- proc.time()[["elapsed"]]
      ans <- tryCatch({
        built <- mde_rebuild_design(clin,cohort,cfg,pcfg,variant)
        fit <- mde_fit_subset(cts,built$design,group,cfg,coefficient)
        c(list(valid=TRUE,reason=""),metrics_of(fit$table,built$receipt),list(warnings=length(fit$warnings)))
      },error=function(e) list(valid=FALSE,reason=conditionMessage(e)))
      ans$elapsed_seconds <- proc.time()[["elapsed"]]-t0
      res$designs[[variant]] <- ans
    }
    tmp <- paste0(path,".partial"); jsonlite::write_json(res,tmp,auto_unbox=TRUE,null="null",na="null",digits=NA)
    if(!file.rename(tmp,path)) stop("checkpoint rename failed")
    cat("CURVE_CHECKPOINT",lv,r,"primary",isTRUE(res$designs$primary$valid),"reduced",isTRUE(res$designs$reduced$valid),"\n")
    res
  }
  Sys.setenv(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
  mde_write_json(list(pid=Sys.getpid(),started_utc=started,workers=workers,status="RUNNING",tasks=nrow(tasks)),file.path(out,"process.json"))
  results <- if(workers==1L) lapply(seq_len(nrow(tasks)),fit_one) else parallel::mclapply(seq_len(nrow(tasks)),fit_one,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  for(res in results) if(inherits(res,"try-error") || is.null(res$designs)) stop("a worker failed; completed checkpoints are preserved")

  # ---- assemble ----
  metric_names <- c("family_size","estimable_genes","n_sig","n_sig_up","n_sig_down","recall_full","recall_full_any_sign","sign_disagreement",
                    "n_sig_with_full_sign","top20_overlap","median_se","median_df_total","typical_mde_alpha_0.05","typical_mde_bh_k100",
                    "residual_df","rank","n_batch_levels")
  rows <- list()
  for(res in results) for(variant in c("primary","reduced")) {
    d <- res$designs[[variant]]
    row <- data.frame(cohort=cohort,level_n_per_arm=res$level_n_per_arm,replicate=res$replicate,design=variant,valid=isTRUE(d$valid),
                      reason=if(is.null(d$reason)) "" else d$reason,stringsAsFactors=FALSE)
    for(mn in metric_names) row[[mn]] <- if(isTRUE(d$valid) && !is.null(d[[mn]])) as.numeric(d[[mn]]) else NA_real_
    row$dropped_constant_nuisance <- if(isTRUE(d$valid)) d$dropped_constant_nuisance else NA_character_
    row$elapsed_seconds <- d$elapsed_seconds
    rows[[length(rows)+1]] <- row
  }
  draws_tab <- do.call(rbind,rows)
  draws_tab <- draws_tab[order(draws_tab$design,draws_tab$level_n_per_arm,draws_tab$replicate),]
  pcts <- sub$summary_percentiles
  summ <- list()
  for(variant in c("primary","reduced")) for(lv in levels) {
    x <- draws_tab[draws_tab$design==variant & draws_tab$level_n_per_arm==lv,]
    v <- x[x$valid,]
    row <- data.frame(cohort=cohort,design=variant,level_n_per_arm=lv,n_donors=2*lv,draws=nrow(x),valid=nrow(v),
                      invalid_reasons=paste(names(table(x$reason[!x$valid])),table(x$reason[!x$valid]),sep="=",collapse=";"),stringsAsFactors=FALSE)
    for(mn in c("n_sig","recall_full","sign_disagreement","top20_overlap","typical_mde_alpha_0.05","typical_mde_bh_k100","median_se","family_size","residual_df")) {
      qv <- if(nrow(v)) quantile(v[[mn]],pcts,type=sub$quantile_type,names=FALSE,na.rm=TRUE) else rep(NA_real_,length(pcts))
      row[[paste0(mn,"_q025")]] <- qv[1]; row[[paste0(mn,"_median")]] <- qv[2]; row[[paste0(mn,"_q975")]] <- qv[3]
    }
    summ[[length(summ)+1]] <- row
  }
  summary_tab <- do.call(rbind,summ)
  mde_write_tsv(draws_tab,file.path(out,"curve_draws.tsv"))
  mde_write_tsv(summary_tab,file.path(out,"curve_summary.tsv"))
  figure <- tryCatch({
    png(file.path(out,"curve.png"),width=1800,height=600,res=150); par(mfrow=c(1,3),mar=c(4.5,4.5,3,1))
    panel <- function(mn,ylab,main,ylim=NULL,ref=NULL) {
      plot(NA,xlim=range(levels)+c(-0.5,0.5),ylim=if(is.null(ylim)) range(c(summary_tab[[paste0(mn,"_q025")]],summary_tab[[paste0(mn,"_q975")]],ref),na.rm=TRUE) else ylim,xlab="donors per phase arm",ylab=ylab,main=main)
      cols <- c(primary="#1b4f72",reduced="#c0392b"); off <- c(primary=-0.12,reduced=0.12)
      for(variant in c("primary","reduced")) { s <- summary_tab[summary_tab$design==variant,]; ok <- s$valid>0
        arrows(s$level_n_per_arm[ok]+off[variant],s[[paste0(mn,"_q025")]][ok],s$level_n_per_arm[ok]+off[variant],s[[paste0(mn,"_q975")]][ok],angle=90,code=3,length=0.03,col=cols[variant])
        lines(s$level_n_per_arm[ok]+off[variant],s[[paste0(mn,"_median")]][ok],type="b",pch=19,col=cols[variant]) }
      if(!is.null(ref)) abline(h=ref,lty=3,col="grey40")
      legend("topleft",legend=c("primary design","reduced design"),col=cols,pch=19,lty=1,bty="n",cex=0.9)
    }
    panel("n_sig","genes at BH q<0.05",paste(cohort,": detected genes"),ref=length(S))
    panel("recall_full","recall of the full-cohort program",paste(cohort,": recall"),ylim=c(0,1))
    panel("typical_mde_alpha_0.05","typical-gene MDE (log2, alpha 0.05)",paste(cohort,": typical MDE"))
    dev.off(); "curve.png"
  },error=function(e){try(dev.off(),silent=TRUE);paste("figure_not_rendered:",conditionMessage(e))})
  for(p in input_paths) if(!identical(mde_sha(p),binding$inputs[[p]])) stop(paste("input changed during the run:",p))
  valid_counts <- lapply(c("primary","reduced"),function(v) setNames(as.list(summary_tab$valid[summary_tab$design==v]),paste0("n",levels)))
  names(valid_counts) <- c("primary","reduced")
  summary <- list(cohort=cohort,coefficient=coefficient,synthetic=synthetic,levels_per_arm=levels,draws_per_level=draws,seed=sub$seed,
                  arm_sizes=list(secretory=length(sec_ix),proliferative=length(pro_ix)),full_cohort=list(n_sig=length(S),family_size=nrow(ft),top20=top_full),
                  identity_check=identity[c("passed","max_abs_diff_logFC","max_abs_diff_t","max_abs_diff_log10_p")],valid_draws=valid_counts,
                  total_fits=2*nrow(tasks),workers=workers,figure=figure)
  mde_write_json(summary,file.path(out,"summary.json"))
  outputs <- setdiff(list.files(out,full.names=TRUE,recursive=FALSE),c(file.path(out,"provenance.json"),file.path(out,"checkpoints"),file.path(out,"process.json")))
  prov <- list(started_utc=started,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),analysis="MDE_PART_B_SUBSAMPLING_CURVE",protocol=mcfg$protocol,
               cohort=cohort,coefficient=coefficient,synthetic=synthetic,data_kind=if(synthetic)"synthetic_test_only" else "public_cohort_positive_control_matrices",
               reference_model=list(path=model_path,sha256=mde_sha(model_path),chain=chain),
               freeze=if(synthetic) NULL else list(path=opts$freeze,sha256=mde_sha(opts$freeze),phase="MDE_ANALYSIS_RULES_FROZEN_PRE_COMPUTATION",artifact_check=freeze_check),
               input_sha256=binding$inputs,binding_sha256=fingerprint,levels_per_arm=levels,draws_per_level=draws,workers=workers,
               identity_check=identity,software=mde_software(),output_sha256=setNames(lapply(outputs,mde_sha),outputs),
               checkpoints=length(list.files(ckdir)),v1_artifacts_modified=FALSE,new_disease_fits=FALSE,
               scientific_scope=if(synthetic)"invented fixture, no biological content" else "descriptive degradation of a known program with donor number under the frozen estimator; not a power calculation for disease effects")
  mde_write_json(prov,file.path(out,"provenance.json"))
  mde_write_json(list(pid=Sys.getpid(),started_utc=started,finished_utc=prov$completed_utc,workers=workers,status="COMPLETED",tasks=nrow(tasks)),file.path(out,"process.json"))
  cat("MDE_CURVE_COMPLETE",cohort,"tasks",nrow(tasks),"valid_primary",sum(draws_tab$valid[draws_tab$design=="primary"]),
      "valid_reduced",sum(draws_tab$valid[draws_tab$design=="reduced"]),"\n")
}
if(sys.nframe()==0L) curve_main(commandArgs(trailingOnly=TRUE))
