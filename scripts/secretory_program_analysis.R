#!/usr/bin/env Rscript
# Secretory-phase program-level disease contrast (docs/SECRETORY_PROGRAM_PROTOCOL.md), one cohort
# per run. Verifies the freeze and every pinned input, reproduces the frozen expression filter,
# TMM factors and design of the phase positive-control fit (identity checks), scores every donor
# on the other cohort's frozen phase program, fits the fixed donor-level models, the sensitivities
# and the gene-level interaction fit with its set test, and writes every result with provenance.
# Real runs require the frozen rules; synthetic runs take a work/synthetic_secretory/ inputs file.
source("scripts/secretory_program_lib.R")
opts <- mde_parse_opts(commandArgs(trailingOnly=TRUE))
synthetic <- identical(opts[["synthetic-test"]],"true")
required <- c("cohort","config","phase-config","sp-config","output-dir",if(synthetic) "inputs" else "freeze")
if(!all(required %in% names(opts))) stop(paste("required:",paste(required,collapse=" ")))
cohort <- opts$cohort; if(!cohort %in% c("discovery","heca")) stop("cohort must be discovery or heca")
cfg <- jsonlite::fromJSON(opts$config,simplifyVector=TRUE)
pcfg <- jsonlite::fromJSON(opts[["phase-config"]],simplifyVector=TRUE)
scfg <- jsonlite::fromJSON(opts[["sp-config"]],simplifyVector=TRUE)
mcfg <- scfg$models; gate_cfg <- scfg$validity_gate; lab <- scfg$evidence_labels
FREEZE_PHASE <- "SECRETORY_PROGRAM_RULES_FROZEN_PRE_ANALYSIS"
started <- format(Sys.time(),tz="UTC",usetz=TRUE)
freeze_check <- NULL; freeze_sha <- NULL; input_hashes <- list()
if(synthetic) {
  if(!startsWith(opts$inputs,"work/synthetic_secretory/")) stop("synthetic mode requires a work/synthetic_secretory/ inputs file")
  inp <- jsonlite::fromJSON(opts$inputs,simplifyVector=TRUE)
  if(!isTRUE(inp$synthetic) || !identical(inp$cohort,cohort)) stop("inputs file is not the synthetic one for this cohort")
  for(k in c("counts","clinical","phase_fit_table","expression_filter","normalization","design","signature_table")) input_hashes[[inp[[k]]]] <- mde_sha(inp[[k]])
} else {
  fz <- mde_verify_freeze(opts$freeze,FREEZE_PHASE); freeze_check <- fz$check; freeze_sha <- mde_sha(opts$freeze)
  if(!identical(fz$freeze$artifacts[["config/secretory_program.json"]],mde_sha(opts[["sp-config"]]))) stop("secretory program config is not the frozen one")
  inp <- scfg$inputs[[cohort]]; other <- scfg$signature$source[[cohort]]
  inp$signature_table <- scfg$inputs[[other]]$phase_fit_table; inp$signature_table_sha256 <- scfg$inputs[[other]]$phase_fit_table_sha256
  for(k in c("counts","clinical","phase_fit_table","expression_filter","normalization","design","signature_table")) {
    p <- inp[[k]]; pinned <- inp[[paste0(k,"_sha256")]]
    if(is.null(p) || is.null(pinned) || !file.exists(p)) stop(paste("missing pinned input:",k))
    actual <- mde_sha(p); if(!identical(actual,pinned)) stop(paste("input differs from the pinned hash:",p))
    input_hashes[[p]] <- actual
  }
}
out <- opts[["output-dir"]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("output directory is not empty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

# ---- inputs and contract ----
counts <- mde_read_counts(inp$counts); clinical <- mde_read_tsv(inp$clinical)
PRO <- pcfg$phase_factor$reference; SEC <- pcfg$phase_factor$case; levels_phase <- c(PRO,SEC)
cell_col <- cfg$cohorts[[cohort]]$cell_column; batch_col <- cfg$cohorts[[cohort]]$batch_column
sub_col <- inp$sub_phase_column; sub_levels <- inp$sub_phase_levels; depth_col <- inp$depth_column
need <- c("donor_id","phase","disease",cfg$design$age_column,cfg$design$adenomyosis_column,batch_col,cell_col,sub_col,depth_col)
if(!all(need %in% names(clinical)) || anyDuplicated(clinical$donor_id)) stop("invalid clinical ledger")
if(!identical(colnames(counts),as.character(clinical$donor_id))) stop("exact donor order required")
if(anyDuplicated(rownames(counts)) || any(!grepl("^ENSG[0-9]{11}$",rownames(counts)))) stop("invalid stable gene IDs")
if(any(!is.finite(counts)) || any(counts<0) || any(colSums(counts)<=0)) stop("counts contract failed")
if(isTRUE(pcfg$cohorts[[cohort]]$integer_counts) && any(counts!=floor(counts))) stop("integer count contract failed")
if(any(!clinical$phase %in% levels_phase) || any(!clinical$disease %in% c("ENDO","non-ENDO"))) stop("unknown phase or disease label")
cells <- as.numeric(clinical[[cell_col]]); if(any(!is.finite(cells)) || any(cells<cfg$design$min_cells)) stop("donor cell threshold failed")
sec_ix <- which(clinical$phase==SEC); pro_ix <- which(clinical$phase==PRO)
if(length(sec_ix)<pcfg$cohorts[[cohort]]$minimum_per_phase_group || length(pro_ix)<pcfg$cohorts[[cohort]]$minimum_per_phase_group) stop("phase group minimum failed")
if(length(unique(clinical$disease[sec_ix]))!=2) stop("both disease groups are required among the secretory donors")
if(any(!clinical[[sub_col]][sec_ix] %in% sub_levels)) stop("secretory sub-phase label outside the configured levels")
group_all <- factor(clinical$phase,levels=levels_phase)
depth <- as.numeric(clinical[[depth_col]]); if(any(!is.finite(depth))) stop("invalid depth column")

# ---- identity checks against the frozen positive-control fit ----
ic <- scfg$identity_checks
fn <- sp_filter_norm(counts,group_all,cfg)
frozen_filter <- mde_read_tsv(inp$expression_filter)
if(!identical(as.character(frozen_filter$gene_id),rownames(counts))) stop("frozen filter gene order differs from the counts")
filter_identical <- identical(unname(fn$keep),as.logical(frozen_filter$expression_pass))
frozen_norm <- mde_read_tsv(inp$normalization)
if(!identical(as.character(frozen_norm$donor_id),colnames(counts))) stop("frozen normalization donor order differs")
rel <- function(a,b) max(abs(a-b)/pmax(abs(b),1e-300))
norm_diff <- rel(fn$norm_factor,as.numeric(frozen_norm$norm_factor)); lib_diff <- rel(fn$filtered_library_size,as.numeric(frozen_norm$filtered_library_size))
norm_identical <- norm_diff<=ic$tolerance_relative && lib_diff<=ic$tolerance_relative
built <- mde_rebuild_design(clinical,cohort,cfg,pcfg,"primary")
frozen_design <- mde_read_tsv(inp$design); D <- as.matrix(frozen_design[,-1,drop=FALSE]); storage.mode(D) <- "double"; rownames(D) <- frozen_design$donor_id
design_identical <- identical(colnames(D),colnames(built$design)) && identical(rownames(D),rownames(built$design)) && max(abs(D-built$design))<=ic$tolerance_design
family <- rownames(counts)[fn$keep]
phase_table <- mde_read_tsv(inp$phase_fit_table)
family_identical <- identical(as.character(phase_table$gene_id),family)
identity <- list(expression_filter_identical=filter_identical,normalization_max_relative_diff=norm_diff,filtered_library_max_relative_diff=lib_diff,
                 normalization_identical=norm_identical,design_identical=design_identical,design_columns=colnames(built$design),
                 own_table_family_identical=family_identical,family_size=length(family),
                 passed=filter_identical && norm_identical && design_identical && family_identical)
mde_write_json(identity,file.path(out,"identity_check.json"))
if(!isTRUE(identity$passed)) stop("identity checks against the frozen positive-control fit failed")
logcpm <- sp_logcpm(fn$y); if(!identical(rownames(logcpm),family)) stop("log-CPM rows differ from the family")

# ---- signature and donor scores ----
signature_table <- mde_read_tsv(inp$signature_table)
sig <- sp_signature(signature_table,family,scfg$signature$q_threshold)
own_sig <- sp_signature(phase_table,family,scfg$signature$q_threshold)
scores <- sp_scores(logcpm,sig,pro_ix); own_scores <- sp_scores(logcpm,own_sig,pro_ix)
donors <- data.frame(donor_id=clinical$donor_id,phase=clinical$phase,sub_phase=clinical[[sub_col]],disease=clinical$disease,
                     age=as.numeric(clinical[[cfg$design$age_column]]),batch=as.character(clinical[[batch_col]]),
                     adenomyosis=as.numeric(clinical[[cfg$design$adenomyosis_column]]),depth=depth,
                     projection=scores$projection,mean_signed_shift=scores$mean_signed_shift,own_program_projection=own_scores$projection,
                     stringsAsFactors=FALSE,check.names=FALSE)
mde_write_tsv(donors,file.path(out,"donor_scores.tsv"))
mde_write_tsv(sig,file.path(out,"signature_genes.tsv"))
grp <- paste(donors$phase,donors$disease)
group_means <- lapply(split(donors$projection,grp),function(v) list(n=length(v),mean=mean(v),sd=if(length(v)>1) sd(v) else NA_real_,median=median(v)))
sec_sd <- sd(donors$projection[sec_ix])

# ---- M0 validity gate ----
m0 <- sp_lm(donors$projection,built$design,"SEC_minus_PRO",mcfg)
auc <- sp_auc(donors$projection,donors$phase==SEC)
gate <- list(model=m0,auc=auc,pass=isTRUE(m0$status=="ESTIMATED" && sign(m0$estimate)==gate_cfg$phase_coefficient_sign && m0$p<gate_cfg$phase_coefficient_p_max && auc>=gate_cfg$phase_separation_auc_min))

# ---- M1 primary and its MDE ----
with_mde <- function(m) { if(identical(m$status,"ESTIMATED")) m$mde_log <- mde_solve(m$se,m$residual_df,scfg$mde$alpha,scfg$mde$power)$mde[1]; m }
d1 <- sp_add_interaction(built$design)
m1 <- with_mde(sp_lm(donors$projection,d1,"SEC_x_ENDO",mcfg))
if(identical(m1$status,"ESTIMATED")) m1$estimate_in_secretory_sd_units <- m1$estimate/sec_sd
# ---- M2 secretory-only and the sub-phase variant ----
clin_sec <- clinical[sec_ix,,drop=FALSE]
d2 <- sp_secretory_design(clin_sec,cfg); m2 <- with_mde(sp_lm(donors$projection[sec_ix],d2,"ENDO_minus_non_ENDO",mcfg))
spc <- sp_subphase_contrasts(clin_sec[[sub_col]],sub_levels)
d2s <- sp_secretory_design(clin_sec,cfg,spc$contrasts); m2s <- with_mde(sp_lm(donors$projection[sec_ix],d2s,"ENDO_minus_non_ENDO",mcfg))
m2s$sub_phase_reference <- spc$reference; m2s$sub_phase_levels_present <- spc$levels_present

# ---- sensitivities of M1 ----
refit_subset <- function(ix,y) {
  tryCatch({
    b <- mde_rebuild_design(clinical[ix,,drop=FALSE],cohort,cfg,pcfg,"primary")
    m <- sp_lm(y[ix],sp_add_interaction(b$design),"SEC_x_ENDO",mcfg); m$dropped_constant_nuisance <- b$receipt$dropped_constant_nuisance; m
  },error=function(e) list(status="NOT_ESTIMABLE",reason=conditionMessage(e),n=length(ix)))
}
sens <- list()
sens$reduced_design <- tryCatch({ b <- mde_rebuild_design(clinical,cohort,cfg,pcfg,"reduced"); with_mde(sp_lm(donors$projection,sp_add_interaction(b$design),"SEC_x_ENDO",mcfg)) },
                                error=function(e) list(status="NOT_ESTIMABLE",reason=conditionMessage(e)))
sens$secondary_score <- with_mde(sp_lm(donors$mean_signed_shift,d1,"SEC_x_ENDO",mcfg))
sens$own_cohort_program <- with_mde(sp_lm(donors$own_program_projection,d1,"SEC_x_ENDO",mcfg)); sens$own_cohort_program$circular <- TRUE
keep_depth <- which(depth>=mcfg$depth_minimum)
sens$depth_minimum <- if(length(keep_depth)==nrow(clinical)) list(status="IDENTICAL_TO_PRIMARY",excluded=0L) else c(with_mde(refit_subset(keep_depth,donors$projection)),list(excluded=nrow(clinical)-length(keep_depth),excluded_donors=clinical$donor_id[-keep_depth]))
keep_adeno <- which(as.numeric(clinical[[cfg$design$adenomyosis_column]])==0)
sens$adenomyosis_excluded <- if(length(keep_adeno)==nrow(clinical)) list(status="IDENTICAL_TO_PRIMARY",excluded=0L) else c(with_mde(refit_subset(keep_adeno,donors$projection)),list(excluded=nrow(clinical)-length(keep_adeno)))
loo <- do.call(rbind,lapply(seq_len(nrow(clinical)),function(i) {
  m <- refit_subset(setdiff(seq_len(nrow(clinical)),i),donors$projection)
  data.frame(left_out=clinical$donor_id[i],phase=clinical$phase[i],disease=clinical$disease[i],status=m$status,
             theta=if(identical(m$status,"ESTIMATED")) m$estimate else NA_real_,se=if(identical(m$status,"ESTIMATED")) m$se else NA_real_,
             p=if(identical(m$status,"ESTIMATED")) m$p else NA_real_,reason=if(identical(m$status,"ESTIMATED")) "" else m$reason,stringsAsFactors=FALSE)
}))
mde_write_tsv(loo,file.path(out,"leave_one_donor_out.tsv"))
loo_summary <- list(fits=nrow(loo),estimated=sum(loo$status=="ESTIMATED"),theta_min=suppressWarnings(min(loo$theta,na.rm=TRUE)),theta_max=suppressWarnings(max(loo$theta,na.rm=TRUE)),
                    max_abs_change=if(identical(m1$status,"ESTIMATED")) suppressWarnings(max(abs(loo$theta-m1$estimate),na.rm=TRUE)) else NA_real_,
                    most_influential=if(identical(m1$status,"ESTIMATED") && any(!is.na(loo$theta))) loo$left_out[which.max(abs(loo$theta-m1$estimate))] else NA_character_,
                    all_same_sign_as_primary=if(identical(m1$status,"ESTIMATED")) all(sign(loo$theta[!is.na(loo$theta)])==sign(m1$estimate)) else NA)

# ---- gene-level interaction fit and set test ----
gl <- mde_fit_subset(counts,d1,group_all,cfg,"SEC_x_ENDO")
if(!identical(as.character(gl$table$gene_id),family)) stop("gene-level family differs from the recomputed family")
vo <- sp_voom(fn$y,d1,cfg); if(!identical(rownames(vo),family)) stop("voom rows differ from the family")
in_fam <- which(sig$in_family); idx <- match(sig$gene_id[in_fam],family)
w_all <- numeric(length(family)); w_all[idx] <- sig$sign[in_fam]
fry <- sp_fry(vo,d1,"SEC_x_ENDO",idx,w_all)
lf <- gl$table$logFC[idx]; b <- sig$beta[in_fam]
gene_level <- list(family_size=gl$family_size,residual_df=gl$table$df_residual[1],estimable_genes=sum(gl$table$fit_estimable),
                   genes_BH_q_lt_0_05=sum(gl$table$adj.P.Val<0.05,na.rm=TRUE),max_abs_t=max(abs(gl$table$t[gl$table$fit_estimable])),
                   signature_genes_in_family=length(idx),fry=fry,
                   interaction_projection_on_signature=sum(b*lf)/sum(b^2),mean_signed_interaction_log2fc=mean(sig$sign[in_fam]*lf),
                   share_signature_genes_interaction_opposes_beta=mean(sign(lf)*sig$sign[in_fam]<0),
                   estimator_warnings=gl$warnings)
mde_write_tsv(gl$table,file.path(out,"interaction_de_all.tsv"))
top <- gl$table[gl$table$fit_estimable,]; top <- top[order(-abs(top$t)),][seq_len(min(20,nrow(top))),c("gene_id","logFC","t","P.Value","adj.P.Val","moderated_SE")]
top$in_signature <- top$gene_id %in% sig$gene_id[in_fam]; top$note <- "descriptive; not a shortlist"
mde_write_tsv(top,file.path(out,"interaction_top20_by_abs_t.tsv"))

# ---- evidence label for this cohort (fixed rules) ----
label <- if(!gate$pass) "NOT_INTERPRETABLE" else if(!identical(m1$status,"ESTIMATED")) "NOT_ESTIMABLE" else if(m1$p<lab$per_cohort_alpha && m1$estimate<0) "attenuation_detected" else if(m1$p<lab$per_cohort_alpha && m1$estimate>0) "enhancement_detected" else "no_detectable_difference"

# ---- figure ----
fig <- tryCatch({
  png(file.path(out,"scores.png"),width=1600,height=650,res=150); par(mfrow=c(1,2),mar=c(7,4.5,3,1))
  glev <- c(paste(PRO,"non-ENDO"),paste(PRO,"ENDO"),paste(SEC,"non-ENDO"),paste(SEC,"ENDO")); gi <- match(grp,glev)
  pch_levels <- c("proliferative",sub_levels); pchs <- c(1,15,16,17)[match(ifelse(donors$phase==PRO,"proliferative",donors$sub_phase),pch_levels)]
  for(k in c("projection","mean_signed_shift")) {
    y <- donors[[k]]; set.seed(1); x <- gi+runif(length(gi),-0.18,0.18)
    plot(x,y,xaxt="n",xlab="",ylab=k,pch=pchs,col=ifelse(donors$disease=="ENDO","firebrick","steelblue"),main=paste(cohort,"-",k),xlim=c(0.5,4.5))
    axis(1,at=1:4,labels=glev,cex.axis=0.75,las=2); abline(h=0,lty=3)
    for(g in 1:4) { v <- y[gi==g]; if(length(v)) segments(g-0.25,mean(v),g+0.25,mean(v),lwd=2) }
  }
  legend("topleft",legend=c("non-ENDO","ENDO",pch_levels),col=c("steelblue","firebrick",rep("black",length(pch_levels))),pch=c(16,16,1,15,16,17),cex=0.7,bty="n")
  dev.off(); "scores.png"
},error=function(e){try(dev.off(),silent=TRUE);paste("figure_not_rendered:",conditionMessage(e))})

# ---- summary and provenance ----
cross <- function(a,b) { d <- as.data.frame(table(a,b),stringsAsFactors=FALSE); names(d) <- c("level","disease","n"); unname(split(d,seq_len(nrow(d)))) }
summary <- list(analysis="SECRETORY_PROGRAM_DISEASE_CONTRAST",cohort=cohort,synthetic=synthetic,
                donors=list(n=nrow(clinical),by_phase_disease=as.list(table(grp)),
                            secretory_sub_phase_by_disease=lapply(cross(donors$sub_phase[sec_ix],donors$disease[sec_ix]),as.list),
                            secretory_batch_by_disease=lapply(cross(donors$batch[sec_ix],donors$disease[sec_ix]),as.list),
                            secretory_age_mean_by_disease=as.list(tapply(donors$age[sec_ix],donors$disease[sec_ix],mean)),
                            secretory_adenomyosis_by_disease=as.list(tapply(donors$adenomyosis[sec_ix],donors$disease[sec_ix],sum)),
                            secretory_depth_min_by_disease=as.list(tapply(depth[sec_ix],donors$disease[sec_ix],min))),
                signature=list(source_cohort=if(synthetic) "synthetic_other" else scfg$signature$source[[cohort]],q_threshold=scfg$signature$q_threshold,
                               genes=nrow(sig),in_family=sum(sig$in_family),up=sum(sig$sign[sig$in_family]>0),down=sum(sig$sign[sig$in_family]<0),
                               median_abs_beta_in_family=median(abs(sig$beta[sig$in_family])),own_program_genes=nrow(own_sig),own_program_in_family=sum(own_sig$in_family)),
                scores=list(group_means_projection=group_means,secretory_projection_sd=sec_sd,
                            correlation_projection_vs_mean_signed_shift=cor(donors$projection,donors$mean_signed_shift)),
                validity_gate=c(gate,list(rules=gate_cfg)),
                M1_primary=m1,M2_secretory=m2,M2_subphase=m2s,sensitivities=sens,leave_one_donor_out=loo_summary,gene_level=gene_level,
                label=label,label_rules=lab$per_cohort,figure=fig,started_utc=started,finished_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
mde_write_json(summary,file.path(out,"summary.json"))
outputs <- c("identity_check.json","donor_scores.tsv","signature_genes.tsv","leave_one_donor_out.tsv","interaction_de_all.tsv","interaction_top20_by_abs_t.tsv","summary.json",if(fig=="scores.png") "scores.png")
prov <- list(analysis="SECRETORY_PROGRAM_DISEASE_CONTRAST",cohort=cohort,synthetic=synthetic,
             freeze=if(synthetic) NULL else list(path=opts$freeze,sha256=freeze_sha,artifact_check=freeze_check),
             config_sha256=list(statistical=mde_sha(opts$config),phase=mde_sha(opts[["phase-config"]]),secretory_program=mde_sha(opts[["sp-config"]])),
             code_sha256=list(lib=mde_sha("scripts/secretory_program_lib.R"),analysis=mde_sha("scripts/secretory_program_analysis.R"),mde_lib=mde_sha("scripts/mde_lib.R")),
             inputs=inp,input_sha256=input_hashes,identity_check=identity,software=mde_software(),
             started_utc=started,finished_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
             output_sha256=setNames(lapply(outputs,function(f) mde_sha(file.path(out,f))),outputs))
mde_write_json(prov,file.path(out,"provenance.json"))
cat("SECRETORY_PROGRAM_DONE",cohort,"label",label,"theta",if(identical(m1$status,"ESTIMATED")) sprintf("%.4f [%.4f, %.4f] p=%.3g",m1$estimate,m1$ci_low,m1$ci_high,m1$p) else m1$status,
    "gate",gate$pass,"auc",sprintf("%.3f",auc),"signature",sum(sig$in_family),"\n")
