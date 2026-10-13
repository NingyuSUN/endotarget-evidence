stopifnot(file.exists("scripts/resample_donor_models.R"))
source("scripts/resample_donor_models.R")
cfg <- jsonlite::fromJSON("config/statistical_analysis.json")
set.seed(117)
n <- 22L; genes <- sprintf("ENSG%011d",seq_len(1500))
clinical <- data.frame(donor_id=paste0("toy",seq_len(n)), disease=c(rep("ENDO",15),rep("non-ENDO",7)), age=sample(20:50,n,TRUE), library_group_Lane=rep(c("A","B"),11), adenomyosis=rep(c(0,0,1),length.out=n), primary_nuclei=100)
counts <- matrix(rnbinom(1500*n,mu=100,size=4),1500,n,dimnames=list(genes,clinical$donor_id))
counts[1500,] <- 0
built <- rebuild_resample_design(clinical,"heca",cfg)
original <- fit_donor_model(counts,built$design,clinical,"heca",cfg)
timing <- system.time(fixed <- fit_fixed_family(counts[original$keep,,drop=FALSE],built$design,cfg))
stopifnot(isTRUE(all.equal(original$table$logFC,fixed$table$logFC,tolerance=1e-12)), isTRUE(all.equal(original$table$t,fixed$table$t,tolerance=1e-12)), isTRUE(all.equal(original$table$adj.P.Val,fixed$table$adj.P.Val,tolerance=1e-12)))
# Repeat individuals are intentional bootstrap weights; no unique-donor eligibility gate.
ix <- c(rep(c(1,2,3),5),rep(c(16,17),length.out=7))
dup <- clinical[ix,]; dup$donor_id <- paste0("draw",seq_len(n))
b <- rebuild_resample_design(dup,"heca",cfg)
stopifnot(b$receipt$unique_original_donors == n) # receipt counts caller IDs; runner records true source IDs
ct <- counts[original$keep,ix,drop=FALSE]; colnames(ct)<-dup$donor_id
f <- fit_fixed_family(ct,b$design,cfg)
stopifnot(nrow(f$table)==sum(original$keep), b$receipt$nominal_residual_df>=5)
# A low count gene retained by original family remains, even if the draw fails expression support.
ct[nrow(ct),] <- 1
f <- fit_fixed_family(ct,b$design,cfg)
stopifnot(tail(f$table$gene_id,1)==tail(rownames(ct),1))
bad <- clinical; bad$age <- as.numeric(bad$disease=="ENDO")
stopifnot(inherits(try(rebuild_resample_design(bad,"heca",cfg),silent=TRUE),"try-error"))
bad <- clinical; bad$library_group_Lane <- bad$disease
stopifnot(inherits(try(rebuild_resample_design(bad,"heca",cfg),silent=TRUE),"try-error"))
constant <- clinical; constant$adenomyosis <- 0
stopifnot("adenomyosis" %in% rebuild_resample_design(constant,"heca",cfg)$receipt$dropped_constant_nuisance)
ledger <- generate_draws(clinical,cfg)
stopifnot(nrow(ledger)==500*22, all(table(ledger$replicate)==22), identical(ledger,generate_draws(clinical,cfg)))
mask <- do.call(rbind,lapply(c("B0","B1","G1"),function(m)do.call(rbind,lapply(c(20,50),function(k)data.frame(method=m,k=k,rank=seq_len(k),gene_id=genes[seq_len(k)],discovery_sign=1,original_evaluable=TRUE)))))
validate_mask(mask)
valid <- fixed$table; valid$logFC<-1; valid$fit_estimable<-TRUE
paired <- evaluate_fixed_masks(valid,mask)
stopifnot(paired$valid, paired$difference==0)
valid$logFC[1] <- 0
stopifnot(!evaluate_fixed_masks(valid,mask)$valid)
valid$logFC[1]<-1; mask$original_evaluable[mask$gene_id==genes[1]]<-FALSE
valid$logFC[1]<-0
stopifnot(evaluate_fixed_masks(valid,mask)$valid)
# Outside primary top20 but inside secondary top50 also invalidates joint replicate.
mask$original_evaluable[]<-TRUE; valid$logFC[40]<-0
stopifnot(!evaluate_fixed_masks(valid,mask)$valid)
stopifnot(bootstrap_summary(c(rep(NA_real_,51),rep(0,449)),500)$interval_suppressed, !bootstrap_summary(c(rep(NA_real_,50),rep(0,450)),500)$interval_suppressed)
stopifnot(identical(bootstrap_summary(seq_len(500),500)$interval,as.numeric(quantile(seq_len(500),c(.025,.975),type=7))))
# A syntactically valid checkpoint edit must never change a resumed paired result.
stopifnot(exists("read_bootstrap_checkpoint"),exists("write_bootstrap_checkpoint"))
checkpoint_dir<-tempfile("checkpoint_",tmpdir="work/synthetic");dir.create(checkpoint_dir,recursive=TRUE)
checkpoint_path<-file.path(checkpoint_dir,"replicate_001.json")
valid$logFC[]<-1;result<-evaluate_fixed_masks(valid,mask);result$replicate<-1L;result$binding_sha256<-"fixture-binding";result$family_size<-nrow(valid)
write_bootstrap_checkpoint(result,checkpoint_path)
stopifnot(read_bootstrap_checkpoint(checkpoint_path,"fixture-binding",1L,mask,nrow(valid))$difference==0)
lines<-readLines(checkpoint_path);writeLines(c(lines," "),checkpoint_path)
stopifnot(inherits(try(read_bootstrap_checkpoint(checkpoint_path,"fixture-binding",1L,mask,nrow(valid)),silent=TRUE),"try-error"))
# Even a freshly hashed bad difference must fail fixed-rate consistency.
result$difference<-.2;write_bootstrap_checkpoint(result,checkpoint_path)
stopifnot(inherits(try(read_bootstrap_checkpoint(checkpoint_path,"fixture-binding",1L,mask,nrow(valid)),silent=TRUE),"try-error"))
result$difference<-0;result$rates$B0_top20$evaluable<-19L;write_bootstrap_checkpoint(result,checkpoint_path)
stopifnot(inherits(try(read_bootstrap_checkpoint(checkpoint_path,"fixture-binding",1L,mask,nrow(valid)),silent=TRUE),"try-error"))
unlink(checkpoint_dir,recursive=TRUE)
cat("RESAMPLING_TESTS_PASS", "toy_fit_elapsed",timing[["elapsed"]],"\n")
