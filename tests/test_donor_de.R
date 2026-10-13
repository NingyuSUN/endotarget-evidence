# Synthetic boundary and independent weighted-regression checks; not scientific validation.
source("scripts/donor_de_lib.R")
suppressPackageStartupMessages(library(jsonlite))
cfg <- fromJSON("config/statistical_analysis.json", simplifyVector=TRUE)
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(paste("FAIL",label))
  checks <<- checks+1L
  cat("PASS",label,"\n")
}
reject <- function(expr,label) {
  failed <- inherits(try(force(expr),silent=TRUE),"try-error")
  check(failed,label)
}
set.seed(42)
ids <- paste0("toy",seq_len(20))
genes <- sprintf("ENSG%011d",seq_len(300))
groups <- rep(c(0,1),each=10)
cts <- matrix(rnbinom(300*20,mu=80,size=10),nrow=300,
              dimnames=list(genes,ids))
cts[1:30,groups==1] <- cts[1:30,groups==1]*4
design <- cbind(intercept=1,ENDO_minus_non_ENDO=groups,
                age_centered_years=rep(seq(-4.5,4.5),2))
rownames(design) <- ids
clinical <- data.frame(donor_id=ids,disease=ifelse(groups==1,"ENDO","non-ENDO"),
                       primary_lineage_cells=200,primary_nuclei=200)
check(validate_donor_inputs(cts,design,clinical,"discovery",cfg)$residual_df==17,
      "valid independent donor design")
bad <- cts[,rev(ids)]
reject(validate_donor_inputs(bad,design,clinical,"discovery",cfg),"donor swap")
bad <- cts; bad[1,1] <- -1
reject(validate_donor_inputs(bad,design,clinical,"discovery",cfg),"negative counts")
bad[1,1] <- Inf
reject(validate_donor_inputs(bad,design,clinical,"discovery",cfg),"nonfinite counts")
bad <- cts; bad[,1] <- 0
reject(validate_donor_inputs(bad,design,clinical,"discovery",cfg),"zero library")
bad <- design; bad <- cbind(bad,confounded_disease=groups)
reject(validate_donor_inputs(cts,bad,clinical,"discovery",cfg),"confounded full design")
bad <- clinical; bad$donor_id[2] <- bad$donor_id[1]
reject(validate_donor_inputs(cts,design,bad,"discovery",cfg),"duplicate person")
bad <- cts; bad[1,1] <- bad[1,1]+0.25
check(validate_donor_inputs(bad,design,clinical,"discovery",cfg)$n_donors==20,
      "fractional discovery estimates accepted unchanged")
reject(validate_donor_inputs(bad,design,clinical,"heca",cfg),"fractional HECA rejected")
bad <- clinical; bad$primary_lineage_cells[1] <- 49
reject(validate_donor_inputs(cts,design,bad,"discovery",cfg),"below primary cell threshold")
fit <- fit_donor_model(cts,design,clinical,"discovery",cfg)
check(fit$filter_audit$minimum_support_donors==10,"fixed group support count")
check(is.character(fit$warnings),"warnings are retained explicitly")
check(nrow(fit$table)==sum(fit$keep),"all filtered genes reported")
check(median(fit$table$logFC[fit$table$gene_id %in% genes[1:30]])>1,
      "known case increase has positive coefficient")
g <- 40
manual <- lm.wfit(design,fit$voom$E[g,],fit$voom$weights[g,])
check(abs(manual$coefficients[2]-fit$table$logFC[g])<1e-10,
      "coefficient matches independent weighted regression")
check(max(abs(p.adjust(fit$table$P.Value,"BH")-fit$table$adj.P.Val))<1e-12,
      "BH family is all externally or internally tested genes")
check(all(fit$table$CI_low<=fit$table$logFC & fit$table$CI_high>=fit$table$logFC),
      "moderated coefficient intervals contain estimates")
check(max(abs(log(prod(fit$normalization$norm_factor))))<1e-10,
      "TMM multiplicative factors centered")
cts2 <- rbind(cts,ENSG99999999999=rep(c(10,0),c(1,19)))
fit2 <- fit_donor_model(cts2,design,clinical,"discovery",cfg)
check(!tail(fit2$keep,1),"one-person low-support gene removed")
flipped_clinical <- clinical
flipped_clinical$disease <- ifelse(clinical$disease=="ENDO","non-ENDO","ENDO")
flipped_design <- design
flipped_design[,"ENDO_minus_non_ENDO"] <- 1-design[,"ENDO_minus_non_ENDO"]
flipped <- fit_donor_model(cts,flipped_design,flipped_clinical,"discovery",cfg)
check(identical(fit$keep,flipped$keep),"filter unchanged by label permutation with fixed group sizes")
check(max(abs(fit$table$logFC+flipped$table$logFC))<1e-10,
      "reversing disease coding reverses all fitted coefficients")
reject(validate_voom_weights(matrix(c(1,0),1)),"zero observation weight fails")
reject(validate_voom_weights(matrix(c(1,Inf),1)),"nonfinite observation weight fails")
check(all(fit$table$df_residual>=cfg$minimum_residual_df),"gene residual support retained")
cat("RESULT",checks,"checks passed\n")
