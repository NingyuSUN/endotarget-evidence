stopifnot(file.exists("scripts/run_clinical_sensitivities.R"))
source("scripts/run_clinical_sensitivities.R")
cfg<-jsonlite::fromJSON("config/statistical_analysis.json")
set.seed(810)
cl<-data.frame(donor_id=paste0("toy",1:22),disease=c(rep("ENDO",12),rep("non-ENDO",10)),age=sample(20:50,22,TRUE),capture_batch=rep(c("A","B"),11),adenomyosis=c(0,rep(1,4),rep(0,17)),strict_histology_sensitivity=!(1:22 %in% c(1,13)),adenomyosis_free_sensitivity=!(1:22 %in% 2:5),primary_lineage_cells=100)
ct<-matrix(rnbinom(1500*22,mu=100,size=4),1500,22,dimnames=list(sprintf("ENSG%011d",1:1500),cl$donor_id))
ids<-rev(cl$donor_id[cl$strict_histology_sensitivity])
design<-cbind(intercept=1,ENDO_minus_non_ENDO=as.numeric(cl$disease[match(ids,cl$donor_id)]=="ENDO"));rownames(design)<-ids
x<-select_sensitivity_inputs(ct,cl,design,"discovery","strict_histology_conflict_exclusion",cfg)
stopifnot(identical(colnames(x$counts),ids),identical(x$clinical$donor_id,ids),identical(rownames(x$design),ids))
# Frozen design may not omit an eligible donor or include a conflict.
stopifnot(inherits(try(select_sensitivity_inputs(ct,cl,design[-1,,drop=FALSE],"discovery","strict_histology_conflict_exclusion",cfg),silent=TRUE),"try-error"))
bad<-ct[,22:1];stopifnot(inherits(try(select_sensitivity_inputs(bad,cl,design,"discovery","strict_histology_conflict_exclusion",cfg),silent=TRUE),"try-error"))
ids<-cl$donor_id[cl$strict_histology_sensitivity & cl$adenomyosis_free_sensitivity]
design<-cbind(intercept=1,ENDO_minus_non_ENDO=as.numeric(cl$disease[match(ids,cl$donor_id)]=="ENDO"));rownames(design)<-ids
x<-select_sensitivity_inputs(ct,cl,design,"discovery","strict_and_adenomyosis_free",cfg)
a<-sensitivity_eligibility(x,cfg)
stopifnot(!a$estimable,a$ENDO==7,a$non_ENDO==9,grepl("minimum",a$reason))
# Synthetic fit sanity uses exactly the frozen estimator after a true subset.
ids<-cl$donor_id[cl$adenomyosis_free_sensitivity]
design<-cbind(intercept=1,ENDO_minus_non_ENDO=as.numeric(cl$disease[match(ids,cl$donor_id)]=="ENDO"));rownames(design)<-ids
x<-select_sensitivity_inputs(ct,cl,design,"discovery","adenomyosis_free",cfg)
stopifnot(sensitivity_eligibility(x,cfg)$estimable)
m<-fit_donor_model(x$counts,x$design,x$clinical,"discovery",cfg)
stopifnot(m$audit$ENDO==8,m$audit$non_ENDO==10,nrow(m$table)>1000,all(is.finite(m$table$logFC)))
cat("CLINICAL_SENSITIVITY_TESTS_PASS\n")
