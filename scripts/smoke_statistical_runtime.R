#!/usr/bin/env Rscript
suppressPackageStartupMessages({
 library(edgeR); library(limma); library(statmod); library(jsonlite)
})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
set.seed(42)
# Synthetic runtime exercise only; these are not disease/cohort data.
counts <- matrix(rnbinom(200L*12L, mu=100, size=10), nrow=200L)
rownames(counts) <- sprintf("TOY_GENE_%03d", seq_len(nrow(counts)))
group <- factor(rep(c("TOY_A","TOY_B"), each=6L))
design <- model.matrix(~group)
y <- DGEList(counts=counts)
y <- calcNormFactors(y, method="TMM")
v <- voom(y, design=design, span=0.5, plot=FALSE)
fit <- eBayes(lmFit(v, design), robust=TRUE)
stopifnot(all(is.finite(v$E)), all(is.finite(v$weights)),
          all(v$weights>0), all(is.finite(fit$t)),
          qr(design)$rank==2L, all(fit$df.residual==10L))
versions <- setNames(lapply(c("edgeR","limma","statmod","jsonlite"),
 function(x) as.character(packageVersion(x))), c("edgeR","limma","statmod","jsonlite"))
write_json(list(status="PASS", purpose="SYNTHETIC_RUNTIME_SMOKE_ONLY",
 biological_validation=FALSE, actual_cohort_effects_computed=FALSE,
 R=R.version.string, packages=versions, toy_genes=200L, toy_samples=12L,
 design_rank=2L, residual_df=10L, seed=42L, TMM="TMM", voom_span=0.5,
 ebayes_robust=TRUE, session_info=capture.output(sessionInfo())),
 args[[1L]], auto_unbox=TRUE, pretty=TRUE)
