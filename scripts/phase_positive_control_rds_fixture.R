#!/usr/bin/env Rscript
# Synthetic Seurat-like fixture for the ported RDS stream parser.
# Base R S4 classes with the dgCMatrix slot layout; no Matrix or Seurat dependency.
# Saved with saveRDS defaults (XDR, serialization v3, gzip), like the GEO object.
# Writes fixture.rds and FIXTURE_EXPECTED.json to the output directory. Synthetic only.
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript phase_positive_control_rds_fixture.R <outdir>")
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261026)
ng <- 12000L
nc <- 1100L
per <- sample(50:200, nc, replace = TRUE)
ilist <- vector("list", nc)
xlist <- vector("list", nc)
for (j in seq_len(nc)) {
  ix <- sort(sample.int(ng, per[j]) - 1L)
  v <- as.double(sample(1:30, length(ix), replace = TRUE))
  frac <- seq_along(v) %% 7L == 0L
  v[frac] <- v[frac] + 0.375
  ilist[[j]] <- ix
  xlist[[j]] <- v
}
i <- as.integer(unlist(ilist))
x <- as.double(unlist(xlist))
p <- as.integer(c(0L, cumsum(per)))
genes <- paste0("gene", seq_len(ng))
cells <- paste0("sample", rep(1:4, length.out = nc), "_", sprintf("CELL%05d", seq_len(nc)))
setClass("dgCMatrix", representation(i = "integer", p = "integer", Dim = "integer", Dimnames = "list", x = "numeric", factors = "list"))
setClass("FixtureAssay", representation(counts = "dgCMatrix", data = "matrix"))
setClass("StreamFixture", representation(assays = "list", meta.data = "data.frame"))
counts <- new("dgCMatrix", i = i, p = p, Dim = c(ng, nc), Dimnames = list(genes, cells), x = x, factors = list())
counts2 <- new("dgCMatrix", i = i, p = p, Dim = c(ng, nc), Dimnames = list(genes, cells), x = x * 2, factors = list())
colsums <- vapply(seq_len(nc), function(j) if (p[j + 1] > p[j]) sum(x[(p[j] + 1):p[j + 1]]) else 0, numeric(1))
rowsums <- numeric(ng)
agg <- tapply(x, i, sum)
rowsums[as.integer(names(agg)) + 1L] <- as.numeric(agg)
refined <- sample(c("eF1", "dS2", "Prv-MYH11", NA), nc, replace = TRUE)
meta <- data.frame(
  sample = factor(sub("_.*", "", cells)),
  nCount_RNA = colsums,
  AnnotationMain = factor(ifelse(is.na(refined), "epithelial", "stromal")),
  AnnotationRefined = refined,
  intcol = seq_len(nc),
  flag = rep(c(TRUE, FALSE, NA), length.out = nc),
  row.names = cells,
  stringsAsFactors = FALSE
)
fx <- new("StreamFixture",
  assays = list(
    RNA = new("FixtureAssay", counts = counts, data = matrix(0, ng, nc)),
    SCT = new("FixtureAssay", counts = counts2, data = matrix(1, ng, nc))
  ),
  meta.data = meta
)
saveRDS(fx, file.path(out, "fixture.rds"))
expected <- list(
  counts_dim = c(ng, nc), counts_nnz = length(x), counts_sum = sum(x),
  fractional_entries = sum(x != floor(x)), colSums = unname(colsums), rowSums = unname(rowsums),
  genes_head = head(genes, 3), cells_head = head(cells, 3),
  meta = list(columns = names(meta), rows = nrow(meta), sample_levels = levels(meta$sample),
              sample_codes_head = as.integer(meta$sample)[1:10], intcol_sum = sum(meta$intcol),
              refined_na = sum(is.na(refined)), refined_head = refined[1:10], flag_head = meta$flag[1:10],
              rownames_head = rownames(meta)[1:3]),
  assay_names = c("RNA", "SCT"), dense_bytes_each = ng * nc * 8,
  r_version = R.version.string
)
write_json(expected, file.path(out, "FIXTURE_EXPECTED.json"), auto_unbox = TRUE, digits = NA, na = "null")
cat("fixture written:", file.path(out, "fixture.rds"), "\n")
