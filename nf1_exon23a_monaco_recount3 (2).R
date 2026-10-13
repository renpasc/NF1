# NF1 exon 23a inclusion (isoform II fraction) in Monaco et al. 2019 healthy PBMCs
# and sorted immune subsets (GEO GSE107011), using recount3 junction counts.
#
# Install once:
#   install.packages("BiocManager")
#   BiocManager::install(c("recount3", "SummarizedExperiment"))
#
# Run in the R console:   source(file.choose(), echo = TRUE)
# NOTE: downloads the project's whole junction file (can be a few hundred MB) and
# may need several GB of RAM while loading. Run on a good connection.

needed  <- c("recount3", "SummarizedExperiment")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install first: BiocManager::install(c(\"",
                          paste(missing, collapse = "\", \""), "\"))")
suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
})

PROJECT     <- "SRP125125"                 # SRA project expected to hold GSE107011
EXPECT_EXON <- c(31252938, 31253000)       # exon 23a hg38, found in the Snaptron run
EXON_LEN    <- 63
MIN_READS   <- 20
NF1_REGION  <- GenomicRanges::GRanges("chr17:31090000-31390000")   # generous NF1 window

# 1. Locate the project and download junction counts ---------------------------------
hp <- available_projects()
pinfo <- hp[hp$project == PROJECT & hp$project_type == "data_sources", ]
if (!nrow(pinfo)) stop(PROJECT, " not found in recount3")
print(pinfo)

if (!exists("rse_jx")) rse_jx <- create_rse(pinfo, type = "jxn")
message("Junctions: ", nrow(rse_jx), "  Samples: ", ncol(rse_jx))

# Confirm this is Monaco: titles should name immune cell types / PBMC
md <- as.data.frame(colData(rse_jx))
title_col <- intersect(c("sra.sample_title", "sra.experiment_title"), names(md))[1]
message("\nFirst sample titles (check these are Monaco PBMC / sorted subsets):")
print(head(md[[title_col]], 15))

# 2. Keep NF1 junctions and find exon 23a ---------------------------------------------
nf1 <- rse_jx[IRanges::overlapsAny(rowRanges(rse_jx), NF1_REGION, type = "within"), ]
gr  <- rowRanges(nf1)
j   <- data.frame(idx = seq_along(gr), start = start(gr), end = end(gr))
tot <- Matrix::rowSums(assay(nf1, 1))

L  <- merge(j, j, by = "start", suffixes = c(".skip", ".left"));  L <- L[L$end.left < L$end.skip, ]
R  <- merge(j, j, by = "end",   suffixes = c(".skip", ".right")); R <- R[R$start.right > R$start.skip, ]
LR <- merge(L, R, by = "idx.skip")
LR <- LR[LR$start.right - LR$end.left - 1 == EXON_LEN, ]
if (!nrow(LR)) stop("No 63-bp cassette exon found in NF1")
LR$exon_start <- LR$end.left + 1
LR$exon_end   <- LR$start.right - 1
LR$cov_left   <- tot[LR$idx.left]; LR$cov_right <- tot[LR$idx.right]; LR$cov_skip <- tot[LR$idx.skip]
LR$support    <- pmin(LR$cov_left, LR$cov_right, LR$cov_skip)
LR <- LR[order(-LR$support), ]
print(head(LR[, c("exon_start", "exon_end", "cov_left", "cov_right", "cov_skip")], 5))
cx <- LR[1, ]
message("Selected cassette exon: chr17:", cx$exon_start, "-", cx$exon_end)
if (cx$exon_start != EXPECT_EXON[1] || cx$exon_end != EXPECT_EXON[2])
  warning("Selected exon differs from the expected exon 23a coordinates - check before using")

# 3. Per-sample isoform II fraction ----------------------------------------------------
cnt  <- as.matrix(assay(nf1, 1)[c(cx$idx.left, cx$idx.right, cx$idx.skip), , drop = FALSE])
incl <- (cnt[1, ] + cnt[2, ]) / 2
skip <- cnt[3, ]
res  <- data.frame(run = colnames(nf1), title = md[[title_col]],
                   left = cnt[1, ], right = cnt[2, ], skip = skip,
                   informative = incl + skip, isoformII_frac = incl / (incl + skip))

# Cell-type label = sample title with donor/replicate suffix stripped (eyeball and adjust)
res$cell_type <- trimws(sub("[_ -]*(S|donor|rep)?[_ -]*[0-9A-Z]{0,6}$", "", res$title, ignore.case = TRUE))

ok <- res[res$informative >= MIN_READS, ]
message("\nSamples with >= ", MIN_READS, " informative reads: ", nrow(ok), " of ", nrow(res))
message("Informative reads per sample:")
print(summary(res$informative))
message("Total reads across all samples - left: ", sum(res$left), "  right: ", sum(res$right),
        "  skip: ", sum(res$skip))

# 4. Pooled (read-weighted) estimate per cell type with Wilson 95% CI -------------------
wilson <- function(x, n) {
  z <- 1.96; p <- x / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))
}
# Pool ALL samples of a cell type (per-sample depth is low; pooling across donors is the point).
by_type <- do.call(rbind, lapply(split(res, res$cell_type), function(g) {
  i <- sum((g$left + g$right) / 2); s <- sum(g$skip); n <- i + s
  ci <- if (n > 0) wilson(i, n) else c(NA, NA)
  deep <- g$isoformII_frac[g$informative >= MIN_READS]
  data.frame(cell_type = g$cell_type[1], n_samples = nrow(g), total_reads = n,
             median_sample = if (length(deep)) round(median(deep), 3) else NA,
             pooled = if (n > 0) round(i / n, 3) else NA,
             ci_low = round(ci[1], 3), ci_high = round(ci[2], 3))
}))
if (is.null(by_type) || !nrow(by_type)) stop("No cell-type groups formed - check the title column")
by_type <- by_type[order(by_type$pooled), ]
print(by_type, row.names = FALSE)

write.csv(res,     "nf1_exon23a_monaco_samples.csv",   row.names = FALSE)
write.csv(by_type, "nf1_exon23a_monaco_celltypes.csv", row.names = FALSE)
message("Wrote nf1_exon23a_monaco_samples.csv and nf1_exon23a_monaco_celltypes.csv in ", getwd())
