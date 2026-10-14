# NF1 exon 23a inclusion (isoform II fraction) in public PBMC RNA-seq via Snaptron
# Run locally in R (needs internet). Install once:
#   install.packages("BiocManager")
#   BiocManager::install(c("snapcount", "SummarizedExperiment"))
#
# In RStudio: open this file and click "Source" (Ctrl+Shift+S), not "Run".
# If `jx` already exists in your session, the (slow) Snaptron query is skipped.

needed  <- c("snapcount", "SummarizedExperiment")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these first, then re-run: BiocManager::install(c(\"",
                          paste(missing, collapse = "\", \""), "\"))")
suppressPackageStartupMessages({
  library(snapcount)
  library(SummarizedExperiment)
})

EXON_LEN  <- 63   # NF1 exon 23a is 63 bp
MIN_READS <- 20   # minimum informative junction reads per sample

# 1. All junctions within NF1 (skipped if jx already in memory) -----------------------
if (!exists("jx")) {
  comps <- names(Compilation)
  COMP_NAME <- if ("srav3h" %in% comps) "srav3h" else if ("srav2" %in% comps) "srav2" else
    stop("Neither srav3h nor srav2 available - update snapcount")
  message("Using compilation: ", COMP_NAME)
  qb <- QueryBuilder(compilation = Compilation[[COMP_NAME]], regions = "NF1")
  jx <- query_jx(qb)
}
message("Junctions: ", nrow(jx), "  Samples: ", ncol(jx))

# 2. Find the 63-bp cassette exon -------------------------------------------------------
# Snaptron coords: start = first intronic base, end = last intronic base
# skip (s,e) + left inclusion (s, a-1) + right inclusion (b+1, e), exon = a..b
gr <- rowRanges(jx)
j  <- data.frame(idx = seq_along(gr), start = start(gr), end = end(gr))
cov <- if ("coverage_sum" %in% names(mcols(gr))) mcols(gr)$coverage_sum else rep(0, length(gr))

L <- merge(j, j, by = "start", suffixes = c(".skip", ".left"))
L <- L[L$end.left < L$end.skip, ]
R <- merge(j, j, by = "end", suffixes = c(".skip", ".right"))
R <- R[R$start.right > R$start.skip, ]
LR <- merge(L, R, by = "idx.skip")
LR <- LR[LR$start.right - LR$end.left - 1 == EXON_LEN, ]
if (!nrow(LR)) stop("No 63-bp cassette exon found among NF1 junctions")

LR$exon_start <- LR$end.left + 1
LR$exon_end   <- LR$start.right - 1
LR$cov_left  <- cov[LR$idx.left]
LR$cov_right <- cov[LR$idx.right]
LR$cov_skip  <- cov[LR$idx.skip]
# A true alternative exon has ALL THREE junctions well used. Rank by the weakest one,
# so artefacts that borrow one busy constitutive junction don't win.
LR$support   <- pmin(LR$cov_left, LR$cov_right, LR$cov_skip)
LR <- LR[order(-LR$support), ]
print(head(LR[, c("exon_start", "exon_end", "cov_left", "cov_right", "cov_skip", "support")], 5))
cx <- LR[1, ]   # expect exon 23a at chr17:31,252,938-31,253,000 (hg38); check against Ensembl
message("Selected cassette exon: chr17:", cx$exon_start, "-", cx$exon_end)

# 3. Counts for just those three junctions (keeps memory small) ------------------------
cnt   <- as.matrix(assay(jx, 1)[c(cx$idx.left, cx$idx.right, cx$idx.skip), , drop = FALSE])
left  <- cnt[1, ]; right <- cnt[2, ]; skip <- cnt[3, ]
incl  <- (left + right) / 2
res <- data.frame(rail_id = colnames(jx), left = left, right = right, skip = skip,
                  informative = incl + skip, isoformII_frac = incl / (incl + skip))

# 4. Keep PBMC samples using metadata text ----------------------------------------
md <- as.data.frame(colData(jx))
chr_cols <- vapply(md, function(x) is.character(x) || is.factor(x), logical(1))
txt <- do.call(paste, c(lapply(md[chr_cols], function(x) ifelse(is.na(x), "", as.character(x))),
                        sep = " | "))

is_pbmc    <- grepl("PBMC|peripheral blood mononuclear", txt, ignore.case = TRUE)
is_exclude <- grepl("leukemia|lymphoma|cancer|tumou?r|HIV|COVID|SARS|sepsis|lupus|SLE|arthritis|stimulat|LPS|infect|patient",
                    txt, ignore.case = TRUE)
res$study    <- if ("study_accession" %in% names(md)) md$study_accession else NA
res$metadata <- txt

pbmc <- res[is_pbmc & !is_exclude & res$informative >= MIN_READS, ]
message("PBMC-labelled samples: ", sum(is_pbmc),
        " | after exclusions + read filter: ", nrow(pbmc))

if (nrow(pbmc)) {
  print(summary(pbmc$isoformII_frac))
  message("Median isoform II fraction: ", round(median(pbmc$isoformII_frac), 3))
}

write.csv(pbmc, "nf1_exon23a_pbmc_psi.csv", row.names = FALSE)
write.csv(res[res$informative >= MIN_READS, ], "nf1_exon23a_all_samples_psi.csv", row.names = FALSE)
message("Wrote nf1_exon23a_pbmc_psi.csv and nf1_exon23a_all_samples_psi.csv in ", getwd())
