# NF1 exon 23a splicing in public RNA-seq data

R scripts that estimate the fraction of NF1 isoform II (exon 23a included) in public RNA-seq data. They provide an independent benchmark for ddPCR measurements of NF1 isoform I and II fractional abundance in peripheral blood mononuclear cells (PBMCs).

The analysis was written for a Master of Diagnostic Genomics thesis (Queensland University of Technology).

## Background

NF1 exon 23a is a 63 bp alternative exon. It adds 21 amino acids to the GAP-related domain of neurofibromin. Transcripts that skip it encode isoform I. Transcripts that include it encode isoform II.

Published estimates of the isoform balance in blood differ. Assunto et al. (2019) reported isoform II as the predominant form in blood. ddPCR in PBMCs for this thesis gave roughly 75% isoform I and 25% isoform II. These scripts use public RNA-seq junction counts as an independent check.

## Method

All three scripts follow the same approach.

1. Download all exon-exon junction counts within NF1.
2. Find every 63 bp cassette exon, defined by two inclusion junctions and one skipping junction.
3. Select the candidate whose weakest junction has the most reads. A true alternative exon uses all three junctions. Ranking by total reads instead can select artefacts that borrow one busy constitutive junction.
4. Check the selected exon against the expected exon 23a coordinates, chr17:31,252,938–31,253,000 (GRCh38).
5. Calculate the isoform II fraction for each sample as `mean(left inclusion, right inclusion) / (mean(left inclusion, right inclusion) + skipping)`.
6. Report per-sample estimates for samples with at least 20 informative reads. Also report a pooled, read-weighted estimate per group with a Wilson 95% confidence interval.

## Scripts

### `nf1_exon23a_gtex_snaptron.R` (main result)

Uses the Snaptron GTEx compilation. Compares whole blood with EBV-transformed lymphocytes. Tries the newer `gtexv2` compilation first and falls back to the original `gtex` compilation if snapcount fails to build the result.

Outputs `nf1_exon23a_gtex_samples.csv` and `nf1_exon23a_gtex_summary.csv`.

### `nf1_exon23a_pbmc_snaptron.R` (exploratory)

Searches the Snaptron SRA compilation for samples whose metadata mentions PBMCs. This did not find healthy PBMCs. The samples that passed the filter were K562 cells, in vitro polarised T helper cells, CD34+ cells undergoing erythroid differentiation and monocyte-derived macrophages. Their methods sections mentioned PBMC isolation. The script is kept to document this.

### `nf1_exon23a_monaco_recount3.R` (exploratory)

Uses recount3 junction counts for SRA project SRP125125. This project holds the Monaco et al. (2019) dataset of total PBMCs and 29 sorted immune cell types from healthy donors (GEO GSE107011). Coverage at exon 23a was too low to use. All 127 samples together gave only 79 informative reads.

## Results

| GTEx tissue | Samples | Samples with ≥20 reads | Isoform II fraction |
|---|---|---|---|
| Whole blood | 456 | 5 | 0.27 pooled (95% CI 0.25–0.30) |
| EBV-transformed lymphocytes | 139 | 87 | 0.40 median (IQR 0.33–0.46), 0.41 pooled (95% CI 0.39–0.42) |

Isoform I was the predominant NF1 transcript in both tissues. The whole blood estimate is close to the ddPCR result in PBMCs. Neither tissue showed isoform II predominance.

### Limitations

- Most whole blood samples have low coverage at exon 23a, so the whole blood estimate relies on pooling.
- Pooled confidence intervals treat reads as independent. They do not account for variation between donors and are likely to understate uncertainty.
- EBV-transformed lymphocytes are activated B cells. They are not equivalent to resting PBMCs.
- Junction counts come from the Rail-RNA aligner used by Snaptron. They are not directly comparable with ddPCR.

## Requirements

R with these Bioconductor packages:

```r
install.packages("BiocManager")
BiocManager::install(c("snapcount", "SummarizedExperiment"))   # Snaptron scripts
BiocManager::install("recount3")                                # Monaco script only
```

The scripts need internet access to the Snaptron and recount3 servers. Institutional networks may block these.

## Usage

Run each script as a whole file so that it stops at the first error.

```r
source("nf1_exon23a_gtex_snaptron.R", echo = TRUE)
```

Output CSVs are written to the R working directory. The Monaco script downloads a large junction file and may need several GB of memory.

## References

Assunto, A., et al. (2019). [Add full citation.]

GTEx Consortium. (2017). Genetic effects on gene expression across human tissues. *Nature, 550*(7675), 204–213. https://doi.org/10.1038/nature24277

Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., et al. (2019). RNA-seq signatures normalized by mRNA abundance allow absolute deconvolution of human immune cell types. *Cell Reports*. [Add volume, pages and DOI.]

Wilks, C., Gaddipati, P., Nellore, A., & Langmead, B. (2018). Snaptron: Querying splicing patterns across tens of thousands of RNA-seq samples. *Bioinformatics, 34*(1), 114–116. https://doi.org/10.1093/bioinformatics/btx547

## Author

Renee Pascouau
