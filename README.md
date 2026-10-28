# OSCC_cytokine_bulkRNAseq
Analysis code for OSCC in bulk transcriptomic analyses
## Workflow

| File | Analysis |
|---|---|
| `01_bulk_RNAseq_preprocessing.sh` | FastQC, MultiQC, fastp, HISAT2, SAMtools, StringTie and prepDE |
| `02_cytokine_bulk_RNAseq.Rmd` | IL-10/TNF-alpha differential expression, batch diagnostics, enrichment and plots |
| `03_TPM_for_CIBERSORTx.Rmd` | Gene-symbol TPM matrices for immune deconvolution |

This repository covers bulk RNA-seq analyses. Clinical survival/machine learning,
single-cell analyses and high-resolution cell-type-specific expression analyses
are outside the scope of these files.

## Inputs

- `name.txt`: the 103 study sample IDs.
- `raw/`: one R1 and one R2 compressed FASTQ per sample, matching
  `<sample_id>_*_R1_001.fastq.gz` and `<sample_id>_*_R2_001.fastq.gz`.
  Concatenate corresponding reads for technical runs of the same library
  before preprocessing. Keep only the resulting read pair in `raw/`.
- `ref/`: the HISAT2 index, `Homo_sapiens.GRCh38.102.gtf`, and `prepDE.py3`.
  Reference paths are configured at the top of script 01.
- `gene_counts_103_samples.csv`: gene IDs in the first column, formatted as
  `Ensembl_ID|gene_symbol`, and raw-count columns named by sample ID.
- `cytokineRNAseq.csv`: combined metadata with `sample_id`, `batch`, `PS`,
  `IL10_roc`, `TNFa_roc`; `ages` is optional for the age diagnostic.
- StringTie `<sample_id>_gene_abund.tab` tables for notebook 03, with gene IDs
  and TPM values. These tables are separate inputs: script 01's StringTie
  command does not request a gene-abundance report with `-A`.

`metadata_columns_template.csv` documents the metadata structure and contains
an artificial example, not a patient record. Sequencing files, expression
matrices, reference files and private clinical tables are not bundled.
Data accessions and access conditions should accompany the manuscript.

## Software

Linux: FastQC, MultiQC, fastp, HISAT2, SAMtools, StringTie, Python 3 and prepDE.
R: rmarkdown/knitr, dplyr/tidyverse, DESeq2, ComplexHeatmap, org.Hs.eg.db,
AnnotationDbi, EnrichmentBrowser, circlize, clusterProfiler, enrichplot, vegan,
ReactomePA, fgsea, msigdbr, limma, ggrepel and grid.
The gene-set code uses the `msigdbr(category=...)` / `entrez_gene` interface.
Use compatible versions; record the software environment used for analysis.

## Execution

From the project root, after preparing inputs and configuring reference paths:

```bash
bash 01_bulk_RNAseq_preprocessing.sh
```

Script 01 writes `mapped/GC.csv` and `mapped/TC.csv`. Notebook 02 reads the
study gene-count matrix as `gene_counts_103_samples.csv`. Gene identifiers and
sample columns must agree when regenerating that input from prepDE.

Run the R notebooks in RStudio with Knit, or:

```r
rmarkdown::render("02_cytokine_bulk_RNAseq.Rmd")
rmarkdown::render("03_TPM_for_CIBERSORTx.Rmd")
```

## Analysis definitions

Genes require at least 10 counts in `round(0.9*n_samples)` samples.
Differential expression uses DESeq2, cytokine reference code 1, and adjusted
P < 0.05. `PS` is a categorical pathologic-stage covariate; `batch2` distinguishes
batch 2 from other batches. Biological high/low labels and concentration
cutoffs belong in the metadata codebook.

TNF-alpha ORA uses an exploratory gene-selection cutoff of adjusted P < 0.6.
GO GSEA uses `pAdjustMethod="none"`; this output is not BH-adjusted GSEA.
Enrichment backgrounds, rankings and collections are specified in notebook 02.
Volcano plotting uses fixed axis limits. Relative fractions and absolute
scores are distinct quantities; only fractions represent proportions.

Record the deconvolution service/version, signature matrix, input scale,
permutations, relative/absolute mode, batch correction and normalization
settings with the exported deconvolution results.

## Citation and reuse

Use the manuscript citation and archived code DOI associated with the version
being analyzed. Consult the repository license for reuse terms; third-party
tools and signature matrices retain their own licenses.
