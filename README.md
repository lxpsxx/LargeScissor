# LargeScissor: Optimized Scissor for Large-Scale Single-Cell Data

### About LargeScissor
The original `Scissor` implementation can become memory-intensive on large-scale single-cell datasets, particularly in the network-regularized workflow, because several steps materialize large dense objects; as a result, Out-of-Memory (OOM) errors may occur when cell numbers become very large.

**LargeScissor** Spent half of the Chinese New Year holiday on LargeScissor, I didn't do anything else, basically just three things: First Seurat V5 compatibility, Second preventing dense gene expression matrices, Third keeping the network matrix strictly sparse. I'm quite ashamed, just did a little bit of tiny work. Thank you everyone!

### Acknowledgements
Sincere thanks to Sun and others, the original authors of Scissor.

LargeScissor is now maintained with AI-assisted updates to preserve compatibility with the evolving upstream Scissor software environment, including Seurat and related tools.

## LargeScissor v1.2

LargeScissor v1.2 is a Scissor-compatible fork focused on three practical goals:

- making large-scale runs more memory-aware
- improving compatibility with modern Seurat objects
- fixing several result-affecting issues in the original workflow

The package keeps the original phenotype-guided cell selection interface while making the implementation safer for current single-cell analysis pipelines.

## Quick Start

Install the current GitHub release with:

```r
remotes::install_github("lxpsxx/LargeScissor@v1.2.0")
```

```r
library(LargeScissor)
library(Seurat)

# sc_dataset: a preprocessed Seurat object with an RNA_snn graph
# bulk_dataset: gene x sample bulk expression matrix
# phenotype: numeric vector, binary phenotype, or survival matrix

sc_dataset <- NormalizeData(sc_dataset, verbose = FALSE)
sc_dataset <- FindVariableFeatures(sc_dataset, verbose = FALSE)
sc_dataset <- ScaleData(sc_dataset, verbose = FALSE)
sc_dataset <- RunPCA(sc_dataset, features = VariableFeatures(sc_dataset), verbose = FALSE)
sc_dataset <- FindNeighbors(sc_dataset, dims = 1:10, verbose = FALSE)

fit <- Scissor(
  bulk_dataset = bulk_dataset,
  sc_dataset = sc_dataset,
  phenotype = phenotype,
  family = "cox",
  alpha = 0.05,
  Save_file = "Scissor_inputs.RData",
  Mthread = TRUE,
  Mcore = 8
)

str(fit$Scissor_pos)
str(fit$Scissor_neg)
head(fit$Coefs)
```

## Highlights

### Large-scale and Seurat v5 support

- Added compatibility with Seurat v5 / Assay5 by reading RNA expression from the `data` layer when appropriate
- Preserved the cell-cell graph as a sparse matrix instead of forcing an early dense conversion
- Reduced avoidable dense operations in the preprocessing path where possible

### Method-consistent fixes

- Corrected network penalty alignment for `gaussian` and `cox`
  - In the original implementation, the R wrapper padded the network matrix with an extra intercept-like dimension even though the solver optimized only the original cell coefficients
  - LargeScissor v1.2 keeps the network penalty on the same cell dimension as the optimized coefficient vector
  - `binomial` still retains one explicit intercept dimension because the logistic solver uses an augmented design matrix

- Corrected binomial phenotype preprocessing
  - Two-level factor input is now converted to `0/1` instead of `1/2`
  - Logical and binary numeric inputs are also normalized to a stable binary response

- Added proper continuous phenotype handling for `gaussian`
  - Continuous phenotypes are now treated as numeric responses rather than grouped categories
  - `tag` is optional for `gaussian` and is used only for summary reporting when appropriate

### API and workflow improvements

- Normalized `family` handling in both `Scissor()` and `reliability.test()`
  - This removes the invalid default-vector behavior inherited from the original package

- Added checkpoint safety for `Load_file`
  - Saved preprocessing inputs now record the original regression family
  - Reloading the same checkpoint with a mismatched `family` now fails early with a clear message

- Added optional multiprocessing
  - `Scissor()` and `reliability.test()` now accept `Mthread` and `Mcore`
  - Cross-validation folds and permutation loops can run through `parallel::mclapply()` on Unix-like systems
  - Unsupported platforms or disabled multiprocessing fall back to serial execution

- Preserved cell names in returned coefficients
  - `Coefs` now keep the column identity from `X`, which makes downstream comparison and visualization more reliable

## Validation

LargeScissor v1.2 was validated with focused installation checks and smoke tests.

- `R CMD INSTALL` completed successfully
- `family` normalization in `Scissor()` and `reliability.test()` was verified
- factor-to-`0/1` conversion for `binomial` phenotypes was verified
- continuous phenotype handling for `gaussian` was verified
- `Load_file` family consistency checking was verified
- serial and parallel smoke tests matched for:
  - `APML1` gaussian
  - `APML1` binomial
  - `APML1` cox
  - `reliability.test()` gaussian

## Known Limitations

LargeScissor v1.2 improves the original package substantially, but it is still an incremental fork rather than a full redesign of the Scissor pipeline.

- Quantile normalization still requires a dense expression block and can remain memory-intensive on very large inputs
- Multiprocessing is process-based rather than thread-based, so large `X` and `network` objects can still increase memory pressure
- Alpha-grid parallelization has not been added in this release
- The original graph-regularized design in Scissor remains a key strength for stabilizing related cells, but on atlas-scale data, where many near-redundant cells accumulate within a much denser SNN graph, the binarized neighborhood structure can encourage smoother coefficient sharing across broader local neighborhoods rather than sharply isolated single-cell selection; in such settings, a metacell-first workflow is often a pragmatic execution strategy before running Scissor or LargeScissor. This is also consistent with the upstream discussions on large matrices and high RAM usage, together with the original author's public suggestion to merge cells into pseudo-cells or metacells before running Scissor:
  - https://github.com/sunduanchen/Scissor/discussions/3
  - https://github.com/sunduanchen/Scissor/discussions/64
  - https://github.com/sunduanchen/Scissor/issues/50#issuecomment-1527472087

## Relationship to Scissor

LargeScissor is built on top of the original Scissor project and should be understood as a compatibility-focused, memory-aware fork rather than a new method.

- Original Scissor repository: https://github.com/sunduanchen/Scissor
- Original tutorial: https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html

## How to Cite

If you use LargeScissor in academic work, please cite the original Scissor paper and the LargeScissor repository as follows:

1. Scissor:
   Sun, D., Guan, X., Moran, A. E., Wu, L. Y., Qian, D. Z., Schedin, P., Dai, M. S., Danilov, A. V., Alumkal, J. J., Adey, A. C., Spellman, P. T., & Xia, Z. (2021). Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nature Biotechnology. https://doi.org/10.1038/s41587-021-01091-3

2. LargeScissor:
   lxpsxx. LargeScissor: Optimized Scissor for Large-Scale Single-Cell Data. GitHub repository, version v1.2.0, 2026. https://github.com/lxpsxx/LargeScissor

## License

LargeScissor is released under the GNU General Public License v3.0, consistent with the upstream Scissor project.
