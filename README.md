# LargeScissor: Optimized Scissor for Large-Scale Single-Cell Data

### About LargeScissor
The network regularization (Omega matrix) code used in the official `Scissor` R package often encounters Out-of-Memory (OOM) errors when processing large-scale single-cell datasets (e.g., >100,000 cells). 

**LargeScissor** Spent half of the Chinese New Year holiday on LargeScissor, I didn't do anything else, basically just three things: First Seurat V5 compatibility, Second preventing dense gene expression matrices, Third keeping the network matrix strictly sparse. I'm quite ashamed, just did a little bit of tiny work. Thank you everyone!

### Installation for LargeScissor
You can install this memory-optimized version directly from GitHub: devtools::install_github('lxpsxx/LargeScissor')

### Cite and acknowledge: https://github.com/sunduanchen/Scissor
Sun, D., Guan, X., Moran, A. E., Wu, L. Y., Qian, D. Z., Schedin, P., Dai, M. S., Danilov, A. V., Alumkal, J. J., Adey, A. C., Spellman, P. T., & Xia, Z. (2022). Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nature biotechnology, 40(4), 527–538. https://doi.org/10.1038/s41587-021-01091-3

Sincere thanks to Sun and others, the original authors of Scissor. The following is the authors' biography page.


# Original Scissor Documentation: Single-Cell Identification of Subpopulations with bulk Sample phenOtype coRrelation #

### Introduction ###
`Scissor` is a novel approach that utilizes the phenotypes, such as disease stage, tumor metastasis, treatment response, and survival outcomes, collected from bulk assays to identify the most highly phenotype-associated cell subpopulations from single-cell data. The workflow of Scissor is shown in the following Figure:

<p align="center">
<img src=Figure_Method.jpg height="702" width="600">
</p>

### News ###
* May, 2021: Scissor version 2.1.0 is updated.  
    + Add utilities for cell level evaludations including correlation check and bootstrap (function: evaluate.cell)
* Feb, 2021: Scissor version 2.0.0 is launched.  
    + Optimize the inputs and outputs in Scissor main function
    + Add utilities for the reliability significance test (function: reliability.test)
* Jun, 2020: Scissor version 1.0.0 is launched.

### Installation ###
* Prerequisites:
Scissor is developed under R (*version >= 3.6.1*). The [Seurat](https://satijalab.org/seurat/) package (*version >= 3.2.0*) is used for loading data and preprocessing.

* Latest version: The latest developmental version of Scissor can be downloaded from GitHub and installed from source by
`devtools::install_github('sunduanchen/Scissor')`

### Manual ###
Please see https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html for details. In the R terminal, please use the command `?Scissor` to access the help documents.

### Examples ###
In our [Scissor Tutorial](https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html), we use several applications on the Lung Adenocarcinoma (LUAD) scRNA-seq cancer cells as examples to show how to execute Scissor in real applications.

### How to cite `Scissor` ###
Please cite the following manuscript:

> *Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data*. Nature Biotechnology (2021). https://doi.org/10.1038/s41587-021-01091-3.    
Duanchen Sun, Xiangnan Guan, Amy E. Moran, Ling-Yun Wu, David Z. Qian, Pepper Schedin, Mu-Shui Dai, Alexey V. Danilov, Joshi J. Alumkal, Andrew C. Adey, Paul T. Spellman and Zheng Xia<br />

### License ###
Scissor is licensed under the GNU General Public License v3.0.

Improvements and new features of Scissor will be updated on a regular basis. Please post on the [GitHub discussion page](https://github.com/sunduanchen/Scissor/discussions) with any questions.
