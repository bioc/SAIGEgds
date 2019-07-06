SAIGEgds: Scalable Implementation of Generalized mixed models using GDS files
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

Scalable and accurate implementation of generalized mixed mode with the support of Genomic Data Structure ([GDS](https://github.com/zhengxwen/SeqArray)) files and highly optimized C++ implementation. It is designed for single variant tests in large-scale phenome-wide association studies (PheWAS) with millions of variants and hundreds of thousands of samples, e.g., [UK Biobank genotype data](https://www.ukbiobank.ac.uk).

The implementation of SAIGEgds is based on the original [SAIGE](https://github.com/weizhouUMICH/SAIGE) R package (v0.29.4.4). All of the calculation with single-precision floating-point numbers in [SAIGE](https://github.com/weizhouUMICH/SAIGE) are replaced by the double-precision calculation in SAIGEgds. SAIGEgds also implements some of the [SPAtest](https://cran.r-project.org/web/packages/SPAtest/index.html) functions in C to speed up the calculation of Saddlepoint Approximation.


## Citation

Zheng X, Davis J.Wade, SAIGEgds -- an efficient statistical tool for large-scale PheWAS with mixed models; (Abstract/Program #XX). Presented at the 69th Annual Meeting of The American Society of Human Genetics (ASHG), Oct 15-19, Houston, US.

Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN, LeFaive J, VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ, Denny JC, Lin M, Hveem K, Kang HM, Abecasis GR, Willer CJ, Lee S. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. *Nat Genet* (2018). Sep;50(9):1335-1341.

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D. SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics* (2017). [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## Installation

* Bioconductor repository:
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SAIGEgds")
```

* Development version from Github:
```R
library("devtools")
install_github("AbbVie-ComputationalGenomics/SAIGEgds")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.


## Examples

```R
library(SeqArray)
library(SAIGEgds)

# open a GDS file
fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
gdsfile <- seqOpen(fn)

# load phenotype
phenofn <- system.file("extdata/pheno.txt.gz", package="SAIGEgds")
pheno <- read.table(phenofn, head=TRUE, as.is=TRUE)
head(pheno)
##   y      x1 x2 sample.id
## 1 0  1.5118  1        s1
## 2 0  0.3898  1        s2
## 3 0 -0.6212  1        s3
## ...

# fit the null model
glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
## SAIGE association analysis:
## Filtering variants:
## [==================================================] 100%, completed (0s)
## Fit the null model: y ~ x1 + x2 + var(GRM)
##     # of samples: 1,000
##     # of variants: 9,976
##     using 1 thread
## ...

# p-value calculation
assoc <- seqAssocGLMM_SPA(gdsfile, glmm, mac=10)
## SAIGE association analysis:
##     # of samples: 1,000
##     # of variants: 10,000
##     p-value threshold for SPA adjustment: 0.05
##     variance ratio for approximation: 0.9410486
## [==================================================] 100%, completed (1s)
## # of variants after filtering MAF/MAC: 9,976
## Done.

head(assoc)
##   id chr pos rs.id ref alt AF.alt AC.alt  num        beta        SE      pval pval.noadj converged
## 1  1   1   1   rs1   1   2 0.0305     61 1000  0.60500665 0.4720839 0.1999950  0.1999950      TRUE
## 2  2   1   2   rs2   1   2 0.0380     76 1000 -0.09606419 0.4101637 0.8148224  0.8148224      TRUE
## 3  3   1   3   rs3   1   2 0.0215     43 1000 -0.55131755 0.5699739 0.3334101  0.3334101      TRUE
## ...

# close the GDS file
seqClose(gdsfile)
```


## Also See

[SeqArray](http://www.bioconductor.org/packages/SeqArray): Data Management of Large-scale Whole-genome Sequence Variant Calls
