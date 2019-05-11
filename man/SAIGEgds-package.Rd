\name{SAIGEgds-package}
\alias{SAIGEgds-package}
\alias{SAIGEgds}
\docType{package}
\title{Scalable and Accurate Implementation of Generalized mixed models using GDS files}
\description{
    Genome-wide association studies are widely used to investigate
the genetic basis of diseases and traits, but they pose many computational
challenges. We developed SNPRelate (R package for multi-core symmetric
multiprocessing computer architectures) to accelerate two key computations
on SNP data: principal component analysis (PCA) and relatedness analysis
using identity-by-descent measures. The kernels of our algorithms are
written in C/C++ and highly optimized.
}
\details{
    \tabular{ll}{
        Package: \tab SAIGEgds \cr
        Type: \tab Package \cr
        License: \tab GPL version 3 \cr
        Depends: \tab gdsfmt (>= 1.0.4) \cr
    }
    The genotypes stored in GDS format can be analyzed by the R functions in
    SNPRelate, which utilize the multi-core feature of machine for a single
    computer.

    Webpage: \url{http://github.com/zhengxwen/SAIGEgds}

    Tutorial: \url{http://corearray.sourceforge.net/tutorials/SNPRelate/}
}
\references{
    Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN, LeFaive J,
VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ, Denny JC, Lin M, Hveem K,
Kang HM, Abecasis GR, Willer CJ, Lee S. Efficiently controlling for case-control
imbalance and sample relatedness in large-scale genetic association studies.
*Nat Genet* (2018). Sep;50(9):1335-1341.

    Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D.
SeqArray -- A storage-efficient high-performance data format for WGS variant calls.
*Bioinformatics* (2017). DOI: 10.1093/bioinformatics/btx145.
}
\author{
    Xiuwen Zheng \email{xiuwen.zheng@abbvie.com},
    Wei Zhou, contributor (the original author of the SAIGE R package)
}

\examples{
# open a GDS file
fn <- system.file("extdata/grm1k_10k_snp.gds", package="SAIGEgds")
gdsfile <- seqOpen(fn)

# load phenotype
phenofn <- system.file("extdata/pheno.txt.gz", package="SAIGEgds")
pheno <- read.table(phenofn, head=TRUE, as.is=TRUE)

# fit the null model
glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)

# p-value calculation
assoc <- seqAssocGLMM_SPA(gdsfile, glmm, mac=4)

head(assoc)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
