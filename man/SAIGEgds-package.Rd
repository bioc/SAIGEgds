\name{SAIGEgds-package}
\alias{SAIGEgds-package}
\alias{SAIGEgds}
\docType{package}
\title{Parallel Computing Toolset for Genome-Wide Association Studies}
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
    Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
    A High-performance Computing Toolset for Relatedness and Principal
    Component Analysis of SNP Data.
    Bioinformatics (2012); doi: 10.1093/bioinformatics/bts610
}
\author{
    Xiuwen Zheng \email{xiuwen.zheng@gmail.com}
}

\examples{
####
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
