\name{seqGetGenoByGeno}
\alias{seqGetGenoByGeno}
\title{Genotype-by-Genotype Interaction Matrix}
\description{
    Calculate the matrix for genotype-by-genotype interaction terms between
two sets of variants, excluding interaction columns without variation.
}

\usage{
seqGetGenoByGeno(gt1, gt2, use_colnm=TRUE)
}

\arguments{
    \item{gt1}{a numeric matrix or Matrix object (samples by variants) for the
        first set of genotypes}
    \item{gt2}{a numeric matrix or Matrix object (samples by variants) for the
        second set of genotypes}
    \item{use_colnm}{if \code{TRUE}, set column names in the output matrix
        using the column names of \code{gt1} and \code{gt2}; if column names
        are missing, default names \code{"a1"}, \code{"a2"}, ... and
        \code{"b1"}, \code{"b2"}, ... are used}
}

\details{
    The function computes element-wise products for all pairs of columns
(one from \code{gt1}, one from \code{gt2}) to form interaction terms.
Interaction columns with no variation (i.e., all zeros) are excluded from
the result. The computation is performed using sparse matrix operations
for efficiency, and the inputs are swapped internally if \code{gt1} has
more columns than \code{gt2}.
}

\value{
    A sparse matrix (\code{dgCMatrix}) where each column represents the
interaction between a variant in \code{gt1} and a variant in \code{gt2}.
Column names are in the format \code{"name1_x_name2"} when
\code{use_colnm=TRUE}.
}

\author{Xiuwen Zheng}

\seealso{
    \code{\link{seqAssocGLMM_SPA}}
}

\examples{
library(Matrix)

# dense matrix input
g1 <- matrix(c(0,1,2,0,1,0), nrow=3)
g2 <- matrix(c(1,0,1,2,1,0), nrow=3)
colnames(g1) <- c("snp1", "snp2")
colnames(g2) <- c("snp3", "snp4")

m <- seqGetGenoByGeno(g1, g2)
m

# sparse matrix input
g1s <- as(g1, "dgCMatrix")
g2s <- as(g2, "dgCMatrix")

m2 <- seqGetGenoByGeno(g1s, g2s)
m2
}

\keyword{GDS}
\keyword{genetics}
