\name{seqFitNullGLMM_SPA}
\alias{seqFitNullGLMM_SPA}
\title{Fit the null model with GRM}
\description{
    Fit the null model in the mixed frame with genetic relationship matrix (GRM).
}

\usage{
seqFitNullGLMM_SPA(formula, data, gdsfile,
    trait.type=c("binary", "quantitative"), maf=0.01, missing.rate=0.01,
    max.num.snp=1000000L, variant.id=NULL, inv.norm=TRUE, X.transform=FALSE,
    tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5, maxiterPCG=500L,
    num.marker=30L, tau.init=c(0,0), traceCVcutoff=1, ratioCVcutoff=1,
    geno.sparse=TRUE, num.thread=1L, model.savefn="", seed=200L, verbose=TRUE)
}

\arguments{
    \item{formula}{an object of class \code{formula} (or one that can be coerced to
          that class), e.g., \code{y ~ x1 + x2}, see \code{\link{lm}}}
    \item{data}{}
    \item{gdsfile}{}
    \item{trait.type}{"binary" for binary outcomes, "quantitative" for continuous
        outcomes}
    \item{maf}{minor allele frequency for imported genotypes, if \code{variant.id=NULL}}
    \item{missing.rate}{threshold of missing rate, if \code{variant.id=NULL}}
    \item{max.num.snp}{the maximum number of SNPs used}
    \item{variant.id}{a list of variant IDs, used to construct GRM}
    \item{inv.norm}{}
    \item{X.transform}{}
    \item{tol}{}
    \item{maxiter}{}
    \item{nrun}{}
    \item{tolPCG}{}
    \item{maxiterPCG}{}
    \item{num.marker}{the number of SNPs used to calculate the variance ratio}
    \item{tau.init}{}
    \item{traceCVcutoff}{}
    \item{ratioCVcutoff}{}
    \item{geno.sparse}{}
    \item{num.thread}{the number of threads}
    \item{model.savefn}{}
    \item{seed}{an integer as a seed for random numbers}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    XXXX
}

\value{
    An object of class \code{\link{gdsn.class}} of the new node.
}

\references{
    \url{http://github.com/zhengxwen/gdsfmt}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAssocGLMM_SPA}}
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


# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{utilities}
