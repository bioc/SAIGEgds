\name{seqFitSparseGRM}
\alias{seqFitSparseGRM}
\alias{seqFitDenseGRM}
\title{Sparse Genetic relationship matrix}
\description{
    Construct sparse genetic relationship matrix (GRM) from the full GRM.
}

\usage{
seqFitSparseGRM(gdsfile, sample.id=NULL, variant.id=NULL, nsnp.sub.random=2000L,
    rel.cutoff=0.125, maf=0.01, missing.rate=0.01, num.thread=1L,
    return.ID=FALSE, verbose=TRUE)
seqFitDenseGRM(gdsfile, sample.id=NULL, variant.id=NULL, maf=0.01,
    missing.rate=0.01, num.thread=1L, use.double=TRUE, return.ID=FALSE,
    verbose=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename or a GDS object for genotypes used
        in the GRM calculation; see details}
    \item{sample.id}{NULL for all samples, or sample IDs in GRM}
    \item{variant.id}{candidate variant IDs used for constructing GRM; or
        \code{NULL} for using all available candidate variants}
    \item{nsnp.sub.random}{the number of SNP markers randomly selected from the
        candidate variants for the initial scan of relatedness;
        if \code{nsnp.sub.random=0}, use all candidate SNPs}
    \item{rel.cutoff}{relatedness threshold for treating two individuals as
        unrelated (check >= rel.cutoff); \code{NaN} or \code{-Inf} for no
        threshold}
    \item{maf}{minor allele frequency threshold for genotypes in \code{gdsfile}
        (checking >= maf), if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{missing.rate}{threshold of missing rate (checking <= missing.rate),
        if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{num.thread}{the number of threads, 1 by default}
    \item{use.double}{\code{TRUE} for using 64-bit double precision to
        calculate GRM; otherwise, to use 32-bit single precision}
    \item{return.ID}{if \code{TRUE}, return the IDs for samples, the full
        variant set and the subset of variants}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    The genetic relationship matrix (GRM) is defined as
$grm_{ij} = avg_l [(g_{il} - 2*p_l)*(g_{jl} - 2*p_l) / 2*p_l*(1 - p_l)]$ for
individuals i,j and locus l, where $g_{il}$ is 0, 1 or 2, and $p_l$ is the
allele frequency at locus l. The missing genotypes are dropped from the
calculation.
}

\value{
    If \code{return.ID=TRUE}, returns a list with \code{sample.id} for
sample IDs, \code{variant.id} for the full set of variants,
\code{variant.sub.id} for the subset of variants, and the GRM matrix.
    Otherwise, it returns a sparse symmetric matrix for GRM, with sample IDs in
\code{colnames()} and \code{rownames()}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}
}

\examples{
library(Matrix)

# open a GDS file
gds_fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
gdsfile <- seqOpen(gds_fn)

seqSetFilter(gdsfile, variant.sel=1:100)
m <- seqFitSparseGRM(gdsfile, rel.cutoff=0.125)

is(m)
nnzero(m)                 # num of non-zero
nnzero(m) / prod(dim(m))  # percentage of non-zero

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
