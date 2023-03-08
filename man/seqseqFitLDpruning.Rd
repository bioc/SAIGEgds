\name{seqFitLDpruning}
\alias{seqFitLDpruning}
\title{LD pruning}
\description{
    Construct LD-pruned SNP sets for genetic relationship matrix (GRM).
}

\usage{
seqFitLDpruning(gdsfile, sample.id=NULL, variant.id=NULL,
    ld.threshold=0.2, maf=0.005, missing.rate=0.005, use.cateMAC=FALSE,
    num.marker=30L, outfn=NULL, parallel=FALSE, verbose=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename or a GDS object for genotypes}
    \item{sample.id}{NULL for all samples, or sample IDs in GRM}
    \item{variant.id}{a list of variant IDs, used to construct GRM}
    \item{ld.threshold}{the LD threshold; see \code{snpgdsLDpruning()}}
    \item{maf}{minor allele frequency threshold for genotypes in \code{gdsfile}
        (checking >= maf), if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{missing.rate}{threshold of missing rate (checking <= missing.rate),
        if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{use.cateMAC}{\code{FALSE}, to use a single global variance ratio;
        \code{TRUE} (equal to
        \code{use.cateMAC=c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)}) for MAC
        categories (0, 1.5), [1.5, 2.5), ... [10.5, 20.5) and [20.5, Inf);
        or a numeric vector (strictly increasing) for unique cut points}
    \item{num.marker}{the number of SNPs used to calculate the variance ratio
        in each MAC category}
    \item{outfn}{x}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), a numeric value for the number of cores, or other value;
        \code{parallel} is passed to the argument \code{cl} in
        \code{\link{seqParallel}}, see \code{\link{seqParallel}} for more
        details}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    Utilizing the sparse structure of genotypes could significantly improve
the computational efficiency of model fitting, but it also increases the memory
usage.
    For more details of SAIGE algorithm, please refer to the SAIGE paper
[Zhou et al. 2018] (see the reference section).
}

\value{
    Returns a sparse matrix for GRM, with sample IDs in \code{colnames()}
and \code{rownames()}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqAssocGLMM_SPA}}
}

\examples{
library(SeqArray)
library(SAIGEgds)

# open a GDS file
gds_fn <- seqExampleFileName("KG_Phase1")

seqFitLDpruning(gds_fn, outfn="grm.gds", use.cateMAC=TRUE)

# delete the temporary file
unlink("grm.gds", force=TRUE)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
