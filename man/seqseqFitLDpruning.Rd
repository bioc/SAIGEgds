\name{seqFitLDpruning}
\alias{seqFitLDpruning}
\title{Linkage disequilibrium pruning}
\description{
    Construct LD-pruned SNP sets for genetic relationship matrix (GRM).
}

\usage{
seqFitLDpruning(gdsfile, sample.id=NULL, variant.id=NULL,
    ld.threshold=0.1, maf=0.01, missing.rate=0.005, autosome.only=TRUE,
    use.cateMAC=TRUE, num.marker=100L, num.total=100000L, save.gdsfn=NULL,
    seed=200L, parallel=FALSE, parallel.multi=NULL, verbose=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS file name, a GDS object for genotypes, or
        a character vector for a list of GDS file names when genotyeps are
        split by chromosomes}
    \item{sample.id}{NULL for all samples, or sample IDs in GRM}
    \item{variant.id}{a list of variant IDs, used to construct GRM}
    \item{ld.threshold}{the LD threshold; see \code{snpgdsLDpruning()}}
    \item{maf}{minor allele frequency threshold for genotypes in \code{gdsfile}
        (checking >= maf), if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{missing.rate}{threshold of missing rate (checking <= missing.rate),
        if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{autosome.only}{\code{TRUE} for using autosomes only}
    \item{use.cateMAC}{\code{FALSE}, to use a single global variance ratio;
        \code{TRUE} (equal to
        \code{use.cateMAC=c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)}) for MAC
        categories (0, 1.5), [1.5, 2.5), ... [10.5, 20.5) and [20.5, Inf);
        or a numeric vector (strictly increasing) for unique cut points}
    \item{num.marker}{the number of SNPs used to calculate the variance ratio
        in each MAC category}
    \item{num.total}{the total number of LD-pruned variants excluding ultra
        rare variants; if the number of variants selected by the process of
        LD pruning is larger than \code{num.total}, the random set is used}
    \item{save.gdsfn}{if a file name is specified, construct a GDS genotype
        file to include all selected variants}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), a numeric value for the number of cores, or other value;
        \code{parallel} is passed to the argument \code{cl} in
        \code{\link{seqParallel}}, see \code{\link{seqParallel}} for more
        details}
    \item{parallel.multi}{only applicable when there are multiple input files
        in \code{gdsfile} and \code{save.gdsfn} is used; multiple GDS files can
        be extracted to create a sub GDS file in parallel via
        \code{parallel.multi}; if \code{NULL}, use \code{parallel} instead;
        \code{parallel.multi} could have less CPU cores than \code{parallel}
        to reduce memory usage}
    \item{seed}{an integer passed to \code{set.seed()} as a seed for random
        numbers, or \code{NULL} for no action}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    This function calls \code{\link[SNPRelate]{snpgdsLDpruning}} in the
SNPRelate package to perform linkage disequilibrium (LD) pruning. When
\code{use.cateMAC} is not \code{FALSE}, the ultra rare variants will be
selected according to the MAC categories, which could be used in the null
model fitting.
}

\value{
    Returns variant IDs or a list of variant IDs for multiple input GDS files.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqAssocGLMM_SPA}},
    \code{\link[SNPRelate]{snpgdsLDpruning}}
}

\examples{
library(SeqArray)
library(SAIGEgds)

# open a GDS file
gds_fn <- seqExampleFileName("KG_Phase1")

seqFitLDpruning(gds_fn, save.gdsfn="grm.gds")

# delete the temporary file
unlink("grm.gds", force=TRUE)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
