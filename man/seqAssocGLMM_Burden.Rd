\name{seqAssocGLMM_Burden}
\alias{seqAssocGLMM_Burden}
\title{Burden tests}
\description{
    Burden p-value calculations using mixed models and the
Saddlepoint approximation method for case-control imbalance.
}

\usage{
seqAssocGLMM_Burden(gdsfile, modobj, units, wbeta=AggrParamBeta, ccimb.adj=TRUE,
    summac=3, dsnode="", res.savefn="", res.compress="LZMA", parallel=FALSE,
    verbose=TRUE, verbose.maf=FALSE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{units}{a list of units of selected variants, with S3 class
        \code{"SeqUnitListClass"} defined in the SeqArray package}
    \item{wbeta}{weights for per-variant effect, using beta distribution
        \code{dbeta()} according to variant's MAF; a length-two vector, or a
        matrix with two rows for multiple beta parameters; by default, using
        beta(1,1) and beta(1,25) both}
    \item{ccimb.adj}{whether adjusting for case-control imbalance or not}
    \item{summac}{a threshold for the sum of minor allele counts (MAC)
        (checking \code{>= summac}); no calculation if the sum of MAC
        \code{< summac}}
    \item{dsnode}{"" for automatically searching the GDS nodes "genotype" and
        "annotation/format/DS", or use a user-defined GDS node in the file}
    \item{res.savefn}{an RData or GDS file name, "" for no saving}
    \item{res.compress}{the compression method for the output file, it should
        be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), a numeric value for the number of cores, or other value;
        \code{parallel} is passed to the argument \code{cl} in
        \code{\link{seqParallel}}, see \code{\link{seqParallel}} for more
        details}
    \item{verbose}{if \code{TRUE}, show information}
    \item{verbose.maf}{if \code{TRUE}, show summary of MAFs in units}
}

\details{
    \code{seqUnitFilterCond()} in the SeqArray package can be used to restrict
the variant sets \code{units} to a range of MAC and/or MAF.
}

\value{
    Return a \code{data.frame} with the following components if not saving to
a file:
    \code{chr}, chromosome;
    \code{start}, the starting position;
    \code{end}, the ending position;
    \code{numvar}, the number of variants in the window;
    \code{maf.avg}, the average of MAFs in the window;
    \code{maf.sd}, the standard deviation of MAFs in the window;
    \code{maf.min}, the minimum of MAFs in the window;
    \code{maf.max}, the maximum of MAFs in the window;
    \code{mac.avg}, the average of MACs in the window;
    \code{mac.sd}, the standard deviation of MACs in the window;
    \code{mac.min}, the minimum of MACs in the window;
    \code{mac.max}, the maximum of MACs in the window;
    \code{summac}, the sum of minor allele counts;
    \code{beta}, beta coefficient, log odds ratio if binary outcomes;
    \code{SE}, standard error for beta coefficient;
    \code{pval}, p-value (adjusted p-value with Saddlepoint approximation
        if \code{ccimb.adj=TRUE} and binary outcomes);
    \item{p.norm}{p-values based on asymptotic normality (could be 0 if it
        is too small, e.g., \code{pnorm(-50) = 0} in R; used for checking only}
    \code{cvg}, whether the SPA algorithm converges or not for adjusted p-value.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAssocGLMM_ACAT_V}},
    \code{\link{seqAssocGLMM_ACAT_O}},
    \code{\link{seqAssocGLMM_SKAT}},
    \code{\link{seqUnitFilterCond}}
}

\examples{
# open a GDS file
fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
gdsfile <- seqOpen(fn)

# load phenotype
phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)
head(pheno)

# fit the null model
glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile, trait.type="binary")

# get a list of variant units for burden tests
units <- seqUnitSlidingWindows(gdsfile, win.size=500, win.shift=250)

assoc <- seqAssocGLMM_Burden(gdsfile, glmm, units)
head(assoc)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
