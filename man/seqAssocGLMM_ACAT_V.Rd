\name{seqAssocGLMM_ACAT_V}
\alias{seqAssocGLMM_ACAT_V}
\title{ACAT-V tests}
\description{
    Set-based ACAT-V tests using mixed models.
}

\usage{
seqAssocGLMM_ACAT_V(gdsfile, modobj, units, maxMAF=0.01, wbeta=AggrParamBeta,
    missing=0.05, ccimb.adj=TRUE, collapse.mac=10, ER.mac=4.5, dsnode="",
    res.savefn="", res.compress="ZIP", parallel=FALSE, verbose=TRUE,
    verbose.maf=FALSE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{units}{a list of units of selected variants, with S3 class
        \code{"SeqUnitListClass"} defined in the SeqArray package}
    \item{maxMAF}{a numeric vector for MAF thresholds, e.g.,
        \code{c(0.01, 0.001)} for rare variants}
    \item{wbeta}{weights for per-variant effect, using beta distribution
        \code{dbeta()} according to variant's MAF; a length-two vector, or a
        matrix with two rows for multiple beta parameters; by default, using
        beta(1,1) and beta(1,25) both}
    \item{missing}{missing threshold for variants (checking <= \code{missing}),
        \code{NaN} for no filter}
    \item{ccimb.adj}{whether adjusting for case-control imbalance or not}
    \item{collapse.mac}{a threshold of minor allele count for collapsing ultra
        rare variants if \code{mac <= collapse.mac}, 10 by default}
    \item{ER.mac}{the maximum sum of MAC for the variants with efficient
        resampling p-values}
    \item{dsnode}{"" for automatically searching the GDS nodes "genotype" and
        "annotation/format/DS", or use a user-defined GDS node in the file}
    \item{res.savefn}{an RData or GDS file name, "" for no saving}
    \item{res.compress}{the compression method for the output file, it should
        be one of LZMA, LZMA_RA, ZIP, ZIP_RA and none}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), a numeric value for the number of cores, or other value;
        \code{parallel} is passed to the argument \code{cl} in
        \code{\link[SeqArray]{seqParallel}}, see
        \code{\link[SeqArray]{seqParallel}} for more details}
    \item{verbose}{if \code{TRUE}, show information}
    \item{verbose.maf}{if \code{TRUE}, show summary of MAFs in units}
}

\details{
    \code{seqUnitFilterCond()} in the SeqArray package can be used to restrict
the variant sets \code{units} to a range of MAC and/or MAF.

    For more details of the ACAT-V method, please refer to the ACAT paper
[Liu et al. 2019] (see the reference section).
}

\value{
    Return a \code{data.frame} with the following components if not saving to
a file:
    \item{chr}{chromosome;}
    \item{start}{the starting position;}
    \item{end}{the ending position;}
    \item{maxMAF}{maximum of MAFs;}
    \item{numvar}{the number of variants in the window; if
        \code{weight="Cauchy"}, it is the number of individual p-values in
        the Cauchy test}
    \item{macmin}{the minimum of MAC in the window;}
    \item{macmed}{the median of MAC in the window;}
    \item{macmax}{the maximum of MAC in the window;}
    \item{summac}{the sum of minor allele counts;}
    \item{weight}{the weighting strategy with the beta distribution parameters}
    \item{n_single}{the number of variants running for single variant tests;}
    \item{n_collapse}{the number of ultra rare variants for collapsing burden
        tests (\code{<= collapse.mac});}
    \item{pval}{ACAT-V p-value;}
}

\references{
    Liu Y., Chen S., Li Z., Morrison A.C., Boerwinkle E., Lin X.
ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant
Analysis in Sequencing Studies. Am J Hum Genetics 104, 410-421 (2019).
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAssocGLMM_ACAT_O}},
    \code{\link{seqAssocGLMM_Burden}},
    \code{\link{seqAssocGLMM_SKAT}},
    \code{\link[SeqArray]{seqUnitFilterCond}}
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

a <- seqAssocGLMM_ACAT_V(gdsfile, glmm, units, maxMAF=c(0.1,0.01))
head(a)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
