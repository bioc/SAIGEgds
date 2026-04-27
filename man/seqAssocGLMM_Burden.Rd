\name{seqAssocGLMM_Burden}
\alias{seqAssocGLMM_Burden}
\title{Burden tests}
\description{
    Set-based burden tests using mixed models.
}

\usage{
seqAssocGLMM_Burden(gdsfile, modobj, units, maxMAF=0.01, wbeta=AggrParamBeta,
    missing=0.05, ccimb.adj=TRUE, ER.mac=4.5, dsnode="",
    geno.model=c("additive", "dominant", "recessive"),
    res.savefn="", res.compress="ZIP", parallel=FALSE, verbose=TRUE,
    verbose.maf=FALSE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{units}{a list of units of selected variants, with S3 class
        \code{"SeqUnitListClass"} defined in the SeqArray package}
    \item{maxMAF}{a numeric vector for MAF thresholds, e.g.,
        \code{c(0.01, 0.001)} for rare variants, or \code{c(0.5, 0.1, 0.01)}
        for both common and rare variants}
    \item{wbeta}{weights for per-variant effect, using beta distribution
        \code{dbeta()} according to variant's MAF; a length-two vector, or a
        matrix with two rows for multiple beta parameters; by default, using
        beta(1,1) and beta(1,25) both}
    \item{missing}{missing threshold for variants (checking <= \code{missing}),
        \code{NaN} for no filter}
    \item{ccimb.adj}{whether adjusting for case-control imbalance or not}
    \item{ER.mac}{the maximum sum of MAC for the variants with efficient
        resampling p-values}
    \item{dsnode}{"" for automatically searching the GDS nodes "genotype" and
        "annotation/format/DS", or use a user-defined GDS node in the file}
    \item{geno.model}{a character string for the genetic model: \code{"additive"}
        (default), \code{"dominant"} or \code{"recessive"}; the dominant model
        recodes dosages as \code{min(G, 1)}, and the recessive model recodes
        as \code{max(G-1, 0)}}
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

    When \code{geno.model="dominant"}, the dosages are recoded as
\code{min(G, 1)} after allele flipping to the minor allele, testing whether
having at least one copy of the minor allele is associated with the trait.
When \code{geno.model="recessive"}, dosages are recoded as \code{max(G-1, 0)},
testing whether having two copies of the minor allele is associated.
MAF/MAC filtering and the reported allele frequencies use the original
(additive) genotype values regardless of the genetic model.
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
    \item{beta}{beta coefficient, log odds ratio if binary outcomes
        (alternative allele vs. reference allele);}
    \item{SE}{standard error for beta coefficient;}
    \item{pval}{final p-value (if the Saddlepoint approximation or efficient
        resampling is applied);}
    \item{method}{\code{"Normal"}: \code{pval} is calculated by assuming
        asymptotic normality; \code{"SPA"}: using the Saddlepoint approximation
        method; \code{"ER"}: using the efficient resampling}
    \item{p.norm}{p-values based on asymptotic normality (could be 0 if it
        is too small, e.g., \code{pnorm(-50) = 0} in R; used for checking only}
    \item{converged}{whether the SPA algorithm converges or not for p-values.}
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAssocGLMM_ACAT_V}}, \code{\link{seqAssocGLMM_ACAT_O}},
    \code{\link{seqAssocGLMM_SKAT}}, \code{\link[SeqArray]{seqUnitFilterCond}}
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

a <- seqAssocGLMM_Burden(gdsfile, glmm, units, maxMAF=c(0.1,0.01))
head(a)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
