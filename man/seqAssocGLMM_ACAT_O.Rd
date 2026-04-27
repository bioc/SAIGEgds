\name{seqAssocGLMM_ACAT_O}
\alias{seqAssocGLMM_ACAT_O}
\alias{seqAssocGLMM_ACAT_O_GT}
\title{ACAT-O tests}
\description{
    Set-based ACAT-O tests using mixed models.
}

\usage{
seqAssocGLMM_ACAT_O(gdsfile, modobj, units, maxMAF=0.01, wbeta=AggrParamBeta,
    missing=0.05, collapse.mac=10, collapse.method=c("max", "sum"),
    ccimb.adj=TRUE, ER.mac=4.5, dsnode="",
    geno.model=c("additive", "dominant", "recessive"),
    res.savefn="", res.compress="ZIP",
    parallel=FALSE, verbose=TRUE, verbose.maf=FALSE)
seqAssocGLMM_ACAT_O_GT(gt_lst, modobj, maxMAF=0.01, wbeta=AggrParamBeta,
    missing=0.05, collapse.mac=10, collapse.method=c("max", "sum"),
    ccimb.adj=TRUE, ER.mac=4.5,
    geno.model=c("additive", "dominant", "recessive"),
    verbose=TRUE, verbose.maf=FALSE)
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
    \item{collapse.mac}{a threshold of minor allele count for collapsing
        ultra rare variants used in ACAT-V and SKAT tests if
        \code{mac <= collapse.mac}, 10 by default}
    \item{collapse.method}{\code{"max"} (by default),
        \code{"sum"} }
    \item{ccimb.adj}{whether adjusting for case-control imbalance}
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
    \item{gt_lst}{a list of numeric genotype matrices (sample by variant)
        or \code{Matrix} objects, where rows correspond to the samples in
        the model and columns to variants; each list element is a unit}
    \item{verbose}{if \code{TRUE}, show information}
    \item{verbose.maf}{if \code{TRUE}, show summary of MAFs in units}
}

\details{
    \code{seqAssocGLMM_ACAT_O()} reads genotype data from a GDS file and
supports parallel computation. \code{seqAssocGLMM_ACAT_O_GT()} takes a list
of in-memory genotype matrices \code{gt_lst} directly, which is useful when
genotypes are already loaded or constructed outside of the GDS framework.
The rows of each matrix must match the samples used in the null model.

    \code{seqUnitFilterCond()} in the SeqArray package can be used to restrict
the variant sets \code{units} to a range of MAC and/or MAF. ACAT-O combines
the p-values from ACAT-V, burden and SKAT together to calculate the final
p-value via Cauchy distribution.

    For more details of the ACAT-O method, please refer to the ACAT paper
[Liu et al. 2019] (see the reference section).

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
    \item{pval}{ACAT-O p-value;}
    \item{n_collapse}{the number of ultra rare variants for collapsing burden
        tests (\code{<= collapse.mac});}
    \item{p.burden}{burden tests;}
    \item{p.skat}{SKAT tests;}
    \item{p.acatv}{ACAT-V tests;}
    \item{burden.beta}{beta coefficient of burden test;}
    \item{burden.se}{SE of burden beta coefficient;}
}

\references{
    Liu Y., Chen S., Li Z., Morrison A.C., Boerwinkle E., Lin X.
ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant
Analysis in Sequencing Studies. Am J Hum Genetics 104, 410-421 (2019).
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAssocGLMM_ACAT_V}},
    \code{\link{seqAssocGLMM_Burden}},
    \code{\link{seqAssocGLMM_SKAT}},
    \code{\link[SeqArray]{seqUnitFilterCond}}
}

\examples{
# open a GDS file
fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
gdsfile <- seqOpen(fn)

sp_grm_fn <- system.file("extdata", "grm1k_10k_sp_grm.rds", package="SAIGEgds")
sp_grm <- readRDS(sp_grm_fn)

# load phenotype
phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)
head(pheno)

# fit the null model
glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile, sp_grm,
    trait.type="binary")

# get a list of variant units for burden tests
units <- seqUnitSlidingWindows(gdsfile, win.size=50, win.shift=250)

a <- seqAssocGLMM_ACAT_O(gdsfile, glmm, units, maxMAF=c(0.1,0.01))
head(a)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
