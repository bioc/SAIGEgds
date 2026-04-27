\name{seqAssocGLMM_SPA}
\alias{seqAssocGLMM_SPA}
\alias{seqAssocGLMM_GT}
\title{P-value calculation}
\description{
    P-value calculations using variance approximation and an adjustment of
Saddlepoint approximation in the mixed model framework.
}

\usage{
seqAssocGLMM_SPA(gdsfile, modobj, maf=NaN, mac=10, missing=0.05, spa=TRUE,
    ER.mac=4.5, dsnode="",
    geno.model=c("additive", "dominant", "recessive"), geno.ploidy=2L,
    res.savefn="", res.compress="ZIP",
    parallel=FALSE, load.balancing=TRUE,
    verbose.pval=c(0, 5e-10, 5e-8, 5e-6, 5e-4, 1), verbose=TRUE)
seqAssocGLMM_GT(gt, modobj, spa=TRUE, ER.mac=4.5,
    geno.model=c("additive", "dominant", "recessive"), geno.ploidy=2L,
    verbose=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{maf}{minor allele frequency threshold (checking >= maf), \code{NaN}
        for no filter}
    \item{mac}{minor allele count threshold (checking >= mac), \code{NaN}
        for no filter}
    \item{missing}{missing threshold for variants (checking <= \code{missing}),
        \code{NaN} for no filter}
    \item{spa}{\code{TRUE} for using the Saddlepoint approximation method for
        adjusting case-control imbalance}
    \item{ER.mac}{the maximum MAC for the variants with efficient resampling
        p-values}
    \item{dsnode}{"" for automatically searching the GDS nodes "genotype" and
        "annotation/format/DS", or use a user-defined GDS node in the file}
    \item{geno.model}{a character string for the genetic model: \code{"additive"}
        (default), \code{"dominant"} or \code{"recessive"}; the dominant model
        recodes dosages as \code{min(G, 1)} (having >= 1 copy of the minor
        allele), and the recessive model recodes as \code{max(G-1, 0)} (having
        2 copies of the minor allele); only applicable when \code{geno.ploidy > 0}}
    \item{geno.ploidy}{specify the ploidy (2 by default); 0 or NA used for
        non-genotype data (e.g., CNV and \code{dsnode} should be specified
        for CNVs meanwhile)}
    \item{res.savefn}{an RData or GDS file name, "" for no saving}
    \item{res.compress}{the compression method for the output file, it should
        be one of ZIP, ZIP_RA, LZMA, LZMA_RA and none; see
        \code{\link[gdsfmt]{compression.gdsn}} for more details}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), a numeric value for the number of cores, or other value;
        \code{parallel} is passed to the argument \code{cl} in
        \code{\link[SeqArray]{seqParallel}}, see
        \code{\link[SeqArray]{seqParallel}} for more details}
    \item{load.balancing}{load balancing for the calculation in parallel
        if \code{TRUE}}
    \item{gt}{a numeric genotype matrix (sample by variant) or a Matrix
        object, where rows correspond to the samples in the model and columns
        to variants}
    \item{verbose.pval}{\code{NULL} or a numeric vector for dividing p-values
        into the intervals and show the corresponding counts}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    \code{seqAssocGLMM_SPA()} reads genotype data from a GDS file and
supports parallel computation. \code{seqAssocGLMM_GT()} takes an in-memory
numeric matrix \code{gt} (sample by variant) directly, which is useful when
genotypes are already loaded or generated outside of the GDS framework.
The rows of \code{gt} must match the samples used in the null model.

    For more details of SAIGE algorithm, please refer to the SAIGE paper
[Zhou et al. 2018] (see the reference section).

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
    \item{id}{variant ID in the GDS file;}
    \item{chr}{chromosome;}
    \item{pos}{position;}
    \item{rs.id}{the RS IDs if it is available in the GDS file as a node
        "annotation/id";}
    \item{ref}{the reference allele;}
    \item{alt}{the alternative allele;}
    \item{AF.alt}{allele frequency for the alternative allele; the minor allele
        frequency is \code{pmin(AF.alt, 1-AF.alt)};}
    \item{mac}{minor allele count; the allele count for the alternative allele
        is \code{ifelse(AF.alt<=0.5, mac, 2*num-mac)};}
    \item{num}{the number of samples with non-missing genotypes;}
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

\references{
    Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN,
LeFaive J, VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ,
Denny JC, Lin M, Hveem K, Kang HM, Abecasis GR, Willer CJ, Lee S.
Efficiently controlling for case-control imbalance and sample relatedness in
large-scale genetic association studies. Nat Genet (2018). Sep;50(9):1335-1341.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqSAIGE_LoadPval}}
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

# p-value calculation
a <- seqAssocGLMM_SPA(gdsfile, glmm, mac=10)
head(a)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
