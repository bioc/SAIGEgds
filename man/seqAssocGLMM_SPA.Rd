\name{seqAssocGLMM_SPA}
\alias{seqAssocGLMM_SPA}
\title{P-value calculation}
\description{
    P-value calculations using variance approximation and an adjustment of
Saddlepoint approximation in the mixed model framework.
}

\usage{
seqAssocGLMM_SPA(gdsfile, modobj, maf=NaN, mac=10, missing=0.05, spa=TRUE,
    dsnode="", geno.ploidy=2L, res.savefn="", res.compress="ZIP",
    parallel=FALSE, load.balancing=TRUE, verbose=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{maf}{minor allele frequency threshold (checking >= maf), \code{NaN}
        for no filter}
    \item{mac}{minor allele count threshold (checking >= mac), \code{NaN}
        for no filter}
    \item{missing}{missing threshold for variants (checking <= missing),
        \code{NaN} for no filter}
    \item{spa}{\code{TRUE} for using the Saddlepoint approximation method for
        adjusting case-control imbalance}
    \item{dsnode}{"" for automatically searching the GDS nodes "genotype" and
        "annotation/format/DS", or use a user-defined GDS node in the file}
    \item{geno.ploidy}{specify the ploidy (2 by default); 0 or NA used for
        non-genotype data (e.g., CNV and \code{dsnode} should be specified
        for CNVs meanwhile)}
    \item{res.savefn}{an RData or GDS file name, "" for no saving}
    \item{res.compress}{the compression method for the output file, it should
        be one of ZIP, ZIP_RA, LZMA, LZMA_RA and none; see
        \code{\link{compression.gdsn}} for more details}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), a numeric value for the number of cores, or other value;
        \code{parallel} is passed to the argument \code{cl} in
        \code{\link{seqParallel}}, see \code{\link{seqParallel}} for more
        details}
    \item{load.balancing}{load balancing for the calculation in parallel
        if \code{TRUE}}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    For more details of SAIGE algorithm, please refer to the SAIGE paper
[Zhou et al. 2018] (see the reference section).
}

\value{
    Return a \code{data.frame} with the following components if not saving to
a file:
    \item{id}{variant ID in the GDS file;}
    \item{chr}{chromosome;}
    \item{pos}{position;}
    \item{rs.id}{the RS IDs if it is available in the GDS file;}
    \item{ref}{the reference allele;}
    \item{alt}{the alternative allele;}
    \item{AF.alt}{allele frequency for the alternative allele; the minor allele
        frequency is \code{pmin(AF.alt, 1-AF.alt)};}
    \item{mac}{minor allele count; the allele count for the alternative allele
        is \code{ifelse(AF.alt<=0.5, mac, 2*num-mac)};}
    \item{num}{the number of samples with non-missing genotypes;}
    \item{beta}{beta coefficient, odds ratio if binary outcomes (alternative
        allele vs. reference allele);}
    \item{SE}{standard error for beta coefficient;}
    \item{pval}{adjusted p-value with the Saddlepoint approximation method;}
    \item{p.norm}{p-values based on asymptotic normality (could be 0 if it
        is too small, e.g., \code{pnorm(-50) = 0} in R; used for checking only}
    \item{converged}{whether the SPA algorithm converges or not for adjusted
        p-values.}
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
    \code{\link{seqAssocGLMM_SPA}}, \code{\link{seqSAIGE_LoadPval}}
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
assoc <- seqAssocGLMM_SPA(gdsfile, glmm, mac=10)
head(assoc)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
