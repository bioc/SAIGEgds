\name{seqAssocGLMM_SPA}
\alias{seqAssocGLMM_SPA}
\title{Fit the null model with GRM}
\description{
    Fit the null model in the mixed frame with genetic relationship matrix (GRM).
}

\usage{
seqAssocGLMM_SPA(gdsfile, modobj, maf=NaN, mac=10, dsnode="", spa.pval=0.05,
    res.savefn="", res.compress=c("LZMA", "LZMA_RA", "ZIP", "ZIP_RA", "none"),
    parallel=FALSE, verbose=TRUE)
}

\arguments{
    \item{gdsfile}{a SeqArray GDS filename, or a GDS object}
    \item{modobj}{an R object for SAIGE model parameters}
    \item{maf}{minor allele frequency cutoff}
    \item{mac}{minor allele count cutoff}
    \item{dsnode}{}
    \item{spa.pval}{the p-value threshold for SPA adjustment}
    \item{res.savefn}{an RData or GDS file name}
    \item{res.compress}{the compression method for the output file}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (multicore
        processing), numeric value or other value; \code{parallel} is passed
        to the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details.}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    XXXX
}

\value{
    Return a \code{data.frame} with the following components:
    \code{id}, variant ID in the GDS file;
    \code{chr}, chromosome;
    \code{pos}, position;
    \code{rs.id}, the RS IDs if it is available in the GDS file;
    \code{ref}, the reference allele;
    \code{alt}, the alternative allele;
    \code{AF.alt}, allele frequency for the alternative allele;
    \code{AC.alt}, allele count for the alternative allele;
    \code{num}, the number of samples with non-missing genotypes;
    \code{beta}, beta coefficient, effect size;
    \code{SE}, standard error for beta coefficient;
    \code{pval}, adjusted p-value with SPA;
    \code{pval.noadj}, p-value without adjustment;
    \code{converged}, whether the SPA algorithm converges or not for adjusted p-value
}

\references{
    Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN, LeFaive J,
    VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ, Denny JC, Lin M, Hveem K,
    Kang HM, Abecasis GR, Willer CJ, Lee S.
    Efficiently controlling for case-control imbalance and sample relatedness in
    large-scale genetic association studies.
    Nat Genet (2018). Sep;50(9):1335-1341.
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

# p-value calculation
assoc <- seqAssocGLMM_SPA(gdsfile, glmm, mac=4)

head(assoc)

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
