\name{seqFitNullGLMM_SPA}
\alias{seqFitNullGLMM_SPA}
\alias{seqRefitNullGLMM}
\title{Fit the null model for the SAIGE mixed model framework}
\description{
    Fit the null model in the mixed model framework with an optional genetic
relationship matrix (GRM).
}

\usage{
seqFitNullGLMM_SPA(formula, data, gdsfile=NULL, grm.mat=NULL,
    trait.type=c("binary", "quantitative", "survival"), event.time=NULL,
    sample.col="sample.id", maf=0.01,
    missing.rate=0.01, max.num.snp=1000000L, variant.id=NULL,
    variant.id.varratio=NULL, nsnp.sub.random=2000L, rel.cutoff=0.125,
    inv.norm=c("residuals", "quant", "none"), use.cateMAC=FALSE,
    cateMAC.inc.maf=TRUE, cateMAC.simu=TRUE, use.offset=FALSE, X.transform=TRUE,
    tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5, maxiterPCG=500L,
    num.marker=30L, tau.init=c(0,0), traceCVcutoff=0.0025, ratioCVcutoff=0.001,
    geno.sparse=TRUE, use.gpu=FALSE, save.packed.geno=FALSE, num.thread=1L,
    model.savefn="", seed=200L, fork.loading, parallel.loading=FALSE,
    verbose=TRUE)
seqRefitNullGLMM(formula, data, model=NULL, event.time=NULL,
    sample.col="sample.id",
    inv.norm=c("residuals", "quant", "none"), use.offset=FALSE,
    X.transform=TRUE, tau.update=FALSE, recalcVR=FALSE,
    tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5, maxiterPCG=500L,
    num.marker=30L, traceCVcutoff=0.0025, ratioCVcutoff=0.001,
    num.thread=1L, seed=200L, use.gpu=FALSE, verbose=TRUE)
}

\arguments{
    \item{formula}{an object of class \code{formula} (or one that can be
        coerced to that class), e.g., \code{y ~ x1 + x2}, see
        \code{\link[stats]{lm}}}
    \item{data}{a data frame for the formulas}
    \item{gdsfile}{a SeqArray GDS filename or a GDS object for genotypes used
        in the GRM calculation; see details}
    \item{grm.mat}{NULL, a user-defined GRM or an object of 'snpgdsGRMClass'
        returned from \code{SNPRelate::snpgdsGRM()}: a dense numeric matrix or
        a sparse symmetric matrix 'sparseMatrix' defined in \code{Matrix};
        \code{colnames} and \code{rownames} should be sample IDs;
        or \code{TRUE}, call \code{seqFitSparseGRM()} to get a sparse GRM;
        see details}
    \item{trait.type}{"binary" for binary outcomes, "quantitative" for
        continuous outcomes, or "survival" for time-to-event outcomes fitted
        with a Cox proportional-hazards frailty model (via the Breslow
        Cox-Poisson equivalence, cf. the GATE method); for "survival" the
        response in \code{formula} is the 0/1 event status (1=event,
        0=censored) and \code{event.time} gives the follow-up time}
    \item{event.time}{only used when \code{trait.type="survival"}: a column
        name in \code{data} giving the non-negative event/censoring time; the
        intercept is removed (absorbed by the baseline hazard) so a GRM (via
        \code{gdsfile} or \code{grm.mat}) is required and \code{X.transform}/
        \code{use.offset} are forced to \code{FALSE}. For
        \code{seqRefitNullGLMM()} the event time is not stored in \code{model}
        and must be re-supplied here (as a column of the new \code{data})}
    \item{sample.col}{the column name of sample IDs corresponding to the GDS
        file}
    \item{maf}{minor allele frequency threshold for genotypes in \code{gdsfile}
        (checking >= maf), if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{missing.rate}{threshold of missing rate (checking <= missing.rate),
        if \code{variant.id=NULL}; \code{NaN} for no filter}
    \item{max.num.snp}{the maximum number of SNPs used, or -1 for no limit}
    \item{variant.id}{a list of variant IDs, considered to be used in GRM;
        or \code{NULL}, all variants in the GDS file to be the candidate
        variants}
    \item{variant.id.varratio}{a list of variant IDs, to be used in the
        variance ratio estimation; or \code{NULL} to use the candidate variants
        according to \code{variant.id}}
    \item{nsnp.sub.random}{used in \code{seqFitSparseGRM()} when
        \code{grm.mat=TRUE}: the number of SNP markers randomly selected from
        the candidate variants for the initial scan of relatedness}
    \item{rel.cutoff}{relatedness threshold for treating two individuals as
        unrelated (check >= rel.cutoff), \code{NaN} or \code{-Inf} for no
        threshold; only applicable when \code{grm.mat=TRUE}}
    \item{use.cateMAC}{\code{FALSE}, to use a single global variance ratio;
        \code{TRUE} (equal to
        \code{use.cateMAC=c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)}) for MAC
        categories (0, 1.5), [1.5, 2.5), ... [10.5, 20.5) and [20.5, Inf);
        or a numeric vector (strictly increasing) for unique cut points}
    \item{cateMAC.inc.maf}{\code{TRUE} to include a MAC cut point according to
        MAF=\code{maf}, or a numeric vector for MAF; only applicable if
        \code{use.cateMAC} is not \code{FALSE}}
    \item{cateMAC.simu}{\code{TRUE} for using simulated independent
        genotypes under Hardy-Weinberg equilibrium if there are not enough
        SNP markers in the MAC category; only applicable if \code{use.cateMAC}
        is not \code{FALSE}}
    \item{inv.norm}{\code{"residuals"} (by default), perform the inverse normal
        transformation on the residuals after fitting the fixed effects;
        \code{"quant"}, perform the inverse normal transformation on the
        outcome variable directly (the same behavior as the SAIGE package);
        \code{"none"}, no transformation.
        If \code{inv.norm=TRUE}, perform \code{inv.norm="residuals"};
        if \code{inv.norm=FALSE}, it is as the same as \code{inv.norm="none"};
        See the reference for more discussion [Sofer et al., 2019]}
    \item{use.offset}{if \code{TRUE}, use the estimated residual offset to
        avoid calculating each fixed effect coefficient}
    \item{X.transform}{if \code{TRUE}, perform QR decomposition on the design
        matrix}
    \item{tol}{overall tolerance for model fitting}
    \item{maxiter}{the maximum number of iterations for model fitting;
        \code{maxiter=0} for using a fixed tau parameter in \code{tau.init}}
    \item{nrun}{the number of random vectors in the trace estimation}
    \item{tolPCG}{tolerance of PCG iterations}
    \item{maxiterPCG}{the maximum number of PCG iterations}
    \item{num.marker}{the number of SNPs used to calculate the variance ratio}
    \item{tau.init}{a 2-length numeric vector, the initial values for variance
        components, tau; for binary traits, the first element is always set
        to 1; if \code{tau.init} is not specified, the second element will be
        0.5 for binary traits}
    \item{traceCVcutoff}{the threshold for coefficient of variation (CV) for
        the trace estimator, and the number of runs for trace estimation will
        be increased until the CV is below the threshold}
    \item{ratioCVcutoff}{the threshold for coefficient of variation (CV) for
        estimating the variance ratio, and the number of randomly selected
        markers will be increased until the CV is below the threshold}
    \item{geno.sparse}{if \code{TRUE} (by default), store the sparse structure
        for genotypes; otherwise, save genotypes in a 2-bit dense matrix;
        see details}
    \item{save.packed.geno}{if \code{TRUE}, save packed genotypes of GRM in the
        output model object; it can be used for refitting the null model with
        different fixed effects}
    \item{num.thread}{the number of threads, 1 by default}
    \item{model.savefn}{the filename of model output, R data file '.rda',
        '.RData', or '.rds'}
    \item{seed}{an integer passed to \code{set.seed()} as a seed for random
        numbers, or \code{NULL} for no action}
    \item{fork.loading}{deprecated, use \code{parallel.loading} instead}
    \item{parallel.loading}{load genotypes via parallel or not; multiple
        processes can reduce loading time of genotypes, but may increase the
        memory usage}
    \item{model}{a fitted null model object of class
        \code{"ClassSAIGE_NullModel"} created by
        \code{seqFitNullGLMM_SPA()} with \code{save.packed.geno=TRUE}}
    \item{tau.update}{if \code{TRUE}, re-estimate the variance components
        tau; otherwise use the tau values from the input model}
    \item{recalcVR}{if \code{TRUE}, recalculate the variance ratio; otherwise
        use the variance ratio from the input model}
    \item{use.gpu}{if \code{TRUE}, attempt to use OpenCL GPU acceleration
        for the GRM cross-product computation; falls back to CPU if no
        compatible GPU or OpenCL runtime is found}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    \code{seqFitNullGLMM_SPA()} fits the null model from scratch using
genotype data from a GDS file and/or a genetic relationship matrix (GRM).
    \code{seqRefitNullGLMM()} refits the null model using packed genotypes
saved from a previous call to \code{seqFitNullGLMM_SPA()} (with
\code{save.packed.geno=TRUE}), allowing different fixed-effect covariates
without reloading genotypes from the GDS file. By default, the variance
components (tau) are kept from the input model (\code{tau.update=FALSE}) and
the variance ratio is recalculated (\code{recalcVR=TRUE}).

    Utilizing the sparse structure of genotypes could significantly improve
the computational efficiency of model fitting, but it also increases the memory
usage.
    For more details of SAIGE algorithm, please refer to the SAIGE paper
[Zhou et al., 2018] (see the reference section).
}

\value{
    Returns a list with the following components:
    \item{coefficients}{the beta coefficients for fixed effects;}
    \item{tau}{a numeric vector of variance components 'Sigma_E' and 'Sigma_G';}
    \item{linear.predictors}{the linear fit on link scale;}
    \item{fitted.values}{fitted values from objects returned by modeling
        functions using \code{glm.fit};}
    \item{residuals}{residuals;}
    \item{converged}{whether the model is fitted or not;}
    \item{obj.noK}{internal use, returned object from the SPAtest package;}
    \item{Sigma_inv}{optional inverse covariance matrix when a GRM is used;}
    \item{chol_inv_X_Sigma}{optional projection matrix component when a GRM is used;}
    \item{var.ratio}{a data.frame with columns 'id' (variant.id), 'maf'
        (minor allele frequency), 'mac' (minor allele count),
        'var1' (the variance of score statistic),
        'var2' (a variance estimate without accounting for estimated random
        effects) and 'ratio' (var1/var2, estimated variance ratio for variance
        approximation);}
    \item{trait.type}{"binary", "quantitative" or "survival";}
    \item{sample.id}{the sample IDs used in the model fitting;}
    \item{variant.id}{the variant IDs used in the model fitting.}
}

\references{
    Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN,
LeFaive J, VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ,
Denny JC, Lin M, Hveem K, Kang HM, Abecasis GR, Willer CJ, Lee S.
Efficiently controlling for case-control imbalance and sample relatedness in
large-scale genetic association studies. Nat Genet (2018). Sep;50(9):1335-1341.

    Zhou W, Zhao Z, Nielsen JB, Fritsche LG, LeFaive J, Taliun SAG, Bi W,
Gabrielsen ME, Daly MJ, Neale BM, Hveem K, Abecasis GR, Willer CJ, et al.
Scalable generalized linear mixed model for region-based association tests
in large biobanks and cohorts. Nat Genet. 2020; 52: 634-9.

    Bi W, Fritsche LG, Mukherjee B, Kim S, Lee S. A fast and accurate method
for genome-wide time-to-event data analysis and its application to UK Biobank.
Am J Hum Genet. 2020; 107(2): 222-233. (the GATE method, used for
\code{trait.type="survival"})
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqAssocGLMM_SPA}}, \code{\link{seqFitSparseGRM}}
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
glmm

# fit a survival null model (Cox frailty via the GATE method):
# the response is the 0/1 event status and 'event.time' the follow-up time
set.seed(100)
pheno$status <- rbinom(nrow(pheno), 1, 0.5)         # event indicator
pheno$atime  <- rexp(nrow(pheno), exp(0.3*pheno$x1))  # event/censoring time
survmod <- seqFitNullGLMM_SPA(status ~ x1 + x2, pheno, gdsfile,
    trait.type="survival", event.time="atime")
survmod

# close the GDS file
seqClose(gdsfile)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
