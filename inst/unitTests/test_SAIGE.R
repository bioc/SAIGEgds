#############################################################
#
# DESCRIPTION: test the SAIGE calculation
#

library(RUnit)
library(SAIGEgds)


# create an R object file for package checking
.create_test_file <- function()
{
	# 'Rounding' was the default in versions prior to R_3.6.0
	# it is used for reproduction of the results created by R (v3.5.2)
	tryCatch(suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")),
		error=function(e) FALSE)

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model for binary outcomes
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
	# save model
	saveRDS(glmm, file="saige_model.rds", compress="xz")
	# p-value calculation
	seqAssocGLMM_SPA(gdsfile, glmm, mac=4, res.savefn="saige_pval.rds")

	# fit the null model for quantitative outcomes
	glmm <- seqFitNullGLMM_SPA(yy ~ x1 + x2, pheno, gdsfile,
		trait.type="quantitative")
	# save model
	saveRDS(glmm, file="saige_model_quant.rds", compress="xz")
	# p-value calculation
	seqAssocGLMM_SPA(gdsfile, glmm, mac=4, res.savefn="saige_pval_quant.rds")
}


test.saige_fit_null_model <- function()
{
	# 'Rounding' was the default in versions prior to R_3.6.0
	# it is used for reproduction of the results created by R (v3.5.2)
	tryCatch(suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")),
		error=function(e) FALSE)

	# load the previous models
	mod1 <- readRDS(system.file("unitTests", "saige_model.rds",
		package="SAIGEgds"))
	mod2 <- readRDS(system.file("unitTests", "saige_model_quant.rds",
		package="SAIGEgds"))

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model, check binary outcomes
	glmm <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile)
	checkEquals(mod1, glmm, "check the SAIGE parameters (binary outcomes)",
		tolerance=1e-4)

	# fit the null model, check quantitative outcomes
	glmm <- seqFitNullGLMM_SPA(yy ~ x1 + x2, pheno, gdsfile,
		trait.type="quantitative")
	checkEquals(mod2, glmm,
		"check the SAIGE parameters (quantitative outcomes)", tolerance=1e-4)
}


test.saige_pval <- function()
{
	# load the previous models and results
	mod1 <- readRDS(system.file("unitTests", "saige_model.rds",
		package="SAIGEgds"))
	mod2 <- readRDS(system.file("unitTests", "saige_model_quant.rds",
		package="SAIGEgds"))
	pval1 <- readRDS(system.file("unitTests", "saige_pval.rds",
		package="SAIGEgds"))
	pval2 <- readRDS(system.file("unitTests", "saige_pval_quant.rds",
		package="SAIGEgds"))

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# p-value calculation, check binary outcomes
	assoc <- seqAssocGLMM_SPA(gdsfile, mod1, mac=4)
	checkEquals(pval1, assoc,
		"check the SAIGE p-value output (binary outcomes)", tolerance=1e-7)

	# p-value calculation, check quantitative outcomes
	assoc <- seqAssocGLMM_SPA(gdsfile, mod2, mac=4)
	checkEquals(pval2, assoc,
		"check the SAIGE p-value output (quantitative outcomes)",
		tolerance=1e-7)
}


test.saige_acta_o <- function()
{
	# load the prefit model
	mod <- readRDS(system.file("unitTests", "saige_model.rds",
		package="SAIGEgds"))

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# get a list of variant units for aggregate tests
	ut <- seqUnitSlidingWindows(gdsfile, win.size=200, win.shift=100)

	# run burden, ACAT-V & ACAT-O
	o <- seqAssocGLMM_spaACAT_O(gdsfile, mod, ut)
	v <- seqAssocGLMM_spaACAT_V(gdsfile, mod, ut)
	b <- seqAssocGLMM_spaBurden(gdsfile, mod, ut)

	# check p-value
	checkEquals(o$pval.b1_1, b$pval.b1_1, "ACAT-O vs Burden, beta(1,1)")
	checkEquals(o$pval.b1_25, b$pval.b1_25, "ACAT-O vs Burden, beta(1,25)")
	checkEquals(o$pval.v1_1, v$pval.v1_1, "ACAT-O vs ACAT-V, beta(1,1)")
	checkEquals(o$pval.v1_25, v$pval.v1_25, "ACAT-O vs ACAT-V, beta(1,25)")
}


test.pACAT <- function()
{
    # R implementation
    ps <- 10^-seq(1, 15, 0.1)
    A1 <- matrix(0, nrow=length(ps), ncol=length(ps))
    for (i in seq_along(ps))
    {
        for (j in seq_along(ps))
        {
            T <- mean(c(tanpi(0.5 - ps[i]), tanpi(0.5 - ps[j])))
            A1[i, j] <- 0.5 - atan(T)/pi
        }
    }

    # C implementation
    A2 <- matrix(0, nrow=length(ps), ncol=length(ps))
    for (i in seq_along(ps))
        for (j in seq_along(ps))
            A2[i, j] <- pACAT(c(ps[i], ps[j]))

    # check
	checkEquals(A1, A2, "R / C implementation of ACAT")
}


test.saige_gpu_fallback <- function()
{
	# This test verifies that use.gpu=TRUE works without error
	# even when no GPU is available (graceful fallback to CPU).

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype
	phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)

	# fit the null model with use.gpu=TRUE (should fall back to CPU gracefully)
	glmm_cpu <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile,
		use.gpu=FALSE, verbose=FALSE)
	glmm_gpu <- seqFitNullGLMM_SPA(y ~ x1 + x2, pheno, gdsfile,
		use.gpu=TRUE, verbose=FALSE)

	# results should match (if GPU was used, within tolerance; if CPU fallback, exact)
	checkEquals(glmm_cpu$tau, glmm_gpu$tau,
		"GPU fallback: tau should match CPU", tolerance=1e-6)
	checkEquals(glmm_cpu$coefficients, glmm_gpu$coefficients,
		"GPU fallback: coefficients should match CPU", tolerance=1e-6)
}


test.saige_survival <- function()
{
	tryCatch(suppressWarnings(RNGkind("Mersenne-Twister", "Inversion",
		"Rounding")), error=function(e) FALSE)

	# open a GDS file
	fn <- system.file("extdata", "grm1k_10k_snp.gds", package="SAIGEgds")
	gdsfile <- seqOpen(fn)
	on.exit(seqClose(gdsfile))

	# load phenotype and simulate a survival outcome (null wrt the SNPs)
	phenofn <- system.file("extdata", "pheno.txt.gz", package="SAIGEgds")
	pheno <- read.table(phenofn, header=TRUE, as.is=TRUE)
	set.seed(100)
	n <- nrow(pheno)
	eta <- 0.3*pheno$x1 - 0.2*pheno$x2
	ftime <- rexp(n, rate=exp(eta)*0.05)
	ctime <- rexp(n, rate=0.03)
	pheno$status <- as.integer(ftime <= ctime)
	pheno$atime  <- pmin(ftime, ctime)

	# fit the survival null model
	glmm <- seqFitNullGLMM_SPA(status ~ x1 + x2, pheno, gdsfile,
		trait.type="survival", event.time="atime", verbose=FALSE)
	checkEquals("survival", glmm$trait.type, "survival trait type")
	checkTrue(glmm$converged, "survival null model converged")
	# martingale residuals sum to ~0
	checkEqualsNumeric(0, sum(glmm$residuals), "martingale residuals sum to 0",
		tolerance=1e-6)
	checkTrue(all(glmm$fitted.values > 0), "all Poisson means are positive")

	# fixed effects should be close to coxph (Breslow ties); the frailty model
	# coefficients differ slightly when the estimated frailty variance > 0
	if (requireNamespace("survival", quietly=TRUE))
	{
		cx <- survival::coxph(survival::Surv(pheno$atime, pheno$status) ~
			pheno$x1 + pheno$x2, ties="breslow")
		checkEqualsNumeric(unname(coef(cx)), unname(glmm$coefficients),
			"survival fixed effects close to coxph(Breslow)", tolerance=0.05)
	}

	# single-variant association with the Poisson saddlepoint approximation
	assoc <- seqAssocGLMM_SPA(gdsfile, glmm, mac=4, verbose=FALSE)
	checkTrue(all(c("beta","SE","pval","p.norm","converged") %in%
		colnames(assoc)), "survival output columns")
	p <- assoc$pval[is.finite(assoc$pval) & assoc$pval > 0]
	checkTrue(all(p >= 0 & p <= 1), "p-values in [0,1]")
	# the SPA should engage for some variants
	checkTrue(any(assoc$method == "SPA"), "Poisson SPA engaged")
	# well calibrated under the null (genomic inflation near 1)
	lambda <- median(qchisq(p, 1, lower.tail=FALSE)) / qchisq(0.5, 1)
	checkTrue(lambda > 0.85 && lambda < 1.15,
		"genomic inflation factor near 1")

	# aggregate tests: burden & ACAT-V support survival, SKAT & ACAT-O reject
	units <- SeqArray::seqUnitSlidingWindows(gdsfile, win.size=5000,
		win.shift=5000)
	bd <- seqAssocGLMM_Burden(gdsfile, glmm, units, verbose=FALSE)
	checkTrue(all(bd$pval[is.finite(bd$pval)] >= 0 &
		bd$pval[is.finite(bd$pval)] <= 1), "survival burden p-values in [0,1]")
	checkTrue(all(c("p.norm","converged") %in% colnames(bd)),
		"survival burden output columns")
	av <- seqAssocGLMM_ACAT_V(gdsfile, glmm, units, verbose=FALSE)
	checkTrue(all(av$pval[is.finite(av$pval)] >= 0 &
		av$pval[is.finite(av$pval)] <= 1), "survival ACAT-V p-values in [0,1]")
	checkException(seqAssocGLMM_SKAT(gdsfile, glmm, units, verbose=FALSE),
		"SKAT should reject survival", silent=TRUE)
	checkException(seqAssocGLMM_ACAT_O(gdsfile, glmm, units, verbose=FALSE),
		"ACAT-O should reject survival", silent=TRUE)

	# refit the survival null model (event.time must be re-supplied)
	glmm2 <- seqFitNullGLMM_SPA(status ~ x1 + x2, pheno, gdsfile,
		trait.type="survival", event.time="atime", save.packed.geno=TRUE,
		verbose=FALSE)
	rf <- seqRefitNullGLMM(status ~ x1 + x2, pheno, glmm2, event.time="atime",
		verbose=FALSE)
	checkEquals("survival", rf$trait.type, "refit survival trait type")
	checkTrue(rf$converged, "refit survival converged")
	# with tau fixed (default), refit reproduces the original coefficients to
	# the inner IRLS tolerance (tol_coef = 1e-4)
	checkEqualsNumeric(unname(glmm2$coefficients), unname(rf$coefficients),
		"refit reproduces survival coefficients", tolerance=1e-3)
	# event.time is required for a survival refit
	checkException(seqRefitNullGLMM(status ~ x1 + x2, pheno, glmm2,
		verbose=FALSE), "refit needs event.time", silent=TRUE)
}
