#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2019        Xiuwen Zheng
# License: GPL-3
# Email: xiuwen.zheng@abbvie.com
#


#######################################################################
# Internal functions
#
.cfunction0 <- function(name)
{
    fn <- function() NULL
    f <- quote(.Call(SEQ_ExternalName0))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.cfunction <- function(name)
{
    fn <- function(x) NULL
    f <- quote(.Call(SEQ_ExternalName1, x))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

.version <- function() .Call(saige_simd_version)


# Internal model checking
.check_saige_model <- function(obj)
{
    stopifnot(is.list(obj))
    for (nm in c("sample.id", "trait.type", "var.ratio"))
    {
        if (!(nm %in% names(obj)))
            stop("'", nm, "' should be stored in the SAIGE model.")
    }
    if (!(obj$trait.type %in% c("binary")))
        stop("'trait.type' should be binary or .")
    invisible()
}


getCoefficients <- function(Yvec, Xmat, wVec, tauVec, maxiterPCG, tolPCG)
{
    .Call('_SAIGE_getCoefficients', Yvec, Xmat, wVec, tauVec, maxiterPCG,
        tolPCG, PACKAGE='SAIGEgds')
}



# Functon to get working vector and fixed & random coefficients
# Run iterations to get converged alpha and eta
Get_Coef <- function(y, X, tau, family, alpha0, eta0, offset, maxiterPCG,
    tolPCG, maxiter, verbose=FALSE)
{
    tol.coef <- 0.1
    mu <- family$linkinv(eta0)
    mu.eta <- family$mu.eta(eta0)
    Y <- eta0 - offset + (y - mu)/mu.eta
    sqrtW <- mu.eta / sqrt(family$variance(mu))
    W <- sqrtW^2

    for(i in 1:maxiter)
    {
        cat("iGet_Coef: ", i, "\n", sep="")
        re.coef <- getCoefficients(Y, X, W, tau, maxiterPCG, tol=tolPCG)
        alpha <- re.coef$alpha
        eta <- re.coef$eta + offset
        if (verbose)
        {
            cat("Tau: "); print(tau)
            cat("Fixed-effect coefficients: "); print(alpha)
        }
        mu <- family$linkinv(eta)
        mu.eta <- family$mu.eta(eta)
        Y <- eta - offset + (y - mu)/mu.eta
        sqrtW <- mu.eta/sqrt(family$variance(mu))
        W <- sqrtW^2
        if (max(abs(alpha-alpha0)/(abs(alpha)+abs(alpha0)+tol.coef)) < tol.coef)
            break
        alpha0 = alpha
    }

    list(Y=Y, alpha=alpha, eta=eta, W=W, cov=re.coef$cov, sqrtW=sqrtW,
        Sigma_iY=re.coef$Sigma_iY, Sigma_iX=re.coef$Sigma_iX, mu=mu)
}



#######################################################################
# Fit the null model
#

seqFitNullGLMM_SPA <- function(formula, data, gdsfile,
    trait.type=c("binary", "quantitative"), invNormalize=FALSE,
    maf=0.01, missing.rate=0.01, max.num.snp=100000L, variant.id=NULL,
    tol=0.02, maxiter=20L, tolPCG=1e-5, maxiterPCG=500L,
    nThreads=1, Cutoff=2,  numMarkers=30, tau.init = c(0,0),
    traceCVcutoff=1, ratioCVcutoff=1, model.save.fn=NA_character_,
    verbose=TRUE)
{
    stopifnot(inherits(formula, "formula"))
    stopifnot(is.data.frame(data))
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    trait.type <- match.arg(trait.type)
    stopifnot(is.character(model.save.fn), length(model.save.fn)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile)
        on.exit(seqClose(gdsfile))
    }

    # variables in the formula
    vars <- all.vars(formula)
    phenovar <- all.vars(formula)[1L]

    # check sample id
    if ("sample.id" %in% vars)
        stop("'sample.id' should not be in the formula.")
    if (!("sample.id" %in% colnames(data)))
        stop("'sample.id' should be one of the columns in 'data'.")
    if (is.factor(data$sample.id))
        stop("'sample.id' should not be a factor variable.")
    if (any(is.na(data$sample.id)))
        stop("'sample.id' should not have any missing value.")
    if (anyDuplicated(data$sample.id))
        stop("'sample.id' in data should be unique.")

    # remove missing values
    data <- data[, c("sample.id", vars)]
    data <- na.omit(data)
    seqResetFilter(gdsfile, sample=TRUE, verbose=FALSE)
    sid <- seqGetData(gdsfile, "sample.id")
    i <- match(sid, data$sample.id)
    i <- i[!is.na(i)]
    data <- data[i, ]
    if (nrow(data) <= 0L)
        stop("No common sample.id between 'data' and the GDS file.")
    seqSetFilter(gdsfile, sample.id=data$sample.id, verbose=FALSE)

    if (is.null(variant.id))
    {
        # filters of maf, mac, missing.rate
        seqSetFilter(gdsfile, sample.id=data$sample.id, verbose=FALSE)
        seqSetFilterCond(gdsfile, maf=maf, missing.rate=missing.rate, verbose=FALSE)
    } else {
        setSetFilter(gdsfile, variant.id=variant.id, verbose=FALSE)
    }

    if (verbose)
    {
        cat("Fit the null model: ")
        print(formula)
        dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)
        cat("    # of samples: ", .pretty(dm$seldim[2L]), "\n", sep="")
        cat("    # of variants: ", .pretty(dm$seldim[3L]), "\n", sep="")
    }

    # 2-bit packed genotypes
    if (verbose) cat("Start reading genotypes")
    PackedGeno <- SeqArray:::.seqGet2bGeno(gdsfile)
    if (verbose) cat(" [done].\n")

    # fit the model
    if (trait.type == "binary")
    {
        # binary outcome
        cat(phenovar, " is a binary trait:\n")

        fit0 <- glm(formula, data=data, family=binomial)
        if (verbose) print(fit0)
        obj.noK <- SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula, data)

        y <- fit0$y
        n <- length(y)
        X <- model.matrix(fit0)
        offset = fit0$offset
        if (is.null(offset)) offset <- rep(0, n)

        family <- fit0$family
        eta <- fit0$linear.predictors
        mu <- fit0$fitted.values
        mu.eta <- family$mu.eta(eta)
        Y <- eta - offset + (y - mu)/mu.eta
        alpha0 <- fit0$coef
        eta0 <- eta
        if (family$family %in% c("binomial", "poisson"))
            tau[1] = fixtau[1] = 1

        # change, use 0.5 as a default value, and use Get_Coef before getAIScore
        q = 1

        if (tauInit[fixtau==0] == 0)
            tau[fixtau==0] = 0.5
        else
            tau[fixtau==0] = tauInit[fixtau==0]
        cat("Inital tau is ", tau, "\n", sep="")
        tau0 <- tau

        re.coef <- Get_Coef(y, X, tau, family, alpha0, eta0,  offset,
            maxiterPCG, tolPCG, maxiter, verbose)
        re <- getAIScore(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY,
            re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG,
            traceCVcutoff)

        tau[2] <- max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)
        if(verbose)
        {
            cat("Variance component estimates: ")
            print(tau)
        }

        for (i in 1:maxiter)
        {
            if (verbose) cat("Iteration", i, ":", tau, "\n")
            alpha0 <- re.coef$alpha
            tau0 = tau
            eta0 = eta
            re.coef <- Get_Coef(y, X, tau, family, alpha0, eta0,  offset,
                maxiterPCG, tolPCG, maxiter, verbose)
            fit = fitglmmaiRPCG(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY,
                re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG,
                tol, traceCVcutoff)

            tau <- as.numeric(fit$tau)
            cov <- re.coef$cov
            alpha <- re.coef$alpha
            eta <- re.coef$eta
            Y <- re.coef$Y
            mu <- re.coef$mu

            if (tau[2] == 0) break
            if (max(abs(tau-tau0)/(abs(tau)+abs(tau0)+tol)) < tol) break
            if (max(tau) > tol^(-2))
            {
                warning("Large variance estimate observed in the iterations, model not converged ...",
                    call.=FALSE, immediate.=TRUE)
                i = maxiter
                break
            }
        }
        if (verbose) cat("Final:" ,tau, "\n")

        re.coef <- Get_Coef(y, X, tau, family, alpha0, eta0,  offset,
            maxiterPCG, tolPCG, maxiter, verbose)
        cov <- re.coef$cov
        alpha <- re.coef$alpha
        eta <- re.coef$eta
        Y <- re.coef$Y
        mu <- re.coef$mu
        converged <- ifelse(i < maxiter, TRUE, FALSE)
        res <- y - mu
        # coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
        coef.alpha <- alpha

        ans <- list(theta=tau, coefficients=coef.alpha, linear.predictors=eta,
  	        fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,
  	        obj.glm.null=fit0, obj.noK=obj.noK, traitType=trait.type,
  	        sample.id=data$sample.id)
    }

    ans
}




#######################################################################
# SAIGE single-variant analysis
#

seqAssocGLMM_SPA <- function(gdsfile, modobj, maf=NaN, mac=NaN,
    parallel=FALSE, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(mac), length(mac)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
        cat("SAIGE association analyses:\n")

    # check model
    if (is.character(modobj))
    {
        stopifnot(length(modobj)==1L)
        modobj <- get(load(modobj))
    }
    .check_saige_model(modobj)

    # save the current filter
    seqSetFilter(gdsfile, action="push", verbose=FALSE)
    on.exit({ seqSetFilter(gdsfile, action="pop", verbose=FALSE) })

    # check sample ID
    seqSetFilter(gdsfile, sample.id=modobj$sample.id, verbose=FALSE)
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs are not available in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")

    if (verbose)
    {
        dm <- SeqArray:::.seldim(gdsfile)
        cat("    # of samples: ", .pretty(dm[2L]), "\n", sep="")
        cat("    # of variants: ", .pretty(dm[3L]), "\n", sep="")
    }

    # initialize the internal model
    y <- unname(modobj$obj.noK$y)
    mu <- unname(modobj$fitted.values)
    X1 <- modobj$obj.noK$X1[ii, ]
    n <- length(ii)
	mobj <- list(
	    maf = maf, mac = mac,
        y = y[ii], mu = mu[ii],
        y_mu = y[ii] - mu[ii],  # y - mu
        mu2 = (mu * (1 - mu))[ii],
        t_XXVX_inv = t(modobj$obj.noK$XXVX_inv[ii, ]),  # K x n_samp (K << n_samp, more efficient)
        XV = modobj$obj.noK$XV[, ii],  # K x n_samp
        t_XVX_inv_XV = t(modobj$obj.noK$XXVX_inv[ii, ] * modobj$obj.noK$V[ii]),  # K x n_samp
        t_X = t(X1),  # K x n_samp
        var.ratio = mean(modobj$var.ratio, na.rm=TRUE),
        buf_coeff = double(nrow(modobj$obj.noK$XV)),
        buf_adj_g = double(n),
        buf_index = integer(n),
        buf_B = double(n),
        buf_g_tilde = double(n),
        buf_tmp = double(ncol(X1))
	)
    mobj$XVX <- t(X1) %*% (X1 * mobj$mu2)  # a matrix: K x K
    mobj$S_a <- colSums(X1 * mobj$y_mu)  # a vector of size K

    if (!is.finite(mobj$var.ratio))
        stop("Invalid variance ratio in the SAIGE model.")
    # initialize internally
    .Call(saige_score_test_init, mobj)

    # scan all (selected) variants
    if (modobj$trait.type == "binary")
    {
        rv <- seqApply(gdsfile, "annotation/format/DS",
            .cfunction("saige_score_test_bin"), as.is="list",
            parallel=parallel, .progress=verbose, .list_dup=FALSE)
    }

    # if any maf/mac filter
    x <- sapply(rv, is.null)
    if (any(x))
    {
        x <- !x
        seqSetFilter(gdsfile, variant.sel=x, action="intersect", verbose=T)
        rv <- rv[x]
    }
    if (verbose)
    {
        cat("# of variants after filtering MAF/MAC: ", .pretty(length(rv)),
            "\n", sep="")
    }

    # output
    data.frame(
        id  = seqGetData(f, "variant.id"),
        chr = seqGetData(f, "chromosome"),
        pos = seqGetData(f, "position"),
        ref = seqGetData(f, "$ref"),
        alt   = seqGetData(f, "$alt"),
        AF   = sapply(rv, `[`, i=1L),
        AC   = sapply(rv, `[`, i=2L),
        num  = as.integer(sapply(rv, `[`, i=3L)),
        beta = sapply(rv, `[`, i=4L),
        SE   = sapply(rv, `[`, i=5L),
        pval = sapply(rv, `[`, i=6L),
        pval.noadj = sapply(rv, `[`, i=7L),
        SPA.converged = sapply(rv, `[`, i=8L)==1,
        stringsAsFactors = FALSE
    )
}
