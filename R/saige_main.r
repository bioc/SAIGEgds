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

SIMD <- function() .Call(saige_simd_version)


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


##suggested by Shawn 01-19-2018
Covariate_Transform <- function(formula, data)
{
  X1<-model.matrix(formula,data=data)
#  X1=X1[,c(2:ncol(X1))] #remove intercept
  formula.frame<-model.frame(formula,data=data)
  Y = model.response(formula.frame, type = "any")
  X_name = colnames(X1)
		
  # First run linear regression to identify multi collinearity 
  out.lm<-lm(Y ~ X1 - 1, data=data)
#  out.lm<-lm(Y ~ X1, data=data)
  idx.na<-which(is.na(out.lm$coef))
  if(length(idx.na)> 0){
	X1<-X1[, -idx.na]
	X_name = X_name[-idx.na]		
        cat("Warning: multi collinearity is detected in covariates! ", X_name[idx.na], " will be excluded in the model\n")
  }
  if(!(1 %in% idx.na)){
    X_name[1] = "minus1"
  }

	
 # QR decomposition
  Xqr = qr(X1)
  X1_Q = qr.Q(Xqr)
  qrr = qr.R(Xqr)
	
  N<-nrow(X1)
	
  # Make square summation=N (so mean=1)
  X1_new<-X1_Q * sqrt(N)	
  Param.transform<-list(qrr=qrr, N=N, X_name = X_name, idx.na=idx.na)
  re<-list(Y =Y, X1 = X1_new, Param.transform=Param.transform)
}


# In case to recover original scale coefficients
# X \beta = Q R \beta = (Q \sqrt(N)) ( R \beta / \sqrt(N))
# So coefficient from fit.new is the same as R \beta / \sqrt(N)
Covariate_Transform_Back<-function(coef, Param.transform)
{	
	#coef<-fit.new$coef; Param.transform=out.transform$Param.transform
	coef1<-coef * sqrt(Param.transform$N)
	coef.org<-solve(Param.transform$qrr, coef1)
	
	names(coef.org)<-Param.transform$X_name
	return(coef.org)
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

    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    n_samp <- dm[2L]
    n_var  <- dm[3L]
    if (verbose)
    {
        cat("Fit the null model: ")
        print(formula)
        cat("    # of samples: ", .pretty(n_samp), "\n", sep="")
        cat("    # of variants: ", .pretty(n_var), "\n", sep="")
    }


	out.transform<-Covariate_Transform(formula, data=data)
	formulaNewList = c("Y ~ ", out.transform$Param.transform$X_name[1])
	if(length(out.transform$Param.transform$X_name) > 1){
	for(i in c(2:length(out.transform$Param.transform$X_name))){
	  formulaNewList = c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
	}
	}
	formulaNewList = paste0(formulaNewList, collapse="")
	formulaNewList = paste0(formulaNewList, "-1")
	formula.new = as.formula(paste0(formulaNewList, collapse=""))
	data.new = data.frame(cbind(out.transform$Y, out.transform$X1))
	colnames(data.new) = c("Y",out.transform$Param.transform$X_name)
	cat("colnames(data.new) is ", colnames(data.new), "\n")
	cat("out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr), "\n")


    # 2-bit packed genotypes
    if (verbose) cat("Start loading genotypes")
    PackedGeno <- SeqArray:::.seqGet2bGeno(gdsfile)
    if (verbose)
    {
        cat(" [done] ")
        print(object.size(PackedGeno))
    }

    # initialize internal variables and buffers
    buf_geno <- double(4*n_var)
    buf_sigma <- double(n_samp)
    .Call(saige_store_geno, PackedGeno, n_samp, buf_geno, buf_sigma)

    # fit the model
    if (trait.type == "binary")
    {
        # binary outcome
        cat(phenovar, " is a binary trait:\n")

        fit0 <- glm(formula.new, data=data.new, family=binomial)
        if (verbose) print(fit0)
        obj.noK <- SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula.new, data.new)

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
        tau <- fixtau <- c(0,0)
        if (family$family %in% c("binomial", "poisson"))
            tau[1] = fixtau[1] = 1

        # change, use 0.5 as a default value, and use Get_Coef before getAIScore
        q = 1

        if (tau.init[fixtau==0] == 0)
            tau[fixtau==0] = 0.5
        else
            tau[fixtau==0] = tau.init[fixtau==0]
        cat("initial tau is ", paste(tau, collapse=", "), "\n", sep="")
        tau0 <- tau

        re.coef <- .Call(SAIGEgds:::saige_get_coeff, y, X, tau, family, alpha0, eta0,
            offset, maxiterPCG, maxiter, tolPCG, verbose)
        re <- .Call(SAIGEgds:::saige_get_AI_score, re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY,
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
            re.coef <- .Call(saige_get_coeff, y, X, tau, family, alpha0, eta0,
                offset, maxiterPCG, maxiter, tolPCG, verbose)
            fit <- fitglmmaiRPCG(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY,
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

        re.coef <- .Call(saige_get_coeff, y, X, tau, family, alpha0, eta0,  offset,
            maxiterPCG, maxiter, tolPCG, verbose)
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
    dsnode="", parallel=FALSE, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(mac), length(mac)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(!is.na(dsnode))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check model
    if (is.character(modobj))
    {
        stopifnot(length(modobj)==1L)
        modobj <- get(load(modobj))
    }
    .check_saige_model(modobj)

    # GDS file
    if (is.character(gdsfile))
    {
        gdsfile <- seqOpen(gdsfile)
        if (verbose)
            cat("Open '", gdsfile, "'\n", sep="")
        on.exit(seqClose(gdsfile))    
    }

    # determine the GDS node for dosages
    if (dsnode == "")
    {
        n <- index.gdsn(gdsfile, "genotype/data", silent=TRUE)
        if (!is.null(n))
        {
            dsnode <- "genotype"
        } else {
            n <- index.gdsn(gdsfile, "annotation/format/DS", silent=TRUE)
            if (!is.null(n))
            {
                dsnode <- "annotation/format/DS"
            } else {
                stop("Dosages should be stored in genotype or annotation/format/DS.")
            }
        }
    }

    if (verbose)
        cat("SAIGE association analyses:\n")

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
    mobj$S_a <- colSums(X1 * mobj$y_mu)    # a vector of size K

    if (!is.finite(mobj$var.ratio))
        stop("Invalid variance ratio in the SAIGE model.")
    # initialize internally
    .Call(saige_score_test_init, mobj)

    # scan all (selected) variants
    if (modobj$trait.type == "binary")
    {
        rv <- seqApply(gdsfile, dsnode, .cfunction("saige_score_test_bin"),
            as.is="list", parallel=parallel, .progress=verbose,
            .list_dup=FALSE)
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
        cat("# of variants after filtering MAF/MAC threshold: ",
            .pretty(length(rv)), "\n", sep="")
    }

    # output
    data.frame(
        id  = seqGetData(f, "variant.id"),
        chr = seqGetData(f, "chromosome"),
        pos = seqGetData(f, "position"),
        ref = seqGetData(f, "$ref"),
        alt   = seqGetData(f, "$alt"),
        AF.alt = sapply(rv, `[`, i=1L),
        AC.alt = sapply(rv, `[`, i=2L),
        num  = as.integer(sapply(rv, `[`, i=3L)),
        beta = sapply(rv, `[`, i=4L),
        SE   = sapply(rv, `[`, i=5L),
        pval = sapply(rv, `[`, i=6L),
        pval.noadj = sapply(rv, `[`, i=7L),
        converged = sapply(rv, `[`, i=8L)==1,
        stringsAsFactors = FALSE
    )
}
