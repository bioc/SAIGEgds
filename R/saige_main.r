#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2019        Xiuwen Zheng (xiuwen.zheng@abbvie.com)
# License: GPL-3
#


#######################################################################
# Internal functions
#

.cfunction <- function(name)
{
    fn <- function(x) NULL
    f <- quote(.Call(SEQ_ExternalName1, x))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.cfunction2 <- function(name)
{
    fn <- function(x, y) NULL
    f <- quote(.Call(SEQ_ExternalName2, x, y))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

SIMD <- function() .Call(saige_simd_version)

.crayon_inverse <- function(s)
{
    if (requireNamespace("crayon", quietly=TRUE))
        s <- crayon::inverse(s)
    s
}

.crayon_underline <- function(s)
{
    if (requireNamespace("crayon", quietly=TRUE))
        s <- crayon::underline(s)
    s
}


# Internal model checking
.check_saige_model <- function(obj)
{
    stopifnot(is.list(obj))
    for (nm in c("sample.id", "trait.type", "var.ratio"))
    {
        if (!(nm %in% names(obj)))
            stop("'", nm, "' should be stored in the SAIGE model.")
    }
    if (!(obj$trait.type %in% c("binary", "quantitative")))
        stop("'trait.type' should be binary or quantitative.")
    invisible()
}



#######################################################################
# Fit the null model
#

seqFitNullGLMM_SPA <- function(formula, data, gdsfile,
    trait.type=c("binary", "quantitative"), maf=0.01, missing.rate=0.01,
    max.num.snp=1000000L, variant.id=NULL, inv.norm=TRUE, X.transform=TRUE,
    tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5, maxiterPCG=500L,
    num.marker=30L, tau.init=c(0,0), traceCVcutoff=1, ratioCVcutoff=1,
    geno.sparse=TRUE, num.thread=1L, model.savefn="", seed=200L, verbose=TRUE)
{
    stopifnot(inherits(formula, "formula"))
    stopifnot(is.data.frame(data))
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    trait.type <- match.arg(trait.type)
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.numeric(max.num.snp), length(max.num.snp)==1L)
    stopifnot(is.logical(inv.norm), length(inv.norm)==1L)
    stopifnot(is.logical(X.transform), length(X.transform)==1L)
    stopifnot(is.numeric(tol), length(tol)==1L)
    stopifnot(is.numeric(maxiter), length(maxiter)==1L)
    stopifnot(is.numeric(nrun), length(nrun)==1L)
    stopifnot(is.numeric(tolPCG), length(tolPCG)==1L)
    stopifnot(is.numeric(maxiterPCG), length(maxiterPCG)==1L)
    stopifnot(is.numeric(num.marker), length(num.marker)==1L)
    stopifnot(is.numeric(tau.init), length(tau.init)==2L)
    stopifnot(is.numeric(traceCVcutoff), length(traceCVcutoff)==1L)
    stopifnot(is.numeric(ratioCVcutoff), length(ratioCVcutoff)==1L)
    stopifnot(is.logical(geno.sparse), length(geno.sparse)==1L)
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.character(model.savefn), length(model.savefn)==1L)
    stopifnot(is.numeric(seed), length(seed)==1L, is.finite(seed))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
        cat(.crayon_inverse("SAIGE association analysis:\n"))
    if (is.character(gdsfile))
    {
        if (verbose)
            cat("Open the genotype file '", gdsfile, "'\n", sep="")
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
        if (verbose)
            cat("Filtering variants:\n")
        seqSetFilterCond(gdsfile, maf=maf, missing.rate=missing.rate,
            parallel=num.thread, .progress=TRUE, verbose=FALSE)
    } else {
        seqSetFilter(gdsfile, variant.id=variant.id, verbose=FALSE)
    }

    # the max mumber of SNPs
    v <- seqGetFilter(gdsfile)$variant.sel
    n <- sum(v, na.rm=TRUE)
    if (n > max.num.snp)
    {
        set.seed(seed)
        seqSetFilter(gdsfile, variant.sel=sample(which(v), max.num.snp),
            verbose=FALSE)
    }

    # get the number of samples / variants
    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    n_samp <- dm[2L]
    n_var  <- dm[3L]
    if (verbose)
    {
        cat("Fit the null model:", format(formula), "+ var(GRM)\n")
        cat("    # of samples: ", .pretty(n_samp), "\n", sep="")
        cat("    # of variants:", .pretty(n_var))
        if (n > max.num.snp)
            cat(" (randomly selected from ", .pretty(n), ")", sep="")
        cat("\n")
    }

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    setThreadOptions(num.thread)
    if (verbose)
    {
        cat("    using ", num.thread, " thread",
            ifelse(num.thread>1L, "s", ""), "\n", sep="")
    }

    if (isTRUE(X.transform))
    {
        if (verbose)
            cat("Transform on the design matrix with QR decomposition:\n")
        X <- model.matrix(formula, data)
        frm <- model.frame(formula, data)
        y <- model.response(frm, type="any")
        # check multi-collinearity
        m <- lm(y ~ X - 1)
        i_na <- which(is.na(m$coefficients))
        if (length(i_na) > 0L)
        {
            X <- X[, -i_na]
            if (verbose)
            {
                cat("    exclude ", length(i_na), " covariates (",
                    paste(colnames(X)[i_na], collapse=", "),
                    ") to avoid multi collinearity.\n", sep="")
            }
        }
        X_name <- colnames(X)
        Xqr = qr(X)  # QR decomposition
        X_new <- qr.Q(Xqr) * sqrt(nrow(X))
        X_qrr <- qr.R(Xqr)
        data <- data.frame(cbind(y, X_new))
        nm <- paste0("x", seq_len(ncol(X_new))-1L)
        colnames(data) <- c("y", nm)
        formula <- as.formula(paste("y ~", paste(nm, collapse=" + "), "-1"))
        if (verbose)
            cat("    new formula: ", format(formula), "\n", sep="")
    }

    # load SNP genotypes
    if (verbose)
        cat("Start loading SNP genotypes:\n")
    if (isTRUE(geno.sparse))
    {
        # sparse genotypes
        buffer <- integer(n_samp + 4L)
        packed.geno <- seqApply(gdsfile, "$dosage_alt", .cfunction2("saige_get_sparse"),
            as.is="list", .useraw=TRUE, .list_dup=FALSE, y=buffer, .progress=verbose)
        rm(buffer)
    } else {
        # 2-bit packed genotypes
        packed.geno <- SeqArray:::.seqGet2bGeno(gdsfile, verbose)
    }
    if (verbose)
    {
        cat("    using ")
        cat(SeqArray:::.pretty_size(as.double(object.size(packed.geno))))
        cat(ifelse(isTRUE(geno.sparse), " (sparse matrix)\n", " (dense matrix)\n"))
    }

    # initialize internal variables and buffers
    buf_std_geno <- double(4L*n_var)
    buf_sigma <- double(n_samp)
    buf_crossprod <- matrix(0.0, nrow=n_samp, ncol=num.thread)
    if (isTRUE(geno.sparse))
    {
        .Call(saige_store_sp_geno, packed.geno, n_samp, buf_std_geno, buf_sigma,
            buf_crossprod)
    } else {
        .Call(saige_store_2b_geno, packed.geno, n_samp, buf_std_geno, buf_sigma,
            buf_crossprod)
    }

    # parameters for fitting the model
    param <- list(
        num.thread = num.thread,
        seed = seed,
        tol = tol, tolPCG = tolPCG,
        maxiter = maxiter, maxiterPCG = maxiterPCG,
        nrun = nrun,
        num.marker = num.marker,
        traceCVcutoff = traceCVcutoff,
        ratioCVcutoff = ratioCVcutoff,
        verbose = verbose
    )

    # fit the model
    if (trait.type == "binary")
    {
        # binary outcome
        cat("Binary outcome: ", phenovar, "\n", sep="")

        # fit the null model
        fit0 <- glm(formula, data=data, family=binomial)
        if (verbose)
        {
            cat("Initial fixed-effect coefficients:\n")
            v <- fit0$coefficients
            names(v) <- c(paste0("    ", names(v)[1L]), names(v)[-1L])
            print(v)
        }
        obj.noK <- SPAtest:::ScoreTest_wSaddleApprox_NULL_Model(formula, data)

        # initial tau
        tau <- fixtau <- c(0,0)
        if (fit0$family$family %in% c("binomial", "poisson"))
            tau[1] <- fixtau[1] <- 1
        if (tau.init[fixtau==0] == 0)
            tau[fixtau==0] <- 0.5
        else
            tau[fixtau==0] <- tau.init[fixtau==0]

        # iterate
        X <- model.matrix(fit0)
        glmm <- .Call(saige_fit_AI_PCG_binary, fit0, X, tau, param)

        # calculate the variance ratio
        if (verbose)
            cat(.crayon_underline("Calculate the average ratio of variances:\n"))
        set.seed(seed)
        var.ratio <- .Call(saige_calc_var_ratio_binary, fit0, glmm, obj.noK,
            param, sample.int(n_var, n_var))
        var.ratio <- var.ratio[order(var.ratio$id), ]
        var.ratio$id <- seqGetData(gdsfile, "variant.id")[var.ratio$id]
        rownames(var.ratio) <- NULL

        # ans$obj.glm.null <- fit0
        glmm$obj.noK <- obj.noK
        glmm$var.ratio <- var.ratio

    } else if (trait.type == "quantitative")    
    {
        # quantitative outcome
        cat("Quantitative outcome: ", phenovar, "\n", sep="")

        # fit the null model
        fit0 <- glm(formula, data=data, family=gaussian)
        if (verbose)
        {
            cat("Initial fixed-effect coefficients:\n")
            v <- fit0$coefficients
            names(v) <- c(paste0("    ", names(v)[1L]), names(v)[-1L])
            print(v)
        }
        obj.noK <- SPAtest:::ScoreTest_wSaddleApprox_NULL_Model_q(formula, data)

      system.time(modglmm<-glmmkin.ai_PCG_Rcpp_Quantitative(plinkFile,fit0, tau = c(0,0), fixtau = c(0,0), maxiter =maxiter, tol = tol, verbose = TRUE, nrun=30, tolPCG = tolPCG, maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, obj.noK=obj.noK, out.transform=out.transform, tauInit=tauInit, memoryChunk = memoryChunk, LOCO=LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec, chromosomeEndIndexVec = chromosomeEndIndexVec, traceCVcutoff = traceCVcutoff))
      save(modglmm, file = modelOut)
      print("step2")



    } else {
        stop("Invalid 'trait.type'.")    
    } 

    if (verbose)
    {
        cat("    ratio avg. is ", mean(var.ratio$ratio),
            ", sd: ", sd(var.ratio$ratio), "\n", sep="")
    }

    # tweak the result
    if (!isTRUE(X.transform))
    {
        names(glmm$coefficients) <- colnames(obj.noK$X1)
    } else {
        coef <- solve(X_qrr, glmm$coefficients * sqrt(nrow(data)))
        names(coef) <- X_name
        glmm$coefficients <- coef
    }
    names(glmm$tau) <- c("Sigma_E", "Sigma_G")
    glmm$trait.type <- trait.type
    glmm$sample.id <- seqGetData(gdsfile, "sample.id")
    glmm$variant.id <- seqGetData(gdsfile, "variant.id")

    if (!is.na(model.savefn) && model.savefn!="")
    {
        cat("Save the model to '", model.savefn, "'\n", sep="")
        save(glmm, file=model.savefn)
    }
    if (verbose)
        cat(.crayon_inverse("Done."), "\n", sep="")

    if (!is.na(model.savefn) && model.savefn!="")
        return(invisible(glmm))
    else
        return(glmm)
}




#######################################################################
# SAIGE single variant analysis
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
            dsnode <- "$dosage_alt"
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
        cat("SAIGE association analysis:\n")

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

    # initialize the internal model parameters
    y <- unname(modobj$obj.noK$y)
    mu <- unname(modobj$fitted.values)
    X1 <- modobj$obj.noK$X1[ii, ]
    n <- length(ii)
    mobj <- list(
        maf = maf, mac = mac,
        tau = modobj$tau,
        y = y[ii], mu = mu[ii],
        y_mu = y[ii] - mu[ii],  # y - mu
        mu2 = (mu * (1 - mu))[ii],
        t_XXVX_inv = t(modobj$obj.noK$XXVX_inv[ii, ]),  # K x n_samp (K << n_samp, more efficient)
        XV = modobj$obj.noK$XV[, ii],  # K x n_samp
        t_XVX_inv_XV = t(modobj$obj.noK$XXVX_inv[ii, ] * modobj$obj.noK$V[ii]),  # K x n_samp
        t_X = t(X1),  # K x n_samp
        var.ratio = mean(modobj$var.ratio$ratio, na.rm=TRUE),
        # buffer
        buf_dosage = double(n),
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
    } else if (modobj$trait.type == "quantitative")    
    {
        rv <- seqApply(gdsfile, dsnode, .cfunction("saige_score_test_quant"),
            as.is="list", parallel=parallel, .progress=verbose,
            .list_dup=FALSE)
    } else {
        stop("Invalid 'modobj$trait.type'.")
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
    ans <- data.frame(
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
        stringsAsFactors = FALSE
    )
    if (modobj$trait.type == "binary")
    {
        ans$pval.noadj <- sapply(rv, `[`, i=7L)
        ans$converged <- as.logical(sapply(rv, `[`, i=8L))
    }

    ans
}
