#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2019-2022    Xiuwen Zheng / AbbVie-ComputationalGenomics
# License: GPL-3
#


#######################################################################
# Fit the null model
#

# check grm.mat
.check_grm_mat <- function(grm.mat, verbose)
{
    # check class
    if (inherits(grm.mat, "snpgdsGRMClass"))
    {
        m <- grm.mat$grm
        colnames(m) <- rownames(m) <- grm.mat$sample.id
        grm.mat <- m
    } else if (is.character(grm.mat))
    {
        if (length(grm.mat) != 1L)
            stop("'grm.mat' should be a file name.")
        if (verbose)
        {
            cat("Load a user-defined GRM:\n")
            .cat("    File: ", grm.mat)
        }
        if (grepl("\\.rds$", grm.mat, ignore.case=TRUE))
        {
            grm.mat <- readRDS(grm.mat)
        } else if (grepl("\\.(rda|RData)$", grm.mat, ignore.case=TRUE))
        {
            grm.mat <- get(load(grm.mat))
        } else
            stop("Unknown format of ", sQuote(grm.mat))
        if (verbose)
        {
            cat("    ")
            if (inherits(grm.mat, "sparseMatrix"))
                cat("Sparse matrix: ")
            else if (is.matrix(grm.mat))
                cat("Dense matrix: ")
            print(object.size(grm.mat))
        }
    }
    if (is.matrix(grm.mat) || inherits(grm.mat, "sparseMatrix"))
    {
        if (nrow(grm.mat) != ncol(grm.mat))
            stop("grm.mat should a numeric square matrix.")
        if ((is.matrix(grm.mat) && !is.double(grm.mat)) ||
            (inherits(grm.mat, "sparseMatrix") && !is.double(grm.mat@x)))
        {
            stop("grm.mat should be a numeric matrix.")
        }
        if (is.null(colnames(grm.mat)) || is.null(rownames(grm.mat)))
            stop("The row and column names of 'grm.mat' should be sample IDs.")
        if (!identical(colnames(grm.mat), rownames(grm.mat)))
            stop("The row and column names of 'grm.mat' should be the same.")
        if (anyDuplicated(colnames(grm.mat)))
            stop("The sample IDs in 'grm.mat' should be unique.")
    } else if (!is.null(grm.mat) && !isTRUE(grm.mat))
        stop("grm.mat should be NULL, a matrix, a file name or TRUE.")
    # check
    d <- diag(grm.mat)  # diag(TRUE)=1
    if (anyNA(d) || any(d <= 0))
        stop("The diagonal of GRM should be all positive.")
    # output
    grm.mat
}


# coerce to dgCMatrix (column-oriented sparse form)
.sp_to_dgCMatrix <- function(m)
{
    stopifnot(inherits(m, "sparseMatrix"))
    # if RsparseMatrix ==> TsparseMatrix first
    if (inherits(m, "RsparseMatrix"))
        m <- as(m, "TsparseMatrix")
    if (inherits(m, "TsparseMatrix") && !inherits(m, "dgTMatrix"))
        m <- as(m, "dgTMatrix")
    if (!inherits(m, "dgCMatrix"))
        m <- as(m, "generalMatrix")
    stopifnot(is(m, "dgCMatrix"))
    m
}


# get diagnoal of the full GRM
.get_grm_diag <- function() .Call(saige_get_grm_diag)


# calculate the variance ratio
.calcVR <- function(gdsfile, seed, fit0, glmm, obj.noK, param, verbose)
{
    var.ratio <- 1
    if (!is.null(gdsfile))
    {
        if (verbose)
        {
            .cat(.crayon_underline(date()))
            .cat(.crayon_inverse("Calculate the average ratio of variances:"))
        }
        set.seed(seed)
        var.ratio <- .Call(saige_calc_var_ratio, fit0, glmm, obj.noK, param)
    }
    glmm$var.ratio <- var.ratio
    glmm
}


# binary outcome
.fit_binary <- function(verbose, X.transform, phenovar, data, formula, param,
    tau.init, gdsfile, grm.mat, seed, n_var)
{
    if (verbose)
    {
        .cat("Binary outcome: ", phenovar)
        if (isTRUE(X.transform))
            y <- data$y
        else
            y <- data[[phenovar]]
        v <- table(y)
        n <- length(v) 
        v <- data.frame(v, as.numeric(prop.table(v)))
        v[, 1L] <- paste0("      ", v[, 1L])
        colnames(v) <- c(phenovar, "Number", "Proportion")
        print(v, row.names=FALSE)
        if (n != 2L)
            stop("The outcome variable has more than 2 categories!")
    }

    # fit the null model
    fit0 <- glm(formula, data=data, family=binomial)
    if (verbose)
    {
        if (is.null(gdsfile) && is.null(grm.mat))
            cat("Fixed-effect coefficients:\n")
        else
            cat("Initial fixed-effect coefficients:\n")
        v <- as.data.frame(t(fit0$coefficients))
        rownames(v) <- "   "
        print(v, width=128L)
    }

    # design matrix
    X <- unname(model.matrix(fit0))
    attr(X, "assign") <- attr(X, "contrasts") <- NULL

    if (!is.null(gdsfile) || !is.null(grm.mat))
    {
        # initial tau
        tau <- fixtau <- c(0, 0)
        if (fit0$family$family %in% c("binomial", "poisson"))
            tau[1] <- fixtau[1] <- 1
        if (sum(tau.init[fixtau==0]) == 0)
            tau[fixtau==0] <- 0.1
        else
            tau[fixtau==0] <- tau.init[fixtau==0]
        # iterate
        glmm <- .Call(saige_fit_AI_PCG, fit0, X, tau, param)
    } else {
        # no random effect
        glmm <- list(
            coefficients = unname(fit0$coefficients),
            tau = c(1, 0),
            linear.predictors = unname(fit0$linear.predictors),
            fitted.values = unname(fit0$fitted.values),
            residuals = unname(fit0$residuals),
            converged = fit0$converged)
    }

    # use updated mu to set obj.noK
    mu <- glmm$fitted.values
    V <- mu*(1-mu)
    XV <- t(X * V)
    XVX_inv <- solve(crossprod(X, X * V))
    XXVX_inv <- X %*% XVX_inv
    obj.noK <- list(y=unname(fit0$y), mu=mu, V=V, X1=X, XV=XV, XXVX_inv=XXVX_inv)
    glmm$obj.noK <- obj.noK

    if (!is.null(grm.mat))
    {
        # get the inverse of the Sigma matrix
        sigma <- grm.mat
        if (!is.null(gdsfile))
            diag(sigma) <- .get_grm_diag()  # for better approximation
        sigma <- glmm$tau[2L] * sigma
        diag(sigma) <- diag(sigma) + 1/V
        colnames(sigma) <- rownames(sigma) <- NULL
        m <- chol2inv(chol(sigma))
        if (inherits(m, "sparseMatrix")) m <- as(m, "symmetricMatrix")
        glmm$Sigma_inv <- obj.noK$Sigma_inv <- m
        # get the part of projection matrix
        X1 <- obj.noK$X1
        s_X1 <- solve(sigma, X1)
        m <- as.matrix(solve(crossprod(X1, s_X1)))
        m <- t(as.matrix(tcrossprod(chol(m), s_X1)))
        dimnames(m) <- NULL
        glmm$chol_inv_X_Sigma <- m
    } else {
        obj.noK$Sigma_inv <- FALSE
    }

    # calculate the variance ratio
    .calcVR(gdsfile, seed, fit0, glmm, obj.noK, param, verbose)
}


# quantitative outcome
.fit_quant <- function(verbose, X.transform, phenovar, data, inv.norm, formula,
    param, tau.init, gdsfile, grm.mat, seed, n_var)
{
    if (verbose)
    {
        .cat("Quantitative outcome: ", phenovar)
        if (isTRUE(X.transform))
            y <- data$y
        else
            y <- data[[phenovar]]
        v <- data.frame(mean=mean(y), sd=sd(y), min=min(y), max=max(y))
        rownames(v) <- "   "
        print(v)
    }

    # inverse normal transformation
    if (inv.norm == "residuals")
    {
        if (isTRUE(X.transform)) phenovar <- "y"
        fit0 <- glm(formula, data=data)
        resid.sd <- sd(fit0$residuals)
        new.y <- .rank_norm(fit0$residuals, s=resid.sd)
        data[[phenovar]] <- new.y
        if (verbose)
        {
            cat("Inverse normal transformation on the residuals with:\n")
            .cat("    standard deviation = ", resid.sd)
        }
    } else if (inv.norm == "quant")
    {
        if (isTRUE(X.transform)) phenovar <- "y"
        y <- data[[phenovar]]
        new.y <- .rank_norm(y)
        data[[phenovar]] <- new.y
        if (verbose)
            cat("Inverse normal transformation on the outcome variable\n")
    }

    # fit the null model
    fit0 <- glm(formula, data=data)
    if (verbose)
    {
        if (is.null(gdsfile) && is.null(grm.mat))
            cat("Fixed-effect coefficients:\n")
        else
            cat("Initial fixed-effect coefficients:\n")
        v <- as.data.frame(t(fit0$coefficients))
        rownames(v) <- "   "
        print(v, width=128L)
    }

    # design matrix
    X <- unname(model.matrix(fit0))
    attr(X, "assign") <- attr(X, "contrasts") <- NULL

    if (!is.null(gdsfile) || !is.null(grm.mat))
    {
        y <- fit0$y
        offset <- fit0$offset
        if (is.null(offset)) offset <- rep(0, length(y))
        eta <- fit0$linear.predictors
        mu <- fit0$fitted.values
        mu.eta <- fit0$family$mu.eta(eta)
        Y <- eta - offset + (y - mu)/mu.eta
        # initial tau
        tau <- tau.init
        if (sum(tau) == 0) tau <- c(0.5, 0.5)
        tau <- var(Y) * tau / sum(tau)
        # iterate
        glmm <- .Call(saige_fit_AI_PCG, fit0, X, tau, param)
    } else {
        # no random effect
        glmm <- list(
            coefficients = unname(fit0$coefficients),
            tau = c(var(fit0$residuals), 0),
            linear.predictors = unname(fit0$linear.predictors),
            fitted.values = unname(fit0$fitted.values),
            residuals = unname(fit0$residuals),
            converged = fit0$converged)
        if (verbose)
            .cat("Estimated tau: (", glmm$tau[1L], ", ", glmm$tau[2L], ")")
    }

    obj.noK <- list(y=unname(fit0$y), mu=glmm$fitted.values,
        V=rep(1, length(fit0$y)),
        X1=X, XV=t(X), XXVX_inv=X %*% solve(crossprod(X)))
    glmm$obj.noK <- obj.noK

    if (!is.null(grm.mat))
    {
        # get the inverse of the Sigma matrix
        sigma <- grm.mat
        if (!is.null(gdsfile))
            diag(sigma) <- .get_grm_diag()  # for better approximation
        sigma <- glmm$tau[2L] * sigma
        diag(sigma) <- diag(sigma) + glmm$tau[1L]
        colnames(sigma) <- rownames(sigma) <- NULL
        m <- chol2inv(chol(sigma))
        if (inherits(m, "sparseMatrix")) m <- as(m, "symmetricMatrix")
        glmm$Sigma_inv <- obj.noK$Sigma_inv <- m
        # get the part of projection matrix
        s_X <- solve(sigma, X)
        m <- as.matrix(solve(crossprod(X, s_X)))
        m <- t(as.matrix(tcrossprod(chol(m), s_X)))
        dimnames(m) <- NULL
        glmm$chol_inv_X_Sigma <- m
    } else {
        obj.noK$Sigma_inv <- FALSE
    }

    # calculate the variance ratio
    .calcVR(gdsfile, seed, fit0, glmm, obj.noK, param, verbose)
}


# check use.cateMAC
.cateMAC_default <- c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)

.check_use_cateMAC <- function(use.cateMAC)
{
    err <- "'use.cateMAC' should be FALSE, TRUE or a numeric vector."
    if (is.logical(use.cateMAC))
    {
        if (length(use.cateMAC) != 1L) stop(err)
        if (is.na(use.cateMAC)) stop(err)
        if (use.cateMAC)
            use.cateMAC <- .cateMAC_default
    } else if (is.numeric(use.cateMAC))
    {
        if (!is.vector(use.cateMAC) || length(use.cateMAC)<=0L)
            stop(err)
        if (anyNA(use.cateMAC))
            stop("'use.cateMAC' should not contain NA or NaN.")
        if (is.unsorted(use.cateMAC, strictly=TRUE))
            stop("'use.cateMAC' should be strictly increasing.")
        if (use.cateMAC[1L] <= 0)
            stop("'use.cateMAC[1]' should be > 0.")
        if (!all(is.finite(use.cateMAC)))
            stop("The numeric values in 'use.cateMAC' should be finite.")
    }
    use.cateMAC
}

# show use.cateMAC
.show_use_cateMAC <- function(use.cateMAC, verbose)
{
    if (verbose)
    {
        cat("MAC categories: ")
        if (is.logical(use.cateMAC))
        {
            if (isFALSE(use.cateMAC))
                cat("none\n")
            else
                .cat(use.cateMAC)
        } else {
            x1 <- c(0, use.cateMAC)
            x2 <- c(use.cateMAC, Inf)
            s <- paste0("[", x1, ",", x2, ")")
            .cat(paste(s, collapse=", "))
        }
    }
    invisible()
}

# simulated genotype packed RAW matrix (mac_low <= MAC < mac_high)
.simu_geno <- function(nsamp, nsnp, mac_low, mac_high)
{
    lv <- ceiling(mac_low)
    hv <- floor(mac_high)
    if (hv == mac_high) hv <- hv - 1L
    if (lv <= hv)
    {
        gmat <- matrix(as.raw(0L), nrow=ceiling(nsamp/4L), ncol=nsnp)
        mac_int <- seq.int(lv, hv)
        for (i in seq_len(nsnp))
        {
            if (length(mac_int) > 1L)
                mac <- sample(mac_int, 1L)
            else
                mac <- mac_int
            j <- sample.int(2L*nsamp, mac)
            # set genotypes according to j
            .Call(saige_set_geno2b_raw, gmat, j, i)
        }
        gmat
    } else {
        stop(sprintf("No integer MAC in [%g, %g).", mac_low, mac_high))
    }
}

# simulated genotype packed RAW matrix (mac_low <= MAC < mac_high)
.get_sparse_geno <- function(gdsfile, nfork, verbose)
{
    # integer genotypes or numeric dosages
    if (exist.gdsn(gdsfile, "genotype/data"))
    {
        nm <- "$dosage_alt"
    } else if (exist.gdsn(gdsfile, "annotation/format/DS/data"))
    {
        nm <- "annotation/format/DS"
        if (verbose) cat("    using 'annotation/format/DS'\n")
    } else {
        stop("'genotype' and 'annotation/format/DS' are not available.")
    }
    # internal buffer
    n_samp <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim[2L]
    buf_b256 <- integer((ceiling(n_samp/256)+1)*3 + 1)
    buf_b1 <- raw(3*ceiling(n_samp/256)*257)
    .Call(saige_init_sparse, n_samp, buf_b256, buf_b1)
    fc <- .cfunction("saige_get_sparse")
    # run
    seqParallel(nfork, gdsfile, FUN=function(f) {
        seqApply(f, nm, fc, as.is="list", .useraw=TRUE, .list_dup=FALSE,
            .progress=nfork==1L && verbose)
    }, .balancing=TRUE, .bl_size=1000L, .bl_progress=verbose)
}

# get information from the output of .get_sparse_geno()
.get_sparse_smallblock_len <- function(x)
{
    stopifnot(is.raw(x))
    .Call(saige_get_sparse_info, x)
}


# fit the null model
seqFitNullGLMM_SPA <- function(formula, data, gdsfile=NULL, grm.mat=NULL,
    trait.type=c("binary", "quantitative"), sample.col="sample.id", maf=0.01,
    missing.rate=0.01, max.num.snp=1000000L, variant.id=NULL,
    variant.id.varratio=NULL, nsnp.sub.random=2000L, rel.cutoff=0.125,
    inv.norm=c("residuals", "quant", "none"), use.cateMAC=FALSE,
    cateMAC.inc.maf=TRUE, cateMAC.simu=FALSE, X.transform=TRUE, tol=0.02,
    maxiter=20L, nrun=30L, tolPCG=1e-5, maxiterPCG=500L, num.marker=30L,
    tau.init=c(0,0), traceCVcutoff=0.0025, ratioCVcutoff=0.001,
    geno.sparse=TRUE, num.thread=1L, model.savefn="", seed=200L,
    fork.loading=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(formula, "formula"))
    stopifnot(is.data.frame(data))
    stopifnot(is.null(gdsfile) | inherits(gdsfile, "SeqVarGDSClass") |
        is.character(gdsfile))
    trait.type <- match.arg(trait.type)
    stopifnot(is.character(sample.col), length(sample.col)==1L, !is.na(sample.col))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.numeric(max.num.snp), length(max.num.snp)==1L)
    stopifnot(is.null(variant.id) | is.vector(variant.id))
    stopifnot(is.null(variant.id.varratio) | is.vector(variant.id.varratio))
    stopifnot(is.numeric(nsnp.sub.random), length(nsnp.sub.random)==1L,
        nsnp.sub.random>=0L)
    stopifnot(is.numeric(rel.cutoff), length(rel.cutoff)==1L)
    if (is.na(rel.cutoff)) rel.cutoff <- -Inf
    if (is.logical(inv.norm))
        inv.norm <- if (isTRUE(inv.norm)) "residuals" else "none"
    inv.norm <- match.arg(inv.norm)
    use.cateMAC <- .check_use_cateMAC(use.cateMAC)
    if (is.logical(cateMAC.inc.maf))
    {
        stopifnot(length(cateMAC.inc.maf) == 1L)
    } else if (is.numeric(cateMAC.inc.maf))
    {
        stopifnot(is.vector(cateMAC.inc.maf))
        if (anyNA(cateMAC.inc.maf))
            stop("'cateMAC.inc.maf' should not include NA/NaN.")
        if (any(cateMAC.inc.maf<=0 | cateMAC.inc.maf>=1))
            stop("'cateMAC.inc.maf' should be between 0 and 1.")
    } else {
        stop("'cateMAC.inc.maf' should be FALSE, TRUE or a numeric vector for MAF.")
    }
    stopifnot(is.logical(cateMAC.simu), length(cateMAC.simu)==1L)
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
    stopifnot(is.logical(fork.loading), length(fork.loading)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (verbose)
    {
        .cat(.crayon_inverse("SAIGE association analysis:"))
        .cat(.crayon_underline(date()))
    }

    # check GRM matrix if specified
    if (isFALSE(grm.mat)) grm.mat <- NULL
    noRE <- is.null(gdsfile) && is.null(grm.mat)  # no random effect
    grm.mat <- .check_grm_mat(grm.mat, verbose)

    # GDS file
    if (is.character(gdsfile))
    {
        if (verbose)
            .cat("Open ", sQuote(gdsfile))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else if (!is.null(gdsfile))
    {
        # save the filter on GDS file
        seqFilterPush(gdsfile)
        on.exit(seqFilterPop(gdsfile))
    }

    # show warnings immediately
    saveopt <- options(warn=1L)
    on.exit(options(warn=saveopt$warn), add=TRUE)
    set.seed(seed)

    # variables in the formula
    s <- as.character(formula)
    formula_str <- paste(s[2L], s[1L], s[3L])
    vars <- all.vars(formula)
    phenovar <- all.vars(formula)[1L]
    y <- data[[phenovar]]
    if (is.null(y))
        stop("There is no '", phenovar, "' in the input data frame.")
    if (!is.factor(y) && !is.numeric(y) && !is.logical(y))
        stop("The response variable should be numeric or a factor.")

    # check sample id
    if (sample.col %in% vars)
        stop(sprintf("'%s' should not be in the formula.", sample.col))
    if (!(sample.col %in% colnames(data)))
        stop(sprintf("'%s' should be one of the columns in 'data'.", sample.col))
    if (is.factor(data[[sample.col]]))
        stop(sprintf("'%s' should not be a factor variable.", sample.col))
    if (any(is.na(data[[sample.col]])))
        stop(sprintf("'%s' should not have any missing value.", sample.col))
    if (anyDuplicated(data[[sample.col]]))
        stop(sprintf("'%s' in data should be unique.", sample.col))

    # remove missing values
    data <- data[, c(sample.col, vars)]
    data <- na.omit(data)
    data <- droplevels(data)
    sid <- NULL
    if (!is.null(gdsfile))
    {
        seqResetFilter(gdsfile, sample=TRUE, verbose=FALSE)
        sid <- seqGetData(gdsfile, "sample.id")
    }
    if (!is.null(grm.mat) && !isTRUE(grm.mat))
    {
        if (is.null(sid))
            sid <- colnames(grm.mat)
        else
            sid <- intersect(sid, colnames(grm.mat))
        if (length(sid) <= 0L)
            stop("'gdsfile' and 'grm.mat' should have shared sample IDs.")
    } else if (is.null(sid))
        sid <- data[[sample.col]]
    i <- match(sid, data[[sample.col]])
    i <- i[!is.na(i)]
    data <- data[i, ]
    if (nrow(data) <= 0L)
        stop("No common sample.id between 'data' and the GDS file.")
    if (!is.null(gdsfile))
        seqSetFilter(gdsfile, sample.id=data[[sample.col]], verbose=FALSE)
    sid <- data[[sample.col]]

    if (!is.null(gdsfile))
    {
        # use gds genotype file
        if (is.null(variant.id))
            seqResetFilter(gdsfile, sample=FALSE, variant=TRUE, verbose=FALSE)
        else
            seqSetFilter(gdsfile, variant.id=variant.id, verbose=verbose)

        # filters of maf, mac, missing.rate
        if (verbose)
            cat("Filtering variants:\n")
        v <- seqGetAF_AC_Missing(gdsfile, minor=TRUE, parallel=num.thread,
            verbose=verbose)
        sel <- (v$miss <= missing.rate) & (v$ac > 0)
        sel[is.na(sel)] <- FALSE
        s1 <- sel & (v$af >= maf)  # variants used in GRM

        # if using the list of variant IDs to estimate variance ratio
        if (!is.null(variant.id.varratio))
        {
            seqFilterPush(gdsfile)
            seqSetFilter(gdsfile, variant.id=variant.id.varratio,
                verbose=FALSE)
            v <- seqGetAF_AC_Missing(gdsfile, minor=TRUE, parallel=num.thread,
                verbose=FALSE)
            sel <- (v$miss <= missing.rate) & (v$ac > 0)
            sel[is.na(sel)] <- FALSE
            if (verbose)
            {
                .cat("# of variants specified in the variance ratio estimation: ",
                    nrow(v))
            }
        }

        # need random markers in variance ratio estimation
        set.seed(seed)
        if (isFALSE(use.cateMAC))
        {
            if (verbose)
                cat("MAC category for estimating variance ratio:\n")
            last <- 20  # mac threshold used in SAIGE
            mac_cat <- Inf
        } else {
            # using MAC categories
            stopifnot(is.numeric(use.cateMAC), is.vector(use.cateMAC))
            if (verbose)
                cat("MAC categories for estimating variance ratios:\n")
            if (isTRUE(cateMAC.inc.maf))
                cateMAC.inc.maf <- ifelse(is.finite(maf), maf, numeric())
            if (is.numeric(cateMAC.inc.maf))
            {
                for (m in sort(cateMAC.inc.maf))
                {
                    a <- 2L * length(sid) * m
                    if (all(abs(a-use.cateMAC) > 1L) && (a > max(use.cateMAC)))
                    {
                        use.cateMAC <- sort(c(use.cateMAC, a))
                        if (verbose)
                        {
                            .cat("    (including a cut point ", a,
                                " according to MAF=", m, ")")
                        }
                    }
                }
            }
            last <- .Machine$double.eps  # so that > 0
            mac_cat <- c(use.cateMAC, Inf)
        }
        rand.packed.geno <- vector("list", length(mac_cat))
        rand.packed.geno.vid <- NULL
        for (k in seq_along(mac_cat))
        {
            mac <- mac_cat[k]
            s2 <- sel & (last <= v$ac) & (v$ac < mac)
            ii <- which(s2)
            if (length(ii) < num.marker)
            {
                if (!isFALSE(use.cateMAC) && isTRUE(cateMAC.simu))
                {
                    if (verbose)
                    {
                        cat(sprintf("    MAC%s, %g):\t",
                            if (last > .Machine$double.eps) paste0("[", last) else "(0",
                            mac))
                    }
                    n <- 5L * num.marker
                    rand.packed.geno[[k]] <- .simu_geno(length(sid), n, last, mac)
                    rand.packed.geno.vid <- c(rand.packed.geno.vid,
                        paste0("simu", k, "_g", seq_len(n)))
                    if (verbose)
                        cat(sprintf("%d+ simulated variants\n", num.marker))
                } else {
                    stop(sprintf("Less variants (n=%d) than %d in MAC[%g, %g)",
                        length(ii), num.marker, last, mac),
                        ", consider using simulated genotypes via 'cateMAC.simu=TRUE'.")
                }
            } else {
                n <- length(ii)
                if (length(ii) >= num.marker*5L)
                {
                    ii <- sample(ii, num.marker*5L)
                } else {
                    ii <- sample(ii, length(ii))
                }
                if (verbose)
                {
                    if (last > .Machine$double.eps)
                        a <- paste0("[", last)
                    else
                        a <- "(0"
                    if (n > num.marker)
                    {
                        .cat(sprintf(
                            "    MAC%s, %g):\t%d+ randomly from %s variants",
                            a, mac, num.marker, .pretty(n)))
                    } else {
                        .cat(sprintf("    MAC%s, %g):\t%d variants", a, mac, n))
                    }
                }
                seqFilterPush(gdsfile)
                seqSetFilter(gdsfile, variant.sel=ii, action="intersect",
                    verbose=FALSE)
                i <- match(ii, sort(ii))
                rand.packed.geno[[k]] <- seqGet2bGeno(gdsfile, verbose=FALSE)[,i]
                rand.packed.geno.vid <- c(rand.packed.geno.vid,
                    seqGetData(gdsfile, "variant.id")[i])
                seqFilterPop(gdsfile)
            }
            last <- mac
        }

        # SNP markers used in GRM
        if (!is.null(variant.id.varratio))
            seqFilterPop(gdsfile)
        seqSetFilter(gdsfile, variant.sel=s1, action="intersect", verbose=FALSE)
        dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
        n_samp <- dm[2L]
        n_var  <- n <- dm[3L]
        if (max.num.snp>0L && n_var>max.num.snp)
        {
            set.seed(seed)
            v <- seqGetData(gdsfile, "$variant_index")
            seqSetFilter(gdsfile, variant.sel=sample(v, max.num.snp),
                warn=FALSE, verbose=FALSE)
            n_var <- as.integer(max.num.snp)
        }

        if (verbose)
        {
            .cat("Fit the null model: ", format(formula),
                ifelse(noRE, "", " + var(GRM)"))
            .cat("    # of samples: ", .pretty(n_samp))
            cat("    # of variants in GRM:", .pretty(n_var))
            if (n > max.num.snp)
                cat(" (randomly selected from ", .pretty(n), ")", sep="")
            cat("\n")
            .cat("    MAF threshold for GRM: >= ", maf)
        }
    } else {
        n_samp <- length(sid)
        n_var  <- NA_integer_
        rand.packed.geno <- NULL
        if (verbose)
        {
            .cat("Fit the null model: ", format(formula),
                ifelse(noRE, "", "+ var(GRM)"))
            .cat("    # of samples: ", .pretty(n_samp))
        }
    }

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    .Call(saige_set_numthread, num.thread)
    if (verbose)
        .cat("    using ", num.thread, " thread", ifelse(num.thread>1L, "s", ""))

    # rearrange grm.mat if needed
    if (isTRUE(grm.mat))
    {
        # need a sparse GRM from the full GRM
        # correct nsnp.sub.random
        if (nsnp.sub.random > 0L)
        {
            nsnp.sub.random <- as.integer(floor(nsnp.sub.random/4) * 4L)
            if (nsnp.sub.random < 1L) nsnp.sub.random <- 4L
        }
        if (verbose)
        {
            cat("Calculating the sparse GRM:\n")
            .cat("    using relatedness threshold: ", rel.cutoff)
        }
        # calculate sparse GRM
        grm.mat <- .fit_calc_sp_grm(gdsfile, nsnp.sub.random, maf,
            missing.rate, rel.cutoff, num.thread, FALSE, FALSE, verbose)
        gc(verbose=FALSE)  # reduce memory usage
        if (verbose) cat("Done (sparse GRM)\n")
    }
    if (!is.null(grm.mat) && !identical(sid, colnames(grm.mat)))
    {
        i <- match(sid, colnames(grm.mat))
        grm.mat <- grm.mat[i, i]
    }

    X <- model.matrix(formula, data)
    if (NCOL(X) <= 1L) X.transform <- FALSE
    if (isTRUE(X.transform))
    {
        if (verbose)
            cat("Transform on the design matrix with QR decomposition:\n")
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
                .cat("    exclude ", length(i_na), " covariates (",
                    paste(colnames(X)[i_na], collapse=", "),
                    ") to avoid multi collinearity.")
            }
        }
        X_name <- colnames(X)
        Xqr <- qr(X)  # QR decomposition
        X_new <- qr.Q(Xqr) * sqrt(nrow(X))
        X_qrr <- qr.R(Xqr)
        data <- data.frame(cbind(y, X_new))
        nm <- paste0("x_", seq_len(ncol(X_new))-1L)
        colnames(data) <- c("y", nm)
        formula <- as.formula(paste("y ~", paste(nm, collapse=" + "), "-1"))
        if (verbose)
            .cat("    new formula: ", format(formula))
    }

    # clear the internal matrix
    .Call(saige_init_fit_grm)
    # internal buffer for diagonal of GRM
    buf_sigma_diag <- double(n_samp)

    # load SNP genotypes
    if (!is.null(gdsfile))
    {
        if (verbose)
        {
            cat("Use dense genetic relationship matrix in the file:\n")
            .cat("    ", gdsfile$filename)
            cat("Loading SNP genotypes from the GDS file:\n")
        }
        nfork <- 1L
        if (SeqArray:::.IsForking(num.thread) && isTRUE(fork.loading))
            nfork <- num.thread
        if (isTRUE(geno.sparse))
        {
            # sparse genotypes
            packed.geno <- .get_sparse_geno(gdsfile, nfork, verbose)
        } else {
            # 2-bit packed genotypes
            packed.geno <- seqGet2bGeno(gdsfile, verbose=verbose)
        }
        if (verbose)
        {
            .cat("    using ",
                .pretty_size(as.double(object.size(packed.geno))),
                " (", ifelse(isTRUE(geno.sparse), "sparse", "dense"),
                " genotype matrix)")
        }

        # initialize internal variables and buffers
        buf_std_geno <- double(4L*n_var)
        buf_crossprod <- matrix(0.0, nrow=n_samp, ncol=num.thread)
        if (isTRUE(geno.sparse))
        {
            .Call(saige_store_sp_geno, packed.geno, rand.packed.geno,
                n_samp, buf_std_geno, buf_sigma_diag, buf_crossprod)
        } else {
            .Call(saige_store_2b_geno, packed.geno, rand.packed.geno,
                n_samp, buf_std_geno, buf_sigma_diag, buf_crossprod)
        }
    }

    # if use GRM matrix
    if (is.matrix(grm.mat))
    {
        # set dense GRM
        .Call(saige_store_dense_grm, n_samp, grm.mat, buf_sigma_diag)
        if (verbose)
        {
            cat("User-defined")
            if (!is.null(gdsfile)) cat(" approximate")
            cat(" genetic relationship matrix:\n")
            cat(sprintf("    %d x %d (dense matrix)\n", n_samp, n_samp))
        }
    } else if (inherits(grm.mat, "sparseMatrix"))
    {
        # column-oriented sparse form
        grm.mat <- .sp_to_dgCMatrix(grm.mat)
        # set sparse matrix
        .Call(saige_store_sparse_grm, n_samp, grm.mat, buf_sigma_diag)
        if (verbose)
        {
            a <- nnzero(grm.mat)
            s <- sprintf("%.3f%%", a/prod(dim(grm.mat))*100)
            if (s=="0.000%") s <- "<0.001%"
            if (!is.null(gdsfile)) cat("Approximate u") else cat("U")
            cat("ser-defined sparse genetic relationship matrix:\n")
            cat(sprintf("    %d x %d, # of nonzero: %d (%s)\n",
                n_samp, n_samp, a, s))
        }
    } else if (is.null(grm.mat) && is.null(gdsfile) && verbose)
        cat("Assuming independent outcomes\n")

    # parameters for fitting the model
    param <- list(
        trait = match(trait.type, .trait_list),
        num.thread = num.thread, seed = seed,
        tol = tol, tolPCG = tolPCG,
        maxiter = maxiter, maxiterPCG = maxiterPCG, no_iteration = FALSE,
        nrun = nrun, num.marker = num.marker,
        traceCVcutoff = traceCVcutoff, ratioCVcutoff = ratioCVcutoff,
        verbose = verbose,
        indent = ""
    )

    tau.init[is.na(tau.init)] <- 0
    tau.init[tau.init < 0] <- 0

    # fit the model
    if (trait.type == "binary")
    {
        # logistic regression
        glmm <- .fit_binary(verbose, X.transform, phenovar, data, formula,
            param, tau.init, gdsfile, grm.mat, seed, n_var)
    } else if (trait.type == "quantitative")
    {
        # linear regression
        glmm <- .fit_quant(verbose, X.transform, phenovar, data, inv.norm,
            formula, param, tau.init, gdsfile, grm.mat, seed, n_var)
    } else
        stop("Invalid 'trait.type'.")

    glmm <- c(list(formula=formula_str), glmm)
    glmm$use.cateMAC <- use.cateMAC
    if (is.data.frame(glmm$var.ratio))
    {
        glmm$var.ratio$id <- rand.packed.geno.vid[glmm$var.ratio$id]
        if (verbose)
        {
            CV <- function(x) sd(x)/(mean(x)*length(x))
            if (isFALSE(use.cateMAC))
            {
                m <- mean(glmm$var.ratio$ratio)
                s <- sd(glmm$var.ratio$ratio)
                cv <- CV(glmm$var.ratio$ratio)
                .cat("    ratio avg: ", m, ", sd: ", s, ", CV: ", cv)
            } else {
                x <- glmm$var.ratio$ratio
                d <- cut(glmm$var.ratio$mac, c(0, use.cateMAC, Inf),
                    right=FALSE, dig.lab=15L)
                m <- suppressWarnings(data.frame(
                    n  = vapply(levels(d), function(s) sum(d==s), 0L),
                    mean = vapply(levels(d), function(s) mean(x[d==s]), 0),
                    sd = vapply(levels(d), function(s) sd(x[d==s]), 0),
                    CV = vapply(levels(d), function(s) CV(x[d==s]), 0),
                    stringsAsFactors=FALSE
                ))
                rownames(m) <- paste("   ", rownames(m))
                cat("    MAC categories:\n")
                print(m)
            }
        }
    }

    # tweak the result
    if (!isTRUE(X.transform))
    {
        names(glmm$coefficients) <- colnames(glmm$obj.noK$X1)
    } else {
        coef <- solve(X_qrr, glmm$coefficients * sqrt(nrow(data)))
        names(coef) <- X_name
        glmm$coefficients <- coef
    }
    names(glmm$tau) <- c("Sigma_E", "Sigma_G")
    glmm$trait.type <- trait.type
    if (!is.null(gdsfile))
    {
        glmm$sample.id <- seqGetData(gdsfile, "sample.id")
        glmm$variant.id <- seqGetData(gdsfile, "variant.id")
    } else {
        glmm$sample.id <- sid
        glmm$variant.id <- NULL
    }
    class(glmm) <- "ClassSAIGE_NullModel"

    if (!is.na(model.savefn) && model.savefn!="")
    {
        .cat("Save the model to ", sQuote(model.savefn))
        if (grepl("\\.(rda|RData)$", model.savefn, ignore.case=TRUE))
        {
            .glmm <- glmm
            save(.glmm, file=model.savefn)
        } else if (grepl("\\.rds$", model.savefn, ignore.case=TRUE))
        {
            saveRDS(glmm, file=model.savefn)
        } else {
            stop("Unknown format of the output file, and it should be RData or RDS.")
        }
    }
    if (verbose)
    {
        .cat(.crayon_underline(date()))
        .cat(.crayon_inverse("Done."))
    }

    if (!is.na(model.savefn) && model.savefn!="")
        return(invisible(glmm))
    else
        return(glmm)
}
