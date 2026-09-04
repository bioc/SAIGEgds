#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2019-2026    Xiuwen Zheng / AbbVie-ComputationalGenomics
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
    # if RsparseMatrix, convert to TsparseMatrix
    if (inherits(m, "RsparseMatrix"))
        m <- as(m, "TsparseMatrix")
    if (inherits(m, "TsparseMatrix"))
        m <- as(m, "CsparseMatrix")
    if (inherits(m, "dsCMatrix"))
        m <- as(m, "generalMatrix")
    if (!inherits(m, "dgCMatrix"))
        m <- as(m, "dgCMatrix")
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
            .cat(.crayon_underline(.tm()))
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
    tau.init, gdsfile, grm.mat, seed, n_varm, calc_vr=TRUE)
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
        rowmore <- nrow(v) > 10L
        if (rowmore) v <- v[1:10, ]
        print(v, row.names=FALSE)
        if (rowmore) cat("        ...\n")
        if (n > 2L)
            stop("The outcome variable has more than 2 categories!")
        else if (n <= 1L)
            stop("The outcome variable has only a unique value!")
    }

    # fit the null model
    fit0 <- glm(formula, data=data, family=binomial, offset=param$covoffset)
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
        if (isTRUE(param$no_iteration))
        {
            tau <- tau.init
        } else {
            tau <- fixtau <- c(0, 0)
            if (fit0$family$family %in% c("binomial", "poisson"))
                tau[1] <- fixtau[1] <- 1
            if (sum(tau.init[fixtau==0]) == 0)
                tau[fixtau==0] <- 0.1
            else
                tau[fixtau==0] <- tau.init[fixtau==0]
        }
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

    # use updated mu to set obj.noK (ScoreTest_NULL_Model)
    if (!is.null(param$Xmat)) X <- param$Xmat
    mu <- glmm$fitted.values
    V <- mu*(1-mu)
    XV <- t(X * V)
    XVX_inv <- solve(crossprod(X, X * V))
    XXVX_inv <- X %*% XVX_inv
    obj.noK <- list(y=unname(fit0$y), mu=mu, V=V, X1=X, XV=XV,
        XXVX_inv=XXVX_inv)
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
    if (calc_vr)
        .calcVR(gdsfile, seed, fit0, glmm, obj.noK, param, verbose)
    else
        glmm
}


# survival (time-to-event) outcome: Cox proportional-hazards frailty model
# fitted via the Breslow Cox<->Poisson equivalence (cf. GATE).
.fit_survival <- function(verbose, X.transform, phenovar, data, formula, param,
    tau.init, gdsfile, grm.mat, seed, n_var, calc_vr=TRUE)
{
    # event status must be 0 (censored) / 1 (event)
    ystatus <- data[[phenovar]]
    uy <- sort(unique(ystatus))
    if (length(uy)!=2L || !all(uy %in% c(0, 1)))
        stop("The event status (response) should be 0 (censored) or 1 (event).")
    if (verbose)
    {
        ev <- sum(ystatus==1)
        .cat("Survival outcome (event status): ", phenovar)
        .cat("    # of events: ", ev, " (",
            sprintf("%.2f%%", 100*ev/length(ystatus)), "), # censored: ",
            sum(ystatus==0))
    }
    if (is.null(gdsfile) && is.null(grm.mat))
        stop("Survival analysis requires a GRM (via 'gdsfile' or 'grm.mat').")

    # initial fixed-effect coefficients via logistic regression
    # (no covariate offset for survival: use.offset is forced FALSE)
    fit0 <- glm(formula, data=data, family=binomial)
    if (verbose)
    {
        cat("Initial fixed-effect coefficients:\n")
        v <- as.data.frame(t(fit0$coefficients))
        rownames(v) <- "   "
        print(v, width=128L)
    }

    # design matrix WITHOUT the intercept (Cox baseline absorbs it)
    cn <- colnames(model.matrix(fit0))
    X <- unname(model.matrix(fit0))
    attr(X, "assign") <- attr(X, "contrasts") <- NULL
    icpt <- which(cn == "(Intercept)")
    if (length(icpt))
    {
        X <- X[, -icpt, drop=FALSE]
        coef0 <- fit0$coefficients[-icpt]
    } else
        coef0 <- fit0$coefficients
    if (NCOL(X) < 1L)
        stop("Survival analysis needs at least one covariate (besides the ",
            "intercept).")
    # coefficients passed to the C++ fitter must match X (no intercept term)
    fit0$coefficients <- coef0

    # initial tau: Sigma_E fixed at 1 (Poisson), estimate Sigma_G
    if (isTRUE(param$no_iteration))
        tau <- tau.init
    else {
        tau <- c(1, 0.1)
        if (sum(tau.init[2L]) != 0) tau[2L] <- tau.init[2L]
    }
    # iterate the null model (Cox-via-Poisson IRLS + AI-REML for tau)
    glmm <- .Call(saige_fit_AI_PCG, fit0, X, tau, param)

    # score-test null object with Poisson variance V = mu.
    # The intercept was dropped for *fitting* (absorbed by the baseline hazard),
    # but the score-test projection must use the intercept-augmented design
    # Xa = [X, 1] (GATE/SAIGE's X1_fg = cbind(X1, 1)). Projecting the genotype
    # on [X, 1] makes the adjusted genotype g_tilde orthogonal (in the V=mu
    # metric) to the baseline-hazard direction too, so that m1 = sum(mu*g_t) = 0.
    # This is required for the correct score variance var2 = sum(mu*g_tilde^2)
    # AND for the Poisson-SPA CGF: sum(y-mu)=0 makes the *score* invariant to a
    # constant shift of g_tilde, but var2 and the CGF are not -- the intercept
    # column removes the baseline direction that estimating Lambda0 accounts for.
    if (!is.null(param$Xmat)) X <- param$Xmat
    mu <- glmm$fitted.values
    Xa <- cbind(X, 1)              # intercept-augmented projection design
    V <- mu
    XV <- t(Xa * V)
    XVX_inv <- solve(crossprod(Xa, Xa * V))
    XXVX_inv <- Xa %*% XVX_inv
    obj.noK <- list(y=unname(fit0$y), mu=mu, V=V, X1=Xa, XV=XV,
        XXVX_inv=XXVX_inv)
    # no Sigma_inv: the frailty correction is captured by the variance ratio
    # (computed from the GRM via PCG in saige_calc_var_ratio).
    obj.noK$Sigma_inv <- FALSE
    glmm$obj.noK <- obj.noK

    # calculate the variance ratio
    if (calc_vr)
        .calcVR(gdsfile, seed, fit0, glmm, obj.noK, param, verbose)
    else
        glmm
}


# quantitative outcome
.fit_quant <- function(verbose, X.transform, phenovar, data, inv.norm, formula,
    param, tau.init, gdsfile, grm.mat, seed, n_var, calc_vr=TRUE)
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
        if (isTRUE(param$no_iteration))
        {
            tau <- tau.init
        } else {
            tau <- tau.init
            if (sum(tau) == 0) tau <- c(0.5, 0.5)
            tau <- var(Y) * tau / sum(tau)
        }
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
    if (calc_vr)
        .calcVR(gdsfile, seed, fit0, glmm, obj.noK, param, verbose)
    else
        glmm
}


# check use.cateMAC
.cateMAC_default <- c(5.5, 10.5, 20.5)

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
        cat("MAC categories:")
        if (is.logical(use.cateMAC))
        {
            if (isFALSE(use.cateMAC))
                cat(" none\n")
            else
                .cat(" ", use.cateMAC)
        } else {
            x1 <- c(0, use.cateMAC)
            x2 <- c(use.cateMAC, Inf)
            s <- paste0("[", x1, ",", x2, ")")
            .cat("\n    ", paste(s, collapse=", "))
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

# get a sparse form of genotypes
.get_sparse_geno <- function(gdsfile, nproc, verbose)
{
    # check
    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile)==1L)
        if (isTRUE(verbose))
            .cat("Open ", sQuote(basename(gdsfile)))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    }
    # integer genotypes or numeric dosages
    if (exist.gdsn(gdsfile, "genotype/data"))
    {
        varnm <- "$dosage_alt2"
    } else if (exist.gdsn(gdsfile, "annotation/format/DS/data"))
    {
        varnm <- "annotation/format/DS"
        if (verbose) cat("    using 'annotation/format/DS'\n")
    } else {
        stop("'genotype' and 'annotation/format/DS' are not available.")
    }
    # initialize
    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    nsamp <- dm[2L]; nvar  <- dm[3L]
    remove(dm)
    # process
    stopifnot(is.numeric(nproc))
    if (nproc > 1L)
    {
        # parallel loading
        cl <- parallel::makeCluster(nproc)
        on.exit(parallel::stopCluster(cl), add=TRUE)
        bs <- ceiling(nvar/100L)  # block size
        seqParallel(cl, gdsfile, FUN=function(gds, varnm)
        {
            gc(FALSE, reset=TRUE)
            seqApply(gds, varnm, .cfunction("saige_get_sparse"),
                as.is="list", .useraw=TRUE, .list_dup=FALSE)
        },
        .initialize=function(proc_id, nsamp)
        {
            # internal buffer
            buf_b256 <- integer((ceiling(nsamp/256L) + 1L)*3L + 1L)
            buf_b1 <- raw(3L*ceiling(nsamp/256L)*257L)
            assign("buf_b256", buf_b256, envir=.PkgEnv)
            assign("buf_b1", buf_b1, envir=.PkgEnv)
            .Call(saige_init_sparse, nsamp, buf_b256, buf_b1)
        },
        .finalize=function(proc_id, nsamp)
        {
            .PkgEnv$buf_b256 <- .PkgEnv$buf_b1 <- NULL
            remove("buf_b256", "buf_b1", envir=.PkgEnv)
        }, .initparam=nsamp, .balancing=TRUE, .bl_size=bs, .bl_progress=verbose,
        varnm=varnm)
    } else {
        # internal buffer
        buf_b256 <- integer((ceiling(nsamp/256L) + 1L)*3L + 1L)
        buf_b1 <- raw(3L*ceiling(nsamp/256L)*257L)
        .Call(saige_init_sparse, nsamp, buf_b256, buf_b1)
        # apply
        seqApply(gdsfile, varnm, .cfunction("saige_get_sparse"),
            as.is="list", .useraw=TRUE, .list_dup=FALSE, .progress=verbose)
    }
}

# get information from the output of .get_sparse_geno()
.get_sparse_smallblock_len <- function(x)
{
    stopifnot(is.raw(x))
    .Call(saige_get_sparse_info, x)
}


# load genotypes and prepare random markers from a GDS file
.load_geno_grm_gds <- function(gdsfile, sid, variant.id, variant.id.varratio,
    use.cateMAC, cateMAC.inc.maf, cateMAC.simu, maf, missing.rate,
    max.num.snp, num.marker, seed, noRE, formula, num.thread, verbose)
{
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
                .cat("# of variants specified in ",
                    "the variance ratio estimation: ", nrow(v))
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
                                if (last > .Machine$double.eps)
                                    paste0("[", last) else "(0",
                                mac))
                    }
                    n <- 5L * num.marker
                    rand.packed.geno[[k]] <-
                        .simu_geno(length(sid), n, last, mac)
                    rand.packed.geno.vid <- c(rand.packed.geno.vid,
                        paste0("simu", k, "_g", seq_len(n)))
                    if (verbose)
                        cat(sprintf("%d+ simulated variants\n", num.marker))
                } else {
                    stop(sprintf("Less variants (n=%d) than %d in MAC[%g, %g)",
                            length(ii), num.marker, last, mac),
                        ", consider using simulated genotypes ",
                        "via 'cateMAC.simu=TRUE'.")
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
                rand.packed.geno[[k]] <-
                    seqGet2bGeno(gdsfile, verbose=FALSE)[,i]
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
        rand.packed.geno.vid <- NULL
        if (verbose)
        {
            .cat("Fit the null model: ", format(formula),
                ifelse(noRE, "", "+ var(GRM)"))
            .cat("    # of samples: ", .pretty(n_samp))
        }
    }

    # output
    list(n_samp=n_samp, n_var=n_var, rand.packed.geno=rand.packed.geno,
        rand.packed.geno.vid=rand.packed.geno.vid, use.cateMAC=use.cateMAC)
}


# fit the null model
seqFitNullGLMM_SPA <- function(formula, data, gdsfile=NULL, grm.mat=NULL,
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
{
    # check
    stopifnot(inherits(formula, "formula"))
    stopifnot(is.data.frame(data))
    stopifnot(is.null(gdsfile) || inherits(gdsfile, "SeqVarGDSClass") ||
            is.character(gdsfile))
    trait.type <- match.arg(trait.type)
    if (trait.type == "survival")
    {
        if (is.null(event.time) || !is.character(event.time) ||
                length(event.time)!=1L || is.na(event.time))
            stop("'event.time' should be a column name in 'data' giving the ",
                "event/censoring time for survival analysis.")
        if (!(event.time %in% colnames(data)))
            stop(sprintf("'%s' should be one of the columns in 'data'.",
                event.time))
        # the Cox baseline hazard absorbs the intercept, which is therefore
        # non-identifiable: disable the QR transform and covariate offset so
        # the intercept can be removed explicitly (as in GATE)
        X.transform <- FALSE
        use.offset <- FALSE
    } else {
        event.time <- NULL
    }
    stopifnot(is.character(sample.col), length(sample.col)==1L,
        !is.na(sample.col))
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
        stop("'cateMAC.inc.maf' should be FALSE, TRUE or ",
            "a numeric vector for MAF.")
    }
    stopifnot(is.logical(cateMAC.simu), length(cateMAC.simu)==1L)
    stopifnot(is.logical(use.offset), length(use.offset)==1L)
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
    stopifnot(is.logical(save.packed.geno), length(save.packed.geno)==1L)
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.character(model.savefn), length(model.savefn)==1L)
    stopifnot(is.numeric(seed), length(seed)==1L, is.finite(seed))
    stopifnot(is.logical(parallel.loading), length(parallel.loading)==1L)
    stopifnot(is.logical(use.gpu), length(use.gpu)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (!missing(fork.loading))
    {
        warning("'fork.loading' is deprecated, ",
            "please use 'parallel.loading' instead.")
    }
    if (verbose)
    {
        .cat(.crayon_inverse("SAIGE association analysis:"))
        .cat(.crayon_underline(.tm()))
    }

    # initialize GPU if requested
    if (isTRUE(use.gpu))
    {
        use.gpu <- .Call(saige_gpu_init, verbose)
        if (!use.gpu && verbose)
            cat("GPU not available, falling back to CPU.\n")
        if (use.gpu)
        {
            on.exit(.Call(saige_gpu_cleanup), add=TRUE)
            geno.sparse <- FALSE  # sparse genotype is not supported for GPU
        }
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
    if (!is.null(seed)) set.seed(seed)

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
    {
        stop(sprintf("'%s' should be one of the columns in 'data'.",
                sample.col))
    }
    if (is.factor(data[[sample.col]]))
        stop(sprintf("'%s' should not be a factor variable.", sample.col))
    if (any(is.na(data[[sample.col]])))
        stop(sprintf("'%s' should not have any missing value.", sample.col))
    if (anyDuplicated(data[[sample.col]]))
        stop(sprintf("'%s' in data should be unique.", sample.col))

    # remove missing values
    if (!is.null(event.time))
    {
        if (event.time %in% vars)
            stop("'event.time' should not be in the formula.")
        data <- data[, c(sample.col, vars, event.time)]
    } else {
        data <- data[, c(sample.col, vars)]
    }
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
    # survival: drop subjects censored before the first event time -- they
    # carry no risk-set information (Lambda0 = mu = 0). GATE removes these
    # before fitting; do it here, before the sample set (GDS filter, GRM subset,
    # n_samp) is finalized, so all downstream buffers stay aligned.
    if (!is.null(event.time))
    {
        st <- data[[phenovar]]
        st <- if (is.factor(st)) as.numeric(as.character(st)) else
            as.numeric(st)
        tm <- data[[event.time]]
        ev1 <- which(st == 1)
        if (length(ev1) > 0L)
        {
            rm_i <- which(st == 0 & tm < min(tm[ev1]))
            if (length(rm_i) > 0L)
            {
                if (verbose)
                    .cat("    survival: removing ", length(rm_i), " subject",
                        if (length(rm_i)>1L) "s" else "",
                        " censored before the first event time")
                data <- data[-rm_i, ]
                if (nrow(data) <= 0L)
                    stop("No samples remain after removing early-censored ",
                        "subjects for survival analysis.")
            }
        }
    }
    if (!is.null(gdsfile))
        seqSetFilter(gdsfile, sample.id=data[[sample.col]], verbose=FALSE)
    sid <- data[[sample.col]]

    # load genotypes for GRM
    v <- .load_geno_grm_gds(gdsfile, sid, variant.id, variant.id.varratio,
        use.cateMAC, cateMAC.inc.maf, cateMAC.simu, maf, missing.rate,
        max.num.snp, num.marker, seed, noRE, formula, num.thread, verbose)
    n_samp <- v$n_samp
    n_var  <- v$n_var
    rand.packed.geno <- v$rand.packed.geno
    rand.packed.geno.vid <- v$rand.packed.geno.vid
    use.cateMAC <- v$use.cateMAC
    remove(v)

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    .Call(saige_set_numthread, num.thread)
    if (verbose)
    {
        .cat("    using ", num.thread, " thread",
            if (num.thread>1L) "s" else "")
    }

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
            missing.rate, rel.cutoff, num.thread, FALSE, use.gpu,
            verbose, verbose)
        gc(verbose=FALSE, reset=TRUE, full=TRUE)  # reduce memory usage
        if (verbose) cat("Done (sparse GRM)\n")
    }
    if (!is.null(grm.mat) && !identical(sid, colnames(grm.mat)))
    {
        i <- match(sid, colnames(grm.mat))
        grm.mat <- grm.mat[i, i]
    }

    # extract the event time aligned to the final sample order (survival),
    # before any X.transform rebuild of 'data' drops the column
    evtime <- NULL
    if (!is.null(event.time))
    {
        evtime <- data[[event.time]]
        if (!is.numeric(evtime))
            stop("'event.time' column should be numeric.")
        if (anyNA(evtime) || any(evtime < 0))
            stop("'event.time' values should be non-negative and non-missing.")
    }

    X <- model.matrix(formula, data)
    if (NCOL(X) <= 1L) use.offset <- X.transform <- FALSE

    # transform to avoid multi-collinearity and improve numeric stability
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

    # estimate the fixed effect coefficients or not
    covoffset <- Xmat <- NULL
    if (isTRUE(use.offset))
    {
        if (verbose)
        {
            cat("    using covariate offset instead of estimating",
                "each fixed effect coefficient\n")
        }
        Xmat <- model.matrix(formula, data=data)
        if (trait.type == "binary")
        {
            mod <- glm(formula, data=data, family=binomial)
        } else {
            mod <- glm(formula, data=data, family=gaussian)
        }
        covoffset <- Xmat[, -1L, drop=F] %*%  mod$coefficients[-1L]
        formula <- as.formula("y ~ 1")
    }	    

    # clear the internal GRM matrix
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
        nproc <- 1L
        if (isTRUE(parallel.loading)) nproc <- num.thread
        if (isTRUE(geno.sparse))
        {
            # sparse genotypes
            packed.geno <- .get_sparse_geno(gdsfile, nproc, verbose)
        } else {
            # 2-bit packed genotypes
            packed.geno <- seqGet2bGeno(gdsfile, parallel=nproc,
                verbose=verbose)
        }
        if (verbose)
        {
            .cat("    using ",
                .pretty_size(as.double(object.size(packed.geno))),
                " (stored in a ",
                ifelse(isTRUE(geno.sparse), "sparse", "dense"), " form)")
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
        covoffset = covoffset, Xmat = Xmat,
        num.thread = num.thread, seed = seed,
        tol = tol, tolPCG = tolPCG,
        maxiter = maxiter, maxiterPCG = maxiterPCG, no_iteration = FALSE,
        nrun = nrun, num.marker = num.marker,
        traceCVcutoff = traceCVcutoff, ratioCVcutoff = ratioCVcutoff,
        verbose = verbose,
        eventTime = evtime,
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
    } else if (trait.type == "survival")
    {
        # Cox proportional-hazards frailty model (Cox-via-Poisson)
        glmm <- .fit_survival(verbose, X.transform, phenovar, data, formula,
            param, tau.init, gdsfile, grm.mat, seed, n_var)
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
    if (!isTRUE(X.transform) || isTRUE(use.offset))
    {
        if (isTRUE(use.offset))
            names(glmm$coefficients) <- "(Offset)"
        else
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
    if (isTRUE(save.packed.geno))
    {
        glmm$packed.geno <- packed.geno
        glmm$rand.packed.geno <- rand.packed.geno
        glmm$rand.packed.geno.vid <- rand.packed.geno.vid
        glmm$grm.mat <- grm.mat
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
            stop("Unknown format of the output file, ",
                "and it should be RData or RDS.")
        }
    }
    if (verbose)
    {
        .cat(.crayon_underline(.tm()))
        .cat(.crayon_inverse("Done."))
    }

    if (!is.na(model.savefn) && model.savefn!="")
        return(invisible(glmm))
    else
        return(glmm)
}


# refit the null model
seqRefitNullGLMM <- function(formula, data, model=NULL, event.time=NULL,
    sample.col="sample.id", inv.norm=c("residuals", "quant", "none"),
    use.offset=FALSE, X.transform=TRUE, tau.update=FALSE, recalcVR=FALSE,
    tol=0.02, maxiter=20L, nrun=30L, tolPCG=1e-5, maxiterPCG=500L,
    num.marker=30L, traceCVcutoff=0.0025, ratioCVcutoff=0.001,
    num.thread=1L, seed=200L, use.gpu=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(formula, "formula"))
    stopifnot(is.data.frame(data))
    stopifnot(inherits(model, "ClassSAIGE_NullModel"))
    if (is.null(model$packed.geno))
        stop("'model' should be created with 'save.packed.geno=TRUE'.")
    stopifnot(is.character(sample.col), length(sample.col)==1L,
        !is.na(sample.col))
    if (is.logical(inv.norm))
        inv.norm <- if (isTRUE(inv.norm)) "residuals" else "none"
    inv.norm <- match.arg(inv.norm)
    stopifnot(is.logical(use.offset), length(use.offset)==1L)
    stopifnot(is.logical(X.transform), length(X.transform)==1L)
    stopifnot(is.logical(tau.update), length(tau.update)==1L)
    stopifnot(is.logical(recalcVR), length(recalcVR)==1L)
    stopifnot(is.numeric(tol), length(tol)==1L)
    stopifnot(is.numeric(maxiter), length(maxiter)==1L)
    stopifnot(is.numeric(nrun), length(nrun)==1L)
    stopifnot(is.numeric(tolPCG), length(tolPCG)==1L)
    stopifnot(is.numeric(maxiterPCG), length(maxiterPCG)==1L)
    stopifnot(is.numeric(num.marker), length(num.marker)==1L)
    stopifnot(is.numeric(traceCVcutoff), length(traceCVcutoff)==1L)
    stopifnot(is.numeric(ratioCVcutoff), length(ratioCVcutoff)==1L)
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.numeric(seed), length(seed)==1L, is.finite(seed))
    stopifnot(is.logical(use.gpu), length(use.gpu)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        .cat(.crayon_inverse("SAIGE association analysis (refitting):"))
        .cat(.crayon_underline(.tm()))
    }

    # initialize GPU if requested
    if (isTRUE(use.gpu))
    {
        use.gpu <- .Call(saige_gpu_init, verbose)
        if (!use.gpu && verbose)
            cat("GPU not available, falling back to CPU.\n")
        if (use.gpu)
            on.exit(.Call(saige_gpu_cleanup), add=TRUE)
    }

    # extract from saved model
    trait.type <- model$trait.type
    packed.geno <- model$packed.geno
    rand.packed.geno <- model$rand.packed.geno
    rand.packed.geno.vid <- model$rand.packed.geno.vid
    grm.mat <- model$grm.mat
    use.cateMAC <- model$use.cateMAC
    geno.sparse <- is.list(packed.geno)

    # survival: the saved model does not carry the event time, so it must be
    # re-supplied via 'event.time' (a column of the new 'data'). The intercept
    # is absorbed by the baseline hazard, so -- as in seqFitNullGLMM_SPA -- the
    # QR transform and covariate offset are disabled for survival.
    if (trait.type == "survival")
    {
        if (is.null(event.time) || !is.character(event.time) ||
                length(event.time)!=1L || is.na(event.time))
            stop("'event.time' should be a column name in 'data' giving the ",
                "event/censoring time when refitting a survival model.")
        if (!(event.time %in% colnames(data)))
            stop(sprintf("'%s' should be one of the columns in 'data'.",
                event.time))
        X.transform <- FALSE
        use.offset <- FALSE
    } else {
        event.time <- NULL
    }

    # show warnings immediately
    saveopt <- options(warn=1L)
    on.exit(options(warn=saveopt$warn), add=TRUE)
    if (!is.null(seed)) set.seed(seed)

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
    {
        stop(sprintf("'%s' should be one of the columns in 'data'.",
                sample.col))
    }
    if (is.factor(data[[sample.col]]))
        stop(sprintf("'%s' should not be a factor variable.", sample.col))
    if (any(is.na(data[[sample.col]])))
        stop(sprintf("'%s' should not have any missing value.", sample.col))
    if (anyDuplicated(data[[sample.col]]))
        stop(sprintf("'%s' in data should be unique.", sample.col))

    # remove missing values and match samples to model
    if (!is.null(event.time))
    {
        if (event.time %in% vars)
            stop("'event.time' should not be in the formula.")
        data <- data[, c(sample.col, vars, event.time)]
    } else {
        data <- data[, c(sample.col, vars)]
    }
    data <- na.omit(data)
    data <- droplevels(data)
    sid <- model$sample.id
    i <- match(sid, data[[sample.col]])
    if (any(is.na(i)))
        stop("All samples in the model must be present in the new data.")
    data <- data[i, ]
    n_samp <- length(sid)

    # event time aligned to the model's sample order (survival). The saved
    # model's sample set already excludes subjects censored before the first
    # event, so no further removal is done here -- the model's sample.id is
    # used verbatim (and the packed genotypes are aligned to it).
    evtime <- NULL
    if (!is.null(event.time))
    {
        evtime <- data[[event.time]]
        if (!is.numeric(evtime))
            stop("'event.time' column should be numeric.")
        if (anyNA(evtime) || any(evtime < 0))
            stop("'event.time' values should be non-negative and non-missing.")
    }

    # number of variants from saved model
    n_var <- length(model$variant.id)

    if (verbose)
    {
        .cat("Refit the null model: ", format(formula))
        .cat("    # of samples: ", .pretty(n_samp))
        .cat("    # of variants: ", .pretty(n_var))
        .cat("    trait type: ", trait.type)
    }

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    .Call(saige_set_numthread, num.thread)
    if (verbose)
    {
        .cat("    using ", num.thread, " thread",
            if (num.thread>1L) "s" else "")
    }

    # rearrange grm.mat if needed
    if (!is.null(grm.mat) && !identical(sid, colnames(grm.mat)))
    {
        i <- match(sid, colnames(grm.mat))
        if (anyNA(i))
            stop("All samples in the model should be present in GRM.")
        if (!all(i == seq_along(i))) grm.mat <- grm.mat[i, i]
    }

    X <- model.matrix(formula, data)
    if (NCOL(X) <= 1L) use.offset <- X.transform <- FALSE

    # transform to avoid multi-collinearity and improve numeric stability
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

    # estimate the fixed effect coefficients or not
    covoffset <- Xmat <- NULL
    if (isTRUE(use.offset))
    {
        if (verbose)
        {
            cat("    using covariate offset instead of estimating",
                "each fixed effect coefficient\n")
        }
        Xmat <- model.matrix(formula, data=data)
        if (trait.type == "binary")
        {
            mod <- glm(formula, data=data, family=binomial)
        } else {
            mod <- glm(formula, data=data, family=gaussian)
        }
        covoffset <- Xmat[, -1L, drop=F] %*%  mod$coefficients[-1L]
        formula <- as.formula("y ~ 1")
    }

    # clear the internal GRM matrix and reinitialize
    .Call(saige_init_fit_grm)
    buf_sigma_diag <- double(n_samp)

    # reload SNP genotypes from the saved model
    if (verbose)
        cat("Reloading SNP genotypes from the saved model:\n")
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
    if (verbose)
    {
        .cat("    using ",
            .pretty_size(as.double(object.size(packed.geno))),
            " (stored in a ",
            ifelse(geno.sparse, "sparse", "dense"), " form)")
    }

    # reload GRM if available
    if (is.matrix(grm.mat))
    {
        .Call(saige_store_dense_grm, n_samp, grm.mat, buf_sigma_diag)
        if (verbose)
        {
            cat("User-defined genetic relationship matrix:\n")
            cat(sprintf("    %d x %d (dense matrix)\n", n_samp, n_samp))
        }
    } else if (inherits(grm.mat, "sparseMatrix"))
    {
        grm.mat <- .sp_to_dgCMatrix(grm.mat)
        .Call(saige_store_sparse_grm, n_samp, grm.mat, buf_sigma_diag)
        if (verbose)
        {
            a <- nnzero(grm.mat)
            s <- sprintf("%.3f%%", a/prod(dim(grm.mat))*100)
            if (s=="0.000%") s <- "<0.001%"
            cat("User-defined sparse genetic relationship matrix:\n")
            cat(sprintf("    %d x %d, # of nonzero: %d (%s)\n",
                n_samp, n_samp, a, s))
        }
    }

    # parameters for fitting the model
    param <- list(
        trait = match(trait.type, .trait_list),
        covoffset = covoffset, Xmat = Xmat,
        num.thread = num.thread, seed = seed,
        tol = tol, tolPCG = tolPCG,
        maxiter = maxiter, maxiterPCG = maxiterPCG,
        no_iteration = !isTRUE(tau.update),
        nrun = nrun, num.marker = num.marker,
        traceCVcutoff = traceCVcutoff, ratioCVcutoff = ratioCVcutoff,
        verbose = verbose,
        eventTime = evtime,
        indent = ""
    )

    # initialize tau for fitting the model
    tau.init <- model$tau
    tau.init[tau.init < 0] <- 0
    # gdsfile flag: TRUE indicates genotypes are loaded (for GRM diagonal)
    has.geno <- if (!is.null(packed.geno)) TRUE else NULL

    # fit the null model (AIREML) and build obj.noK
    if (trait.type == "binary")
    {
        glmm <- .fit_binary(verbose, X.transform, phenovar, data, formula,
            param, tau.init, has.geno, grm.mat, seed, n_var,
            calc_vr=isTRUE(recalcVR))
    } else if (trait.type == "quantitative")
    {
        glmm <- .fit_quant(verbose, X.transform, phenovar, data, inv.norm,
            formula, param, tau.init, has.geno, grm.mat, seed, n_var,
            calc_vr=isTRUE(recalcVR))
    } else if (trait.type == "survival")
    {
        glmm <- .fit_survival(verbose, X.transform, phenovar, data, formula,
            param, tau.init, has.geno, grm.mat, seed, n_var,
            calc_vr=isTRUE(recalcVR))
    } else
        stop("Invalid 'trait.type'.")

    # variance ratio
    if (!isTRUE(recalcVR))
        glmm$var.ratio <- model$var.ratio

    glmm <- c(list(formula=formula_str), glmm)
    glmm$use.cateMAC <- use.cateMAC

    # tweak the result
    if (!isTRUE(X.transform) || isTRUE(use.offset))
    {
        if (isTRUE(use.offset))
            names(glmm$coefficients) <- "(Offset)"
        else
            names(glmm$coefficients) <- colnames(glmm$obj.noK$X1)
    } else {
        coef <- solve(X_qrr, glmm$coefficients * sqrt(nrow(data)))
        names(coef) <- X_name
        glmm$coefficients <- coef
    }
    names(glmm$tau) <- c("Sigma_E", "Sigma_G")
    glmm$trait.type <- trait.type
    glmm$sample.id <- sid
    glmm$variant.id <- model$variant.id
    class(glmm) <- "ClassSAIGE_NullModel"

    if (verbose)
    {
        .cat(.crayon_underline(.tm()))
        .cat(.crayon_inverse("Done."))
    }

    glmm
}
