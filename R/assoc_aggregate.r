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
# Beta parameters in aggregate tests
#

AggrParamBeta <- structure(c(1,1,1,25), dim=c(2L,2L),
    dimnames=list(c("shape1", "shape2"), c("b_1_1", "b_1_25")))

.check_wbeta <- function(wbeta)
{
    stopifnot(is.numeric(wbeta), is.finite(wbeta))
    err <- "'wbeta' should be a length-two vector or a matrix with two rows."
    if (is.vector(wbeta))
    {
        if (length(wbeta) != 2L) stop(err)
    } else if (is.matrix(wbeta))
    {
        if (NCOL(wbeta) <= 0L) stop(err)
    }
    invisible()
}

.show_wbeta <- function(wbeta, verbose)
{
    if (verbose)
    {
        v <- apply(wbeta, 2L, function(x)
                paste0("beta(", x[1L], ",", x[2L], ")"))
        .cat("    variant weights: ", paste(v, collapse=", "))
    }
    invisible()
}

.set_check_sample_id <- function(gdsfile, modobj)
{
    seqSetFilter(gdsfile, sample.id=modobj$sample.id, warn=FALSE, verbose=FALSE)
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs are not available in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")
    ii
}

.set_show_units <- function(gdsfile, mobj, units, spa.pval, var.ratio, verbose)
{
    seqSetFilter(gdsfile, variant.sel=unlist(units$index), warn=FALSE,
        verbose=FALSE)
    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    if (dm[2L] <= 0)
        stop("No sample in the genotypic data set!")
    if (dm[3L] <= 0)
        stop("No variant in the genotypic data set!")
    v <- lengths(units$index, use.names=FALSE)
    sz_wmax <- max(v)
    if (verbose)
    {
        .cat("    # of samples: ", .pretty(dm[2L]))
        .cat("    # of variants in total: ", .pretty(dm[3L]))
        .cat("    trait: ", mobj$trait.type)
        .cat("Input sets of variants:")
        .cat("    # of units: ", .pretty(length(units$index)))
        .cat("    avg # of variants per unit: ", mean(v))
        .cat("    min # of variants in a unit: ", min(v))
        .cat("    max # of variants in a unit: ", sz_wmax)
        .cat("    sd  # of variants in a unit: ", sd(v))
        .cat("Parameters:")
        if (length(var.ratio) == 1L)
        {
            .cat("    variance ratio for approximation: ", var.ratio)
        } else {
            cat("    variance ratio for approximation (MAC categories):\n")
            v <- var.ratio
            attr(v, "cateMAC") <- NULL
            print(v, width=1024L)
        }
    }
    sz_wmax
}

.show_maf <- function(gdsfile, parallel)
{
    cat("Calculating minor allele frequencies (MAF):\n")
    maf <- seqAlleleFreq(gdsfile, minor=TRUE, verbose=TRUE, parallel=parallel)
    cat(sprintf(
        "    MAF: avg (%.5f), min (%.5f), max (%.5f), sd (%.5f)\n",
        mean(maf, na.rm=TRUE), min(maf, na.rm=TRUE),
        max(maf, na.rm=TRUE), sd(maf, na.rm=TRUE)))
    invisible()
}

.mapply <- function(mat_lst, irow)
{
    unlist(lapply(mat_lst, function(m) m[irow, ]))
}

.ncol <- function(x) if (length(d <- dim(x)) > 1L) d[2L] else 0L

.aggr_ret_obj <- function(units, obj, wbeta)
{
    # data.frame
    ans <- units$desp
    if (!is.null(ans$numvar)) ans$numvar <- NULL  # remove numvar
    ns <- vapply(obj, .ncol, 0L)
    ans <- lapply(ans, function(v) rep(v, times=ns))
    ans <- as.data.frame(ans)
    # set-based variables
    ans$maxMAF <- .mapply(obj, 1L)
    ans$numvar <- as.integer(.mapply(obj, 2L))
    ans$macmin <- .mapply(obj, 3L)
    ans$macmed <- .mapply(obj, 4L)
    ans$macmax <- .mapply(obj, 5L)
    ans$summac <- .mapply(obj, 6L)
    # weight beta
    v <- as.integer(.mapply(obj, 7L))
    attr(v, "levels") <- c(sprintf("(%g,%g)", wbeta[1L,], wbeta[2L,]), "Cauchy")
    attr(v, "class") <- "factor"
    ans$weight <- v
    # output
    ans
}

.aggr_ret_gds <- function(outf, gdsfile, units, obj, wbeta, Add)
{
    # add attributes for file format
    put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT_SET")
    put.attr.gdsn(outf$root, "Version",
        paste0("SAIGEgds_", packageVersion("SAIGEgds")))
    # add sample IDs
    Add("sample.id", seqGetData(gdsfile, "sample.id"))
    # data.frame
    ans <- units$desp
    if (!is.null(ans$numvar)) ans$numvar <- NULL  # remove numvar
    ns <- vapply(obj, .ncol, 0L)
    ans <- lapply(ans, function(v) rep(v, times=ns))
    for (nm in names(ans)) Add(nm, ans[[nm]])
    # set-based variables
    Add("maxMAF", .mapply(obj, 1L))
    Add("numvar", as.integer(.mapply(obj, 2L)))
    Add("macmin", .mapply(obj, 3L))
    Add("macmed", .mapply(obj, 4L))
    Add("macmax", .mapply(obj, 5L))
    Add("summac", .mapply(obj, 6L))
    # weight beta
    v <- as.integer(.mapply(obj, 7L))
    attr(v, "levels") <- c(sprintf("(%g,%g)", wbeta[1L,], wbeta[2L,]), "Cauchy")
    attr(v, "class") <- "factor"
    Add("weight", v)
    invisible()
}


#######################################################################
# SAIGE burden tests
#

# .check_modobj() plus a survival veto, for SKAT / ACAT-O. The GATE
# method (Bi et al. 2020) defines the Cox-via-Poisson score test + saddlepoint
# approximation only for single-variant statistics; this carries over to the
# burden test (a weighted collapse tested as one score) and to ACAT-V (a Cauchy
# combination of single-variant p-values), but NOT to the SKAT variance-
# component test, which has no published Poisson/Cox analogue. ACAT-O combines
# Burden + ACAT-V + SKAT, so it is unavailable for survival as well.
.check_modobj_no_surv <- function(modobj, verbose)
{
    modobj <- .check_modobj(modobj, verbose)
    if (!is.null(modobj$trait.type) && modobj$trait.type == "survival")
        stop("SKAT and ACAT-O do not support trait.type='survival': the GATE ",
            "method defines the score/saddlepoint test only for single-variant ",
            "and burden-style collapsing, and the SKAT variance-component test ",
            "has no Poisson/Cox analogue. Use seqAssocGLMM_Burden() or ",
            "seqAssocGLMM_ACAT_V() for survival outcomes.")
    modobj
}

.scan_fc_aggr_burden <- function(g, maxMAF)
    .Call(saige_burden_test_pval, g, maxMAF)

seqAssocGLMM_Burden <- function(gdsfile, modobj, units, maxMAF=0.01,
    wbeta=AggrParamBeta, missing=0.05, ccimb.adj=TRUE, ER.mac=4.5, dsnode="",
    geno.model=c("additive", "dominant", "recessive"),
    res.savefn="", res.compress="ZIP", parallel=FALSE,
    verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    stopifnot(is.numeric(maxMAF), 0<maxMAF & maxMAF<=1)
    .check_wbeta(wbeta)
    stopifnot(is.numeric(missing), length(missing)==1L)
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.numeric(ER.mac), length(ER.mac)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    geno.model <- match.arg(geno.model)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)
    if (verbose)
        .cat(.crayon_inverse("SAIGE burden analysis:"))

    # check model (survival supported: burden = single collapsed score test)
    modobj <- .check_modobj(modobj, verbose)
    var.ratio <- .get_var_ratio(modobj)
    if (!length(maxMAF)) maxMAF <- 1
    maxMAF <- sort(maxMAF, decreasing=TRUE)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest
    if (!isTRUE(ccimb.adj)) spa.pval <- -1
    if (is.na(ER.mac)) ER.mac <- 0

    # GDS file
    if (is.character(gdsfile))
    {
        if (verbose)
            .cat("    open ", sQuote(gdsfile))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        # save the filter on GDS file
        seqFilterPush(gdsfile)
        on.exit(seqFilterPop(gdsfile))
    }

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode, TRUE)

    # check sample ID
    ii <- .set_check_sample_id(gdsfile, modobj)

    # set variant filter and show summary
    sz_wmax <- .set_show_units(gdsfile, modobj, units, spa.pval, var.ratio,
        verbose)

    # show beta weights
    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    .show_wbeta(wbeta, verbose)
    if (verbose)
    {
        .cat("    MAF threshold", ifelse(length(maxMAF)>1L, "s", ""), ": ",
            paste(maxMAF, collapse=", "))
        .cat("    missing proportion threshold: ", .pretty_lt_eq(missing))
        .cat("    genetic model: ", geno.model)
    }

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # show allele frequencies
    if (verbose && isTRUE(verbose.maf))
        .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, missing, spa.pval, ER.mac,
        var.ratio, 2L, modobj$Sigma_inv, modobj$chol_inv_X_Sigma, maxMAF,
        wbeta, sz_wmax, geno.model=geno.model)

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(.load_lib)
                .PkgEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() { .PkgEnv$mobj <- NULL })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=.scan_fc_aggr_burden,
        as.is="list", parallel=parallel, .useraw=NA, .progress=verbose,
        maxMAF=maxMAF)

    # check
    if (length(rv) != length(units$index))
        stop("seqUnitApply() returns a vector of wrong length.")

    # output to a GDS file?
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE))
    {
        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=res.compress[1L], closezip=TRUE)
        # create a GDS file
        outf <- createfn.gds(res.savefn)
        on.exit(closefn.gds(outf), add=TRUE)
        .aggr_ret_gds(outf, gdsfile, units, rv, wbeta, Add)
        Add("beta", .mapply(rv, 8L))
        Add("SE", .mapply(rv, 9L))
        Add("pval", .mapply(rv, 10L))
        Add("method", .pval_method(.mapply(rv, 11L)))
        if (modobj$trait.type %in% c("binary", "survival"))
        {
            Add("p.norm", .mapply(rv, 12L))
            Add("converged", .mapply(rv, 13L)==1L)
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output nothing
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv, wbeta)
        ans$beta <- .mapply(rv, 8L)
        ans$SE <- .mapply(rv, 9L)
        ans$pval <- .mapply(rv, 10L)
        ans$method <- .pval_method(.mapply(rv, 11L))
        if (modobj$trait.type %in% c("binary", "survival"))
        {
            ans$p.norm <- .mapply(rv, 12L)
            ans$converged <- .mapply(rv, 13L)==1L
        }
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}



#######################################################################
# SAIGE SKAT tests
#

# initialize loading SKAT functions
.skat_init_func <- function()
{
    .Call("saige_init_skat_pkg", PACKAGE="SAIGEgds")
}

# get p-value from mixed Chi square distribution
#   (try CompQuadForm first, then saddle point method, suggested by UW DCC)
# used in saige_main.cpp
.skat_eig_chiq <- function(Q, eigval)
{
    # try integration method
    v <- suppressWarnings(CompQuadForm::davies(Q, eigval, acc=1e-9))
    if ((v$ifault > 0L) || (v$Qq < 1e3*.Machine$double.eps) || (v$Qq > 1))
    {
        # try saddlepoint method
        survey:::saddle(Q, eigval)
    } else {
        v$Qq
    }
}

.scan_fc_aggr_skat <- function(g) .Call(saige_skat_test_pval, g)

# SKAT tests
seqAssocGLMM_SKAT <- function(gdsfile, modobj, units, maxMAF=0.01,
    wbeta=AggrParamBeta, missing=0.05, collapse.mac=10,
    collapse.method=c("max", "sum"), ccimb.adj=TRUE, ER.mac=4.5, dsnode="",
    geno.model=c("additive", "dominant", "recessive"),
    res.savefn="", res.compress="ZIP", parallel=FALSE,
    verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    stopifnot(is.numeric(maxMAF), 0<maxMAF & maxMAF<=1)
    .check_wbeta(wbeta)
    stopifnot(is.numeric(missing), length(missing)==1L)
    stopifnot(is.numeric(collapse.mac), length(collapse.mac)==1L,
        is.finite(collapse.mac))
    collapse.method <- match.arg(collapse.method)
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.numeric(ER.mac), length(ER.mac)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    geno.model <- match.arg(geno.model)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    # check packages
    pkg_cqf <- suppressPackageStartupMessages(
        requireNamespace("CompQuadForm", quietly=TRUE))
    pkg_svy <- suppressPackageStartupMessages(
        requireNamespace("survey", quietly=TRUE))
    if (!pkg_cqf || !pkg_svy)
        stop("The packages 'CompQuadForm' and 'survey' should be installed.")

    if (verbose)
        .cat(.crayon_inverse("SAIGE SKAT analysis:"))

    # check model
    modobj <- .check_modobj_no_surv(modobj, verbose)
    if (is.null(modobj$Sigma_inv) || is.null(modobj$chol_inv_X_Sigma))
    {
        stop("A (sparse) genetic relationship matrix 'grm.mat' should be ",
            "specified in seqFitNullGLMM_SPA(), when the null model is built ",
            "for SKAT.")
    }

    # variance ratio
    var.ratio <- .get_var_ratio(modobj)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest
    if (!isTRUE(ccimb.adj)) spa.pval <- -1
    if (is.na(ER.mac)) ER.mac <- 0

    # GDS file
    if (is.character(gdsfile))
    {
        if (verbose)
            .cat("    open ", sQuote(gdsfile))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        # save the filter on GDS file
        seqSetFilter(gdsfile, action="push", verbose=FALSE)
        on.exit(seqSetFilter(gdsfile, action="pop", verbose=FALSE))
    }

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode, TRUE)

    # check sample ID
    ii <- .set_check_sample_id(gdsfile, modobj)

    # set variant filter and show summary
    sz_wmax <- .set_show_units(gdsfile, modobj, units, spa.pval, var.ratio,
        verbose)
    if (verbose)
    {
        .cat("    MAC threshold for collapsing ultra rare variants: <= ",
            sprintf("%.15g", collapse.mac))
        if (is.finite(collapse.mac) && collapse.mac>0)
        {
            if (collapse.method == "max")
                cat("        max: maximum of dosages of rare variants\n")
            else if (collapse.method == "sum")
                cat("        sum: sum up rare genotypes\n")
        }
        if (modobj$trait.type == "binary")
        {
            if (isTRUE(ccimb.adj))
                cat("    accounting for case-control imbalance\n")
            else
                cat("    not accounting for case-control imbalance\n")
        }
    }

    # show beta weights
    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    .show_wbeta(wbeta, verbose)
    if (verbose)
    {
        .cat("    MAF threshold", ifelse(length(maxMAF)>1L, "s", ""), ": ",
            paste(maxMAF, collapse=", "))
        .cat("    missing proportion threshold: ", .pretty_lt_eq(missing))
        .cat("    genetic model: ", geno.model)
    }

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, missing, spa.pval, ER.mac,
        var.ratio, 2L, modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        maxMAF, wbeta, sz_wmax, collapse_method=collapse.method,
        geno.model=geno.model)
    # load package(s)
    .load_skat(FALSE)

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
        # initialize SKAT
        mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
        .Call(saige_skat_test_init, mobj$Sigma_inv_cg, mobj$t_XVX_inv_XV,
            mobj$Si_X, mobj$XVX_inv_XV_X_Si_X)
        # finalize
        on.exit(.Call(saige_skat_test_done), add=TRUE)
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(.load_lib)
                mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
                .PkgEnv$mobj <- mobj
                # initialize SKAT
                .Call(saige_score_test_init, mobj)
                .Call(saige_skat_test_init, mobj$Sigma_inv_cg,
                    mobj$t_XVX_inv_XV, mobj$Si_X, mobj$XVX_inv_XV_X_Si_X)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() {
                    .PkgEnv$mobj <- NULL
                    .Call(saige_skat_test_done)
                })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=.scan_fc_aggr_skat, as.is="list",
        parallel=parallel, .useraw=NA, .progress=verbose)
    # check
    if (length(rv) != length(units$index))
        stop("seqUnitApply() returns a vector of wrong length.")

    # output to a GDS file?
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE))
    {
        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        cm <- res.compress[1L]
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # create a GDS file
        outf <- createfn.gds(res.savefn)
        on.exit(closefn.gds(outf), add=TRUE)
        .aggr_ret_gds(outf, gdsfile, units, rv, wbeta, Add)
        Add("n_collapse", as.integer(.mapply(rv, 8L)))
        Add("g_ncol", as.integer(.mapply(rv, 9L)))
        Add("g_minMAC", .mapply(rv, 10L))
        Add("pval", .mapply(rv, 11L))
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv, wbeta)
        ans$n_collapse <- as.integer(.mapply(rv, 8L))
        ans$g_ncol <- as.integer(.mapply(rv, 9L))
        ans$g_minMAC <- .mapply(rv, 10L)
        ans$pval <- .mapply(rv, 11L)
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}



#######################################################################
# SAIGE ACAT-V tests
#

.scan_fc_aggr_acatv <- function(g) .Call(saige_acatv_test_pval, g)

seqAssocGLMM_ACAT_V <- function(gdsfile, modobj, units, maxMAF=0.01,
    wbeta=AggrParamBeta, missing=0.05, ccimb.adj=TRUE, collapse.mac=10,
    ER.mac=4.5, dsnode="",
    geno.model=c("additive", "dominant", "recessive"),
    res.savefn="", res.compress="ZIP", parallel=FALSE,
    verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    stopifnot(is.numeric(maxMAF), 0<maxMAF & maxMAF<=1)
    .check_wbeta(wbeta)
    stopifnot(is.numeric(missing), length(missing)==1L)
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.numeric(collapse.mac), length(collapse.mac)==1L,
        is.finite(collapse.mac))
    stopifnot(is.numeric(ER.mac), length(ER.mac)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    geno.model <- match.arg(geno.model)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)
    if (verbose)
        .cat(.crayon_inverse("SAIGE ACAT-V analysis:"))

    # check model (survival supported: ACAT-V = Cauchy combination of p-values)
    modobj <- .check_modobj(modobj, verbose)
    var.ratio <- .get_var_ratio(modobj)
    if (!length(maxMAF)) maxMAF <- 1
    maxMAF <- sort(maxMAF, decreasing=TRUE)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest
    if (!isTRUE(ccimb.adj)) spa.pval <- -1
    if (is.na(ER.mac)) ER.mac <- 0

    # GDS file
    if (is.character(gdsfile))
    {
        if (verbose)
            .cat("    open ", sQuote(gdsfile))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        # save the filter on GDS file
        seqFilterPush(gdsfile)
        on.exit(seqFilterPop(gdsfile))
    }

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode, TRUE)

    # check sample ID
    ii <- .set_check_sample_id(gdsfile, modobj)

    # set variant filter and show summary
    sz_wmax <- .set_show_units(gdsfile, modobj, units, spa.pval, var.ratio,
        verbose)
    if (verbose)
    {
        .cat("    MAC threshold for collapsing ultra rare variants: <= ",
            sprintf("%.15g", collapse.mac))
    }

    # show beta weights
    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    .show_wbeta(wbeta, verbose)
    if (verbose)
    {
        .cat("    MAF threshold", ifelse(length(maxMAF)>1L, "s", ""), ": ",
            paste(maxMAF, collapse=", "))
        .cat("    missing proportion threshold: ", .pretty_lt_eq(missing))
        .cat("    genetic model: ", geno.model)
    }

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, missing, spa.pval, ER.mac,
        var.ratio, 2L, modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        maxMAF, wbeta, sz_wmax, collapse.mac, geno.model=geno.model)

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(.load_lib)
                .PkgEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() { .PkgEnv$mobj <- NULL })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=.scan_fc_aggr_acatv, as.is="list",
        parallel=parallel, .useraw=NA, .progress=verbose)
    # check
    if (length(rv) != length(units$index))
        stop("seqUnitApply() returns a vector of wrong length.")

    # output to a GDS file?
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE))
    {
        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        cm <- res.compress[1L]
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # create a GDS file
        outf <- createfn.gds(res.savefn)
        on.exit(closefn.gds(outf), add=TRUE)
        .aggr_ret_gds(outf, gdsfile, units, rv, wbeta, Add)
        Add("n_single", as.integer(.mapply(rv, 8L)))
        Add("n_collapse", as.integer(.mapply(rv, 9L)))
        Add("pval", .mapply(rv, 10L))
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv, wbeta)
        ans$n_single <- as.integer(.mapply(rv, 8L))
        ans$n_collapse <- as.integer(.mapply(rv, 9L))
        ans$pval <- .mapply(rv, 10L)
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}



#######################################################################
# SAIGE ACAT-O tests
#

.scan_fc_aggr_acato <- function(g) .Call(saige_acato_test_pval, g)

seqAssocGLMM_ACAT_O <- function(gdsfile, modobj, units, maxMAF=0.01,
    wbeta=AggrParamBeta, missing=0.05, collapse.mac=10,
    collapse.method=c("max", "sum"), ccimb.adj=TRUE, ER.mac=4.5, dsnode="",
    geno.model=c("additive", "dominant", "recessive"),
    res.savefn="", res.compress="ZIP", parallel=FALSE,
    verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    stopifnot(is.numeric(maxMAF), 0<maxMAF & maxMAF<=1)
    .check_wbeta(wbeta)
    stopifnot(is.numeric(missing), length(missing)==1L)
    stopifnot(is.numeric(collapse.mac), length(collapse.mac)==1L,
        is.finite(collapse.mac))
    collapse.method <- match.arg(collapse.method)
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.numeric(ER.mac), length(ER.mac)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    geno.model <- match.arg(geno.model)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    # check packages
    pkg_cqf <- suppressPackageStartupMessages(
        requireNamespace("CompQuadForm", quietly=TRUE))
    pkg_svy <- suppressPackageStartupMessages(
        requireNamespace("survey", quietly=TRUE))
    if (!pkg_cqf || !pkg_svy)
        stop("The packages 'CompQuadForm' and 'survey' should be installed.")

    if (verbose)
        .cat(.crayon_inverse("SAIGE ACAT-O analysis:"))

    # check model
    modobj <- .check_modobj_no_surv(modobj, verbose)
    if (is.null(modobj$Sigma_inv) || is.null(modobj$chol_inv_X_Sigma))
    {
        stop("A (sparse) genetic relationship matrix 'grm.mat' should be ",
            "specified in seqFitNullGLMM_SPA(), when the null model is built ",
            "for SKAT.")
    }

    # variance ratio
    var.ratio <- .get_var_ratio(modobj)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest
    if (!isTRUE(ccimb.adj)) spa.pval <- -1
    if (is.na(ER.mac)) ER.mac <- 0

    # GDS file
    if (is.character(gdsfile))
    {
        if (verbose)
            .cat("    open ", sQuote(gdsfile))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    } else {
        # save the filter on GDS file
        seqSetFilter(gdsfile, action="push", verbose=FALSE)
        on.exit(seqSetFilter(gdsfile, action="pop", verbose=FALSE))
    }

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode, TRUE)

    # check sample ID
    ii <- .set_check_sample_id(gdsfile, modobj)

    # set variant filter and show summary
    sz_wmax <- .set_show_units(gdsfile, modobj, units, spa.pval, var.ratio,
        verbose)
    if (verbose)
    {
        .cat("    MAC threshold for collapsing ultra rare variants in ",
            "ACAT-V and SKAT: <= ", sprintf("%.15g", collapse.mac))
        cat("    ACAT-O p-values combine Burden, ACAT-V and SKAT\n")
        if (modobj$trait.type == "binary")
        {
            if (isTRUE(ccimb.adj))
                cat("    accounting for case-control imbalance\n")
            else
                cat("    not accounting for case-control imbalance\n")
        }
    }

    # show beta weights
    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    .show_wbeta(wbeta, verbose)
    if (verbose)
    {
        .cat("    MAF threshold", ifelse(length(maxMAF)>1L, "s", ""), ": ",
            paste(maxMAF, collapse=", "))
        .cat("    missing proportion threshold: ", .pretty_lt_eq(missing))
        .cat("    genetic model: ", geno.model)
    }

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, missing, spa.pval, ER.mac,
        var.ratio, 2L, modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        maxMAF, wbeta, sz_wmax, collapse.mac, geno.model=geno.model)
    i <- match(collapse.method, c("max", "sum"))
    if (is.na(i)) stop("Internal error in 'collapse.method'.")
    mobj$collapse.method <- i
    # load package(s)
    .load_skat(FALSE)

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
        # initialize SKAT
        mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
        .Call(saige_skat_test_init, mobj$Sigma_inv_cg, mobj$t_XVX_inv_XV,
            mobj$Si_X, mobj$XVX_inv_XV_X_Si_X)
        # finalize
        on.exit(.Call(saige_skat_test_done), add=TRUE)
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(.load_lib)
                mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
                .PkgEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
                if (!is.null(mobj$Sigma_inv))
                {
                    # initialize SKAT
                    .Call(saige_skat_test_init, mobj$Sigma_inv_cg,
                        mobj$t_XVX_inv_XV, mobj$Si_X, mobj$XVX_inv_XV_X_Si_X)
                } else {
                    .Call(saige_skat_test_reset)
                }
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() {
                    .PkgEnv$mobj <- NULL
                    .Call(saige_skat_test_done)
                })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=.scan_fc_aggr_acato, as.is="list",
        parallel=parallel, .useraw=NA, .progress=verbose)
    # check
    if (length(rv) != length(units$index))
        stop("seqUnitApply() returns a vector of wrong length.")

    # output to a GDS file?
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE))
    {
        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        cm <- res.compress[1L]
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE)
        # create a GDS file
        outf <- createfn.gds(res.savefn)
        on.exit(closefn.gds(outf), add=TRUE)
        .aggr_ret_gds(outf, gdsfile, units, rv, wbeta, Add)
        Add("n_collapse", as.integer(.mapply(rv, 8L)))
        Add("pval", .mapply(rv, 9L))
        Add("p.burden", .mapply(rv, 10L))
        Add("p.skat", .mapply(rv, 11L))
        Add("p.acatv", .mapply(rv, 12L))
        Add("burden.beta", .mapply(rv, 13L))
        Add("burden.se", .mapply(rv, 14L))
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv, wbeta)
        ans$n_collapse <- as.integer(.mapply(rv, 8L))
        ans$pval <- .mapply(rv, 9L)
        ans$p.burden <- .mapply(rv, 10L)
        ans$p.skat <- .mapply(rv, 11L)
        ans$p.acatv <- .mapply(rv, 12L)
        ans$burden.beta <- .mapply(rv, 13L)
        ans$burden.se <- .mapply(rv, 14L)
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}



#######################################################################
# SAIGE ACAT-O tests on in-memory genotype matrices
#

seqAssocGLMM_ACAT_O_GT <- function(gt_lst, modobj, maxMAF=0.01,
    wbeta=AggrParamBeta, missing=0.05, collapse.mac=10,
    collapse.method=c("max", "sum"), ccimb.adj=TRUE, ER.mac=4.5,
    geno.model=c("additive", "dominant", "recessive"),
    verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(is.list(gt_lst), length(gt_lst) > 0L)
    stopifnot(is.numeric(maxMAF), 0<maxMAF & maxMAF<=1)
    .check_wbeta(wbeta)
    stopifnot(is.numeric(missing), length(missing)==1L)
    stopifnot(is.numeric(collapse.mac), length(collapse.mac)==1L,
        is.finite(collapse.mac))
    collapse.method <- match.arg(collapse.method)
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.numeric(ER.mac), length(ER.mac)==1L)
    geno.model <- match.arg(geno.model)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    # check packages
    pkg_cqf <- suppressPackageStartupMessages(
        requireNamespace("CompQuadForm", quietly=TRUE))
    pkg_svy <- suppressPackageStartupMessages(
        requireNamespace("survey", quietly=TRUE))
    if (!pkg_cqf || !pkg_svy)
        stop("The packages 'CompQuadForm' and 'survey' should be installed.")

    if (verbose)
        .cat(.crayon_inverse("SAIGE ACAT-O analysis (in-memory):"))

    # check model
    modobj <- .check_modobj_no_surv(modobj, verbose)
    if (is.null(modobj$Sigma_inv) || is.null(modobj$chol_inv_X_Sigma))
    {
        stop("A (sparse) genetic relationship matrix 'grm.mat' should be ",
            "specified in seqFitNullGLMM_SPA(), when the null model is built ",
            "for SKAT.")
    }

    # check gt_lst: each element should be a numeric matrix (sample x variant)
    n <- length(modobj$sample.id)
    nv <- integer(length(gt_lst))
    for (j in seq_along(gt_lst))
    {
        g <- gt_lst[[j]]
        if (!(is.numeric(g) && is.matrix(g)) && !inherits(g, "Matrix"))
            stop("gt_lst[[", j, "]] should be a numeric matrix or Matrix.")
        if (nrow(g) != n)
            stop("nrow(gt_lst[[", j, "]]) should be ", n, ".")
        nv[j] <- ncol(g)
    }
    sz_wmax <- max(nv)

    # variance ratio
    var.ratio <- .get_var_ratio(modobj)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest
    if (!isTRUE(ccimb.adj)) spa.pval <- -1
    if (is.na(ER.mac)) ER.mac <- 0

    if (verbose)
    {
        .cat("    # of samples: ", .pretty(n))
        .cat("    # of units: ", .pretty(length(gt_lst)))
        if (length(nv) > 1L)
        {
            .cat("    avg # of variants per unit: ", mean(nv))
            .cat("    min # of variants in a unit: ", min(nv))
            .cat("    max # of variants in a unit: ", sz_wmax)
            .cat("    sd  # of variants in a unit: ", sd(nv))
        } else {
            .cat("    # of variants per unit: ", sz_wmax)
        }
        .cat("    trait: ", modobj$trait.type)
        .cat("    MAC threshold for collapsing ultra rare variants in ",
            "ACAT-V and SKAT: <= ", sprintf("%.15g", collapse.mac))
        cat("    ACAT-O p-values combine Burden, ACAT-V and SKAT\n")
        if (modobj$trait.type == "binary")
        {
            if (isTRUE(ccimb.adj))
                cat("    accounting for case-control imbalance\n")
            else
                cat("    not accounting for case-control imbalance\n")
        }
    }

    # show beta weights
    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    .show_wbeta(wbeta, verbose)
    if (verbose)
    {
        .cat("    MAF threshold", ifelse(length(maxMAF)>1L, "s", ""), ": ",
            paste(maxMAF, collapse=", "))
        .cat("    missing proportion threshold: ", .pretty_lt_eq(missing))
        .cat("    genetic model: ", geno.model)
    }

    # initialize the internal model parameters
    ii <- TRUE
    mobj <- .init_nullmod(modobj, ii, 0, 0, missing, spa.pval, ER.mac,
        var.ratio, 2L, modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        maxMAF, wbeta, sz_wmax, collapse.mac, geno.model=geno.model)
    i <- match(collapse.method, c("max", "sum"))
    if (is.na(i)) stop("Internal error in 'collapse.method'.")
    mobj$collapse.method <- i
    # load package(s)
    .load_skat(FALSE)

    # initialize internally
    .Call(saige_score_test_init, mobj)
    # initialize SKAT
    mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
    .Call(saige_skat_test_init, mobj$Sigma_inv_cg, mobj$t_XVX_inv_XV,
        mobj$Si_X, mobj$XVX_inv_XV_X_Si_X)
    # finalize
    on.exit(.Call(saige_skat_test_done), add=TRUE)

    # scan all units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- vector("list", length(gt_lst))
    for (j in seq_along(gt_lst))
    {
        # convert numeric matrix to dgCMatrix for the C code
        g <- gt_lst[[j]]
        if (!inherits(g, "dgCMatrix"))
            g <- methods::as(gt_lst[[j]], "CsparseMatrix")
        rv[[j]] <- .Call(saige_acato_test_pval, g)
    }

    # output
    ns <- vapply(rv, .ncol, 0L)
    # build id column
    if (!is.null(names(gt_lst)))
    {
        id <- rep(names(gt_lst), times=ns)
    } else {
        id <- rep(seq_along(gt_lst), times=ns)
    }
    ans <- data.frame(id=id, stringsAsFactors=FALSE)
    ans$maxMAF <- .mapply(rv, 1L)
    ans$numvar <- as.integer(.mapply(rv, 2L))
    ans$macmin <- .mapply(rv, 3L)
    ans$macmed <- .mapply(rv, 4L)
    ans$macmax <- .mapply(rv, 5L)
    ans$summac <- .mapply(rv, 6L)
    # weight beta
    v <- as.integer(.mapply(rv, 7L))
    attr(v, "levels") <- c(sprintf("(%g,%g)", wbeta[1L,], wbeta[2L,]), "Cauchy")
    attr(v, "class") <- "factor"
    ans$weight <- v
    ans$n_collapse <- as.integer(.mapply(rv, 8L))
    ans$pval <- .mapply(rv, 9L)
    ans$p.burden <- .mapply(rv, 10L)
    ans$p.skat <- .mapply(rv, 11L)
    ans$p.acatv <- .mapply(rv, 12L)
    ans$burden.beta <- .mapply(rv, 13L)
    ans$burden.se <- .mapply(rv, 14L)
    ans
}
