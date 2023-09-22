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

.aggr_ret_obj <- function(units, obj)
{
    ans <- units$desp  # it is a data.frame
    ans$numvar  <- as.integer(vapply(obj, `[`, 0, i=1L))
    ans$maf.avg <- vapply(obj, `[`, 0, i=2L)
    ans$maf.sd  <- vapply(obj, `[`, 0, i=3L)
    ans$maf.min <- vapply(obj, `[`, 0, i=4L)
    ans$maf.max <- vapply(obj, `[`, 0, i=5L)
    ans$mac.avg <- vapply(obj, `[`, 0, i=6L)
    ans$mac.sd  <- vapply(obj, `[`, 0, i=7L)
    ans$mac.min <- vapply(obj, `[`, 0, i=8L)
    ans$mac.max <- vapply(obj, `[`, 0, i=9L)
    ans
}

.aggr_ret_gds <- function(outf, gdsfile, units, obj, Add)
{
    put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT_SET")
    put.attr.gdsn(outf$root, "Version",
        paste0("SAIGEgds_", packageVersion("SAIGEgds")))
    # add sample IDs
    Add("sample.id", seqGetData(gdsfile, "sample.id"))
    # write summary variant data
    for (nm in names(units$desp))
        Add(nm, units$desp[[nm]])
    Add("numvar", as.integer(vapply(obj, `[`, 0, i=1L)))
    Add("maf.avg", vapply(obj, `[`, 0, i=2L))
    Add("maf.sd",  vapply(obj, `[`, 0, i=3L))
    Add("maf.min", vapply(obj, `[`, 0, i=4L))
    Add("maf.max", vapply(obj, `[`, 0, i=5L))
    Add("mac.avg", vapply(obj, `[`, 0, i=6L))
    Add("mac.sd",  vapply(obj, `[`, 0, i=7L))
    Add("mac.min", vapply(obj, `[`, 0, i=8L))
    Add("mac.max", vapply(obj, `[`, 0, i=9L))
    invisible()
}


#######################################################################
# SAIGE burden tests
#

seqAssocGLMM_Burden <- function(gdsfile, modobj, units, wbeta=AggrParamBeta,
    ccimb.adj=TRUE, summac=3, dsnode="", res.savefn="", res.compress="LZMA",
    parallel=FALSE, verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    .check_wbeta(wbeta)
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.numeric(summac), length(summac)==1L, is.finite(summac))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    if (verbose)
        .cat(.crayon_inverse("SAIGE burden analysis:"))

    # check model
    modobj <- .check_modobj(modobj, verbose)
    var.ratio <- .get_var_ratio(modobj)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest
    if (!isTRUE(ccimb.adj)) spa.pval <- -1

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
    wb_colnm <- sprintf("%g_%g", wbeta[1L,], wbeta[2L,])
    .show_wbeta(wbeta, verbose)

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
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio, 2L,
        modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        summac, wbeta, sz_wmax)

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
                .packageEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() { .packageEnv$mobj <- NULL })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=function(x) .Call(saige_burden_test_pval, x), as.is="list",
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
        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=res.compress[1L], closezip=TRUE)
        # create a GDS file
        outf <- createfn.gds(res.savefn)
        on.exit(closefn.gds(outf), add=TRUE)
        .aggr_ret_gds(outf, gdsfile, units, rv, Add)
        Add("summac", vapply(rv, `[`, 0, i=10L))
        # write p-values
        k <- 10L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            Add(paste0("beta", s), vapply(rv, `[`, 0, i=k+1L))
            Add(paste0("SE", s),   vapply(rv, `[`, 0, i=k+2L))
            Add(paste0("pval", s), vapply(rv, `[`, 0, i=k+3L))
            k <- k + 3L
            if (modobj$trait.type == "binary")
            {
                Add(paste0("p.norm", s), vapply(rv, `[`, 0, i=k+1L))
                Add(paste0("cvg", s), as.logical(vapply(rv, `[`, 0, i=k+2L)))
                k <- k + 2L
            }
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output nothing
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv)
        ans[["summac"]] <- vapply(rv, `[`, 0, i=10L)
        k <- 10L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            ans[[paste0("beta", s)]] <- vapply(rv, `[`, 0, i=k+1L)
            ans[[paste0("SE", s)]]   <- vapply(rv, `[`, 0, i=k+2L)
            ans[[paste0("pval", s)]] <- vapply(rv, `[`, 0, i=k+3L)
            k <- k + 3L
            if (modobj$trait.type == "binary")
            {
                ans[[paste0("p.norm", s)]] <- vapply(rv, `[`, 0, i=k+1L)
                ans[[paste0("cvg", s)]] <- as.logical(vapply(rv, `[`, 0, i=k+2L))
                k <- k + 2L
            }
        }
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}



#######################################################################
# SAIGE SKAT tests
#

# get p-value from mixed Chi square distribution
#   (try CompQuadForm first, then saddle point method, suggested by UW DCC)
.skat_eig_chiq <- function(Q, eigval)
{
    # try integration method
    # print(head(eigval)); print(summary(eigval))
    v <- suppressWarnings(CompQuadForm::davies(Q, eigval, acc=1e-9))
    if ((v$ifault > 0L) || (v$Qq < 1e3*.Machine$double.eps) || (v$Qq > 1))
    {
        # try saddlepoint method
        survey:::saddle(Q, eigval)
    } else {
        v$Qq
    }
}

# SKAT tests
seqAssocGLMM_SKAT <- function(gdsfile, modobj, units, wbeta=AggrParamBeta,
    collapse.mac=10, collapse.method=c("PA", "PA_int", "SumG"), dsnode="",
    ccimb.adj=TRUE, res.savefn="", res.compress="LZMA", parallel=FALSE,
    verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    .check_wbeta(wbeta)
    stopifnot(is.numeric(collapse.mac), length(collapse.mac)==1L,
        is.finite(collapse.mac))
    collapse.method <- match.arg(collapse.method)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    # check packages
    pkg_cqf <- suppressPackageStartupMessages(requireNamespace("CompQuadForm",
        quietly=TRUE))
    pkg_svy <- suppressPackageStartupMessages(requireNamespace("survey",
        quietly=TRUE))
    if (!pkg_cqf || !pkg_svy)
        stop("The packages 'CompQuadForm' and 'survey' should be installed.")

    if (verbose)
        .cat(.crayon_inverse("SAIGE SKAT analysis:"))

    # check model
    modobj <- .check_modobj(modobj, verbose)
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
            if (collapse.method == "PA")
                cat("        PA: presence or absence using dosage 1, 2 or (max if < 0.5)\n")
            else if (collapse.method == "PA_int")
                cat("        PA_int: presence or absence using integer dosages (0, 1, 2)\n")
            else if (collapse.method == "SumG")
                cat("        SumG: sum up rare genotypes\n")
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
    wb_colnm <- sprintf("%g_%g", wbeta[1L,], wbeta[2L,])
    .show_wbeta(wbeta, verbose)

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio, 2L,
        modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        0, wbeta, sz_wmax, skat_mac=collapse.mac)
    i <- match(collapse.method, c("PA", "PA_int", "SumG"))
    if (is.na(i)) stop("Internal error in 'collapse.method'.")
    mobj$collapse.method <- i

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
        # initialize SKAT
        mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
        .Call(saige_skat_test_init, mobj$Sigma_inv_cg, mobj$t_XVX_inv_XV,
            mobj$Si_X, mobj$XVX_inv_XV_X_Si_X, mobj$collapse.method)
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
                .packageEnv$mobj <- mobj
                # initialize SKAT
                .Call(saige_score_test_init, mobj)
                .Call(saige_skat_test_init, mobj$Sigma_inv_cg, mobj$t_XVX_inv_XV,
                    mobj$Si_X, mobj$XVX_inv_XV_X_Si_X, mobj$collapse.method)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() {
                    .packageEnv$mobj <- NULL
                    .Call(saige_skat_test_done)
                })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=function(x) .Call(saige_skat_test_pval, x), as.is="list",
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
        .aggr_ret_gds(outf, gdsfile, units, rv, Add)
        Add("n_collapse", as.integer(vapply(rv, `[`, 0, i=10L)))
        Add("g_ncol", as.integer(vapply(rv, `[`, 0, i=11L)))
        Add("g_minMAC", vapply(rv, `[`, 0, i=12L))
        # write p-values
        st <- 12L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            Add(paste0("pval", s),  vapply(rv, `[`, 0, i=st+1L))
            st <- st + 1L
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv)
        ans$n_collapse <- as.integer(vapply(rv, `[`, 0, i=10L))
        ans$g_ncol <- as.integer(vapply(rv, `[`, 0, i=11L))
        ans$g_minMAC <- vapply(rv, `[`, 0, i=12L)
        st <- 12L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            ans[[paste0("pval", s)]]  <- vapply(rv, `[`, 0, i=st+1L)
            st <- st + 1L
        }
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}



#######################################################################
# SAIGE ACAT-V tests
#

seqAssocGLMM_ACAT_V <- function(gdsfile, modobj, units, wbeta=AggrParamBeta,
    ccimb.adj=TRUE, collapse.mac=10, burden.summac=3, dsnode="", res.savefn="",
    res.compress="LZMA", parallel=FALSE, verbose=TRUE, verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    .check_wbeta(wbeta)
    stopifnot(is.logical(ccimb.adj), length(ccimb.adj)==1L)
    stopifnot(is.numeric(collapse.mac), length(collapse.mac)==1L,
        is.finite(collapse.mac))
    stopifnot(is.numeric(burden.summac), length(burden.summac)==1L,
        is.finite(burden.summac))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    if (verbose)
        .cat(.crayon_inverse("SAIGE ACAT-V analysis:"))

    # check model
    modobj <- .check_modobj(modobj, verbose)
    var.ratio <- .get_var_ratio(modobj)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest
    if (!isTRUE(ccimb.adj)) spa.pval <- -1

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
    wb_colnm <- sprintf("%g_%g", wbeta[1L,], wbeta[2L,])
    .show_wbeta(wbeta, verbose)

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio, 2L,
        modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        burden.summac, wbeta, sz_wmax, collapse.mac)

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
                .packageEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() { .packageEnv$mobj <- NULL })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=function(x) .Call(saige_acatv_test_pval, x), as.is="list",
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
        .aggr_ret_gds(outf, gdsfile, units, rv, Add)
        Add("n_single", as.integer(sapply(rv, `[`, i=9L)))
        Add("n_collapse", as.integer(sapply(rv, `[`, i=10L)))
        # write p-values
        st <- 10L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            Add(paste0("pval", s),  sapply(rv, `[`, i=st+1L))
            Add(paste0("p.med", s), sapply(rv, `[`, i=st+2L))
            Add(paste0("p.min", s), sapply(rv, `[`, i=st+3L))
            Add(paste0("p.max", s), sapply(rv, `[`, i=st+4L))
            st <- st + 4L
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv)
        ans$n.single <- as.integer(sapply(rv, `[`, i=10L))
        ans$n.burden <- as.integer(sapply(rv, `[`, i=11L))
        st <- 11L
        for (i in seq_len(ncol(wbeta)))
        {
            s <- ""
            if (length(wb_colnm) > 1L)
                s <- paste0(".", wb_colnm[i])
            ans[[paste0("pval", s)]]  <- sapply(rv, `[`, i=st+1L)
            ans[[paste0("p.med", s)]] <- sapply(rv, `[`, i=st+2L)
            ans[[paste0("p.min", s)]] <- sapply(rv, `[`, i=st+3L)
            ans[[paste0("p.max", s)]] <- sapply(rv, `[`, i=st+4L)
            st <- st + 4L
        }
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}



#######################################################################
# SAIGE ACAT-O tests
#

seqAssocGLMM_ACAT_O <- function(gdsfile, modobj, units, wbeta=AggrParamBeta,
    acatv.collapse.mac=10, skat.collapse.mac=10,
    skat.collapse.method=c("PA", "PA_int", "SumG"), burden.summac=3, dsnode="",
    res.savefn="", res.compress="LZMA", parallel=FALSE, verbose=TRUE,
    verbose.maf=FALSE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(inherits(units, "SeqUnitListClass"))
    .check_wbeta(wbeta)
    stopifnot(is.numeric(acatv.collapse.mac), length(acatv.collapse.mac)==1L,
        is.finite(acatv.collapse.mac))
    stopifnot(is.numeric(skat.collapse.mac), length(skat.collapse.mac)==1L,
        is.finite(skat.collapse.mac))
    skat.collapse.method <- match.arg(skat.collapse.method)
    stopifnot(is.numeric(burden.summac), length(burden.summac)==1L,
        is.finite(burden.summac))
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.maf), length(verbose.maf)==1L)

    if (verbose)
        .cat(.crayon_inverse("SAIGE ACAT-O analysis:"))

    # check model
    modobj <- .check_modobj(modobj, verbose)
    var.ratio <- .get_var_ratio(modobj)
    spa.pval <- NaN  # according to Cutoff=2 in SPAtest

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

    has_skat <- !is.null(modobj$Sigma_inv)
    if (has_skat)
    {
        # check packages
        pkg_cqf <- suppressPackageStartupMessages(requireNamespace(
            "CompQuadForm", quietly=TRUE))
        pkg_svy <- suppressPackageStartupMessages(requireNamespace(
            "survey", quietly=TRUE))
        if (!pkg_cqf || !pkg_svy)
            stop("The packages 'CompQuadForm' and 'survey' should be installed.")
    }
    if (verbose)
    {
        .cat("    MAC threshold for collapsing ultra rare variants for ACAT-V: <= ",
            sprintf("%.15g", acatv.collapse.mac))
        if (has_skat)
        {
            .cat("    MAC threshold for collapsing ultra rare variants for SKAT: <= ",
                sprintf("%.15g", skat.collapse.mac))
            cat("    ACAT-O p-values combine Burden, ACAT-V and SKAT\n")
        } else {
            cat("    ACAT-O p-values combine Burden and ACAT-V\n")
        }
    }

    # show beta weights
    if (!is.matrix(wbeta))
        wbeta <- matrix(wbeta, nrow=2L)
    wb_colnm <- sprintf("%g_%g", wbeta[1L,], wbeta[2L,])
    .show_wbeta(wbeta, verbose)

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # get allele frequencies
    if (verbose && isTRUE(verbose.maf)) .show_maf(gdsfile, parallel)

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, 0, 0, 1, spa.pval, var.ratio, 2L,
        modobj$Sigma_inv, modobj$chol_inv_X_Sigma,
        burden.summac, wbeta, sz_wmax, acatv.collapse.mac, skat.collapse.mac)
    i <- match(skat.collapse.method, c("PA", "PA_int", "SumG"))
    if (is.na(i)) stop("Internal error in 'skat.collapse.method'.")
    mobj$collapse.method <- i

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
        if (has_skat)
        {
            # initialize SKAT
            mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
            .Call(saige_skat_test_init, mobj$Sigma_inv_cg, mobj$t_XVX_inv_XV,
                mobj$Si_X, mobj$XVX_inv_XV_X_Si_X, mobj$collapse.method)
            # finalize
            on.exit(.Call(saige_skat_test_done), add=TRUE)
        } else {
            .Call(saige_skat_test_reset)
        }
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize
        seqParallel(parallel, NULL, split="none", .combine="none",
            FUN = function(mobj) {
                eval(.load_lib)
                mobj$Sigma_inv_cg <- .sp_to_dgCMatrix(mobj$Sigma_inv)
                .packageEnv$mobj <- mobj
                .Call(saige_score_test_init, mobj)
                if (!is.null(mobj$Sigma_inv))
                {
                    # initialize SKAT
                    .Call(saige_skat_test_init, mobj$Sigma_inv_cg, mobj$t_XVX_inv_XV,
                        mobj$Si_X, mobj$XVX_inv_XV_X_Si_X, mobj$collapse.method)
                } else {
                    .Call(saige_skat_test_reset)
                }
            }, mobj=mobj)
        # finalize
        on.exit({
            seqParallel(parallel, NULL, split="none", .combine="none",
                FUN = function() {
                    .packageEnv$mobj <- NULL
                    .Call(saige_skat_test_done)
                })
        }, add=TRUE)
    }

    # scan all variant units
    if (verbose)
        cat("Calculating p-values:\n")
    rv <- seqUnitApply(gdsfile, units, dsnode,
        FUN=function(x) .Call(saige_acato_test_pval, x), as.is="list",
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
        .aggr_ret_gds(outf, gdsfile, units, rv, Add)
        Add("pval", vapply(rv, `[`, 0, i=10L))
        nn <- ifelse(has_skat, 3L, 2L)
        # write p-values
        for (i in seq_len(ncol(wbeta)))
        {
            s <- wb_colnm[i]
            Add(paste0("pval.b", s), vapply(rv, `[`, 0, i=11L+(i-1L)*nn))
            Add(paste0("pval.v", s), vapply(rv, `[`, 0, i=12L+(i-1L)*nn))
            if (has_skat)
                Add(paste0("pval.s", s), vapply(rv, `[`, 0, i=13L+(i-1L)*nn))
        }
        if (verbose) cat(.crayon_inverse("Done.\n"))
        # output
        invisible()

    } else {
        # output
        ans <- .aggr_ret_obj(units, rv)
        ans$pval <- vapply(rv, `[`, 0, i=10L)
        nn <- ifelse(has_skat, 3L, 2L)
        for (i in seq_len(ncol(wbeta)))
        {
            s <- wb_colnm[i]
            ans[[paste0("pval.b", s)]] <- vapply(rv, `[`, 0, i=11L+(i-1L)*nn)
            ans[[paste0("pval.v", s)]] <- vapply(rv, `[`, 0, i=12L+(i-1L)*nn)
            if (has_skat)
                ans[[paste0("pval.s", s)]] <- vapply(rv, `[`, 0, i=13L+(i-1L)*nn)
        }
        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)
    }
}
