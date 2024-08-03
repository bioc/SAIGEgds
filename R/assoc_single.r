#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2019-2024    Xiuwen Zheng / AbbVie-ComputationalGenomics
# License: GPL-3
#


#######################################################################
# Internal R functions

.load_lib <- quote(
    suppressPackageStartupMessages(library("SAIGEgds", quietly=TRUE)))

.trait_list <- c("quantitative", "binary")


# Internal model initialization
.init_nullmod <- function(modobj, ii, maf, mac, missing, spa.pval, var.ratio,
    geno.ploidy, Sigma_inv, chol_inv_X_Sigma,
    summac=NaN, wbeta=double(), num_wbuf=0L, acatv_mac=10, skat_mac=10)
{
    # check
    if (!is.numeric(var.ratio) || anyNA(var.ratio))
        stop("Invalid variance ratio in the SAIGE model.")
    if (!is.null(Sigma_inv))
    {
        if (!is.matrix(Sigma_inv) && !inherits(Sigma_inv, "dsCMatrix"))
            stop("Sigma_inv should be a dense matrix or a sparse matrix 'dsCMatrix'.")
    }
    if (!is.null(chol_inv_X_Sigma) && !is.matrix(chol_inv_X_Sigma))
        stop("chol_inv_X_Sigma should be a dense matrix.")

    if (!identical(class(modobj$obj.noK), "list"))
    {
        if (identical(class(modobj$obj.noK), "SA_NULL"))
        {
            warning(
                "For more accurate model building, the null model should be built using SAIGEgds>=v1.9.1!",
                call.=FALSE, immediate.=TRUE)
        } else
            stop("Unknown model object.")
    }

    # initialize the internal model parameters
    y <- unname(modobj$obj.noK$y)[ii]
    mu <- unname(modobj$fitted.values)[ii]
    X1 <- modobj$obj.noK$X1[ii,, drop=FALSE]
    t_X1 <- t(X1)
    V <- modobj$obj.noK$V[ii]
    n <- length(ii)  # Num of samples
    K <- ncol(X1)    # Num of fixed-effect coefficients
    cateMAC <- attr(var.ratio, "cateMAC")
    mobj <- list(
        trait = match(modobj$trait.type, .trait_list),
        maf = maf, mac = mac, missing = missing, spa.pval = spa.pval,
        geno.ploidy = geno.ploidy,
        tau = modobj$tau,
        y = y, mu = mu, y_mu = y - mu,
        mu2 = mu * (1 - mu),
        t_XXVX_inv = t(modobj$obj.noK$XXVX_inv[ii,, drop=FALSE]),  # K x n_samp (K << n_samp, more efficient)
        XV = modobj$obj.noK$XV[, ii, drop=FALSE],  # K x n_samp
        t_XVX_inv_XV = t(modobj$obj.noK$XXVX_inv[ii,, drop=FALSE] * V),  # K x n_samp
        t_X = t_X1,  # K x n_samp
        var.ratio = var.ratio,
        vr_sqrt = sqrt(var.ratio),
        vr_inv_sqrt = 1/sqrt(var.ratio),
        vr_cateMAC = if (isFALSE(cateMAC)) double() else cateMAC,
        Sigma_inv = Sigma_inv, chol_inv_X_Sigma = chol_inv_X_Sigma,
        # buffer
        buf_dosage = double(n),
        buf_coeff = double(2L*K),
        buf_adj_g = double(n),
        buf_index = integer(n),
        buf_B = double(n),
        buf_g_tilde = double(n),
        buf_X1 = double(K),
        buf_spa = double(2L*n)
    )

    if (is.na(mobj$trait))
        stop("Invalid 'modobj$trait.type'.")

    # additional design matrix
    if (modobj$trait.type == "binary")
    {
        mobj$XVX <- t(X1) %*% (X1 * mobj$mu2)  # a matrix: K x K
        mobj$S_a <- colSums(X1 * mobj$y_mu)    # a vector of size K
    } else if (modobj$trait.type == "quantitative")
    {
        mobj$XVX <- t(X1) %*% X1               # a matrix: K x K
        mobj$S_a <- colSums(X1 * mobj$y_mu)    # a vector of size K
    } else
        stop("Invalid 'modobj$trait.type': ", modobj$trait.type, ".")

    if (!is.null(Sigma_inv))
    {
        mobj$Si_X <- as.matrix(t(Sigma_inv %*% X1))
        mobj$XVX_inv_XV_X_Si_X <-
            as.matrix(t_X1 %*% t(mobj$Si_X) %*% mobj$t_XVX_inv_XV)
    }

    # additional for aggregate tests
    mobj$summac <- summac
    mobj$acatv_mac <- acatv_mac
    mobj$skat_mac <- skat_mac
    mobj$buf_wbeta <- as.double(wbeta)
    mobj$num_unitsz <- num_wbuf
    mobj$buf_unitsz <- matrix(0.0, nrow=num_wbuf, ncol=7L)

    # output
    mobj
}

.dsnode <- function(gdsfile, nm, aggregate=FALSE)
{
    if (nm == "")
    {
        if (!aggregate)
        {
            n <- index.gdsn(gdsfile, "genotype/data", silent=TRUE)
            if (!is.null(n))
            {
                nm <- "$dosage_alt2"
            } else {
                nm <- getOption("seqarray.node_ds", "annotation/format/DS")
                n <- index.gdsn(gdsfile, nm, silent=TRUE)
                if (is.null(n))
                    stop("Dosages should be stored in genotype or annotation/format/DS.")
            }
        } else {
            nm <- "$dosage_sp2"
        }
    }
    nm
}

# check compression option
.check_compress <- function(res.compress)
{
    if (!is.character(res.compress))
        stop("'res.compress' should be character.")
    if (length(res.compress) != 1L)
        stop("'res.compress' should be character vector of length 1.")
    cp <- c("ZIP", "ZIP_RA", "LZMA", "LZMA_RA", "none")
    if (!(res.compress %in% cp))
    {
        stop("`res.compress` should be one of ",
            paste(paste(cp[-length(cp)], collapse=", "), "and", cp[length(cp)]),
            ".")
    }
}

# get the variance ratio estimate
.get_var_ratio <- function(modobj)
{
    # categorical MACs
    cateMAC <- modobj$use.cateMAC
    if (is.null(cateMAC)) cateMAC <- FALSE
    cateMAC <- .check_use_cateMAC(cateMAC)
    # variance ratio(s)
    if (is.numeric(modobj$var.ratio))
    {
        vr <- modobj$var.ratio[1L]
        if (!isFALSE(cateMAC))
            stop("'use.cateMAC' should be FALSE for a single variance ratio!")
    } else if (is.data.frame(modobj$var.ratio))
    {
        if (isFALSE(cateMAC))
        {
            vr <- mean(modobj$var.ratio$ratio, na.rm=TRUE)
        } else {
            if (!is.numeric(cateMAC))
                stop("'use.cateMAC' should be a numeric vector for multiple variance ratios.")
            x <- modobj$var.ratio$ratio
            d <- cut(modobj$var.ratio$mac, c(0, cateMAC, Inf),
                right=FALSE, dig.lab=15L)
            vr <- suppressWarnings(
                vapply(levels(d), function(s) mean(x[d==s]), 0))
            if (length(vr) != length(cateMAC)+1L)
            {
                stop("'use.cateMAC' should be a ", length(vr)-1L,
                    "-length vector!")
            }
            if (anyNA(vr))
            {
                print(vr)
                stop("variance ratios should not contain NaN.")
            }
        }
    } else
        stop("Invalid variance ratio in the model!")
    if (any(!is.finite(vr)))
        stop("The variance ratio should be a finite number.")
    attr(vr, "cateMAC") <- cateMAC
    vr
}

# save to R object file
.save_R_obj <- function(obj, res.compress, res.savefn, verbose)
{
    isfn <- !is.na(res.savefn) && res.savefn!=""
    if (isfn)
    {
        cm <- switch(res.compress, LZMA="xz", LZMA_RA="xz", ZIP="gzip",
            ZIP_RA="gzip", TRUE)
        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        if (grepl("\\.(rda|RData)$", res.savefn, ignore.case=TRUE))
        {
            .res <- obj
            save(.res, file=res.savefn, compress=cm)
        } else if (grepl("\\.rds$", res.savefn, ignore.case=TRUE))
        {
            saveRDS(obj, file=res.savefn, compress=cm)
        } else {
            stop("Unknown format of the output file, and it should be RData, RDS or gds.")
        }
        if (verbose) .cat(.crayon_inverse("Done."))
        invisible()
    } else {
        if (verbose) .cat(.crayon_inverse("Done."))
        obj
    }
}



#######################################################################
# SAIGE single variant analysis
#

seqAssocGLMM_SPA <- function(gdsfile, modobj, maf=NaN, mac=10, missing=0.05,
    spa=TRUE, dsnode="", geno.ploidy=2L, res.savefn="", res.compress="ZIP",
    parallel=FALSE, load.balancing=TRUE, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(mac), length(mac)==1L)
    stopifnot(is.numeric(missing), length(missing)==1L)
    stopifnot(is.logical(spa), length(spa)==1L)
    stopifnot(is.character(dsnode), length(dsnode)==1L, !is.na(dsnode))
    stopifnot(is.numeric(geno.ploidy) | is.na(geno.ploidy),
        length(geno.ploidy)==1L)
    if (is.numeric(geno.ploidy) && !is.na(geno.ploidy))
        stopifnot(geno.ploidy >= 0L)
    stopifnot(is.character(res.savefn), length(res.savefn)==1L)
    .check_compress(res.compress)
    stopifnot(is.logical(load.balancing), length(load.balancing)==1L,
        !is.na(load.balancing))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
        cat(.crayon_inverse("SAIGE association analysis:\n"))

    # check model
    modobj <- .check_modobj(modobj, verbose)
    var.ratio <- .get_var_ratio(modobj)
    spa.pval <- if (isTRUE(spa)) NaN else -1

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

    # show warnings immediately
    saveopt <- options(warn=1L)
    on.exit(options(warn=saveopt$warn), add=TRUE)

    # determine the GDS node for dosages
    dsnode <- .dsnode(gdsfile, dsnode)

    # check sample ID
    seqSetFilter(gdsfile, sample.id=modobj$sample.id, warn=FALSE, verbose=FALSE)
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs are not available in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")

    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    nVariant <- dm[3L]
    if (verbose)
    {
        .cat("    trait type: ", modobj$trait.type)
        .cat("    # of samples: ", .pretty(dm[2L]))
        .cat("    # of variants: ", .pretty(dm[3L]))
        if (!is.na(geno.ploidy) && (geno.ploidy>0L))
        {
            .cat("    MAF threshold: ", .pretty_gt_eq(maf))
            .cat("    MAC threshold: ", .pretty_gt_eq(mac))
            if (geno.ploidy != 2L)
                .cat("    ploidy: ", geno.ploidy)
        } else {
            cat("    not genotypic data\n")
        }
        .cat("    missing proportion threshold: ", .pretty_lt_eq(missing))
        v <- var.ratio
        if (isFALSE(attr(v, "cateMAC")))
        {
            .cat("    variance ratio for approximation: ", v)
        } else {
            cat("    variance ratios for approximation (MAC categories):\n")
            attr(v, "cateMAC") <- NULL
            print(v, width=1024L)
        }
        if (is.matrix(modobj$Sigma_inv))
            cat("    using a dense matrix of Sigma (covariance)\n")
        else if (inherits(modobj$Sigma_inv, "sparseMatrix"))
            cat("    using a sparse matrix of Sigma (covariance)\n")
        if (!isTRUE(spa) && (modobj$trait.type=="binary"))
            cat("    no adjustment for case-control imbalance\n")
    }

    if (dm[2L] <= 0) stop("No sample in the genotypic data set!")
    if (dm[3L] <= 0) stop("No variant in the genotypic data set!")

    # initialize the internal model parameters
    mobj <- .init_nullmod(modobj, ii, maf, mac, missing, spa.pval, var.ratio,
        geno.ploidy, modobj$Sigma_inv, modobj$chol_inv_X_Sigma)

    # update parallel object
    njobs <- SeqArray:::.NumParallel(parallel)
    parallel <- SeqArray:::.McoreParallel(parallel)
    is_fork <- SeqArray:::.IsForking(parallel)  # is forking or not?
    if (verbose)
        .cat("    # of processes: ", njobs)

    # initialize internally
    if (njobs<=1L || is_fork)
    {
        # forking, no need to distribute model parameters
        .Call(saige_score_test_init, mobj)
        initfun <- finalfun <- NULL
    } else {
        # pass the model parameters to each process
        if (verbose)
            cat("Distribute the model parameters to the", njobs, "processes\n")
        # initialize child processes internally
        initfun <- function(proc_id, mobj)
        {
            eval(.load_lib)
            .packageEnv$modobj <- mobj
            .Call(saige_score_test_init, mobj)
        }
        # clear when exit
        finalfun <- function(proc_id, param)
        {
            .packageEnv$modobj <- NULL
            remove(modobj, envir=.packageEnv)
        }
    }

    # reset memory GC to reduce memory usage (especially in forked processes)
    gc(verbose=FALSE, reset=TRUE, full=TRUE)

    # defined functions
    isfn <- !is.na(res.savefn) && res.savefn!=""
    isfn_gds <- isfn && grepl("\\.gds$", res.savefn, ignore.case=TRUE)
    if (!isfn_gds)
    {
        combine_fun <- "unlist"
        scan_fun <- function(f, dsnode, pverbose)
        {
            seqApply(f, dsnode, .cfunction("saige_score_test_pval"),
                as.is="list", parallel=FALSE,
                .progress=pverbose && SeqArray:::process_index==1L,
                .list_dup=FALSE, .useraw=NA)
        }
    } else {
        # output to a GDS file
        Append <- function(nm, val)
            append.gdsn(index.gdsn(outf, nm), val)
        combine_fun <- function(v)
        {
            Append("id", v$id)
            if (!is.na(geno.ploidy) && (geno.ploidy>0L))
            {
                Append("AF.alt", vapply(v$rv, `[`, 0, i=1L))
                Append("mac", vapply(v$rv, `[`, 0, i=2L))
            } else {
                Append("mean", vapply(v$rv, `[`, 0, i=1L))
                Append("nnzero", vapply(v$rv, `[`, 0, i=2L))
            }
            Append("num", vapply(v$rv, `[`, 0, i=3L))
            Append("beta", vapply(v$rv, `[`, 0, i=4L))
            Append("SE", vapply(v$rv, `[`, 0, i=5L))
            Append("pval", vapply(v$rv, `[`, 0, i=6L))
            if (modobj$trait.type == "binary")
            {
                Append("p.norm", vapply(v$rv, `[`, 0, i=7L))
                Append("converged", vapply(v$rv, `[`, 0, i=8L)==1)
            }
            NULL
        }
        scan_fun <- function(f, dsnode, pverbose)
        {
            v<- seqApply(f, dsnode, .cfunction("saige_score_test_pval"),
                as.is="list", parallel=FALSE,
                .progress=pverbose && SeqArray:::process_index==1L,
                .list_dup=FALSE, .useraw=NA)
            i <- seqGetData(f, "$variant_index")
            x <- !vapply(v, is.null, FALSE)
            list(id=i[x], rv=v[x])
        }

        if (verbose)
            .cat("Save to ", sQuote(res.savefn), " ...")
        # create a GDS file
        cm <- res.compress[1L]
        outf <- createfn.gds(res.savefn)
        on.exit({ if (!is.null(outf)) closefn.gds(outf) }, add=TRUE)
        put.attr.gdsn(outf$root, "FileFormat", "SAIGE_OUTPUT")
        put.attr.gdsn(outf$root, "Version",
            paste0("SAIGEgds_", packageVersion("SAIGEgds")))

        # add function
        Add <- function(varnm, val)
            add.gdsn(outf, varnm, val, compress=cm, closezip=TRUE, replace=TRUE)
        # add sample IDs
        Add("sample.id", seqGetData(gdsfile, "sample.id"))
        add.gdsn(outf, "id", integer())
        if (!is.na(geno.ploidy) && (geno.ploidy>0L))
        {
            add.gdsn(outf, "AF.alt", double())
            add.gdsn(outf, "mac", double())
        } else {
            add.gdsn(outf, "mean", double())
            add.gdsn(outf, "nnzero", integer())
        }
        add.gdsn(outf, "num", integer())
        add.gdsn(outf, "beta", double())
        add.gdsn(outf, "SE", double())
        add.gdsn(outf, "pval", double())
        if (modobj$trait.type == "binary")
        {
            add.gdsn(outf, "p.norm", double())
            add.gdsn(outf, "converged", logical())
        }
    }

    # scan variants
    rv <- seqParallel(parallel, gdsfile, split="by.variant",
        .initialize=initfun, .finalize=finalfun, .initparam=mobj,
        .combine=combine_fun, FUN=scan_fun,
        .balancing=load.balancing, .bl_size=50000L, .bl_progress=verbose,
        dsnode=dsnode,
        pverbose=verbose && ((njobs==1L) || !load.balancing))

    # output to a GDS file?
    if (!isfn_gds)
    {
        # if any maf/mac filter
        if (length(rv) != nVariant)
            stop("seqParallel() returns a vector of wrong length.")
        x <- vapply(rv, is.null, FALSE)
        if (any(x))
        {
            x <- !x
            seqSetFilter(gdsfile, variant.sel=x, action="intersect", verbose=FALSE)
            rv <- rv[x]
        }
        if (verbose)
        {
            cat("# of variants after filtering by ")
            if (!is.na(geno.ploidy) && (geno.ploidy>0L))
                cat("MAF, MAC and ")
            .cat("missing thresholds: ", .pretty(length(rv)))
        }

        # output
        ans <- data.frame(
            id  = seqGetData(gdsfile, "variant.id"),
            chr = seqGetData(gdsfile, "chromosome"),
            pos = seqGetData(gdsfile, "position"),
            stringsAsFactors = FALSE
        )
        # add RS IDs if possible
        if (!is.null(index.gdsn(gdsfile, "annotation/id", silent=TRUE)))
            ans$rs.id <- seqGetData(gdsfile, "annotation/id")
        ans$ref <- seqGetData(gdsfile, "$ref")
        ans$alt <- seqGetData(gdsfile, "$alt")
        if (!is.na(geno.ploidy) && (geno.ploidy>0L))
        {
            ans$AF.alt <- vapply(rv, `[`, 0, i=1L)
            ans$mac <- vapply(rv, `[`, 0, i=2L)
        } else {
            ans$mean <- vapply(rv, `[`, 0, i=1L)
            ans$nnzero <- vapply(rv, `[`, 0, i=2L)
        }
        ans$num  <- as.integer(vapply(rv, `[`, 0, i=3L))
        ans$beta <- vapply(rv, `[`, 0, i=4L)
        ans$SE   <- vapply(rv, `[`, 0, i=5L)
        ans$pval <- vapply(rv, `[`, 0, i=6L)
        if (modobj$trait.type == "binary")
        {
            ans$p.norm <- vapply(rv, `[`, 0, i=7L)
            ans$converged <- vapply(rv, `[`, 0, i=8L)==1
        }

        # save file?
        .save_R_obj(ans, res.compress, res.savefn, verbose)

    } else {
        # output to a GDS file
        # write data
        i <- read.gdsn(index.gdsn(outf, "id"))
        if (verbose)
        {
            cat("# of variants after filtering by ")
            if (!is.na(geno.ploidy) && (geno.ploidy>0L))
                cat("MAF, MAC and ")
            .cat("missing thresholds: ", .pretty(length(i)))
        }
        seqSetFilter(gdsfile, variant.sel=i, warn=FALSE, verbose=FALSE)
        nd_id <- .write_gds(outf, "id", gdsfile, "variant.id", cm)
        nd_chr <- .write_gds(outf, "chr", gdsfile, "chromosome", cm)
        moveto.gdsn(nd_chr, nd_id)
        nd_pos <- .write_gds(outf, "pos", gdsfile, "position", cm)
        moveto.gdsn(nd_pos, nd_chr)
        nd <- nd_pos
        # rs.id
        if (!is.null(index.gdsn(gdsfile, "annotation/id", silent=TRUE)))
        {
            nd <- .write_gds(outf, "rs.id", gdsfile, "annotation/id", cm)
            moveto.gdsn(nd, nd_pos)
        }
        # ref and alt alleles
        nd_ref <- Add("ref", seqGetData(gdsfile, "$ref"))
        moveto.gdsn(nd_ref, nd)
        nd_alt <- Add("alt", seqGetData(gdsfile, "$alt"))
        moveto.gdsn(nd_alt, nd_ref)
        # sub variables
        i <- order(i)
        for (nm in c("AF.alt", "mac", "mean", "nnzero", "num", "beta", "SE",
            "pval", "p.norm", "converged"))
        {
            nd <- index.gdsn(outf, nm, silent=TRUE)
            if (!is.null(nd)) Add(nm, read.gdsn(nd)[i])
        }

        # close the GDS file
        closefn.gds(outf); outf <- NULL
        cleanup.gds(res.savefn, verbose=FALSE)
        if (verbose) .cat(.crayon_inverse("Done."))
        invisible()
    }
}
