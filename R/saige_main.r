#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Copyright (C) 2019        Xiuwen Zheng
# License: GPL-3
# Email: zhengxwen@gmail.com
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


#######################################################################
# Internal model checking
#
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


#######################################################################
# Open a SNP GDS file
#

seqAssocGMMAT_SPA <- function(gdsfile, modobj, maf=NaN, mac=NaN,
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
	mobj <- list(
	    maf = maf, mac = mac,
        y = y[ii], mu = mu[ii],
        y_mu = y[ii] - mu[ii],  # y - mu
        mu2 = (mu * (1 - mu))[ii],
        t_XXVX_inv = t(modobj$obj.noK$XXVX_inv[ii, ]),  # K x n_samp (K << n_samp, more efficient)
        XV = modobj$obj.noK$XV[, ii],  # K x n_samp
        var.ratio = mean(modobj$var.ratio, na.rm=TRUE),
        buf1 = double(nrow(modobj$obj.noK$XV)),
        buf2 = double(length(y))
	)
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
        cat("# of variants after MAF/MAC filtering: ", .pretty(length(rv)),
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
        converged  = sapply(rv, `[`, i=8L)==1,
        stringsAsFactors = FALSE
    )
}
