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
# Internal C functions
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

.cfunction2 <- function(name)
{
    fn <- function(x, y) NULL
    f <- quote(.Call(SEQ_ExternalName2, x, y))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.cfunction3 <- function(name)
{
    fn <- function(x, y, z) NULL
    f <- quote(.Call(SEQ_ExternalName3, x, y, z))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.cfunction4 <- function(name)
{
    fn <- function(w, x, y, z) NULL
    f <- quote(.Call(SEQ_ExternalName4, w, x, y, z))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}

.cfunction5 <- function(name)
{
    fn <- function(v, w, x, y, z) NULL
    f <- quote(.Call(SEQ_ExternalName5, v, w, x, y, z))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SAIGEgds")$address
    body(fn) <- f
    fn
}



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

seqAssocGMMAT_SPA <- function(gdsfile, modobj, parallel=FALSE, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check model
    if (is.character(modobj))
    {
        stopifnot(length(modobj)==1L)
        modobj <- get(load(modobj))
    }
    .check_saige_model(modobj)

    # check sample ID
    seqSetFilter(gdsfile, sample.id=modobj$sample.id, verbose=FALSE)
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("Some of sample IDs do not exist in the GDS file.")
    ii <- match(sid, modobj$sample.id)
    if (any(is.na(ii)))
        stop("Sample IDs do not match.")

    # initialize
    mu <- unname(modobj$fitted.values)
	mobj <- list(
        y     = unname(modobj$obj.noK$y[ii]),
        mu.a  = mu[ii],
        mu2.a = (mu * (1 - mu))[ii],
        XXVX_inv = modobj$obj.noK$XXVX_inv[ii, ],
        XV = modobj$obj.noK$XV[, ii],
        var.ratio = mean(modobj$var.ratio, na.rm=TRUE)
	)
    if (!is.finite(mobj$var.ratio))
        stop("Invalid variance ratio in the SAIGE model.")
    .Call(saige_score_test_init, mobj)

    # scan
    if (modobj$trait.type == "binary")
    {
        rv <- seqApply(gdsfile, "annotation/format/DS",
            .cfunction("saige_score_test_bin"), as.is="list",
            parallel=parallel, .progress=verbose, .list_dup=FALSE)
    }

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
        stringsAsFactors = FALSE
    )
}
