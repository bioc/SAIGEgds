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

seqAssocGMMAT_SPA <- function(gdsfile, modobj, verbose=TRUE)
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
    seqSetFilter(gdsfile, sample.id=modobj$sample.id)
    sid <- seqGetData(gdsfile, "sample.id")
    if (length(sid) != length(modobj$sample.id))
        stop("")
    ii <- match()

    # initialize
    mu <- as.vector(modobj$fitted.values)
	mobj <- list(
        y     = modobj$obj.noK$y,
        mu.a  = mu,
        mu2.a = mu * (1 - mu),
        XXVX_inv = modobj$obj.noK$XXVX_inv,
        XV = modobj$obj.noK$XV,
        var.ratio = mean(modobj$var.ratio, na.rm=TRUE)
	)
    .Call(saige_score_test_init, mobj)

    # 
    if (modobj$trait.type == "binary")
    {
        rv <- seqApply(gdsfile, "annotation/format/DS",
            .cfunction("saige_score_test_bin"),
            as.is="list", .progress=verbose, .list_dup=FALSE)
    }

	return(rv)


    invisible()
}
