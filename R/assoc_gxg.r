#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2026    Xiuwen Zheng / AbbVie-ComputationalGenomics
# License: GPL-3
#


#######################################################################
# Calculate the matrix for the genotype-by-genotype interaction terms
#

seqGetGenoByGeno <- function(gt1, gt2, use_colnm=TRUE)
{
    # check
    stopifnot(is.matrix(gt1) | inherits(gt1, "Matrix"),
        is.matrix(gt2) | inherits(gt2, "Matrix"))
    stopifnot(nrow(gt1) == nrow(gt2))
    stopifnot(ncol(gt1) > 0L, ncol(gt2) > 0L)
    stopifnot(is.logical(use_colnm), length(use_colnm)==1L)
    # colnames
    if (use_colnm)
    {
        cn1 <- colnames(gt1)
        if (is.null(cn1)) cn1 <- paste0("a", seq_len(ncol(gt1)))
        cn2 <- colnames(gt2)
        if (is.null(cn2))
        {
            cn2 <- paste0("b", seq_len(ncol(gt2)))
            colnames(gt2) <- cn2
        }
    }
    # swap if gt1 has more columns than gt2 for better efficiency
    if (ncol(gt1) > ncol(gt2))
    {
        tmp <- gt1; gt1 <- gt2; gt2 <- tmp
        if (use_colnm)
        {
            tmp <- cn1; cn1 <- cn2; cn2 <- tmp
        }
    }
    # for-loop
    lst <- lapply(seq_len(ncol(gt1)), function(i)
    {
        m <- drop0(gt1[, i] * gt2)
        # identify columns with variation
        n0 <- diff(m@p) == 0L  # non-empty columns
        if (any(n0)) m <- m[, !n0, drop=FALSE]
        # set column names
        if (use_colnm & ncol(m) > 0L)
            colnames(m) <- paste(cn1[i], colnames(m), sep="_x_")
        # return
        if (ncol(m) > 0L) m else NULL
    })
    # merge
    do.call(cbind, lst)
}



#######################################################################
# Calculate the matrix for the genotype-by-genotype interaction terms
#



