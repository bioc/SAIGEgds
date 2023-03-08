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


# Package-wide variable
.packageEnv <- new.env()


#######################################################################
# Internal functions
#

.set_option <- function(always_fastSPA=TRUE, use_avx=c("avx512f","avx2","no"),
    verbose=TRUE)
{
    use_avx <- match.arg(use_avx)
    use_avx <- match(use_avx, c("avx512f","avx2","no"))
	.Call(saige_set_option, always_fastSPA, use_avx, verbose)
	invisible()
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

.cat <- function(...) cat(..., "\n", sep="")

.load_pkg <- quote({
    library(Rcpp, quietly=TRUE)
    library(SAIGEgds, quietly=TRUE)
})

.pretty <- function(x)
{
    prettyNum(x, big.mark=",", scientific=FALSE)
}

.pretty_size <- function(x)
{
    if (x >= 1024^4)
        sprintf("%.1fT", x/1024^4)
    else if (x >= 1024^3)
        sprintf("%.1fG", x/1024^3)
    else if (x >= 1024^2)
        sprintf("%.1fM", x/1024^2)
    else if (x >= 1024)
        sprintf("%.1fK", x/1024)
    else
        sprintf("%g bytes", x)
}

.pretty_gt_eq <- function(x)
{
    if (is.na(x)) "no" else sprintf(">= %.15g", x)
}

.pretty_lt_eq <- function(x)
{
    if (is.na(x)) "no" else sprintf("<= %.15g", x)
}

SIMD <- function()
{
    .Call(saige_simd_version)
}

.rank_norm <- function(y, m=0, s=1)
{
    if (anyNA(y))
        stop("NA is not allowed in the inverse normal transform.")
    qnorm((rank(y) - 0.5)/length(y), mean=m, sd=s)
}

.crayon_inverse <- function(s)
{
    if (getOption("gds.crayon", TRUE) && requireNamespace("crayon", quietly=TRUE))
        s <- crayon::inverse(s)
    s
}

.crayon_underline <- function(s)
{
    if (getOption("gds.crayon", TRUE) && requireNamespace("crayon", quietly=TRUE))
        s <- crayon::underline(s)
    s
}

# Write to GDS file
.write_gds <- function(out.gds, out.nm, in.gds, in.nm, cm)
{
    i <- index.gdsn(in.gds, in.nm)
    n <- add.gdsn(out.gds, out.nm, storage=i, compress=cm)
    seqApply(in.gds, in.nm, `c`, as.is=n)
    readmode.gdsn(n)
    invisible()
}

# Check null model
.check_modobj <- function(modobj, verbose)
{
    if (is.character(modobj))
    {
        stopifnot(length(modobj)==1L)
        if (verbose)
            .cat("    load the null model from ", sQuote(modobj))
        if (grepl("\\.(rda|RData)$", modobj, ignore.case=TRUE))
        {
            modobj <- get(load(modobj))
        } else if (grepl("\\.rds$", modobj, ignore.case=TRUE))
        {
            modobj <- readRDS(modobj)
        } else
            stop("It should be an RData, rda or rds file.")
    }
    stopifnot(inherits(modobj, "ClassSAIGE_NullModel"))
    modobj
}

# Show the distribution
.show_outcome <- function(trait.type, y, phenovar=NULL)
{
    if (trait.type == "binary")
    {
        # binary outcome
        if (!is.null(phenovar))
            .cat("Binary outcome: ", phenovar)
        v <- table(y)
        n <- length(v) 
        v <- data.frame(v, as.numeric(prop.table(v)))
        v[, 1L] <- paste0("      ", v[, 1L])
        colnames(v) <- c(phenovar, "Number", "Proportion")
        print(v, row.names=FALSE)
        if (n != 2L)
            stop("The outcome variable has more than 2 categories!")
    } else if (trait.type == "quantitative")
    {
        # quantitative outcome
        if (!is.null(phenovar))
            .cat("Quantitative outcome: ", phenovar)
        v <- data.frame(mean=mean(y), sd=sd(y), min=min(y), max=max(y))
        rownames(v) <- "   "
        print(v)
    }
}


# S3 print method
print.ClassSAIGE_NullModel <- function(x, ...) str(x)


#######################################################################
# p-value from Cauchy-based ACAT combination method
#

pACAT <- function(p, w=NULL)
{
    .Call(saige_acat_p, p, w)
}

pACAT2 <- function(p, maf, wbeta=c(1,25))
{
    stopifnot(is.numeric(wbeta), length(wbeta)==2L)
    stopifnot(length(p) == length(maf))
    w <- dbeta(maf, wbeta[1L], wbeta[2L])
    .Call(saige_acat_p, p, w * w * maf * (1-maf))
}



#######################################################################
# Load the association p-values in a GDS file
#

seqSAIGE_LoadPval <- function(fn, varnm=NULL, index=NULL, verbose=TRUE)
{
    # check
    stopifnot(is.character(fn), length(fn)>0L, all(!is.na(fn)))
    stopifnot(is.null(varnm) || is.character(varnm))
    stopifnot(is.null(index) || is.numeric(index) || is.logical(index))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    fn <- unique(fn)
    if (length(fn) == 1L)
    {
        if (verbose)
            .cat("Loading ", sQuote(fn), " ...")
        if (grepl("\\.gds$", fn, ignore.case=TRUE))
        {
            f <- openfn.gds(fn)
            on.exit(closefn.gds(f))
            fm <- get.attr.gdsn(f$root)$FileFormat[1L]
            if (fm %in% c("SAIGE_OUTPUT", "SAIGE_OUTPUT_SET"))
            {
                if (is.null(varnm))
                    varnm <- ls.gdsn(f$root)
                varnm <- setdiff(varnm, "sample.id")
                rv <- list()
                for (nm in varnm)
                    rv[[nm]] <- readex.gdsn(index.gdsn(f, nm), index)
                rv <- as.data.frame(rv, stringsAsFactors=FALSE)
            } else {
                stop("FileFormat should be 'SAIGE_OUTPUT' or 'SAIGE_OUTPUT_BURDEN'.")
            }
        } else if (grepl("\\.(rda|RData)$", fn, ignore.case=TRUE))
        {
            rv <- get(load(fn))
            if (!is.null(varnm)) rv <- rv[, varnm]
            if (!is.null(index)) rv <- rv[index, ]
        } else if (grepl("\\.rds$", fn, ignore.case=TRUE))
        {
            rv <- readRDS(fn)
            if (!is.null(varnm)) rv <- rv[, varnm]
            if (!is.null(index)) rv <- rv[index, ]
        } else
            stop("Unknown format, should be RData, RDS or gds.")
    } else {
        if (!is.null(index))
            stop("'index' should be NULL for multiple input files.")
        rv <- lapply(fn, function(nm)
            seqSAIGE_LoadPval(nm, varnm, verbose=verbose))
        if (verbose) cat("Merging ...\n")
        rv <- do.call(rbind, rv)
        if (verbose)
            cat(sprintf("Done: %d rows merged\n", NROW(rv)))
    }
    rv
}



#######################################################################
# Heritability estimation
#

glmmHeritability <- function(modobj, adjust=TRUE)
{
    # check
    stopifnot(is.logical(adjust), length(adjust)==1L)
    modobj <- .check_modobj(modobj, FALSE)

    if (modobj$trait.type == "binary")
    {
        tau <- modobj$tau[2L]
        r <- 1
        if (isTRUE(adjust))
        {
            y <- unname(modobj$obj.noK$y)
            p <- sum(y==1) / length(y)
            # coefficients are based on Supplementary Table 7 (Zhou et al. 2018)
            r <- 2.970 + 0.372*log10(p)
        }
        h <- tau / (pi*pi/3 + tau) * r
    } else if (modobj$trait.type == "quantitative")
    {
        h <- modobj$tau[2L] / sum(modobj$tau)
    } else
        stop("Invalid 'modobj$trait.type'.")

    unname(h)
}
