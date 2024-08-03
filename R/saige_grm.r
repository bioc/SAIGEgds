#######################################################################
#
# Package name: SAIGEgds
#
# Description:
#     Scalable and accurate implementation of generalized mixed models
# using GDS files
#
# Copyright (C) 2022-2024    Xiuwen Zheng / AbbVie-ComputationalGenomics
# License: GPL-3
#


# Select LD-pruned variants for genetic relationship matrix
seqFitLDpruning <- function(gdsfile, sample.id=NULL, variant.id=NULL,
    ld.threshold=0.1, maf=0.01, missing.rate=0.005, autosome.only=TRUE,
    use.cateMAC=TRUE, num.marker=100L, num.total=100000L, save.gdsfn=NULL,
    seed=200L, parallel=FALSE, parallel.multi=NULL, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.numeric(ld.threshold), length(ld.threshold)==1L)
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.logical(autosome.only), length(autosome.only)==1L)
    use.cateMAC <- .check_use_cateMAC(use.cateMAC)
    stopifnot(is.numeric(num.marker), length(num.marker)==1L, num.marker>0L)
    stopifnot(is.numeric(num.total), length(num.total)==1L, num.total>0L)
    if (is.na(num.total)) num.total <- Inf
    stopifnot(is.null(save.gdsfn) | is.character(save.gdsfn))
    if (is.character(save.gdsfn))
        stopifnot(!anyNA(save.gdsfn), length(save.gdsfn)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check packages
    if (!suppressPackageStartupMessages(requireNamespace("SNPRelate")))
        stop("The package 'SNPRelate' should be installed.")

    stopifnot(is.numeric(seed), length(seed)==1L, is.finite(seed))
    if (!is.null(seed)) set.seed(seed)
    if (is.null(parallel.multi)) parallel.multi <- parallel

    is_multi_fn <- is.character(gdsfile) && length(gdsfile)>1L
    if (verbose && !is_multi_fn)
        cat(.crayon_inverse("Perform linkage disequilibrium (LD) pruning:\n"))

    # GDS file(s)
    if (is_multi_fn)
    {
        if (anyDuplicated(gdsfile))
            stop("'gdsfile' has duplicated file name(s).")
        if (anyNA(gdsfile))
            stop("'gdsfile' has NA.")
        num.marker.each <- ceiling(num.marker / length(gdsfile))
        # for-each file name
        ans <- lapply(seq_along(gdsfile), function(i)
        {
            if (verbose && i>1L)
                .cat(paste(rep("=", 72L), collapse=""))
            # reset memory GC to reduce memory usage
            gc(verbose=FALSE, reset=TRUE, full=TRUE)
            # run LD pruning
            v <- seqFitLDpruning(gdsfile[i], sample.id=sample.id,
                ld.threshold=ld.threshold, maf=maf, missing.rate=missing.rate,
                autosome.only=autosome.only, use.cateMAC=use.cateMAC,
                num.marker=num.marker.each, num.total=Inf,
                parallel=parallel, verbose=verbose)
            if (is.null(names(v))) names(v) <- rep("", length(v))
            v
        })
        # check num.total
        ns <- vapply(ans, function(v) sum(names(v)==""), 0L)
        if (sum(ns) > num.total)
        {
            d <- data.frame(
                g = rep(1:length(ns), times=ns),
                i = unlist(lapply(ans, function(v) which(names(v)==""))))
            if (verbose)
            {
                .cat(paste(rep("=", 72L), collapse=""))
                .cat("Randomly select ", num.total, " from ", nrow(d),
                    " LD-pruned variants.")
            }
            d <- d[sort(sample.int(nrow(d), num.total)), ]
            ans <- lapply(seq_along(ns), function(j)
            {
                ii <- ans[[j]]
                jj <- d$i[d$g==j]
                c(ii[names(ii)!=""], ii[jj])
            })
        }
        names(ans) <- gdsfile
        # merging process
        if (!is.null(save.gdsfn))
        {
            if (verbose)
                .cat(.crayon_inverse(paste(c("\n", rep("=",72L)), collapse="")))
            # reset memory GC to reduce memory usage
            gc(verbose=FALSE, reset=TRUE, full=TRUE)
            save.gdsfn.part <- paste0(save.gdsfn, "_part", seq_along(gdsfile))
            on.exit(unlink(save.gdsfn.part, force=TRUE))

            # in parallel if possible
            seqParApply(parallel.multi, seq_along(gdsfile),
                FUN = function(i, snplst, gdsfn, sampid, verbose)
            {
                outfn <- paste0(gdsfn, "_part", i)
                if (verbose)
                {
                    if (sink.number() > 0) sink()  # TODO
                    cat(paste0("Export to ", sQuote(basename(outfn)),
                        " [start at ", date(), "]\n"))
                }
                f <- seqOpen(names(snplst)[i])
                on.exit(seqClose(f))
                seqSetFilter(f, sample.id=sampid, variant.id=snplst[[i]],
                    warn=FALSE, verbose=FALSE)
                seqExport(f, outfn, info.var=character(), fmt.var=character(),
                    samp.var=character(), verbose=FALSE)
                NULL  # no output
            }, snplst=ans, gdsfn=save.gdsfn, sampid=sample.id, verbose=verbose)

            # merge files
            if (verbose)
                .cat(.crayon_inverse(paste(rep("=",72L), collapse="")))
            seqMerge(save.gdsfn.part, save.gdsfn, verbose=verbose)
        }
        # output
        return(invisible(ans))

    } else if (is.character(gdsfile))
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

    nthread <- 1L
    if (is.numeric(parallel)) nthread <- parallel

    # run LD pruning using SNPRelate
    snpset <- SNPRelate::snpgdsLDpruning(gdsfile, sample.id, variant.id,
        autosome.only=autosome.only, maf=maf, missing.rate=missing.rate,
        ld.threshold=ld.threshold, num.thread=nthread, verbose=verbose)
    nchrom <- length(snpset)
    snpset <- unlist(snpset, use.names=FALSE)
    if (length(snpset) > num.total)
    {
        if (verbose)
        {
            .cat("Randomly select ", num.total, " from ",
                length(snpset), " variants.")
        }
        snpset <- snpset[sort(sample.int(length(snpset), num.total))]
    }

    # need variants for MAC categories
    if (!isFALSE(use.cateMAC))
    {
        if (verbose)
            cat("Select ultra rare variants ...\n")
        .show_use_cateMAC(use.cateMAC, verbose)
        seqSetFilter(gdsfile, sample.id=sample.id, verbose=FALSE)
        seqSetFilter(gdsfile, variant.id=variant.id, verbose=FALSE)
        ii <- seqGetData(gdsfile, "$variant_index")
        v <- seqGetAF_AC_Missing(gdsfile, minor=TRUE, parallel=parallel,
            verbose=verbose)
        mac <- v$ac
        # missing rates
        if (is.finite(missing.rate))
        {
            x <- v$miss <= missing.rate
            x[is.na(x)] <- FALSE
            ii <- ii[x]; mac <- mac[x]
            if (verbose)
            {
                .cat("    excluding ", .pretty(sum(!x)),
                    " variants according missing rates")
            }
        }
        remove(v)  # remove unused v
        # MAC categories
        var_ii <- NULL; nm <- NULL
        low <- 0
        for (u in use.cateMAC)
        {
            s <- paste0("[", low, ",", u, ")")
            x <- mac < u
            x[is.na(x)] <- FALSE
            if (any(x))
            {
                iix <- ii[x]
                n <- min(num.marker, length(iix))
                var_ii <- c(var_ii, sort(sample(iix, n)))
                nm <- c(nm, rep(s, n))
                ii <- ii[!x]; mac <- mac[!x]
            }
            if (verbose)
            {
                .cat("    ", s, "\trandom ", num.marker, " from ", sum(x),
                    " variants")
            }
            low <- u
        }
        # finally
        seqSetFilter(gdsfile, variant.sel=var_ii, warn=FALSE, verbose=FALSE)
        snpset <- unique(c(seqGetData(gdsfile, "variant.id"), snpset))
        names(snpset) <- c(nm, rep("", length(snpset)-length(nm)))
    }

    # output to a GDS file?
    if (is.character(save.gdsfn))
    {
        seqSetFilter(gdsfile, variant.id=snpset, warn=FALSE, verbose=FALSE)
        seqExport(gdsfile, save.gdsfn, info.var=character(),
            fmt.var=character(), samp.var=character(), verbose=verbose)
    } else {
        if (verbose && nchrom>1L)
            .cat("Done (", length(snpset), " variants in total).")
    }

    invisible(snpset)
}



#######################################################################
# Calculation of GRM

.simd_calc_cpuinfo <- function() .Call(saige_simd_sp_grm)

.fit_calc_sp_grm <- function(gdsfile, nsnp.sub.random, maf, missing.rate,
    rel.cutoff, num.thread, return.ID, verbose, verbose_progress)
{
    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    nsamp <- dm[2L]
    nvar <- nvar_tot <- dm[3L]
    if ((nsnp.sub.random > 0L) && (nvar > nsnp.sub.random) && (nvar > 1L))
    {
        # random indices
        ii <- sort(sample.int(nvar, nsnp.sub.random))
        nvar <- nsnp.sub.random
    }
    if (verbose)
    {
        cat("Genotype matrix:\n")
        .cat("    # of samples: ", .pretty(nsamp))
        if (nvar == nvar_tot)
        {
            .cat("    # of variants: ", .pretty(nvar))
        } else {
            .cat("    # of variants: ", .pretty(nvar),
                " randomly selected from the full set (",
                .pretty(nvar_tot), ")")
        }
        .cat("    MAF threshold: >= ", maf)
        .cat("    missing genotype threshold: <= ", missing.rate)
        .cat("    # of threads: ", num.thread)
        .cat("    relatedness threshold: ", rel.cutoff)
        .cat("    using CPU capability: ", .simd_calc_cpuinfo())
    }
    if (nsamp <= 0L || nvar <= 0L)
        stop("No selected sample or variant.")

    # fill 2-bit genotypes
    if (verbose_progress)
    {
        if (nvar < nvar_tot)
        {
            cat(">>>>  First step (sparse GRM)  <<<<\n")
            cat("Loading the genotypes (the full SNP set):\n")
        } else
            cat("Loading the genotypes:\n")
    }
    nr <- ceiling(nvar_tot / 4L)  # in bytes
    ext_nb <- ceiling(nr/4L)*4L - nr  # 32-bit aligned
    g_pack <- seqGet2bGeno(gdsfile, samp_by_var=FALSE, ext_nbyte=ext_nb,
        verbose=verbose)
    g_pack2 <- g_pack
    if (verbose_progress)
    {
        cat("    ", nrow(g_pack), " x ", ncol(g_pack), ": ", sep="")
        print(object.size(g_pack))
    }
    # randomly selected SNPs
    if (nvar < nvar_tot)
    {
        # rearrange, the randomly selected SNPs are stored before non-selection
        j <- rep(1L, nrow(g_pack)*4L)
        j[ii] <- 0L
        j <- order(j) - 1L  # ZERO-based
        buf <- raw(nrow(g_pack))
        .Call(saige_grm_sp_reraw, g_pack, j, buf)
        remove(buf)
        # get a new RAW matrix
        nr <- ceiling(nvar/4L)   # in bytes
        nr <- ceiling(nr/4L)*4L  # 32-bit aligned
        if (nr < nrow(g_pack))
            g_pack2 <- g_pack[seq_len(nr), ]
    }

    # initialize
    g_lookup <- matrix(NaN, nrow=8L, ncol=nrow(g_pack)*4L)
    bl_size <- 256L
    n <- ceiling(nsamp / bl_size)
    n <- n*(n+1L)/2L    # total number of blocks

    # calculation
    if (verbose_progress)
    {
        if (nvar < nvar_tot)
            .cat("Calculating GRM from the reduced SNP set (m=", nvar, "):")
        else
            cat("Calculating GRM:\n")
    }
    prog_func <- SeqArray:::.seqProgForward
    prog <- if (verbose_progress) SeqArray:::.seqProgress(n, 1L) else NULL
    v <- .Call(saige_grm_sp_calc, nvar, g_pack2, g_lookup, rel.cutoff,
        bl_size, prog, prog_func)
    remove(prog)

    # using the full set for non-zero entries
    n <- length(v$i)
    if (verbose_progress)
    {
        .cat("# of non-zero entries in GRM: ", .pretty(n))
        cat(">>>>  Second step (sparse GRM)  <<<<\n")
        .cat("Calculating the non-zero entries using the full SNP set (m=",
            nvar_tot, "):")
    }
    bl_size <- 1024L
    n <- ceiling(n/bl_size)
    prog <- if (verbose_progress) SeqArray:::.seqProgress(n, 1L) else NULL
    # v$x will be updated
    v$x <- .Call(saige_grm_sp_calc_ijx, v$i, v$j, nvar_tot,
        g_pack, g_lookup, bl_size, prog, prog_func)
    remove(prog)

    # output
    b <- any(v$x < rel.cutoff)
    if (is.na(b)) b <- FALSE
    if (b)
    {
        z <- v$x >= rel.cutoff
        m <- sparseMatrix(i=v$i[z], j=v$j[z], x=v$x[z], dims=c(nsamp, nsamp),
            symmetric=TRUE, index1=FALSE)
    } else {
        m <- sparseMatrix(i=v$i, j=v$j, x=v$x, dims=c(nsamp, nsamp),
            symmetric=TRUE, index1=FALSE)
    }

    samp_id <- seqGetData(gdsfile, "sample.id")
    colnames(m) <- rownames(m) <- samp_id
    if (isTRUE(return.ID))
    {
        m <- list(sample.id = samp_id,
            variant.id = seqGetData(gdsfile, "variant.id"))
        if (nvar < nvar_tot)
            m$variant.sub.id <- m$variant.id[ii]
        m$grm <- m
    }
    m
}

# create a sparse genetic relationship matrix
seqFitSparseGRM <- function(gdsfile, sample.id=NULL, variant.id=NULL,
    nsnp.sub.random=2000L, rel.cutoff=0.125, maf=0.01, missing.rate=0.005,
    num.thread=1L, return.ID=FALSE, seed=200L, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.numeric(nsnp.sub.random), length(nsnp.sub.random)==1L,
        nsnp.sub.random>=0L)
    stopifnot(is.numeric(rel.cutoff), length(rel.cutoff)==1L)
    if (is.na(rel.cutoff)) rel.cutoff <- -Inf
    stopifnot(is.numeric(maf), length(maf) %in% c(1L,2L))
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.logical(return.ID), length(return.ID)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
        .cat(.crayon_inverse("Genetic Relationship Matrix (GRM):"))
    if (nsnp.sub.random > 0L)
    {
        nsnp.sub.random <- as.integer(floor(nsnp.sub.random/4) * 4L)
        if (nsnp.sub.random < 1L) nsnp.sub.random <- 4L
    }
    if (!is.null(seed)) set.seed(seed)

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

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    .Call(saige_set_numthread, num.thread)

    # set sample & variant filters
    if (!is.null(sample.id))
        seqSetFilter(gdsfile, sample.id=sample.id, verbose=FALSE)
    if (!is.null(variant.id))
        seqSetFilter(gdsfile, variant.id=variant.id, verbose=FALSE)
    if (verbose)
        cat("Filtering variants:\n")
    seqSetFilterCond(gdsfile, maf=maf, missing.rate=missing.rate,
         parallel=num.thread, .progress=verbose, verbose=FALSE)

    # calculating ...
    m <- .fit_calc_sp_grm(gdsfile, nsnp.sub.random, maf, missing.rate,
        rel.cutoff, num.thread, return.ID, verbose, verbose)
    # output
    if (verbose) .cat(.crayon_inverse("Done."))
    m
}


# create a dense genetic relationship matrix
seqFitDenseGRM <- function(gdsfile, sample.id=NULL, variant.id=NULL,
    maf=0.01, missing.rate=0.005, num.thread=1L, use.double=TRUE,
    return.ID=FALSE, verbose=TRUE)
{
    stopifnot(inherits(gdsfile, "SeqVarGDSClass") | is.character(gdsfile))
    stopifnot(is.numeric(maf), length(maf) %in% c(1L,2L))
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.logical(return.ID), length(return.ID)==1L)
    stopifnot(is.logical(use.double), length(use.double)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (verbose)
        .cat(.crayon_inverse("Genetic Relationship Matrix (GRM):"))

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

    # set the number of internal threads
    if (is.na(num.thread) || num.thread < 1L)
        num.thread <- 1L
    .Call(saige_set_numthread, num.thread)

    # set sample & variant filters
    if (!is.null(sample.id))
        seqSetFilter(gdsfile, sample.id=sample.id, verbose=FALSE)
    if (!is.null(variant.id))
        seqSetFilter(gdsfile, variant.id=variant.id, verbose=FALSE)
    if (verbose)
        cat("Filtering variants:\n")
    seqSetFilterCond(gdsfile, maf=maf, missing.rate=missing.rate,
         parallel=num.thread, .progress=verbose, verbose=FALSE)
    dm <- seqSummary(gdsfile, "genotype", verbose=FALSE)$seldim
    nsamp <- dm[2L]
    nvar <- dm[3L]
    if (verbose)
    {
        cat("Genotype matrix:\n")
        .cat("    # of samples: ", .pretty(nsamp))
        .cat("    # of variants: ", .pretty(nvar))
        .cat("    MAF threshold: >= ", maf)
        .cat("    missing genotype threshold: <= ", missing.rate)
        .cat("    # of threads: ", num.thread)
        .cat("    using CPU capability: ", .simd_calc_cpuinfo())
    }
    if (nsamp <= 0L || nvar <= 0L)
        stop("No selected sample or variant.")

    # fill 2-bit genotypes
    if (verbose)
        cat("Loading the genotypes:\n")
    nr <- ceiling(nvar / 4L)  # in bytes
    ext_nb <- ceiling(nr/4L)*4L - nr  # 32-bit aligned
    g_pack <- seqGet2bGeno(gdsfile, samp_by_var=FALSE, ext_nbyte=ext_nb,
        verbose=verbose)
    if (verbose)
    {
        cat("    ", nrow(g_pack), " x ", ncol(g_pack), ": ", sep="")
        print(object.size(g_pack))
    }

    # initialize
    g_lookup <- matrix(NaN, nrow=8L, ncol=nrow(g_pack)*4L)
    bl_size <- 256L
    n <- ceiling(nsamp / bl_size)
    n <- n*(n+1L)/2L    # total number of blocks

    # calculation
    if (verbose)
        cat("Calculating dense GRM:\n")
    prog_func <- SeqArray:::.seqProgForward
    prog <- if (verbose) SeqArray:::.seqProgress(n, 1L) else NULL
    m <- .Call(saige_grm_ds_calc, nvar, g_pack, g_lookup, use.double, bl_size,
        prog, prog_func)
    remove(prog)

    # output
    samp_id <- seqGetData(gdsfile, "sample.id")
    colnames(m) <- rownames(m) <- samp_id
    if (isTRUE(return.ID))
    {
        m <- list(sample.id = samp_id,
            variant.id = seqGetData(gdsfile, "variant.id"),
            grm = m)
    }
    if (verbose) .cat(.crayon_inverse("Done."))
    m
}
