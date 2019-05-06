\name{seqSAIGE_LoadPval}
\alias{seqSAIGE_LoadPval}
\title{Load the association results}
\description{
    Load the association results from a RData or GDS file.
}

\usage{
seqSAIGE_LoadPval(fn, varnm=NULL, verbose=TRUE)
}

\arguments{
    \item{fn}{RData or GDS file names, merging datasets if multiple files}
    \item{varnm}{NULL, or a character vector to include the column names}
    \item{verbose}{if \code{TRUE}, show information}
}

\value{
    Return a \code{data.frame} including p-values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqAssocGLMM_SPA}}
}

\examples{
fn <- system.file("unitTests/saige_pval.rda", package="SAIGEgds")

pval <- seqSAIGE_LoadPval(fn)
head(pval)
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
