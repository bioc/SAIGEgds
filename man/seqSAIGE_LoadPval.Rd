\name{seqSAIGE_LoadPval}
\alias{seqSAIGE_LoadPval}
\title{Load the association results}
\description{
    Load the association results from a RData or GDS file.
}

\usage{
seqSAIGE_LoadPval(fn, varnm=NULL)
}

\arguments{
    \item{fn}{a RData or GDS file name}
    \item{varnm}{NULL, or a character vector to include the column names}
}

\value{
    Return a \code{data.frame} including p-values.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqFitNullGLMM_SPA}}, \code{\link{seqAssocGLMM_SPA}}
}

\keyword{GDS}
\keyword{genetics}
\keyword{association}
