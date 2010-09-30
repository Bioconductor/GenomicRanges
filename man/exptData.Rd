\name{exptData}

\alias{exptData}
\alias{exptData<-}

\alias{rowData}
\alias{rowData<-}

\alias{colData}
\alias{colData<-}

\alias{assays}
\alias{assays<-}

\alias{assay}
\alias{assay<-}

\title{Accessing SummarizedExperiment instances}

\description{

  The generics documented on this page get and set data in the
  \code{\linkS4class{SummarizedExperiment}-class}. See that page for
  additional detail.

}

\usage{

exptData(x, ...)
exptData(x, ...) <- value

rowData(x, ...)
rowData(x, ...) <- value

colData(x, ...)
colData(x, ...) <- value

assays(x, ...)
assays(x, ...) <- value

assay(x, i, ...)
assay(x, i, ...) <- value

}

\arguments{
  \item{x}{The instance from which data is to be extracted.}

  \item{value}{The value to be assigned to the instance.}

  \item{i}{A subscript indicating which assay element to access or
    replace.}

  \item{...}{Additional arguments, currently ignored.}

}

\seealso{
  \code{\linkS4class{SummarizedExperiment}-class}
}

\examples{
  if (interactive())
    ?`SummarizedExperiment-class`
}

\keyword{methods}
