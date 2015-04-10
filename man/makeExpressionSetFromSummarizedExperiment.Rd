\name{makeExpressionSetFromSummarizedExperiment}

\alias{makeExpressionSetFromSummarizedExperiment}

\title{Make an ExpressionSet object from a SummarizedExperiment}

\description{
  \code{makeExpressionSetFromSummarizedExperiment} takes a \link[Biobase]{ExpressionSet}
  object as input and a Range mapping function that maps the features to
  ranges.  It then returns a \link{SummarizedExperiment} object.
}

\usage{
makeExpressionSetFromSummarizedExperiment(from,
                                          mapFun=naiveRangeMapper)
}

\arguments{
  \item{from}{
    An \link[Biobase]{ExpressionSet} object.
  }
  \item{mapFun}{
    A function which takes an \link[Biobase]{ExpressionSet} object and
    returns a \link{GRanges}, or \link{GRangesList} object which
    corresponds to the genomic ranges used in the ExpressionSet.  Defaults to
    \link{naiveRangeMapper}.
  }
}
\value{
  A \link{SummarizedExperiment} object that corresponds to the input.
}

\seealso{
  \link{naiveRangeMapper}, \link{probeRangeMapper}
}
