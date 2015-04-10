\name{makeSummarizedExperimentFromExpressionSet}

\alias{makeSummarizedExperimentFromExpressionSet}

\title{Make an SummarizedExperiment object from an ExpressionSet}

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
    corresponds to the genomic ranges used in the ExpressionSet.  The
    \link[base]{rownames} of the returned \link{GRanges} are used to match the
    \link[Biobase]{featureNames} of the \link[Biobase]{ExpressionSet}.
    \link{naiveRangeMapper}.
  }
}

\value{
  A \link{SummarizedExperiment} object that corresponds to the input.
}

\seealso{
  \link{naiveRangeMapper}, \link{probeRangeMapper}, \link{geneRangeMapper}
}
\examples{

  library(Biobase)
  data(sample.ExpressionSet, package = "Biobase")

  # using the naive coercion
  makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)

  # using probe range mapper
  makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet, probeRangeMapper)

  # using the gene range mapper
  makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet,
                                            geneRangeMapper("TxDb.Hsapiens.UCSC.hg19.knownGene"))
}
