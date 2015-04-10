\name{mappingFunctions}

\alias{mappingFunctions}
\alias{naiveRangeMapper}
\alias{probeRangeMapper}
\alias{geneRangeMapper}

\title{Mapping Functions used in \link{makeSummarizedExperimentFromExpressionSet}}

\description{
  These mapping functions are used in the \code{mapFun} argument to
  \link{makeSummarizedExperimentFromExpressionSet}.
}

\usage{
  naiveRangeMapper(from)

  probeRangeMapper(from)

  geneRangeMapper(txDbPackage, key = "ENTREZID")
}

\arguments{
  \item{from}{
    An \link[Biobase]{ExpressionSet} object.
  }
  \item{txDbPackage}{
    A character string with the Transcript Database to use for the mapping.
  }
  \item{key}{
    A character string with the Gene key to use for the mapping.
  }
}

\value{
  A \link{GRanges} object, with the \code{rownames} corresponding to the
  \link[Biobase]{featureNames} of the \link[Biobase]{ExpressionSet} object.
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
