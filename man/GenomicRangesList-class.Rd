\name{GenomicRangesList-class}
\docType{class}

% Classes

\alias{class:GenomicRangesList}
\alias{GenomicRangesList-class}
\alias{GenomicRangesList}

\alias{class:GenomicRanges_OR_GenomicRangesList}
\alias{GenomicRanges_OR_GenomicRangesList-class}
\alias{GenomicRanges_OR_GenomicRangesList}

\alias{class:SimpleGenomicRangesList}
\alias{SimpleGenomicRangesList-class}
\alias{SimpleGenomicRangesList}

\alias{class:CompressedGenomicRangesList}
\alias{CompressedGenomicRangesList-class}
\alias{CompressedGenomicRangesList}

% Methods

\alias{updateObject,GenomicRangesList-method}
\alias{length,GenomicRangesList-method}
\alias{names,GenomicRangesList-method}
\alias{names<-,GenomicRangesList-method}
\alias{extractROWS,GenomicRangesList,ANY-method}
\alias{getListElement,GenomicRangesList-method}
\alias{unlist,GenomicRangesList-method}

\alias{show,GenomicRangesList-method}
\alias{seqnames,GenomicRangesList-method}
\alias{strand,GenomicRangesList-method}
\alias{elementMetadata,GenomicRangesList-method}
\alias{seqinfo,GenomicRangesList-method}
\alias{seqinfo,CompressedGenomicRangesList-method}

\alias{seqnames<-,CompressedGenomicRangesList-method}
\alias{strand<-,CompressedGenomicRangesList,ANY-method}
\alias{strand<-,CompressedGenomicRangesList,character-method}
\alias{elementMetadata<-,CompressedGenomicRangesList-method}
\alias{seqinfo<-,CompressedGenomicRangesList-method}
\alias{start<-,CompressedGenomicRangesList-method}
\alias{end<-,CompressedGenomicRangesList-method}
\alias{width<-,CompressedGenomicRangesList-method}

\alias{score,GenomicRangesList-method}
\alias{score<-,GenomicRangesList-method}

\title{GenomicRangesList objects}

\description{
  The GenomicRangesList \emph{virtual} class is a general container for
  storing a list of \link{GenomicRanges} objects.

  Most users are probably more interested in the \link{GRangesList}
  container, a GenomicRangesList derivative for storing a list of
  \link{GRanges} objects.
}

\details{
  The place of GenomicRangesList in the \emph{Vector class hierarchy}:
  \preformatted{
                               Vector
                                 ^
                                 |
                                List
                                 ^
                                 |
                             RangesList
                              ^     ^
                             /       \
                            /         \
                           /           \
                          /             \
                         /               \
                        /                 \
         IntegerRangesList               GenomicRangesList
                ^                                 ^
                |                                 |
           IRangesList                       GRangesList
             ^     ^                           ^     ^
            /       \                         /       \
           /         \                       /         \
          /           \                     /           \
   SimpleIRangesList   \             SimpleGRangesList   \
            CompressedIRangesList             CompressedGRangesList
  }
  Note that the \emph{Vector class hierarchy} has many more classes.
  In particular \link[S4Vectors]{Vector}, \link[S4Vectors]{List},
  \link[IRanges]{IRangesList}, and \link[IRanges]{IntegerRangesList}
  have other subclasses not shown here.
}

\author{H. Pagès and M. Lawrence}

\seealso{
  \itemize{
    \item \link{GRangesList} objects.

    \item \link{GenomicRanges} and \link{GRanges} objects.
  }
}

\keyword{classes}
\keyword{methods}
