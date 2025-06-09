\name{SurfaceColloc}
\docType{data}
\alias{SurfaceColloc}

\title{

  A small data set of surface collocations from the English Wikipedia

}

\description{

  This data set demonstrates how co-occurrence and marginal frequencies can be provided for collocation analysis with \code{\link{am.score}}.
  It contains surface co-occurrence counts for 7 English nouns as nodes and 7 selected collocates. The counts are based on a collocational span of two tokens to the left and right of the node (L2/R2) in the WP500 corpus.
  Marginal frequencies for the nodes are overall corpus frequencies of the nouns, so expected co-occurrence frequency needs to be adjusted with the total span size of 4 tokens.
  
}

\usage{

SurfaceColloc

}

\format{

  A list with the following components:
  
  \describe{
    \item{\code{cooc}:}{
      A data frame with 34 rows and the following columns:
      \itemize{
        \item{\code{w1}:} node word (noun)
        \item{\code{w2}:} collocate
        \item{\code{f}:} co-occurrence frequency within L2/R2 span
      }
    }
    \item{\code{f1}:}{
      Labelled integer vector of length 7 specifying the marginal frequencies of the node nouns.
    }
    \item{\code{f2}:}{
      Labelled integer vector of length 7 specifying the marginal frequencies of the collocates.
    }
    \item{\code{N}:}{
      Sample size, i.e. the total number of tokens in the WP500 corpus.
    }
  }

}

%%\references{
%%}

\examples{
head(SurfaceColloc$cooc, 10)
SurfaceColloc$f1
SurfaceColloc$f2
SurfaceColloc$N
}

\seealso{
  \code{\link{am.score}}
}

\author{Stephanie Evert (\url{https://purl.org/stephanie.evert})}
\keyword{datasets}
