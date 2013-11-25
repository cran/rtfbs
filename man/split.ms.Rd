\name{split.ms}
\alias{split.ms}
\title{Split sequences}
\usage{
  \method{split}{ms} (x, f, drop = FALSE, ...)
}
\arguments{
  \item{x}{Multiple sequences object}

  \item{f}{Numeric window size, or Features object used to
  determine where to split sequences}

  \item{drop}{Currently not used (for S3 compatibility)}

  \item{...}{Currently not used (for S3 compatibility)}
}
\value{
  MS object, containing the split sequences
}
\description{
  Split sequences in MS object at given locations, or based
  on a numeric window size. If splitting based on a numeric
  window size, (i.e. every 500 bases) specify the maximum
  number of bases per sequence in the \code{f} parameter.
  If splitting based on specific locations, use
  read.feat(FullPath) to read in a BED file specifying the
  locations, and use the resulting object in the \code{f}
  parameter.
}
\seealso{
  \code{\link{read.ms}}, or read.feat from package rphast
  for more about Features object
}

