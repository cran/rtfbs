\name{is.pointer.ms}
\alias{is.pointer.ms}
\title{Data in R or C}
\usage{
  is.pointer.ms(x)
}
\arguments{
  \item{x}{MS object containing at least one sequence}
}
\value{
  TRUE if data in MS is stored in C, otherwise FALSE
  indicating it is stored in R
}
\description{
  Determine if data in MS object is stored in R (by value)
  or C (by reference). If data is stored in C, it can be
  copied to R using the function
  \code{\link{from.pointer.ms}}. Conversely, if the data is
  stored in R, it can be copied to C using the function
  \code{\link{as.pointer.ms}}.
}

