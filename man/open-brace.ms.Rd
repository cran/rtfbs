% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfbs.R
\name{[.ms}
\alias{[.ms}
\title{Extract, replace, reorder MS}
\usage{
\method{[}{ms}(x, rows)
}
\arguments{
\item{x}{An object of type \code{ms}}

\item{rows}{A numeric vector of sequence indices,
character vector (containing sequence name), or
logical vector (containing sequences to keep).  If logical vector it
will be recycled as necessary to the same length as \code{nrow.ms(x)}.}
}
\value{
An MS object sampled from x as indicated by rows.
}
\description{
Treat multiple sequences as a array where each row
corresponds to a sequence for one species.
}
\details{
The bracket notation can return a set of sequences,
or re-order rows.
}
\note{
This function will not alter the value of x even if it is stored as
a pointer to a C structure.
}
\author{
Nick Peterson
}
\keyword{ms}

