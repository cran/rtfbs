\name{summary.ms}
\alias{summary.ms}
\title{MS Summary}
\usage{summary.ms(object, ..., print.all=(length.ms(object) < 15 &&
    sum(lengths.ms(object)) < 500))
}
\description{Prints a short description of an Multiple Sequence (MS) object.
Omits names, offsets, and bases by default, but these can be printed 
using the print.all argument.}
\keyword{ms}
\seealso{\code{\link{print.ms}}}
\author{Nick Peterson}
\arguments{\item{object}{MS object}
\item{...}{Not used (exists for S3 compatibility)}
\item{print.all}{whether to suppress printing of the bases, offsets, and 
names.  If \code{TRUE} prints sequences, offsets, and names of all sequences.}
}

