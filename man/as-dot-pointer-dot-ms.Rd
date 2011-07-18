\name{as.pointer.ms}
\alias{as.pointer.ms}
\title{MS To Pointer}
\usage{as.pointer.ms(src)
}
\description{RTFBS can store MS objects in R's memory, or C's memory for
efficency reasons.  This function transforms an MS object
whose data is stored in memory within R into an MS object
whose data is stored in C's memory.  The result is a new MS object in R
whose data is stored in C memory.
Copying an MS object from R into C is performed
automatically by RTFBS when it needs to run C code.  The user
will probably not need to call this function.}
\value{an MS object containing only a pointer to an object
created in C which contains the data.}
\seealso{\code{\link{ms}} for details on MS storage options.}
\keyword{ms}
\author{Nick Peterson}
\arguments{\item{src}{An MS object stored by value in R}
}

