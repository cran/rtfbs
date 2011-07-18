\name{ms}
\alias{ms}
\title{Multiple Sequence (MS) Objects}
\usage{ms(seqs, names=NULL, offsets=NULL, pointer.only=FALSE)
}
\description{Creates a new Multiple Sequences (MS) object to hold given sequences.}
\details{Make a new multiple sequence (MS) object given a vector of
character strings.  They can be optionally annotated with sample names.

The number of elements in names (if provided) must match the number 
of elements in seqs.

An alphabet (valid non-missing characters) of "ACGT" is automatically 
assumed for all sequences that RTFBS operates on.

About storing objects as pointers:
If \code{pointer.only==FALSE}, the MS object will be stored in R and can be
viewed and modified by base R code as well as RTFBS functions.
Setting \code{pointer.only=TRUE} will cause the object to be stored by
reference, as an external pointer to an object created by C code.  This
may be necessary to improve performance, but the object can then only
be viewed/manipulated via RTFBS functions.  Furthermore, if an object
is stored as a pointer, then its value is liable to be changed when
passed as an argument to a function.}
\value{An ms object.  These are stored in an array-like format,
so that they can be subsetted with the [] operator.}
\seealso{Functions for accessing/viewing ms objects:
\code{\link{sequences.ms}, \link{offsets.ms}, \link{length.ms},
\link{lengths.ms}, \link{names.ms}, \link{print.ms},
\link{write.ms}, \link{[.ms}, \link{as.pointer.ms}, \link{from.pointer.ms},
\link{is.pointer.ms}}}
\keyword{ms}
\author{Nick Peterson}
\arguments{\item{seqs}{A character vector containing sequences, one per sample}
\item{names}{A character vector identifying the sample name for each
sequence.  If \code{NULL}, use "seq1", "seq2", ...}
\item{offsets}{List of integers giving the offset for each sequences from the
start of its unsplit sequence.  If \code{NULL}, offsets of zero are assumed
for each sequence.}
\item{pointer.only}{a boolean indicating whether returned alignment object
should be stored by reference (see Details)}
}

