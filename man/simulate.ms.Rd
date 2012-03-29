\name{simulate.ms}
\alias{simulate.ms}
\title{Generate sequence from Markov Model}
\usage{
  \method{simulate}{ms} (object, nsim, seed = NULL,
    pointer.only = FALSE, ...)
}
\arguments{
  \item{object}{Markov Model \code{\link{build.mm}}}

  \item{nsim}{Length of the sequence to simulate}

  \item{seed}{A random number seed.  Either \code{NULL}
  (the default; do not re-seed random number generator), or
  an integer to be sent to set.seed.}

  \item{pointer.only}{If \code{TRUE}, keep sequence data
  stored in a C structure, otherwise it is automatically
  copied into an R object.}

  \item{...}{Not used; for S3 compatibility}
}
\value{
  MS object containing a single sequence with nsim bases.
}
\description{
  Simulate a single sequence based from a Markov Model.
  These are referred to as simulated sequences and used
  compute the background rates and False Discovery Rates.
}
\examples{
require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
ms <- read.ms(seqFile);
# Build a 3rd order Markov Model to represent the sequences
#   in the MS object "ms".  The Model will be a list of
#   matrices  corrisponding in size to the order of the 
#   Markov Model
mm <- build.mm(ms, 3);
# Generate a sequence 1000 bases long using the supplied
#   Markov Model and random numbers
v <- simulate.ms(mm, 1000);
}
\seealso{
  \code{\link{build.mm}} for details on Markov models,
  \code{\link{ms}} for details on MS objects
}

