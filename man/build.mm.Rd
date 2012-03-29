\name{build.mm}
\alias{build.mm}
\title{Build Markov Model to represent sequences in an MS object}
\usage{
  build.mm(ms, order, pseudoCount = 0,
    considerReverse = FALSE)
}
\arguments{
  \item{ms}{Sequence used to build Markov Model}

  \item{order}{Order of Markov model to build; ie, the
  number of preceding bases to consider when calculating
  the probability of a given base.}

  \item{pseudoCount}{(Optional) Integer added to the number
  of observed cases of each possible sequence pattern.}

  \item{considerReverse}{(Optional) Logical value.  If
  \code{TRUE}, considers reverse complement frequencies in
  addition to forward strand frequencies when building
  model.}
}
\value{
  A list of matrices, each representing a markov model from
  order 0, 1, ..., order.  Each matrix gives the
  probability of observing a particular base (column) given
  the preceding bases (row).
}
\description{
  Build a Markov Model of user specified order to represent
  sequences in an MS object.
}
\examples{
require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
ms <- read.ms(seqFile)
# Build a 3rd order Markov Model to represent the sequences
#   in the MS object "ms".  The Model will be a list of
#   matrices  corrisponding in size to the order of the 
#   Markov Model
mm <- build.mm(ms, 3);
# Print the list of 4 Markov Matrices that make up the 
#   Markov Model 
print(mm)
}
\seealso{
  \code{\link{read.ms} \link{split.ms} \link{groupByGC.ms}}
}

