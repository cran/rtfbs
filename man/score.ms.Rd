\name{score.ms}
\alias{score.ms}
\title{Score sequences against a PWM}
\usage{
  score.ms(ms, pwm, mm, conservative = TRUE, threshold = 0,
    strand = "best")
}
\arguments{
  \item{ms}{MS object containing at least one sequence}

  \item{pwm}{Position Weight matrix representing
  transcription factor motif}

  \item{mm}{Markov Model associated with given sequences,
  which represents the null model}

  \item{conservative}{(Logical value) If TRUE, sequences
  containing N's are given a log likelihood of negative
  infinity under the PWM model. If FALSE, any 'N'
  encountered does not contributes to the score.}

  \item{threshold}{(Numeric value) Only sites with scores
  above this threshold are returned (default = 0)}

  \item{strand}{One of "best", "both", "+", or "-"
  specifying which strand(s) to return results for. If
  "both" search for binding sites in both directions,
  return all results found. If "best" search for binding
  sites in both directions, but for each N-mer, return the
  maximum score over either strand. If "+" look only on the
  forward strand, and if "-" look only on the reverse
  strand.}
}
\value{
  Scores and locations for possible binding sites returned
  as a feature object.
}
\description{
  Score all potential binding sites in an MS object. If a
  PWM has N rows, then score every observed N-mer in the MS
  object. The score is given by the log likelihood of the
  N-mer given the PWM, minus the log likelihood of the
  N-mer under the Markov model specified by mm. By default,
  only potential binding sites with scores > 0 are
  returned, but this can be modified with the
  \code{threshold} argument.
}
\note{
  If a PWM file contains multiple PWMs, then read.pwm will
  return a list of PWMs.  This function takes a single PWM.
}
\examples{
require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
ms <- read.ms(seqFile);
pwmFile <- "pwm.meme"
unzip(exampleArchive, pwmFile)
# Read in Position Weight Matrix (PWM) from MEME file from
#  the examples into a Matrix object
pwm <- read.pwm(pwmFile)
# Build a 3rd order Markov Model to represent the sequences
#   in the MS object "ms".  The Model will be a list of
#   matrices  corrisponding in size to the order of the 
#   Markov Model
mm <- build.mm(ms, 3);
# Match the PWM against the sequences provided to find
#   possible transcription factor binding sites.  A 
#   Features object is returned, containing the location
#   of each possible binding site and an associated score.
#   Sites with a negative score are not returned unless 
#   we set threshold=-Inf as a parameter.
score.ms(ms, pwm, mm)
}
\seealso{
  \code{\link{read.ms} \link{split.ms} \link{groupByGC.ms}
  \link{build.mm} \link{read.pwm}}
}

