\name{calc.fdr}
\alias{calc.fdr}
\title{Calculate FDR}
\usage{calc.fdr(realSeqs, realSeqsScores, simSeqs, simSeqsScores, interval=0.01)
}
\description{Calculate False Discovery Rate (FDR) of possible binding sites.
This function uses two sets of scores, realSeqsScores and simSeqsScores.  
realSeqsScores are scores for the sequences being scanned for binding sites.
simSeqsScores are scores for the simulated sequence.  The simulated 
sequences and simSeqsScores must be made using the same Markov Model as the realSeqsScores.}
\seealso{\code{ \link{score.ms}}}
\value{Data.Frame with two columns 'score' and 'FDR' mapping a single score to a single FDR.  Data frame is sorted by score if any exist.}
\note{realSeqsScores and simSeqsScores are both objects returned by
score.ms; the same arguments (threshold, conservative, strand) should be
used in both calls to score.ms or FDR will not be valid.}
\arguments{\item{realSeqs}{MS object containing non-simulated sequences}
\item{realSeqsScores}{Feat object obtained from scoring realSeqs}
\item{simSeqs}{MS object containing simulated sequences}
\item{simSeqsScores}{Feat object obtained from scoring simSeqs}
\item{interval}{Float specifying distance between steps at which the FDR will be calculated (lower is better). If NULL, calculate FDR for each unique score.}
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
cs <- score.ms(ms, pwm, mm, threshold=-2)
# Generate a sequence 1000 bases long using the supplied
#   Markov Model and random numbers
v <- simulate.ms(mm, 100000)
# Match the PWM against the sequences provided to find
#   possible transcription factor binding sites.  A 
#   Features object is returned, containing the location
#   of each possible binding site and an associated score.
#   Sites with a negative score are not returned unless 
#   we set threshold=-Inf as a parameter. Any identified
#   binding sites from simulated data are false positives
#   and used to calculate False Discovery Rate
xs <- score.ms(v, pwm, mm, threshold=-2)
# Calculate the False Discovery Rate for each possible
#   binding site in the Features object CS.  Return
#   a mapping between each binding site score and the
#   associated FDR.
fdr <- calc.fdr(ms, cs, v, xs)
# Print the Data.Frame containing the FDR/Score mapping
fdr
}
