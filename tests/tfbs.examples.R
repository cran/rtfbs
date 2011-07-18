library(rtfbs)

#' read.ms
require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
seqs <- read.ms(seqFile)
# Print the number of sequences in the MS object and whether
#   stored in C or R
print(seqs)

###############################################

#' read.pwm
require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
pwmFile <- "pwm.meme"
unzip(exampleArchive, pwmFile)
# Read in Position Weight Matrix (PWM) from MEME file from
#  the examples into a Matrix object
pwm <- read.pwm(pwmFile)
# Print PWM as an R matrix
print(pwm)

################################################

#' groupByGC.ms
require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
seqFile <- "input.fas"
unzip(exampleArchive, seqFile)
# Read in FASTA file "input.fas" from the examples into an 
#   MS (multiple sequences) object
seqs <- read.ms(seqFile)
# Group sequences from the "seqs" MS object based on each 
#   sequences's GC content into 4 new MS objects, one for
#   each GC content range
groups <- groupByGC.ms(seqs, 4)
sapply(groups, length)
sapply(groups, function(x) {range(gcContent.ms(x))})

################################################

#' gcContent.ms
require("rtfbs")
seqs <- ms(seqs=c("AAAA", "ACACACAC", "CGCCG", "ACGTACGTACGT", "CGGGGGGGGGG"),
           paste("fake", 1:5, sep=""))
gcContent.ms(seqs)


################################################

#' build.mm
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

#################################################

#' simulate.ms
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


##################################################

#' score.ms
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


################################################

#' calc.fdr
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

###############################################

#' plot.fdr
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
cs <- score.ms(ms, pwm, mm)
# Generate a sequence 1000 bases long using the supplied
#   Markov Model and random numbers
v <- simulate.ms(mm, 10000)
# Match the PWM against the sequences provided to find
#   possible transcription factor binding sites.  A 
#   Features object is returned, containing the location
#   of each possible binding site and an associated score.
#   Sites with a negative score are not returned unless 
#   we set threshold=-Inf as a parameter. Any identified
#   binding sites from simulated data are false positives
#   and used to calculate False Discovery Rate
xs <- score.ms(v, pwm, mm)
# Calculate the False Discovery Rate for each possible
#   binding site in the Features object CS.  Return
#   a mapping between each binding site score and the
#   associated FDR.
fdr <- calc.fdr(ms, cs, v, xs)
# Plot the False Discovery Rate v.s. score for one or
#   more groups.  To plot multiple FDR/Score mapping 
#   data frames on a single plot, simply supply a list
#   of FDR/Score data.frames 
plot.fdr(fdr)

##############################################

#' output.sites
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
cs <- score.ms(ms, pwm, mm)
# Generate a sequence 1000 bases long using the supplied
#   Markov Model and random numbers
v <- simulate.ms(mm, 10000)
# Match the PWM against the sequences provided to find
#   possible transcription factor binding sites.  A 
#   Features object is returned, containing the location
#   of each possible binding site and an associated score.
#   Sites with a negative score are not returned unless 
#   we set threshold=-Inf as a parameter. Any identified
#   binding sites from simulated data are false positives
#   and used to calculate False Discovery Rate
xs <- score.ms(v, pwm, mm)
# Calculate the False Discovery Rate for each possible
#   binding site in the Features object CS.  Return
#   a mapping between each binding site score and the
#   associated FDR.
fdr <- calc.fdr(ms, cs, v, xs)
# Output identified transcription factor binding sites 
#   below a 0.5 FDR threshold
output.sites(cs, NULL, fdr, 0.5)
# OR
# Output identified transcription factor binding sites
#   above a 5.2 score threshold 
output.sites(cs, 5.2)


################################################

rm(list = ls())
gc()
