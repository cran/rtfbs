\name{groupByGC.ms}
\alias{groupByGC.ms}
\title{Group sequences by GC}
\usage{groupByGC.ms(ms, ngroups)
}
\description{Group sequences in an MS object by their GC content.}
\seealso{\code{\link{read.ms}}}
\value{List of MS objects, where element i represents the i'th quantile
according to GC content.}
\arguments{\item{ms}{MS object, containing at least ngroups sequences.}
\item{ngroups}{Number of quantiles to group sequences into.}
}

\examples{
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
}
