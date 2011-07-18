\name{gcContent.ms}
\alias{gcContent.ms}
\title{Get GC content of each sequence in an MS object...}
\usage{gcContent.ms(ms)
}
\description{Get GC content of each sequence in an MS object}
\value{numeric vector containing GC content of each sequence in ms.
GC content is computed as (# of GC bases)/(# of ACGT bases).}
\arguments{\item{ms}{MS object}
}

\examples{
require("rtfbs")
seqs <- ms(seqs=c("AAAA", "ACACACAC", "CGCCG", "ACGTACGTACGT", "CGGGGGGGGGG"),
           paste("fake", 1:5, sep=""))
gcContent.ms(seqs)
}
