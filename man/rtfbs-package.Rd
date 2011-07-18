\name{rtfbs-package}
\alias{rtfbs-package}
\alias{rtfbs}
\docType{package}
\title{
R Transcription Factor Binding Site identification tool
}
\description{
RTFBS identifies and scores possible Transcription Factor
Binding Sites and allows for FDR analysis and pruning.  It supports
splitting of sequences based on size or a specified GFF, grouping
by GC content, and specification of Markov model order.  The heavy
lifting is done in C while all results are made available via R.
}
\details{
\tabular{ll}{
Copyright: \tab The code in src/pcre is Copyright (c) 1997-2010 University\cr
\tab of Cambridge.  All other code is Copyright (c) 2002-2011 University\cr
\tab of California, Cornell University.\cr
Package: \tab rtfbs\cr
License: \tab BSD\cr
Version: \tab 0.1.0\cr
URL: \tab http://compgen.bscb.cornell.edu/rtfbs\cr
Date: \tab 2011-07-18\cr
Depends: \tab stats, rphast\cr
Collate: \tab 'tfbs.R'\cr
Built: \tab R 2.13.1; x86_64-pc-linux-gnu; 2011-07-23 18:14:59 UTC; unix\cr
}

Index:
\preformatted{
[.ms                    Extract, replace, reorder MS...
as.pointer.ms           MS To Pointer
build.mm                Build Markov Model to represent sequences in an
                        MS object
calc.fdr                Calculate FDR
concat.ms               Concat MS
from.pointer.ms         MS From Pointer
gcContent.ms            Get GC content of each sequence in an MS
                        object...
groupByGC.ms            Group sequences by GC
is.pointer.ms           Data in R or C
length.ms               Length of MS object
lengths.ms              MS sequence lengths
ms                      Multiple Sequence (MS) Objects
names.ms                MS Sequence Names
offsets.ms              Get index offsets
output.sites            Threshold possible binding sites by Score or
                        FDR
plot.fdr                Plot FDR
print.ms                Printing MS objects
read.mm                 Read Markov Model from file
read.ms                 Reading in sequences from file
read.pwm                Read PWM object
score.ms                Score sequences against a PWM
sequences.ms            Get sequences
simulate.ms             Generate sequence from Markov Model
split.ms                Split sequences
summary.ms              MS Summary
write.mm                Write Markov Model to file
write.ms                Writing MS Object to FASTA file
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{vignette} \tab rtfbs (source)\cr
}
}
\author{
Nicholas Peterson, Andre Martins, Melissa Hubisz, and Adam
Siepel

Maintainer: RTFBS Group <rtfbs-help-l@cornell.edu>
}
\keyword{ package }
