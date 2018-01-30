\name{rtfbs-package}
\alias{rtfbs-package}
\alias{rtfbs}
\docType{package}
\title{
Transcription Factor Binding Site Identification Tool
}
\description{
Identifies and scores possible Transcription Factor
    Binding Sites and allows for FDR analysis and pruning.  It supports
    splitting of sequences based on size or a specified GFF, grouping
    by G+C content, and specification of Markov model order.  The heavy
    lifting is done in C while all results are made available via R.
}
\details{
\tabular{ll}{
Copyright: \tab Copyright (c) 2002-2018 University of California, Cornell\cr
\tab University, Cold Spring Harbor Laboratory.\cr
Package: \tab rtfbs\cr
License: \tab BSD_3_clause + file LICENSE\cr
Version: \tab 0.3.9\cr
URL: \tab http://compgen.cshl.edu/rtfbs\cr
Date: \tab 2018-01-30\cr
Imports: \tab \cr
\tab rphast,\cr
\tab methods\cr
Depends: \tab \cr
\tab stats\cr
Collate: \tab \cr
\tab 'tfbs.R'\cr
RoxygenNote: \tab 5.0.1\cr
Built: \tab R 3.1.2; x86_64-unknown-linux-gnu; 2018-01-30 20:25:08 UTC; unix\cr
}

Index:
\preformatted{
[.ms                    Extract, replace, reorder MS
as.pointer.ms           MS To Pointer
build.mm                Build Markov Model to represent sequences in an
                        MS object
calc.fdr                Calculate FDR
concat.ms               Concat MS
from.pointer.ms         MS From Pointer
gcContent.ms            Get GC content of each sequence in an MS object
groupByGC.ms            Group sequences by GC
is.pointer.ms           Data in R or C
label.matrix            Name PWM & MM rows and columns
length.ms               Length of MS object
lengths.ms              MS sequence lengths
makeFdrPlot             Plot FDR
ms                      Multiple Sequence (MS) Objects
names.ms                MS Sequence Names
offsets.ms              Get index offsets
output.sites            Threshold possible binding sites by Score or
                        FDR
print.ms                Prints an MS (multiple sequence) object. #'
                        @title Printing MS objects
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
}
\author{
Nicholas Peterson, Andre Martins, Melissa Hubisz, and Adam
    Siepel

Maintainer: Ritika Ramani<rramani@cshl.edu>
}
\keyword{ package }
