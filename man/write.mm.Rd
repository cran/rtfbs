\name{write.mm}
\alias{write.mm}
\title{Write Markov Model to file}
\usage{
  write.mm(x, file = NULL, ...)
}
\arguments{
  \item{x}{Markov Model (list of Markov Matrices, one for
  each order)}

  \item{file}{Full path of where to save file, if NULL
  print to screen}

  \item{...}{Not Used (for S3 compatability)}
}
\description{
  Write markov model into the format used by MEME's
  fasta-get-markov tool
}

