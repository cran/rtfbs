\name{offsets.ms}
\alias{offsets.ms}
\title{Get index offsets}
\usage{
  offsets.ms(ms)
}
\arguments{
  \item{ms}{An MS object}
}
\value{
  List of integers, each indicating the offset at which
  this sequence starts compared to the unsplit sequence it
  came from
}
\description{
  Return list of index offsets for each sequence in an MS
  object
}

