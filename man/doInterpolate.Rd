% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFunctions.R
\name{doInterpolate}
\alias{doInterpolate}
\title{calculate the linear interpolation for each consecutive set of empty blocks}
\usage{
doInterpolate(dataDF, dateDF, timeframe, baselineDate)
}
\arguments{
\item{dataDF}{the data frame with the measurements records}

\item{dateDF}{corresponding data frame with date records}

\item{timeframe}{length of blocks in days}

\item{baselineDate}{date of the baseline visit}
}
\value{
data frame with (linearly) interpolated values
}
\description{
calculate the linear interpolation for each consecutive set of empty blocks
}
\keyword{internal}
