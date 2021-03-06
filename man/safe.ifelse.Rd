% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFunctions.R
\name{safe.ifelse}
\alias{safe.ifelse}
\title{Type safe if.else implementation}
\usage{
safe.ifelse(cond, yes, no)
}
\arguments{
\item{cond}{the condition}

\item{yes}{what happens if yes condition is met}

\item{no}{what happens if no condition is not met}
}
\value{
type safe value
}
\description{
Type safe if.else implementation that make sure that the output has the same class as the yes condition part.
This is important for dates for examples.
}
\examples{
myDate1 <- NA
myDate2 <- as.Date('2000-01-31', '\%Y-\%m-\%d')
safe.ifelse(is.na(myDate1),myDate2,myDate1)
class(safe.ifelse(is.na(myDate1),myDate2,myDate1))

ifelse(is.na(myDate1),myDate2,myDate1)
class(ifelse(is.na(myDate1),myDate2,myDate1))
}
