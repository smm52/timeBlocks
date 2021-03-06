% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/responseScore.R
\name{calculateResponseScore}
\alias{calculateResponseScore}
\title{Calculates a score based on principle component analysis}
\usage{
calculateResponseScore(
  responseDf,
  idColumn = "PATIENT",
  blockColumn = "block",
  sexColumn = NA,
  responseColumns,
  confounderColumns = NA
)
}
\arguments{
\item{responseDf}{the data frame that contains the time blocks and the biomarker data per time block}

\item{idColumn}{name of ID column}

\item{blockColumn}{name of time block column}

\item{sexColumn}{name of column that contains sex information (optional). If provided, the score will be adjusted for sex.}

\item{confounderColumns}{names of columns that contain confounders (optional). If provided, the score will be adjusted for these confounders using linear regression}

\item{responseColumn}{names of the columns containing the patient block value for a biomarker (e.g. vector of string names). The score will be based only on these columns.}
}
\value{
the response data frame with an added column phenoScore that shows the computed score
}
\description{
The function calcuculates a score for biomarker data represented as time blocks. The score is calculated via PCA.
}
\examples{
\dontrun{
results <- read_tsv('../data/responseSpreadsheet.txt')

resultsCor <- calculateResponseScore(results,responseColumns = c('sbp','dbp','aerAll','acrAll','egfr'), sexColumn = 'SEX')
}
}
