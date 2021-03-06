% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisations.R
\name{visualiseOnePatientBinary}
\alias{visualiseOnePatientBinary}
\title{Creates a simple visualisation for serial and binary prescription data for one patient}
\usage{
visualiseOnePatientBinary(
  df,
  patient,
  idColumn = "PATIENT",
  blockColumn = "block",
  dataColumn = "measure",
  dateColumn = "date.measure",
  prescriptionColumn = "prescription"
)
}
\arguments{
\item{df}{the response data frame. It needs to contain binary prescription information.}

\item{patient}{The PATIENT ID number}

\item{idColumn}{name of ID column: default is PATIENT}

\item{blockColumn}{name of the block column: default is block}

\item{dataColumn}{name of the data column: default is measure}

\item{dateColumn}{name of date column: default is date.measure. This column has to be of class Date}

\item{prescriptionColumn}{name of the prescription column: default is prescription (it has to be binary)}
}
\value{
a visualisation
}
\description{
Creates a simple visualisation for serial and binary prescription data for one patient
}
\examples{
\dontrun{
visualiseOnePatientBinary(results,patient = 10090, idColumn = 'PATIENT', blockColumn = 'block',
              dataColumn = 'sbp', dateColumn = 'date.sbp', prescriptionColumn = 'anyRAASBinary')
}
}
