#' Checks the format of the data set quick visualisation
#'
#' @param df the response data frame. It needs to contain binary prescription information.
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is date.measure. This column has to be of class Date
#' @param dataColumn name of the data column: default is measure
#' @param blockColumn name of the block column: default is block
#' @param prescriptionColumn name of the prescription column: default is prescription (it has to be binary)
#' @return NULL if there are format problems, otherwise a data frame that has a PATIENT, VISIT, MEASURE, BLOCK and PRESCRIPTION column
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkVizFormat <- function(df, idColumn = 'PATIENT', blockColumn = 'block', dataColumn = 'measure', dateColumn = 'date.measure', prescriptionColumn = 'prescription') {
  result <- NULL
  if (class(df %>% extract2(dateColumn)) == "Date" & all(na.omit(df %>% extract2(prescriptionColumn)) %in% 0:1)) {
    names(df)[which(names(df) == idColumn)] <- "PATIENT"
    names(df)[which(names(df) == dateColumn)] <- "VISIT"
    names(df)[which(names(df) == dataColumn)] <- "MEASURE"
    names(df)[which(names(df) == blockColumn)] <- "BLOCK"
    names(df)[which(names(df) == prescriptionColumn)] <- "PRESCRIPTION"
    result <- df
  }
  result
}

#' Creates a simple visualisation for serial and binary prescription data for one patient
#'
#' @param df the response data frame. It needs to contain binary prescription information.
#' @param patient The PATIENT ID number
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is date.measure. This column has to be of class Date
#' @param dataColumn name of the data column: default is measure
#' @param blockColumn name of the block column: default is block
#' @param prescriptionColumn name of the prescription column: default is prescription (it has to be binary)
#' @return a visualisation
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 alpha
#' @importFrom ggplot2 theme_bw
#' @examples
#' \dontrun{
#' visualiseOnePatientBinary(results,patient = 10090, idColumn = 'PATIENT', blockColumn = 'block',
#'               dataColumn = 'sbp', dateColumn = 'date.sbp', prescriptionColumn = 'anyRAASBinary')
#' }
visualiseOnePatientBinary <- function(df, patient, idColumn = 'PATIENT', blockColumn = 'block', dataColumn = 'measure', dateColumn = 'date.measure', prescriptionColumn = 'prescription'){
  df <- checkVizFormat(df, idColumn = idColumn, blockColumn = blockColumn, dataColumn = dataColumn, dateColumn = dateColumn, prescriptionColumn = prescriptionColumn)
  if(is.null(df)){
    stop('Result data frame is incorrectly formatted')
  }
  df <- df %>% select(PATIENT, BLOCK, MEASURE, VISIT,PRESCRIPTION) %>% filter(PATIENT == patient)
  if(nrow(df) == 0){
    stop('Patient not found')
  }
  p <- ggplot( df, aes(x = BLOCK, y = MEASURE)) +
    geom_point(aes(shape = factor(ifelse(is.na(VISIT),1,2), levels = c(1,2))), size = 3) +
    geom_rect(aes(x = NULL, y = NULL, xmin = as.numeric(BLOCK) - 0.5, xmax = as.numeric(BLOCK) + 0.5, fill = factor(ifelse(is.na(PRESCRIPTION),2,PRESCRIPTION), levels = c(0,1,2))), ymin = -Inf, ymax = +Inf, data = df) +
    scale_shape_manual(values=c(19, 1), labels=c("estimated", "measured"), name = 'Biomarker', drop = FALSE) +
    geom_line() +
    scale_fill_manual(values = alpha(c("#1b9e77", "#d95f02","#7570b3"), .3), labels=c("No prescriptions", prescriptionColumn, "n/a"), name = 'Prescription', drop = FALSE) +
    theme_bw()
  plot(p)
}
