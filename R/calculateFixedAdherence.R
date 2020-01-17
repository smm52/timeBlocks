#' Calculates medication adherence based on a fixed, specified refill time period
#'
#' @param serialDf the data frame with the serial data
#' @param startDates dataframe with ID and a date to start the adherence calculation (optional)
#' @param atcCode regular expression for the class of medication ATC code (default: C09 which stands for any RAAS)
#' @param refillPeriod length of the refill period (default 90 days)
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the column with the ATC codes: default is ATC
#' @return time block data and dates in wide format
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate

calculateFixedAdherence <- function(serialDf, startDates = NA, atcCode = "C09", refillPeriod = 90, 
                                    idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC") {
  # check data frames
  if (nrow(serialDf) == 0 ) {
    stop("Serial prescription data is empty")
  } else {
    serialDf <- serialDf[ , colSums(is.na(serialDf)) == 0]
    if (nrow(serialDf) == 0 ) {
      stop("Serial prescription data is empty after removing missing data")
    }
  }
  
  # check block length
  if (refillPeriod <= 0) {
    stop("The refill period cannot be 0 or negative")
  }
  
  # check baseline data
  if(!is.null(nrow(startDates))){
    startDate <- checkBaselineFormat(startDates, idColumn = idColumn, dateColumn = dateColumn)
    if (is.null(startDate)) {
      stop("incorrect format of start dates")
    }
  }
  
  
  # check serial prescription data
  serialDf <- checkBinaryPrescriptionFormat(serialDf, idColumn = idColumn, dateColumn = dateColumn, atcColumn = atcColumn)
  if (is.null(serialDf)) {
    stop("incorrect format of serial prescription input data")
  }
  # restrict serial data to patients in startDates (if available)
  if(!is.null(nrow(startDates))){
    serialDf <- serialDf %>% filter(PATIENT %in% startDates$PATIENT)
    if (nrow(serialDf) == 0) {
      stop("Serial prescription data does not contain any data for the patients.")
    }
  }
  
  # restrict serial to ATC code under investigation
  serialDf <- serialDf %>% filter(grepl(atcCode, ATC))
  if (nrow(serialDf) == 0) { 
    stop('Serial data contains no prescriptions with the specific ATC code expression.') 
  }
  
  #if start dates are available, start looking one refill period back to see whether the patient was medicated at start
  if(!is.null(nrow(startDates))){
    startDates <- startDates %>%
      mutate(VISIT = VISIT - refillPeriod)
  }
  
  #create results data frame
  resultsDf <- data.frame(PATIENT=character(),adherence=double(),firstPrescription=character(),lastPrescription=character(),
                          lastCoverDate=character(), numPrescriptions=numeric(), stringsAsFactors = F)
  
  allPatients <- unique(serialDf$PATIENT)
  #go through all patients
  for(i in 1:length(allPatients)){
    myPatient <- allPatients[i]
    resultsDf[i,'PATIENT']  <- myPatient
    myPrescriptions <- serialDf %>%
      filter(PATIENT == myPatient)
    #restrict prescription for this patient to start at a certain time (optional)
    if(!is.null(nrow(startDates))){
      myStartDate <- (startDates %>% filter(PATIENT == myPatient))$VISIT
      myPrescriptions <- myPrescriptions %>%
        filter(VISIT >= myStartDate)
    }
    #if there are 0 prescriptions, there is 0 adherence, if there is 1 prescription, the adhrence period is the refill period
    if(nrow(myPrescriptions) > 1){
      resultsDf[i,'firstPrescription'] <- as.character(as.Date(min(myPrescriptions$VISIT)))
      resultsDf[i,'lastPrescription'] <- as.character(as.Date(max(myPrescriptions$VISIT)))
      resultsDf[i,'lastCoverDate'] <- as.character(as.Date(max(myPrescriptions$VISIT)) + refillPeriod)
      resultsDf[i,'numPrescriptions'] <- nrow(myPrescriptions)
      daysToCover <- as.numeric(max(myPrescriptions$VISIT) - min(myPrescriptions$VISIT)) + refillPeriod
      #go through the precsriptions in a sorted way starting from the earliest
      myPrescriptions <- myPrescriptions[order(myPrescriptions$VISIT),]
      uncoveredDays <- 0
      for(j in 2:nrow(myPrescriptions)){
        pOld <- myPrescriptions[(j-1),'VISIT'][[1]]
        pNew <- myPrescriptions[j,'VISIT'][[1]]
        if(as.numeric(pNew - pOld) > refillPeriod){
          uncoveredDays <- uncoveredDays + as.numeric((pNew - pOld)) - refillPeriod
        }  
      }
      resultsDf[i,'adherence'] <- (daysToCover - uncoveredDays)/daysToCover
    } else if(nrow(myPrescriptions) == 0){
      resultsDf[i,'adherence'] <- 0
      resultsDf[i,'firstPrescription'] <- NA
      resultsDf[i,'lastPrescription'] <- NA
      resultsDf[i,'lastCoverDate'] <- NA
      resultsDf[i,'numPrescriptions'] <- 0
    } else if(nrow(myPrescriptions) == 1){
      resultsDf[i,'adherence'] <- 1
      resultsDf[i,'firstPrescription'] <- as.character(as.Date(min(myPrescriptions$VISIT)))
      resultsDf[i,'lastPrescription'] <- as.character(as.Date(min(myPrescriptions$VISIT)))
      resultsDf[i,'lastCoverDate'] <- as.character(as.Date(min(myPrescriptions$VISIT)) + refillPeriod)
      resultsDf[i,'numPrescriptions'] <- 1
    }
  }
  resultsDf$firstPrescription <- as.Date(resultsDf$firstPrescription)
  resultsDf$lastPrescription <- as.Date(resultsDf$lastPrescription)
  resultsDf$lastCoverDate <- as.Date(resultsDf$lastCoverDate)
  resultsDf
}

