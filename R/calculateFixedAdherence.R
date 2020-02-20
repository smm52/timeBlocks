#' Calculates medication adherence based on a fixed, specified refill time period
#'
#' @param serialDf the data frame with the serial data
#' @param startDates dataframe with ID and a date to start the adherence calculation (optional)
#' @param endDates dataframe with ID and a date to end the adherence calculation (optional)
#' @param atcCode regular expression for the class of medication ATC code (default: C09 which stands for any RAAS)
#' @param refillPeriod length of the refill period (default 90 days)
#' @param stopPeriod number of days after which the medication is considered to have stopped (optional)
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the column with the ATC codes: default is ATC
#' @param doPCAScore flag whether a PCA score from the medication data should be calculated
#' @return time block data and dates in wide format
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr select_if
#' @importFrom dplyr starts_with
#' @importFrom dplyr ends_with

calculateFixedAdherence <- function(serialDf, startDates = NA, endDates = NA, atcCode = "C09", refillPeriod = 90, stopPeriod = NA,
                                    idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC", doPCAScore = F) {
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
  
  if(!is.na(stopPeriod)){
    if(stopPeriod <= refillPeriod){
      stop("Stop period needs to be bigger than the refill period")
    }
  }
  
  # check baseline data
  if(!is.null(nrow(startDates))){
    startDate <- checkBaselineFormat(startDates, idColumn = idColumn, dateColumn = dateColumn)
    if (is.null(startDate)) {
      stop("incorrect format of start dates. Check column class. No missing values allowed.")
    }
  }
  
  #chek end date data
  if(!is.null(nrow(endDates))){
    endDate <- checkBaselineFormat(endDates, idColumn = idColumn, dateColumn = dateColumn)
    if (is.null(endDate)) {
      stop("incorrect format of end dates. Check column class. No missing values allowed.")
    }
  }
  
  
  
  # check serial prescription data
  serialDf <- checkBinaryPrescriptionFormat(serialDf, idColumn = idColumn, dateColumn = dateColumn, atcColumn = atcColumn)
  if (is.null(serialDf)) {
    stop("incorrect format of serial prescription input data. Check column class. No missing values allowed.")
  }
  
  # check start and end date patients are the same
  if(!is.null(nrow(startDates)) & !is.null(nrow(endDates))){
    if(!all(startDates$PATIENT %in% endDates$PATIENT) & !all(endDates$PATIENT %in% startDates$PATIENT)){
      stop("start and end dates must be specified for the same patients")  
    }  
  }
  
  # restrict serial data to patients in startDates (if available)
  if(!is.null(nrow(startDates))){
    serialDf <- serialDf %>% filter(PATIENT %in% startDates$PATIENT)
    if (nrow(serialDf) == 0) {
      stop("Serial prescription data does not contain any data for the patients.")
    }
  }
  
  # restrict serial data to patients in endDates (if available)
  if(!is.null(nrow(endDates))){
    serialDf <- serialDf %>% filter(PATIENT %in% endDates$PATIENT)
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
                          lastCoverDate=character(), numPrescriptions=numeric(), startStops=numeric(), stoppedYears = double(),
                          coveredYears = double(), stringsAsFactors = F)
  
  #if there are start dates, use the start date patients to iterate over to find adherence, otherwise use all patients in the prescriptions
  if(!is.null(nrow(startDates))){
    allPatients <- unique(startDates$PATIENT)   
  } else {
    allPatients <- unique(serialDf$PATIENT)  
  }
  #go through all patients
  for(i in 1:length(allPatients)){
    myPatient <- allPatients[i]
    resultsDf[i,'PATIENT']  <- myPatient
    myPrescriptions <- serialDf %>%
      filter(PATIENT == myPatient)
    #restrict prescription for this patient to start at a certain time (optional)
    if(nrow(myPrescriptions) > 0){
      if(!is.null(nrow(startDates))){
        myStartDate <- (startDates %>% filter(PATIENT == myPatient))$VISIT
        myPrescriptions <- myPrescriptions %>%
          filter(VISIT >= myStartDate)
      }
      if(!is.null(nrow(endDates))){
        myEndDate <- (endDates %>% filter(PATIENT == myPatient))$VISIT
        myPrescriptions <- myPrescriptions %>%
          filter(VISIT <= myEndDate)
      }
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
      stoppedDays <- 0
      numStops <- 0
      for(j in 2:nrow(myPrescriptions)){
        pOld <- myPrescriptions[(j-1),'VISIT'][[1]]
        pNew <- myPrescriptions[j,'VISIT'][[1]]
        if(as.numeric(pNew - pOld) > refillPeriod){
          if(is.na(stopPeriod)){
            uncoveredDays <- uncoveredDays + as.numeric((pNew - pOld)) - refillPeriod
          } else if(as.numeric(pNew - pOld) > (stopPeriod + refillPeriod)){
            stoppedDays <- stoppedDays + as.numeric((pNew - pOld))
            numStops <- numStops + 1
            daysToCover <- daysToCover - as.numeric((pNew - pOld)) 
          } else {
            uncoveredDays <- uncoveredDays + as.numeric((pNew - pOld)) - refillPeriod  
          }
        }  
      }
      resultsDf[i,'adherence'] <- (daysToCover - uncoveredDays)/daysToCover
      resultsDf[i,'coveredYears'] <- (daysToCover - uncoveredDays)/365.25
      resultsDf[i,'startStops'] <- numStops
      resultsDf[i,'stoppedYears'] <- stoppedDays/365.25
    } else if(nrow(myPrescriptions) == 0){
      resultsDf[i,'adherence'] <- 0
      resultsDf[i,'firstPrescription'] <- NA
      resultsDf[i,'lastPrescription'] <- NA
      resultsDf[i,'lastCoverDate'] <- NA
      resultsDf[i,'numPrescriptions'] <- 0
      resultsDf[i,'startStops'] <- NA
      resultsDf[i,'stoppedYears'] <- NA
      resultsDf[i,'coveredYears'] <- 0
    } else if(nrow(myPrescriptions) == 1){
      resultsDf[i,'adherence'] <- 1
      resultsDf[i,'firstPrescription'] <- as.character(as.Date(min(myPrescriptions$VISIT)))
      resultsDf[i,'lastPrescription'] <- as.character(as.Date(min(myPrescriptions$VISIT)))
      resultsDf[i,'lastCoverDate'] <- as.character(as.Date(min(myPrescriptions$VISIT)) + refillPeriod)
      resultsDf[i,'numPrescriptions'] <- 1
      resultsDf[i,'startStops'] <- 0
      resultsDf[i,'stoppedYears'] <- 0
      resultsDf[i,'coveredYears'] <- refillPeriod/365.25
    }
  }
  resultsDf$firstPrescription <- as.Date(resultsDf$firstPrescription)
  resultsDf$lastPrescription <- as.Date(resultsDf$lastPrescription)
  resultsDf$lastCoverDate <- as.Date(resultsDf$lastCoverDate)
  
  if(!is.null(nrow(endDates))){
    resultsDf <- resultsDf %>%
      merge(endDates %>% select(PATIENT, VISIT))
    resultsDf <- resultsDf %>%
      mutate(gapEndFollowUpYears = safe.ifelse(lastCoverDate < VISIT,as.numeric(VISIT - lastCoverDate)/365.25,0)) %>%
      select(-VISIT)
  }
  
  if(doPCAScore){
    scoring <- resultsDf %>% 
      select(-ends_with('Prescription'), -ends_with('Date'), -numPrescriptions) 
    scoring[is.na(scoring)] <- 0
    s <- prcomp(scoring %>% select(-PATIENT, -starts_with('ATC')) %>% select_if(~ length(unique(.)) > 1), center = TRUE,scale. = TRUE)
    s <- scale(s$x[,1])
    resultsDf <- cbind(resultsDf,s) %>%
      as.data.frame()
    names(resultsDf)[which(names(resultsDf) == 's')] <- 'score'
  }
  
  resultsDf
}

#' Calculates medication adherence based on a fixed, specified refill time period
#'
#' @param serialDf the data frame with the serial data
#' @param startDates dataframe with ID and a date to start the adherence calculation (optional)
#' @param endDates dataframe with ID and a date to end the adherence calculation (optional)
#' @param atcCode list of regular expressions for the classes of medication ATC code (default: c('C09', 'C10'))
#' @param refillPeriod length of the refill period (default 90 days)
#' @param stopPeriod number of days after which the medication is considered to have stopped (optional)
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the column with the ATC codes: default is ATC
#' @param doPCAScore flag whether a PCA score from the medication data should be calculated
#' @param combinedScore flag whether a combined PCA scores for all ATC codes in the list should be calculated
#' @return time block data and dates in wide format
#' @export
#' @importFrom magrittr %>%
#' @importFrom plyr ldply
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select
#' @importFrom dplyr select_if
#' @importFrom dplyr starts_with
#' @importFrom dplyr ends_with

calculateFixedAdherenceList <- function(serialDf, startDates = NA, endDates = NA, atcCodeList = c("C09", "C10"), refillPeriod = 90, stopPeriod = NA,
                                    idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC", doPCAScore = F, combinedScore = T) {
  results <- atcCodeList %>%
    lapply(listAllAdherences, serialDf, startDates, endDates, refillPeriod, stopPeriod, idColumn, dateColumn, atcColumn) %>%
    ldply()
  
  if(doPCAScore){
    scoring <- results %>% 
      select(-ends_with('Prescription'), -ends_with('Date'), -numPrescriptions) 
    scoring <- scoring %>%
      pivot_wider(id_cols = PATIENT, names_from = ATC, 
                  values_from = names(scoring)[which(names(scoring) != 'PATIENT')])
    scoring[is.na(scoring)] <- 0
    if(combinedScore){
      s <- prcomp(scoring %>% select(-PATIENT, -starts_with('ATC')) %>% select_if(~ length(unique(.)) > 1), center = TRUE,scale. = TRUE)
      s <- scale(s$x[,1])
      results <- cbind(results,s) %>%
        as.data.frame()
      names(results)[which(names(results) == 's')] <- 'score'
    } else {
      results[,'score'] <- NA
      codes <- unique(results$ATC)
      for(i in 1:length(codes)){
        s <- prcomp(scoring %>% select(ends_with(paste0('_',codes[i]))) %>% select_if(~ length(unique(.)) > 1), center = TRUE,scale. = TRUE)
        s <- scale(s$x[,1])
        results[which(results$ATC == codes[i]),'score'] <- s
      }
    }
  }
  
  results
  
}

#' Helper function to call individual (one ATC code) adherence calculation for a list of ATC codes
#'
#' @param atcCode regular expression for the class of medication ATC code
#' @param serialDf the data frame with the serial data
#' @param startDates dataframe with ID and a date to start the adherence calculation (optional)
#' @param endDates dataframe with ID and a date to end the adherence calculation (optional)
#' @param refillPeriod length of the refill period (default 90 days)
#' @param stopPeriod number of days after which the medication is considered to have stopped (optional)
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the column with the ATC codes: default is ATC
#' @return time block data and dates in wide format
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
listAllAdherences <- function(atcCode, serialDf, startDates, endDates, refillPeriod, stopPeriod, idColumn, dateColumn, atcColumn){
  results <- calculateFixedAdherence(serialDf = serialDf, startDates = startDates, endDates = endDates, atcCode = atcCode, 
                                     refillPeriod = refillPeriod, stopPeriod =  stopPeriod, idColumn = idColumn, dateColumn = dateColumn,
                                     atcColumn = atcColumn, doPCAScore = F) %>%
    mutate(ATC = atcCode)
  results
  
}
