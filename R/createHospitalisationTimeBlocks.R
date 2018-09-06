#' Calculates block data for hospitalisations
#'
#' @param serialDf the data frame with the hospital events (one row per event)
#' @param baselineInfo dataframe with ID and baseline date
#' @param studyStartDate start Date of the study (has to be of class Date)
#' @param studyEndDate end Date of the study (has to be of class Date)
#' @param blockDays length of a time block in days (default 365.25 days)
#' @param blDayDiff time in days that should be searched for a measurement around the baseline visit (default 90)
#' @param idColumn name of ID column: default is PATIENT
#' @param startDateColumn name of start date column (day of hospitalisation or outpatient visit): default is VISIT_START. This column has to be of class Date
#' @param endDateColumn name of end date column (release date): default is VISIT_END. This column has to be of class Date
#' @param countNights flag indicating whether to count the days in hospital or the night: default TRUE, count nights (if nights are counted, an outpatient stay has length 0, otherwise 1)
#' @param blDateColumn name of the column containing the baseline date (date has to be of class Date): default is VISIT
#' @param longFormat flag to indicate whether results should be in long or wide format. Wide has some additional columns (default: TRUE)
#' @return time block data and dates in wide format
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom plyr ldply
#' @importFrom dplyr select
#' @importFrom dplyr rowwise
createHospitalTimeBlocks <- function(serialDf, baselineInfo, studyStartDate, studyEndDate,
                                     blockDays = 365.25, blDayDiff = 90, idColumn = "PATIENT", startDateColumn = "VISIT_START", 
                                     endDateColumn = "VISIT_END", countNights = TRUE, blDateColumn = "VISIT", longFormat = TRUE) {
  # check data frames
  if (nrow(serialDf) == 0 | nrow(baselineInfo) == 0) {
    stop("Serial data or baseline data is empty")
  }
  # check study dates
  if (class(studyStartDate) != "Date" | class(studyEndDate) != "Date") {
    stop("Study start and end date have to be of class Date")
  }
  if (studyEndDate <= studyStartDate) {
    stop("Start date after end date")
  }
  # check block length
  if (blockDays < 30) {
    stop("A block has to be at least 30 days")
  }
  # check days around baseline
  if (blDayDiff < 0 | blDayDiff >= blockDays) {
    stop("The number days around the baseline visit to check for a measurement has to be positive and it has to be smaller than a block.")
  }
  # restrict serial data to baseline patients
  serialDf <- serialDf %>% filter(PATIENT %in% baselineInfo$PATIENT)
  
  # check format and find baseline hospital stay in serial data
  baselineInfo <- findBaselineHospitalStay(serialDf, baselineDates = baselineInfo, blDayDiff = blDayDiff, idColumn = idColumn,
                                           startDateColumn = startDateColumn,  endDateColumn = endDateColumn, blDateColumn = blDateColumn, countNights = countNights)
  # restrict serial data to study period
  serialDf <- serialDf %>% filter(studyStartDate <= VISIT_END & studyEndDate >= VISIT_START)
  
  allRecords <- list()
  dataName <- ifelse(countNights == TRUE, "HOSPITAL.NIGHTS", "HOSPITAL.DAYS")
  # go through all baseline patients
  for (i in 1:nrow(baselineInfo)) {
    myPatient <- baselineInfo[i, ]
    
    #date of the baseline study visit
    myBaselineDate <- myPatient$VISIT
    
    #initialise lastMeasure with baseline measure
    lastMeasure <- ifelse(countNights == TRUE, 
                          myPatient$BL_HOSPITAL_NIGHTS,
                          myPatients$BL_HOSPITAL_DAYS)
    minContainerStart <- myBaselineDate
    
    # initial lists and set baseline (measure0) date and value
    longitudonalRecord <- list()
    longitudonalRecordNames <- list()
    longitudonalRecord[1] <- lastMeasure
    longitudonalRecordNames[1] <- paste0(paste0(dataName,'_'), 0)
    
    # look at a specific patient
    if(nrow(serialDf) > 0){
      patientSerialDf <- serialDf %>% filter(PATIENT == myPatient$PATIENT)
    } else {
      patientSerialDf <- serialDf
    }
    
    # how many containers do we need for this individual
    patientContainers <- ceiling(as.numeric(studyEndDate - myBaselineDate)/blockDays)
    if (patientContainers <= 0) {
      stop("Study end date problem")
    }
    
    # for each patient the containers start at baseline
    containerStart <- myBaselineDate
    
    # fill each container
    for (j in 1:patientContainers) {
      containerEnd <- containerStart + blockDays
      if(nrow(patientSerialDf) > 0){
        thisBlock <- patientSerialDf %>% 
          filter(containerStart <= VISIT_END & containerEnd >= VISIT_START)
        
        # is there any data in the block
        if (nrow(thisBlock) > 0) {
          thisBlock <- thisBlock %>% 
            mutate(VISIT_END = safe.ifelse(VISIT_END > containerEnd,containerEnd,VISIT_END))
          thisBlock <- thisBlock%>% 
            mutate(VISIT_START = safe.ifelse(VISIT_START < containerStart,containerStart,VISIT_START))
          thisBlock <- thisBlock %>%
            rowwise() %>%
            mutate(hospitalTime = 
                     safe.ifelse(countNights == TRUE,
                                 as.numeric(VISIT_END - VISIT_START),
                                 as.numeric(VISIT_END - VISIT_START) + 1))
          
          longitudonalRecord[j + 1] <- sum(thisBlock$hospitalTime)
        } else {
          # no data in the block
          longitudonalRecord[j + 1] <- 0
        }
      } else {
        #patient has never been hospitalised
        longitudonalRecord[j + 1] <- 0  
      }
      
      # prepare variables for next block
      longitudonalRecordNames[j + 1] <- paste0(paste0(dataName,'_'), j)
      containerStart <- containerEnd
    }
    
    
    dbRecord <- longitudonalRecord %>% unlist() %>% t() %>% as.data.frame()
    names(dbRecord) <- longitudonalRecordNames %>% unlist()
    dbRecord <- dbRecord %>% mutate(PATIENT = myPatient$PATIENT)
    allRecords[[i]] <- dbRecord
  }
  resultDf <- allRecords %>% ldply()
  
  #create wide form
  if(!longFormat){
    sortedHeadings <- c("PATIENT", c(0:maxContainers) %>% lapply(function(x,
                                                                          myName) {
      paste0(paste0(myName,'_'), x)
    }, dataName)) %>% unlist()
    resultDf <- resultDf %>% select(sortedHeadings)
  } else { #create long form
    resultDf <- makeLong(resultDf, dataName = dataName, idColumn = idColumn)
  }
  resultDf
  
}
