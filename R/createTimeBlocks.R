#' Calculates the time block data
#'
#' @param serialDf the data frame with the serial data
#' @param baselineInfo dataframe with ID and baseline date
#' @param studyStartDate start Date of the study (has to be of class Date)
#' @param studyEndDate end Date of the study (has to be of class Date)
#' @param blockDays length of a time block in days (default 365.25 days)
#' @param interpolate shall data for empty blocks be interpolated (default TRUE)
#' @param measureIsFactor when FALSE (default) data will be interpolated linearly, otherwise, empty blocks will take the value of the last preceding non-empty block
#' @param blockSummaryStats how to summarise the data in each block. Currently available: min, max, median, IQR
#' @param blDayDiff time in days that should be searched for a measurement around the baseline visit (default 90)
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param dataColumn name of the data column: default is MEASURE
#' @param dataName set a user specified name for the data column in the output: default is measure
#' @param longFormat flag to indicate whether results should be in long or wide format. Wide has some additional columns (default: TRUE)
#' @return time block data and dates in wide format
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom plyr ldply
#' @examples
#' \dontrun{
#' baselinePatients <- read_tsv('../data/ourBaselinePatients.txt', guess_max = 3582) %>%
#'     select(PATIENT, VISIT = VISIT_TL)
#' bmi <- read_tsv('../data/serial_bmi.txt')
#' resBmiWide <- createTimeBlocks(bmi, baselinePatients, studyStartDate = as.Date('1995-01-01', '%Y-%m-%d'),
#'         studyEndDate = as.Date('2015-12-31', '%Y-%m-%d'), blDayDiff = 90, dataName = "bmi", longFormat = FALSE)
#' resBmiTRUE <- createTimeBlocks(bmi, baselinePatients, studyStartDate = as.Date('1995-01-01', '%Y-%m-%d'),
#'         studyEndDate = as.Date('2015-12-31', '%Y-%m-%d'), blDayDiff = 90, dataName = "bmi", longFormat = TRUE)
#' }
createTimeBlocks <- function(serialDf, baselineInfo, studyStartDate, studyEndDate,
                             blockDays = 365.25, interpolate = TRUE, measureIsFactor = FALSE,
                             blockSummaryStats = "median", blDayDiff = 90,
                             idColumn = "PATIENT", dateColumn = "VISIT", dataColumn = "MEASURE",
                             dataName = "measure", longFormat = TRUE) {
    # did we implement the blockSummaryStats already
    if (!(blockSummaryStats %in% c("median", "min", "max", "IQR","mean","sd","cv","geom_CV","quartile_cd"))) {
        stop(paste0(blockSummaryStats, ": has not been implemented as a block summary statistics yet"))
    }
    #currently an underscore is a special character
    if(grepl('_',dataName)){
      stop('The current implementation does not support underscores (_) in data names')
    }
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
    # no interpolation with 'iqr' stats
    if (blockSummaryStats == "IQR" & interpolate == TRUE) {
        stop("Interpolation for IQR not implemented")
    }

    # restrict serial data to baseline patients
    serialDf <- serialDf %>% filter(PATIENT %in% baselineInfo$PATIENT)
    # check whether there are patients
    if (nrow(serialDf) == 0) {
        stop("Serial data does not contain any of the baseline patients.")
    }

    # check format and find baseline measure and date in serial data
    baselineInfo <- findClosestToBaseline(serialDf, baselineDates = baselineInfo, blDayDiff = blDayDiff, idColumn = idColumn,
        dateColumn = dateColumn, dataColumn = dataColumn)
    # restrict serial data to study period
    serialDf <- serialDf %>% filter(VISIT <= studyEndDate & VISIT >= studyStartDate)
    # check whether there are patients
    if (nrow(serialDf) == 0) {
        stop("Serial data does not contain any data for the study period.")
    }


    allRecords <- list()
    maxContainers <- 0

    # go through all baseline patients
    for (i in 1:nrow(baselineInfo)) {
        myPatient <- baselineInfo[i, ]

        #date of the baseline study visit (might not have a measure)
        myBaselineDate <- myPatient$VISIT

        #initialise lastMeasure with baseline measure
        lastMeasure <- myPatient$MEASURE
        if (is.na(lastMeasure)) {
            minContainerStart <- NA
        } else {
            minContainerStart <- myBaselineDate
        }
        # we cannot calculate a few statistics based on one measurement
        if (blockSummaryStats %in% c("IQR","sd","cv","geom_CV","quartile_cd")) {
            lastMeasure <- NA
        }

        # initial lists and set baseline (measure0) date and value
        longitudonalRecord <- list()
        longitudonalRecordNames <- list()
        longitudonalRecordDates <- list()
        longitudonalRecord[1] <- lastMeasure
        longitudonalRecordNames[1] <- paste0(paste0(dataName,'_'), 0)
        longitudonalRecordDates[1] <- as.character(myBaselineDate)
        timeWithRecords <- 0

        # look at a specific patient
        patientSerialDf <- serialDf %>% filter(PATIENT == myPatient$PATIENT)

        # do we have serial data for this patient?
        if (nrow(patientSerialDf) > 0) {

            # last known measurement
            lastFollowUpTime <- patientSerialDf[[which.max(patientSerialDf$VISIT), "VISIT"]]

            # how many containers do we need for this individual
            patientContainers <- ceiling(as.numeric(lastFollowUpTime - myBaselineDate)/blockDays)

            # record the current maximum of containers needed
            maxContainers <- max(maxContainers, patientContainers)

            # for each patient the containers start at baseline
            containerStart <- myBaselineDate
            finalRecord <- NA

            # fill each container if there are any
            if (patientContainers > 0) {
                for (j in 1:patientContainers) {
                  thisMeasure <- lastMeasure
                  containerEnd <- containerStart + blockDays
                  thisBlock <- patientSerialDf %>% filter(VISIT > containerStart & VISIT <= containerEnd)

                  # is there any data in the block
                  if (nrow(thisBlock) > 0) {
                    minContainerStart <- min(minContainerStart, containerStart, na.rm = TRUE)

                    # calculate the block summary stats
                    if (blockSummaryStats == "max") {
                      indexMax <- which.max(thisBlock$MEASURE)
                      thisMeasure <- thisBlock[indexMax, "MEASURE"]
                    } else if (blockSummaryStats == "min") {
                      indexMax <- which.min(thisBlock$MEASURE)
                      thisMeasure <- thisBlock[indexMax, "MEASURE"]
                    } else if (blockSummaryStats == "median") {
                      thisMeasure <- median(thisBlock$MEASURE, na.rm = TRUE)
                    } else if (blockSummaryStats == "iqr") {
                      thisMeasure <- IQR(thisBlock$MEASURE, na.rm = TRUE)
                    } else if (blockSummaryStats == "mean") {
                      thisMeasure <- mean(thisBlock$MEASURE, na.rm = TRUE)
                    } else if (blockSummaryStats == "sd") {
                      thisMeasure <- sd(thisBlock$MEASURE, na.rm = TRUE)
                    } else if (blockSummaryStats == "cv") {
                      thisMeasure <- sd(thisBlock$MEASURE, na.rm = TRUE) / abs(mean(thisBlock$MEASURE, na.rm = TRUE))
                    } else if (blockSummaryStats == "geom_cv") {
                      thisMeasure <- sqrt(exp((sd(logb(thisBlock$MEASURE), na.rm = TRUE))^2) - 1)
                    } else if (blockSummaryStats == "quartile_cd") {
                      thisMeasure <- IQR(thisBlock$MEASURE, na.rm = TRUE)/(quantile(thisBlock$MEASURE, na.rm = TRUE)[[4]] + quantile(thisBlock$MEASURE, na.rm = TRUE)[[2]])
                    }
                    longitudonalRecord[j + 1] <- thisMeasure
                    finalRecord <- thisMeasure

                    # calculate the block date
                    if (blockSummaryStats %in% c("min", "max")) {
                      thisBlock <- thisBlock %>% filter(MEASURE == thisMeasure)
                      indexMin <- which.min(thisBlock$VISIT)
                    } else {
                      # date with measure closest to stats
                      measureCloseTo <- median(thisBlock$MEASURE, na.rm = TRUE)
                      indexMin <- which(abs(thisBlock$MEASURE - measureCloseTo) == min(abs(thisBlock$MEASURE - measureCloseTo),
                        na.rm = TRUE))
                      indexMin <- indexMin[[1]]
                    }
                    longitudonalRecordDates[j + 1] <- format(thisBlock[[indexMin, "VISIT"]],'%Y-%m-%d')
                  } else {
                    # no data in the block
                    if (is.na(minContainerStart)) {
                      longitudonalRecord[j + 1] <- NA
                    } else {
                      if (interpolate == TRUE & measureIsFactor == FALSE & blockSummaryStats != "iqr") {
                        longitudonalRecord[j + 1] <- -1
                      } else {
                        if (blockSummaryStats == "iqr" | interpolate == FALSE) {
                          longitudonalRecord[j + 1] <- NA
                        } else {
                          longitudonalRecord[j + 1] <- lastMeasure
                        }
                      }
                    }
                    longitudonalRecordDates[j + 1] <- NA
                  }

                  # prepare variables for next block
                  longitudonalRecordNames[j + 1] <- paste0(paste0(dataName,'_'), j)
                  lastMeasure <- thisMeasure
                  containerStart <- containerEnd
                }
                timeWithRecords <- as.numeric(lastFollowUpTime - minContainerStart)/365.25
            }
        }

        dbRecord <- longitudonalRecord %>% unlist() %>% t() %>% as.data.frame()
        names(dbRecord) <- longitudonalRecordNames %>% unlist()
        dateRecord <- longitudonalRecordDates %>% unlist() %>% t() %>% as.data.frame()
        names(dateRecord) <- longitudonalRecordNames %>% lapply(function(x) {
            paste0("date.", x)
        }) %>% unlist()

        # check whether we have to a linear interpolation
        if (interpolate == TRUE & measureIsFactor == FALSE & -1 %in% dbRecord[1, ]) {
            dbRecord <- doInterpolate(dbRecord, dateRecord, blockDays, myBaselineDate)
        }
        if (timeWithRecords > 0) {
            if (blockSummaryStats == "min") {
                dbRecord <- dbRecord %>% mutate(lowestRecord = min(unlist(dbRecord[1, ]), na.rm = TRUE))
            } else if (blockSummaryStats == "max") {
                dbRecord <- dbRecord %>% mutate(highestRecord = max(unlist(dbRecord[1, ]), na.rm = TRUE))
            } else if (blockSummaryStats == "median") {
                dbRecord <- dbRecord %>% mutate(medianRecord = median(unlist(dbRecord[1, ]), na.rm = TRUE))
            } else if (blockSummaryStats %in%  c("IQR","sd","cv","geom_cv","quartile_cd")) {
                dbRecord <- dbRecord %>% mutate(iqrRecord = median(unlist(dbRecord[1, ]), na.rm = TRUE))
            } else if (blockSummaryStats == "mean") {
              dbRecord <- dbRecord %>% mutate(meanRecord = mean(unlist(dbRecord[1, ]), na.rm = TRUE))
            }
        }
        dbRecord <- dbRecord %>% mutate(finalRecord = finalRecord)
        dbRecord <- dbRecord %>% mutate(timeWithRecords = timeWithRecords)
        dbRecord <- cbind(dbRecord, dateRecord) %>% as.data.frame()
        dbRecord <- dbRecord %>% mutate(PATIENT = myPatient$PATIENT)
        allRecords[[i]] <- dbRecord
    }
    resultDf <- allRecords %>% ldply()

    #create wide form
    if(!longFormat){
      if (maxContainers > 0) {
          extremeValue <- c()
          if (blockSummaryStats == "min") {
              extremeValue <- "lowestRecord"
          } else if (blockSummaryStats == "max") {
              extremeValue <- "highestRecord"
          } else if (blockSummaryStats == "median") {
              extremeValue <- "medianRecord"
         } else if (blockSummaryStats == "IQR") {
              extremeValue <- "iqrRecord"
         } else if (blockSummaryStats == "mean") {
              extremeValue <- "meanRecord"
         } else if (blockSummaryStats == "cv") {
              extremeValue <- "coeffVariationRecord"
         } else if (blockSummaryStats == "sd") {
              extremeValue <- "sdRecord"
         } else if (blockSummaryStats == "geom_cv") {
           extremeValue <- "geomCoeffVariationRecord"
         } else if (blockSummaryStats == "quartile_cd") {
           extremeValue <- "quartileCoeffDispersionRecord"
         }
          sortedHeadings <- c("PATIENT", "timeWithRecords", extremeValue[1], "finalRecord", c(0:maxContainers) %>% lapply(function(x,
            myName) {
            paste0(paste0(myName,'_'), x)
          }, dataName), c(0:maxContainers) %>% lapply(function(x, myName) {
            paste0("date.", paste0(paste0(myName,'_'), x))
          }, dataName)) %>% unlist()
      } else {
          sortedHeadings <- c("PATIENT", "timeWithRecords", paste0("date", paste0(".", paste0(dataName, 0))))
      }
      resultDf <- resultDf %>% select(sortedHeadings)
    } else { #create long form
      resultDf <- makeLong(resultDf, dataName = dataName, idColumn = idColumn)
    }
    resultDf
}

