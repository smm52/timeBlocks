#' Calculates the time block data for prescriptions
#'
#' @param serialDf the data frame with the serial data
#' @param baselineInfo dataframe with ID and baseline date
#' @param studyStartDate start Date of the study (has to be of class Date)
#' @param studyEndDate end Date of the study (has to be of class Date)
#' @param atcCode regular expression for the class of medication ATC code (default: C09 which stands for any RAAS)
#' @param blockDays length of a time block in days (default 365.25 days)
#' @param useDoses flag indicating whether presence of prescription events or sum of doses should be used in a block (default FALSE)
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the column with the ATC codes: default is ATC
#' @param doseColumn name of the column with the dosage information: default is cDDD
#' @param dataName set a user specified name for the data column in the output: default is prescription
#' @return time block data and dates in wide format
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_at
#' @importFrom plyr ldply
#' @examples
#' \dontrun{
#' baselinePatients <- read_tsv('../data/ourBaselinePatients.txt', guess_max = 3582) %>%
#'     select(PATIENT, VISIT = VISIT_TL)
#' allFimea <- read_tsv('../data/all_fimea_kela.txt')
#' allFimea <- allFimea %>% mutate(VISIT = ymd(allFimea$ostopv))
#' allFimea <- allFimea %>% mutate(PATIENT = fd)
#' allFimea <- allFimea %>% filter(!is.na(all_ddd)) #check these as they should not exists
#' medResLong <- createMedicationTimeBlocks(allFimea, baselinePatients,
#'    studyStartDate = as.Date('1995-01-01', '%Y-%m-%d'),
#'    studyEndDate = as.Date('2015-12-31', '%Y-%m-%d'),
#'    atcCode = 'C09', blockDays = 365.25, useDoses = FALSE, idColumn = "PATIENT", dateColumn = "VISIT",
#'    atcColumn = "ATC_FIMEA", dddColumn = 'all_ddd', dataName = "anyRAAS", longFormat = TRUE,
#'    noUseAtBaseline = FALSE)
#' medResWide <- createMedicationTimeBlocks(allFimea, baselinePatients,
#'    studyStartDate = as.Date('1995-01-01', '%Y-%m-%d'),
#'    studyEndDate = as.Date('2015-12-31', '%Y-%m-%d'),
#'    atcCode = 'C09', blockDays = 365.25, useDoses = FALSE, idColumn = "PATIENT", dateColumn = "VISIT",
#'    atcColumn = "ATC_FIMEA", dddColumn = 'all_ddd', dataName = "anyRAAS", longFormat = FALSE,
#'    noUseAtBaseline = FALSE)
#' }
createMedicationTimeBlocks <- function(serialDf, baselineInfo, studyStartDate, studyEndDate, atcCode = "C09", blockDays = 365.25, useDoses = FALSE,
                                       idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC", dddColumn = "cDDD", dataName = "prescription", longFormat = TRUE, noUseAtBaseline = FALSE) {
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
  # check baseline data
  baselineInfo <- checkBaselineFormat(baselineInfo, idColumn = idColumn, dateColumn = dateColumn)
  if (is.null(baselineInfo)) {
    stop("incorrect format of baseline input data")
  }
  # check serial prescription data
  serialDf <- checkPrescriptionFormat(serialDf, idColumn = idColumn, dateColumn = dateColumn, atcColumn = atcColumn, dddColumn = dddColumn)
  if (is.null(serialDf)) {
    stop("incorrect format of serial prescription input data")
  }
  # restrict serial data to baseline patients
  serialDf <- serialDf %>% filter(PATIENT %in% baselineInfo$PATIENT)
  if (nrow(serialDf) == 0) {
    stop("Serial data does not contain any data for the patients.")
  }

  # for each patient we need to know when they started or ended taking any medication to calculate follow-up times and to know whether a block
  # will be NA (because no information avaialble) or 0 (no prescription with atcCode but other)
  prescritptionDataAvailable <- serialDf %>% group_by(PATIENT) %>% summarise_at(vars(VISIT), funs(min, max), na.rm = TRUE)

  # restrict serial to ATC code under investigation
  serialDf <- serialDf %>% filter(grepl(atcCode, ATC))
  # if (nrow(serialDf) == 0) { stop('Serial data contains no prescriptions with the specific ATC code expression.') }

  # create baseline information
  baselineInfo <- findPrescriptionstAtBaseline(prescriptions = serialDf, baselineDates = baselineInfo, anyPrescriptions = prescritptionDataAvailable,
                                               blockDays = blockDays, useDoses = useDoses)

  # restrict serial data and remove all prescriptions outside the study period
  serialDf <- serialDf %>% filter(VISIT <= studyEndDate & VISIT >= studyStartDate)
  # check whether there are patients if (nrow(serialDf) == 0) { stop('Serial data does not contain any data for the study period.') }


  # initialise
  allRecords <- list()
  maxContainers <- 0
  maxContainersBeforeBl <- 0

  # go through all individual drug records
  for (i in 1:nrow(baselineInfo)) {
    myPatient <- baselineInfo[i, ]
    myBaselineDate <- myPatient$VISIT

    # how many consecutive blocks immediately after baseline are NA
    emptyBlocksAfterBl <- 0
    longitudonalRecord <- list()
    longitudonalRecordNames <- list()
    finalRecord <- NA
    timeContainers <- 0
    timeContainersBeforeBl <- 0
    lengthFollowUp <- 0

    # put baseline in
    longitudonalRecord[1] <- myPatient$BL_MEASURE
    longitudonalRecordNames[1] <- paste0(paste0(dataName, "_"), 0)

    # did the patient take any medication
    myAllPrescriptions <- prescritptionDataAvailable %>% filter(PATIENT == myPatient$PATIENT)

    if (nrow(myAllPrescriptions) > 0) {
      anyPrescriptionStart <- myAllPrescriptions$min
      anyPrescriptionEnd <- myAllPrescriptions$max
      lastFollowUpTime <- min(studyEndDate, anyPrescriptionEnd)
      firstFollowUpTime <- max(studyStartDate, anyPrescriptionStart)
      lengthStudyRecords <- as.numeric(lastFollowUpTime - firstFollowUpTime)/365.25
      if (lastFollowUpTime > myBaselineDate) {
        timeContainers <- ceiling(as.numeric(lastFollowUpTime - myBaselineDate)/blockDays)
        lengthFollowUp <- as.numeric(lastFollowUpTime - max(myBaselineDate, firstFollowUpTime))/365.25
      }
      if (myBaselineDate > firstFollowUpTime) {
        timeContainersBeforeBl <- ceiling(as.numeric(myBaselineDate - firstFollowUpTime)/blockDays)
      }
      maxContainers <- max(maxContainers, timeContainers)
      maxContainersBeforeBl <- max(maxContainersBeforeBl, timeContainersBeforeBl)
      containerStart <- myBaselineDate

      # did this specific patient have any prescriptions of the chosen atcCode
      myPrescriptions <- serialDf %>% filter(PATIENT == myPatient$PATIENT)

      # go through all timeContianers after baseline
      if (timeContainers > 0) {
        for (j in 1:timeContainers) {
          containerEnd <- containerStart + blockDays
          myContrainerPrescriptions <- myPrescriptions %>% filter(VISIT > containerStart & VISIT <= containerEnd)

          # were their prescriptions for the patient in the timeContainer
          if (nrow(myContrainerPrescriptions) > 0) {
            # for binary containers mark them with 1 otherwise sum up the dose
            if (!useDoses) {
              longitudonalRecord[j + 1] <- 1
            } else {
              longitudonalRecord[j + 1] <- sum(myContrainerPrescriptions$cDDD, na.rm = TRUE)
            }

          } else {
            longitudonalRecord[j + 1] <- safe.ifelse(containerStart <= myAllPrescriptions$max & containerEnd >= myAllPrescriptions$min,
                                                     0, NA)
            if (is.na(longitudonalRecord[j + 1]) & emptyBlocksAfterBl == j - 1) {
              emptyBlocksAfterBl <- j
            }
          }

          # add the container name
          longitudonalRecordNames[j + 1] <- paste0(paste0(dataName, "_"), j)

          # remember what happened at the end of follow-up
          finalRecord <- longitudonalRecord[[j + 1]]

          # initialise next container
          containerStart <- containerEnd
        }
      }

      # go through all timeContianers before baseline
      if (timeContainersBeforeBl > 0) {
        containerEnd <- myBaselineDate
        for (j in 1:timeContainersBeforeBl) {
          containerStart <- containerEnd - blockDays
          index <- timeContainers + j
          myContrainerPrescriptions <- myPrescriptions %>% filter(VISIT > containerStart & VISIT <= containerEnd)

          # were their prescriptions for the patient in the timeContainer
          if (nrow(myContrainerPrescriptions) > 0) {
            # for binary containers mark them with 1 otherwise sum up the dose
            if (!useDoses) {
              longitudonalRecord[index + 1] <- 1
            } else {
              longitudonalRecord[index + 1] <- sum(myContrainerPrescriptions$cDDD, na.rm = TRUE)
            }

          } else {
            longitudonalRecord[index + 1] <- safe.ifelse(containerStart <= myAllPrescriptions$max & containerEnd >= myAllPrescriptions$min,
                                                         0, NA)
          }

          # add the container name
          longitudonalRecordNames[index + 1] <- paste0(paste0(dataName, "_"), as.character(j * -1))

          # initialise next container
          containerEnd <- containerStart
        }
      }
    }

    # create one dataframe with all the data
    dbRecord <- longitudonalRecord %>% unlist() %>% t() %>% as.data.frame()
    names(dbRecord) <- longitudonalRecordNames %>% unlist()

    # add patient information and follow-up time
    dbRecord <- dbRecord %>% mutate(PATIENT = myPatient$PATIENT)
    dbRecord <- dbRecord %>% mutate(timeWithFollowUpRecords = lengthFollowUp)
    dbRecord <- dbRecord %>% mutate(totalTimeWithRecords = lengthStudyRecords)
    dbRecord <- dbRecord %>% mutate(finalRecord = finalRecord)

    # set empty baseline to no use if there is no data before baseline and individuals have actually been followed up
    if (noUseAtBaseline & is.na(myPatient$BL_MEASURE) & dbRecord$timeWithFollowUpRecords > 0 & dbRecord$timeWithFollowUpRecords == dbRecord$totalTimeWithRecords) {
      dbRecord <- noUseAtMissingBaseline(dbRecord, blockLength = emptyBlocksAfterBl, dataName = dataName)
    }

    allRecords[[i]] <- dbRecord
  }
  resultDf <- allRecords %>% ldply()

  # create wide form
  if (!longFormat) {
    # create a sort order for the headings
    if (maxContainers > 0 | maxContainersBeforeBl > 0) {
      sortedHeadings <- c("PATIENT", "totalTimeWithRecords", "timeWithFollowUpRecords", "finalRecord", c((-1 * maxContainersBeforeBl):maxContainers) %>%
                            lapply(function(x, myName) {
                              paste0(paste0(myName, "_"), x)
                            }, dataName)) %>% unlist()
    } else {
      sortedHeadings <- c("PATIENT", "totalTimeWithRecords", "timeWithFollowUpRecords", "finalRecord", paste0(paste0(dataName, "_"), 0))
    }

    # re-order the dataframe according to the headings
    resultDf <- resultDf %>% select(sortedHeadings)
  } else {
    # create long form
    resultDf <- makeLong(resultDf, dataName = dataName, idColumn = idColumn)
  }
  resultDf

}

