#' Type safe if.else implementation
#'
#' Type safe if.else implementation that make sure that the output has the same class as the yes condition part.
#' This is important for dates for examples.
#'
#' @param cond the condition
#' @param yes what happens if yes condition is met
#' @param no what happens if no condition is not met
#' @return type safe value
#' @export
#' @examples
#' myDate1 <- NA
#' myDate2 <- as.Date('2000-01-31', '%Y-%m-%d')
#' safe.ifelse(is.na(myDate1),myDate2,myDate1)
#' class(safe.ifelse(is.na(myDate1),myDate2,myDate1))
#'
#' ifelse(is.na(myDate1),myDate2,myDate1)
#' class(ifelse(is.na(myDate1),myDate2,myDate1))
safe.ifelse <- function(cond, yes, no) {
  class.y <- class(yes)
  X <- ifelse(cond, yes, no)
  class(X) <- class.y
  return(X)
}

#' Checks whether a column is in a data frame and all its rows are not NA
#'
#' @param columnName name of the column
#' @param df the data frame to check
#' @return TRUE if the  column passes the test
#' @keywords internal
checkStructureColumn <- function(columnName,df) {
  myReturn <- FALSE
  if (!is.na(columnName) & is.character(columnName)) {
    if (any(names(df) == columnName)) {
      if (length(which(is.na(df[, columnName]) == TRUE)) == 0) {
        myReturn <- TRUE
      }
    }
  }
  myReturn
}

#' Checks the format of the data set for the time block building
#'
#' We need an ID, data and a date column. the date column has to be of class Date.
#'
#' @param df the data frame to check
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param dataColumn name of the data column: default is MEASURE
#' @return NULL if there are format problems, otherwise a data frame that has a PATIENT, VISIT and MEASURE column
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkBlockFormat <- function(df, idColumn = "PATIENT", dateColumn = "VISIT", dataColumn = "MEASURE") {
  result <- NULL
  if (class(df %>% extract2(dateColumn)) == "Date") {
    if (checkStructureColumn(idColumn, df) & checkStructureColumn(dateColumn, df) & checkStructureColumn(dataColumn, df)) {
      names(df)[which(names(df) == idColumn)] <- "PATIENT"
      names(df)[which(names(df) == dateColumn)] <- "VISIT"
      names(df)[which(names(df) == dataColumn)] <- "MEASURE"
      result <- df
    }
  }
  result
}

#' Checks the format of the data set for the baseline dates
#'
#' We need an ID and a date column. the date column has to be of class Date.
#'
#' @param df the data frame to check
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @return  NULL if there are format problems, otherwise a data frame that has a PATIENT andVISIT column
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkBaselineFormat <- function(df, idColumn = "PATIENT", dateColumn = "VISIT") {
  result <- NULL
  if (class(df %>% extract2(dateColumn)) == "Date") {
    if (checkStructureColumn(idColumn, df) & checkStructureColumn(dateColumn, df)) {
      names(df)[which(names(df) == idColumn)] <- "PATIENT"
      names(df)[which(names(df) == dateColumn)] <- "VISIT"
      result <- df
    }
  }
  result
}

#' Checks the format of the hospital data
#'
#' We need an ID and two date columns. The date columns have to be of class Date.
#'
#' @param df the data frame to check
#' @param idColumn name of ID column: default is PATIENT
#' @param startDateColumn name of date column: default is VISIT_START. This column has to be of class Date
#' @param endDateColumn name of date column: default is VISIT_END. This column has to be of class Date
#' @return  NULL if there are format problems, otherwise a data frame that has a PATIENT andVISIT column
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkBaselineHospitalFormat <- function(df, idColumn = "PATIENT", startDateColumn = "VISIT_START", endDateColumn = "VISIT_END") {
  result <- NULL
  if (class(df %>% extract2(startDateColumn)) == "Date" & class(df %>% extract2(endDateColumn)) == "Date") {
    if (checkStructureColumn(idColumn, df) & checkStructureColumn(startDateColumn, df) & checkStructureColumn(endDateColumn, df)) {
      names(df)[which(names(df) == idColumn)] <- "PATIENT"
      names(df)[which(names(df) == startDateColumn)] <- "VISIT_START"
      names(df)[which(names(df) == endDateColumn)] <- "VISIT_END"
      result <- df
    }
  }
  result
}

#' Checks the format of the data set for the prescriptions
#'
#' We need an ID, date column, ATC and DDD column. The date column has to be of class Date.
#'
#' @param df the data frame to check
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the ATC column: default is ATC
#' @param dddColumn name of the DDD column: default is DDD
#' @return NULL if there are format problems, otherwise a data frame that has a PATIENT, VISIT, ATC and DDD column
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkPrescriptionFormat <- function(df, idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC", dddColumn = "cDDD") {
  result <- NULL
  if (class(df %>% extract2(dateColumn)) == "Date") {
    if (checkStructureColumn(idColumn, df) & checkStructureColumn(dateColumn, df) & checkStructureColumn(atcColumn, df) & checkStructureColumn(dddColumn, df)) {
      names(df)[which(names(df) == idColumn)] <- "PATIENT"
      names(df)[which(names(df) == dateColumn)] <- "VISIT"
      names(df)[which(names(df) == atcColumn)] <- "ATC"
      names(df)[which(names(df) == dddColumn)] <- "cDDD"
      result <- df
    }
  }
  result
}

#' Checks the format of the data set for the prescriptions (without DDDs)
#'
#' We need an ID, date column and ATC. The date column has to be of class Date.
#'
#' @param df the data frame to check
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the ATC column: default is ATC
#' @return NULL if there are format problems, otherwise a data frame that has a PATIENT, VISIT, ATC 
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkBinaryPrescriptionFormat <- function(df, idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC") {
  result <- NULL
  if (class(df %>% extract2(dateColumn)) == "Date") {
    if (checkStructureColumn(idColumn, df) & checkStructureColumn(dateColumn, df) & checkStructureColumn(atcColumn, df)) {
      names(df)[which(names(df) == idColumn)] <- "PATIENT"
      names(df)[which(names(df) == dateColumn)] <- "VISIT"
      names(df)[which(names(df) == atcColumn)] <- "ATC"
      result <- df
    }
  }
  result
}



#' Finds the measurement date that is closest to baseline
#'
#' Finds the measurement date that is closest to baseline in the given time period around baseline.
#'
#' @param df the data frame to check
#' @param dayDiff time around baseline (default: 90 days, will look three months before and after baseline)
#' @inheritParams checkBlockFormat
#' @return data frame with baseline measure or NA if no measurement was close to baseline
#' @export
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
#' @examples
#' \dontrun{
#' baselinePatients <- read_tsv('../data/ourBaselinePatients.txt', guess_max = 3582) %>%
#'     select(PATIENT, VISIT = VISIT_TL)
#' bmi <- read_tsv('../data/serial_bmi.txt') # contains serial and baseline data
#' baselinePatients <- findClosestToBaseline(bmi, baselinePatients, blDayDiff = 90)
#' }
findClosestToBaseline <- function(serialDf, baselineDates, blDayDiff = 90, idColumn = "PATIENT", dateColumn = "VISIT", dataColumn = "MEASURE") {
  serialDf <- checkBlockFormat(serialDf, idColumn = idColumn, dateColumn = dateColumn, dataColumn = dataColumn)
  if (is.null(serialDf)) {
    stop("incorrect format of serial input data")
  }
  baselineDates <- checkBaselineFormat(baselineDates, idColumn = idColumn, dateColumn = dateColumn)
  if (is.null(baselineDates)) {
    stop("incorrect format of baseline input data")
  }
  baselineDates <- baselineDates %>% mutate(BL_MEASURE = NA, BL_MEASURE_DATE = NA)
  for (i in 1:nrow(baselineDates)) {
    df <- serialDf[which(serialDf$PATIENT == baselineDates[[i, "PATIENT"]]), ]
    if (nrow(df) > 0) {
      index <- which.min(abs(as.numeric(baselineDates[[i, "VISIT"]] - df$VISIT)))
      if (abs(as.numeric(baselineDates[[i, "VISIT"]] - df[[index, "VISIT"]])) <= blDayDiff) {
        baselineDates[i, "MEASURE"] <- df[[index, "MEASURE"]]
        baselineDates[i, "BL_MEASURE_DATE"] <- format(df[[index, "VISIT"]], "%Y-%m-%d")
      }
    }
  }
  baselineDates$BL_MEASURE_DATE <- as.Date(baselineDates$BL_MEASURE_DATE)
  baselineDates
}

#' Finds the length of a hospital stay at the baseline block
#'
#' Finds the length of hospital stays around baseline
#'
#' @param df the data frame to check
#' @param dayDiff time around baseline (default: 90 days, will look three months before and after baseline)
#' @inheritParams checkBaselineHospitalFormat
#' @param countNights whether to count nights or days in hopsital (if night, than an outpatient stay has length 0, otherwise 1)
#' @return data frame with baseline measure or NA if no measurement was close to baseline
#' @export
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
#' @importFrom dplyr rowwise
#' @importFrom dplyr mutate
findBaselineHospitalStay <- function(serialDf, baselineDates = baselineInfo, 
                                     blDayDiff = blDayDiff, idColumn = "PATIENT", 
                                     startDateColumn = "VISIT_START",  
                                     endDateColumn = "VISIT_END", blDateColumn = "VISIT",
                                     countNights = TRUE){
  serialDf <- serialDf %>% 
    checkBaselineHospitalFormat(idColumn = idColumn, startDateColumn = startDateColumn, 
                                endDateColumn = endDateColumn)
  if (is.null(serialDf)) {
    stop("incorrect format of serial input data")
  }
  baselineDates <- baselineDates %>%
    checkBaselineFormat(idColumn = idColumn, dateColumn = blDateColumn)
  if (is.null(baselineDates)) {
    stop("incorrect format of baseline input data")
  }
  if(countNights == TRUE){
    baselineDates <- baselineDates %>% mutate(BL_HOSPITAL_NIGHTS = 0)
  } else {
    baselineDates <- baselineDates %>% mutate(BL_HOSPITAL_DAYS = 0)  
  }
  if(nrow(serialDf) > 0){
    for (i in 1:nrow(baselineDates)) {
      df <- serialDf[which(serialDf$PATIENT == baselineDates[[i, "PATIENT"]]), ]
      startDate <- as.Date(baselineDates[[i, "VISIT"]]) - blDayDiff
      endDate <- as.Date(baselineDates[[i, "VISIT"]]) + blDayDiff
      df <- df %>% filter(startDate <= VISIT_END & endDate >= VISIT_START)
      if (nrow(df) > 0) {
        df <- df %>% 
          mutate(VISIT_END = safe.ifelse(VISIT_END > endDate,endDate,VISIT_END))
        df <- df %>% 
          mutate(VISIT_START = safe.ifelse(VISIT_START < startDate,startDate,VISIT_START))
        df <- df %>%
          rowwise() %>%
          mutate(hospitalTime = 
                   safe.ifelse(countNights == TRUE,
                               as.numeric(VISIT_END - VISIT_START),
                               as.numeric(VISIT_END - VISIT_START) + 1))
        if(countNights == TRUE){
          baselineDates[i,'BL_HOSPITAL_NIGHTS'] <- sum(df$hospitalTime)
        } else {
          baselineDates[i,'BL_HOSPITAL_DAYS'] <- sum(df$hospitalTime)
        }
      }
    }
  }
  baselineDates
}

#' Finds the prescription information around the baseline visit
#'
#' Finds the prescription information around the baseline visit by looking half a block on either side of baseline.
#' The information is either returned in binary form or in cumulative DDDs depending on the useDoses flag
#'
#' @param serialDf the data frame with all the prescription data
#' @param baselineDates the dataframe with all patient IDs and baseline visit dates
#' @param blockDays length of a block. Baseline info will be searched for the block of length blockDays where the baseline visit sits right in the middle
#' @param usesDoses flag indicating whether to use the presence of a prescription or the cumulate DDDs
#' @return data frame with baseline prescription information
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
findPrescriptionstAtBaseline <- function(prescriptions, baselineDates, anyPrescriptions, blockDays, useDoses) {
  # set the period to look at either side of baseline for prescriptions to half a block. This way the DDDs (if used) will be comparable to a
  # block's DDDs
  blDiff <- blockDays/2
  baselineDates <- baselineDates %>% mutate(BL_MEASURE = 0)
  for (i in 1:nrow(baselineDates)) {
    myRecord <- baselineDates[i, ]
    myPatient <- myRecord$PATIENT
    myBaselineDate <- myRecord$VISIT
    df <- prescriptions %>% filter(PATIENT == myPatient)
    anyPres <- anyPrescriptions %>% filter(PATIENT == myPatient)
    df <- df %>% filter(VISIT >= myBaselineDate - blDiff & VISIT <= myBaselineDate + blDiff)
    anyPres <- anyPres %>% filter(myBaselineDate - blDiff <= max & myBaselineDate + blDiff >= min)
    if (nrow(df) > 0) {
      if (!useDoses) {
        baselineDates[i, "BL_MEASURE"] <- 1
      } else {
        baselineDates[i, "BL_MEASURE"] <- sum(df$cDDD, na.rm = TRUE)
      }
    } else if (nrow(anyPres) == 0) {
      baselineDates[i, "BL_MEASURE"] <- NA
    }
  }
  baselineDates
}

#' calculate the linear interpolation from one measurement to the next available one
#'
#' @param toInterpolate index position of the start of the interpolation
#' @param dataDF the data frame with the measurements records
#' @param dateDF corresponding data frame with date records
#' @param timeframe length of blocks in days
#' @param baselineDate the baseline date
#' @return data frame with (linearly) interpolated values
#' @keywords internal
calculateInterpolation <- function(toInterpolate, dataDF, dateDF, timeframe, baselineDate) {
  myStart <- as.Date(dateDF[[toInterpolate[1] - 1]])
  myEnd <- as.Date(dateDF[[toInterpolate[length(toInterpolate)] + 1]])
  startMeasure <- dataDF[toInterpolate[1] - 1]
  endMeasure <- dataDF[toInterpolate[length(toInterpolate)] + 1]
  dateDiff <- myEnd - myStart
  myDailyFactor <- (endMeasure - startMeasure)/as.numeric(dateDiff)
  for (i in 1:length(toInterpolate)) {
    midDate <- baselineDate + ((toInterpolate[i] - 2) * timeframe) + (0.5 * timeframe)
    myFactor <- as.numeric(midDate - myStart) * myDailyFactor
    dataDF[toInterpolate[i]] <- startMeasure + myFactor
  }
  dataDF
}

#' calculate the linear interpolation for each consecutive set of empty blocks
#'
#' @param dataDF the data frame with the measurements records
#' @param dateDF corresponding data frame with date records
#' @param timeframe length of blocks in days
#' @param baselineDate date of the baseline visit
#' @return data frame with (linearly) interpolated values
#' @keywords internal
doInterpolate <- function(dataDF, dateDF, timeframe, baselineDate) {
  # if baseline date is empty, use the first available date
  if (is.na(baselineDate)) {
    baselineDate <- apply(dateDF, 1, min, na.rm = TRUE)
  }
  # consecutive set of empty blocks
  toInterpolate <- which(apply(dataDF, 2, function(r) any(r == -1)))
  inARow <- split(toInterpolate, cumsum(c(TRUE, diff(toInterpolate) != 1)))
  for (i in 1:length(inARow)) {
    myIndex <- inARow[[i]]
    dataDF <- calculateInterpolation(myIndex, dataDF, dateDF, timeframe, baselineDate)
  }
  dataDF
}

#' Transform the time block results from wide to long format
#'
#' @param df the resulting wide format time block data
#' @param dataName the name of the data column in the output
#' @param idColumn the name of the ID column
#' @return data frame in long format
#' @keywords internal
#' @importFrom tidyr separate
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom dplyr select
#' @importFrom dplyr mutate_at
#' @importFrom dplyr starts_with
#' @importFrom dplyr arrange
#' @importFrom lubridate ymd
#' @importFrom magrittr %>%
makeLong <- function(df, dataName, idColumn) {
  df1 <- df %>% select(c(starts_with(dataName), starts_with(idColumn)))
  df1 <- gather(df1, key, value, -PATIENT) %>% separate(key, sep = "_", into = c(dataName, "block")) %>% spread(dataName, value)
  
  # check whether there are date columns
  if (ncol(df %>% select(starts_with("date"))) > 0) {
    df2 <- df %>% select(c(starts_with("date"), starts_with(idColumn))) %>% mutate_at(vars(starts_with("date")), as.character)
    df2 <- gather(df2, key, value, -PATIENT) %>% separate(key, sep = "_", into = c("date", "block")) %>% spread(date, value) %>% mutate_at(vars(paste0('date.',dataName)), funs(ymd))
    df <- merge(df1, df2, by = c("PATIENT", "block"), all = TRUE)
  } else {
    df <- df1
  }
  df$block <- as.integer(df$block)
  df <- df %>% arrange(PATIENT, block)
  df
}

#' Set missing medication info at baseline and immediately after to 0 which means no use
#'
#' Use this function carefully as it assumes that missing data means no prescritpions
#'
#' @param df wide format result data frame
#' @param blockLength length of consecutive blocks that are empty after baseline
#' @param dataName name of the data columns
#' @return data frame with no medication use added at baseline and immediately after if missing
#' @keywords internal
noUseAtMissingBaseline <- function(df, blockLength, dataName) {
  for (i in 0:blockLength) {
    if (is.na(df[1, paste0(paste0(dataName, "_"), i)])) {
      df[1, paste0(paste0(dataName, "_"), i)] <- 0
    }
  }
  df
}

#' Checks whether a column is in a data frame
#'
#' @param columnName name of the column
#' @param df the data frame to check
#' @return TRUE if the  column passes the test
#' @keywords internal
checkStructureColumnResponse <- function(columnName,df) {
  myReturn <- FALSE
  if (!is.na(columnName) & is.character(columnName)) {
    if (any(names(df) == columnName)) {
      myReturn <- TRUE
    }
  }
  myReturn
}

#' Checks the format whether the format for the data set is a time block based
#'
#'
#' @param df the data frame to check
#' @param idColumn name of ID column
#' @param blockColumn name of time block column
#' @param sexColumn name of the sex column (optional)
#' @param responseColumn name of the column containing the patient block value for a biomarker
#' @return NULL if there are format problems
#' @keywords internal
#' @importFrom magrittr %>%
checkResponseFormat <- function(df, idColumn, blockColumn, sexColumn, responseColumns) {
  result <- NULL
  if (length(responseColumns) > 0) {
    if (checkStructureColumn(idColumn, df) && checkStructureColumn(blockColumn, df) && (is.na(sexColumn) || (checkStructureColumn(sexColumn, df) && length(unique(df[,sexColumn])) == 2))){
      if(all(responseColumns %>% lapply(checkStructureColumnResponse,df) %>% unlist())){
        names(df)[which(names(df) == idColumn)] <- "PATIENT"
        names(df)[which(names(df) == blockColumn)] <- "BLOCK"
        if(!is.na(sexColumn)){
          names(df)[which(names(df) == sexColumn)] <- "SEX"
        }
        result <- df
      }
    }
  }
  result
}

#' Plot the score of an individual patient for all blocks (with binary medication data)
#'
#'
#' @param scoreData data frame containing at least the PATIENT, BLOCK, phenoScore and binaryRAAS
#' @param patient patient number to plot
#' @return a ggplot element that can be plotted with the plot function
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_grid
plotIndividualScore <- function(scoreData, patient){
  scoreData <-  scoreData %>% filter(!(is.na(phenoScore)))
  p <- ggplot(scoreData %>% filter(PATIENT == patient), aes(as.factor(BLOCK), phenoScore, colour = as.factor(anyRAASBinary))) + 
    geom_line(aes(group = PATIENT)) + 
    facet_grid(PATIENT ~ .) +
    scale_colour_discrete(name  ="RAAS",
                          breaks=c("0", "1"),
                          labels=c("No", "Yes"), 
                          drop = FALSE)
  p
}
