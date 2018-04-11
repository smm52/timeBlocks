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
safe.ifelse <- function(cond, yes, no) {
    class.y <- class(yes)
    X <- ifelse(cond, yes, no)
    class(X) <- class.y
    return(X)
}

#' Checks whether a column is in a data frame and all its rows are not NA
#'
#' @param df the data frame to check
#' @param columnName name of the column
#' @return TRUE if the  column passes the test
#' @keywords internal
checkStructureColumn <- function(df, columnName) {
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
#' @return NA if there are format problems, otherwise a data frame that has a PATIENT, VISIT and MEASURE column
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkBlockFormat <- function(df, idColumn = "PATIENT", dateColumn = "VISIT", dataColumn = "MEASURE") {
    result <- NULL
    if (class(df %>% extract2(dateColumn)) == "Date") {
        if (checkStructureColumn(df, idColumn) & checkStructureColumn(df, dateColumn) & checkStructureColumn(df, dataColumn)) {
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
#' @return  NA if there are format problems, otherwise a data frame that has a PATIENT andVISIT column
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom magrittr extract2
checkBaselineFormat <- function(df, idColumn = "PATIENT", dateColumn = "VISIT") {
    result <- NULL
    if (class(df %>% extract2(dateColumn)) == "Date") {
        if (checkStructureColumn(df, idColumn) & checkStructureColumn(df, dateColumn)) {
            names(df)[which(names(df) == idColumn)] <- "PATIENT"
            names(df)[which(names(df) == dateColumn)] <- "VISIT"
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
findClosestToBaseline <- function(serialDf, baselineDates, blDayDiff = 90, idColumn = 'PATIENT', dateColumn = 'VISIT', dataColumn = 'MEASURE') {
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
                baselineDates[i, "BL_MEASURE_DATE"] <- format(df[[index, "VISIT"]],'%Y-%m-%d')
            }
        }
    }
    baselineDates$BL_MEASURE_DATE <- as.Date(baselineDates$BL_MEASURE_DATE)
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


