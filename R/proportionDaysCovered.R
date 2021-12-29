#' Calculates adherence to drug treatments based on the proportion of days covered method. 
#'
#' @param serialDf the data frame with prescription data
#' @param startDates a data frame containing the start dates of the study for each patient
#' @param endDates a data frame containing the end dates of the study for each patient
#' @param atcCode a vector containing regular expressions, each encoding for one component/drug class of the treatment
#' @param refillPeriod length of a prescription refill period in days (default 90 days)
#' @param idColumn name of ID column: default is PATIENT
#' @param dateColumn name of date column: default is VISIT. This column has to be of class Date
#' @param atcColumn name of the column with the ATC codes: default is ATC
#' @param treatmentBreakDays a vector containing the number of days (one entry for each drug class) after which the treatment is considered discontinued (default: no breaks applied)
#' @param createGraphs flag indicating whether graphs should be produced (default: FALSE)
#' @return adherence rates for the full prescription period and between start and end dates
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom lubridate interval
#' @importFrom tidyr pivot_longer
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_x_date
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 theme
#' @examples
#' \dontrun{
#'  dfStart <- read_tsv('/home/ad/home/s/stefmutt/projects/former/cadGRS/data/bl_all_new.txt') %>%
#'    select(PATIENT, VISIT)
#'  dfEnd <- dfStart %>%
#'    mutate(VISIT = as.Date('2015-12-31'))
#'  
#'  kela <- read_tsv('/home/ad/home/s/stefmutt/projects/former/cadGRS/data/kela_all.txt') %>%
#'    filter(!is.na(ATC)) %>%
#'    select(-all_ddd) %>%
#'    filter(PATIENT %in% dfStart$PATIENT)
#'  
#'  adherences <- pdc_treatment(serialDf = kela, startDates = dfStart, endDates = dfEnd, atcCode = c('^C09', '^C10'), refillPeriod = 90, 
#'                           idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC", createGraphs = T, 
#'                           treatmentBreakDays = c(181,181))
#' }
pdc_treatment <- function(serialDf, startDates, endDates, atcCode = c(), refillPeriod = 90, 
                          idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC", treatmentBreakDays = c(),
                          createGraphs = F) {
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
  
  # check start date data
  if (nrow(startDates) == 0 ) {
    stop("Start date data is empty")
  } else {
    startDate <- checkBaselineFormat(startDates, idColumn = idColumn, dateColumn = dateColumn)
    if (is.null(startDate)) {
      stop("incorrect format of start dates. Check column class. No missing values allowed.")
    }
  }
  
  #check end date data
  if (nrow(endDates) == 0 ) {
    stop("End date data is empty")
  } else {
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
  if(!all(startDates$PATIENT %in% endDates$PATIENT) & !all(endDates$PATIENT %in% startDates$PATIENT)){
    stop("start and end dates must be specified for the same patients")  
  }
  
  #check that there no duplicates
  if(length(which(duplicated(startDates$PATIENT))) > 0){
    stop('Remove duplicated individuals in start dates')
  }
  if(length(which(duplicated(endDates$PATIENT))) > 0){
    stop('Remove duplicated individuals in end dates')
  }
  
  #check that the underscore is not used in patient IDs
  if(any(grepl('_', startDates$PATIENT))){
    stop('Please remove all underscores from patient IDs!')
  }
  
  # restrict serial data to patients in startDates (if available)
  serialDf <- serialDf %>% filter(PATIENT %in% startDates$PATIENT)
  if (nrow(serialDf) == 0) {
    stop("Serial prescription data does not contain any data for the patients.")
  }
  
  # restrict serial data to patients in endDates
  serialDf <- serialDf %>% filter(PATIENT %in% endDates$PATIENT)
  if (nrow(serialDf) == 0) {
    stop("Serial prescription data does not contain any data for the patients.")
  }
  
  # check that list of ATC codes is not emtpy
  if(length(atcCode) == 0){
    stop("List of ATC codes missing.")
  }
  
  # generate one regular expression for all medications
  allCodesRegExp <- paste(atcCode, sep = '', collapse = ')|(')
  allCodesRegExp <- paste0('(', allCodesRegExp, ')')
  
  # restrict serial to ATC codes under investigation
  serialDf <- serialDf %>% 
    filter(grepl(allCodesRegExp, ATC))
  if (nrow(serialDf) == 0) { 
    stop('Serial data contains no prescriptions with the specific ATC codes expression.') 
  }
  
  # if no treatment breaks given set them to 0, otherwise check that list has same length than ATC codes
  if(length(treatmentBreakDays) == 0){
    for(i in 1:length(atcCode)){
      treatmentBreakDays[i] <- NA
    }
  } else if(length(atcCode) != length(treatmentBreakDays)){
    stop("You need to provide treament break days for each medication")
  }
  
  # create medication data for each treatment subclass
  treatmentRegimes <- list()
  firstPrescription <- min(serialDf$VISIT)
  lastPrescription <- max(serialDf$VISIT)
  mySubclass <- 1
  for(i in 1: length(atcCode)){
    subDf <- serialDf %>%
      filter(grepl(atcCode[i],ATC))
    treatmentRegimes[[i]] <- subDf %>%
      mutate(TREATMENT_SUBCLASS = mySubclass) %>%
      mutate(ROW_ID = paste0(PATIENT, '_', TREATMENT_SUBCLASS)) %>%
      mutate(COL_ID = interval(firstPrescription,VISIT)/days(1)) %>%
      dplyr::arrange(ROW_ID, COL_ID) %>%
      select(ROW_ID, COL_ID) %>%
      mutate(COL_ID = paste0(COL_ID))
    # ensure that the first non-empty medication is subclass 1
    if(nrow(subDf) > 0){
      mySubclass <- mySubclass + 1
    }
  }
  # only keep the treatment regimes that have prescriptions
  hasPrescriptions <- sapply(treatmentRegimes,function(x){nrow(x)  != 0})
  treatmentRegimes <- treatmentRegimes[hasPrescriptions]
  atcCode <- atcCode[hasPrescriptions]
  treatmentBreakDays <- treatmentBreakDays[hasPrescriptions]
  
  #add days after first prescription to start dates and end dates
  startDates <- startDates %>%
    mutate(DAYS = interval(firstPrescription,VISIT)/days(1))
  endDates <- endDates %>%
    mutate(DAYS = interval(firstPrescription,VISIT)/days(1))
  
  #create results matrix
  allPatients <- unique(startDates$PATIENT) 
  myRowNames <- paste0(rep(allPatients,each = length(atcCode)), '_', seq(1:length(atcCode)))
  daysBetweenFirstLast <- (interval(firstPrescription, lastPrescription)/days(1)) + refillPeriod
  myColNames <- c(c(0,seq(0:daysBetweenFirstLast)))
  treatmentTable <- matrix(NA, nrow = length(myRowNames), ncol = length(myColNames), dimnames = list(myRowNames, myColNames))
  
  print('----------------------creating individual treatment table')
  # fill in
  print('------------------------------fill in medicated days')
  for(i in 1:length(treatmentRegimes)){
    myData <- treatmentRegimes[[i]] %>%
      as.matrix()
    oldPatient <- ''
    myBreak <- treatmentBreakDays[i]
    treatmentStart <- daysBetweenFirstLast + 2
    treatmentEnd <- 1
    colEnd_corrected <- 1
    oldRowNumber <- ''
    for(j in 1:nrow(myData)){
      rowNumber <- myData[j,'ROW_ID']
      newPatient <-  sub("_.*", "", rowNumber)
      colNumber <- as.numeric(myData[j,'COL_ID']) + 1
      if(newPatient != oldPatient | is.na(treatmentTable[rowNumber,colNumber])){
        offset <- 0
      } else {
        offset <- colEnd_corrected - colNumber + 1
      }
      colNumber <-colNumber +  offset
      if(colNumber < (daysBetweenFirstLast + 2)){
        colEnd <- colNumber - 1 + refillPeriod
        colEnd_corrected <- min(colEnd, (daysBetweenFirstLast + 2))
        alreadyCovered <- sum(treatmentTable[rowNumber, colNumber:colEnd_corrected], na.rm = T)
        if(alreadyCovered == 0){
          offset <- 0
        } else {
          offset <- offset + alreadyCovered
          colEnd_corrected <- min((colEnd_corrected + offset), (daysBetweenFirstLast + 2))
        }
        treatmentTable[rowNumber, colNumber:colEnd_corrected] <- 1
        oldPatient <- sub("_.*", "", rowNumber)
        oldRowNumber <- rowNumber
      } 
    }
  }
  
  #fill in breaks and uncovered days
  print('------------------------------fill in uncovered and break days')
  myRownames <- rownames(treatmentTable)
  treatmentTable <- apply(as.matrix(rownames(treatmentTable)),1,addUncovered,treatmentTable, treatmentBreakDays) %>%
    t()
  rownames(treatmentTable) <- myRownames
  
  #fill in hospital days
  
  #calculate adherence
  print('----------------------calculate adherences')
  myResults <- calculateAdherences(treatmentTable, startDates, endDates, refillPeriod)
  
  #create graphs
  if(createGraphs){
    print('----------------------create graphs')
    allPatientsWithAdherence <- unique((myResults %>% filter(!is.na(adherenceFullTime)))$PATIENT)
    treatmentTable <- treatmentTable %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      mutate(PATIENT = sub("_.*", "", rowname),
             treatmentCode = sub("^.*_", "", rowname)) %>%
      select(-rowname) %>%
      filter(PATIENT %in% allPatientsWithAdherence)
    treatmentTable <- treatmentTable %>%
      pivot_longer(col = !one_of('PATIENT', 'treatmentCode'), names_to = 'days', values_to = 'covered', values_drop_na = F)
    treatmentTable  <- treatmentTable %>% 
      mutate(dateDay = as.Date(firstPrescription) + as.numeric(days))
    l <- 1
    myPlots <- list()
    plotsPerPage <- 20
    pagesNeeded <- ceiling(length(allPatientsWithAdherence)/plotsPerPage)
    pdf('treatmentGraphs.pdf', width=21, height=27)
    for(n in 1:length(allPatientsWithAdherence)){
      patientTreatments <- treatmentTable %>%
        filter(PATIENT == allPatientsWithAdherence[n]) %>% 
        filter(!is.na(covered)) %>%
        mutate(covered = factor(covered, levels = c(0,1), labels = c('no', 'yes')))
      patientStart <- startDates[which(startDates$PATIENT == allPatientsWithAdherence[n]),'VISIT'][[1]]
      patientEnd <- endDates[which(endDates$PATIENT == allPatientsWithAdherence[n]),'VISIT'][[1]]
      myXmin <- min(min(patientTreatments$dateDay), patientStart)
      myXmax <- max(max(patientTreatments$dateDay), patientEnd)
      fa <- NA
      if(length(atcCode) > 1){
        fa <- myResults %>%
          filter(treatmentCode == 'treatment' & PATIENT == allPatientsWithAdherence[n])
      } else {
        fa <- myResults %>%
          filter(treatmentCode == 1 & PATIENT == allPatientsWithAdherence[n])  
      }
      myPlots[[l]] <- ggplot(patientTreatments, aes(dateDay, treatmentCode, colour = covered)) + 
        geom_point(shape = 15) +
        geom_vline(xintercept = patientStart, linetype = 'dotdash') +
        geom_vline(xintercept = patientEnd, linetype = 'dotdash') +
        ggtitle(paste0('Adherences for ', allPatientsWithAdherence[n],'\n full-time adherence ', fa[1,'adherenceFullTime'],
                       '% \n start to end adherence ', fa[1,'adherenceStartToEnd'], '%')) +
        xlab('timeline') +
        scale_y_discrete(name = 'treatment', labels = atcCode) +
        scale_color_manual(values = c("#999999", "#E69F00"), drop = F) +
        scale_x_date(date_minor_breaks = "1 year", date_labels = "%Y", limits = c(myXmin - 10,myXmax + 10)) +
        theme_classic() +
        theme(plot.title = element_text(size=18, face="bold"), axis.text = element_text(size = 12)) 
      if (l %% plotsPerPage == 0) { ## print plotsPerPage plots on a page
        print(paste0('------------------------------printing page ',ceiling(n/plotsPerPage),' of ',pagesNeeded))
        do.call(grid.arrange,  myPlots)
        myPlots = list() # reset plot 
        l = 0 # reset index
      }
      l = l + 1
    }
    if (length(myPlots) != 0) { 
      print(paste0('------------------------------printing page ',ceiling(n/plotsPerPage),' of ',pagesNeeded))
      do.call(grid.arrange,  myPlots)
    }
    dev.off()
  }
  
  myResults
}

#' Adds uncovered and breaks days into the treatment of one patient and one drug class
#'
#' @param treatmentName the row ID that defines the patient and drug class
#' @param treatmentTable the complete treatment table
#' @param myBreakDays the break days for all drug classes
#' @return uncovered and breaks days for one patient and one drug class
addUncovered <- function(treatmentName, treatmentTable, myBreakDays){
  myTreatmentRegime <- as.numeric(sub("^.*_", "", treatmentName))
  myBreak<- myBreakDays[myTreatmentRegime]
  treatmentRow <- treatmentTable[treatmentName,]
  
  #skip unmedicated
  if(length(which(!is.na(treatmentRow))) > 0){
    treatmentStart <- which.min(treatmentRow == 1)
    treatmentEnd <- max(which(treatmentRow == 1))
    allStops <- which(is.na(treatmentRow)) 
    if(length(allStops) > 0){
      allStops <- allStops[allStops > treatmentStart & allStops < treatmentEnd]
      if(length(allStops) > 0){
        if(!is.na(myBreak)){
          individualBreaks <- allStops
          breakIndex <- 1
          if(length(individualBreaks) > 1){
            for(m in 1:(length(individualBreaks) - 1)){
              keep <- individualBreaks[m]
              individualBreaks[m] <- breakIndex
              if((individualBreaks[m+1] - 1) != keep){
                breakIndex <- breakIndex+1
              }
            }
          }
          individualBreaks[length(individualBreaks)] <- breakIndex
          result <- rle(individualBreaks) 
          tooLong <- which(result$lengths >= myBreak)
          if(length(tooLong) > 0){
            endBreak <- names(result$values[tooLong])
            revResult <- rle(rev(individualBreaks))
            tooLongRev <- which(revResult$lengths >= myBreak)
            startBreak <- rev(names(revResult$values[tooLongRev]))
            for(k in 1:length(tooLong)){
              allStops <- allStops[which(!between(allStops,
                                                  as.numeric(startBreak[k]) + 1 + myBreak, 
                                                  as.numeric(endBreak[k]) + 1))]
            }
          }  
        }
        treatmentRow[allStops] <- 0
      }
    }
  }
  treatmentRow
}

#' Controls the adherence calculation for the whole treatment table
#'
#' @param treatmentTable the complete treatment table matrix
#' @param startDates a data frame containing the start dates of the study for each patient
#' @param endDates a data frame containing the end dates of the study for each patient
#' @param refillPeriod length of a prescription refill period in days
#' @importFrom plyr ldply
#' @return adherences for the entire treatment table as a data frame
calculateAdherences <- function(treatmentTable, startDates, endDates, refillPeriod){
  res <- apply(as.matrix(rownames(treatmentTable)), 1, calculateMedicationSpecificAdherence, treatmentTable, startDates, 
               endDates, refillPeriod) %>% ldply()
  res
}

#' Calculate the adherence for each patient and all their drug classes
#'
#' @param treatmentRowNames the row ID from the treatment table identifying one patent and one drug class
#' @param treatment the complete treatment table matrix
#' @param startDates a data frame containing the start dates of the study for each patient
#' @param endDates a data frame containing the end dates of the study for each patient
#' @param refillPeriod length of a prescription refill period in days
#' @return adherences for a patient and drug class
calculateMedicationSpecificAdherence <- function(treatmentRowNames, treatment, startDates, endDates, refillPeriod){
  #result dataframe
  result <- data.frame(PATIENT=character(), treatmentCode = integer(), adherenceFullTime = double(), 
                       treatmentDaysFullTime = integer(), adherenceStartToEnd = double(), treatmentDaysStartEnd = integer(), 
                       stringsAsFactors = F)
  myPatient <- sub("_.*", "", treatmentRowNames)
  myTreatmentRegime <- sub("^.*_", "", treatmentRowNames)
  result[1,'PATIENT'] <- myPatient
  result[1,'treatmentCode'] <- myTreatmentRegime
  treatmentRow <- treatment[treatmentRowNames,]
  myStart <- startDates[which(startDates$PATIENT == sub("_.*", "", treatmentRowNames)),'DAYS'][[1]]
  myEnd <- endDates[which(endDates$PATIENT == sub("_.*", "", treatmentRowNames)),'DAYS'][[1]]
  
  #overall adherence during the full follow-up
  result[1,'treatmentDaysFullTime'] <-  length(which(treatmentRow %in% c(0,1)))
  if(result[1,'treatmentDaysFullTime'] > 0) {
    result[1,'adherenceFullTime'] <- calculateOneAdherenceValue(treatmentRow, refillPeriod)
    
    #adherence between start and date date
    treatmentRow <- treatmentRow[which(as.numeric(names(treatmentRow)) >= myStart & as.numeric(names(treatmentRow)) <= myEnd)]
    result[1,'treatmentDaysStartEnd'] <-  length(which(treatmentRow %in% c(0,1)))
    if(result[1,'treatmentDaysStartEnd'] > 0){
      result[1,'adherenceStartToEnd'] <-  calculateOneAdherenceValue(treatmentRow, refillPeriod)
    }
  } else {
    result[1,'treatmentDaysStartEnd'] <- 0  
  }
  
  if(myTreatmentRegime == '1' & length(which(grepl(myPatient,rownames(treatment)))) > 1){
    result[2,'PATIENT'] <- myPatient
    result[2,'treatmentCode'] <- 'treatment'
    myTreatment <- treatment[which(grepl(myPatient,rownames(treatment))),]
    combinedRow <- apply(myTreatment,2,combineRows)
    result[2,'treatmentDaysFullTime'] <-  length(which(combinedRow %in% c(0,1)))
    if(result[2,'treatmentDaysFullTime'] > 0){
      result[2,'adherenceFullTime'] <- calculateOneAdherenceValue(combinedRow, refillPeriod)
      combinedRow <-combinedRow[which(as.numeric(names(combinedRow)) >= myStart & as.numeric(names(combinedRow)) <= myEnd)]
      result[2,'treatmentDaysStartEnd'] <-  length(which(combinedRow %in% c(0,1)))
      if(result[2,'treatmentDaysStartEnd'] > 0){
        result[2,'adherenceStartToEnd'] <-  calculateOneAdherenceValue(combinedRow, refillPeriod)
      }
    } else {
      result[2,'treatmentDaysStartEnd'] <- 0
    }
  }
  
  result
}

#' Combine the information from all drug classes into one row (logical AND)
#'
#' @param myValue 0, 1 or NA from different drug classes on the same date
#' @return combined row
combineRows <- function(myValues){
  newValue <- NA
  if(any(!is.na(myValues))){
    if(length(which(myValues == 0)) > 0){
      newValue <- 0
    } else {
      newValue <- 1
    }
  }
  newValue
}

#' Calculate one adherence rate
#'
#' @param treatmentRow one row from the treatment table identifying one patient and one drug class
#' @param refillPeriod length of a prescription refill period in days
#' @return adherence
calculateOneAdherenceValue <- function(treatmentRow, refillPeriod){
  treatmentDays <- length(which(treatmentRow %in% c(0,1)))
  #if there were 0 prescriptions (treatmentDays == 0) or just one (treatmentDays == refillPeriod), do not calculate an adherence
  if(treatmentDays <= refillPeriod){
    adherence <- NA
  } else {
    adherence <- round(100 * sum(treatmentRow, na.rm = T) / treatmentDays, digits = 2)
  }
  adherence
}

