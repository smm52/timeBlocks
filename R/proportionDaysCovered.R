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
#' @param absenceDays a data frame containing start dates and end dates of absences for each patient. This time will be removed from the calculation. The first day should be stored in a column called START and the final in one called END. (optional)
#' @param createGraphs flag indicating whether graphs should be produced (default: FALSE)
#' @param savePrescriptionTable flag indicating whether the whole prescription table should be saved in a file (default: FALSE)
#' @return adherence rates for the full prescription period and between start and end dates
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr rowwise
#' @importFrom dplyr bind_rows
#' @importFrom lubridate interval
#' @importFrom lubridate days
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
#'                           treatmentBreakDays = c(181,181), absenceDays = NULL)
#' }
pdc_treatment <- function(serialDf, startDates, endDates, atcCode = c(), refillPeriod = 90, 
                          idColumn = "PATIENT", dateColumn = "VISIT", atcColumn = "ATC", 
                          treatmentBreakDays = c(), absenceDays = NULL, createGraphs = F, savePrescriptionTable = F) {
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
  
  # check absence data
  if(!is.null(absenceDays)){
    if(nrow(absenceDays) == 0){
      stop("absence data is empty.")
    }
    if(length(which(names(absenceDays) == idColumn)) == 0) {
      stop("incorrect ID in absence date")
    }
    names(absenceDays)[which(names(absenceDays) == idColumn)] <- "PATIENT"
    if(!any(c('START', 'END') %in% names(absenceDays))){
      stop("The first day of a absence should be in a column called START and the final in one called END.")
    }
    if(!(class(absenceDays$START) == 'Date' & class(absenceDays$END) == 'Date')){
      stop("The columns START and END have to be of class Date.")
    }
    # remove patients not in startDates from absence
    absenceDays <- absenceDays %>%
      filter(PATIENT %in% startDates$PATIENT)
    if(nrow(absenceDays) == 0){
      stop("absence data is empty for patients in start dates.")
    }
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
      mutate(COL_ID = as.numeric(interval(firstPrescription,VISIT)/days(1))) %>%
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
  
  # remove absence time before first and after last prescription and change the presentation to days after first prescription
  if(!is.null(absenceDays)){
    absenceDays <- absenceDays %>%
      filter(END > firstPrescription) %>%
      filter(START < (lastPrescription + refillPeriod)) %>%
      filter(END > START) %>%
      mutate(START = safe.ifelse(START < firstPrescription, firstPrescription, START)) %>%
      mutate(END = safe.ifelse(END > (lastPrescription + refillPeriod), (lastPrescription + refillPeriod), END)) %>%
      mutate(START_DAY = interval(firstPrescription,START)/days(1)) %>%
      mutate(END_DAY = interval(firstPrescription,END)/days(1)) %>%
      select(-START, -END)
  }
  
  #create results matrix
  allPatients <- unique(startDates$PATIENT) 
  myRowNames <- paste0(rep(allPatients,each = length(atcCode)), '_', seq(1:length(atcCode)))
  daysBetweenFirstLast <- (interval(firstPrescription, lastPrescription)/days(1)) + refillPeriod
  myColNames <- c(0,seq(1:(daysBetweenFirstLast - 1)))
  treatmentTable <- matrix(NA, nrow = length(myRowNames), ncol = length(myColNames), dimnames = list(myRowNames, myColNames))
  
  print('----------------------creating individual treatment table')
  # fill in
  #go through all treatment regimes
  for(i in 1:length(treatmentRegimes)){
    myTreatment <- treatmentRegimes[[i]]
    myBreakDays <- treatmentBreakDays[i]
    tableRowNames <- myRowNames[which(as.numeric(sub("^[^_]*_", "", myRowNames)) == i)]
    #go through all patients with the chosen treatment
    for(j in 1:length(tableRowNames)){
      myID = tableRowNames[j]
      #nothing to do if the patient does not have the medication
      if(myID %in% myTreatment$ROW_ID){
        patientTreatment <- myTreatment %>%
          filter(ROW_ID == myID)
        myStart <- as.numeric(patientTreatment[1,'COL_ID'][[1]])
        for(k in 1:nrow(patientTreatment)){
          currentDay <- as.numeric(patientTreatment[k,'COL_ID'][[1]])
          currentPeriod <- seq(from = currentDay, length.out = refillPeriod)
          #check if some of the days have already been covered
          extraDays <- sum(treatmentTable[myID,paste(currentPeriod)], na.rm = T)
          #put the current prescription period in the treatment table
          treatmentTable[myID, paste(currentPeriod)] <- 1
          myEnd <- max(currentPeriod)
          #resolve extraDays
          while(extraDays > 0){
            newStart <-max(currentPeriod) + 1
            currentPeriod <- seq(from = newStart, length.out = extraDays)
            #if extra days go beyond or reach the end of the prescription period end 
            if(max(currentPeriod) >= max(myColNames)){
              extraDays <- 0
              if(newStart <= max(myColNames)){
                currentPeriod <- seq(from = newStart, to = max(myColNames)) 
                treatmentTable[myID, paste(currentPeriod)] <- 1
                myEnd <- max(currentPeriod)
              }
            } else {
              extraDays <-sum(treatmentTable[myID,paste(currentPeriod)], na.rm = T)   
              treatmentTable[myID, paste(currentPeriod)] <- 1
              myEnd <- max(currentPeriod)
            }
          }
        }
        #fill in unmedicated days
        inBetween <- seq(from = myStart, to = myEnd)
        
        #remove break days from unmedicated
        if(!is.na(myBreakDays)){
          vec <- treatmentTable[myID,paste(inBetween)]
          rl <- rle(is.na(vec))
          i1 <- rl$lengths > myBreakDays & rl$values
          lst <- split(vec, rep(cumsum(c(TRUE, i1[-length(i1)])), rl$lengths)) 
          unmedicated <- lapply(lst,function(x){ if(is.na(tail(x,1))){x <- x[1:((max(which(!is.na(x))) + myBreakDays))]}; x}) %>%
            unlist()
          names(unmedicated) <- sub("^[^.]*.", "", names(unmedicated))
          unmedicated <- which(is.na(unmedicated))
          treatmentTable[myID,paste(names(unmedicated))] <- 0
        } else {
          #when no breaks, set everything between start and end to 0 which is not 1 (e.g. NA)
          unmedicated <- which(is.na(treatmentTable[myID,paste(inBetween)]))  
          treatmentTable[myID,unmedicated] <- 0
        }

        #fill in absence days
        if(!is.null(absenceDays)){
          #check again as potentially all absence dates might have been outside the prescription time
          if(nrow(absenceDays) > 0){
            myDays <- absenceDays %>%
              filter(PATIENT == gsub("_.*","",myID)) %>%
              dplyr::arrange(START_DAY)
            if(nrow(myDays) > 0){
              for(m in 1:nrow(myDays)){
                tryCatch({
                  theseDays <- myDays[m,]
                  if(theseDays$START_DAY <= myEnd & theseDays$END_DAY >= myStart){
                    thisPeriod <- seq(from = theseDays$START_DAY, length.out = theseDays$END_DAY - theseDays$START_DAY + 1)
                    #check if some of the days have already been covered, add more extra days after the absence stay
                    extraDays <- sum(treatmentTable[myID,paste(thisPeriod)], na.rm = T)
                    #check that absence data does not overlap
                    if(length(which(treatmentTable[myID,paste(thisPeriod)] == 2)) > 0){
                      stop(paste('overlapping absence periods for ID ', gsub("_.*","",myID), '. No absences added for this ID.'))
                    }
                    #put the current period in the treatment table
                    treatmentTable[myID, paste(thisPeriod)] <- 2 #set to 2 to differentiate from NA days in the plotting, for adherence calc. needs to be NA again  
                    while(extraDays > 0){
                      newStart <-max(thisPeriod) + 1
                      thisPeriod <- seq(from = newStart, length.out = extraDays)
                      #if extra days go beyond or reach the end of the prescription period end 
                      if(max(thisPeriod) >= max(myColNames)){
                        extraDays <- 0
                        if(newStart <= max(myColNames)){
                          thisPeriod <- seq(from = newStart, to = max(myColNames)) 
                          treatmentTable[myID, paste(thisPeriod)] <- 1
                        }
                      } else {
                        extraDays <-sum(treatmentTable[myID,paste(thisPeriod)], na.rm = T)   
                        treatmentTable[myID, paste(thisPeriod)] <- 1
                      }
                    }
                  }
                }, error = function(e){
                  print(paste("ERROR :",conditionMessage(e)))
                })
              }
            }
          }
        }
      }
    }
  }
  #save treatment table
  if(savePrescriptionTable){
    print('----------------------saving prescription table')
    write.csv(treatmentTable, file = 'prescriptionTable.csv')  
  }
  
  #calculate adherence
  print('----------------------calculate adherences for individual medications')
  myResults <- calculateAdherences(treatmentTable, startDates, endDates, refillPeriod)
  if(length(atcCode) > 1){
    print('----------------------calculate adherences for full treatment')
    #combine rows for each patient ID with treatment code 0
    for(i in 1:length(allPatients)){
      newRow <- apply(treatmentTable[which(grepl(allPatients[i],rownames(treatmentTable))),],2,combineRows)
      treatmentTable <- rbind(treatmentTable, newRow)
      rownames(treatmentTable)[which(rownames(treatmentTable) == 'newRow')] <- paste0(allPatients[i],'_0')
    }
    myResultsTreatment <- calculateAdherences(treatmentTable[which(grepl('_0$',rownames(treatmentTable))),], 
                                              startDates, endDates, refillPeriod)
    myResults <- bind_rows(myResults, myResultsTreatment) %>%
      dplyr::arrange(PATIENT, treatmentCode) %>%
      mutate(treatmentCode = as.numeric(treatmentCode)) %>%
      rowwise() %>%
      mutate(treatment = ifelse(treatmentCode == 0,'polyPharmacy',atcCode[treatmentCode]))
  }
  
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
    myLabels <- c('no', 'yes')
    myLevels <- c(0,1)
    myColours <- c("#999999", "#E69F00")
    if(any(treatmentTable == 2)){
      myLevels <- c(0,1,2)  
      myLabels <- c('no', 'yes', 'absence')
      myColours <- c("#999999", "#E69F00", "#CC79A7")
    }
    pdf('treatmentGraphs.pdf', width=21, height=27)
    for(n in 1:length(allPatientsWithAdherence)){
      patientTreatments <- treatmentTable %>%
        filter(PATIENT == allPatientsWithAdherence[n]) %>% 
        filter(!is.na(covered)) %>%
        mutate(covered = factor(covered, levels = myLevels, labels = myLabels))
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
        geom_point(shape = 15, size = 0.1) +
        geom_vline(xintercept = patientStart, linetype = 'dotdash') +
        geom_vline(xintercept = patientEnd, linetype = 'dotdash') +
        ggtitle(paste0('Adherences for ', allPatientsWithAdherence[n],'\n full-time adherence ', fa[1,'adherenceFullTime'],
                       '% \n start to end adherence ', fa[1,'adherenceStartToEnd'], '%')) +
        xlab('timeline') +
        scale_y_discrete(name = 'treatment', labels = atcCode) +
        scale_color_manual(values = myColours, drop = F) +
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
                       treatmentDaysFullTime = integer(), adherenceStartEnd = double(), treatmentDaysStartEnd = integer(), 
                       stringsAsFactors = F)
  myPatient <- sub("_.*", "", treatmentRowNames)
  myTreatmentRegime <- sub("^.*_", "", treatmentRowNames)
  result[1,'PATIENT'] <- myPatient
  result[1,'treatmentCode'] <- sub(".*_", "", treatmentRowNames)
  treatmentRow <- treatment[treatmentRowNames,]
  #ignore hospital days to calculate the adherence
  treatmentRow[which(treatmentRow == 2)] <- NA
  myStart <- startDates[which(startDates$PATIENT == sub("_.*", "", treatmentRowNames)),'DAYS'][[1]]
  myEnd <- endDates[which(endDates$PATIENT == sub("_.*", "", treatmentRowNames)),'DAYS'][[1]]
  
  #overall adherence during the full follow-up
  result[1,'treatmentDaysFullTime'] <-  length(which(treatmentRow %in% c(0,1)))
  if(result[1,'treatmentDaysFullTime'] > 0) {
    result[1,'adherenceFullTime'] <- calculateOneAdherenceValue(treatmentRow, refillPeriod)
    
    #adherence between start and date date
    treatmentRow <- treatmentRow[paste(seq(from = myStart, to = myEnd))]
    result[1,'treatmentDaysStartEnd'] <-  length(which(treatmentRow %in% c(0,1)))
    if(result[1,'treatmentDaysStartEnd'] > 0){
      result[1,'adherenceStartEnd'] <-  calculateOneAdherenceValue(treatmentRow, refillPeriod)
    }
  } else {
    result[1,'treatmentDaysStartEnd'] <- 0  
  }
  
  result
}

#' Combine the information from all drug classes into one row (logical AND)
#'
#' @param myValue 0, 1, 2 or NA from different drug classes on the same date
#' @return combined row
combineRows <- function(myValues){
  #if the value is 2 from hospitalizations, ignore (e.g. set to NA)
  myValues[which(!myValues %in% c(0,1))] <- NA
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
  #if there were 0 prescriptions (treatmentDays == 0) or just the equivalent number of treatment days to one pres.
  #(treatmentDays == refillPeriod), do not calculate an adherence
  if(treatmentDays <= refillPeriod){
    adherence <- NA
  } else {
    adherence <- round(100 * sum(treatmentRow, na.rm = T) / treatmentDays, digits = 2)
  }
  adherence
}


