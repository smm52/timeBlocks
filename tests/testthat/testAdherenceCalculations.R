meds <- read_tsv('../examples/medications.csv')
starts <- read_tsv('../examples/starts.csv')
ends <- read_tsv('../examples/ends.csv')
abs <- read_tsv('../examples/absences.csv')

test_that("one medication, no breaks, no absence", {
  res <- pdc_treatment(serialDf = meds  %>% filter(PATIENT %in% c(1,2,3,4,5,6)),
                              startDates = starts  %>% filter(PATIENT %in% c(1,2,3,4,5,6)), 
                              endDates = ends  %>% filter(PATIENT %in% c(1,2,3,4,5,6)), 
                              atcCode = c('C09'), 
                              refillPeriod = 7, 
                              treatmentBreakDays = c(), 
                              absenceDays = NULL)
  
  #days
  expect_equal((res %>% filter(PATIENT == 1))$treatmentDaysFullTime,14) 
  expect_equal((res %>% filter(PATIENT == 1))$treatmentDaysStartEnd,14) 
  expect_equal((res %>% filter(PATIENT == 4))$treatmentDaysFullTime,28) 
  expect_equal((res %>% filter(PATIENT == 4))$treatmentDaysStartEnd,28) 
  expect_equal((res %>% filter(PATIENT == 5))$treatmentDaysFullTime,0) 
  expect_equal((res %>% filter(PATIENT == 5))$treatmentDaysStartEnd,0) 
  expect_equal((res %>% filter(PATIENT == 2))$treatmentDaysFullTime,12) 
  expect_equal((res %>% filter(PATIENT == 2))$treatmentDaysStartEnd,7) 
  expect_equal((res %>% filter(PATIENT == 3))$treatmentDaysFullTime,21) 
  expect_equal((res %>% filter(PATIENT == 3))$treatmentDaysStartEnd,21) 
  expect_equal((res %>% filter(PATIENT == 6))$treatmentDaysFullTime,14) 
  expect_equal((res %>% filter(PATIENT == 6))$treatmentDaysStartEnd,7) 
  
  #adherence
  expect_true(is.na((res %>% filter(PATIENT == 5))$adherenceFullTime)) 
  expect_true(is.na((res %>% filter(PATIENT == 5))$adherenceStartEnd))
  expect_equal((res %>% filter(PATIENT == 1))$adherenceFullTime,100) 
  expect_equal((res %>% filter(PATIENT == 1))$adherenceStartEnd,100) 
  expect_equal((res %>% filter(PATIENT == 2))$adherenceFullTime,100) 
  expect_true(is.na((res %>% filter(PATIENT == 2))$adherenceStartEnd)) 
  expect_equal((res %>% filter(PATIENT == 3))$adherenceFullTime,100) 
  expect_equal((res %>% filter(PATIENT == 3))$adherenceStartEnd,100) 
  expect_equal((res %>% filter(PATIENT == 4))$adherenceFullTime,50) 
  expect_equal((res %>% filter(PATIENT == 4))$adherenceStartEnd,50)
  expect_equal((res %>% filter(PATIENT == 6))$adherenceFullTime,100) 
  expect_true(is.na((res %>% filter(PATIENT == 6))$adherenceStartEnd)) 
  
  #check that no extra treatment rows are added as this is only one medication
  expect_equal(nrow(res),6)
  
})

test_that("one medication with breaks, no absence", {
  res <- pdc_treatment(serialDf = meds %>% filter(PATIENT %in% c(4)),
                        startDates = starts  %>% filter(PATIENT %in% c(4)), 
                        endDates = ends  %>% filter(PATIENT %in% c(4)), 
                        atcCode = c('C09'), 
                        refillPeriod = 7, 
                        treatmentBreakDays = c(3), 
                        absenceDays = NULL)
  expect_equal((res %>% filter(PATIENT == 4))$adherenceFullTime,round(100 * 14/17, 2)) 
  expect_equal((res %>% filter(PATIENT == 4))$adherenceStartEnd,round(100 * 14/17, 2)) 
  
  res <- pdc_treatment(serialDf = meds %>% filter(PATIENT %in% c(7)),
                        startDates = starts  %>% filter(PATIENT %in% c(7)), 
                        endDates = ends  %>% filter(PATIENT %in% c(7)), 
                        atcCode = c('C09'), 
                        refillPeriod = 3,
                        treatmentBreakDays = c(2), 
                        absenceDays = NULL)
  
  expect_equal((res %>% filter(PATIENT == 7))$adherenceFullTime,round(100 * 15/22, 2)) 
  expect_equal((res %>% filter(PATIENT == 7))$adherenceStartEnd,round(100 * 15/22, 2)) 
  expect_equal((res %>% filter(PATIENT == 7))$treatmentDaysFullTime,22)
  expect_equal((res %>% filter(PATIENT == 7))$treatmentDaysFullTime,22)
})
  
  
test_that("one medication with breaks and absences", {
  
  res <- pdc_treatment(serialDf = meds %>% filter(PATIENT %in% c(3,4,6,7)),
                        startDates = starts  %>% filter(PATIENT %in% c(3,4,6,7)), 
                        endDates = ends  %>% filter(PATIENT %in% c(3,4,6,7)), 
                        atcCode = c('C09'), 
                        refillPeriod = 3,
                        treatmentBreakDays = c(2), 
                        absenceDays = abs)
  
  expect_equal((res %>% filter(PATIENT == 6))$adherenceFullTime,round(100 * 6/7, 2)) 
  expect_equal((res %>% filter(PATIENT == 6))$adherenceStartEnd,round(100 * 6/7, 2)) 
  expect_equal((res %>% filter(PATIENT == 6))$treatmentDaysFullTime,7) 
  expect_equal((res %>% filter(PATIENT == 6))$treatmentDaysStartEnd,7) 
  expect_equal((res %>% filter(PATIENT == 4))$adherenceFullTime,100) 
  expect_equal((res %>% filter(PATIENT == 4))$adherenceStartEnd,100) 
  expect_equal((res %>% filter(PATIENT == 4))$treatmentDaysFullTime,6) 
  expect_equal((res %>% filter(PATIENT == 4))$treatmentDaysStartEnd,6)
  expect_equal((res %>% filter(PATIENT == 7))$adherenceFullTime,round(100 * 12/15, 2)) 
  expect_equal((res %>% filter(PATIENT == 7))$adherenceStartEnd,round(100 * 12/15, 2)) 
  expect_equal((res %>% filter(PATIENT == 7))$treatmentDaysFullTime,15) 
  expect_equal((res %>% filter(PATIENT == 7))$treatmentDaysStartEnd,15) 
})

test_that("polypharmacy with breaks and absences", {
  
  res <- pdc_treatment(serialDf = meds %>% filter(PATIENT %in% c(7,8)),
                        startDates = starts  %>% filter(PATIENT %in% c(7,8)), 
                        endDates = ends  %>% filter(PATIENT %in% c(7,8)), 
                        atcCode = c('C09','C10'), 
                        refillPeriod = 3,
                        treatmentBreakDays = c(2,2), 
                        absenceDays = abs)
  
  expect_equal((res %>% filter(PATIENT == 7 & treatment == 'polyPharmacy'))$adherenceFullTime,round(100 * 12/15, 2)) 
  expect_equal((res %>% filter(PATIENT == 7 & treatment == 'polyPharmacy'))$adherenceStartEnd,round(100 * 12/15, 2)) 
  expect_equal((res %>% filter(PATIENT == 7 & treatment == 'polyPharmacy'))$treatmentDaysFullTime,15) 
  expect_equal((res %>% filter(PATIENT == 7 & treatment == 'polyPharmacy'))$treatmentDaysStartEnd,15) 
  expect_equal((res %>% filter(PATIENT == 8 & treatment == 'polyPharmacy'))$adherenceFullTime,round(100 * 11/16, 2)) 
  expect_equal((res %>% filter(PATIENT == 8 & treatment == 'polyPharmacy'))$adherenceStartEnd,round(100 * 11/16, 2)) 
  expect_equal((res %>% filter(PATIENT == 8 & treatment == 'polyPharmacy'))$treatmentDaysFullTime,16) 
  expect_equal((res %>% filter(PATIENT == 8 & treatment == 'polyPharmacy'))$treatmentDaysStartEnd,16) 
})

test_that("polypharmacy some individuals with breaks and absences", {
  
  res <- pdc_treatment(serialDf = meds %>% filter(PATIENT %in% c(4,8)),
                       startDates = starts  %>% filter(PATIENT %in% c(4,8)), 
                       endDates = ends  %>% filter(PATIENT %in% c(4,8)), 
                       atcCode = c('C09','C10'), 
                       refillPeriod = 3,
                       treatmentBreakDays = c(2,2), 
                       absenceDays = abs)
  
  expect_equal((res %>% filter(PATIENT == 4 & treatment == 'polyPharmacy'))$adherenceFullTime,100) 
  expect_equal((res %>% filter(PATIENT == 4 & treatment == 'polyPharmacy'))$adherenceStartEnd,100) 
  expect_equal((res %>% filter(PATIENT == 8 & treatment == 'polyPharmacy'))$adherenceFullTime,round(100 * 11/16, 2)) 
  expect_equal((res %>% filter(PATIENT == 8 & treatment == 'polyPharmacy'))$adherenceStartEnd,round(100 * 11/16, 2)) 
})


