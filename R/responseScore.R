#' Calculates a score based on principle component analysis
#'
#' The function calcuculates a score for biomarker data represented as time blocks. The score is calculated via PCA.
#'
#' @param responseDf the data frame that contains the time blocks and the biomarker data per time block
#' @param idColumn name of ID column
#' @param blockColumn name of time block column
#' @param sexColumn name of column that contains sex information (optional). If provided, the score will be adjusted for sex.
#' @param responseColumn names of the columns containing the patient block value for a biomarker (e.g. vector of string names). The score will be based only on these columns.
#' @param confounderColumns names of columns that contain confounders (optional). If provided, the score will be adjusted for these confounders using linear regression
#' @return the response data frame with an added column phenoScore that shows the computed score
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @examples
#' \dontrun{
#' results <- read_tsv('../data/responseSpreadsheet.txt')
#'
#' resultsCor <- calculateResponseScore(results,responseColumns = c('sbp','dbp','aerAll','acrAll','egfr'), sexColumn = 'SEX')
#' }
calculateResponseScore <- function(responseDf, idColumn = 'PATIENT', blockColumn = 'block', sexColumn = NA, responseColumns, confounderColumns = NA){
  if(length(responseColumns) < 2 ){
    stop('Provide at least 2 response columns')
  }
  responseDf <- checkResponseFormat(responseDf,idColumn = idColumn, blockColumn = blockColumn, sexColumn = sexColumn, responseColumns = responseColumns)
  if(is.null(responseDf)){
    stop('Incorrect response format')
  }
  
  if(is.na(sexColumn)){
    pc <- calculatePCA(responseDf = responseDf, responseColumns = responseColumns)
    responseDf <- calculateScore(responseDf = responseDf, responseColumns = responseColumns, firstPrincipalComponent = pc)
  } else {
   pc1 <- calculatePCA(responseDf = responseDf %>% filter(SEX == unique(responseDf$SEX)[1]), responseColumns = responseColumns)
   pc2 <- calculatePCA(responseDf = responseDf %>% filter(SEX == unique(responseDf$SEX)[2]), responseColumns = responseColumns)
   if(any(sign(pc1) != sign(pc2))){
     if(all(sign(pc1) != sign(pc2))){
      pc2 <- pc2 * (-1) 
     } else {
       stop('PCA for men and women have inconsistent directions')
     }
   }
   responseDf1 <- calculateScore(responseDf = responseDf %>% filter(SEX == unique(responseDf$SEX)[1]), responseColumns = responseColumns, firstPrincipalComponent = pc1)
   responseDf2 <- calculateScore(responseDf = responseDf %>% filter(SEX == unique(responseDf$SEX)[2]), responseColumns = responseColumns, firstPrincipalComponent = pc2)
   responseDf1$phenoScore <- scale(responseDf1$phenoScore)[,1]
   responseDf2$phenoScore <- scale(responseDf2$phenoScore)[,1]
   responseDf <- rbind(responseDf1,responseDf2) %>% as.data.frame()
   responseDf <- responseDf[order(responseDf$PATIENT, responseDf$BLOCK),]
  }
  responseDf <- responseDf %>% adjustScore(confounderColumns = confounderColumns, responseColumns = responseColumns)
  responseDf
}

#' Perform the principle component analysis using a correlation matrix
#'
#' @param responseDf the data frame that contains the time blocks and the biomarker data per time block
#' @param responseColumn names of the columns containing the patient block value for a biomarker (e.g. vector of string names). The score will be based only on these columns.
#' @return the first principal component
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
calculatePCA <- function(responseDf, responseColumns){
  #remove rows that only have NAs for biomarkers
  responseNoMissingAll <- responseDf[(is.na(responseDf[,responseColumns])) %>% apply(1,function(x,rl){length(which(x == TRUE)) < rl},length(responseColumns)),]
  blockData <- responseNoMissingAll %>% select(c('BLOCK',responseColumns))  
  
  #only select baseline blocks for PCA
  #selectedBlocks <- blockData %>% filter(BLOCK == 0) %>% select(-BLOCK)
  
  #select all blocks for PCA
  selectedBlocks <- blockData %>% select(-BLOCK)

  #standardise with overall mean and sd
  blockData <- blockData %>% select(-BLOCK)
  blockData <- blockData %>% lapply(scale) %>% as.data.frame()
  
  #use the cor way to calculate PCs
  principalComponents <- NULL
  #calculate correlation/covariance
  inputBaselineData <- cor(selectedBlocks, use = 'pairwise.complete.obs', method = 'spearman')
   
  #eigenvalues
  eigenData <- eigen(inputBaselineData)
  #if the matrix is not positive semi-definite, add some weight to the diagonal
  while(any(eigenData$values < 0)){
    diag(inputBaselineData) <- diag(inputBaselineData) + 0.01
    eigenData <- eigen(inputBaselineData)
  }
  principalComponents <- eigenData$vector
  if(is.null(principalComponents)){
    stop('PCA failed. Check your input data. The method impute is currently not implemented.')
  }
  
  #only look at the first PC
  principalComponents <- principalComponents[,1]
  
  principalComponents
}

#' Get the score from the first PC
#'
#' @param responseDf the data frame that contains the time blocks and the biomarker data per time block
#' @param responseColumn names of the columns containing the patient block value for a biomarker (e.g. vector of string names). The score will be based only on these columns.
#' @param firstPrincipalComponent a vector containing the first PC
#' @return the response data frame with an added column phenoScore that shows the computed score
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_all
calculateScore <- function(responseDf, responseColumns, firstPrincipalComponent){
  #remove rows that only have NAs for biomarkers
  responseNoMissingAll <- responseDf[(is.na(responseDf[,responseColumns])) %>% apply(1,function(x,rl){length(which(x == TRUE)) < rl},length(responseColumns)),]
  blockData <- responseNoMissingAll %>% select(c('BLOCK',responseColumns))  
  
  #standardise with overall mean and sd
  blockData <- blockData %>% select(-BLOCK)
  blockData <- blockData %>% lapply(scale) %>% as.data.frame()
  
  #set NAs to 0 in order to compute score (otherwise missing biomarker will lead to missing score)
  blockData <- blockData %>% mutate_all(funs(replace(., is.na(.), 0)))
  
  #score is a linear combination of the first component
  blockData <- blockData %>% mutate(phenoScore = data.matrix(blockData) %*% data.matrix(firstPrincipalComponent))
  
  responseNoMissingAll <- cbind(responseNoMissingAll,blockData$phenoScore) %>% as.data.frame()
  responseNoMissingAll <- responseNoMissingAll %>% select(PATIENT,BLOCK,phenoScore = `blockData$phenoScore`)
  responseDf <- merge(responseDf,responseNoMissingAll, by = c('PATIENT','BLOCK'), all = TRUE)
  responseDf  
}

#' Get the score from the first PC
#'
#' @param responseDf the data frame that contains the time blocks and the scores and the confounders per time block
#' @param scoreColumn name of the column that contains the score
#' @param confounderColumns names of the columns that contain confounding information
#' @param responseColumn names of the columns containing the patient block value for a biomarker (e.g. vector of string names). The score will be based only on these columns.
#' @return the response data frame with an adjusted score
#' @keywords internal
adjustScore <- function(responseDf, scoreColumn = 'phenoScore', confounderColumns, responseColumns){
  if(!any(is.na(confounderColumns))){
    myFormula <- reformulate(termlabels = c(0,confounderColumns), response = scoreColumn)
    fit <- lm(myFormula, data=responseDf, na.action = na.exclude)
    responseDf[,scoreColumn] <- ifelse(is.na(residuals(fit)),responseDf[,scoreColumn],residuals(fit))
  }
  responseDf
}



