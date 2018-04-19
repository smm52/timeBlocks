#' Calculates a score based on principle component analysis
#'
#' The function calcuculates a score for biomarker data represented as time blocks. The score is calculated via PCA.
#'
#' @param responseDf the data frame that contains the time blocks and the biomarker data per time block
#' @param idColumn name of ID column
#' @param blockColumn name of time block column
#' @param responseColumn names of the columns containing the patient block value for a biomarker (e.g. vector of string names). The score will be based only on these columns.
#' @param method string that indiciates how the PCA should be calculated. Values are cor for correlation matrix, cov for covariance matrix and impute to impute missing values via k-nearest neighbour and the use base R PCA (currently not implemented)
#' @return the response data frame with an added column phenoScore that shows the computed score
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_all
#' @examples
#' \dontrun{
#' results <- read_tsv('../data/responseSpreadsheet.txt')
#'
#' resultsCov <- calculateResponseScore(results,responseColumns = c('sbp','dbp','aerAll','acrAll','egfr'), method = 'cov')
#' resultsCor <- calculateResponseScore(results,responseColumns = c('sbp','dbp','aerAll','acrAll','egfr'), method = 'cor')
#' }
calculateResponseScore <- function(responseDf, idColumn = 'PATIENT', blockColumn = 'block', responseColumns, method = 'cor'){
  if(length(responseColumns) < 1 ){
    stop('No response columns provided')
  }
  if(!(method %in% c('cor','cov','impute'))){
    stop('Principal component are either calculated via the (spearman) correlation (cor) or covariance (cov) matrix through their eigenvectors or via imputing missing data nd base R (impute).')
  }
  responseDf <- checkResponseFormat(responseDf,idColumn = idColumn, blockColumn = blockColumn, responseColumns = responseColumns)
  if(is.null(responseDf)){
    stop('Incorrect response format')
  }

  #remove rows that only have NAs for biomarkers
  responseNoMissingAll <- responseDf[(is.na(responseDf[,responseColumns])) %>% apply(1,function(x,rl){length(which(x == TRUE)) < rl},length(responseColumns)),]
  blockData <- responseNoMissingAll %>% select(c('BLOCK',responseColumns))
  #only select baseline blocks for PCA
  baselineBlocks <- blockData %>% filter(BLOCK == 0) %>% select(-BLOCK)
  baselineBlocks <- baselineBlocks[,responseColumns]
  #remember mean and sd
  baselineMeans <- baselineBlocks %>% apply(2,mean,na.rm = TRUE)
  baselineSDs <- baselineBlocks %>% apply(2,sd,na.rm = TRUE)

  #standardise with the baseline mean and sd
  blockData <- blockData %>% select(-BLOCK)
  blockData[,responseColumns] <-
    (blockData[,responseColumns] - t(replicate(nrow(blockData),baselineMeans))) / (t(replicate(nrow(blockData),baselineSDs)))

  #use the cor/cov way to calculate PCs or impute missing data and then use base R on complete data
  principalComponents <- NULL
  if(method %in% c('cor','cov')){
    #calculate correlation/covariance
    if(method == 'cor'){
      inputBaselineData <- cor(baselineBlocks, use = 'pairwise.complete.obs', method = 'spearman')
    } else if(method == 'cov'){
      inputBaselineData <- cov(baselineBlocks, use = 'pairwise.complete.obs', method = 'pearson') #cov and pairwise.complete requires pearson
    }
    #eigenvalues
    eigenData <- eigen(inputBaselineData)
    #if the matrix is not positive semi-definite, add some weight to the diagonal
    while(any(eigenData$values < 0)){
      diag(inputBaselineData) <- diag(inputBaselineData) + 0.01
      eigenData <- eigen(inputBaselineData)
    }
    principalComponents <- eigenData$vector
  } else if( method == 'impute'){
    principalComponents <- NULL
  }
  if(is.null(principalComponents)){
    stop('PCA failed. Check your input data. The method impute is currently not implemented.')
  }

  #only look at the first PC
  principalComponents <- principalComponents[,1]

  #set NAs to 0 in order to compute score (otherwise missing biomarker will lead to missing score)
  blockData <- blockData %>% mutate_all(funs(replace(., is.na(.), 0)))

  #score is a linear combination of the first component
  blockData <- blockData %>% mutate(phenoScore = data.matrix(blockData) %*% data.matrix(principalComponents))

  responseNoMissingAll <- cbind(responseNoMissingAll,blockData$phenoScore) %>% as.data.frame()
  responseNoMissingAll <- responseNoMissingAll %>% select(PATIENT,BLOCK,phenoScore = `blockData$phenoScore`)
  responseDf <- merge(responseDf,responseNoMissingAll, by = c('PATIENT','BLOCK'), all = TRUE)
  responseDf
}
