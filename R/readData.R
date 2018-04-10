#' Read in a dataset
#'
#' Read in a tab-delimited dataset that also has a description file available with column types.
#'
#' The description file is a direct downlad from the BCOS table descritpions
#'
#' If only the fileName is given, the function behaves like read_tsv from readr
#'
#' @param fileName path to the tab-delimited data file.
#' @param descName path to the description file (optional)
#' @param doRemove removes PATIENT entries that start with !, DSP, TH or CC
#' @param guessRow when no description file is given, how many rows should be used to guess the column type (default: 1000)
#' @return The dataset as a data frame.
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom readr problems
#' @importFrom dplyr select
#' @importFrom dplyr filter

readData <- function(fileName, descName = NA, doRemove = TRUE, guessRows = NA) {
  if (!is.na(descName)) {
    df <- read_tsv(fileName, n_max = 1)
    dfDesc <- read_tsv(descName) %>% select(QUESTIONID, TYPE)
    if (length(which(names(df) == dfDesc$QUESTIONID)) != nrow(dfDesc)) {
      stop("Column headings in data file and in description file have different length")
    }
    myColTypes <- dfDesc$TYPE %>% as.vector()
    myColTypes <- gsub("ALT", "i", myColTypes)
    myColTypes <- gsub("NUM", "d", myColTypes)
    myColTypes <- gsub("DATE", "D", myColTypes)
    myColTypes <- gsub("TEXT(.)*", "c", myColTypes)
    df <- read_tsv(fileName, col_types = paste(myColTypes, collapse = ""))
  } else {
    df <- read_tsv(fileName, guess_max = ifelse(is.na(guessRows), 1000, guessRows))
  }
  if (nrow(problems(df)) > 0) {
    print(problems(df))
    stop("Problems while reading data.")
  }
  if (doRemove == TRUE) {
    df <- df %>% filter(substr(PATIENT, 1, 1) != "!")
    df <- df %>% filter(substr(PATIENT, 1, 3) != "DSP")
    df <- df %>% filter(substr(PATIENT, 1, 2) != "TH")
    df <- df %>% filter(substr(PATIENT, 1, 2) != "CC")
  }
  df
}
