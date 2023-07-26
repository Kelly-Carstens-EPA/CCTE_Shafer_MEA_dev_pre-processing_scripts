#' Read the cytotoxicity data for a given cyto_type (AB or LDH) from an input data chunk for a single plate
#'
#' @param cyto_data data chunk from a file containing the chemical, concentration, and blank-corrected values for a single assay and plate
#' @param cyto_type type of cytotoxicity data in the input cyto_data (either AB or LDH)
#' @param Plate.SN plate serial number of plate represented in cyto_data
#' @param date culture data of data represented in cyto_data
#'
#' @return
#' @export
#'
#' @details The chemical, concentration, and blank-corrected data values are read. The following process is used to identify the indicies in the input data chunk that 
#' contain the desired data:
#' 
#' * Chemical data = 
#' 1) Search for word “Chemical”
#' 2) Search for word “Row” in same column as “Chemical”
#' 3) Take all data that is 1 – 6 rows below the word “Row” and in the same column as “Chemical” to + 8 rows out
#' 
#' * Concentration data = 
#' Same as previous, except replace “Chemical” with “Concentration mM”
#' 
#' * Blank-corrected values = 
#' 1) Same as previous, except “Chemical” is replaced with one of the following: c("Corrected for Blank", "Corrected Optical Denisty 490 nm", "Corrected Optical Density 490 nm", "Corrected Fluorescence")
#' 
#' 
#' Note that any negative blank-corrected values are set to 0.
#'
#' @examples
createCytoData <- function(cyto_data,cyto_type,Plate.SN = NULL, date = NULL) {
  
  # compound map
  compound_col <- which(cyto_data == "Chemical", arr.ind = T)[1,"col"]
  Row_row <- which(cyto_data[, compound_col] == "Row")[1] # find first the occurence of the phrase "Row" in the same column as "Chemical"
  compoundmap <- cyto_data[(Row_row + 1):(Row_row + 6), compound_col:(compound_col+8)] # get the treatment labels under "Row"
  colnames(compoundmap) <- cyto_data[Row_row, compound_col:(compound_col+8)]
  compoundmap[, setdiff(names(compoundmap),"Row")] <- sapply(compoundmap[, setdiff(names(compoundmap),"Row")], as.character) # each col of the df is read as a list element
  setDT(compoundmap)
  compoundmap <- melt(compoundmap, id.vars = "Row", variable.name = "coli", value.name = "treatment", variable.factor = FALSE)
  
  # concentrations
  conc_col <- which(cyto_data == "Concentration mM", arr.ind = T)[1,"col"] # "Concentration [micro]M" is read from excel into R as "Concentration mM"
  Row_row <- which(cyto_data[, conc_col] == "Row")[1] # find first the occurence of the phrase "Row" in the same column as "Concentration mM"
  concmap <- cyto_data[(Row_row + 1):(Row_row + 6), conc_col:(conc_col+8)] # get the conc labels under "Row"
  colnames(concmap) <- cyto_data[Row_row, conc_col:(conc_col+8)]
  concmap[, setdiff(names(concmap),"Row")] <- sapply(concmap[, setdiff(names(concmap),"Row")], as.numeric) # each col of the df is read as a list element
  setDT(concmap)
  concmap <- melt(concmap, id.vars = "Row", variable.name = "coli", value.name = "conc", variable.factor = FALSE)
  
  # checking that there were no issues with data offset in sheet
  if(any(is.na(compoundmap) | is.na(concmap))) {
    print(compoundmap)
    print(concmap)
    stop(paste("NA's found in well ID data for",Plate.SN,cyto_type))
  }
  
  # Get desired values
  tagPhrases = c("Corrected for Blank", "Corrected Optical Denisty 490 nm", "Corrected Optical Density 490 nm", "Corrected Fluorescence")
  value_index <- matrix(data = NA_real_, nrow = 0, ncol = 2) # initialize value_index as empty
  i = 1
  while (nrow(value_index)==0) {
    if (i > length(tagPhrases)) {
      stop("No corrected for blank-corrected data found for",Plate.SN,cyto_type,
           "\nUpdate the tagPhrases in createCytoData() to include the tagPhrase needed to identify the blank-corrected values for the current file")
      print(cyto_data)
    }
    else {
      value_index <- which(cyto_data == tagPhrases[i], arr.ind = TRUE)
      i = i+1
    }
  }
  value_row <- value_index[1,"row"] # in case there were multiple occurences, we want the first one
  value_col <- value_index[1,"col"]
  # find the first occurence of the phrase "Row" in the same column as tagPhrase
  Row_row <- which(cyto_data[value_row:nrow(cyto_data), value_col] == "Row")[1] + value_row - 1
  valuemap <- cyto_data[(Row_row + 1):(Row_row + 6), value_col:(value_col+8)]
  colnames(valuemap) <- cyto_data[Row_row, value_col:(value_col+8)]
  valuemap[, setdiff(names(valuemap),"Row")] <- sapply(valuemap[, setdiff(names(valuemap),"Row")], as.numeric) # each col of the df is read as a list element
  setDT(valuemap)
  valuemap <- melt(valuemap, id.vars = "Row", variable.name = "coli", value.name = "rval", variable.factor = FALSE)
  
  # moving this check to a separate function
  # if (any(is.na(valuemap))) {
  #   stop("valuemap contains NAs")
  # }
  
  # if any blank-corrected values are negative, then we should set these to zero
  valuemap[rval < 0, rval := 0.0]
  
  # merge data together!
  longdat <- Reduce(merge, x = list(compoundmap, concmap, valuemap))
  if (nrow(longdat) != 48) {
    print(compoundmap)
    print(concmap)
    print(valuemap)
    stop(paste("Did not merge correctly with 48 rows for",Plate.SN,cyto_type))
  }
  
  # get numerical rowi from character Row
  longdat[, rowi := match(Row, LETTERS)]
  
  if (length(date) != 1 || is.null(date) || is.na(date)) {
    stop(paste("Culture date cannot be determined from file name for",Plate.SN,cyto_type))
  }
  
  # Add info columns
  longdat$acsn <- cyto_type
  longdat$Plate.SN <- Plate.SN
  longdat$date <- date
  longdat[, coli := as.numeric(coli)]
  
  # reorder columns
  longdat <- longdat[,c("date","Plate.SN","treatment","rowi","coli","conc","rval","acsn")]
  return(longdat)
}
