createCytoData <- function(cyto_data,cyto_type,Plate.SN = NULL, srcname = NULL, date = NULL) {
  
  cat(Plate.SN,cyto_type,"\n")
  
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
    stop("NA's found in well ID data")
  }
  
  # Get desired values
  tagPhrases = c("Corrected for Blank", "Corrected Optical Denisty 490 nm", "Corrected Optical Density 490 nm", "Corrected Fluorescence", )
  value_index <- matrix(data = NA_real_, nrow = 0, ncol = 2) # initialize value_index as empty
  i = 1
  while (nrow(value_index)==0) {
    if (i > length(tagPhrases)) {
      cat("no corrected for blank-corrected data found for",Plate.SN,cyto_type,"\n")
      print(cyto_data)
      value_row <- readline(prompt = "Enter the Row number in table above corresponding to well A1 of the blank-corrected values: ")
      value_col <- readline(prompt = "Enter the Column number in table above corresponding to well A1 of the blank-corrected values : ")
      value_index <- array(c((as.numeric(value_row)-1), (as.numeric(value_col)-1)), dim  = c(1,2), dimnames = list(NULL,c("row","col"))) # Back up to the Row and Col id info
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
    stop("Did not merge correctly with 48 rows")
  }
  
  # get numerical rowi from character Row
  longdat[, rowi := match(Row, LETTERS)]
  
  if (length(Plate.SN) != 1 || is.null(Plate.SN) || is.na(Plate.SN)) {
    stop("Plate cannot be determined from file name. ")
  }
  if (length(date) != 1 || is.null(date) || is.na(date)) {
    stop("Culture date cannot be determined from file name. ")
  }
  
  # Add info columns
  longdat[, src_acsn := cyto_type]
  longdat$srcf = srcname
  longdat$Plate.SN = Plate.SN
  longdat$date = date
  longdat[, coli := as.numeric(coli)]
  
  # reorder columns
  longdat <- longdat[,c("date","Plate.SN","treatment","rowi","coli","conc","rval","srcf","src_acsn")]
  return(longdat)
}