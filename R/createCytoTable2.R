# splits the cyto data for each plate, assay, then calls createCytoData()
createCytoTable2 <- function(cytoFile, cyto_type) {
  
  # get the date, file name
  srcname <- basename(cytoFile)
  srcname <- sub(" ","_",srcname) # sometimes file contains space instead of _, esp in Summary files
  namesplit <- strsplit(srcname, split ="_")[[1]]
  date <- grep(pattern="[0-9]{8}",namesplit, val = T) # looks for 8-digit string in namesplit
  
  # Determine which tab in cytoFile contains the data for the given cyto_type
  ABnames <- c("Alamar Blue", "AB", "CTB", "CellTiter Blue")
  LDHnames <- c("LDH", "Total LDH")
  sheetNames <- switch(cyto_type,
                     AB = ABnames,
                     LDH = LDHnames)
  cytoFile_sheets <- getSheetNames(cytoFile)
  sheetName_for_cyto_type <- intersect(sheetNames, getSheetNames(cytoFile_sheets))
  if (length(sheetName_for_cyto_type) != 1) {
    stop('Not sure which sheet in ',basename(cytoFile),' to read for ',cyto_type,' data.\nAvailable sheets: ',paste0(cytoFile_sheets,collapse = ", "))
  }
  cyto_data_all <- read.xlsx(cytoFile, sheet = sheetName_for_cyto_type, colNames = FALSE)
  
  # if the input file was a "Summary" file, should see the phrase "Chemical" only once
  if (length(which(cyto_data_all == "Chemical")) == 1) {
    cyto_data_list <- list(cyto_data_all)
  } else{
    # "Calculations" file with 3 plates
    # split the data wherever there is a new occurrence of the word "Chemical" in the first column
    # Only want the first 3 - because these should correspond to the first 3 plates
    # Any occurences of "chemical" after that are probably other calculations
    cyto_plate_indicies <- which(cyto_data_all[,1] == "Chemical")
    cyto_plate_indicies <- cyto_plate_indicies[1:min(3,length(cyto_plate_indicies))]
    cyto_data_list <- lapply(cyto_plate_indicies, function(i) cyto_data_all[i:(i+9),])
  }
  
  # Loop through each data chunk in cyto_data_list
  # Read the data with createCytoData and put into table
  file_dat <- data.table()
  
  for (cyto_data in cyto_data_list) {
    
    # get Plate.SN
    plateindex <- which(cyto_data == "Plate", arr.ind = T)
    if (nrow(plateindex) == 0) {
      # then get plate from file name
      Plate.SN <-  grep("-", namesplit, val = T)
      Plate.SN <- sub("Summary\\.xlsx","",Plate.SN) # get rid of file name ending if present
    } else {
      Plate.SN <- paste("MW",cyto_data[plateindex[1,"row"],(plateindex[1,"col"]+1)], sep = "")
    }
    
    if (all(is.na(cyto_data))) {
      cat("**No data found for ",basename(cytoFile),", ",cyto_type,", reading from sheet ",sheetName_for_cyto_type,", Plate ",Plate.SN,"\n")
    } else {
      cyto_longdat <- createCytoData(cyto_data, cyto_type, Plate.SN, srcname, date, masterChemFiles)
    }
    file_dat <- rbind(file_dat, cyto_longdat)
  }
  return(file_dat)
}
