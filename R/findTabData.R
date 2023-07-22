# uses the list of possible tab names to get the AB or LDH data from excel sheet
findTabData <- function(sourcefile, assay = c("AB", "LDH")) {
  
  ABnames <- c("Alamar Blue", "AB", "CTB", "CellTiter Blue")
  LDHnames <- c("LDH", "Total LDH")
  
  tabNames <- switch(assay,
                     AB = ABnames,
                     LDH = LDHnames)
  
  tabName <- intersect(tabNames, getSheetNames(sourcefile))
  if (length(tabName) != 1) {
    tabName <- readline(prompt = paste0("Enter name of tab in ", basename(sourcefile)," for ",assay," data, or 'skip' to not include any ",assay," data for this plate: "))
  }
  if (tabName == "skip" ) {
    my_data <- data.frame()
  } else {
    my_data <- read.xlsx(sourcefile, sheet = tabName, colNames = FALSE)
  }
  return(my_data)
}