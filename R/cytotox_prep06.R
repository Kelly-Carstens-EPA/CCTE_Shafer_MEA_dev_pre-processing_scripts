# This script will extract the blank-corrected flurescense values (instead of the percent of control)
# Script to process cytotoxicity data to prepare for tcpl
# Output is a long file with all the necessary columns for tcpl (except no spid, just treatment column)
# This script will extract the cytotoxicity data (LDH and Alamar Blue) from the input files
# Example input files:
# - "ON_20160720_MW1139-19_Summary.xlsx" - these contain data for 1 plate per sheet (LDH and Alamar Blue)
# - "20171011_ON G8_2 Calculations.xlsx" - these contain data for 3 plates per sheet (LDH and Alamar Blue)

###################################################################################
# USER INPUT
###################################################################################
# ALL OF THESE ARE NOW CONTROL IN THE FUNCTION
# 'run_cytotox_functions()' below

# set the location for the output file
# basepath <- main.output.dir

# set the name of the output file
# filename = paste0(project_name,"_cytotox_longfile.csv")

# check if the filename already exists first?
# remake_all <- TRUE

# append to existing data?
# append <- FALSE
###################################################################################
# END USER INPUT
###################################################################################

require(openxlsx)
require(data.table)

###################### FUNCTIONS

# function to  look up index of value in a data frame
#  there probably is a more efficient, data table way to do this...
returnindex = function(value, mydata) {
  for (i in 1:nrow(mydata)) {
    for (j in 1:ncol(mydata)) {
      if (is.na(mydata[i,j])) {
        next
      }
      if (strcmp(as.character(mydata[i,j]),as.character(value))) {
        return(c(i,j))
      }
    }
  }
  # print("could not find index in data frame")
  return(NULL)
}

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


# splits the cyto data for each plate, assay, then calls createCytoData()
createCytoTable2 <- function(sourcefile, cyto_type, masterChemFiles = c()) {

  # get the date, file name
  srcname <- basename(sourcefile)
  srcname <- sub(" ","_",srcname) # sometimes file contains space instead of _, esp in Summary files
  namesplit <- strsplit(srcname, split ="_")[[1]]
  date <- grep(pattern="[0-9]{8}",namesplit, val = T) # looks for 8-digit string in namesplit

  cyto_data_all <- findTabData(sourcefile, cyto_type)
  if (nrow(cyto_data_all) == 0) {
    # is user elected to "skip" this plate/cyto_type
    return(cyto_data_all)
  }

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
      cat(cyto_type,"data not found\n")
    } else {
      cyto_longdat <- createCytoData(cyto_data, cyto_type, Plate.SN, srcname, date, masterChemFiles)
    }
    file_dat <- rbind(file_dat, cyto_longdat)
  }
  return(file_dat)
}

createCytoData <- function(sourcedata,cyto_type,Plate.SN = NULL, srcname = NULL, date = NULL, masterChemFiles = c()) {

  cat(Plate.SN,cyto_type,"\n")

  # compound map
  compound_col <- which(sourcedata == "Chemical", arr.ind = T)[1,"col"]
  Row_row <- which(sourcedata[, compound_col] == "Row")[1] # find first the occurence of the phrase "Row" in the same column as "Chemical"
  compoundmap <- sourcedata[(Row_row + 1):(Row_row + 6), compound_col:(compound_col+8)] # get the treatment labels under "Row"
  colnames(compoundmap) <- sourcedata[Row_row, compound_col:(compound_col+8)]
  compoundmap[, setdiff(names(compoundmap),"Row")] <- sapply(compoundmap[, setdiff(names(compoundmap),"Row")], as.character) # each col of the df is read as a list element
  setDT(compoundmap)
  compoundmap <- melt(compoundmap, id.vars = "Row", variable.name = "coli", value.name = "treatment", variable.factor = FALSE)

  # concentrations
  conc_col <- which(sourcedata == "Concentration mM", arr.ind = T)[1,"col"] # "Concentration [micro]M" is read from excel into R as "Concentration mM"
  Row_row <- which(sourcedata[, conc_col] == "Row")[1] # find first the occurence of the phrase "Row" in the same column as "Concentration mM"
  concmap <- sourcedata[(Row_row + 1):(Row_row + 6), conc_col:(conc_col+8)] # get the conc labels under "Row"
  colnames(concmap) <- sourcedata[Row_row, conc_col:(conc_col+8)]
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
  tagPhrases = c("Corrected for Blank", "Corrected Optical Denisty 490 nm", "Corrected Fluorescence")
  value_index <- matrix(data = NA_real_, nrow = 0, ncol = 2) # initialize value_index as empty
  i = 1
  while (nrow(value_index)==0) {
    if (i > length(tagPhrases)) {
      cat("no corrected for blank-corrected data found for",Plate.SN,cyto_type,"\n")
      print(sourcedata)
      value_row <- readline(prompt = "Enter the Row number in table above corresponding to well A1 of the blank-corrected values: ")
      value_col <- readline(prompt = "Enter the Column number in table above corresponding to well A1 of the blank-corrected values : ")
      value_index <- array(c((as.numeric(value_row)-1), (as.numeric(value_col)-1)), dim  = c(1,2), dimnames = list(NULL,c("row","col"))) # Back up to the Row and Col id info
    }
    else {
      value_index <- which(sourcedata == tagPhrases[i], arr.ind = TRUE)
      i = i+1
    }
  }
  value_row <- value_index[1,"row"] # in case there were multiple occurences, we want the first one
  value_col <- value_index[1,"col"]
  # find the first occurence of the phrase "Row" in the same column as tagPhrase
  Row_row <- which(sourcedata[value_row:nrow(sourcedata), value_col] == "Row")[1] + value_row - 1
  valuemap <- sourcedata[(Row_row + 1):(Row_row + 6), value_col:(value_col+8)]
  colnames(valuemap) <- sourcedata[Row_row, value_col:(value_col+8)]
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
    cat("Plate cannot be determined from file name. ")
    Plate.SN = readline("Enter plate sn: ")
  }
  if (length(date) != 1 || is.null(date) || is.na(date)) {
    cat("Culture date cannot be determined from file name. ")
    date = readline("Enter culture date (as yyyymmdd): ")
  }

  # determine correct assay component source name
  longdat[, src_acsn := cyto_type]
  longdat$srcf = srcname
  longdat$Plate.SN = Plate.SN
  longdat$date = date
  longdat[, coli := as.numeric(coli)]

  # if provided, replace the treatment names with the names in the master chemical lists
  if (length(masterChemFiles) != 0) {
    # Get the masterChemFile with the same plate number
    masterChemFile = grep(pattern = paste0(date,"_",Plate.SN,"[_ ]"), masterChemFiles, value = T)
    if (length(masterChemFile) == 0) {
      # sometimes the master chem file does not include "MW" in the file name
      masterChemFile <- grep(pattern = paste0(date,"_",sub("MW","",Plate.SN),"[_ ]"), masterChemFiles, value = T)
    }
    # sometimes, maestroExperimentLog does not contain the plate.SN in the file name. Match by plate folder and date in file name instead
    if (length(masterChemFile) == 0) {
      masterChemFile <- Filter(function(mcf) grepl(sub("MW( )*","",Plate.SN),mcf) && grepl(date,basename(mcf)), masterChemFiles)
    }

    # If still no match found
    if (length(masterChemFile) != 1) {
      warning(paste("master chem file match not found for",Plate.SN,sep = " "))
    }
    else {
      # removing this 4/20/21 -> want to keep trt names from Calc file as a secondary check
      # but, still want to check from matching plate name, so keeping steps above
      # masterChemData = as.data.table(read.csv(masterChemFile, stringsAsFactors = FALSE))
      # masterChemData[, `:=`(coli = as.numeric(sub("[[:alpha:]]","",Well)), rowi = match(sub("[[:digit:]]","",Well), LETTERS),
      #                       date = as.character(Experiment.Date))]
      #
      # # only replacing the treatment names for now, might add conc's in the future
      # longdat[, treatment := NULL] # remove current treatment column
      # longdat <- merge(longdat, masterChemData[, .(date,Plate.SN,rowi,coli,treatment = Treatment)], by = c("date","Plate.SN","rowi","coli"))
    }
  }

  # reorder columns
  longdat <- longdat[,c("date","Plate.SN","treatment","rowi","coli","conc","rval","srcf","src_acsn")]
  return(longdat)
}


wllq_updates_cytotox <- function(longdat, basepath = NULL, get_files_under_basepath = TRUE) {
  if (get_files_under_basepath) wllq_info <- as.data.table(read.csv(file.path(basepath, "wells_with_well_quality_zero.csv"),
                                                                    colClasses = c(rep("character",4),"numeric",rep("character",2)), stringsAsFactors = FALSE))
  else wllq_info <- as.data.table(read.csv(choose.files(default = basepath, caption = "Select well quality csv table"), stringsAsFactors = FALSE))

  longdat[, `:=`(date = as.character(date), rowi = as.numeric(rowi))]

  # checking for typo's under "affected_endpoints"
  if(nrow(wllq_info[!grepl("(mea)|(CTB)|(LDH)",affected_endpoints)])>0) {
    cat("The following rows don't match any of the expected affected_endpoints (mea,CTB,LDH):\n")
    print(wllq_info[!grepl("(mea)|(CTB)|(LDH)",affected_endpoints)])
    stop("Update wells_with_well_quality_zero.csv")
  }

  # expand the rows where well == "all" to include every well in plate
  for(table_row in which(wllq_info$well == "all")) {
    full_plate <- wllq_info[table_row, .(date, Plate.SN, DIV, well = paste0(unlist(lapply(LETTERS[1:6],rep,times=8)),rep(1:8,6)), wllq, wllq_notes, affected_endpoints)]
    wllq_info <- rbind(wllq_info, full_plate) # rbind the new rows to the end of wllq_info
  }
  wllq_info <- wllq_info[well != "all"]

  # add rowi and coli to wllq_info
  wllq_info[, well := paste0(toupper(substring(well,1,1)), substring(well,2,2))] # ensure that all well ID's start with a capital letter
  wllq_info[, `:=`(coli = as.numeric(sub("[[:alpha:]]","",well)), rowi = match(sub("[[:digit:]]","",well), LETTERS))]

  # Check for any rows where there may have been a typo in wllq_info
  unmatched_wells <- merge(longdat[, .(date, Plate.SN, rowi, coli, treatment, conc, srcf)],
                           wllq_info[grepl("(LDH)|(CTB)|(AB)",affected_endpoints), .(date, Plate.SN, rowi, coli, wllq, wllq_notes, affected_endpoints)], all.y = T)[is.na(srcf)]
  if(nrow(unmatched_wells) > 0) {
    cat("The following rows from wells_with_well_quality_zero.csv did not match any rows in longdata:\n")
    print(unmatched_wells)
    cat("\nSummary of longdata:\n")
    print(longdata[, .(plates = paste0(sort(unique(Plate.SN)),collapse=",")), by = "date"][order(date)])
    stop("Update wells_with_well_quality_zero.csv")
  }

  # transfrom wllq_info into longdat
  ldh_wllq <- wllq_info[grepl("LDH",affected_endpoints), .(wllq = min(wllq),
                                                           wllq_notes = paste0(unique(wllq_notes), collapse = "; ")),
                        by = .(date, Plate.SN, rowi, coli)] # collapsing in case of multiple wllq notes for a given well for LDH
  ldh_wllq[, src_acsn := grep("LDH",unique(longdat$src_acsn),val = T)]
  ctb_wllq <- wllq_info[grepl("(CTB)|(AB)",affected_endpoints), .(wllq = min(wllq),
                                                                  wllq_notes = paste0(unique(wllq_notes), collapse = "; ")),
                        by = .(date, Plate.SN, rowi, coli)] # collapsing in case of multiple wllq notes for a given well for
  ctb_wllq[, src_acsn := grep("AB",unique(longdat$src_acsn),val = T)]

  # set the wllq in longdat
  longdat <- merge(longdat, rbind(ldh_wllq, ctb_wllq), all.x = TRUE)
  longdat[is.na(wllq), `:=`(wllq = 1, wllq_notes = "")]

  # flag NA rval's where there is no wllq note
  na_rvals <- longdat[is.na(rval) & wllq == 1]
  if (nrow(na_rvals) > 0) {
    cat("The following rval's are missing in sourcedata, but wllq==1\n")
    print(na_rvals)
    resp <- readline(prompt = "Do you wish to set wllq=0 for all of these wells? (Only do this if you know that these values should be NA) (y/n): ")
    if (resp %in% c("y","Y","yes","Yes")) {
      longdat[is.na(rval) & wllq == 1, `:=`(wllq = 0, wllq_notes = "rval is NA; ")]
    }
    else {
      stop("Update wllq_info table")
    }
  }

  # Label the wllq as "wllq_by_well", because further wllq adjustments may be made once merge in "well_quality_notes_per_culture_treatment_cndx.xlsx"
  # (based on treatment and cdx) after have verified meta data
  setnames(longdat, old = 'wllq', new = 'wllq_by_well')
  setnames(longdat, old = 'wllq_notes', new = 'wllq_notes_by_well')

  # summary of wllq updates
  cat("Wllq summary:\n")
  print(longdat[, .N, by = c("src_acsn","wllq_by_well","wllq_notes_by_well")][order(src_acsn, wllq_by_well, wllq_notes_by_well)])

  return(longdat)
}


run_cytotox_functions <- function(basepath, get_files_from_log = TRUE, filename = "cytotox.csv", copy_maestro_exp_log_treatments = TRUE, append = FALSE) {

  cat("\nStarting cytotoxicity data collection...\n")
  cat("Any negative blank-corrected values will be set to 0.\n")

  # set up dir and output_file
  if (!dir.exists(file.path(basepath, "output"))) dir.create(file.path(basepath, "output"))
  # might be a better way to do this... I just don't want to have to pass it thru 4 versions of 3 function calls
  assign("output_file", value = file.path(basepath, "output",filename), envir = .GlobalEnv)

  # if(!append && file.exists(output_file)) {
  #   return(paste0(basename(output_file), " already exits."))
  # }

  # get the source files, either from log file or by selecting
  if(get_files_from_log) {
    cytoFiles <- readFilesLog(basepath, files_type = "Calculations")
    cytoFiles <- c(cytoFiles, readFilesLog(basepath, files_type = "Summary"))
    if(copy_maestro_exp_log_treatments) masterChemFiles <- readFilesLog(basepath, files_type = "MaestroExperimentLog")
    else masterChemFiles <- c()
  }
  else {
    cytoFiles <- choose.files(caption = paste("Select all Summary and/or Calculations files containing cytotoxicity data for 1-3 plates per sheet",sep=" "))
    if(copy_maestro_exp_log_treatments) masterChemFiles <- choose.files(caption = "(optional) Select all Master Chemical Lists Files fot these plates")
    else masterChemFiles <- c()
  }

  if (length(masterChemFiles) == 0) {
    print("no master chem lists available, using treatment names from input files")
  }

  # if (append) {
  #   cyto_dat <- as.data.table(read.csv(output_file, stringsAsFactors = F))
  #   setdiff(basename(cytoFiles), unique(cyto_dat$srcf)) # but ah, I have replaced all " " with "_"
  #   # and, what if I renamed the summary file,a dn that is why I am re-running? #then need tomanually delete...
  #
  #   # other method:
  #   completed_dates <- unique(cyto_dat$date)
  #   cytoFiles <- Filter(function(filei) {
  #     datei <- strsplit(basename(filei), split="_")[[1]][2]
  #     !(datei %in% completed_dates) }, cytoFiles)
  #   # how to make sure the entire date was completed though, for all plates?
  # }

  if (length(cytoFiles) == 0) {
    cat('No Summary or Calculations files found\n')
    return(0)
  }

  # run the functions for each file
  longdat <- data.table()
  for (i in 1:length(cytoFiles)) {

    if (grepl("(Calculations)|(Summary)",basename(cytoFiles[i]))) {
      AB_dat <- createCytoTable2(cytoFiles[i], cyto_type = "AB", masterChemFiles)
      LDH_dat <- createCytoTable2(cytoFiles[i], cyto_type = "LDH", masterChemFiles)
    }
    else {
      cat(paste("can't tell if",cytoFiles[i],"is 'Summary' file or 'Calculations' file\n"))
    }
    longdat <- rbind(longdat, AB_dat, LDH_dat)
    rm(list = c("AB_dat","LDH_dat"))
  }

  # set the wllq
  longdat <- wllq_updates_cytotox(longdat, basepath, get_files_under_basepath = get_files_from_log)
  write.table(longdat, file = output_file, sep = ",", row.names = FALSE, col.names = !append, append = append)

  cat(basename(output_file),"is ready\n")

}

###################### END FUNCTIONS

# example exectue:
# run_cytotox_functions(basepath = getwd(), get_files_from_log = FALSE)
# will allow you to select files






