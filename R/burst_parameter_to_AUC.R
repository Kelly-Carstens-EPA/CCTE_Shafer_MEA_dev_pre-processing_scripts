## MEA neuronal ontogeny AUC calculations
## Adapted from Chris Frank - February-March 2016
## Adapated from Mahmoud 2019 
## Last edited by Amy Carpenter August 2020

##########################################################
# USER INPUTS
##########################################################
# set location where output file should be created
basepath <- main.output.dir
# Set output file name
filename <- paste0(project_name,"_AUC.csv")
get_files_under_basepath <- TRUE
# set the DIVs that should be included
use_divs <- c(5,7,9,12) # note that a point of DIV = 2 value = 0 will be added for every endpoint regardless
interpolate_diff_divs <- TRUE
##########################################################
# END USER INPUTS
##########################################################

require(data.table)
if(!dir.exists(file.path(basepath,"output"))) dir.create(file.path(basepath,"output"))

cat("\nStarting AUC preparations...\n")
if(get_files_under_basepath) {
    # get files with the 16 parameters for all culture dates
    parameter_data_files <- list.files(path = file.path(basepath, "prepared_data"), pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
    cat("Got",length(parameter_data_files),"parameter csv files from",file.path(basepath, "prepared_data"),"\n")
    
    # get files with mi data for all culture dates
    mi_data_files <- list.files(path = file.path(basepath, "All_MI"), pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
    cat("Got",length(mi_data_files),"MI csv files from",file.path(basepath, "All_MI"),"\n")
    
} else {
  # get files with the 16 parameters for all culture dates
  parameter_data_files <- choose.files(default = basepath, caption = "Select all files containing activity parameter values to calculate AUC")
  cat("Got",length(parameter_data_files),"parameter data files.\n")
  # get files with mi data for all culture dates
  mi_data_files <- choose.files(default = basepath, caption = "Select all files containing MI data to calculate AUC")
  cat("Got",length(mi_data_files),"MI data files.\n")
}

options(digits = 6)

# read in the data -------------------------------------------------------------------
parameter_data <- c()
for (i in 1:length(parameter_data_files)) {
  parameter_datai <- read.csv(parameter_data_files[i], stringsAsFactors = F)
  parameter_data <- rbind(parameter_data, parameter_datai)
  rm(parameter_datai)
}
mi_data <- c()
for (i in 1:length(mi_data_files)) {
  mi_datai <- read.csv(mi_data_files[i], stringsAsFactors = F)
  mi_data <- rbind(mi_data, mi_datai)
  rm(mi_datai)
}

parameter_data <- unique(parameter_data) # eliminate duplicate data
mi_data <- unique(mi_data) # eliminate duplicate data
all_data <- merge(parameter_data, mi_data, by = c("date","Plate.SN","DIV","well","trt","dose","units","file.name")) # merge two data frames on common columns (plateID, date, well, treatment, source file name, etc.)
if (nrow(all_data) < nrow(parameter_data) | nrow(all_data) < nrow(mi_data)) {
  stop(paste0("Some rows of parameter_data and mi_data don't match or are missing\n"))
}
setDT(all_data)
all_data[, trt := as.character(trt)] # I keep forgetting how to do this more efficiently!!
all_data[, date := as.character(date)]


# update wllq -------------------------------------------------------------------------
cat("\nUpdating wllq...\n")
wllq_info <- as.data.table(read.csv(file.path(basepath, "wells_with_well_quality_zero.csv"), colClasses = c(rep("character",4),"numeric",rep("character",2)), stringsAsFactors = F))

# checking for typo's under "affected_endpoints"
if(nrow(wllq_info[!grepl("(mea)|(CTB)|(LDH)",affected_endpoints)])>0) {
  cat("The following rows don't match any of the expected affected_endpoints (mea,CTB,LDH):\n")
  print(wllq_info[!grepl("(mea)|(CTB)|(LDH)",affected_endpoints)])
  warning("Update wells_with_well_quality_zero.csv")
}

# expand the rows where well == "all" to include every well in plate
for(table_row in which(wllq_info$well == "all")) {
  full_plate <- wllq_info[table_row, .(date, Plate.SN, DIV, well = paste0(unlist(lapply(LETTERS[1:6],rep,times=8)),rep(1:8,6)), wllq, wllq_notes, affected_endpoints)]
  wllq_info <- rbind(wllq_info, full_plate) # rbind the new rows to the end of wllq_info
}
wllq_info <- wllq_info[well != "all"]

# check for any rows in wllq_info that dont' match a wells in all_data (Which may indicate a typo in wllq_info)
wllq_info[, well := paste0(toupper(substring(well,1,1)), substring(well,2,2))] # ensure that all well ID's start with a capital letter
unmatched_wells <- merge(all_data[, .(date, Plate.SN, well, trt, dose, file.name)], wllq_info[, .(date, Plate.SN, well, wllq, wllq_notes)], all.y = T)[is.na(file.name)]
if(nrow(unmatched_wells) > 0) {
  cat("The following rows from wells_with_well_quality_zero.csv did not match any rows in all_data:\n")
  print(unmatched_wells)
  cat("\nSummary of all_data:\n")
  print(all_data[, .(plates = paste0(sort(unique(Plate.SN)),collapse=",")), by = "date"][order(date)])
  warning("Update wells_with_well_quality_zero.csv")
}

# check for any DIV that do not match all_data
unmatched_DIV <- merge(all_data[, .(date, Plate.SN, well, trt, dose, file.name, DIV)], wllq_info[!is.na(suppressWarnings(as.numeric(DIV))), .(date, Plate.SN, well, DIV = as.numeric(DIV), wllq, wllq_notes)], 
      all.y = T)[is.na(file.name)]
if(nrow(unmatched_DIV) > 0) {
  cat("The following wells and DIV rows from wells_with_well_quality_zero.csv did not match the DIV for the corresponding plates in all_data:\n")
  print(unmatched_DIV)
  cat("\nSummary of all_data:\n")
  print(all_data[Plate.SN %in% unmatched_DIV$Plate.SN, .(DIVs = paste0(sort(unique(DIV)),collapse=",")), by = c("date","Plate.SN")][order(date,Plate.SN)])
  warning("Update wells_with_well_quality_zero.csv")
}

# make sure there is just 1 entry for every date-plate-well-affected_endpoints
wllq_info <- wllq_info[, .(wllq = min(wllq), wllq_notes = paste0(unique(wllq_notes), collapse= '; ')), by = .(date, Plate.SN, DIV, well, affected_endpoints)]

# initializing values
all_data[, `:=`(wllq = 1, wllq_notes = "")] 

# Update wllq for the wells where wllq==0 for all DIV
all_data <- merge(all_data, wllq_info[grepl("mea",affected_endpoints) & DIV == "all", .(date, Plate.SN, well, wllq, wllq_notes)], by = c("date","Plate.SN","well"),
                  suffixes = c("",".wllq_update"), all.x = TRUE)
all_data[!is.na(wllq.wllq_update), `:=`(wllq = wllq.wllq_update, wllq_notes = paste0(wllq_notes.wllq_update, "; "))]
all_data <- all_data[, .SD, .SDcols = names(all_data)[!grepl("wllq_update",names(all_data))]]

# next, remove any data where wllq==0 for specific DIV
wllq_info[, DIV := suppressWarnings(as.numeric(DIV))]
all_data[, full_id := paste(date, Plate.SN, well, DIV, sep = "_")] # hopefully I will find a better way to do this in the future
wllq_info[, full_id := paste(date, Plate.SN, well, DIV, sep = "_")]
cat("The following data rows will be REMOVED from the AUC analysis because wllq==0 and DIV is not 'all':\n")
# (Only removing rows where the wllq is not already 0 because of soem affect on all DIV)
print(merge(all_data[wllq == 1, .(date, Plate.SN, well, DIV, full_id, trt, dose, file.name, meanfiringrate)], 
            wllq_info[grepl("mea",affected_endpoints) & !is.na(DIV) & wllq == 0], 
            by = c("date","Plate.SN","well","DIV","full_id")))
all_data <- all_data[!(wllq == 1 & full_id %in% wllq_info[grepl("mea",affected_endpoints) & !is.na(DIV) & wllq == 0, unique(full_id)])]
all_data[, full_id := NULL]
# setting wllq to 0 instead
# all_data <- merge(all_data, wllq_info[grepl("mea",affected_endpoints) & !is.na(DIV), .(date, Plate.SN, DIV, well, wllq, wllq_notes)], by = c("date","Plate.SN","DIV","well"),
#       suffixes = c("",".wllq_update"), all = TRUE)
# all_data[wllq.wllq_update == 0, `:=`(wllq = wllq.wllq_update, wllq_notes = paste0(wllq_notes, wllq_notes.wllq_update, "; "))]
# all_data <- all_data[, .SD, .SDcols = names(all_data)[!grepl("wllq_update",names(all_data))]]

# Label the wllq as "wllq_by_well", because further wllq adjustments may be made once merge in "well_quality_notes_per_culture_treatment_cndx.xlsx"
# (based on treatment and cdx) after have verified meta data 
setnames(all_data, old = 'wllq', new = 'wllq_by_well')
setnames(all_data, old = 'wllq_notes', new = 'wllq_notes_by_well')

# print summary of wllq updates
cat("Wllq update summary:\n")
print(all_data[, .(number_of_well_recordings = .N, wllq_by_well = paste0(sort(unique(wllq_by_well)),collapse=",")), by = c("wllq_notes_by_well")][order(-wllq_by_well), .(wllq_by_well, wllq_notes_by_well, number_of_well_recordings)])
print(all_data[wllq_by_well == 0, .(date, Plate.SN, well, trt, dose, DIV, wllq_by_well, wllq_notes_by_well, meanfiringrate, nAE)])
rm(wllq_info)


# check for any DIV other than use_DIVS -----------------------------------------------
cat("\nChecking for any DIV other than",use_divs,"...\n")
all_data$date_plate <- paste0(all_data$date, "_", all_data$Plate.SN)
allDIV <- unique(all_data$DIV)
diffDIV <- setdiff(allDIV, use_divs)

if (length(diffDIV) > 0) {
  plates_with_diff_div <- unique(all_data[DIV %in% diffDIV, Plate.SN])
  cat(paste0("Recordings from DIV ",paste0(diffDIV,collapse = ","), " are found in ",paste0(plates_with_diff_div,collapse=",")),"\n")
  
  if(!interpolate_diff_divs) {
    warning(paste0("\nData from ",paste0(diffDIV,collapse = ",")," will be excluded."))
    all_data <- subset(all_data, DIV %in% use_divs)
  }
  else {
    # Interpolate the standard DIV values from the existing reocrdings
    update_plates <- all_data[DIV %in% diffDIV, unique(date_plate)]
    
    for (date_platei in update_plates) {
      cat(date_platei,"\n")
      dat <- all_data[date_plate == date_platei]
      all_data <- all_data[date_plate != date_platei]
      
      plate.DIVs <- unique(dat$DIV)
      # remove.DIV <- setdiff(plate.DIVs, use_divs)
      add.DIVs <- setdiff(use_divs, plate.DIVs)
      
      for (add.DIV in add.DIVs) {
        dat <- linearInterpolateDIV(dat, new.DIV = add.DIV)
      }
      
      # add updated plate data back to all_data
      dat <- dat[DIV %in% use_divs] # remove the non-standard DIVs
      all_data <- rbind(all_data, dat)
      rm(dat)
    }
  }
}


# check if ea plate has a recording for each of use_divs ------------------------
# original method, assuming every plate had same missing DIVs
# date_plates <- unique(all_data$date_plate)
# for (date_platei in date_plates) {
#   plate_divs <- all_data[date_plate == date_platei, sort(unique(DIV))]
#   missing_divs <- setdiff(use_divs, plate_divs)
#   if (length(missing_divs) > 0) {
#     cat(paste0("There is no data for ",sub("_"," ",date_platei), " DIV ",paste0(missing_divs,collapse=","),"\n"))
#     
#     if (length(missing_divs) > 1) {
#       warning(paste0("No data will be used from ",date_platei))
#       all_data <- all_data[date_plate != date_platei]
#     }
#     else {
#       cat("Values will be estimated from corresponding wells in other plates in same culture.\n")
#       # Generate values for missing DIV by the median of other plates from this DIV
#       all_data <- estimate_missing_DIV(dat = all_data, date_platei, missing_divs)
#     }
#   }
# }
# all_data <- all_data[, date_plate := NULL]
cat("\nChecking that every plate has a recording on DIV",use_divs,"...\n")
all_data[, well_id := paste(date, Plate.SN, well, sep = "_")]
wells_missing_div <- all_data[, .(DIV_flag = ifelse(length(setdiff(use_divs, unique(DIV)))>0, 1, 0),
                                  missing_DIV = list(setdiff(use_divs, unique(DIV)))), by = c("date","Plate.SN","date_plate","well","well_id")][DIV_flag == 1]
if (nrow(wells_missing_div) == 0) {
  cat("Every well has data from DIVs",use_divs,"\n")
} else {
  check_plates <- wells_missing_div[, unique(date_plate)]
  for(date_platei in check_plates) {
    missing_divs <- wells_missing_div[date_plate == date_platei, unique(unlist(missing_DIV))]
    cat(date_platei,'is missing some recordings.\n')
    
    # loop through each missing_div on this plate (will check if multiple DIV estimated afterwards)
    # That way I still have estimated values, even if wllq set to 0
    for (add.DIV in missing_divs) {
      # Generate values for missing DIV by the median of other plates from this DIV
      all_data <- estimate_missing_DIV(dat = all_data, date_platei, add.DIV)
    }
  }
  
  # If more than 1 DIV value had to be estimated for a given well, set wllq==0
  cat("\nSetting wllq to 0 for plates with estimated values for more than 1 DIV...\n")
  check_multiple_missing <- function(wllq_notes_vector) {
    if (sum(grepl("estimated as median from corresponding wells",wllq_notes_vector)) > 1)
      paste0(wllq_notes_vector, "Multiple recordings missing; ")
    else
      # no changes needed for wllqnotes
      wllq_notes_vector
  }
  all_data[, wllq_notes_by_well := lapply(.SD, check_multiple_missing), .SDcols = "wllq_notes_by_well", by = "well_id"]
  all_data[grepl("Multiple recordings missing",wllq_notes_by_well), wllq_by_well := 0] # could I merge this step with above?
  cat(all_data[grepl("Multiple recordings missing",wllq_notes_by_well), length(unique(Plate.SN))],'plates affected.\n')
}
all_data <- all_data[, c("date_plate","well_id") := NULL]


# Removing BIC data
#bis_rows <- grep("12_01_", all_data$file.name, fixed=TRUE) #index all Bicuculline-treated wells
#all_data <- all_data[- bis_rows,] #Remove all bic-treated wells

# rename a few columns for readability and consistency
setnames(x = all_data, old = c("trt","Mutual.Information"), new = c("treatment","mi"))

# going back to a data frame, for compatiblity with the functions below
all_data <- as.data.frame(all_data)
# Convert doses to character to ensure that no minor discrepancies in the decimal values of the doses from different files
# interferes with the groupings in the "split" function below
# all_data$dose <- sprintf("%.5f", all_data$dose) # old method - prints 5 decimal places, so may lose some significant figures
all_data$dose <- as.character(signif(all_data$dose, 4)) # 1 more sig fig than is used in TCPL, to ensure more than enough precision

# save a snapshot of the combined prepared data, with the added/interpolated rows where DIV where missing
write.csv(all_data, file = file.path(basepath,"output",paste0(project_name,"_parameters_by_DIV.csv")), row.names = FALSE)

#Replace all NAs with zeros for AUC calculations - This may be undesirable for MEA parameters that are derived from other parameters.
all_data[is.na(all_data)] <- 0

## Split data frame by individual wells over time (interaction term speeds this up greatly for larger datasets)
# all_data_split <- split(all_data, by = c("date","Plate.SN","well","trt","dose"), drop = TRUE, sorted = TRUE) # want to verify that this is identical in the future
all_data_split <- split(all_data, interaction(all_data$date, all_data$Plate.SN, all_data$well, all_data$treatment, all_data$dose, drop=TRUE)) # Split data into bins of single well across time (DIV
if(sum(unlist(lapply(all_data_split, nrow)) != 4) > 0) {
  stop("Some chunks in all_data_split do not have exactly 4 data rows (1 from each DIV)\n(Remove or modify this error if not testing exactly 4 DIV)")
}

#*****************************************************************************
#*                              FUNCTION                                     *
#*****************************************************************************

## Function to calculate area under the curve (AUC) for each ontogeny parameter for each bin (each experiment)
calc_auc <- function(all_data_split, sqrt=FALSE) {
  require(pracma) #pracma package has trapz function that computes AUC based on trapezoidal geometry (no curve fitting)
  
  endpoint_cols <- c('meanfiringrate','burst.per.min','mean.isis','per.spikes.in.burst','mean.dur','mean.IBIs','nAE','nABE',
                     'ns.n','ns.peak.m','ns.durn.m','ns.percent.of.spikes.in.ns','ns.mean.insis','ns.durn.sd','ns.mean.spikes.in.ns','r','mi')
  
  out <- lapply(1:length(all_data_split), function(i) {
    
    all_data_split[[i]] <- all_data_split[[i]][order(all_data_split[[i]][,"DIV"]),]  # Make sure order of rows follows DIV time
    
    date <- all_data_split[[i]]$date[1]
    Plate.SN <- all_data_split[[i]]$Plate.SN[1]
    well <- all_data_split[[i]]$well[1]
    treatment <- all_data_split[[i]]$treatment[1]
    dose <- all_data_split[[i]]$dose[1]
    units <- all_data_split[[i]]$units[1]
    wllq_by_well <- min(all_data_split[[i]]$wllq_by_well) # if any DIV recording for this well still has wllq_by_well==0, don't include that well
    wllq_notes_by_well <- paste0(unique(all_data_split[[i]]$wllq_notes_by_well),collapse="")
    
    for (j in endpoint_cols) {
      param_name <- paste(j, "_auc", sep="") # create auc variable name
      assign(param_name, round(trapz(append(all_data_split[[i]][,"DIV"], 2, after=0), append(all_data_split[[i]][,j], 0, after=0)),6), inherits=TRUE) # calculate auc, assign to variable name
    }
    
    # put vector of AUC values together
    c("date" = date, "Plate.SN" = as.character(Plate.SN), "well" = as.character(well), "treatment" = as.character(treatment), "dose" = dose, "units" = as.character(units), 
      "wllq_by_well" = wllq_by_well, "wllq_notes_by_well" = wllq_notes_by_well, sapply(paste(endpoint_cols,"auc",sep = "_"), get, USE.NAMES = T))
  })
  
  ##FIX!!! - (amy): not sure what needs to be fixed here
  sum_table <- as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE) # Re-form data frame
  # setnames(sum_table, old = paste0(assay_component_map$create_burst_ont_Data_endpoint,"_auc"), new = assay_component_map$tcpl_acsn)
  sum_table[,paste0(endpoint_cols,"_auc")] <- lapply(sum_table[,paste0(endpoint_cols,"_auc")], as.numeric)
  #sum_table[sum_table[,"dose"]==sprintf("%.5f",0),"treatment"] <- ControlTreatmentName # control treatment name will be updated in tcpl_MEA_dev_AUC
  
  # I don't think this is needed anymore (Amy 8/17/2020)
  if (sqrt==TRUE){
    sum_table <- cbind(sum_table[,1:6], sqrt(sum_table[,7:25]))
  }
  
  sum_table
}

#*****************************************************************************
#*                             END FUNCTIONS                                 *
#*****************************************************************************

sum_table <- calc_auc(all_data_split = all_data_split)

write.csv(sum_table, file.path(basepath, "output", filename), row.names = FALSE)

rm(list = c("all_data","all_data_split","parameter_data","mi_data","sum_table"))
cat("\n", filename," is ready\n",sep = "")
