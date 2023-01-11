update_wllq_by_well_DIV <- function(longdat, basepath = NULL, get_files_under_basepath = TRUE) {
  
  cat("\nUpdating wllq...\n")
  
  # Get wllq file
  if (get_files_under_basepath) {
    
    wllq.info.file <- file.path(basepath, "well_quality_notes_by_well_DIV.xlsx")
    
    # For backwards compatibility:
    if(!file.exists(wllq.info.file)) wllq.info.file <- file.path(basepath, "wells_with_well_quality_zero.csv")
  }
  else {
    wllq.info.file <- choose.files(default = basepath, caption = "Select well quality csv table")
  }
  
  # Read file
  if(grepl('\\.xlsx',wllq.info.file)) {
    wllq_info <- as.data.table(openxlsx::read.xlsx(wllq.info.file, sheet = 1))
    wllq_info[, c(setdiff(names(wllq_info),c('DIV','wllq'))) := lapply(.SD, as.character), .SDcols = c(setdiff(names(wllq_info),c('DIV','wllq')))]
    wllq_info[, c('DIV','wllq') := lapply(.SD, as.numeric), .SDcols =  c('DIV','wllq')]
  }
  else {
    # (for backwards compatibility)
    wllq_info <- as.data.table(read.csv(wllq.info.file, 
                                        colClasses = c(rep("character",4),"numeric",rep("character",2)), stringsAsFactors = FALSE))
  }
  
  # add rowi and coli to wllq_info
  wllq_info[, well := toupper(well)] # ensure that all well ID's start with a capital letter
  wllq_info[, `:=`(coli = as.numeric(stringi::stri_extract(well, regex = '[0-9]+')), 
                   rowi = match(stringi::stri_extract(well, regex = '[A-Za-z]+'), LETTERS))]
  
  # expand the rows where well == "all" to include every well in plate
  # (Currently assuming that plates are 6x8, but would like to make this more flexible in the future)
  for(table_row in which(wllq_info$well == "all")) {
    full_plate <- wllq_info[table_row, .(date, Plate.SN, DIV, rowi = rep(1:6, each = 8), coli = rep(1:8, times = 6), wllq, wllq_notes, affected_endpoints)]
    wllq_info <- rbind(wllq_info, full_plate) # rbind the new rows to the end of wllq_info
  }
  wllq_info <- wllq_info[well != "all"]
  
  
  # check for any DIV that do not match data
  unmatched_DIV <- merge(longdat[, .(date, Plate.SN, well, trt, dose, file.name, DIV)], wllq_info[!is.na(suppressWarnings(as.numeric(DIV))), .(date, Plate.SN, well, DIV = as.numeric(DIV), wllq, wllq_notes)], 
                         all.y = T)[is.na(file.name)]
  if(nrow(unmatched_DIV) > 0) {
    cat("The following wells and DIV rows from wells_with_well_quality_zero.csv did not match the DIV for the corresponding plates in longdat:\n")
    print(unmatched_DIV)
    cat("\nSummary of longdat:\n")
    print(longdat[Plate.SN %in% unmatched_DIV$Plate.SN, .(DIVs = paste0(sort(unique(DIV)),collapse=",")), by = c("date","Plate.SN")][order(date,Plate.SN)])
    warning("Update wells_with_well_quality_zero.csv")
  }
  
  # Expand rows where DIV == "all"
  # Where DIV  == "all", get the DIV present in data associated with the given date, Plate, well, etc.
  # (I could just assume DIV 5, 7, 9, and 12, but i want to be more flexible in case that changes)
  wllq_info.div.all <- wllq_info[grepl("mea",affected_endpoints) & DIV == "all", .(date, Plate.SN, rowi, coli, wllq, wllq_notes)]
  setkey(longdat, date, Plate.SN, rowi, coli, wllq, wllq_notes)
  divs.present <- longdat[J(wllq_info.div.all), unique(.SD), .SDcols = c('date', 'Plate.SN', 'rowi','coli', 'DIV')]
  wllq_info.div.all <- merge(wllq_info.div.all, divs.present, by = c('date', 'Plate.SN', 'well', 'rowi','coli'))
  wllq_info <- wllq_info[DIV != "all"]
  wllq_info <- rbind(wllq_info, wllq_info.div.all)
  

  # checking for typo's or missing "affected_endpoints"
  if(nrow(wllq_info[!grepl("(mea)|(CTB)|(LDH)",affected_endpoints)])>0) {
    cat("The following rows don't match any of the expected affected_endpoints (mea,CTB,LDH):\n")
    print(wllq_info[!grepl("(mea)|(CTB)|(LDH)",affected_endpoints)])
    stop(paste0("Update ",wllq.info.file))
  }
  
  # Expand rows for each "affected endpoint"
  affected_endpoints_list <- unique(unlist(stringi::stri_split(wllq_info$affected_endpoints, regex = '[,; ]+')))
  wllq_info <- rbindlist(lapply(affected_endpoints_list,
                                function(endpoint_typei) cbind(wllq_info[grepl(endpoint_typei,affected_endpoints), 
                                                                         .SD, .SDcols = setdiff(names(wllq_info),'affected_endpoints')],
                                                               'endpoint_type' = endpoint_typei)))
  
  # Get the minimum wllq and collapse across notes for every date, plate, well, endpoint_type combination
  wllq_info <- wllq_info[, .(wllq = min(wllq),
                             wllq_notes = paste0(unique(wllq_notes), collapse = "; ")),
                         by = .(date, Plate.SN, DIV, rowi, coli, wllq_ref, endpoint_type)]

  # Prepare columns in longdat
  longdat[, `:=`(date = as.character(date), rowi = as.numeric(rowi))]
  longdat[grepl('LDH',acsn), endpoint_type := 'LDH']
  longdat[grepl('AB',acsn), endpoint_type := 'CTB'] # this is the abbreviation I'm using in the wllq tables... could change to AB
  longdat[is.na(endpoint_type), endpoint_type := 'mea']
  
  # Check for any rows where the date, plate, rowi, coli, affected endpoints is not in longdat
  unmatched_wells <- merge(longdat[, .(date, Plate.SN, rowi, coli, treatment, conc, srcf, endpoint_type)], 
                           wllq_info[, .(date, Plate.SN, rowi, coli, wllq, wllq_notes, endpoint_type)], all.y = T)[is.na(srcf)]
  if(nrow(unmatched_wells) > 0) {
    cat(paste0("The following rows from ",wllq.info.file," did not match any rows in data:\n"))
    print(unmatched_wells)
    cat("\nSummary of longdata:\n")
    print(longdata[, .(plates = paste0(sort(unique(Plate.SN)),collapse=",")), by = "date"][order(date)])
    stop(paste0("Update ",wllq.info.file))
  }
  
  # Merge wllq_info to longdat, define wllq
  longdat <- merge(longdat, wllq_info, all.x = TRUE)
  longdat[is.na(wllq), `:=`(wllq = 1, wllq_notes = "")]
  
  # flag NA rval's where there is no wllq note
  na_rvals <- longdat[is.na(rval) & wllq == 1]
  if (nrow(na_rvals) > 0) {
    cat("The following rval's are missing in data, but wllq==1\n")
    print(na_rvals)
    resp <- readline(prompt = "Do you wish to set wllq=0 for all of these wells? (Only do this if you know that these values should be NA) (y/n): ")
    if (resp %in% c("y","Y","yes","Yes")) {
      longdat[is.na(rval) & wllq == 1, `:=`(wllq = 0, wllq_notes = "rval is NA; ")]
    }
    else {
      stop("Update wllq_info table")
    }
  }
  
  # Label the wllq as "wllq_by_well", because further wllq adjustments may be made once merge in "well_quality_notes_by_culture_treatment_cndx.xlsx"
  # (based on treatment and cndx) after have verified meta data 
  setnames(longdat, old = c('wllq','wllq_notes','wllq_ref'), new = paste0(c('wllq','wllq_notes','wllq_ref'),'_by_well'))
  
  # summary of wllq updates
  cat("Wllq summary:\n")
  print(longdat[, .N, by = c("src_acsn","wllq_by_well","wllq_notes_by_well")][order(src_acsn, wllq_by_well, wllq_notes_by_well)])
  
  return(longdat)
}


# Goal:
# - remove data points where is MEA and is a particular DIV (so can be linearly interpolated
# - add wllq notes and wllq to points where have a note

# question:
# Why even apply the wllq at this stage, if it isn't final until after the run_me?
# I could just apply the wllq to the DIV situation...
# might not be a bad idea...


# Game plan:
# - in burst parameter to AUC, just apply the wllq to remove certian DIV
# - in run_me, apply the wllq_notes to everything.
# -> where there was a note aobut a particular missing DIV, add that in as a wllq note (even though as now been interpolated)

# RESU MEHERE!

# Next steps:
# - take what you need from below into burst_parameter_to_AUC to just remove certain DIV (I think below is a workign draft)
# - check if wllq is referenced in other burst paramter to AUC body or sub-functions
# - edit above function to add wllq note where certain DIV were removed in burst parameter to AUC. But don't set wllq to 0, because the idea is that I already interpolated a new DIV...? Hmm, could get trickey, think about this
# -> document in mea nfa notebook that you've made this update!!
# - Remove wllq update from cytotox prep, and any reference to wllq in tcpl_prep function?
# - Finalize this funciton, do a test run, then and reference in run_me
# - Implement addition of new secondary wllq table to run_me as well, possibly borrow from this function.
# - do a minim test-run from start to fin with wllq


# Get wllq file
if (get_files_under_basepath) {
  
  wllq.info.file <- file.path(basepath, "well_quality_notes_by_well_DIV.xlsx")
  
  # For backwards compatibility:
  if(!file.exists(wllq.info.file)) wllq.info.file <- file.path(basepath, "wells_with_well_quality_zero.csv")
} else {
  wllq.info.file <- choose.files(default = basepath, caption = "Select well quality csv table")
}

# Read file
if(grepl('\\.xlsx',wllq.info.file)) {
  wllq_info <- as.data.table(openxlsx::read.xlsx(wllq.info.file, sheet = 1))
  wllq_info[, c(setdiff(names(wllq_info),c('DIV','wllq'))) := lapply(.SD, as.character), .SDcols = c(setdiff(names(wllq_info),c('DIV','wllq')))]
  wllq_info[, c('DIV','wllq') := lapply(.SD, as.numeric), .SDcols =  c('DIV','wllq')]
} else {
  # (for backwards compatibility)
  wllq_info <- as.data.table(read.csv(wllq.info.file, 
                                      colClasses = c(rep("character",4),"numeric",rep("character",2)), stringsAsFactors = FALSE))
}

# Just concerned about MEA data here
wllq_info <- wllq_info[grepl('mea',tolower(affected_endpoints))]

# expand the rows where well == "all" to include every well in plate
# (Currently assuming that plates are 6x8, but would like to make this more flexible in the future)
# expand the rows where well == "all" to include every well in plate
for(table_row in which(wllq_info$well == "all")) {
  full_plate <- wllq_info[table_row, .(date, Plate.SN, DIV, well = paste0(unlist(lapply(LETTERS[1:6],rep,times=8)),rep(1:8,6)), wllq, wllq_notes, affected_endpoints)]
  wllq_info <- rbind(wllq_info, full_plate) # rbind the new rows to the end of wllq_info
}
wllq_info <- wllq_info[well != "all"]

# Remove rows where DIV == "all" (this info will be merged in later)
wllq_info <- wllq_info[DIV != 'all']
wllq_info[, DIV := as.numeric(DIV)]
stopifnot(nrow(wllq_info[is.na(DIV)]) == 0)

# check for any DIV that do not match data
unmatched_DIV <- merge(longdat[, .(date, Plate.SN, well, trt, dose, file.name, DIV)], wllq_info[, .(date, Plate.SN, well, DIV, wllq, wllq_notes)], 
                       all.y = T)[is.na(file.name)]
if(nrow(unmatched_DIV) > 0) {
  cat("The following wells and DIV rows from wells_with_well_quality_zero.csv did not match the DIV for the corresponding plates in longdat:\n")
  print(unmatched_DIV)
  cat("\nSummary of longdat:\n")
  print(longdat[Plate.SN %in% unmatched_DIV$Plate.SN, .(DIVs = paste0(sort(unique(DIV)),collapse=",")), by = c("date","Plate.SN")][order(date,Plate.SN)])
  warning(paste0("Update ",wllq.info.file))
}

# Get indicies of MEA data where wllq==0 for specific DIV
well.div.to.remove <- wllq_info[wllq == 0, unique(.SD), .SDcols = c('date','Plate.SN','well','DIV')]
# If multiple DIV affected for a given well... then we'd want to just set the wllq to 0, not remove
stopifnot(nrow(well.div.to.remove[, .N, by = .(date, Plate.SN, well)][N > 1]) == 0)

# Merge with longdat
longdat[, in_dat := 1]
well.div.to.remove[, in_well.div.to.remove := 1]
longdat <- merge(longdat, well.div.to.remove, by = c('date','Plate.SN','well','DIV'), all = T)

# Check for any date, Plate.SN, well, DIV combinations in wllq_info not in longdat (likley indicates typo in wllq_info)
if(nrow(longdat[is.na(in_dat)]) > 0) {
  print(longdat[is.na(in_dat), .N, by = .(date, Plate.SN, well, DIV, in_dat, in_well.div.to.remove)])
  stop('Above date, plate, well, DIV combinations in wllq_info are not in data')
}

# Remove where wllq == 0 for specific DIV
cat("The following data rows will be REMOVED from the AUC analysis because wllq==0 and DIV is not 'all':\n")
print(longdat[!is.na(in_well.div.to.remove), .(date, Plate.SN, well, DIV, full_id, trt, dose, file.name, meanfiringrate)])
longdat <- longdat[!is.na(in_well.div.to.remove)]

# Note: 
# Previous version: Only removed rows where the wllq is not already 0 because of soem affect on all DIV
# Now: removing rows where wllq is 0 for an individual DIV, even if there is another wllq note affecting all DIV

# Remove added columns
longdat[, c('in_dat_','in_well.div.to.remove') := NULL]
