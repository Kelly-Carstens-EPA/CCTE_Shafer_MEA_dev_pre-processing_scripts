tcpl_MEA_dev_AUC <- function(basepath, project_name, 
                             AUCsourcefilename = file.path(basepath, "output", paste0(project_name, "_AUC.csv")), 
                             DIVsourcefilename = file.path(basepath, "output", paste0(project_name, "_parameters_by_DIV.csv")), 
                             cytotox_filename = file.path(basepath, "output", paste0(project_name, "_cytotox.csv")),
                             assay_component_map_filename = file.path(dirname(basepath), "mea_nfa_component_name_map.csv"))
{
  
  require(data.table)
  
  # get DIV data and melt
  DIV_data <- fread(DIVsourcefilename)
  idcols <- c("date","Plate.SN","well","treatment","dose","units","wllq_by_well","wllq_notes_by_well")
  endpoint_cols <- setdiff(names(DIV_data),c(idcols,"DIV","file.name"))
  DIV_data[, (endpoint_cols) := lapply(.SD, as.numeric), .SDcols = endpoint_cols]
  DIV_data <- melt(DIV_data, id.vars = c(idcols,"DIV"), measure.vars = endpoint_cols, variable.name = "src_acsn",
                value.name = "rval", variable.factor = FALSE)
  DIV_data[, `:=`(src_acsn = paste0(src_acsn,"_DIV",DIV),
                  srcf = basename(DIVsourcefilename))]
  DIV_data[, DIV := NULL]
  
  # get AUC data and melt
  AUC <- fread(AUCsourcefilename)
  AUC <- melt(AUC, id.vars = idcols, measure.vars = paste0(endpoint_cols,"_auc"), variable.name = "src_acsn",
                value.name = "rval", variable.factor = FALSE)
  AUC[, srcf := basename(AUCsourcefilename)]
  
  # rbind DIV and AUC data
  longdat <- rbind(DIV_data, AUC)
  longdat[, `:=`(coli = as.numeric(sub("[[:alpha:]]","",well)), rowi = match(sub("[[:digit:]]","",well), LETTERS))]
  longdat[, c("well") := list(NULL)]
  setnames(longdat, old = "dose", new = "conc")
  rm(list = c("DIV_data","AUC"))
  
  # add cytotox data
  cytotox_data <- fread(cytotox_filename)
  if(length(setdiff(names(longdat),names(cytotox_data))) > 0) cat('cytotox data does not have ',paste0(setdiff(names(longdat),names(cytotox_data)),collapse=","),'. Will fill with NA\n', sep ='')
  if(length(setdiff(names(cytotox_data),names(longdat))) > 0) cat('MEA data does not have ',paste0(setdiff(names(cytotox_data),names(longdat)),collapse=","),'. Will fill with NA\n', sep = '')
  longdat <- rbind(longdat, cytotox_data, fill = T)
  rm(list = c("cytotox_data"))
  longdat[, treatment := as.character(treatment)] # sometimes the treatment is read as an integer instead of a char
  
  # replace the src_acsn with the TCPL acsn
  assay_component_map <- as.data.table(read.csv(assay_component_map_filename, stringsAsFactors = FALSE))
  longdat <- merge(longdat, assay_component_map, by = c("src_acsn"), all.x = T)
  if (any(is.na(unique(longdat$acsn)))) {
    print(longdat[is.na(acsn), unique(src_acsn)])
    stop(paste0("The above src_acsn's are not found in ",assay_component_map_filename))
  }
  
  # Define apid
  longdat[, apid := paste(date, Plate.SN, sep = "_")]
  
  # Assign wllt
  longdat[, wllt := "t"]
  longdat[ conc == 0, wllt := "n"]
  
  # get the desired columns, in the desired order
  longdat <- longdat[, .SD, .SDcols = intersect(c('apid', 'rowi', 'coli', 'treatment', 'conc', 'wllq_by_well', 'wllq_notes_by_well', 'wllt', 'rval', 'acsn', 'srcf', 'units'), names(longdat))]
  # longdat may or may not include "units"
  
  cat("long-format data is ready.\n")
  return(longdat)
}


# additional functions to prepare the data for TCPL mc0 format

update_control_well_treatment <- function(dat, control_compound, culture_date = c(), plates = c(), control_rowi) {
  apids <- Reduce(f = union, x = lapply(culture_date, function(x) grep(x, unique(dat$apid), val = T)))
  if (length(plates) > 0) {
    apids <- Reduce(f = union, x = lapply(plates, function(x) grep(x, apids, val = T)))
  }
  cat("Control treatment will be updated to",control_compound,"for the following wells:\n")
  print(dat[wllt == "n" & apid %in% apids & rowi %in% control_rowi, unique(.SD), .SDcols = c("apid","treatment","rowi","coli","conc")][order(apid,rowi,coli)])
  dat[wllt == "n" & apid %in% apids & rowi %in% control_rowi, treatment := control_compound]
}

update_treatment_names <- function(date, root_output_dir, project_name) {
  trt_name_map <- as.data.table(read.csv(file.path(root_output_dir, "supplemental_mea_treatment_name_map.csv"), stringsAsFactors = F))
  trt_name_map <- trt_name_map[dataset == project_name, .(mea_treatment_name, updated_treatment_name)]
  unused_trt_names <- setdiff(unique(trt_name_map$mea_treatment_name), unique(dat$treatment))
  if(length(unused_trt_names)> 0 ){
    cat("Some expected mea treatment names in 'supplemental_mea_treatment_name_map.csv' are not in the input data table:", unused_trt_names, "\n", sep = "\n")
  }
  dat <- merge(dat, trt_name_map, by.x = "treatment", by.y = "mea_treatment_name", all.x = T)
  dat[is.na(updated_treatment_name), updated_treatment_name := treatment] # for compound names that do not need to be updated
  setnames(dat, old = c("treatment","updated_treatment_name"), new = c("mea_treatment_name","treatment"))
  return(dat)
}

check_and_assign_spids <- function(dat, spidmap) {
  if (length(setdiff(c("treatment","spid"), names(spidmap))) > 0) {
    stop("The following columns are not found in spidmap: ",paste0(setdiff(c("treatment","spid"), names(spidmap)),collapse =","))
  }
  if (spidmap[!is.na(spid) & treatment %in% unique(dat$treatment), .(length(unique(spid))), by = "treatment"][,any(V1 !=1)]) {
    stop(paste0("The following treatments map to multiple spids: ",
                spidmap[!is.na(spid) & treatment %in% unique(dat$treatment), .(length(unique(spid))), by = "treatment"][V1 != 1,paste0(treatment,collapse=", ")]))
  }
  dat <- merge(dat, spidmap[, .(spid, treatment)], by = "treatment", all.x = TRUE)
  dat[wllt != 't' & is.na(spid), spid := treatment]
  if (dat[wllt == "t", any(is.na(spid))]) {
    cat("The following treatments don't have a corresponding spid in the spidmap:\n")
    print(dat[wllt == "t" & is.na(spid), unique(treatment)])
    stop("Adjust input dat before continuing")
  }
  return(dat)
}
