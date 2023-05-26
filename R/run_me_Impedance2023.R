rm(list=ls()) # clear environment
graphics.off() # clear plot history
# ------------------------------------------------------------------------ #
# USER INPUT
# ------------------------------------------------------------------------ #
dataset_title <- "Impedance2023" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- FALSE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName <- "DMSO" # all compounds other than those listed below should have this vehicle control
different_vehicleControlCompounds = c("Deltamethrin","Glyphosate","Lindane","Tributyltin Chloride", "Rotenone", "Thiamethoxam") # e.g. c("Sodium Orthovanadate", "Amphetamine")
# Enter the names of the vehicle controls as they correspond to the compounds in the previous list
different_vehicleControls = c("EtOH:DMSO","H2O","EtOH", "DMSO", "DMSO", "DMSO") # e.g. c("Water", "Water")

spidmap_file <- ""
spid_sheet <- ""

project.dir <- "L:/Lab/NHEERL_MEA/Project - Impedance/Impedance Paper/prep for tcpl pipeline/NFA/Spike list" # project main folder (where will look for README files)
scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R"
root_output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the dataset_title folder will be created

update_concs_without_prompt <- FALSE
get_new_files_log <- FALSE
# ------------------------------------------------------------------------ #
# END USER INPUT
# ------------------------------------------------------------------------ #

library(data.table)
library(openxlsx)
library(RMySQL)


# Run main steps ----------------------------------------------------------

# create a summary log file and store the
if(save_notes_graphs) {
  sink(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_run_log_",as.character.Date(Sys.Date()),".txt")))
  cat("Output from the script run_me_",dataset_title,".R\n",sep="")
  cat("Date Ran:",as.character.Date(Sys.Date()),"\n")
  cat(R.version.string,"\n")
  cat("USER INPUT settings:\n")
  print(sapply(ls(), get, envir = .GlobalEnv))
  graphics.off()
  pdf(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")), width = 10, height = 8)
}


# > Scan for readme's that might affect dosing, wllq ----------------------

txt.files <- list.files(project.dir, pattern = '\\.txt', recursive = T, full.names = T)
readmes <- txt.files[grepl('read( )*me',tolower(txt.files))]
for (readme in readmes) {
  cat(dirname(readme),'\n')
  cat(scan(readme, what = character(), sep = '\n', quiet = T), sep = '\n')
  cat('\n')
}


# > Automated way to get files (optional) --------------------------------------------
source(file.path(scripts.dir, 'gather_files_functions.R'))

#modified to get neurostats files 
get_NFA_standard_structue <- function(culture.folderi) {
  cyto.files <- list.files(path = culture.folderi, pattern = '(Calculations)|(Summary)', full.names = T)
  plate.dirs <- list.files(path = culture.folderi, pattern = '^(MW)*[0-9\\-]{7}$', full.names = T)
  mfiles <- list.files(path = culture.folderi, pattern = 'MaestroExperimentLog', full.names = T, recursive = T)
  slists <- list.files(path = culture.folderi, pattern = '_spike_list\\.csv', full.names = T, recursive = T)
  neurostats <- list.files(path = culture.folderi, pattern = '_Compiler\\.csv', full.names = T, recursive = T )
  filesi <- c(cyto.files, mfiles, slists, neurostats)
  return(filesi)
}

culture.folders <- list.files(path = project.dir, pattern = "[0-9]{8}", full.names = T, recursive = F)
all.files <- c(culture.folders)
for (culture.folderi in culture.folders) {
  all.files <- c(all.files, get_NFA_standard_structue(culture.folderi))

  
  
  
  }

#Basic cleaning
all.files <- Filter(function(filei) !grepl('\\Q~$\\E',filei), all.files)
all.files <- Filter(function(filei) !grepl('deprecated',tolower(filei)), all.files)
length(all.files)
# should be # groups * (1 Calc file + 3 plates * (4 DIV + 1 maestro exp log file))

# Write to log file
writeLogFile(all.files, output.dir = file.path(root_output_dir, dataset_title), dataset_title, files_type = '')
# Writing 288 files to DNT_NTP2021_files_log_2022-04-19.txt ...
# [1] "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/DNT_NTP2021/DNT_NTP2021_files_log_2022-04-19.txt is ready."
#rm(all.files, culture.folders, readmes)


# > run the main steps ----------------------------------------------------

# Can enter "a" to select files, then hit cancel, to see summary of files selected
source(file.path(scripts.dir, 'source_steps.R'))

##neurostats files added in during first step
#created h5 files, calculated mfr..., calculated mutual information
## need to add in impedance data
neurostats <- Filter(function(filei) !grepl('*spike_list',filei), all.files)
neurostats <- Filter(function(filei) !grepl('*quality',filei), neurostats)
neurostats <- Filter(function(filei) !grepl('*Log',filei), neurostats)
neurostats <- Filter(function(filei) !grepl('*Calculation',filei), neurostats)
# can use function from acute template to extract impedance data from neurostatistic compiler
#Switched script path to acute path, remember to switch it back later
##sourcing every acute scripts because I don't want to go in and check which one to use
start.dir <- "L:/Lab/NHEERL_MEA/Project - Impedance/Impedance Paper/prep for tcpl pipeline/NFA/Impedance"
scripts.dir <-"L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_acute_for_tcpl/mea-acute-neural-stats-to-mc0-scripts"
scripts <- list.files(path = "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_acute_for_tcpl/mea-acute-neural-stats-to-mc0-scripts", pattern = "\\.R$", full.names = T, recursive = F)
sapply(scripts, source)
# read in the data from well.averages.rowi index
#for (i in neurostats) {
#combine multiple DIVs together in a single table
####imp2 is working perfectly## order in imp might be wrong, if wanna write a loop please use imp2 steps
groupbyplate <- list(1:6, 7:12, 19:24, 25:30, 31:36)
#there is a corrupted file DIV 9 for plate MW507214B, will skip over that single file
for (j in groupbyplate) {
imp.plate <- data.frame()
for (i in j) {
  file_scan <- scan(file = neurostats[i], what = character(), sep = "\n", blank.lines.skip = F, quiet=T) # empty lines will be just ""
  file_col1 <- sapply(file_scan, function(x) strsplit(x, split = ",")[[1]][1], USE.NAMES = F) # empty lines will be NA
  file_col2 <- sapply(file_scan, function(x) strsplit(x, split = ",")[[1]][2], USE.NAMES = F) # if nothing in second col, will be NA
  # get the index of the tag phrase 'Well Averages'
  well.averages.rowi <- grep("[Ww]ell [Aa]verages", file_col1)
   # find the next blank line after well averages index
  next.blank.row.dist <- which(is.na(file_col1[well.averages.rowi:length(file_col1)]) | file_col1[well.averages.rowi:length(file_col1)] == "")[1]
 imp <- as.data.table(read.table(neurostats[i], sep = ",", header = F, skip = (well.averages.rowi - 1), nrows = (next.blank.row.dist-1),
                                stringsAsFactors = F))
imp <- imp[c(1, nrow(imp)), 1:49] ##trim down in between, just need last and first row, no header 
#need to transpose into column/row, then change well name into rowi# coli#
imp <- transpose(imp)
imp[1,1] <- "well"#change "Well Averages" into "well"
imp[1,2] <- "wt.mean.resistance"
imp <- janitor::row_to_names(imp, 1)
#remove blank
imp$blank <- "12.107"
imp[ , 2] <- apply(imp[ , 2,drop=F], 2,function(x) as.numeric(as.character(x)))
imp[ , 3] <- apply(imp[ , 3,drop=F], 2,function(x) as.numeric(as.character(x)))
imp[, wt.mean.resistance := imp$wt.mean.resistance - imp$blank]
imp[, blank:= NULL]

#assigning date&PlateSN to match with treatment later
srcname <- basename(neurostats[i]) #"iNFA_20211103_MW502622B_00_(000)(000)_Compiler.csv"
srcname <- sub(" ", "_", srcname)
imp[, file.name := srcname]
srcname <- strsplit(srcname, split = "_")[[1]]
imp[, date := srcname[2]]
imp[, Plate.SN := srcname[3]]
imp[, DIV := as.numeric(srcname[4])]
##apid <- paste(srcname[2], "_", srcname[3])
#assign each well to treatment, then rearrange the columns
#find matching maestroexperimentallog
output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/Impedance2023"
masterChemFiles <- readLogFile(output_dir, files_type = "MaestroExperimentLog")
masterChemFile <- grep(pattern = paste0(srcname[2],"[_ ]"), masterChemFiles, value = T)
masterChemFile <- grep(pattern = paste0(srcname[3],"[_ ]"), masterChemFile, value = T)
#cant use chem.info.2 coz it is for h5 files
#assuming A1-F8 are arranged in order (should be a way of reordering before adding, but they should be in same order for both experimentlog and neurostats)
trt.info <- as.data.table(read.table(masterChemFile, sep = ",", header = T, skip = 0, nrows = 49, stringsAsFactors = F))
imp[, trt := trt.info$Treatment]
imp[, dose :=trt.info$Dose]
imp[, units :=trt.info$Units]
#rearrange columns to match NMI.csv
imp <- imp[,c(4, 5, 6, 1, 7, 8, 9, 3, 2)]
imp.plate <- rbind(imp.plate, imp)
#set the NAs to 0
imp.plate[is.na(imp.plate)] <- 0

#write files
imp.dir<-paste(output_dir, "/IMP",sep="")
dir.create(imp.dir)
setwd(imp.dir)
impfilename <- paste("IMP",srcname[2],srcname[3],sep="_")
write.csv(imp.plate, file = paste(impfilename, ".csv",sep =""), row.names= FALSE )
}
}

setwd(output_dir)

# >run AUC for imp data
#### read in imp data from IMP folder
#### then remove all functions that remove DIV 0 and 2
impfiles <- list.files(paste(file.path(root_output_dir, dataset_title), "/IMP", sep = ""))
impfiles <- paste(paste(file.path(root_output_dir, dataset_title), "/IMP/", sep = ""), impfiles[1], sep = "")
imp <- as.data.table(read.table(impfiles, sep = ",", header = T, skip = 0, stringsAsFactors = F))
imp <- imp[DIV %in% "9"]
impfiles <- list.files(paste(file.path(root_output_dir, dataset_title), "/IMP", sep = ""))
impfiles <- paste(paste(file.path(root_output_dir, dataset_title), "/IMP/", sep = ""), impfiles[2], sep = "")
imp2 <- as.data.table(read.table(impfiles, sep = ",", header = T, skip = 0, stringsAsFactors = F))
imp <- rbind(imp, imp2[DIV %in% "9"])
impfiles <- list.files(paste(file.path(root_output_dir, dataset_title), "/IMP", sep = ""))
impfiles <- paste(paste(file.path(root_output_dir, dataset_title), "/IMP/", sep = ""), impfiles[3], sep = "")
imp507 <- as.data.table(read.table(impfiles, sep = ",", header = T, skip = 0, stringsAsFactors = F))
#no need interpolate
#estimate missing DIV 9 for MW507214B based on values from 2 other plates in the same culture date 20211103_MW502622B,20211103_MW502622D
# calculate the median parameter value for each trt and dose (by trt, with all control wells grouped)
est <-data.frame(matrix(ncol = 4, nrow = 0))
colnames(est) <- c("well", "trt", "dose", "wt.mean.resistance")
treat <- c("Deltamethrin", "Thiamethoxam", "Lindane", "Tributyltin Chloride", "Rotenone", "Glyphosate")
sol <- c("DMSO", "DMSO", "EtOH:DMSO", "H2O", "EtOH", "DMSO")
sol.well <- c("A2", "B2", "C2", "D2", "E2", "F2")
dos <- c("0.03","0.1","0.3", "1","3","10","30")
for(i in 1:6) {
  
  for (j in 1:7) {
  
  est[nrow(est)+1, ] <- c(well = imp507[imp507$trt == treat[i] & imp507$dose == dos[j], ][1]$well, trt = treat[i], dose = dos[j], wt.mean.resistance = median(imp[imp$trt == treat[i] & imp$dose == dos[j] , ]$wt.mean.resistance))
  }
  }
for (k in 1:6) {
  est[nrow(est)+1, ] <- c(well = sol.well[k], trt = sol[k], dose = "0.00", wt.mean.resistance = median(imp[imp$trt == sol[k], ]$wt.mean.resistance))
}
est[ , 3] <- apply(est[ , 3,drop=F], 1,function(x) as.numeric(as.character(x)))
est$date <- imp507$date[1]
est$Plate.SN <- imp507$Plate.SN[1]
est$DIV <- "9"
est$units <- imp507$units[1]
est$file.name <- "median_at_DIV9_in_corresponding_wells_of_20211103_MW502622B,20211103_MW502622D"
est <- est[,c(5, 6, 7, 1, 2, 3, 8, 9, 4)]
imp507 <- rbind(imp507, est)
imp.dir<- paste(file.path(root_output_dir, dataset_title), "/IMP", sep = "")
setwd(imp.dir)
impfilename <- paste("IMP", imp507$date[1], imp507$Plate.SN[1] ,sep="_")
write.csv(imp507, file = paste(impfilename, ".csv",sep =""), row.names= FALSE )

#read in imp data and combine with all other parameters from parameters_by_DIV
##copied and modified from burst parameters to AUC r script
impfiles <- list.files(path = file.path(root_output_dir, dataset_title, "IMP"), pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
imp_data <- c()
for (i in 1:length(impfiles)) {
  imp_datai <- read.csv(impfiles[i], stringsAsFactors = F)
  imp_data <- rbind(imp_data, imp_datai)
  rm(imp_datai)
}
imp_data <- unique(imp_data) # eliminate duplicate data
#update well quality? -later
imp_data$dose <- as.character(signif(imp_data$dose, 4))
## Split data frame by individual wells over time (interaction term speeds this up greatly for larger datasets)
# all_data_split <- split(all_data, by = c("date","Plate.SN","well","trt","dose"), drop = TRUE, sorted = TRUE) # want to verify that this is identical in the future
imp_data_split <- split(imp_data, interaction(imp_data$Plate.SN, imp_data$well, imp_data$trt, imp_data$dose, drop=TRUE)) # Split data into bins of single well across
if(sum(unlist(lapply(imp_data_split, nrow)) != 6) > 0) {
  stop("Some chunks in imp_data_split do not have exactly 4 data rows (1 from each DIV)\n(Remove or modify this error if not testing exactly 4 DIV)")
}
#run AUC function
  require(pracma) #pracma package has trapz function that computes AUC based on trapezoidal geometry (no curve fitting)
  
  endpoint_cols <- "wt.mean.resistance"
  
  out <- lapply(1:length(imp_data_split), function(i) {

    imp_data_split[[i]] <- imp_data_split[[i]][order(imp_data_split[[i]][,"DIV"]),]  # Make sure order of rows follows DIV time
    
    date <- imp_data_split[[i]]$date[1]
    Plate.SN <- imp_data_split[[i]]$Plate.SN[1]
    well <- imp_data_split[[i]]$well[1]
    treatment <- imp_data_split[[i]]$trt[1]
    dose <- imp_data_split[[i]]$dose[1]
    units <- imp_data_split[[i]]$units[1]
   
       for (j in endpoint_cols) {
      param_name <- paste(j, "_auc", sep="") # create auc variable name
      #assign(param_name, round(trapz(append(imp_data_split[[i]][,"DIV"], 2, after=0), append(imp_data_split[[i]][,j], 0, after=0)),6), inherits=TRUE) # calculate auc, assign to variable name
      assign(param_name, round(trapz(imp_data_split[[i]][,"DIV"], imp_data_split[[i]][,j]),6), inherits=TRUE) # calculate auc, assign to variable name
    }
    
    # put vector of AUC values together
    c("date" = date, "Plate.SN" = as.character(Plate.SN), "well" = as.character(well), "treatment" = as.character(treatment), "dose" = dose, "units" = as.character(units), sapply(paste(endpoint_cols,"auc",sep = "_"), get, USE.NAMES = T))
  })
  
  ##FIX!!! - (amy): not sure what needs to be fixed here
  sum_table <- as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE) # Re-form data frame
  # setnames(sum_table, old = paste0(assay_component_map$create_burst_ont_Data_endpoint,"_auc"), new = assay_component_map$tcpl_acsn)
  #sum_table[,paste0(endpoint_cols,"_auc")] <- lapply(sum_table[,paste0(endpoint_cols,"_auc")], as.numeric)# no idea why this step convert everything to one number
  #sum_table[sum_table[,"dose"]==sprintf("%.5f",0),"treatment"] <- ControlTreatmentName # control treatment name will be updated in tcpl_MEA_dev_AUC

#save a copy of imp_auc  
sum_table$wllq_by_well <-"1"
sum_table$wllq_notes_by_well <- ""
sum_table[sum_table$Plate.SN == "MW507214B", ]$wllq_notes_by_well <- "DIV9 estimated as median from corresponding wells of 20211103_MW502622B,20211103_MW502622D;"
sum_table <- sum_table[,c(1, 2, 3, 4, 5, 6, 8, 9, 7)]
filename <- paste0(dataset_title, "_IMP_AUC.csv", sep = "")
write.csv(sum_table, file.path(root_output_dir, dataset_title, "output", filename), row.names = FALSE)

#save a copy of imp raw by DIV
colnames(imp_data)[5] <- "treatment" #change column nam to match para_mi_data
imp_data$wllq_by_well <-"1"
imp_data$wllq_notes_by_well <- ""
imp_data[imp_data$Plate.SN == "MW507214B" & imp_data$DIV == "9", ]$wllq_notes_by_well <- "DIV9 estimated as median from corresponding wells of 20211103_MW502622B,20211103_MW502622D;"
filename <- paste0(dataset_title, "_IMP_by_DIV.csv", sep = "")
write.csv(imp_data, file.path(root_output_dir, dataset_title, "output", filename), row.names = FALSE)

rm(imp_data,imp, imp2,imp_data_split, out, sum_table)

# > run tcpl_MEA_dev_AUC --------------------------------------------------
#source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
#dat <- tcpl_MEA_dev_AUC(basepath = file.path(root_output_dir,dataset_title), dataset_title)
# below is a modified copy of tcpl_MEA_dev_AUC

 AUCsourcefilename = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title, "_AUC.csv"))
 DIVsourcefilename = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title, "_parameters_by_DIV.csv")) 
 cytotox_filename = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title, "_cytotox.csv"))
 imp_rdata_filename =  file.path(root_output_dir, dataset_title, "output", paste0(dataset_title, "_IMP_by_DIV.csv"))
 imp_AUC_filename = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title, "_IMP_AUC.csv"))
 assay_component_map_filename = file.path(root_output_dir, "mea_nfa_component_name_map.csv")
 

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

#add imp_raw and imp_auc data
imp_raw <- fread(imp_rdata_filename)
imp_raw <- melt(imp_raw, id.vars = c(idcols, "DIV"), measure.vars = "wt.mean.resistance", variable.name = "src_acsn",
                value.name = "rval", variable.factor = FALSE)
imp_raw[, `:=`(src_acsn = paste0(src_acsn,"_DIV",DIV),
                srcf = basename(imp_rdata_filename))]
imp_raw[, DIV := NULL]
imp_raw[, srcf := basename(imp_rdata_filename)]

imp_auc <-  fread(imp_AUC_filename)
imp_auc <- melt(imp_auc, id.vars = idcols, measure.vars = "wt.mean.resistance_auc", variable.name = "src_acsn",
                value.name = "rval", variable.factor = FALSE)
imp_auc[, srcf := basename(imp_AUC_filename)]

# rbind DIV and AUC and IMP data
longdat <- rbind(DIV_data, AUC, imp_raw, imp_auc)
longdat[, `:=`(coli = as.numeric(sub("[[:alpha:]]","",well)), rowi = match(sub("[[:digit:]]","",well), LETTERS))]
longdat[, c("well") := list(NULL)]
setnames(longdat, old = "dose", new = "conc")
rm(list = c("DIV_data","AUC", "imp_auc", 'imp_raw'))

# add cytotox data
cytotox_data <- fread(cytotox_filename)
if(length(setdiff(names(longdat),names(cytotox_data))) > 0) cat('cytotox data does not have ',paste0(setdiff(names(longdat),names(cytotox_data)),collapse=","),'. Will fill with NA\n', sep ='')
if(length(setdiff(names(cytotox_data),names(longdat))) > 0) cat('MEA data does not have ',paste0(setdiff(names(cytotox_data),names(longdat)),collapse=","),'. Will fill with NA\n', sep = '')
longdat <- rbind(longdat, cytotox_data, fill = T)
rm(list = c("cytotox_data"))
longdat[, treatment := as.character(treatment)] # sometimes the treatment is read as an integer instead of a char

# replace the src_acsn with the TCPL acsn
#impedance has no tcpl acsn, will make up for imp
assay_component_map <- as.data.frame(read.table(assay_component_map_filename, sep = ",", header = T, skip = 0, stringsAsFactors = FALSE))
DIVs <- c("DIV0", "DIV2", "DIV5", "DIV7", "DIV9", "DIV12")
for (i in 1:6) {
assay_component_map[nrow(assay_component_map)+1, ] <- c(src_acsn = paste("wt.mean.resistance_", DIVs[i],sep = ""), acsn = paste("CCTE_Shafer_MEA_dev_Impedance_", DIVs[i], sep = ""))

}
assay_component_map[nrow(assay_component_map)+1, ] <- c(src_acsn = "wt.mean.resistance_auc", acsn = "CCTE_Shafer_MEA_dev_Impedance")

longdat <- merge(longdat, assay_component_map, by = c("src_acsn"), all.x = T)
if (any(is.na(unique(longdat$acsn)))) {
  print(longdat[is.na(acsn), unique(src_acsn)])
  stop(paste0("The above src_acsn's are not found in ",assay_component_map_filename))
}
longdat <- setDT(longdat)

# Define apid
longdat[, apid := paste(date, Plate.SN, sep = "_")]

# Assign wllt
longdat[, wllt := "t"]
longdat[ conc == 0, wllt := "n"]

# get the desired columns, in the desired order
longdat <- longdat[, .SD, .SDcols = intersect(c('apid', 'rowi', 'coli', 'treatment', 'conc', 'wllq_by_well', 'wllq_notes_by_well', 'wllt', 'rval', 'acsn', 'srcf', 'units'), names(longdat))]
# longdat may or may not include "units"
cat("long-format data is ready.\n")

dat <- longdat



# additional functions to prepare the data for TCPL mc0 format

#update_control_well_treatment <- function(dat, control_compound, culture_date = c(), plates = c(), control_rowi) {
 # apids <- Reduce(f = union, x = lapply(culture_date, function(x) grep(x, unique(dat$apid), val = T)))
  #if (length(plates) > 0) {
  #  apids <- Reduce(f = union, x = lapply(plates, function(x) grep(x, apids, val = T)))
  #}
  #cat("Control treatment will be updated to",control_compound,"for the following wells:\n")
  #print(dat[wllt == "n" & apid %in% apids & rowi, unique(.SD), .SDcols = c("apid","treatment","rowi","coli","conc")][order(apid,rowi,coli)])
  #dat[wllt == "n" & apid %in% apids & rowi, treatment := control_compound]
#}
#Cytotox control treatment is not specified
#need to manually update
#suggest merging metadata from maestroexperimentallog with plate.SN,well to the AB&LDH data to avoid doing this backward

# Updated treatment label for solvent control wells ------------------------------------
dat[, treatment_srcf := treatment]

# Note: Often the treatment name in the control wells is the same as the chemical treatment used in the same row
# But we know that these wells are controls because the concentration is 0
# Need to updated treatment name to the solvent control used.
dat[wllt == 'n', .N, by = .(treatment)] # wllt == 'n' determined where conc == 0

#treatment    N
#1:         Deltamethrin   12
#2:         Thiamethoxam   12
#3: Tributyltin Chloride   12
#4:             Rotenone   12
#5:           Glyphosate   12
#6:              Lindane   12
#7:            EtOH:DMSO  552
#8:                 DMSO 1656
#9:                  H2O  552
#10:                 EtOH  552

# Should all of these treatments be default_ControlTreatmentName?
# If so, use below:

# Can manually update other wells where control treatment is not the default, or use teh function below
# dat <- update_control_well_treatment(longdat, control_compound = "H20",culture_date = "20211103")
for (i in 1:6) {
dat[wllt == "n"& treatment == different_vehicleControlCompounds[i], ]$treatment <- different_vehicleControls[i]
}

#treatment    N
#1: EtOH:DMSO  564
#2:      DMSO 1692
#3:       H2O  564
#4:      EtOH  564

# Set the control well concentration. Adjust as needed
dat[wllt == "n", conc := 0.001]
dat[ , spid := treatment]
dat <- dat %>% replace_na(list(units = "uM")) 

# Assign sample ID's -------------------------------------------------------------
#spidmap <- as.data.table(read.xlsx(spidmap_file, sheet = spid_sheet))
#head(spidmap)
#unique(spidmap$Concentration_Unit) # all mM?
#unique(dat$units) # confirm these are all uM (this taken from maestroexperiment log file)
#setnames(spidmap, old = c(trt_col, spid_col), new = c("treatment","spid"))
# for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
#spidmap[, expected_stock_conc := 20] # initialize expected_stock_conc. Usually this is 20mM. Change as needed.
# update expected_stock_conc for individual compouunds where needed
# for example,
# spidmap[treatment %in% c("2,2',4,4',5,5'-Hexabromodiphenyl ether","Dibenz(a,h)anthracene"), expected_stock_conc := 10.0]
#spidmap[, treatment := as.character(treatment)]
#head(spidmap[, .(treatment, spid, expected_stock_conc)])

# Add additional spidmap's if needed and rbind into 1 spidmap

# check if every treatment name from the mea data maps to a unique sample in spidmap
#setdiff(dat$treatment, spidmap$treatment) # checkign for missed treatments
#spidmap[treatment %in% unique(dat$treatment), .N, by = .(treatment)][N > 1] # checking for treatments that match multiple spid
# if there is not a 1-to-1 correspondence, update treatment names in "supplemental_mea_treatment_name_map.csv"

# update treatment names with entries in "supplemental_mea_treatment_name_map.csv" corresponding to dataset
# (treatment -> "mea_treatment_name", "updated_treatment_name" column will match "PREFERRED_NAME"
#dat <- update_treatment_names(dat, root_output_dir, dataset_title)

# assign spids
#dat <- check_and_assign_spids(dat, spidmap)


# Confirm Conc's ----------------------------------------------------------------
# confirm that the conc's collected from master chem lists and Calc files match
# and that the correct concentration-corrections has been done for each compound
dat[, conc_srcf := conc] # save the original conc's in a column

# check if there are multiple conc's assigned to the same well (usually occurs if there are differences between master chem file and calc file)
# Note: in TCPL mc1, the conc's are set to dat[ , conc := signif(conc, 3)]. So it's okay for us to round here.
dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1]
# if any, standardize those before continuing.
problem_comps <- dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1, unique(treatment)]
problem_comps
# character(0)
save(dat, file = paste(root_output_dir, dataset_title, "output", "Impedance2023_longfile.RData", sep = "/"))

# finally, run this:
source(file.path(scripts.dir, 'confirm_concs.R'))
con <- dbConnect(drv = RMySQL::MySQL(), user = "", pass = "", dbname='',host = "")
dat <- confirm_concs(dat, spidmap, con, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30), update_concs_without_prompt = update_concs_without_prompt)
dbDisconnect(con)


# FINAL DATA CHECKS -------------------------------------------------------------
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Check for the expected number of technical replicates
dat[wllt == 't', .(length(unique(paste0(apid,rowi,coli)))), by = .(spid, conc)][V1 != 3]
# do you except these cases to have more than or less than 3 replicates?
# Were some samples repeated, and only certain repeats meant to be included?
#### Rotenone and Tributyltin Chloride are repeated on a lower concentration range, but it did not show up in the check


# Any other plots or things to check?

# save dat and graphs
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
