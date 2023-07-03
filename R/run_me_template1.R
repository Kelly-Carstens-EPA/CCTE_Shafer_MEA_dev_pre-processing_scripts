rm(list=ls()) # clear environment
graphics.off() # clear plot history
# ------------------------------------------------------------------------ #
# USER INPUT
# ------------------------------------------------------------------------ #
project_name <- "testpipeline2020" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- FALSE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName <- "DMSO" # all compounds other than those listed below should have this vehicle control

spidmap_file <- ""
spid_sheet <- ""

project.dir <- "" # project main folder (where will look for README files)
scripts.dir <- "" # update to the folder where the scripts are located
root_output_dir <- "" # where the project_name folder will be created

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
  sink(file = file.path(root_output_dir, project_name, paste0(project_name,"_run_log_",as.character.Date(Sys.Date()),".txt")))
  cat("Output from the script run_me_",project_name,".R\n",sep="")
  cat("Date Ran:",as.character.Date(Sys.Date()),"\n")
  cat(R.version.string,"\n")
  cat("USER INPUT settings:\n")
  print(sapply(ls(), get, envir = .GlobalEnv))
  graphics.off()
  pdf(file = file.path(root_output_dir, project_name, paste0(project_name,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")), width = 10, height = 8)
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

get_NFA_standard_structue <- function(culture.folderi) {
  cyto.files <- list.files(path = culture.folderi, pattern = '(Calculations)|(Summary)', full.names = T)
  plate.dirs <- list.files(path = culture.folderi, pattern = '^(MW)*[0-9\\-]{7}$', full.names = T)
  mfiles <- list.files(path = culture.folderi, pattern = 'MaestroExperimentLog', full.names = T, recursive = T)
  slists <- list.files(path = culture.folderi, pattern = '_spike_list\\.csv', full.names = T, recursive = T)
  filesi <- c(cyto.files, mfiles, slists)
  return(filesi)
}

culture.folders <- list.files(path = project.dir, pattern = "[0-9]{8}", full.names = T, recursive = F)
all.files <- c()
for (culture.folderi in culture.folders) {
  all.files <- c(all.files, get_NFA_standard_structue(culture.folderi))
}

#Basic cleaning
all.files <- Filter(function(filei) !grepl('\\Q~$\\E',filei), all.files)
all.files <- Filter(function(filei) !grepl('deprecated',tolower(filei)), all.files)
length(all.files)
# should be # groups * (1 Calc file + 3 plates * (4 DIV + 1 maestro exp log file))

# Write to log file
writeLogFile(all.files, output.dir = file.path(root_output_dir, project_name), project_name, files_type = '')
# Writing 288 files to DNT_NTP2021_files_log_2022-04-19.txt ...
# [1] "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/DNT_NTP2021/DNT_NTP2021_files_log_2022-04-19.txt is ready."
rm(all.files, culture.folders, readmes)


# > run the main steps ----------------------------------------------------

# Can enter "a" to select files, then hit cancel, to see summary of files selected
source(file.path(scripts.dir, 'source_steps.R'))



# > run tcpl_MEA_dev_AUC --------------------------------------------------

source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
dat <- tcpl_MEA_dev_AUC(basepath = file.path(root_output_dir,project_name), project_name)


# Updated treatment label for solvent control wells ------------------------------------
dat[, treatment_srcf := treatment]

# Note: Often the treatment name in the control wells is the same as the chemical treatment used in the same row
# But we know that these wells are controls because the concentration is 0
# Need to updated treatment name to the solvent control used.
dat[wllt == 'n', .N, by = .(treatment)] # wllt == 'n' determined where conc == 0
# Should all of these treatments be default_ControlTreatmentName?
# If so, use below:
# dat[wllt == "n", treatment := default_ControlTreatmentName]
# Can manually update other wells where control treatment is not the default, or use teh function below
# dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "")

# Set the control well concentration. Adjust as needed
dat[wllt == "n", conc := 0.001]


# Assign sample ID's -------------------------------------------------------------
spidmap <- as.data.table(read.xlsx(spidmap_file, sheet = spid_sheet))
head(spidmap)
unique(spidmap$Concentration_Unit) # all mM?
unique(dat$units) # confirm these are all uM (this taken from maestroexperiment log file)
setnames(spidmap, old = c(trt_col, spid_col), new = c("treatment","spid"))
# for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
spidmap[, expected_stock_conc := 20] # initialize expected_stock_conc. Usually this is 20mM. Change as needed.
# update expected_stock_conc for individual compouunds where needed
# for example,
# spidmap[treatment %in% c("2,2',4,4',5,5'-Hexabromodiphenyl ether","Dibenz(a,h)anthracene"), expected_stock_conc := 10.0]
spidmap[, treatment := as.character(treatment)]
head(spidmap[, .(treatment, spid, expected_stock_conc)])

# Add additional spidmap's if needed and rbind into 1 spidmap

# check if every treatment name from the mea data maps to a unique sample in spidmap
setdiff(dat$treatment, spidmap$treatment) # checkign for missed treatments
spidmap[treatment %in% unique(dat$treatment), .N, by = .(treatment)][N > 1] # checking for treatments that match multiple spid
# if there is not a 1-to-1 correspondence, update treatment names in "supplemental_mea_treatment_name_map.csv"

# update treatment names with entries in "supplemental_mea_treatment_name_map.csv" corresponding to dataset
# (treatment -> "mea_treatment_name", "updated_treatment_name" column will match "PREFERRED_NAME"
dat <- update_treatment_names(dat, root_output_dir, project_name)

# assign spids
dat <- check_and_assign_spids(dat, spidmap)


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

# Any other plots or things to check?

# save dat and graphs
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, project_name, "output", paste0(project_name,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
