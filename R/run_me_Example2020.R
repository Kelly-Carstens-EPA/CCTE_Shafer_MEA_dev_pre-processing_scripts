rm(list=ls()) # clear environment
graphics.off() # clear plot history
####################################################################################
# USER INPUT
###################################################################################
dataset_title <- "Example2020" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- FALSE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName = "DMSO" # usually DMSO. all compounds other than those listed below should have this vehicle control
# Enter the names of the compounds as they appear in the MEA data that have a vehicle control other than the default
different_vehicleControlCompounds = c() # e.g. c("Sodium Orthovanadate", "Amphetamine")
# Enter the names of the vehicle controls as they correspond to the compounds in the previous list
different_vehicleControls = c() # e.g. c("Water", "Water")

spidmap_file <- "L:/Lab/NHEERL_MEA/Project PFAS 2018/EPA_9238_EPA-Shafer_75_20180511_key_MW Waste Calculations.xlsx"
spid_sheet <- "Worksheet1 (2)"

scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R"
root_output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the dataset_title folder will be created
###################################################################################
# END USER INPUT
###################################################################################

library(data.table)
library(openxlsx)

# create a summary log file and store the 
if(save_notes_graphs) {
  sink(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_run_log_",as.character.Date(Sys.Date()),".txt")))
  cat("Output from the script run_me_",dataset_title,".R\n",sep="")
  cat("Date:",as.character.Date(Sys.Date()),"\n")
  cat("USER INPUT settings:\n")
  print(sapply(ls(), get, envir = .GlobalEnv))
  graphics.off()
  pdf(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")))
}

# run the main steps
source(file.path(scripts.dir, 'source_steps.R'))

# run tcpl_MEA_dev_AUC
source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
dat <- tcpl_MEA_dev_AUC(basepath = file.path(root_output_dir,dataset_title), dataset_title)

# change untreated wells to Control Treatment ------------------------------------
dat[wllt == "n", treatment := default_ControlTreatmentName]
# dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "")


# Assign SPIDs ------------------------------------------------------------------
spidmap <- as.data.table(read.xlsx(spidmap_file, sheet = spid_sheet))
head(spidmap)
unique(spidmap$Concentration_Unit) # all mM?
setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
spidmap[, expected_stock_conc := 30] # initialize expected_stock_conc. Usually this is 20mM. Change as needed.
# update expected_stock_conc for individual compouunds where needed 
# for example, 
# spidmap[treatment %in% c("2,2',4,4',5,5'-Hexabromodiphenyl ether","Dibenz(a,h)anthracene"), expected_stock_conc := 10.0]
spidmap[, treatment := as.character(treatment)]
spidmap[, stock_conc := as.numeric(stock_conc)]
head(spidmap[, .(treatment, spid, stock_conc, expected_stock_conc)])

# # rename any compounds, if needed
# auc <- fread(file.path(root_output_dir,dataset_title, "output", paste0(dataset_title, "_AUC.csv")))
# cyto <- fread(file.path(root_output_dir,dataset_title, "output", paste0(dataset_title, "_cytotox.csv")))
# auc[treatment == "Dibenz[a,c] anthracene", treatment := "Dibenz[a,c]anthracene"]
# cyto[treatment == "Dibenz[a,c] anthracene", treatment := "Dibenz[a,c]anthracene"]
# write.csv(auc, file.path(root_output_dir,dataset_title, "output", paste0(dataset_title, "_AUC.csv")), row.names = FALSE)
# write.csv(cyto, file.path(root_output_dir,dataset_title, "output", paste0(dataset_title, "_cytotox.csv")), row.names = FALSE)
# rm(list = c("auc","cyto"))

# assign spids
dat <- check_and_assign_spids(dat, spidmap)


# Confirm Conc's ----------------------------------------------------------------
# confirm that the conc's collected from master chem lists and Calc files match
# and that the correct concentration-corrections has been done for each compound
source(file.path(scripts.dir, 'confirm_concs.R'))
dat <- confirm_concs(dat, spidmap, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30))


# FINAL DATA CHECKS
# this section is to confirm that the data has been processed correctly
dat <- read.csv(file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.csv")))
setDT(dat)
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Any other plots or things to check?

rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
