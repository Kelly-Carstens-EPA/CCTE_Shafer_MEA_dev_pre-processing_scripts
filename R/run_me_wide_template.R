rm(list=ls()) # clear environment
graphics.off() # clear plot history
###################################################################################
# USER INPUT
###################################################################################
project_name <- "testpipeline2020" # the name for the current dataset, e.g. "name2020" (this should match the name of the output folder)
pause_between_steps <- TRUE # leave this as TRUE for the first run
save_notes_graphs <- FALSE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

# default_ControlTreatmentName <- "DMSO" # all compounds other than those listed below should have this vehicle control

project.dir <- "" # project main folder (where will look for README files)
scripts.dir <- "" # update to the folder where the scripts are located
root.output.dir <- "" # where the project_name folder will be created
###################################################################################
# END USER INPUT
###################################################################################

library(data.table)

# create a summary log file and store the 
if(save_notes_graphs) {
  sink(file = file.path(root.output.dir, project_name, paste0(project_name,"_run_log_",as.character.Date(Sys.Date()),".txt")))
  cat("Output from the script run_me_",project_name,".R\n",sep="")
  cat("Date Ran:",as.character.Date(Sys.Date()),"\n")
  cat(R.version.string,"\n")
  cat("USER INPUT settings:\n")
  print(sapply(ls(), get, envir = .GlobalEnv))
  graphics.off()
  pdf(file = file.path(root.output.dir, project_name, paste0(project_name,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")), width = 10, height = 8)
}

# Scan for readme's that might affect dosing, wllq
txt.files <- list.files(project.dir, pattern = '\\.txt', recursive = T, full.names = T)
readmes <- txt.files[grepl('read( )*me',tolower(txt.files))]
for (readme in readmes) {
  cat(dirname(readme),'\n')
  cat(scan(readme, what = character(), sep = '\n', quiet = T), sep = '\n')
  cat('\n')
}

# run the main steps
source(file.path(scripts.dir, 'source_steps.R'))

# Normalize the AUC values to controls
# Values will be normalized to the median of the controls on every unique Plate.SN and date
# Control wells are assumed to be where the dose is 0
source(file.path(scripts.dir, 'normalize_auc_summary.R'))
auc_table <- read.csv(file.path(root.output.dir, project_name, "output", paste0(project_name,"_AUC.csv")), stringsAsFactors = F)

# normalize in down direction (controls at 100%, down response goes to 0)
auc_table_normalized_dn <- auc_summary(auc_table, direction = 'down')
write.csv(auc_table_normalized_dn, file = file.path(root.output.dir, project_name, "output", paste0(project_name,"_AUC_normalized_down.csv")), row.names = FALSE)

# normalize in up direction (controls at 100%, up response goes to 0)
auc_table_normalized_up <- auc_summary(auc_table, direction = 'up')
write.csv(auc_table_normalized_up, file = file.path(root.output.dir, project_name, "output", paste0(project_name,"_AUC_normalized_up.csv")), row.names = FALSE)

# Skipping:
# - melt data and merge cytotox
# - update control well treatment and conc
# - assign sample IDs
# - confirming conc's consistent between master chem lists and calc files
# - make conc's match stock conc provided in spid file


# FINAL DATA CHECKS -------------------------------------------------------------
# this section is to confirm that the data has been processed correctly
rm(list = setdiff(ls(),c(grep('auc',ls(),val = T),'save_notes_graphs','scripts.dir','root.output.dir','project_name')))
source(file.path(scripts.dir, 'dataset_checks.R'))

# Run dataset checks for each table
dataset_checks_wide(auc_table)
dataset_checks_wide(auc_table_normalized_dn, normalized = T, direction = 'Down')
dataset_checks_wide(auc_table_normalized_up, normalized = T, direction = 'Up')

# confirm that the ranges for the normalized values appear correct

# Do any other desired checks here:

# Save output
if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
