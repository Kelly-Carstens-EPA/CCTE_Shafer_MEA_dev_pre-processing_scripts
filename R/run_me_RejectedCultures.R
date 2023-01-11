rm(list=ls())
graphics.off()
###################################################################################
# USER INPUT
###################################################################################
dataset_title <- "RejectedCultures" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- TRUE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName = "Control" # usually DMSO. I don't want to search all of this, so I will just leave it generic

scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R"
root_output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the dataset_title folder will be created

spidmap_file <- ""
spid_sheet <- ""

update_concs_without_prompt <- TRUE
###################################################################################
# END USER INPUT
###################################################################################

library(data.table)
library(openxlsx)

# create a summary log file and store the 
if(save_notes_graphs) {
  sink(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_run_log_",as.character.Date(Sys.Date()),".txt")))
  cat("Output from the script run_me_",dataset_title,".R\n",sep="")
  cat("Date Ran:",as.character.Date(Sys.Date()),"\n")
  cat(R.version.string,"\n")
  cat("USER INPUT settings:\n")
  print(sapply(ls(), get, envir = .GlobalEnv))
  graphics.off()
  pdf(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")), height = 8, width = 10)
}

# run the main steps
source(file.path(scripts.dir, 'source_steps.R'))

# run tcpl_MEA_dev_AUC
source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
dat <- tcpl_MEA_dev_AUC(basepath = file.path(root_output_dir,dataset_title), dataset_title)

# change untreated wells to Control Treatment ------------------------------------
dat[wllt == "n", treatment := default_ControlTreatmentName]
# update other control wells as needed, e.g.
# dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "20190904", plates = paste0("MW69-381",7:9), control_rowi = which(LETTERS[1:6] %in% c("E","F")))
# leaving teh conc as 0, since this won't go into tcpl regardless


# Assign SPIDs ------------------------------------------------------------------
# skipping this, since this data won't be added to tcpl
# spidmap <- as.data.table(read.xlsx(spidmap_file1, sheet = spid_sheet1))
# head(spidmap)
# unique(spidmap$ALIQUOT_CONC_UNIT) # all mM
# unique(spidmap$ALIQUOT_SOLVENT) # DMSO
# setnames(spidmap, old = c("preferred_name","ALIQUOT_CONC", "EPA_SAMPLE_ID"), new = c("treatment","stock_conc","spid"))
# # for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
# spidmap[, treatment := as.character(treatment)]
# head(spidmap[, .(treatment, spid, stock_conc)])
# 
# # update names in "treatment" col to match "PREFERRED_NAME" in spidmap, set original treatments col to "mea_treatment_name"
# # dat <- update_treatment_names(dat, root_output_dir, dataset_title)
# # don't need to worry about treatment names for this data set
# 
# # assign spids
# dat <- check_and_assign_spids(dat, spidmap)


# Confirm Conc's ----------------------------------------------------------------
# confirm that the conc's collected from master chem lists and Calc files match
# and that the correct concentration-corrections has been done for each compound

# check if there are multiple conc's assigned to the same well (usually occurs if there are differences between master chem file and calc file)
# Note: in TCPL mc1, the conc's are set to dat[ , conc := signif(conc, 3)]. So it's okay for us to round here.
dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1]
# Empty data.table (0 rows and 5 cols): treatment,apid,rowi,coli,num_unique_concs_in_well
# if any, standardize those before continuing.

# finally, run this - again, skipping this for now
# source(file.path(scripts.dir, 'confirm_concs.R'))
# dat <- confirm_concs(dat, spidmap, expected_target_concs = c(0.03,0.1,0.3,1,3,10,20), update_concs_without_prompt = update_concs_without_prompt)


# FINAL DATA CHECKS
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# save the data and graphs
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # save the plots to the pdf device
}

cat("\nDone!\n")
