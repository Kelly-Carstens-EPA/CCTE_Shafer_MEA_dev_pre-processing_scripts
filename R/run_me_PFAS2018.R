rm(list=ls())
graphics.off()
###################################################################################
# USER INPUT
###################################################################################
dataset_title <- "PFAS2018" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- TRUE # can say no if you are in debugging phase, but do do a final run with saving notes adn graphs

default_ControlTreatmentName = "DMSO" # all compounds other than those listed below should have this vehicle control

spidmap_file <- "L:/Lab/NHEERL_MEA/Project PFAS 2018/EPA_9238_EPA-Shafer_75_20180511_key_MW Waste Calculations.xlsx"
spid_sheet <- "Worksheet1 (2)"

scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R"
root_output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the dataset_title folder will be created

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
  cat("Date:",as.character.Date(Sys.Date()),"\n")
  cat("USER INPUT settings:\n")
  print(sapply(ls(), get, envir = .GlobalEnv))
  graphics.off()
  pdf(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")))
}

# source the ultimate function!
source(file.path(scripts.dir, 'source_steps.R'))

# run tcpl_MEA_dev_AUC
source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
dat <- tcpl_MEA_dev_AUC(basepath = file.path(root_output_dir,dataset_title), dataset_title)


# change untreated wells to Control Treatment ------------------------------------
dat[wllt == "n", treatment := default_ControlTreatmentName]
# dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "")
dat[wllt == "n", conc := 0.001]

# prepare spidmap
spidmap <- as.data.table(read.xlsx(spidmap_file, sheet = spid_sheet))
head(spidmap)
unique(spidmap$Unit) # all mM
setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
spidmap[, expected_stock_conc := 30] # initialize expected_stock_conc. Usually this is 20mM. Change as needed.
# update expected_stock_conc for individual compouunds where needed 
# for example, 
# spidmap[treatment %in% c("2,2',4,4',5,5'-Hexabromodiphenyl ether","Dibenz(a,h)anthracene"), expected_stock_conc := 10.0]
spidmap[, treatment := as.character(treatment)]
spidmap[, stock_conc := as.numeric(stock_conc)]
head(spidmap[, .(treatment, spid, stock_conc, expected_stock_conc)])
spidmap <- spidmap[, .(treatment, spid, stock_conc, expected_stock_conc)]

spidmap2 <- as.data.table(read.xlsx(file.path(root_output_dir,"Sample IDs","Shafer_sample_info_to_register_20201110_afc.xlsx"), sheet = 1))
spidmap2 <- spidmap2[dataset == dataset_title]
spidmap2[, `:=`(expected_stock_conc = stock_conc)]
setnames(spidmap2, old = c("SPID","compound"), new = c("spid","treatment"))
usecols <- c("spid","treatment","stock_conc","expected_stock_conc")
spidmap <- rbind(spidmap[, ..usecols], spidmap2[, ..usecols])

# # just confirming that we have all of the spids
# setdiff(parameter_data$trt, spidmap$treatment) # "Acetaminophen" "Bisphenol A"   "Loperamide"
# dat[treatment == "Loperamide", treatment := "4-(4-Chlorophenyl)-4-hydroxy-N,N-dimethyl-alpha,alpha-diphenylpiperidine-1-butyramide monohydrochloride"]

# add "updated_treatment_name" columne to match "PREFERRED_NAME" in spidmap, reaname treatment to "mea_treatment_name"
dat <- update_treatment_names(dat, root_output_dir, dataset_title)

# assign spids
dat <- check_and_assign_spids(dat, spidmap)


# Confirm Conc's ----------------------------------------------------------------
# confirm that the conc's collected from master chem lists and Calc files match
# and that the correct concentration-corrections has been done for each compound

# check if there are multiple conc's assigned to the same well (usually occurs if there are differences between master chem file and calc file)
# Note: in TCPL mc1, the conc's are set to dat[ , conc := signif(conc, 3)]. So it's okay for us to round here.
dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1]
# if any, standardize those before continuing.
problem_comps <- dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1, unique(treatment)]
problem_comps
# character(0)

# finally, run this:
source(file.path(scripts.dir, 'confirm_concs.R'))
dat <- confirm_concs(dat, spidmap, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30), update_concs_without_prompt = update_concs_without_prompt)


# FINAL DATA CHECKS
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Any other plots or things to check?
dat[acsn == "CCTE_Shafer_MEA_dev_inter_network_spike_interval_mean_DIV5", unique(rval)] # NA
dat[acsn == "CCTE_Shafer_MEA_dev_network_spike_duration_std_DIV5", unique(rval)] # NA
dat[acsn == "CCTE_Shafer_MEA_dev_network_spike_duration_std_DIV5", min(rval, na.rm = T)] 
# [1] Inf 
# Warning message:
#   In min(rval, na.rm = T) : no non-missing arguments to min; returning Inf
# yep, this is warnign the warning is coming from
# okay, this is reasonable network spike endpoints on for DIV 5

# checking data that got dropped previously
dat[grepl("MW1208-6",apid), .N, by = "acsn"] # 48 pts for all 87 acsn, including each DIV
dat[grepl("MW1208-7",apid) & grepl("mutual_information",acsn), summary(rval)]
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000013 0.001425 0.002784 0.008230
dat[grepl("MW1208-8",apid) & grepl("mutual_information_norm_DIV5",acsn), summary(rval)] # all 0's, but present.
dat[grepl("MW1208-8",apid) & grepl("DIV5",acsn), summary(rval)] # p resent
dat[grepl("MW1208-8",apid) & grepl("DIV5",acsn), length(unique(acsn))] # 17
dat[grepl("MW1230-53",apid) & grepl("DIV9",acsn), .N, by = "acsn"] # 48 pts for each acsn!

source(file.path(root_output_dir,"supplemental_scripts","view_replicates_by_DIV.R"))
div_dat <- as.data.table(read.csv(file.path(root_output_dir,dataset_title,"output",paste0(dataset_title,"_parameters_by_DIV.csv"))))
view_replicates_by_DIV(div_dat, "Acetaminophen", dosei = 0.3, title_msg = "verifying 1230-53 estimated DIV9")
view_replicates_by_DIV(div_dat, "Bisphenol A", dosei = 10, title_msg = "verifying 1230-53 estimated DIV9")
# nAE looks great! MFR could be better, but it's okay

# checking out these treatments with cytotoxicity "outliers"
plotdat <- dat[grepl("LDH",acsn)]
stripchart(rval ~ signif(log10(conc),3), plotdat[treatment == "4-(4-Chlorophenyl)-4-hydroxy-N,N-dimethyl-alpha,alpha-diphenylpiperidine-1-butyramide monohydrochloride"],
           vertical = T, pch = 1, main = paste0(unique(plotdat$acsn)," PFAS2018 Loperamide Dose-response"))
abline(h = plotdat[apid == "20181114_MW1234-49" & wllt == "n" & wllq == 1, median(rval)])

stripchart(rval ~ signif(log10(conc),3), plotdat[treatment == "1475813"],
           vertical = T, pch = 1, main = paste0(unique(plotdat$acsn)," PFAS2018 1475813 Dose-response"))
abline(h = plotdat[apid == "20181114_MW1234-25" & wllt == "n" & wllq == 1, median(rval)])

stripchart(rval ~ signif(log10(conc),3), plotdat[treatment == "1475815"],
           vertical = T, pch = 1, main = paste0(unique(plotdat$acsn)," PFAS2018 1475815 Dose-response"))
abline(h = plotdat[apid == "20181017_MW1207-43" & wllt == "n" & wllq == 1, median(rval)])

stripchart(rval ~ signif(log10(conc),3), plotdat[treatment == "1475824"],
           vertical = T, pch = 1, main = paste0(unique(plotdat$acsn)," PFAS2018 1475824 Dose-response"))
abline(h = plotdat[apid == "20181017_MW1208-4" & wllt == "n" & wllq == 1, median(rval)])

stripchart(rval ~ signif(log10(conc),3), dat[grepl("AB",acsn) & treatment == "1475815"],
           vertical = T, pch = 1, main = paste0("AB"," PFAS2018 1475815 Dose-response"))
abline(h = dat[grepl("AB",acsn) & apid == "20181017_MW1207-43" & wllt == "n" & wllq == 1, median(rval)])

# confirming that the conc's have been udpated for 20181017
dat[grepl("20181017",apid) & conc < 0.1, .N, by = .(wllt, conc, coli, apid)]
# wllt  conc coli               apid   N
# 1:    t 0.030    1 20181017_MW1207-43 522
# 2:    t 0.030    1 20181017_MW1207-44 522
# 3:    t 0.030    1  20181017_MW1208-1 522
# 4:    t 0.030    1  20181017_MW1208-2 522
# 5:    t 0.030    2  20181017_MW1208-3 522
# 6:    t 0.030    1  20181017_MW1208-4 522
# 7:    n 0.001    2 20181017_MW1207-43 522
# 8:    n 0.001    2 20181017_MW1207-44 522
# 9:    n 0.001    2  20181017_MW1208-1 522
# 10:    n 0.001    2  20181017_MW1208-2 522
# 11:    n 0.001    1  20181017_MW1208-3 522
# 12:    n 0.001    2  20181017_MW1208-4 522

# save dat and graphs
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")

# # manual fix for 20181017
# update_files <- c(list.files(file.path(root_output_dir,dataset_title,"prepared_data"),pattern = "ont_data_summary_ABEfilt_20181017_", full.names = T),
#                   list.files(file.path(root_output_dir,dataset_title,"All_MI"),pattern = "NMI_20181017_", full.names = T))
# update_files <- update_files[!grepl("MW1208-3",update_files)]
# 
# for (filei in update_files) {
#   cat(basename(filei),"\n")
#   pdat <- as.data.table(read.csv(filei, stringsAsFactors = F))
#   pdat[grepl("2",well), dose := 0] # controls are in col 2
#   pdat[grepl("1",well), dose := 0.03] # 0.03 is in col 1
#   write.csv(pdat, filei,row.names = F)
# }
