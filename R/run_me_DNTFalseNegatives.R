rm(list=ls()) # clear environment
graphics.off() # clear plot history
###################################################################################
# USER INPUT
###################################################################################
dataset_title <- "DNTFalseNegatives" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- TRUE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName <- "DMSO" # all compounds other than those listed below should have this vehicle control

spidmap_file <- ""
spid_sheet <- ""

scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R" # update to the folder where the scripts are located
root_output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the dataset_title folder will be created

update_concs_without_prompt <- FALSE
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
  pdf(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")), width = 10, height = 8)
}


# run the main steps
source(file.path(scripts.dir, 'source_steps.R'))

# run tcpl_MEA_dev_AUC
source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
dat <- tcpl_MEA_dev_AUC(basepath = file.path(root_output_dir,dataset_title), dataset_title)


# change untreated wells to Control Treatment ------------------------------------
dat[, .N, by = .(treatment)]
dat[wllt == 'n' & treatment != 'DMSO', .N, by = .(srcf, conc)]
# srcf conc_srcf  N
# 1: 20210818_NFA_False_Negatice_Repeats__Calculations.xlsx         0 36
# Only not labelled as DMSO in the cytotox data. But I know should be solvent control bc treatment conc is 0
dat[wllt == "n", treatment := default_ControlTreatmentName]
# Manually update other wells where control treatment is not the default, or use teh function below
# dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "")

# Set the control well concentration. Adjust as needed
dat[wllt == "n", conc := 0.001]


# Assign sample ID's -------------------------------------------------------------
dat[, .N, by = .(treatment)]
#                treatment    N
# 1:                 6-PPD 1827
# 2:                  DMSO 1566
# 3:         6-PPD Quinone 1827
# 4:              Caffeine 1827
# 5: 5,5-Diphenylhydantoin 1827
# 6:         Dexamethasone 1827
# 7:                 Maneb 1827
# I don't have any spids for these compounds
# but the compound names are consistent! So we'll just stick with these for now
dat[, spid := treatment]

# Will run this later when have spids
# spidmap <- as.data.table(read.xlsx(spidmap_file, sheet = spid_sheet))
# head(spidmap)
# unique(spidmap$Concentration_Unit) # all mM?
# unique(dat$units) # confirm these are all uM (this taken from maestroexperiment log file)
# setnames(spidmap, old = c(trt_col, spid_col), new = c("treatment","spid"))
# # for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
# spidmap[, expected_stock_conc := 20] # initialize expected_stock_conc. Usually this is 20mM. Change as needed.
# # update expected_stock_conc for individual compouunds where needed 
# # for example, 
# # spidmap[treatment %in% c("2,2',4,4',5,5'-Hexabromodiphenyl ether","Dibenz(a,h)anthracene"), expected_stock_conc := 10.0]
# spidmap[, treatment := as.character(treatment)]
# head(spidmap[, .(treatment, spid, expected_stock_conc)])
# 
# # Add additional spidmap's if needed and rbind into 1 spidmap
# 
# # check if every treatment name from the mea data maps to a unique sample in spidmap
# setdiff(dat$treatment, spidmap$treatment) # checkign for missed treatments
# spidmap[treatment %in% unique(dat$treatment), .N, by = .(treatment)][N > 1] # checking for treatments that match multiple spid
# # if there is not a 1-to-1 correspondence, update treatment names in "supplemental_mea_treatment_name_map.csv"
# 
# # update treatment names with entries in "supplemental_mea_treatment_name_map.csv" corresponding to dataset
# # (treatment -> "mea_treatment_name", "updated_treatment_name" column will match "PREFERRED_NAME"
# dat <- update_treatment_names(dat, root_output_dir, dataset_title)
# 
# # assign spids
# dat <- check_and_assign_spids(dat, spidmap)


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

# Since these spids haven't been registered, just check out the conc's
for (treatmenti in unique(dat$treatment)) {
  cat(treatmenti, '\n')
  print(dat[treatment == treatmenti, table(conc)])
}
# 6-PPD 
# conc
# 0.03  0.1  0.3    1    3   10   30 
# 261  261  261  261  261  261  261 
# DMSO 
# conc
# 0.001 
# 1566 
# 6-PPD Quinone 
# conc
# 0.00267  0.0089  0.0267   0.089   0.267    0.89    2.67 
# 261     261     261     261     261     261     261 
# Caffeine 
# conc
# 0.1 0.3   1   3  10  30 100 
# 261 261 261 261 261 261 261 
# 5,5-Diphenylhydantoin 
# conc
# 1    3   10   30  100  300 1000 
# 261  261  261  261  261  261  261 
# Dexamethasone 
# conc
# 0.1 0.3   1   3  10  30 100 
# 261 261 261 261 261 261 261 
# Maneb 
# conc
# 0.01 0.03  0.1  0.3    1    3   10 
# 261  261  261  261  261  261  261 
# Looks good & consistent!

# # finally, run this:
# source(file.path(scripts.dir, 'confirm_concs.R'))
# con <- dbConnect(drv = RMySQL::MySQL(), user = "", pass = "", dbname='',host = "")
# dat <- confirm_concs(dat, spidmap, con, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30), update_concs_without_prompt = update_concs_without_prompt)
# dbDisconnect(con)


# > Concentration units ---------------------------------------------------

dat[, .N, by = .(treatment, units)]
# treatment units    N
# 1:                 6-PPD  <NA>   42
# 2:                  DMSO  <NA>   36
# 3:         6-PPD Quinone  <NA>   42
# 4:              Caffeine  <NA>   42
# 5: 5,5-Diphenylhydantoin  <NA>   42
# 6:         Dexamethasone  <NA>   42
# 7:                 Maneb  <NA>   42
# 8:                 6-PPD    uM 1785
# 9:                  DMSO    uM 1275
# 10:         6-PPD Quinone µg/ml 1785
# 11:                  DMSO µg/ml  255
# 12:              Caffeine    uM 1785
# 13: 5,5-Diphenylhydantoin    uM 1785
# 14:         Dexamethasone    uM 1785
# 15:                 Maneb    uM 1785

# Update DMSO
dat[treatment == 'DMSO', units := 'fraction']

# Fill where is NA
dat[is.na(units), .N, by = .(srcf)]
# srcf   N
# 1: 20210818_NFA_False_Negatice_Repeats__Calculations.xlsx 288
# units not defined in calculations files
# But I know same concentrations used in cyto assay as in NFA (bc same source wells!)
# So I will just update with the units from Meastro log file for MEA endpoints
dat[, .(length(unique(units[!is.na(units)]))), by = .(treatment, spid)][V1 > 1] # empty, good
dat[, units := unique(units[!is.na(units)]), by = .(treatment, spid)]

# Convert to uM
dat[treatment == '6-PPD Quinone', mol_weight_g_per_mol := 298.38] # see notebook "RE: MW for 6-ppd quinone"
dat[units == 'uM', conc_in_uM := conc]
dat[units == 'µg/ml', conc_in_uM := conc*1000/mol_weight_g_per_mol] # updated from ug/mL on Apr 25, since that is different than what is recorded under "units"


# FINAL DATA CHECKS -------------------------------------------------------------
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Check for the expected number of technical replicates
dat[wllt == 't', .(length(unique(paste0(apid,rowi,coli)))), by = .(spid, conc)][V1 != 3]
# do you except these cases to have more than or less than 3 replicates?
# Were some samples repeated, and only certain repeats meant to be included?

# Any other plots or things to check?


# Other checks ------------------------------------------------------------

# > Confirming DIV 12 estimated only for affected wells on plate ----------

# Of the AUC endpoints, which wells from this plate relied on estimated values?
dat[!grepl('DIV',acsn) & !grepl('(LDH)|(AB)',acsn) & grepl('estimated',wllq_notes), .N, by = .(apid, rowi, coli)][order(rowi, coli)]
#                  apid rowi coli  N
# 1: 20210818_MW75-9207    5    1 17
# 2: 20210818_MW75-9207    6    1 17
# 3: 20210818_MW75-9207    6    7 17
# 4: 20210818_MW75-9207    6    8 17
# cool, just these 4 wells

# And the DIV12?
dat[grepl('DIV12',acsn) & grepl('estimated',wllq_notes), .N, by = .(apid, rowi, coli)][order(rowi, coli)]
# apid rowi coli  N
# 1: 20210818_MW75-9207    5    1 17
# 2: 20210818_MW75-9207    6    1 17
# 3: 20210818_MW75-9207    6    7 17
# 4: 20210818_MW75-9207    6    8 17



# > Add'l wllq notes ------------------------------------------------------

# I'm guessing 5% CO2 note got dropped, bc only afffected DIV12?
dat[apid %in% c('20210818_MW75-9206','20210818_MW75-9207'), .N, by = .(wllq_notes)]
# ya, not present
# Just add this note for MEA AUC adn DIV12 endpoints
# bu tdon't change wllq
dat[apid %in% c('20210818_MW75-9206','20210818_MW75-9207') & !(grepl('(LDH)|(AB)',acsn) | grepl('DIV[579]',acsn)), `:=`(wllq_notes = paste0(wllq_notes, '; 5% CO2 ran out during reading'))]

# Add note for precipitate
dat[grepl('precipitate',wllq_notes)]

# Load HCI data with wllq notes info
load(file.path(root_output_dir, '../pre-process_hci_for_tcpl/projects/output_data/DNTFalseNegatives_data.RData'))
hci.dat[grepl('precipitate',wllq_notes), unique(wllq_notes)]
# [1] "neuron count low, group was repeated; Compound observered to precipitate out of solution at this conc in cortical synap 20220112"        
# [2] "neuron count low, group was repeated; Compound observered to precipitate out of solution at this conc"                                   
# [3] "NA; Compound observered to precipitate out of solution at this conc in cortical synap 20220112"                                          
# [4] "well quality set to 0 in xlsx plate map file; Compound observered to precipitate out of solution at this conc in cortical synap 20220112"
# There was 1 assay that was atypical. I'm going to go with just the precipitate notes from cortical synap 20220112
precipitate.observations.tb <- hci.dat[grepl('Compound observered to precipitate out of solution at this conc in cortical synap 20220112',wllq_notes),
                                       .(min_conc_precipitate_observed = min(conc_in_uM, na.rm = T)), by = .(spid)]
precipitate.observations.tb
#                     spid min_conc_precipitate_observed
# 1:                 6-PPD                      9.000000
# 2:         6-PPD Quinone                      2.664388
# 3: 5,5-Diphenylhydantoin                    1000.000000

# Confirm the pseudo-spid names from hci match the acute
setdiff(precipitate.observations.tb$spid, dat$spid) # "empty

# In just the MEA NFA, precipitate also observed at 100, 300, and 1000 uM for 5,5-Diphenylhydantoin
precipitate.observations.tb[spid == '5,5-Diphenylhydantoin', min_conc_precipitate_observed := 100]

# Apply to dat
dat <- merge(dat, precipitate.observations.tb, by = c('spid'), all.x = T)
dat[conc_in_uM >= min_conc_precipitate_observed, .N, by = .(spid, conc_in_uM)]
#                     spid conc_in_uM   N
# 1: 5,5-Diphenylhydantoin  100.00000 261
# 2: 5,5-Diphenylhydantoin  300.00000 261
# 3: 5,5-Diphenylhydantoin 1000.00000 261
# 4:                 6-PPD   10.00000 261
# 5:                 6-PPD   30.00000 261
# 6:         6-PPD Quinone    2.98277 261
# 7:         6-PPD Quinone    8.94832 261
dat[conc_in_uM >= min_conc_precipitate_observed, wllq_notes := paste0(wllq_notes, '; Compound observered to precipitate out of solution at this conc in cortical synap 20220112 or MEA NFA 20210818')]
dat[grepl('precipitate',wllq_notes), .N, by = .(spid, conc_in_uM)]
# looks good!!
dat[, min_conc_precipitate_observed := NULL]
rm(hci.dat)


# > Confirm all wllq assignments ------------------------------------------

# Confirm I identified the correct wells for the precipitate based on compound/conc
dat[grepl('precipitate',wllq_notes), .N, by = .(treatment, wllq, conc_in_uM)]
# treatment wllq conc_in_uM   N
# 1: 5,5-Diphenylhydantoin    1  100.00000 261
# 2: 5,5-Diphenylhydantoin    1  300.00000 259
# 3: 5,5-Diphenylhydantoin    1 1000.00000 259
# 4: 5,5-Diphenylhydantoin    0  300.00000   2
# 5: 5,5-Diphenylhydantoin    0 1000.00000   2
# 6:                 6-PPD    1   10.00000 261
# 7:                 6-PPD    1   30.00000 261
# 8:         6-PPD Quinone    1    2.98277 261
# 9:         6-PPD Quinone    1    8.94832 261


# Confirm all looks okay
library(stringi)
dat[grepl('DIV',acsn), endpoint_type := 'DIV']
dat[grepl('(LDH)|(AB)',acsn), endpoint_type := 'cytotox']
dat[is.na(endpoint_type), endpoint_type := 'AUC']
dat[wllq == 0, .(length(unique(acsn))), by = .(apid, rowi, coli, wllq_notes)]
View(dat[wllq == 0, .(num_AUC_endpoints = length(unique(acsn[endpoint_type == 'AUC'])),
                      num_DIV_endpoints = length(unique(acsn[endpoint_type == 'DIV'])),
                      num_cyto_endpoints = length(unique(acsn[endpoint_type == 'cytotox']))), by = .(apid, rowi, coli, wllq_notes)][order(apid, rowi, coli)])
# The "Contamination" note is only applying to the LDH and CTB endpoints
# because I estimated the missing DIV12 endpoint for the DIV12 and AUC endpionts,
# as is usually done when there is a signal missing DIV


# > View where CO2 turned off, confirm these okay -------------------------

# decided these okay to use
# See investigations/CO2_off_DIV12


# save dat and graphs
dat[, endpoint_type := NULL]
dat[, CO2_off_DIV12 := NULL]
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
