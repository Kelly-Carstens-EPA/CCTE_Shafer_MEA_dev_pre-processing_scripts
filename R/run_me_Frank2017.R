rm(list=ls()) # clear environment
graphics.off() # clear plot history
###################################################################################
# USER INPUT
###################################################################################
dataset_title <- "Frank2017" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probably leave this as true when you first run
save_notes_graphs <- TRUE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName <- "DMSO" # usually DMSO. all compounds other than those listed below should have this vehicle control

spidmap_file <- ""
spid_sheet <- ""

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
# using table1 from published data to match the solvent control in each well to the solvent used for the treatment in the corresponding rows
# Since the control in corresponding row will still be labelled with the treatment name, just with conc==0,
# I can updated the control treatment name by the current treatment name
table1 <- as.data.table(read.xlsx(file.path(root_output_dir,dataset_title,"Table1_CompoundList_dtxsid_updated.xlsx"),sheet = 1))
setnames(table1, old = c("CAS.No.","Compound.name.(abbreviation)"), new = c("CASRN","treatment"))
table1 <- table1[!is.na(treatment)]
table1[treatment == "Acetaminophen", Solvent.used := "Water"] # fixing this, based on lab notebook 20141203

# first, checkign that treatment names are consistent
dat[, .(length(unique(treatment))), by = .(apid, rowi, coli)][V1 != 1] # empty, that's good
setdiff(dat[wllt == "t", unique(treatment)], unique(table1$treatment)) # these will be addressed in updated_treatment_names

# update names in "treatment" col to match names in table1, set original treatments col to "mea_treatment_name"
dat <- update_treatment_names(dat, root_output_dir, dataset_title)

# merge in the solvent info table
dat <- merge(dat, table1[, .(treatment, Solvent.used, CASRN)], by = "treatment", all.x = T)
dat[, any(is.na(Solvent.used))] # FALSE
dat[wllt == "n", `:=`(treatment = Solvent.used,
                      CASRN = "")]
dat[wllt == "n", .N, by = .(treatment)]
# treatment     N
# 1:         DMSO 22431
# 2:        Water  5217
# 3: DMSO/Ethanol   261
# 4:      Ethanol   261
dat[acsn == "CCTE_Shafer_MEA_dev_correlation_coefficient_mean", .(sum(wllt == "n")), by = .(apid)][V1 != 6] # empty, every plate has 6 control wells
dat[is.na(treatment), .N] # 0
dat[, Solvent.used := NULL]

# define the control conc
dat[wllt == "n", conc := 0.001]


# Assign SPIDs ------------------------------------------------------------------
use_spid_files <- file.path(root_output_dir,"Sample IDs",c("EPA_ES202_EPA-Shafer_103_20191218_key.xlsx","EPA_ES204_EPA-Shafer_12_20200117_key.xlsx"))
spidmap <- rbindlist(lapply(use_spid_files, function(filei) {
  tb <- as.data.table(read.xlsx(filei,sheet=1))
  tb$filename <- basename(filei)
  return(tb) }
  ))
head(spidmap)
unique(spidmap$ALIQUOT_CONCENTRATION_UNIT) # all mM? - yes
setnames(spidmap, old = c("PREFERRED_NAME", "ALIQUOT_CONCENTRATION", "EPA_SAMPLE_ID"), new = c("treatment","stock_conc","spid"))
# get rid of the data rows we don't need, to eliminate multiple spids for some treatments
spidmap[treatment == "Chlorpyrifos oxon"] # 5 spids, all from EPA_ES202_EPA-Shafer_103_20191218_key.xlsx
spidmap <- spidmap[!(treatment == "Chlorpyrifos oxon" & ALIQUOT_WELL_ID %in% c(5,6,7,9))] # removing all but the 4th spid (aliquote well ID 8)
spidmap[treatment == "Dexamethasone"] # 2 spids, all from EPA_ES202_EPA-Shafer_103_20191218_key.xlsx
spidmap <- spidmap[!(treatment == "Dexamethasone" & ALIQUOT_WELL_ID == 24)] # we are using the second spid, removign the first
spidmap[treatment == "Methotrexate"] # 2 spids, all from EPA_ES202_EPA-Shafer_103_20191218_key.xlsx
spidmap <- spidmap[!(treatment == "Methotrexate" & ALIQUOT_WELL_ID == 74)] # we are using the second spid, removign the first
spidmap <- spidmap[!(treatment %in% c("Acetaminophen","Glyphosate"))] # these spids we be taken from Shafer_sample_info_to_register_20201110

# get spids from NTP 91 list for TPP and TCEP
spidmap2 <- as.data.table(read.xlsx(file.path(root_output_dir,"Sample IDs","Copy of NTP91_Compounds_4NHEERL_MEA_dev_cg.xlsx"), sheet = "NeuroTox 91 Cmpds"))
setnames(spidmap2, old = c("Chemical.Name","CAS","Conc..(mM)","SPID"), new = c("treatment","CASRN","stock_conc","spid"))
spidmap2 <- spidmap2[treatment %in% c("Triphenyl phosphate","Tris(2-chloroethyl) phosphate")]

# Shafer_42 list for Sodium orthovanadate
spidmap3 <- as.data.table(read.xlsx(file.path(root_output_dir,"Sample IDs","EPA_ES203_EPA-Shafer_42_20200110_key.xlsx"), sheet = 1))
setnames(spidmap3, old = c("PREFERRED_NAME","ALIQUOT_CONCENTRATION","EPA_SAMPLE_ID"), new = c("treatment","stock_conc","spid"))
spidmap3 <- spidmap3[treatment == "Sodium orthovanadate"]

# merge the files
spidmap <- rbindlist(list(spidmap, spidmap2, spidmap3), fill = T)
spidmap[, treatment := as.character(treatment)]
spidmap[, stock_conc := as.numeric(stock_conc)]
spidmap[, expected_stock_conc := ifelse(is.na(TARGET_CONCENTRATION),20,TARGET_CONCENTRATION)] # initialize expected_stock_conc. Usually this is 20mM. Change as needed.
head(spidmap[, .(treatment, spid, stock_conc, expected_stock_conc)])

# the rest of teh control compounds
spidmap4 <- as.data.table(read.xlsx(file.path(root_output_dir, "Sample IDs","Shafer_sample_info_to_register_20201110_afc.xlsx"), sheet = 1))
spidmap4 <- spidmap4[dataset == dataset_title]
spidmap4[, `:=`(expected_stock_conc = stock_conc)]
setnames(spidmap4, old = c("SPID","compound","CAS"), new = c("spid","treatment","CASRN"))
spidmap <- rbind(spidmap[,.(treatment, spid, stock_conc, expected_stock_conc, CASRN)], 
                 spidmap4[, .(treatment, spid, stock_conc, expected_stock_conc, CASRN)])

# confirm unique casn for each treatment name
spidmap[, .(length(unique(treatment))), by = .(CASRN)][V1 != 1] # empty

# see which compounds still have multiple spids. Remove the ones we do not want to use
spidmap[, .(length(unique(spid))), by = .(treatment)][V1 != 1]
# treatment V1
# 1:      Glufosinate-P  2
# 2:         Omeprazole  3
# 3: Boric acid (H3BO3)  2
# Glufosinate-P - neither of these included in current dataset, remove
dat[grepl("mepraz",treatment)] # empty
dat[grepl("[Bb]or",treatment)] # empty
spidmap <- spidmap[!(treatment %in% c("Glufosinate-P","Omeprazole","Boric acid (H3BO3)"))] # none of these included in current data set

# next, check that you can map the compounds by casrn 
setdiff(dat[wllt == "t", unique(CASRN)], unique(spidmap$CASRN)) # empty 

# trick check_and_assign_spids function below to map by casrn instead of treatment name
spidmap[, treatment_name := treatment] # put the updated treatment names in another dummy column
spidmap[, treatment := NULL]
setnames(spidmap, old = c("CASRN"), new = c("treatment"))
dat[, treatment_name := treatment] # put the updated treatment names in another dummy column
dat[, treatment := NULL]
setnames(dat, old = c("CASRN"), new = c("treatment"))
dat[wllt == "n", treatment := treatment_name]

# get rid of duplicate rows (for some compounds, Aliquot well ID might be different, but spid and stock_conc are the same)
spidmap <- unique(spidmap)
spidmap[, .N, by = .(treatment)][N != 1, unique(treatment)] # empty

# assign spids
dat <- check_and_assign_spids(dat, spidmap)
# dat2 <- merge(dat, spidmap[, .(treatment_name, spid)], by = c("spid"), suffixes = c(".dat",".spidmap"), all.x = T)
# unique(dat2[treatment_name.dat != treatment_name.spidmap, .(treatment_name.dat, treatment_name.spidmap)]) # i'm going to leave this as it is
# rm(dat2)
setnames(dat, old = c("treatment","treatment_name"), new = c("CASRN","treatment"))
dat[is.na(spid), .N] # 0
dat[, CASRN := NULL]


# Confirm Conc's ----------------------------------------------------------------
# confirm that the conc's collected from master chem lists and Calc files match
# and that the correct concentration-corrections has been done for each compound

# check if there are multiple conc's assigned to the same well (usually occurs if there are differences between master chem file and calc file)
# Note: in TCPL mc1, the conc's are set to dat[ , conc := signif(conc, 3)]. So it's okay for us to round here.
dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1]
# if any, standardize those before continuing.
problem_comps <- dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1, unique(treatment)]
problem_comps
# character(0). yayay!!

spidmap[, org_stock_conc := stock_conc]

# finally, run this:
source(file.path(scripts.dir, 'confirm_concs.R'))
dat <- confirm_concs(dat, spidmap, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30), update_concs_without_prompt = update_concs_without_prompt)

# I'm a little wary about TPP and TCEP.. will ask Tim again
dat[grepl("TPHP",treatment), unique(apid)]
# [1] "20150128_MW1008-37" "20150128_MW1008-41" "20150128_MW1008-42" "20151125_MW1086-39" "20151125_MW1086-41" "20151125_MW1088-3" 
dat[grepl("TCEP",treatment), unique(apid)]
# [1] "20150128_MW1008-37" "20150128_MW1008-41" "20150128_MW1008-42"

spidmap[org_stock_conc != stock_conc] # empty, cool


# FINAL DATA CHECKS
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Any other plots or things to check?

# Things I noted that I want to check
dat[apid == "20141029_MW1008-41", .N, by = .(acsn)] # looks good
dat[apid == "20141029_MW1008-41" & grepl("DIV5",acsn), .N, by = .(acsn)]
stripchart(rval ~ apid, dat[grepl("CCTE_Shafer_MEA_dev_firing_rate_mean_DIV5",acsn)], vertical = T, method = "jitter", pch = 1, las = 2, cex.axis = 0.6)
dat[treatment %in% c("Saccharin","Spiroxamine","Emamectin benzoate"), .(paste0(sort(unique(conc)),collapse=",")), by = .(treatment)]
# treatment                                                  V1
# 1:          Saccharin                              0.03,0.1,0.3,1,3,10,30
# 2:        Spiroxamine       3e-04,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30
# 3: Emamectin benzoate 1e-04,3e-04,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,3
dat[grepl("20151125",apid), .(paste0(sort(unique(conc)),collapse=",")), by = .(treatment)]
# treatment                                       V1
# 1:                                            DMSO                                    0.001
# 2:                      Triphenyl phosphate (TPHP) 0.0298,0.0995,0.298,0.995,2.98,9.95,29.8
# 3:                                       Saccharin                   0.03,0.1,0.3,1,3,10,30
# 4: Tris (1,3-dichloro-2-propyl) phosphate (TDCIPP)                   0.03,0.1,0.3,1,3,10,30
# 5:                   Tetrabromobisphenol A (TBBPA)                   0.03,0.1,0.3,1,3,10,30
# 6:                                     Spiroxamine      3e-04,0.001,0.003,0.01,0.03,0.1,0.3
# 7:                              Emamectin benzoate    1e-04,3e-04,0.001,0.003,0.01,0.03,0.1
dat[grepl("20151125",apid) & acsn == "CCTE_Shafer_MEA_dev_burst_rate", .N, by = .(treatment)] # 21 for all, 18 for DMSO

# confirm treatments agree with previous data
previous_dat <- as.data.table(read.csv("L:/Lab/NHEERL_MEA/Frank 86 tcpl prep/Intermediate_Output_Round2/Frank86_longfile_withspids_20200501.csv",stringsAsFactors = F))
dat[, plate.SN := sub("^.*_","",apid)]
setnames(previous_dat, old = "apid", new = "plate.SN")
previous_dat <- previous_dat[!grepl("[abc]",plate.SN)] # I did some funky renaming of plate.SN's where a plate was repeated... not goign to mess with that for now
dat_summary <- unique(dat[, .(spid, plate.SN, apid, rowi, coli, wllt)])
pdat_summary <- unique(previous_dat[, .(spid, plate.SN, rowi, coli, wllt)])
test <- merge(dat_summary, pdat_summary, by = c("plate.SN","rowi","coli"), suffixes = c(".new",".prev"))
test[spid.new != spid.prev]
# plate.SN rowi coli     spid.new               apid wllt.new    spid.prev wllt.prev
# 1: MW1045-02    4    2 DMSO/Ethanol 20141231_MW1045-02        n DMSO/ethanol         n
# 2: MW1045-06    4    2 DMSO/Ethanol 20141231_MW1045-06        n DMSO/ethanol         n
# 3: MW1045-09    4    2 DMSO/Ethanol 20141231_MW1045-09        n DMSO/ethanol         n
# yay! That is okay
test[wllt.new != wllt.prev] # empty
dat[, plate.SN := NULL]
rm(previous_dat)

# 20151125
dat[acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12" & rval < 50/60 & wllt == "n", .N, by = "apid"][order(apid)]
# apid N
# 1: 20140910_MW1008-39 2
# 2: 20140924_MW1007-81 1
# 3: 20141231_MW1045-09 1
# 4: 20150128_MW1008-41 1
# 5: 20150805_MW1038-36 2
# 6: 20150805_MW1040-12 1
# 7:  20151125_MW1088-3 3
# 8: 20160601_MW1062-44 1
# okay! So MW1088-3 is the only plate on the table with 3 or more control wells below 50 spikes/min at DIV 12
par(oma = c(2,0,0,0))
stripchart(rval*60 ~ apid, dat[wllq == 1 & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12" & wllt == "n"], vertical = T, method = "jitter", pch = 1, las = 2, cex.axis = 0.6)
stripchart(V1*60 ~ apid, dat[wllq == 1 & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12" & wllt == "n", .(median(rval)),by=.(apid)], vertical = T, method = "jitter", pch = 19, las = 2, cex.axis = 0.6,add = T)
abline(h = 50)
title(main = paste0("Mean Firing Rate in Controls at DIV 12\nfor all apid in ",dataset_title))
legend(x = "topleft", legend = c("control well","median of control wells"), pch = c(1,19), col = "black", bg = "transparent")
par(oma = c(0,0,0,0))
dat[apid == "20151125_MW1088-3" & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12" & wllt == "n", median(rval)*60] # 56.0184

# MW1007-70 no master chem file match
dat[apid == "20140924_MW1007-70", unique(treatment), by = .(rowi)] # looks good, no inconsistencies introduced by acsn

# confirming all wllq stuff
dat[, acsn_cat := ""]
dat[grepl("LDH",acsn), acsn_cat := "LDH"]
dat[grepl("AB",acsn), acsn_cat := "AB"]
dat[acsn_cat == "", acsn_cat := "MEA"]
dat[wllq == 0, .N, by = .(apid, acsn_cat)][order(apid)]
# apid acsn_cat    N
# 1: 20141015_MW1007-107       AB    2
# 2: 20141015_MW1007-107      LDH    2
# 3: 20141015_MW1007-107      MEA  170
# 4: 20141015_MW1007-108       AB    2
# 5: 20141015_MW1007-108      LDH    2
# 6: 20141015_MW1007-108      MEA  170
# 7:  20150128_MW1008-41       AB    2
# 8:  20150128_MW1008-41      LDH    2
# 9:  20150128_MW1008-41      MEA  170
# 10:  20150805_MW1038-36       AB   48
# 11:  20150805_MW1038-36      LDH   48
# 12:  20150805_MW1038-36      MEA 4080
# 13:  20150805_MW1040-11      LDH   48
# 14:   20151125_MW1088-3      LDH   48
dat[, acsn_cat := NULL]

# save dat and graphs
setkey(dat, NULL) # remove all the attributes I inadvertently added to decrease object size
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
