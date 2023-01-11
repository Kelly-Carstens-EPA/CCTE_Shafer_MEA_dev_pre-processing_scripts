rm(list=ls())
graphics.off()
###################################################################################
# USER INPUT
###################################################################################
dataset_title <- "Brown2014" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- TRUE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName <- "DMSO" # usually DMSO. all compounds other than those listed below should have this vehicle control
scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R"
root_output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the dataset_title folder will be created

spidmap_file <- file.path(root_output_dir,"Sample IDs","EPA_ES203_EPA-Shafer_42_20200110_key.xlsx")
spid_sheet <- 1

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
  pdf(file = file.path(root_output_dir, dataset_title, paste0(dataset_title,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")), width = 10, height = 8)
}

# run the main steps
source(file.path(scripts.dir, 'source_steps.R'))

# run tcpl_MEA_dev_AUC
source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
dat <- tcpl_MEA_dev_AUC(basepath = file.path(root_output_dir,dataset_title), dataset_title)

# confirming consistent treatment names
dat[, .N, by = .(treatment)] # yep, no inconsistent naming.


# change untreated wells to Control Treatment ------------------------------------
dat[wllt == "n", treatment := default_ControlTreatmentName]
# Acetaminophen, Domic Acid, and Sodium Orthovandate should be "water" for all cultures
water_trts <- c("Acetaminophen","Domoic Acid","Sodium Orthovanadate","Glyphosate")

# 20140205
dat[grepl("^20140205",apid) & treatment %in% water_trts, .(rowi = unique(rowi)), by = "treatment"]
# treatment rowi
# 1:        Acetaminophen    1
# 2:          Domoic Acid    3
# 3: Sodium Orthovanadate    6
water_trt_rows <- dat[grepl("^20140205",apid) & treatment %in% water_trts, unique(rowi)]
dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "20140205", control_rowi = water_trt_rows)

# 20140212
dat[grepl("^20140212",apid) & treatment %in% water_trts, .(rowi = unique(rowi)), by = "treatment"]
# treatment rowi
# 1:        Acetaminophen    1
# 2:          Domoic Acid    3
# 3: Sodium Orthovanadate    6
water_trt_rows <- dat[grepl("^20140212",apid) & treatment %in% water_trts, unique(rowi)]
dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "20140212", control_rowi = water_trt_rows)
# Control treatment will be updated to Water for the following wells:
#   apid treatment rowi coli conc
# 1: 20140212_MW1007-27      DMSO    1    1    0
# 2: 20140212_MW1007-27      DMSO    1    8    0
# 3: 20140212_MW1007-27      DMSO    3    1    0
# 4: 20140212_MW1007-27      DMSO    3    8    0
# 5: 20140212_MW1007-27      DMSO    6    1    0
# 6: 20140212_MW1007-27      DMSO    6    8    0

# 20140402
dat[grepl("^20140402",apid) & treatment %in% water_trts, .(rowi = unique(rowi)), by = "treatment"]
# treatment rowi
# 1:        Acetaminophen    1
# 2:          Domoic Acid    3
# 3: Sodium Orthovanadate    6
water_trt_rows <- dat[grepl("^20140402",apid) & treatment %in% water_trts, unique(rowi)]
dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "20140402", control_rowi = water_trt_rows)

# 20140423
dat[grepl("^20140423",apid) & treatment %in% water_trts, .(rowi = unique(rowi)), by = "treatment"]
# treatment rowi
# 1:        Acetaminophen    1
# 2:          Domoic Acid    3
# 3: Sodium Orthovanadate    6
water_trt_rows <- dat[grepl("^20140423",apid) & treatment %in% water_trts, unique(rowi)]
dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "20140423", control_rowi = water_trt_rows)

# 20140716
dat[grepl("^20140716",apid) & treatment %in% water_trts, .(rowi = unique(rowi)), by = "treatment"]
# treatment rowi
# 1:           Glyphosate    1
# 2:           Glyphosate    4
# 3: Sodium Orthovanadate    3
# 4: Sodium Orthovanadate    6
water_trt_rows <- dat[grepl("^20140716",apid) & treatment %in% water_trts, unique(rowi)]
dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "20140716", control_rowi = water_trt_rows)

# 20140730
dat[grepl("^20140730",apid) & treatment %in% water_trts, .(rowi = unique(rowi)), by = "treatment"]
# treatment rowi
# 1:           Glyphosate    2
# 2:           Glyphosate    5
# 3: Sodium Orthovanadate    3
# 4: Sodium Orthovanadate    6
water_trt_rows <- dat[grepl("^20140730",apid) & treatment %in% water_trts, unique(rowi)]
dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "20140730", control_rowi = water_trt_rows)

# from lab notebook:
# "accidentally used sterile water as solvent control for Glyphosate" (20140730)
dat[grepl("^20140730",apid) & treatment == "Glyphosate", .N, by = "rowi"]
# rowi   N
# 1:    2 602
# 2:    5 602
dat[grepl("^20140730",apid) & wllt == "n" & rowi %in% c(2,5), treatment := "Water"]

# we are not including the data for Glyphosate, as per the published data file
dat[treatment == "Glyphosate", .N, by = .(rowi, apid)]
# rowi                apid   N
# 1:    1  20140716_MW1007-26 602
# 2:    4  20140716_MW1007-26 602
# 3:    2 20140730_MW1007-104 602
# 4:    5 20140730_MW1007-104 602
dat <- dat[treatment != "Glyphosate"]

dat[wllt == "n", .N, by = c("treatment","rowi")]
# treatment rowi   N
# 1:     Water    1 688
# 2:      DMSO    2 860
# 3:     Water    3 688
# 4:      DMSO    4 946
# 5:      DMSO    5 860
# 6:     Water    6 688
# 7:      DMSO    1 258
# 8:      DMSO    3 172
# 9:      DMSO    6 172
# 10:     Water    2  86
# 11:     Water    5  86
dat[wllt == "n", conc := 0.001] # All wells have 0.1% solvent

# Assign SPIDs ------------------------------------------------------------------
spidmap <- as.data.table(read.xlsx(spidmap_file, sheet = spid_sheet))
head(spidmap)
unique(spidmap$TARGET_CONCENTRATION_UNIT) # mM
setnames(spidmap, old = c("PREFERRED_NAME","ALIQUOT_CONCENTRATION","EPA_SAMPLE_ID","TARGET_CONCENTRATION"), new = c("treatment","stock_conc","spid","expected_stock_conc"))
spidmap <- spidmap[treatment == "Sodium orthovanadate"]
spidmap[, treatment := as.character(treatment)]
spidmap[, stock_conc := as.numeric(stock_conc)]
spidmap

spidmap2 <- as.data.table(read.xlsx(file.path(root_output_dir, "Sample IDs","Shafer_sample_info_to_register_20201110_afc.xlsx"), sheet = 1))
spidmap2 <- spidmap2[dataset == dataset_title]
spidmap2[, `:=`(expected_stock_conc = stock_conc)]
setnames(spidmap2, old = c("SPID","compound"), new = c("spid","treatment"))

# for Bis 1, must differentiate the treatments for each spid by date tested
spid_date_map <- spidmap2[, .(date_range_tested = unlist(strsplit(dates_tested,split=" - "))), by = .(spid, treatment)]
spidmap2 <- merge(spidmap2, spid_date_map[treatment == "Bisindolylmaleimide 1"], all.x = T)
spidmap2[treatment == "Bisindolylmaleimide 1", treatment := paste0(treatment,"_",date_range_tested)]
usecols <- c("spid","treatment","stock_conc","expected_stock_conc")
spidmap <- rbind(spidmap[, ..usecols], spidmap2[, ..usecols])

# Loperamide spid comes from Shafer_103 (HTP_LOG)
spidmap3 <- as.data.table(read.xlsx(file.path(root_output_dir, "Sample IDs","EPA_ES202_EPA-Shafer_103_20191218_key.xlsx"), sheet = 1))
setnames(spidmap3, old = c("PREFERRED_NAME","ALIQUOT_CONCENTRATION","EPA_SAMPLE_ID","TARGET_CONCENTRATION"), new = c("treatment","stock_conc","spid","expected_stock_conc"))
spidmap3 <- spidmap3[treatment == "4-(4-Chlorophenyl)-4-hydroxy-N,N-dimethyl-alpha,alpha-diphenylpiperidine-1-butyramide monohydrochloride"]
spidmap <- rbind(spidmap[, ..usecols], spidmap3[, ..usecols])

# rename any compounds, if needed
dat <- update_treatment_names(dat, root_output_dir, dataset_title)
dat[, .N, by = .(treatment, mea_treatment_name)]

# update the treatment name to include the culture date for Bis 1
dat[treatment == "Bisindolylmaleimide 1", treatment := paste0(treatment,"_",sub("_.*$","",apid))]

spidmap[treatment %in% unique(dat$treatment)][order(treatment)]
# spid                                                                                               treatment stock_conc expected_stock_conc
# 1:      EX000411 4-(4-Chlorophenyl)-4-hydroxy-N,N-dimethyl-alpha,alpha-diphenylpiperidine-1-butyramide monohydrochloride         20                  20
# 2: MEA20201109A6                                                                                           Acetaminophen        100                 100
# 3:      EX000475                                                                          Bisindolylmaleimide 1_20140205         10                  10
# 4:      EX000475                                                                          Bisindolylmaleimide 1_20140212         10                  10
# 5: MEA20201109A1                                                                          Bisindolylmaleimide 1_20140402         10                  10
# 6: MEA20201109A1                                                                          Bisindolylmaleimide 1_20140423         10                  10
# 7: MEA20201109A2                                                                          Bisindolylmaleimide 1_20140716         10                  10
# 8: MEA20201109A2                                                                          Bisindolylmaleimide 1_20140730         10                  10
# 9:      EX000498                                                                                              Mevastatin         30                  30
# 10:      EX000499                                                                                    Sodium orthovanadate        100                 100

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
# "Acetaminophen"        "L-Domoic acid"        "Sodium orthovanadate"
dat[treatment %in% problem_comps, .(paste0(sort(unique(signif(conc,3))),collapse=",")), by = .(treatment)]

# get a summary of the treatments by srcf
summary_dat <- dat[treatment %in% problem_comps, .(conc_shown = unique(conc)), by = .(apid, rowi, coli, treatment, srcf)]
summary_dat[, conc_round := signif(conc_shown, 1)]
summary_dat[, conc_source := ""]
summary_dat[grepl("(Calc)|(Summary)",srcf), conc_source := "Calc"]
summary_dat[grepl("AUC",srcf), conc_source := "AUC"]
summary_dat[grepl("DIV",srcf), conc_source := "DIV"]
summary_wide <- dcast(summary_dat, apid + treatment + rowi + coli ~ conc_source, value.var = "conc_shown")

summary_wide[round(AUC,2) != round(Calc,2)]
summary_wide[AUC != Calc, unique(apid)]
# [1] "20140212_MW1007-27" "20140402_MW1007-27" "20140423_MW1007-38"
# I already verified all of the conc's in the maestro exp log file
# So I am going to take those conc's over the conc's in teh Calc/Summary data

summary_wide[AUC != DIV] # just confirming these are all the same

summary_wide[, use_conc := AUC]
dat <- merge(dat, summary_wide[, .(apid, treatment, rowi, coli, use_conc)], all.x = T, by = c("apid","treatment","rowi","coli"))
dat[!is.na(use_conc), conc := use_conc]
dat[, use_conc := NULL]

# just confirming this, since it was flipped in master chem file originally
dat[grepl("20140716",apid) & grepl("Bis",treatment), .(unique(conc)), by = .(rowi, coli)]
# good, B8 is 10uM, B7 is not present, so must be control
# but wait... E7 is 10 and E8 is 10... E8 should be a control well.
# I will fix this in prepared data file, AUC, and the source h5file.
dat[grepl("20140716",apid) & grepl("Bis",treatment) & grepl("AB",acsn), .(unique(conc)), by = .(rowi, coli)] # the summary file and cytotox data is correct though
# this is fixed now

# confirming this is empty now
dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1]

# confirm conc ranges included with the stated conc ranges in paper
dat[wllq == 1, .(concs_tested = paste0(sort(unique(conc)),collapse=",")), by = .(treatment)][order(treatment)]
# treatment                concs_tested
# 1: 4-(4-Chlorophenyl)-4-hydroxy-N,N-dimethyl-alpha,alpha-diphenylpiperidine-1-butyramide monohydrochloride           0.1,0.3,1,3,10,30
# 2:                                                                                           Acetaminophen           0.1,0.3,1,3,10,30
# 3:                                                                          Bisindolylmaleimide 1_20140205         0.03,0.1,0.3,1,3,10
# 4:                                                                          Bisindolylmaleimide 1_20140212         0.03,0.1,0.3,1,3,10
# 5:                                                                          Bisindolylmaleimide 1_20140402         0.03,0.1,0.3,1,3,10
# 6:                                                                          Bisindolylmaleimide 1_20140423         0.03,0.1,0.3,1,3,10
# 7:                                                                          Bisindolylmaleimide 1_20140716         0.03,0.1,0.3,1,3,10
# 8:                                                                          Bisindolylmaleimide 1_20140730         0.03,0.1,0.3,1,3,10
# 9:                                                                                                    DMSO                       0.001
# 10:                                                                                             Domoic Acid   0.003,0.01,0.03,0.1,0.3,1
# 11:                                                                                              Mevastatin           0.1,0.3,1,3,10,30
# 12:                                                                                    Sodium orthovanadate 0.01,0.03,0.1,0.3,1,3,10,30
# 13:                                                                                                   Water                       0.001
# this agrees with paper

# finally, run this:
# source(file.path(scripts.dir, 'confirm_concs.R'))
# dat <- confirm_concs(dat, spidmap, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30), update_concs_without_prompt = update_concs_without_prompt)
# I am skipping this step because the stock_conc's in the db are...
# - for Bis 1, Mevastatin, and Domoic Acid, the stkc's in the db are all listed as 100, but I think that what I have in spidmap file is more correct
# - for Loperamide, and Sodium orthovanadate, the stkc's in db are 20, 100, and 100, resp - so no conc-correction would be needed
# - Acetamionphen - stock conc is purely provided by me in spidmap file
# con <- dbConnect(drv = RMySQL::MySQL(), user = Sys.getenv('INVITRODB_USER_RO'), pass = Sys.getenv('INVITRODB_PASS_RO'), dbname='invitrodb',host = Sys.getenv('INVITRODB_HOST_RO'))
# query_term <- paste0("SELECT * FROM sample WHERE spid IN('",paste(spidmap[!is.na(spid),unique(spid)],collapse="','",sep=""),"');")
# sample_info <- dbGetQuery(con, query_term)
# dbDisconnect(con)
# spidmap <- merge(spidmap, sample_info, by = "spid", all.x = T)
# spidmap[, .(treatment, spid, expected_stock_conc, stkc)]
# treatment                                                                                                           spid expected_stock_conc stkc
# 1: 4-(4-Chlorophenyl)-4-hydroxy-N,N-dimethyl-alpha,alpha-diphenylpiperidine-1-butyramide monohydrochloride      EX000411                  20   20
# 2:                                                                          Bisindolylmaleimide 1_20140205      EX000475                  10  100
# 3:                                                                          Bisindolylmaleimide 1_20140212      EX000475                  10  100
# 4:                                                                                             Domoic Acid      EX000487                  10  100
# 5:                                                                                              Mevastatin      EX000498                  30  100
# 6:                                                                                    Sodium orthovanadate      EX000499                 100  100
# 7:                                                                          Bisindolylmaleimide 1_20140402 MEA20201109A1                  10   NA
# 8:                                                                          Bisindolylmaleimide 1_20140423 MEA20201109A1                  10   NA
# 9:                                                                          Bisindolylmaleimide 1_20140716 MEA20201109A2                  10   NA
# 10:                                                                          Bisindolylmaleimide 1_20140730 MEA20201109A2                  10   NA
# 11:                                                                                           Acetaminophen MEA20201109A6                 100   NA
rm(list=c("summary_dat","summary_wide"))


# FINAL DATA CHECKS
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Any other plots or things to check?
dat[wllt == "n", .N, by = .(apid, acsn)][N != 12, unique(apid)]
# "20140716_MW1007-26"  "20140730_MW1007-104"
dat[apid %in% c("20140716_MW1007-26", "20140730_MW1007-104") & wllt == "n", .N, by = .(apid, acsn)][N != 8]
# Empty data.table (0 rows and 3 cols): apid,acsn,N

# correlation plot with the published data
pub_dat <- as.data.table(read.csv(file.path(root_output_dir,dataset_title,"Published Downloaded Data","Final_Data_Set_SA1_DNT_Paper1 (2)(updated)_CF.csv"),stringsAsFactors = F))
div_dat <- as.data.table(read.csv(file.path(root_output_dir,dataset_title,"output",paste0(dataset_title,"_parameters_by_DIV.csv")),stringsAsFactors = F))
comb_dat <- merge(pub_dat, div_dat, by = c("date","Plate.SN","well","DIV"), suffixes = c(".org",".new"), all = T)

# check basic well ID stuff
comb_dat[treatment != trt] # empty
comb_dat[dose.new != dose.org] # just Acetaminophen from well A6 and A7 20140212. I switched these 2 doses according to the lab notebook, looks like they did not.
comb_dat[is.na(date)]
setdiff(div_dat$date, pub_dat$date) # empty
setdiff(pub_dat$date, div_dat$date) # empty

# compare rval's
plot(comb_dat[, .(meanfiringrate.new, meanfiringrate.org)], main ="Correlation plot of the Published MFR values vs current, for all DIV")
abline(0,1)
# not too bad. I am really looking mroe for anomalies than for a general pattern of difference -
# I know that there have been slight changes in the script, and that's okay.
# more concerned if the wrong data values got mapped to some wells
comb_dat[abs(meanfiringrate.new - meanfiringrate.org) > 1] # 2 instances

plot(comb_dat[, .(nAE.new, nAE.org)], main = "Correlation plot of the Published nAE values vs current, for all DIV")
abline(0,1)
# not too bad. I am really looking mroe for anomalies than for a general pattern of difference -
# I know that there have been slight changes in the script, and that's okay.
# more concerned if the wrong data values got mapped to some wells
comb_dat[abs(nAE.new - nAE.org) > 2, .(date, Plate.SN, well, DIV, trt, dose.org, dose.new, nAE.org, nAE.new)] # 2 instances

# I think I just want to "feel" that the difference is due to just spike list file chopping
comb_dat[abs(meanfiringrate.new - meanfiringrate.org) > 1, .(date, Plate.SN, well, DIV, trt, dose.org, dose.new, meanfiringrate.org, meanfiringrate.new)]
# date  Plate.SN well DIV        trt dose.org dose.new meanfiringrate.org meanfiringrate.new
# 1: 20140205 MW1007-26   E6   9 Mevastatin       10       10           6.558445           4.603398
# 2: 20140212 MW1007-27   E5   9 Mevastatin        3        3           7.555179           3.642779

# library(rhdf5)
# new_h5 <- h5read(list.files(file.path(root_output_dir,dataset_title,"h5Files"), pattern = "20140212_MW1007-27_09", full.names= T), name = "/")
# str(new_h5)
# is.ordered(new_h5$spikes) # FALSE
# beg_time <- min(new_h5$spikes)
# end_time <- max(new_h5$spikes)
# nspikes <- length(new_h5$spikes)
# mfr_by_channel <- new_h5$sCount / (end_time - beg_time)
# grep("E5",new_h5$names, val = T) # "E5_12" "E5_14" "E5_21" "E5_31"
# (welli_channels <- mfr_by_channel[grepl("E5",new_h5$names)])
# # 7.124624908 0.132818858 0.002213648 0.081904962
# mean(welli_channels[welli_channels > 5/60]) # 3.642779. Yay, that's what is in the actual data!!
# 
# org_h5 <- h5read("L:/Lab/NHEERL_MEA/PIP3 - Project/Data/Specific Aim 1/20140212 Ontogeny/h5Files/ON_20140212_MW1007-27_DIV09_001.h5", name = "/")
# str(org_h5)
# is.ordered(org_h5$spikes) # FALSE
# beg_time <- min(org_h5$spikes)
# end_time <- max(org_h5$spikes)
# mfr_by_channel <- org_h5$sCount / (end_time - beg_time)
# grep("E5",org_h5$names, val = T) # "E5_12" "E5_14" "E5_21" "E5_31" "E5_42"
# (welli_channels <- mfr_by_channel[grepl("E5",org_h5$names)])
# # 7.555179373 0.057554838 0.006640943 0.070836724 0.004427295
# mean(welli_channels[welli_channels > 5/60]) # 7.555179. Yep!! this is the value they got before
# 
# # 5/60 = 0.08333333
# # what really drives the differences is well E5_14
# which(new_h5$names == "E5_14") # 446
# spike_index_start <- sum(new_h5$sCount[1:445])
# range(new_h5$spikes[spike_index_start:(spike_index_start + new_h5$sCount[446])]) # 8.39064 899.96120
# 
# which(org_h5$names == "E5_14") # 443
# spike_index_start <- sum(org_h5$sCount[1:442])
# range(org_h5$spikes[spike_index_start:(spike_index_start + org_h5$sCount[443])]) # 216.8342 1067.3945
# # huh, so they actually just shifted by about 160 seconds...
# min(org_h5$spikes) # 164.0084
# max(org_h5$spikes) # 1067.495
# org_h5$sCount[443] # 52
# new_h5$sCount[446] # 120
# 
# # huh, so I guess the difference is in the recording slice used.
# # But again, this is the biggest anamoly that we see
# # for most wells, there is hardly any difference.
# # I am content to more forward.
# rm(list = c("new_h5","org_h5","comb_dat","pub_dat","div_dat"))

# quick check of DIV 9 estimated values
par(mfrow = c(4, 4))
wells <- div_dat[date == "20140730", unique(well)]
for (welli in wells) {
  plotdat <- div_dat[date == "20140730" & well == welli, .(DIV, meanfiringrate)][order(DIV)]
  plot(plotdat, type = "o", ylim = div_dat[date == "20140730", range(meanfiringrate)])
  title(main = paste0(div_dat[date == "20140730" & well == welli, unique(treatment)]," ",
                      div_dat[date == "20140730" & well == welli, unique(dose)],"uM",
                      "\n", welli, " 20140730 DIV 9 estimation"))
}

par(mfrow = c(4, 4))
wells <- div_dat[date == "20140716", unique(well)]
for (welli in wells) {
  plotdat <- div_dat[date == "20140716" & well == welli, .(DIV, meanfiringrate)][order(DIV)]
  plot(plotdat, type = "o", ylim = div_dat[date == "20140716", range(meanfiringrate)])
  title(main = paste0(div_dat[date == "20140716" & well == welli, unique(treatment)]," ",
                      div_dat[date == "20140716" & well == welli, unique(dose)],"uM",
                      "\n", welli, " 20140716 DIV 9 comparison"))
}

# save dat and graphs
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")


# # EXTRA DATA CHECKS:
# # determining which copy of the spike list files to use, and if they are different at all
# filei <- "L:/Lab/NHEERL_MEA/PIP3 - Project/Data/Specific Aim 1/Regenerated spikelist files SA1 compounds/20140423/ON_20140423_MW1007-38_02_00(000)_Spike Detector (8 x STD)(000)_spike_list.csv"
# data.raw<-read.csv(filei,header=F,colClasses=c("NULL", "NULL", "character","character","character")) # make electrode column char, not factor
# data.info <- read.csv(filei,header=F,colClasses=c("character","character"), nrows = 100)
# head(data.raw)
# data.info
# data.raw # this is it??
# # V3        V4                   V5
# # 1      Time (s) Electrode        Amplitude(mV)
# # 2    56.2676800     A8_23 0.028872283887509825
# # 3   126.3688800     A8_23 0.030162792900087505
# # 4    730.171200     C2_24 0.021529744499527877
# # 5    730.171200     C2_44 0.023524730659420935
# # 6   730.1712800     C2_33 0.023111509512467321
# 
# # well, that was just DIV 2
# filei <- "L:/Lab/NHEERL_MEA/PIP3 - Project/Data/Specific Aim 1/Regenerated spikelist files SA1 compounds/20140423/ON_20140423_MW1007-38_05_00(000)_Spike Detector (8 x STD)(000)_spike_list.csv"
# data.raw<-read.csv(filei,header=F,colClasses=c("NULL", "NULL", "character","character","character")) # make electrode column char, not factor
# data.info <- read.csv(filei,header=F,colClasses=c("character","character"), nrows = 100)
# data.info
# data.raw # okay, this looks good
# data.raw_0 <- data.raw
# 
# # comparing with the other file
# filei <- "L:/Lab/NHEERL_MEA/PIP3 - Project/Data/Specific Aim 1/Regenerated spikelist files SA1 compounds/20140423_1/ON_20140423_MW1007-38_05_00(000)_Spike Detector (8 x STD)(000)_spike_list.csv"
# data.raw<-read.csv(filei,header=F,colClasses=c("NULL", "NULL", "character","character","character")) # make electrode column char, not factor
# data.info <- read.csv(filei,header=F,colClasses=c("character","character"), nrows = 100)
# data.info
# data.raw # okay, this looks good
# data.raw_1 <- data.raw
# 
# # comparison
# nrow(data.raw_0)
# # [1] 9914
# nrow(data.raw_1)
# # [1] 9843
# tail(data.raw_0) # lots of empty rows
# tail(data.raw_1)
# data.raw_0 <- data.raw_0[data.raw_0$V3 != "",] # chopping off several rows at the end
# tail(data.raw_0) # okay, still have well summary info here. snippet:
# # 9841 913.8309600   D3_13 0.061697795920477749
# # 9842  914.167600   A3_43 0.028635573871101324
# # 9843  914.215600   F2_22 0.074398127943828149
# # 9845          A2      A3                   A4
# # 9846          A7      A8                   B1
# # 9847          B4      B5                   B6
# # 9848          C1      C2                   C3
# data.raw_0 <- data.raw_0[1:9843,]
# all.equal(data.raw_0, data.raw_1)
# # TRUE!
# # okay!
# # so now the question is... which do I use? Do I need to check every other file now?
# 
