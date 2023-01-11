rm(list=ls()) # clear environment
graphics.off() # clear plot history
###################################################################################
# USER INPUT
###################################################################################
dataset_title <- "PFAS2019" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- TRUE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName <- "DMSO" # all compounds other than those listed below should have this vehicle control

spidmap_file <- ""
spid_sheet <- ""

scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R" # updated to the folder where the scripts are located
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

# check if trt names consistent between calc and master chem list annotations
dat[, treatment_srcf := treatment]
check.well.trts <- dat[, .(trt_labs_per_well = length(unique(treatment))), by = .(apid, rowi, coli)][trt_labs_per_well > 1, .(apid, rowi, coli)]
setkey(dat, apid, rowi, coli)

# Which cases is this due to a control well?
dat[J(check.well.trts)][, .N, by = .(conc, wllt)]
# conc wllt     N
# 1:    0    n 10179
# all deal with control wells
dat[J(check.well.trts)][, .N, by = .(treatment)]

auc.trt.tb <- dat[J(check.well.trts)][srcf == 'PFAS2019_AUC.csv', unique(.SD), .SDcols = c('apid','treatment','conc','wllt','rowi','coli','srcf')]
calc.trt.tb <- dat[J(check.well.trts)][!(srcf %in% paste0(dataset_title,c('_AUC.csv','_parameters_by_DIV.csv'))), unique(.SD), .SDcols = c('apid','treatment','conc','wllt','rowi','coli','srcf')]
comp.trts <- merge(auc.trt.tb, calc.trt.tb, by = c('apid','rowi','coli'), all = T, suffixes = c('.mcl','.calc'))
comp.trts[conc.mcl != conc.calc] # empty
comp.trts[wllt.mcl != wllt.calc] # empty
View(comp.trts[, .(apid, rowi, coli, treatment.mcl, treatment.calc, conc = conc.mcl, srcf.calc)])
# Default to the treatment in MCL (this lists DMSO where conc is 0)
# Exceptions:
# - 20210127_MW75-562 - use treatment.calc

# Update treatments ------------------------------
# Default to the treatment listed in MCL file
dat[apid != '20210127_MW75-5620', treatment := unique(treatment_srcf[srcf == 'PFAS2019_AUC.csv']), by = .(apid, rowi, coli)]
# For this plate, use treatments listed in calc file
dat[apid == '20210127_MW75-5620', treatment := unique(treatment_srcf[srcf == '20210127_NFA_Pos_Ctrl+Chlorpyrifos oxon_ Calculations.xlsx']), by = .(apid, rowi, coli)]

# Wherever Media was tested, the well labelled 'DMSO' in same row is actually just Media
dat[, media_rowi := unique(rowi[treatment == 'Media']), by = .(apid)]
dat[rowi == media_rowi & treatment != 'Media', .N, by = .(treatment, coli, rowi, srcf, apid)]
#    treatment coli rowi                                                       srcf               apid  N
# 1:      DMSO    2    3 20210127_NFA_Pos_Ctrl+Chlorpyrifos oxon_ Calculations.xlsx 20210127_MW75-5620  2
# 2:      DMSO    2    3                             PFAS2019_parameters_by_DIV.csv 20210127_MW75-5620 68
# 3:      DMSO    2    3                                           PFAS2019_AUC.csv 20210127_MW75-5620 17
# 4:      DMSO    2    6               20210331_NFA_PFAS_Group_8_ Calculations.xlsx 20210331_MW75-8102  2
# 5:      DMSO    2    6                             PFAS2019_parameters_by_DIV.csv 20210331_MW75-8102 68
# 6:      DMSO    2    6                                           PFAS2019_AUC.csv 20210331_MW75-8102 17
# 7:      DMSO    2    4               20210331_NFA_PFAS_Group_8_ Calculations.xlsx 20210331_MW75-8103  2
# 8:      DMSO    2    4                             PFAS2019_parameters_by_DIV.csv 20210331_MW75-8103 68
# 9:      DMSO    2    4                                           PFAS2019_AUC.csv 20210331_MW75-8103 17
# 10:      DMSO    2    5               20210331_NFA_PFAS_Group_8_ Calculations.xlsx 20210331_MW75-8104  2
# 11:      DMSO    2    5                             PFAS2019_parameters_by_DIV.csv 20210331_MW75-8104 68
# 12:      DMSO    2    5                                           PFAS2019_AUC.csv 20210331_MW75-8104 17
# have confirmed these wells should be Media for 20210127_MW75-5620, 
# will confirm teh others with Seline
dat[rowi == media_rowi, treatment := 'Media'] # replace DMSO with Media

# View changes
dat[rowi == media_rowi, .N, by = .(treatment, apid, rowi, coli)]
# all Media, from 4 plates*8 cols = 32 cases
dat[, media_rowi := NULL]
dat[treatment != treatment_srcf, .N, by = .(apid, rowi, coli, treatment, treatment_srcf, conc, wllt)]
# 138 cases

# change untreated wells to Control Treatment ------------------------------------
# dat[wllt == "n", treatment := default_ControlTreatmentName]
# Manually update other wells where control treatment is not the default, or use teh function below
# dat <- update_control_well_treatment(dat, control_compound = "Water",culture_date = "")

# Confirm have appropriate treatment label wherever wllt == 'n'
dat[wllt == 'n', .N, by = .(conc, treatment)]
#    conc treatment    N
# 1:    0      DMSO 9570
# 2:    0     Media 2784
# 3:    0  3360 A04  261
# 4:    0  3360 A07  261
# 5:    0  3360 A11  261
# 6:    0  3360 B03  261
# 7:    0  3360 B04  261
# 8:    0  3360 B12  261
# 9:    0  3360 C09  261
# 10:    0  3360 C11  261
# 11:    0  3360 D01  261
# 12:    0  3360 D03  261
# 13:    0  3360 D07  261
# 14:    0  3360 D12  261
# is the solvent treatment not present in the calc file for these plates?
dat[wllt == 'n' & grepl('Calc',srcf) & !(treatment %in% c('DMSO','Media')), .N, by = .(apid)]
# apid  N
# 1: 20210303_MW75-5813 12
# 2: 20210303_MW75-5814 12
# 3: 20210303_MW75-5815 12
# 4: 20210303_MW75-5816 12
# 5: 20210303_MW75-5817 12
# 6: 20210303_MW75-5818 12
# For 20210303 G1 (MW75-5813, MW75-5814, MW75-5815), I confirmed in '20210303_NFA_PFAS_Group 1_ Calculations.xlsx' sheet "Dosing Prep"
# that DMSO was placed in every well in col 2 of dosing plate
# For 20210303 G2 (MW75-5816, MW75-5817, MW75-5818), I confirmed in '20210303_NFA_PFAS_Group 2_ Calculations.xlsx' sheet "Dosing Prep"
# that DMSO was placed in every well in col 2 of dosing plate
update.apid.controls <- dat[wllt == 'n' & grepl('Calc',srcf) & !(treatment %in% c('DMSO','Media')), .N, by = .(apid)][, unique(apid)]
dat[wllt == 'n' & conc == 0 & apid %in% update.apid.controls, treatment := 'DMSO']
dat[wllt == 'n', .N, by = .(conc, treatment)]
# conc treatment     N
# 1:    0      DMSO 12702
# 2:    0     Media  2784
# All DMSO and Media now

# DMSO only should be wllt n here (not Media)
dat[treatment == 'Media', .N, by = .(wllt)] # all currently n
dat[treatment == 'Media', `:=`(wllt = 'b')]
# For plates with only 5 DMSO controls -> keep the Media well in the same column as the other controls as wllt == 'n'
update.apids <- dat[wllt == 'n', .(length(unique(paste0(rowi,coli)))), by = .(apid)][V1 < 6, unique(apid)]
dat[apid %in% update.apids & wllt == 'n', .(dmso_control_col = unique(coli)), by = .(apid)]
#                  apid dmso_control_col
# 1: 20210127_MW75-5620                2
# 2: 20210331_MW75-8102                2
# 3: 20210331_MW75-8103                2
# 4: 20210331_MW75-8104                2
# 3 for all
dat[apid %in% update.apids & coli == 2 & treatment == 'Media', `:=`(wllt = 'n')]
dat[wllt == 'n', .(length(unique(paste0(rowi,coli)))), by = .(apid)]
# all have exactly 6 controls now

# Any unexpected cases where conc == 0?
dat[conc %in% c(0) | is.na(conc), .N, by = .(treatment, conc, wllt)]
# treatment conc wllt     N
# 1:      DMSO    0    n 12702
# 2:     Media    0    b  2436
# 3:     Media    0    n   348
# all as expected

# update anythign else?
dat[, .N, by = .(wllt, treatment)]

# Set the control well concentration. Adjust as needed
dat[treatment == 'DMSO', conc := 0.001]
dat[treatment == 'Media', conc := NA_real_]



# assign sample ID's -------------------------------------------------------------
# This is how I obtained the SPIDs for the HCI SPS Data. We'll see if this covers the same chem
# (and just confirm with someone that these do indeed correspond t the samples tested)
spidmap <- as.data.table(read.xlsx("L:/Lab/NHEERL_Mundy/Project - PFAS 2019/Supporting Doc/EPA_27864_EPA-Shafer_134_20191001_key.xlsx", sheet = 1))
spidmap[, short_rackplate := sub("SRACK0","",RACKPLATE_BARCODE)]
unique(spidmap$short_rackplate) # "3360" "3361"
spidmap2 <- as.data.table(read.xlsx("L:/Lab/NHEERL_Mundy/Project - PFAS 2019/Supporting Doc/EPA_29885_EPA-Shafer_36_20191112_key.xlsx", sheet = 1))
spidmap2[, short_rackplate := sub("SRACK0","",RACKPLATE_BARCODE)]
spidmap2[, expected_stock_conc := CONCENTRATION]
usecols <- intersect(names(spidmap), names(spidmap2))
spidmap <- rbind(spidmap[,usecols,with=F], spidmap2[,usecols,with=F])
spidmap[, treatment := paste(short_rackplate, WELL_POSITION, sep = " ")]
unique(spidmap$CONCENTRATION_UNIT) # all mM
spidmap[, expected_stock_conc := 30] # perhaps expected isn't the word... but known when doing dilutions adn calculating conc
spidmap3 <- as.data.table(read.xlsx(file.path(root_output_dir,'Sample IDs/EPA_ES209_EPA-Shafer_4_20210504_key.xlsx'), sheet = 1))
spidmap3[, treatment := PREFERRED_NAME]
# I'm not sure why Chris entered 65 as the Target concentration for these 4 samples...
# based on conc's tested and calc files, I am confident that the expected concentrations are as follows:
spidmap3[treatment %in% c('Bisphenol A','Acetaminophen'), expected_stock_conc := 30]
# For Chlorpyrifos and Chlorpyrifos oxon - looking at the calc files, looks like actual and expected stock conc's were 30
# But, these are registered in the sample table as 100mM for the actual conc.
# I want expected conc == stock conc, because the people who made the diltuions made the stock conc equal to the expected conc.
# So I'm going to enter 100 for the expected stock con
spidmap3[treatment %in% c('Chlorpyrifos','Chlorpyrifos oxon'), expected_stock_conc := 100]
spidmap3[, stock_conc := ALIQUOT_CONCENTRATION] # stock_conc not in invitrodb yet for some of these
usecols <- c(intersect(names(spidmap), names(spidmap3)))
spidmap <- rbind(spidmap[,usecols,with=F], spidmap3[,c(usecols,'stock_conc'),with=F], fill = T)
head(spidmap)
setnames(spidmap, old = c('EPA_SAMPLE_ID'), new = c("spid"))
# for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
# update expected_stock_conc for individual compouunds where needed 
# for example, 
# spidmap[treatment %in% c("2,2',4,4',5,5'-Hexabromodiphenyl ether","Dibenz(a,h)anthracene"), expected_stock_conc := 10.0]
spidmap[, treatment := as.character(treatment)]
head(spidmap[, .(treatment, spid, expected_stock_conc)])

# check if every treatment name from the mea data maps to a unique sample in spidmap
setdiff(dat[wllt == 't',unique(treatment)], spidmap$treatment) # checkign for missed treatments "Acetomenophin"     "Chlorpyrifos Oxon" -> these will be renamed
spidmap[treatment %in% unique(dat$treatment), .N, by = .(treatment)][N > 1] # checking for treatments that match multiple spid
# if there is not a 1-to-1 correspondence, update treatment names in "supplemental_mea_treatment_name_map.csv"
# empty

# update treatment names with entries in "supplemental_mea_treatment_name_map.csv" corresponding to dataset
# (treatment -> "mea_treatment_name", "updated_treatment_name" column will match "PREFERRED_NAME"
dat <- update_treatment_names(dat, root_output_dir, dataset_title)

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
# 3361 F01" "3612 H02" "3612 H03"
dat[treatment %in% problem_comps, .N, by = .(treatment, rowi, coli, conc)][order(treatment, rowi, coli, conc)]
# Looking at calc files -> it looks like more decimal places were saved in the Calc file concs than in the mcl file conc's
# If round to 1 sig fig, all conc's are the same
dat[, .(length(unique(signif(conc,1)))), by = .(apid, rowi, coli)][V1 > 1]
# empty

# Just use conc from Calculations files
dat[treatment %in% problem_comps, conc := signif(unique(conc_srcf[grepl('[Cc]alc',srcf)]), 3), by = .(apid, rowi, coli)]
dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1]
# empty


# finally, run this:
source(file.path(scripts.dir, 'confirm_concs.R'))
dat <- confirm_concs(dat, spidmap, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30), update_concs_without_prompt = update_concs_without_prompt)
dat[is.na(conc), .N, by = .(treatment, spid)] # just media


# FINAL DATA CHECKS -------------------------------------------------------------
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Any other plots or things to check?
# for first culture, only include CO data, not the other 2 cultures, right?
dat[, date := sub('_.*$','',apid)]
dat[date == '20210127', .N, by = .(apid, treatment)]
# apid         treatment    N
# 1: 20210127_MW75-5620     Acetaminophen  609
# 2: 20210127_MW75-5620       Bisphenol A  609
# 3: 20210127_MW75-5620 Chlorpyrifos oxon 1827
# 4: 20210127_MW75-5620              DMSO  435
# 5: 20210127_MW75-5620             Media  696
dat[, .(length(unique(apid))), by = .(date)][order(date)]
#        date V1
# 1: 20210127  1
# 2: 20210217  6
# 3: 20210303  6
# 4: 20210317  6
# 5: 20210331  6

# what's with the NA rval's not getting set to 0?
#also, why so many?
dat[is.na(rval) & !(grepl('DIV[579]',acsn)), .N, by = .(acsn)]
# acsn   N
# 1:            CCTE_Shafer_MEA_dev_interburst_interval_mean_DIV12 122
# 2:                 CCTE_Shafer_MEA_dev_burst_duration_mean_DIV12 122
# 3:       CCTE_Shafer_MEA_dev_per_burst_interspike_interval_DIV12 122
# 4:                 CCTE_Shafer_MEA_dev_spike_duration_mean_DIV12 124
# 5:          CCTE_Shafer_MEA_dev_network_spike_duration_std_DIV12 130
# 6:   CCTE_Shafer_MEA_dev_inter_network_spike_interval_mean_DIV12 130
# 7: CCTE_Shafer_MEA_dev_per_network_spike_spike_number_mean_DIV12 124
# 8:                  CCTE_Shafer_MEA_dev_network_spike_peak_DIV12 124
# these are all bursting or network spike related, so can be NA
# NOTE THat wllq set to 0 where rval is NA in the save_lvl0_snapshot() (get_latest_dat.R)

# Making sure i have included everyting...
pfas.hits.tb <- as.data.table(read.xlsx('L:/Lab/NHEERL_MEA/Project PFAS 2019/MEA NFA/PFAS Hits chemicals.xlsx', sheet = 'Hits'))
which(is.na(pfas.hits.tb$Number)) # 43 101 102 103 104 105 106 107 108 109 110 111 112
setdiff(pfas.hits.tb[1:42, EPA_SAMPLE_ID], unique(dat$spid))
# empty
setdiff(unique(dat$spid), pfas.hits.tb[1:42, EPA_SAMPLE_ID])
# [1] "EX000530" "EX000531" "EX000596" "EX000595" "DMSO"     "Media"   
# where wllt != 't' or not a pfas compound, this okay

# OLD:
# # plates that don't have 6 control wells...
# dat[wllt == 'n', .(length(unique(paste0(rowi,coli)))), by = .(apid)][V1 != 6]
# # apid V1
# # 1: 20210127_MW75-5620  5
# # 2: 20210331_MW75-8102  5
# # 3: 20210331_MW75-8103  5
# # 4: 20210331_MW75-8104  5
# check.apids <- dat[wllt == 'n', .(length(unique(paste0(rowi,coli)))), by = .(apid)][V1 != 6, unique(apid)]
# dat[apid %in% check.apids, .N, by = .(wllt, treatment, apid)]
# # ah, yes, these are teh plates that contain Media instead of DMSO in one of the rows on each plate

# Some of the mean firing rate bvals look rather high...
par(mar = c(10,3,3,2))
plotdat <- dat[grepl('firing_rate_mean$',acsn)]
plotdat$apid <- factor(plotdat$apid, levels = sort(unique(plotdat$apid)), ordered = T)
stripchart(rval ~ apid, plotdat[wllt == 'n' & wllq == 1], vertical = T, las = 2, pch = 1)
stripchart(bval ~ apid, plotdat[, .(bval = median(rval[wllt == 'n' & wllq == 1])), by = .(apid)], vertical = T, las = 2, pch = 19, add = T)
title(main = 'Mean Firing Rate AUC bval by apid')

# check when I had included the plate from 20210127, before I realized I should set wllq to 0
# # first plate, nAE
# # let's check out the development
# plotdat <- dat[apid == '20210127_MW75-5620' & wllt == 'n' & grepl('active_electrodes_number_DIV',acsn)]
# plotdat[, DIV := as.numeric(sub('CCTE_Shafer_MEA_dev_active_electrodes_number_DIV','',acsn))]
# plotdat[, well := paste0(apid,rowi,coli)]
# plotdat <- plotdat[order(DIV)]
# plot(rval ~ DIV, plotdat, pch = 19)
# for (welli in unique(plotdat$well)) {
#   points(rval ~ DIV, plotdat[well == welli & wllq == 1], type = 'l')  
# }
# # woah, the highest has only 4 AE on DIV 12!! woah now!
# # I wonder if I messed this up... usuallly kathleen/seline catch this kind of thing
# 
# # what's normal for DIV 12 in control wells?
# dat[wllt == 'n' & wllq == 1 & grepl('active_electrodes_number_DIV12',acsn), .N, by = .(rval)][order(-N)]
# dat[wllt == 'n' & wllq == 1 & grepl('active_electrodes_number_DIV12',acsn), .N, by = .(nAE_above_10 = rval >= 10)]
# #    nAE_above_10   N
# # 1:        FALSE  13
# # 2:         TRUE 133
# # want to get a higher level view of this...
# stripchart(rval*60 ~ apid, dat[grepl('firing_rate_mean_DIV12',acsn) & wllt == 'n' & wllq == 1], vertical = T, las = 2, pch = 1)
# abline(h = 50)
# stripchart(rval ~ apid, dat[grepl('active_electrodes_number_DIV12',acsn) & wllt == 'n' & wllq == 1], vertical = T, las = 2, pch = 1, method = 'jitter')
# abline(h = 10, lty = 'dashed')
# 
# # What are these AB DMSO values that are so low?
# dat[wllt == 'n' & grepl('AB',acsn) & rval < 10000, .(treatment, apid, rowi, coli, rval)]
# # treatment               apid rowi coli rval
# # 1:      DMSO 20210127_MW75-5620    1    2 4808
# # 2:      DMSO 20210127_MW75-5620    2    2 3082
# # 3:      DMSO 20210127_MW75-5620    4    2 3044
# # 4:      DMSO 20210127_MW75-5620    5    2 2614
# # 5:      DMSO 20210127_MW75-5620    6    2 3450
# # yep, these all from the first plate again
# dat[wllt == 'n' & grepl('LDH',acsn) & rval < 0.2, .(treatment, apid, rowi, coli, rval)]
# #    treatment               apid rowi coli       rval
# # 1:      DMSO 20210127_MW75-5620    1    2 0.04423333
# # 2:      DMSO 20210127_MW75-5620    2    2 0.09703333
# # 3:      DMSO 20210127_MW75-5620    4    2 0.04553333
# # 4:      DMSO 20210127_MW75-5620    5    2 0.05153333
# # 5:      DMSO 20210127_MW75-5620    6    2 0.04673333

dat[, date := NULL]

# Additional checks
dat[spid == 'Media', .N, by = .(wllt, conc)]
# wllt conc    N
# 1:    b   NA 2436
# 2:    n   NA  348
dat[spid == "Media" & wllt != 'b', .N, by = .(apid, rowi, coli)]
# all col 2, from 4 plates
dat[wllt == 'n', .(length(unique(paste(rowi,coli)))), by = .(apid)]
# all exactly 6
dat[wllt == 't', .(length(unique(paste(rowi,coli)))), by = .(apid)]
# all 42, except for last culture 35 = 48 - 6 controls - 7 other Media

# first culture
dat[wllq == 0, .N, by = .(apid)]
# apid    N
# 1: 20210331_MW75-8101   87
# 2: 20210127_MW75-5620 4176
dat[apid == '20210127_MW75-5620', .N, by = .(wllq)] # all 0!

# save dat and graphs
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
