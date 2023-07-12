rm(list=ls()) # clear environment
graphics.off() # clear plot history
###################################################################################
# USER INPUT
###################################################################################
project_name <- "SPS_PFAS2019" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- TRUE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

copy_maestro_exp_log_treatments <- FALSE # For cytotox data, keep treatment names from Calc/Summary files, or overwrite with maestro exp log treatment names?
default_ControlTreatmentName = "DMSO" # all compounds other than those listed below should have this vehicle control

spidmap_file <- ""
spid_sheet <- ""

scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R"
root.output.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the project_name folder will be created

update_concs_without_prompt <- TRUE
###################################################################################
# END USER INPUT
###################################################################################

library(data.table)
library(openxlsx)

# create a summary log file and store the 
if(save_notes_graphs) {
  sink(file = file.path(root.output.dir, project_name, paste0(project_name,"_run_log_",as.character.Date(Sys.Date()),".txt")), append = F)
  cat("Output from the script run_me_",project_name,".R\n",sep="")
  cat("Date Ran:",as.character.Date(Sys.Date()),"\n")
  cat(R.version.string,"\n")
  cat("USER INPUT settings:\n")
  print(sapply(ls(), get, envir = .GlobalEnv))
  graphics.off()
  pdf(file = file.path(root.output.dir, project_name, paste0(project_name,"_summary_plots_",as.character.Date(Sys.Date()),".pdf")), width = 10, height = 8)
}


# run the main steps
source(file.path(scripts.dir, 'source_steps.R'))


# run tcpl_MEA_dev_AUC
source(file.path(scripts.dir, 'tcpl_MEA_dev_AUC.R'))
dat <- tcpl_MEA_dev_AUC(basepath = file.path(root.output.dir,project_name), project_name)


# change untreated wells to Control Treatment ------------------------------------
dat[, conc_srcf := conc] # save the original conc's in a column

# wllt ==n wherever conc==0 currently
dat[wllt == 'n', unique(treatment)] # "DMSO"     "Media"    "3360 G12". Cool, most already labelled!
dat[treatment == '3360 G12', unique(.SD), .SDcols = c('apid','rowi','coli','srcf','conc')]
# apid rowi coli                                             srcf conc
# 1: 20201209_MW71-7211    1    6 20201209_NFA_PFAS_Group_2_SPS__Calculations.xlsx   30
# 2: 20201209_MW71-7213    2    6 20201209_NFA_PFAS_Group_2_SPS__Calculations.xlsx   30
# 3: 20201209_MW71-7214    3    6 20201209_NFA_PFAS_Group_2_SPS__Calculations.xlsx   30
# 4: 20201209_MW71-7211    1    6               SPS_PFAS2019_parameters_by_DIV.csv    0
# 5: 20201209_MW71-7213    2    6               SPS_PFAS2019_parameters_by_DIV.csv    0
# 6: 20201209_MW71-7214    3    6               SPS_PFAS2019_parameters_by_DIV.csv    0
# 7: 20201209_MW71-7211    1    6                             SPS_PFAS2019_AUC.csv    0
# 8: 20201209_MW71-7213    2    6                             SPS_PFAS2019_AUC.csv    0
# 9: 20201209_MW71-7214    3    6                             SPS_PFAS2019_AUC.csv    0
# From Seline 1/6 - this conc should be 30. Is not updated in maestroexplog
dat[treatment == '3360 G12' & conc == 0,`:=`(conc = 30, wllt = 't')]
dat[wllt == 'n', unique(treatment)] #  "DMSO"  "Media"
dat[treatment %in% c('DMSO','Media'), unique(wllt)] # "n" for all

# From Seline 1/6: they couldn't test 5 compounds because of low stock availability
# See Readme's from the 5 compounds. 1 from G1, 4 form G4
# currently these 5 labelled as "Media" in mea data, but as treated in calc files
dat[treatment %in% c('3360 C09','3612 G02','3612 H01','3612 H02','3612 H03'), `:=`(treatment ='Media', conc = 0)]

# Check for any disagreements betweecompare treatment names in AUC and Calc data
(qry_wells <- dat[, .(length(unique(treatment))), by = .(apid, rowi, coli)][V1 > 1, .(apid, rowi, coli)])
setkey(dat, apid, rowi, coli)
dat[J(qry_wells), .(unique(treatment)), by = .(apid, srcf, rowi, coli, conc)]
# empty

# Set the control well concentration. Adjust as needed
dat[wllt == "n", conc := 0.001] # not goign to bother finding this right now

# What is the concentration in all Bisphenol wells?
dat[treatment == 'Bisphenol', .N, by = .(conc, srcf)] # all say 30
dat[treatment == 'Bisphenol', unique(wllt)]
# [1] "t" # good, all are wllt t

# Setting Media to 'b' for blank
dat[treatment == 'Media', `:=`(wllt = 'b', conc = NA_real_)]


# assign sample ID's -------------------------------------------------------------
spidmap <- as.data.table(read.xlsx("../../Sample IDs/EPA_27864_EPA-Shafer_134_20191001_key.xlsx", sheet = 1))
spidmap2 <- as.data.table(read.xlsx("../../Sample IDs/EPA_29885_EPA-Shafer_36_20191112_key.xlsx", sheet = 1))
usecols <- intersect(names(spidmap), names(spidmap2))
spidmap <- rbind(spidmap[, ..usecols], spidmap2[, ..usecols])
spidmap[, short_rackplate := sub("SRACK0","",RACKPLATE_BARCODE)]
unique(spidmap$short_rackplate) # "3360" "3361"  "3612"
spidmap[, Compound.Name := paste(short_rackplate, WELL_POSITION, sep = " ")]
spidmap3 <- as.data.table(read.xlsx('../../Sample IDs/EPA_Shafer1_20210609.xlsx', sheet = 1)) # adding spidamp for bisphenol A
spidmap3[, Compound.Name := PREFERRED_NAME] # this sample of BPA will match by name
spidmap3[, EPA_SAMPLE_ID := BLINDED_SAMPLE_ID] # no epa sample id in this sheet, but I know it should be teh same as this col based on the url and email from Chris
spidmap3[, `:=`(CONCENTRATION = TARGET_CONCENTRATION,
                CONCENTRATION_UNIT = TARGET_CONCENTRATION_UNIT)]
usecols <- intersect(names(spidmap), names(spidmap3))
spidmap <- rbind(spidmap, spidmap3[, ..usecols], fill = T)
spidmap[, unique(CONCENTRATION_UNIT)] # mM for all
unique(spidmap$CONCENTRATION) # 30  5 20 10 - Since Kathleen adn Theresa have the spidmap files, I am assuming that they based their dilutions off of the concentrations listed here
setnames(spidmap, old = c('Compound.Name', 'EPA_SAMPLE_ID','CONCENTRATION'), new = c("treatment","spid","expected_stock_conc"))
spidmap[, treatment := as.character(treatment)]
head(spidmap[, .(treatment, spid)])

# check if every treatment name from the mea data maps to a unique sample in spidmap
setdiff(dat$treatment, spidmap$treatment) # "DMSO"      "Media"     "Bisphenol", that's fine, Bipshenol will be renamed below
spidmap[treatment %in% unique(dat$treatment), .N, by = .(treatment)][N > 1] # checking for treatments that match multiple spid - empty
# if there is not a 1-to-1 correspondence, update treatment names in "supplemental_mea_treatment_name_map.csv"

# update treatment names with entries in "supplemental_mea_treatment_name_map.csv" corresponding to dataset
# (treatment -> "mea_treatment_name", "updated_treatment_name" column will match "PREFERRED_NAME"
dat <- update_treatment_names(dat, root.output.dir, project_name)

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
problem_comps # empty
dat[treatment %in% problem_comps, unique(conc), by = .(spid, treatment, srcf)]

# confirmation that there is just 1 conc per spid, for SPS
dat[, .(length(unique(conc))), by = .(spid)][V1 != 1]
# empty, good

# finally, run this:
source(file.path(scripts.dir, 'confirm_concs.R'))
dat <- confirm_concs(dat, spidmap, expected_target_concs = c(30), update_concs_without_prompt = update_concs_without_prompt)


# FINAL DATA CHECKS -------------------------------------------------------------
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Any other plots or things to check?
dat[, well_id := paste(apid, rowi, coli, sep = "_")]
dat[, .(rep_count = length(unique(well_id))), by = .(treatment)][rep_count != 3]
# treatment rep_count
# 1: Bisphenol A        10
# 2:      DMSO        60
# 3:     Media        17
# cool, only control wells
dat[, well_id := NULL]

dat[is.na(rval) & !grepl('DIV',acsn)]
# empty

# Check if non-LDH points from the plate with LDH wllq == 0 look okay
ldh_bad_plate <- dat[grepl('LDH',acsn) & wllq == 0, unique(apid)]
ldh_bad_plate_trts <- dat[apid == ldh_bad_plate & grepl('firing_rate_mean_DIV12',acsn) & wllt == 't' & treatment != 'Bisphenol A', unique(treatment)]
ldh_bad_plate_and_replicates <- dat[treatment %in% ldh_bad_plate_trts, unique(apid)]
library(ggplot2)
for (acsni in c('CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12','CCTE_Shafer_MEA_dev_firing_rate_mean',
                'CCTE_Shafer_MEA_dev_active_electrodes_number_DIV12','CCTE_Shafer_MEA_dev_active_electrodes_number',
                'CCTE_Shafer_MEA_dev_burst_rate_DIV12','CCTE_Shafer_MEA_dev_burst_rate',
                'CCTE_Shafer_MEA_dev_network_spike_peak_DIV12','CCTE_Shafer_MEA_dev_network_spike_peak')) {
  plotdat <- dat[apid %in% ldh_bad_plate_and_replicates & acsn == acsni]
  plotdat[, med_rval := median(rval), by = .(spid, treatment)]
  setkey(plotdat, med_rval)
  plotdat$treatment <- factor(plotdat$treatment, levels = unique(plotdat$treatment), ordered = T)
  p <- ggplot(data = plotdat, mapping = aes(x = treatment, y = rval))+
          geom_point(mapping = aes(color = factor(apid)), shape = 1, stroke = 2, size = 2, alpha = 0.75)+
          scale_color_manual(breaks = c(setdiff(ldh_bad_plate_and_replicates, ldh_bad_plate), ldh_bad_plate),
                       values = c('black','blue','red'), name = 'apid')+
          ggtitle(label = paste0('Comparison of Replicates for Compounds Tested on Apid with Very High LDH Values (',ldh_bad_plate,')\n',acsni))+
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(p)
}

# Check if the plate that was temporarily misdosed looks okay compared to other replicates
check_apid <- dat[grepl('[Mm]isdose',wllq_notes), unique(apid)]
check_apid_trts <- dat[apid == check_apid &  wllt == 't' & rowi %in% c(1,3), unique(treatment)]
for (acsni in c('CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12','CCTE_Shafer_MEA_dev_firing_rate_mean',
                'CCTE_Shafer_MEA_dev_active_electrodes_number_DIV12','CCTE_Shafer_MEA_dev_active_electrodes_number',
                'CCTE_Shafer_MEA_dev_burst_rate_DIV12','CCTE_Shafer_MEA_dev_burst_rate',
                'CCTE_Shafer_MEA_dev_network_spike_number_DIV12','CCTE_Shafer_MEA_dev_network_spike_number')) {
  plotdat <- dat[treatment %in% check_apid_trts & acsn == acsni]
  plotdat[, med_rval := median(rval), by = .(spid, treatment)]
  setkey(plotdat, med_rval)
  plotdat$treatment <- factor(plotdat$treatment, levels = unique(plotdat$treatment), ordered = T)
  plotdat$apid <- factor(plotdat$apid, levels = c(setdiff(unique(plotdat$apid), check_apid), check_apid), ordered = T)
  p <- ggplot(data = plotdat, mapping = aes(x = treatment, y = rval))+
    geom_point(mapping = aes(color = factor(apid)), shape = 1, stroke = 2, size = 2, alpha = 0.75)+
    scale_color_manual(breaks = levels(plotdat),
                       values = c('black','blue','red'), name = 'apid')+
    ggtitle(label = paste0('Comparison of Replicates for Compounds Tested on Apid Momentarily Misdosed on DIV 5\n(',check_apid,', rows 1&3)\n',acsni))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(p)
}

# confirmation that all data values and data present are the same as in Jan 2021
cdat <- dat
# load('../../SPS_PFAS2019/output/output_2021-01-06/SPS_PFAS2019_longfile.RData')
load('L:/Lab/NHEERL_MEA/Project PFAS 2019/MEA NFA/SPS Hit Call Analysis/dat_SPS_PFAS2019_2021-01-06.RData')
setdiff(dat$acnm, cdat$acsn)
# character(0)
setdiff(cdat$acsn, dat$acnm) # ah, this file does not include the DIV 5-9 points
# the only thing that should be different is teh wllq and wllq_notes
testdat <- cdat[acsn %in% dat$acnm]
nrow(testdat) == nrow(dat) # TRUE
testdat_spid_counts <- testdat[, .N, by = .(spid)]
dat_spid_counts <- dat[, .N, by = .(spid)]
setdiff(testdat$spid, dat$spid) # Bisphenol A, of course!
dat[spid == 'Bisphenol', spid := 'Bisphenol A']
testdat_spid_counts <- testdat[, .N, by = .(spid)]
dat_spid_counts <- dat[, .N, by = .(spid)]
all.equal(testdat_spid_counts, dat_spid_counts, ignore.row.order = T)
# TRUE!

setnames(dat, old = 'acnm', new = 'acsn')
usecols <- c('spid','apid','rowi','coli','srcf','acsn','rval','wllt')
all.equal(testdat[, ..usecols], dat[, ..usecols], ignore.row.order = T, ignore.col.order = T, check.attributes = F)
# TRUE
# but, something in conc's and treatment's not matching...
check.dat <- merge(testdat, dat, by = usecols, all = T)
check.dat[is.na(treatment.x)] # empty
check.dat[is.na(treatment.y)] # empty -> cool, so it looks like nothing was dropped
check.dat[conc.x != conc.y] # empty...
check.dat[conc.x != conc.y | xor(is.na(conc.x), is.na(conc.y)), .N, by = .(conc.x, conc.y, spid)] # empty...
# conc.x conc.y spid    N
# 1:  0.001     NA DMSO 5220
# ah, the DMSO conc's!
check.dat[treatment.x != treatment.y, .N, by = .(treatment.x, treatment.y)]
#    treatment.x treatment.y   N
# 1: Bisphenol A   Bisphenol 870
# I know I updated this, so this is okay
rm(list = c('dat','testdat'))
dat <- cdat


# save dat and graphs
setkey(dat, NULL)
save(dat, file = file.path(root.output.dir, project_name, "output", paste0(project_name,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
