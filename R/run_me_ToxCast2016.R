rm(list=ls()) # clear environment
graphics.off() # clear plot history
###################################################################################
# USER INPUT
###################################################################################
dataset_title <- "ToxCast2016" # the name for the current dataset, e.g. "name2020" (this should match the name of the folder under 'pre-process_mea_nfa_for_tcpl', e.g. 'Frank2017' or 'ToxCast2016')
pause_between_steps <- TRUE # probs want to be true when you first run
save_notes_graphs <- TRUE # Do this after have run thru once, to save a log of the steps. Set pause_between_steps to FALSE if saving notes and graphs for speed

default_ControlTreatmentName <- "DMSO" # usually DMSO. all compounds other than those listed below should have this vehicle control

scripts.dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl/nfa-spike-list-to-mc0-r-scripts/R"
root_output_dir <- "L:/Lab/NHEERL_MEA/Carpenter_Amy/pre-process_mea_nfa_for_tcpl" # where the dataset_title folder will be created

spidmap_file1 <- file.path(root_output_dir,"Sample IDs","EPA_12088_EPA-Shafer_96misc_75ul_20160826_key.xlsx")
spid_sheet1 <- 1
spidmap_file2 <- file.path(root_output_dir,"Sample IDs","EPA_11024_TShafer_384ph2_75ul_13May2015.xlsx")
spid_sheet2 <- 1
spidmap_file3 <- file.path(root_output_dir,"Sample IDs","Copy of NTP91_Compounds_4NHEERL_MEA_dev_cg.xlsx")
spid_sheet3 <- "NeuroTox 91 Cmpds"

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
# since I know the conc of DMSO, I am going to input that as well
dat[treatment == "DMSO", conc := 0.001] # 0.1% DMSO by volume


# Assign SPIDs ------------------------------------------------------------------
spidmap1 <- as.data.table(read.xlsx(spidmap_file1, sheet = spid_sheet1))
head(spidmap1)
unique(spidmap1$ALIQUOT_CONC_UNIT) # all mM
unique(spidmap1$ALIQUOT_SOLVENT) # DMSO
setnames(spidmap1, old = c("preferred_name","ALIQUOT_CONC", "EPA_SAMPLE_ID"), new = c("treatment","stock_conc","spid"))
# for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
spidmap1[, treatment := as.character(treatment)]
head(spidmap1[, .(treatment, spid, stock_conc)])

# add the second spidmap
spidmap2 <- as.data.table(read.xlsx(spidmap_file2, sheet = spid_sheet2))
head(spidmap2)
unique(spidmap2$ALIQUOT_CONC_UNIT) # all mM?
spidmap2[ALIQUOT_CONC_UNIT != "mM"] # C10-21 sulfonic acids phenyl esters. we aren't using this compound right now, so no worries 
unique(spidmap2$ALIQUOT_SOLVENT) # DMSO
setnames(spidmap2, old = c("dsstox_preferred_name","ALIQUOT_CONC", "EPA_SAMPLE_ID"), new = c("treatment","stock_conc","spid"))
# for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
spidmap2[, treatment := as.character(treatment)]
spidmap2 <- spidmap2[treatment %in% c("Clotrimazole","1H,1H,2H,2H-Perfluorooctyl iodide","Perfluoroundecanoic acid")]
head(spidmap2[, .(treatment, spid, stock_conc)])

# third spidmap -> Just for Valinomycin from NTP compounds
spidmap3 <- as.data.table(read.xlsx(spidmap_file3, sheet = spid_sheet3))
head(spidmap3)
setnames(spidmap3, old = c("Chemical.Name","Conc..(mM)", "SPID"), new = c("treatment","stock_conc","spid"))
# for example, setnames(spidmap, old = c("Aliquot_Vial_Barcode", "Concentration", "EPA_Sample_ID"), new = c("treatment","stock_conc","spid"))
spidmap3[, treatment := as.character(treatment)]
spidmap3 <- spidmap3[treatment == "Valinomycin"] # don't want other compounds, which could have diff spids
head(spidmap3[, .(treatment, spid, stock_conc)])
spidmap <- rbind(spidmap1[, .(treatment, spid, stock_conc)], spidmap2[, .(treatment, spid, stock_conc)], spidmap3[, .(treatment, spid, stock_conc)])
# spidmap[, length(unique(spid)), by = "treatment"][V1 != 1]
# spidmap[, length(unique(treatment)), by = "spid"][V1 != 1]

# finalize spidmap
spidmap[, expected_stock_conc := 20] # initialize expected_stock_conc. Usually this is 20mM. Change as needed.
# update expected_stock_conc for individual compouunds where needed 
# for example, 
# spidmap[treatment %in% c("2,2',4,4',5,5'-Hexabromodiphenyl ether","Dibenz(a,h)anthracene"), expected_stock_conc := 10.0]
spidmap[, stock_conc := as.numeric(stock_conc)]
head(spidmap[, .(treatment, spid, stock_conc, expected_stock_conc)])

# update names in "treatment" col to match "PREFERRED_NAME" in spidmap, set original treatments col to "mea_treatment_name"
dat <- update_treatment_names(dat, root_output_dir, dataset_title)

# assign spids
dat <- check_and_assign_spids(dat, spidmap)


# Confirm Conc's ----------------------------------------------------------------
# confirm that the conc's collected from master chem lists and Calc files match
# and that the correct concentration-corrections has been done for each compound

# check if there are multiple conc's assigned to the same well (usually occurs if there are differences between master chem file and calc file)
# Note: in TCPL mc1, the conc's are set to dat[ , conc := signif(conc, 3)]. So it's okay for us to round here.
dat[, .(num_unique_concs_in_well = length(unique(signif(conc,3)))), by = .(treatment, apid, rowi, coli)][num_unique_concs_in_well > 1]
# Empty data.table (0 rows and 5 cols): treatment,apid,rowi,coli,num_unique_concs_in_well
# if any, standardize those before continuing.

# con <- dbConnect(drv = RMySQL::MySQL(), user = Sys.getenv('INVITRODB_USER_RO'), pass = Sys.getenv('INVITRODB_PASS_RO'), dbname='invitrodb',host = Sys.getenv('INVITRODB_HOST_RO'))
# query_term <- paste0("SELECT * FROM sample WHERE spid IN('",paste(spidmap[!is.na(spid),unique(spid)],collapse="','",sep=""),"');")
# sample_info <- dbGetQuery(con, query_term)
# dbDisconnect(con)
# sample_info <- merge(sample_info, spidmap, by = c("spid"), suffixes = c(".db",".file"))
# setDT(sample_info)
# sample_info[signif(stkc,3) != signif(stock_conc,3)] # 10 compounds... huh
# spid  chid    stkc stkc_unit tested_conc_unit               treatment stock_conc expected_stock_conc
# 1: TP0001649B02 34695 10.0000        mM               uM                Mancozeb         20                  20
# 2: TP0001649B03 34187 10.0000        mM               uM               Tamoxifen         20                  20
# 3: TP0001649D07 22991  5.0000        mM               uM            Erythromycin         20                  20
# 4: TP0001649D10 20501 14.5000        mM               uM Methadone hydrochloride         20                  20
# 5: TP0001649E02 44175  5.2175     mg/ml             mg/l          Clove leaf oil         20                  20
# 6: TP0001649E06 20827 19.9000        mM               uM            Methoxychlor         20                  20
# 7: TP0001649E12 47525 10.0000        mM               uM      Pravastatin sodium         20                  20
# 8: TP0001649F10 32572 19.9000        mM               uM             Prallethrin         20                  20
# 9: TP0001649F11 20243 19.7000        mM               uM                  Captan         20                  20
# 10: TP0001649G07 47344 14.9000        mM               uM     Cariporide mesylate         20                  20
# confirmed with Ann Richard - the stock conc's in file are just the target conc's. 
# Invitrodb stkc is correct. Since the true stkc's are not present on the L drive, I can assume
# that whoever prepared these dilutions assumed that the stock conc was 20
# Additionally, for 1H,1H,2H,2H-Perfluorooctyl iodide with stkc = 5mM -> okay to apply conc-correction with expected conc at 20mM 
# see OneNote for notes

# finally, run this:
source(file.path(scripts.dir, 'confirm_concs.R'))
dat <- confirm_concs(dat, spidmap, expected_target_concs = c(0.03,0.1,0.3,1,3,10,20), update_concs_without_prompt = update_concs_without_prompt)


# FINAL DATA CHECKS
# this section is to confirm that the data has been processed correctly
source(file.path(scripts.dir, 'dataset_checks.R'))
dataset_checks(dat)

# Any other plots or things to check?
dat[grepl("MW1147-4",apid), .N, by = .(ifelse(grepl("DIV",acsn),"DIV","AUC"))]

# save the data and graphs
setkey(dat, NULL)
save(dat, file = file.path(root_output_dir, dataset_title, "output", paste0(dataset_title,"_longfile.RData")))
rm(dat)

if(save_notes_graphs) {
  sink() # close the txt log file
  graphics.off() # clear the plot history
}

cat("\nDone!\n")
