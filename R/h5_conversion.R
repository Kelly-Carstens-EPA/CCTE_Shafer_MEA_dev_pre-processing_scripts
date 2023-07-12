# convert_spkListFiles_to_h5_MetaData.R
# Diana Hall
#load required package
library(sjemea)
library(rhdf5)
library(meadq)

#get files
analysis <- readFilesLog(main.output.dir, files_type = "spike_list")
spkListFiles <- sort(analysis)

###################################################################################
# USER INPUT
###################################################################################
# set location where you want the "h5Files" folder to be created
basepath <- main.output.dir
# # or, if you want the h5Files folder to be created next the spike list files folder, use the following line
# basepath = dirname(dirname(spkListFiles[1]))
# If an h5Files folder already exists in  this location, the h5Files created will overwrite what is there
# remake_all <- TRUE # re-create existing h5files from spike list files?
###################################################################################
# END USER INPUT
###################################################################################

# create h5Files directory
suppressWarnings( dir.create(paste(basepath,'/h5Files',sep='') ) )
h5.dir<-paste(basepath, "/h5Files",sep="")

# Determine which of the spike list files have already been run
if (!remake_all) {
  existing_h5Files <- list.files(h5.dir, pattern = "\\.h5")
  existing_h5Files <- sub("\\.h5","",existing_h5Files)
  # Find the h5 file names of these spkListFiles,
  # using same h5file naming structure as here and in spike_list_functions.R
  spkListFiles_h5names <- sapply(spkListFiles, function(filei) 
    h5file_namei <- sub("_spike_list.*$","",basename(filei))
    )
  # get only the subset of spkListFiles that are not already in h5.dir
  use_files <- spkListFiles_h5names[!(spkListFiles_h5names %in% existing_h5Files)]
  spkListFiles <- names(use_files) # the "names" of use_files contain the unchanged file.path's
  rm(list = c("existing_h5Files","spkListFiles_h5names","use_files"))
}

#get master chemical lists
masterChemFiles <- readFilesLog(main.output.dir, files_type = "MaestroExperimentLog")

if (length(spkListFiles)/4 > length(masterChemFiles)) {
  cat("Only",length(masterChemFiles), "master experiment log file selected for",length(spkListFiles),"spike list files.\n")
  repeat {
    resp <- readline(prompt = "Continue anyway? (y/n): ")
    if (resp %in% c("y","n")) break
  }
  if (resp == "n") stop()
}

L=length(spkListFiles)

#convert .mapTimestamps to .h5 files
for (i in seq_along(spkListFiles)){
  
  # Find the plate number of current spike list file
  spikefilename <- basename(spkListFiles[i])
  date_plate <- paste(strsplit(spikefilename, split = "_")[[1]][2:3],collapse = "_")
  
  # Get the masterChemFile with the date_plate
  masterChemFile <- grep(pattern = paste0(date_plate,"[_ ]"), masterChemFiles, value = T)
  if (length(masterChemFile) == 0) {
    # sometimes the master chem file does not include "MW" in the file name
    masterChemFile <- grep(pattern = paste0("_",sub("MW","",date_plate),"[_ ]"), masterChemFiles, value = T)
  }
  # sometimes, maestroExperimentLog does not contain the plate.SN in the file name. Match by plate folder and date in file name instead
  if (length(masterChemFile) == 0) {
    masterChemFile <- Filter(function(mcf) grepl(sub("^.*_MW( )*","",date_plate),mcf) && grepl(sub("_.*$","",date_plate),basename(mcf)), masterChemFiles)
  }
  
  # If still no match found, or there are multiple matches, throw an error
  if (length(masterChemFile) != 1) {
    stop(paste("master chem file match not found for",spikefilename,sep = " "))
  } 

  #load data from all files in folder (listed above)
  title<-strsplit(basename(spkListFiles[i]), "[.]")[[1]][1]
  #get plate chemical info for each file in the list
  plate.chem.info<-chem.info.2(spkListFiles[i], masterChemFile)
  #check to see if there's data
  matching.log.file<-!(length(plate.chem.info)==0)
  
  if(length(plate.chem.info)==0) {
    cat('Meta data in',masterChemFile,'does not match',basename(spkListFiles[i]),'.\n')
    cat('Spike list file will not be converted to h5. Look at Maestro Experiment Log and spike list files and resolved any namining mismatch.\n')
    warning(paste0('Meta data in',masterChemFile,'does not match',basename(spkListFiles[i]),'. No h5 file created for this spike list file.'))
    next
  }
  
  #make h5 files that contain chemical info 
  axion.spkList.to.h5(title, spkListFiles[i], plate.chem.info)
  
}
