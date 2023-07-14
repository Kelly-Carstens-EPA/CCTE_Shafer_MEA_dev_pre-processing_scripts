# convert_spkListFiles_to_h5_MetaData.R
# Created by Diana Hall
# Adapted by Amy Carpenter

h5_conversion <- function(project.output.dir, 
                          files.log.output.dir = project.output.dir,
                          remake_all = FALSE, 
                          check_nwells_per_plate = 48,
                          recording_duration_sec = 900.00) {
  
  #load required package
  library(sjemea)
  library(rhdf5)
  library(meadq)
  
  #get files
  spkListFiles <- sort(readFilesLog(project.output.dir, files_type = "spike_list"))
  
  # define h5Files output directory
  h5.dir <- file.path(project.output.dir, "h5Files")
  if(!(dir.exists(h5.dir))) dir.create(h5.dir)
  
  # Determine which of the spike list files have already been converted to h5
  if (!remake_all) {
    existing_h5Files <- list.files(h5.dir, pattern = "\\.h5", full.names = FALSE)
    # Generate names of the h5 file names of the spkListFiles
    spkListFiles_h5names <- sub("_spike_list.csv$",".h5",basename(spkListFiles))
    # get the subset of spkListFiles that have not already been converted to h5Files
    spkListFiles <- spkListFiles[!(spkListFiles_h5names %in% existing_h5Files)]
    rm(list = c("existing_h5Files","spkListFiles_h5names"))
  }
  
  #get master chemical lists
  masterChemFiles <- readFilesLog(files.log.output.dir, files_type = "MaestroExperimentLog")
  
  if (length(spkListFiles)/4 > length(masterChemFiles)) {
    warning(paste("Only",length(masterChemFiles), 
                  "MaestroExperimentLog files found in files log for",
                  length(spkListFiles),
                  "spike list files.\n"))
  }
  
  # convert spike list files to .h5 files
  # (historically used to convert .mapTimestamps to h5 files)
  for (i in seq_along(spkListFiles)){
    
    # Find the date and plate number of current spike list file
    spikefilename <- basename(spkListFiles[i])
    date_plate <- paste(strsplit(spikefilename, split = "_")[[1]][2:3],collapse = "_")
    
    # Get the masterChemFile with the same date_plate (followed by _ or " ") as the spike list file
    masterChemFile <- grep(pattern = paste0(date_plate,"[_ ]"), masterChemFiles, value = T)
    # Check if any matching master chem file was found
    # If none matched, the "MW" from prefix of the plate number may be missing
    # try removing it
    if (length(masterChemFile) == 0) {
      masterChemFile <- grep(pattern = paste0(sub("MW","",date_plate),"[_ ]"), masterChemFiles, value = T)
    }
    
    # If still no matching master chem file found, or there are multiple matches, throw an error
    if (length(masterChemFile) != 1) {
      stop(paste("MaestroExperimentLog file match not found for",spikefilename,sep = " "))
    } 
    
    # get plate chemical info for each file in the list
    plate.chem.info<-chem.info.2(spkListFiles[i], masterChemFile)
    
    # check to see if the expected number of rows meta data there's data
    if(length(plate.chem.info$well) != check_nwells_per_plate) {
      warning(paste(check_nwells_per_plate,' rows not found in master chem file that have meta data that matches the spike list file:\n',
                    basename(masterChemFile),'\n',
                    basename(spkListFiles[i]),
                    '\nNo h5 will be created for this spike list file.\nLook at MaestroExperimentLog body and spike list file name and resolve any mismatches'))
      next
    }
    
    # make h5 files that contain chemical info 
    axion.spkList.to.h5(spkListFiles[i], plate.chem.info, h5.dir, 
                        recording_duration_sec = recording_duration_sec)
    
  }
  
  # indicate successful completion
  cat('h5files are ready under',h5.dir,'\n')
}
