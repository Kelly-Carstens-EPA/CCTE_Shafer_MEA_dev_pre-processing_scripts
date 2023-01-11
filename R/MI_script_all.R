# These functions allow you to select the h5 files
# then calls the functions to calculate the NMI

library(rhdf5) # this is the library from Bio Conductor
library(data.table)
require('pracma')
library('compiler')
require('gtools')

create_nmi_files <- function(h5Files, mi.dir, remake_all = TRUE) {
  
  # get a list of the date-plate combinations selected
  date_plates <- unique( lapply(h5Files, function(x) strsplit(basename(x), split = "_")[[1]][2:3]) )
  
  for (date_plate in date_plates) {
    
    # select the files corresponding to the current plate
    date <- date_plate[1]
    plate <- date_plate[2]
    plate_files <- h5Files[grepl(pattern = paste0(date,"_",plate,"_"), basename(h5Files))]
    file_split<-split(plate_files,1:length(plate_files)) # legacy thing, haven't had the gumption to change to vector yet
    
    file_name <- paste0(mi.dir,"/NMI_",date,"_",plate,".csv")
    if(!remake_all && file.exists(file_name)) {
      cat(basename(file_name),"already exists.\n")
      next
    }
    cat("Starting", date, plate, "at", as.character.Date(Sys.time()), "\n")
    

    #Begin MI analysis
    parsed_data = load.spikedata_final(333,file_split)
    MI_output<-list()
    MI_output<-nmi_wrapper(parsed_data)
    rm(parsed_data)

    # save the result as csv file
    fwrite(MI_output, file = file_name,col.names=TRUE, append =FALSE, row.names=FALSE, sep=",")
    rm(MI_output)
    
    cat("\tCompleted",basename(file_name),"at",as.character.Date(Sys.time()),"\n")
    
  }
  print(paste("MI files are ready in folder",mi.dir))
}


run_mi_functions <- function(basepath, get_h5Files_under_basepath = TRUE, remake_all = TRUE) {
  
  cat("\nStarting Normalized Mutual Information Calculations...\n")
  if (get_h5Files_under_basepath) {
    h5Files <- list.files(path = file.path(basepath, "h5Files"), pattern = "\\.h5$", full.names = TRUE, recursive = FALSE)
  }
  else {
    h5Files<-choose.files(caption="Select .h5 Files for NMI calculation")
    cat("Got",length(h5Files),"h5Files.\n")
  }
  h5Files <- sort(h5Files)
  
  # create All_MI directory
  if(basepath == "use_h5Files_dir") basepath <- dirname(dirname(h5Files[1]))
  mi.dir<-paste(basepath, "/All_MI",sep="")
  if(!(dir.exists(mi.dir))) dir.create(mi.dir)
  
  # calc the NMI for each h5File and save in csv!
  create_nmi_files(h5Files, mi.dir, remake_all = remake_all)
}

# Execute the function on it's own:
# run_mi_functions(basepath = "use_h5Files_dir", get_h5Files_under_basepath = FALSE)
