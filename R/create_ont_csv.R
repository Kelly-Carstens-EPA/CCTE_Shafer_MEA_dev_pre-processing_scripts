# create_ont_csvFile.R
# Diana Hall
# 3-14-2014
# purpose: to create a burst data automatically from package

create_ont_csv<-function(project.output.dir, remake_all = TRUE, get_h5Files_under_project.output.dir = TRUE, save.rdata = FALSE, param.file = NULL, AEfile = FALSE){  

  #load necessary packages
  library(sjemea)
  library(rhdf5)
  # library(lattice) # no longer needed
  # library(tcltk) # tk_choose.files() is intolerant of spaces in file names
  library(meadq)
  
  cat("\nStarting parameter calculations...\n")
  if (get_h5Files_under_project.output.dir) {
    h5Files <- list.files(path = file.path(project.output.dir, "h5Files"), pattern = "\\.h5$", full.names = TRUE, recursive = FALSE)
  }
  else {
    h5Files<-sort(choose.files(caption="Select .h5 Files") )
  }
  
  #create directories
  assign("h5.dir", dirname(h5Files[1]), envir = .GlobalEnv )
  assign("prepared.dir", paste(project.output.dir, "prepared_data", sep="/"), envir = .GlobalEnv )
  if(!dir.exists(prepared.dir)) dir.create( prepared.dir )

  # output file names
  csv.filename.AEfilt <- paste( prepared.dir, "/ont_data_summary_AEfilt",sep="")
  csv.filename.ABEfilt <- paste( prepared.dir, "/ont_data_summary_ABEfilt",sep=""  )
  
  if ( is.null( param.file ) ){
    data('chgv_parameters')
  } else {
    if ( grepl(x=basename( param.file) , pattern=".rda") ){
      load(param.file)
    } else{
      source( param.file, local=TRUE  ) 
    }
  }

  create_burst_ont_Data(h5Files=h5Files, csv.filename.AEfilt, csv.filename.ABEfilt, 
                        save.rdata=save.rdata, AEfile=AEfile, remake_all = remake_all)
  
  cat('Prepared data files are ready under',prepared.dir,'\n')
  
}
