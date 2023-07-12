# create_ont_csvFile.R
# Diana Hall
# 3-14-2014
# purpose: to create a burst data automatically from package

#load necessary packages
library(sjemea)
library(rhdf5)
# library(lattice) # no longer needed
# library(tcltk) # tk_choose.files() is intolerant of spaces in file names
library(meadq)

create_ont_csv<-function(basepath, get_h5Files_under_basepath = TRUE, save.rdata = FALSE, param.file = NULL, AEfile = FALSE, remake_all = TRUE){  

  cat("\nStarting parameter calculations...\n")
  if (get_h5Files_under_basepath) {
    h5Files <- list.files(path = file.path(basepath, "h5Files"), pattern = "\\.h5$", full.names = TRUE, recursive = FALSE)
  }
  else {
    h5Files<-sort(choose.files(caption="Select .h5 Files") )
  }
  
  #create directories
  assign("h5.dir", dirname(h5Files[1]), envir = .GlobalEnv )
  if(basepath == "use_h5Files_dir") basepath <- dirname(h5.dir)
  assign("prepared.dir", paste(basepath, "prepared_data", sep="/"), envir = .GlobalEnv )
  if(!dir.exists(prepared.dir)) dir.create( prepared.dir )

  # output file names
  assign( "csv.filename.AEfilt",paste( prepared.dir, "/ont_data_summary_AEfilt",sep=""),
          envir = .GlobalEnv )
  assign( "csv.filename.ABEfilt",paste( prepared.dir, "/ont_data_summary_ABEfilt",sep=""  ),
          envir = .GlobalEnv )
  
 
  if ( is.null( param.file ) ){
    data('chgv_parameters')
  } else {
    if ( grepl(x=basename( param.file) , pattern=".rda") ){
      load(param.file)
    } else{
      source( param.file, local=TRUE  ) 
    }
  }

  create_burst_ont_Data(h5Files=h5Files, save.rdata=save.rdata, AEfile=AEfile, remake_all = remake_all)

} # end of create_ont_csv.R


# Execute the function
# create_ont_csv()
