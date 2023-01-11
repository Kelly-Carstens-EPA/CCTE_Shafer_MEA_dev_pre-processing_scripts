# script to gather the mea nfa spike list files, and other files

selectInputFiles <- function(start.dir, output.dir, dataset_title,append = F, files_type = "") {
  
  if (append) {
    file_names <- readLogFile(output.dir, files_type)
  }
  else {
    file_names <- c()
  }
  
  # get starting folder, to initialize starting screen
  culture.dirs <- list.dirs(path = start.dir, recursive = F)
  previousfolder <- culture.dirs[1]
  
  repeat {
    
    add.files <- choose.files(default = previousfolder, caption = paste0("Select all ",files_type," analysis files for ",dataset_title))
    
    # loop breaks when user hits cancel
    if (length(add.files) == 0) {
      break
    }
    
    file_names <- c(file_names, add.files)
    previousfolder <- dirname(tail(add.files,n=1))
  }
  # just in case any files were selected twice
  file_names <- unique(file_names)
  
  writeLogFile(file_names, output.dir, dataset_title, files_type = files_type)
}


writeLogFile <- function(file_names, output.dir, dataset_title, files_type) {
  
  # alphabetize/numerically sort file names so that files from same folders will be printed together
  file_names <- sort(file_names)
  
  # create log file name
  log_file <- file.path(output.dir, paste0(dataset_title,"_files_log_",as.character.Date(Sys.Date()),".txt"))
  cat("Writing",length(file_names),"files to",basename(log_file),"...\n")
  
  # create the log file
  sink(file = log_file, append = F)
  cat(paste0(dataset_title," files used for MEA NFA pre-processing for TCPL\n"))
  cat("Created with the script gather_files-functions.R\n")
  cat("Date created: ")
  cat(as.character.Date(Sys.time()))
  cat("\nEvery line ending in '.csv' or '.xlsx' or '.xls' will be read as an input file")
  
  all_dirs <- unique(dirname(dirname(file_names)))
  common_folders <- Reduce(intersect, strsplit(all_dirs, split = "/|\\\\"))
  common_dir <- Reduce(file.path, common_folders)
  common_dir <- sub("\\(","\\\\(",common_dir)
  common_dir <- sub("\\)","\\\\)",common_dir)
  # or, just pass start.dir to writeFileLog, then sub start.dir with "" in for loop belows
  
  cat("\n\nMain directory:",common_dir,"\n")
  cat(paste0("Collected ",length(file_names)," files from ",length(all_dirs)," sub-directories.\n"))
  
  for (diri in all_dirs) {
    cat("\n",sub(paste0(common_dir,"/"),"",diri),"\n",sep = "")
    diri_files <- file_names[dirname(dirname(file_names)) == diri]
    cat(diri_files, sep = "\n")
  }
  
  closeAllConnections()
  
  print(paste0(log_file," is ready."))
  
}


readLogFile <- function(output.dir, files_type = "", files_log = "") {
  
  require(data.table)
  if (!files_type %in% c("","spike_list","MaestroExperimentLog","Calculations","Summary")) {
    stop(paste0("files_type must be one of 'spike_list','Calculations','Summary','MaestroExperimentLog' or empty"))
  }
  
  if (files_log == "") {
    # read the data from the most recent files_log in output.dir
    files_logs <- list.files(path = output.dir, pattern = "_files_log_.*\\.txt", recursive = F, full.names = T)
    files_log <- files_logs[order(basename(files_logs), decreasing = T)[1]] # get the most recent log file
  }
  else if (dirname(files_log) == ".") {
    # make sure files_log is a full.name
    files_log <- file.path(output.dir, files_log)
  }
  
  # send the name of files_log to the environ where read_files was called so that the chosen files_log can be documented
  assign("files_log",files_log, envir = parent.frame())
  
  cat("\nReading from ",basename(files_log),"...",sep="")
  
  # function to read the files from the log file
  files_table <- read.table(files_log, sep = "\n", stringsAsFactors = F, col.names = c("col1"), blank.lines.skip = T)
  setDT(files_table)
  allfiles <- files_table[grepl("(\\.csv$)|(\\.xlsx$)|(\\.xls$)",col1),c(col1)]
  allfiles <- unique(allfiles) # in case any overmatching in grepl statement in writeLogFile resulted in duplicated entries
  
  # isolate the file type that we want
  allfiles <- allfiles[grepl(pattern = files_type, basename(allfiles))]
  
  cat("\nGot",length(allfiles), files_type, "files.\n",sep=" ")
  return(allfiles)
}
