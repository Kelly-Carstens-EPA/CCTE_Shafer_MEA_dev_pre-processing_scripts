writeFilesLog <- function(file_names, files.log.dir, project_name) {
  
  # alphabetize/numerically sort file names so that files from same folders will be printed together
  file_names <- sort(file_names)
  
  # create log file name
  log_file <- file.path(files.log.dir, paste0(project_name,"_files_log.txt"))
  cat("Writing",length(file_names),"files to",basename(log_file),"...\n")
  
  # create the log file
  sink(file = log_file, append = F)
  cat(paste0(project_name," files used for MEA NFA pre-processing for TCPL\n"))
  cat("Created with the script writeFilesLog.R\n")
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