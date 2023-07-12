readFilesLog <- function(files.log.dir, files_type = "", files_log = "") {
  
  require(data.table)
  if (!files_type %in% c("","spike_list","MaestroExperimentLog","Calculations","Summary")) {
    stop(paste0("files_type must be one of 'spike_list','Calculations','Summary','MaestroExperimentLog' or empty"))
  }
  
  if (files_log == "") {
    files_log <- list.files(files.log.dir, pattern = '_files_log', full.names = T)
    stopifnot(length(files_log) == 1)
  }
  else if (dirname(files_log) == ".") {
    # If files_log does not contain path, add the file path
    files_log <- file.path(files.log.dir, files_log)
  }
  
  
  cat("\nReading from ",basename(files_log),"...",sep="")
  
  # function to read the files from the log file
  files_table <- read.table(files_log, sep = "\n", stringsAsFactors = F, col.names = c("col1"), blank.lines.skip = T)
  setDT(files_table)
  allfiles <- files_table[grepl("(\\.csv$)|(\\.xlsx$)|(\\.xls$)",col1),c(col1)]
  allfiles <- unique(allfiles) # in case any overmatching in grepl statement in writeFilesLog resulted in duplicated entries
  
  # isolate the file type that we want
  allfiles <- allfiles[grepl(pattern = files_type, basename(allfiles))]
  
  cat("\nGot",length(allfiles), files_type, "files.\n",sep=" ")
  return(allfiles)
}