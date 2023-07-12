# script to gather the mea nfa spike list files, and other files

selectInputFiles <- function(start.dir, output.dir, project_name,append = F, files_type = "") {
  
  if (append) {
    file_names <- readFilesLog(output.dir, files_type)
  }
  else {
    file_names <- c()
  }
  
  # get starting folder, to initialize starting screen
  culture.dirs <- list.dirs(path = start.dir, recursive = F)
  previousfolder <- culture.dirs[1]
  
  repeat {
    
    add.files <- choose.files(default = previousfolder, caption = paste0("Select all ",files_type," analysis files for ",project_name))
    
    # loop breaks when user hits cancel
    if (length(add.files) == 0) {
      break
    }
    
    file_names <- c(file_names, add.files)
    previousfolder <- dirname(tail(add.files,n=1))
  }
  # just in case any files were selected twice
  file_names <- unique(file_names)
  
  writeFilesLog(file_names, output.dir, project_name, files_type = files_type)
}





