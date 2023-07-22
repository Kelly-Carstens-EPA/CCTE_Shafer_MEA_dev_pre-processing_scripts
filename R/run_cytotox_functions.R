#' Main function to extract Alamar Blue and LDH data from Calculations or Summary .xlsx files
#' 
#'
#' @param project.output.dir folder where the output folder should be saved and where the files_log is located
#' @param project_name name of the project (used to name the output file)
#' @param get_files_from_log binary where to read files from files_log (TRUE) or from input vector of cytoFiles (FALSE, see next argument)
#' @param cytoFiles (optional) vector of full paths to .xlsx files containing blank-corrected LDH and AB data. File names should contain the word "Calculations" or "Summary"
#' @param wllq.tb.by.well.file (optional) full path to table containing well quality definitions and/or notes for specific plates, wells, and assays. 
#' Notes in the wllq.tb.by.well containing "AB" or "LDH" in the affected_assays are used here. If NULL, wllq_by_well will be set to 1 for all data rows
#'
#' @return
#' @export
#' 
#' @details This script will extract the blank-corrected fluorescence and optical density values for the LDH and Alamar Blue assays as part of the NFA.
#' Note that the script gets the blank-corrected values, NOT the pre-normalized percent of control values. Normalization is performed in level 3 of the ToxCast Pipeline.
#' Example input files:
#'   - "ON_20160720_MW1139-19_Summary.xlsx" - these contain data for 1 plate per sheet (LDH and Alamar Blue)
#'   - "20171011_ON G8_2 Calculations.xlsx" - these contain data for 3 plates per sheet (LDH and Alamar Blue)
#'
#' @examples
run_cytotox_functions <- function(project.output.dir, project_name, get_files_from_log = TRUE, cytoFiles = NULL, wllq.tb.by.well.file = NULL) {
  
  require(openxlsx)
  require(data.table)
  
  cat("\nStarting cytotoxicity data collection...\n")
  cat("Note: Any negative blank-corrected values will be set to 0.\n")
  
  # get the Calculations/Summary data files, either from files_log or provided cytoFiles
  if (get_files_from_log) {
    cytoFiles <- readFilesLog(project.output.dir, files_type = "Calculations")
    cytoFiles <- c(cytoFiles, readFilesLog(project.output.dir, files_type = "Summary"))
  }
  
  if (length(cytoFiles) == 0) {
    cat('No Summary or Calculations files found\n')
    return(0)
  }
  
  # run the functions for each file
  all_cyto_data <- data.table()
  for (i in 1:length(cytoFiles)) {
    
    if (grepl("(Calculations)|(Summary)",basename(cytoFiles[i]))) {
      AB_dat <- createCytoTable2(cytoFiles[i], cyto_type = "AB")
      LDH_dat <- createCytoTable2(cytoFiles[i], cyto_type = "LDH")
    }
    else {
      cat(paste("can't tell if",cytoFiles[i],"is 'Summary' file or 'Calculations' file\n"))
    }
    all_cyto_data <- rbind(all_cyto_data, AB_dat, LDH_dat)
    rm(list = c("AB_dat","LDH_dat"))
  }
  
  # Update the wllq
  if(!is.null(wllq.tb.by.well.file)) {
    all_cyto_data[, assay := src_acsn]
    all_cyto_data[, `:=`(date = as.character(date), rowi = as.numeric(rowi))]
    
    all_cyto_data <- add_wllq_by_well(all_cyto_data, wllq.tb.by.well.file, num_rows_per_plate = 6, num_columns_per_plate = 8)
    # remove columns that were created for or by the add_wllq_by_well() and are not needed going forward
    all_cyto_data[, intersect(c('assay'),names(all_cyto_data)) := NULL]
  } else{
    # Default to wllq == 1
    cat('No well quality table given - will default to wllq_by_well = 1 for all rows\n')
    all_cyto_data[, `:=`(wllq_by_well = 1,
                         wllq_notes_by_well = NA_character_,
                         wllq_ref_by_well = NA_character_)]
  }
  
  # flag NA rval's where wllq == 1
  na_rvals <- all_cyto_data[is.na(rval) & wllq_by_well == 1]
  if (nrow(na_rvals) > 0) {
    cat("The following rval's are NA, but wllq_by_well == 1. Wllq_by_well will be set to 0\n")
    print(na_rvals)
    all_cyto_data[is.na(rval) & wllq_by_well == 1, `:=`(wllq_by_well = 0,
                                                        wllq_notes_by_well = paste0("rval is NA; ",wllq_notes_by_well[!is.na(wllq_notes_by_well)]),
                                                        wllq_ref_by_well = paste0("rval is NA; ",wllq_ref_by_well[!is.na(wllq_ref_by_well)]))]
    warning('Some rvals were NA; Wllq_by_well set to 0 for affected wells/assays.\n')
  }
  
  # Print summary of wllq updates
  cat("\nWllq summary:\n")
  print(all_cyto_data[, .N, by = c("src_acsn","wllq_by_well","wllq_notes_by_well")][order(src_acsn, wllq_by_well, wllq_notes_by_well)])
  
  # Save the file
  if (!dir.exists(file.path(project.output.dir, "output"))) dir.create(file.path(project.output.dir, "output"))
  output_file <- file.path(project.output.dir, "output",paste0(project_name,'_cytotox.csv'))
  write.csv(all_cyto_data, file = output_file, row.names = FALSE)
  
  cat(file.path('output',basename(output_file)),"is ready\n")
  
}
