#' Read in wide-format AUC, DIV, and cytotoxicity data files and transform into 1 long-format data table
#'
#' @param project.output.dir folder where the AUC, DIV, and cytotoxicity data files are located
#' @param project_name name of the project
#' @param assay_component_map_filename file path to the csv file containing a map of the parameter value names as they appear in the AUC, DIV, and cytotoxicity data files to the corresponding assay component names registered in invitroDB.
#' Note that this function will append 'DIV#" to the parameter values in the input DIV table in order to match with the names in the assay_component_map_filename.
#' @param auc_data_file (optional) file path to the csv file containing the AUC parameter values (created by prepare_parameter_values.R). Defaults to the AUC data file in the project.out.dir/output folder.
#' @param div_data_file (optional) file path to the csv file containing the parameter values by DIV (created by parameter_values_to_AUC.R). Defaults to the _parameters_by_DIV data file in the project.out.dir/output folder.
#' @param cytotox_data_file (optional) file path to the csv file containing the AUC parameter values (created by run_cytotox_functions.R). Defaults to the _cytotox data file in the project.out.dir/output folder.
#' @param id.columns vector of column names in auc_data_file, div_data_file, and cytotox_data_file that should be interpreted as ID columns (i.e., all columns other than the parameters). 
#' This should include all columns that appear as an ID column in any of the 3 data files (any extra columns that do not apply to a particular data file will be ignored)
#'
#' @return
#' @export
#'
#' @examples
tcpl_MEA_dev_AUC <- function(project.output.dir, 
                             project_name, 
                             assay_component_map_filename,
                             auc_data_file = file.path(project.output.dir, "output", paste0(project_name, "_AUC.csv")), 
                             div_data_file = file.path(project.output.dir, "output", paste0(project_name, "_parameters_by_DIV.csv")), 
                             cytotox_data_file = file.path(project.output.dir, "output", paste0(project_name, "_cytotox.csv")),
                             id.columns = c("date","Plate.SN","DIV","well","trt","dose","units","file.name","wllq_by_well","wllq_notes_by_well","wllq_ref_by_well"))
{
  
  require(data.table)
  
  # get DIV data and melt
  div_data_wide <- as.data.table(read.csv(div_data_file))
  measure.columns <- setdiff(names(div_data_wide),id.columns)
  div_data_wide[, c(measure.columns) := lapply(.SD, as.numeric), .SDcols = measure.columns]
  div_data <- data.table::melt(div_data_wide, 
                               id.vars = intersect(names(div_data_wide),id.columns),
                               variable.name = 'src_acsn',
                               value.name = 'rval',
                               variable.factor = FALSE)
  
  # Exclude any dummy-data that was added in prepare_parameter_values() for the AUC calculation (usually at DIV 2)
  div_data <- div_data[!grepl('not recorded - dummy data at DIV', wllq_notes_by_well, fixed = T) & !grepl('not recorded - dummy data at DIV', file.name, fixed = T)]
  
  # Append the DIV to the source assay component source name (src_acsn) and remove the DIV column
  div_data[, `:=`(src_acsn = paste0(src_acsn,"_DIV",DIV),
                  srcf = basename(div_data_file))]
  div_data[, DIV := NULL]
  div_data[, file.name := NULL] # this column is DIV-specific and is not present in the AUC data
  
  # get auc_data data and melt
  auc_data_wide <- as.data.table(read.csv(auc_data_file))
  measure.columns <- setdiff(names(auc_data_wide),id.columns)
  auc_data_wide[, c(measure.columns) := lapply(.SD, as.numeric), .SDcols = measure.columns]
  auc_data <- data.table::melt(auc_data_wide, 
                               id.vars = intersect(names(auc_data_wide),id.columns),
                               variable.name = 'src_acsn',
                               value.name = 'rval',
                               variable.factor = FALSE)
  auc_data[, `:=`(srcf = basename(auc_data_file))]
  
  # rbind div_data and auc_data and clean some columns
  longdat <- rbind(div_data, auc_data)
  longdat[, `:=`(coli = as.numeric(sub("[[:alpha:]]","",well)), rowi = match(sub("[[:digit:]]","",well), LETTERS))]
  longdat[, well := NULL]
  setnames(longdat, old = "dose", new = "conc")
  rm(list = c("div_data","auc_data"))
  
  # Rename "trt" column as "treatment" (for historical continuity)
  # and convert to character because sometimes the treatment is read as an integer instead of a char
  setnames(longdat, old = 'trt', new = 'treatment', skip_absent = T)
  longdat[, treatment := as.character(treatment)] 
  
  # add cytotox data
  cytotox_data <- as.data.table(read.csv(cytotox_data_file))
  
  # Check alignment of names
  columns.to.keep <- c('apid', 'rowi', 'coli', 'treatment', 'conc', 'wllq_by_well', 'wllq_notes_by_well', 'wllt', 'rval', 'acsn', 'srcf', 'units')
  columns.in.mea.data.not.in.cytotox <- intersect(setdiff(names(longdat),names(cytotox_data)), columns.to.keep)
  if(length(columns.in.mea.data.not.in.cytotox) > 0) 
    warning('cytotox data does not have columns ',paste0(columns.in.mea.data.not.in.cytotox,collapse=","),
            '. Will fill with NA\n', sep ='')
  columns.in.cytotox.not.mea.data <- intersect(setdiff(names(cytotox_data),names(longdat)), columns.to.keep)
  if(length(columns.in.cytotox.not.mea.data) > 0) 
    warning('DIV & AUC data do not have columns ',paste0(columns.in.cytotox.not.mea.data,collapse=","),
            '. Will fill with NA\n', sep = '')
  longdat <- rbind(longdat, cytotox_data, fill = T)
  rm(list = c("cytotox_data"))
  
  # replace the src_acsn with the TCPL acsn
  assay_component_map <- as.data.table(read.csv(assay_component_map_filename, stringsAsFactors = FALSE))
  longdat <- merge(longdat, assay_component_map, by = c("src_acsn"), all.x = T)
  if (any(is.na(unique(longdat$acsn)))) {
    print(longdat[is.na(acsn), unique(src_acsn)])
    stop(paste0("The above src_acsn's are not found in ",assay_component_map_filename))
  }
  
  # Define apid
  longdat[, apid := paste(date, Plate.SN, sep = "_")]
  
  # get the desired columns, in the desired order
  longdat <- longdat[, .SD, .SDcols = intersect(columns.to.keep, names(longdat))]
  # longdat may or may not include "units"
  
  return(longdat)
}
