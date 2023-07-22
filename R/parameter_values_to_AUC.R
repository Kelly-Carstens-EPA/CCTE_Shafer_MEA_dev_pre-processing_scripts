#' Calculate Area-Under-the-Curve of parameter values vs DIV in each well
#' 
#' Original from Chris Frank - February-March 2016
#' Adapted by Mahmoud 2019 
#' Last edited by Amy Carpenter July 2023
#' Legacy script name: burst_parameter_to_AUC.R
#'
#' @param project.output.dir folder where the output folder should be saved
#' @param project_name name of the project (used to name the output file)
#' @param div_data_file full file path to single, wide .csv files containing the parameter data by DIV (created by prepare_parameter_values())
#' @param id.columns vector of column names in div_data_file that should be interpreted as ID columns (i.e., all columns other than the parameters)
#' @param expected_DIVs vector of the DIVs that should be present for all plates in div_data_file (defaults to 2, 5, 7, 9, 12). Used to check for missing DIV.
#' @param dose_column_name (optional) name of column in id.columns that contains the numeric dose. Will be rounded to 5 significant figures and converted to character to ensure consistency in representation of numeric value.
#'
#' @return
#' @export
#'
#' @examples
parameter_values_to_AUC <- function(project.output.dir, 
                                    project_name, 
                                    div_data_file,
                                    id.columns = c("date","Plate.SN","DIV","well","trt","dose","units","file.name","wllq_by_well","wllq_notes_by_well","wllq_ref_by_well"),
                                    expected_DIVs = c(2,5,7,9,12),
                                    dose_column_name = NULL) {
  
  # Read in the div_data_file & melt
  div_data_wide <- as.data.table(read.csv(div_data_file))
  # For data.table::melt, all measure.vars must be same type (numeric)
  measure.columns <- setdiff(names(div_data_wide),id.columns)
  div_data_wide[, c(measure.columns) := lapply(.SD, as.numeric), .SDcols = measure.columns]
  div_data <- data.table::melt(div_data_wide, id.vars = id.columns,
                               variable.name = 'parameter',
                               value.name = 'value_by_DIV',
                               variable.factor = FALSE)
  
  # Exclude columns that distinguish the DIV and wllq-related columns from id.columns
  # (these will be collapsed by well at the AUC-step)
  id.columns.sans.div <- setdiff(id.columns, c('DIV','file.name'))
  id.columns.sans.div <- id.columns.sans.div[!grepl('wllq',id.columns.sans.div)]
  
  
  # Convert the dose column to character
  # This is to ensure that no minor discrepancies in the decimal values 
  # from different input data files interfere with the groupings 
  # when group using "by" to get the AUC below
  # (Note that TCPL uses 3 sig figs for the concentration. 5 are saved here)
  if (!is.null(dose_column_name)) {
    div_data[, c(dose_column_name) := as.character(signif(get(dose_column_name), 5)) ]
  }
  
  # Fill all NA values with zero --------------------------------------------
  
  # Replace all NAs with zeros for AUC calculations 
  # **This may be undesirable for MEA parameters that are derived from other parameters.
  # But has been the standard practice for the past several years
  div_data[is.na(value_by_DIV), value_by_DIV := 0]
  
  
  # Calculate the AUC -------------------------------------------------------
  
  # Sort the rows by DIV
  div_data[, DIV := as.numeric(DIV)] # make sure DIV is numeric, not character
  div_data <- div_data[order(DIV)]
  
  # Exclude the wllq_note regarding intentionally imputed 'dummy data' before calculate AUC
  # (see prepare_parameter_values())
  div_data[, wllq_notes_by_well := sub('not recorded - dummy data at DIV 2[;]*','',wllq_notes_by_well)]
  div_data[wllq_notes_by_well == '', wllq_notes_by_well := NA_character_]
  
  # Calculate the trapezoidal area under the curve
  auc_data <- div_data[, .(
    # ** THE AUC STEP **
    value_AUC = trapz(x = DIV,
                      y = value_by_DIV),
    
    # Collapse wllq notes across all DIV
    wllq_by_well = min(wllq_by_well, na.rm = T),
    wllq_notes_by_well = paste0(unique(wllq_notes_by_well[!is.na(wllq_notes_by_well)]), collapse = "; "),
    wllq_ref_by_well = paste0(unique(wllq_ref_by_well[!is.na(wllq_ref_by_well)]), collapse = "; "),
    DIVs_with_wllq_1 = paste0(sort(unique(DIV[wllq_by_well == 1])),collapse = ",")
  ),
  by = c(id.columns.sans.div,'parameter')]
  
  # Clean up
  auc_data[, value_AUC := round(value_AUC, digits = 6)] # for consistency with historical methods
  auc_data[, parameter := paste0(parameter,'_auc')]
  
  
  # Missing DIV check ---------------------------------------------------
  
  # Check if any of the expected DIV are missing in wells with wllq_by_well == 1
  
  expected_DIVs_str <- paste0(sort(unique(expected_DIVs)),collapse = ",")
  
  if(nrow(auc_data[DIVs_with_wllq_1 != expected_DIVs_str & wllq_by_well == 1]) > 0) {
    
    cat('Some plates/wells appear to not have data (with wllq == 1) for all of the expected DIV (',expected_DIVs_str,'):\n')
    print(auc_data[DIVs_with_wllq_1 != expected_DIVs_str & wllq_by_well == 1, 
                   .(num_wells= length(unique(paste0(date,Plate.SN,well))),
                     expected_DIVs = expected_DIVs_str), by = .(DIVs_with_wllq_1)])
    
    cat('Wllq will be set to 0 for the AUC values for wells/plates with missing DIVs\n', sep = '')
  }
  
  # Update well quality where some DIV are missing
  auc_data[DIVs_with_wllq_1 != expected_DIVs_str, 
           `:=`(wllq_by_well = 0,
                wllq_notes_by_well = paste0('only these DIVs are present and have wllq=1: ',DIVs_with_wllq_1,'; ',
                                            wllq_notes_by_well),
                wllq_ref_by_well = paste0('wllq set to 0 bc some DIV are missing; ',
                                          wllq_ref_by_well))]
  auc_data[, DIVs_with_wllq_1 := NULL]
  
  
  # Save the output ---------------------------------------------------------
  
  # Make auc_data wide and save 
  auc_data_wide <- data.table::dcast(auc_data, ... ~ parameter, value.var = 'value_AUC')
  write.csv(auc_data_wide,
            file.path(project.output.dir, "output", paste0(project_name,"_AUC.csv")),
            row.names = FALSE)
  
  cat("\n", paste0('output/',project_name,"_AUC.csv") ," is ready\n",sep = "")
  
}

