#' Clean tables of parameter values by DIV and prepare for AUC calculation 
#' 
#' This function reads the wide parameter data files in `div_data_files`, merges them, defines the wllq_by_well, and performs other optional checks and clean up. 
#' The resulting single, wide data file containing all parameter values and DIV is saved as a csv in a folder named "output".
#'
#' @param project.output.dir folder where the output folder should be saved
#' @param project_name name of the project (used to name the output file)
#' @param div_data_files full file paths to wide .csv files containing the parameter data by DIV (with parameter values as column name(s))
#' @param id.columns vector of column names in div_data_files that should be interpreted as ID columns (i.e., all columns other than the parameters)
#' @param wllq.tb.by.well.file (optional) full path to table containing well quality definitions and/or notes for specific plates, wells, and DIVs. If NULL, wllq_by_well will be set to 1 for all data rows
#' @param num_rows_per_plate (optional) number of rows in an MEA plate. Used to expand data rows where well = 'all' in the wllq.tb.by.well. Defaults to 6 for a standard MEA 48-well plate.
#' @param num_columns_per_plate (optional) number of columns in an MEA plate. Used to expand data rows where well = 'all' in the wllq.tb.by.well. Defaults to 8 for a standard MEA 48-well plate.
#' @param expected_DIVs vector of the DIVs that should be present for all plates in div_data_files (defaults to 5, 7, 9, 12). Used to expands data rows where DIV = 'all' in wllq.tb.by.well and to check for non-standard DIV (see interpolate_stnd_DIVs)
#' @param add_DIV2_values_of_0 whether to add dummy-data values of 0 at DIV 2 (standard practice for the Shafer Lab). If TRUE, the function will add 2 to expected_DIVs
#' @param interpolate_stnd_DIVs whether to interpolate standard DIVs (defined by expected_DIVs) from non-standard DIVs. Only applies if non-standard DIVs are present in the div_data_files
#'
#' @return
#' @export
#'
#' @examples
prepare_parameter_values <- function(project.output.dir, 
                                     project_name, 
                                     div_data_files,
                                     id.columns = c("date","Plate.SN","DIV","well","trt","dose","units","file.name"),
                                     wllq.tb.by.well.file = NULL,
                                     num_rows_per_plate = 6,
                                     num_columns_per_plate = 8,
                                     expected_DIVs = c(5,7,9,12),
                                     add_DIV2_values_of_0 = TRUE,
                                     interpolate_stnd_DIVs = TRUE) {
  
  # pracma package has trapz function that computes AUC based on trapezoidal geometry 
  require(pracma) 
  require(data.table)
  options(digits = 6) # for consistency with older methods
  
  # Check that "DIV" is included in the id.columns
  # (DIV is used as the x-axis for the AUC calculation)
  stopifnot('DIV' %in% id.columns)
  
  # Get an additional vector that includes non-DIV specific columns
  id.columns.sans.div <- setdiff(id.columns, c('DIV','file.name'))
  
  # Set output folder
  if(!dir.exists(file.path(project.output.dir,"output"))) dir.create(file.path(project.output.dir,"output"))
  
  
  # Read files & melt -------------------------------------------------------
  
  div_data <- data.table()
  for (i in 1:length(div_data_files)) {
    data_wide_i <- as.data.table(read.csv(div_data_files[i], stringsAsFactors = F))
    # For data.table::melt, all measure.vars must be same type (numeric)
    measure.columns <- setdiff(names(data_wide_i),id.columns)
    data_wide_i[, c(measure.columns) := lapply(.SD, as.numeric), .SDcols = measure.columns]
    data_long_i <- data.table::melt(data_wide_i, id.vars = id.columns,
                                    variable.name = 'parameter',
                                    value.name = 'value_by_DIV',
                                    variable.factor = FALSE)
    div_data <- rbind(div_data, data_long_i)
  }
  div_data[, value_by_DIV := as.numeric(value_by_DIV)]
  
  # Update well quality by DIV ----------------------------------------------
  
  if(!is.null(wllq.tb.by.well.file)) {
    div_data[, assay := 'NFA']
    div_data <- add_wllq_by_well(div_data, wllq.tb.by.well.file, num_rows_per_plate, num_columns_per_plate, all_DIVs = expected_DIVs)
    # remove columns that were created for or by the add_wllq_by_well() and are not needed going forward
    div_data[, intersect(c('assay','rowi','coli'),names(div_data)) := NULL]
  } else{
    # Default to wllq == 1
    cat('No well quality table given - will default to wllq_by_well = 1 for all rows\n')
    div_data[, `:=`(wllq_by_well = 1,
                    wllq_notes_by_well = NA_character_,
                    wllq_ref_by_well = NA_character_)]
  }
  
  
  # Add dummy points at DIV 2 -----------------------------------------------------
  
  # For consistency with old methods in Shafer lab, 
  # add a point of 0 at DIV == 2 for every parameter value
  if (add_DIV2_values_of_0) {
    div2_fill_0_table <- unique(div_data, by = c(id.columns.sans.div,'parameter'))
    div2_fill_0_table[, `:=`(DIV = 2,
                             value_by_DIV = 0)]
    div2_fill_0_table[, file.name := 'not recorded - dummy data at DIV 2']
    div2_fill_0_table[, wllq_by_well := 1]
    div2_fill_0_table[, wllq_notes_by_well := 'not recorded - dummy data at DIV 2']
    div2_fill_0_table[, wllq_ref_by_well := NA_character_]
    div_data <- rbind(div_data, div2_fill_0_table)
    rm(div2_fill_0_table)
    expected_DIVs <- union(2, expected_DIVs)
  } 
  
  
  # Interpolate values from non-standard DIV to standard DIV ---------------------------------
  
  # NOTE - this step serves to interpolate values at standard DIV from
  # non-standard DIV recordings (NOT to estimate missing DIV if there are < 4 recordings).
  # For example, recordings usually made on DIV 5, 7, 9, and 12
  # but due to e.g., extreme weather on DIV 9, recordings were made on DIV 5, 7, 10, and 12 for a given plate
  # In this case, we want to interpolate what the value would have been on DIV 9
  # using the data from DIVs 7 and 10
  
  # Cycle through the plates that have non-standard DIV one by one
  div_data[, date_plate := paste0(date,'_',Plate.SN)]
  plates_with_non_standard_DIV <- div_data[!DIV %in% expected_DIVs, unique(date_plate)]
  
  if (length(plates_with_non_standard_DIV) > 0) {
    
    if (interpolate_stnd_DIVs) {
      
      for (date_platei in plates_with_non_standard_DIV) {
        cat(date_platei,"\n")
        
        # Identify DIV to add for current plate
        add.DIVs <- div_data[date_plate == date_platei, setdiff(expected_DIVs, DIV)]
        
        # For each DIV to add, run the linearInterpolate function
        for (add.DIV in add.DIVs) {
          add.dat <- linear_interpolate_DIV(dati = div_data[date_plate == date_platei], DIV.to.interpolate = add.DIV)
          div_data <- rbind(div_data, add.dat)
        }
        
        # Remove the non-standard DIVs
        div_data <- div_data[DIV %in% expected_DIVs]
      }
      
    } else {
      cat('The following non-standard DIVs are present in the div_data:\n')
      print(div_data[!DIV %in% expected_DIVs, .N, by = .(DIV)])
      warning(paste0("\nData from some plates have non-standard DIVs."))
    }
  }
  
  # remove created column
  div_data[, date_plate := NULL]
  
  
  # Clean up -----------------------------------------------------------
  
  # rename a few items for consistency with previous code
  div_data[parameter == 'Mutual.Information', parameter := 'mi']
  
  # Confirm that have the same number of parameter values for every well and DIV
  # (i.e., confirm that have both MI and prepared data files for each plate)
  params.per.well.div <- div_data[, .(num_params = .N), by = c(id.columns)]
  max.num.params.per.well.div <- params.per.well.div[, max(num_params)]
  rows.with.missing.params <- params.per.well.div[num_params < max.num.params.per.well.div]
  if(nrow(rows.with.missing.params) > 0) {
    cat('\nThe following wells/DIVs appear to be missing some parameter values:\n\n')
    print(rows.with.missing.params)
    stop('Some parameter values may be missing for some wells/DIVs')
  }  
  
  
  # Save file with all data by DIV ------------------------------------------
  
  # Save a snapshot of the combined prepared data, 
  # with the added/interpolated rows where DIV where missing
  div_data_wide <- data.table::dcast(div_data, ... ~ parameter, value.var = 'value_by_DIV')
  div_data_wide[, DIV := as.numeric(DIV)]
  setkeyv(div_data_wide, id.columns) # order by ID columns
  write.csv(div_data_wide, file = file.path(project.output.dir,"output",paste0(project_name,"_parameters_by_DIV.csv")), row.names = FALSE)
  cat(paste0('\noutput/',project_name,"_parameters_by_DIV.csv"),'is ready\n')
}