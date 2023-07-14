prepare_parameter_values <- function(project.output.dir, 
                                     project_name, 
                                     div_data_file_folders,
                                     id.columns = c("date","Plate.SN","DIV","well","trt","dose","units","file.name"),
                                     use_DIVs = c(5,7,9,12),
                                     interpolate_diff_divs = TRUE,
                                     add_DIV2_values_of_0 = TRUE,
                                     wllq.tb.by.well.file = NULL,
                                     num_rows_per_plate = 6,
                                     num_columns_per_plate = 8) {
  
  #pracma package has trapz function that computes AUC based on trapezoidal geometry 
  # without fitting a curve
  require(pracma) 
  require(data.table)
  options(digits = 6) # for consistency with older methods
  
  cat("\nStarting AUC preparations...\n")
  
  # Check that "DIV" is including in the id.columns
  # (DIV is used as the x-axis for the AUC calculation)
  stopifnot('DIV' %in% id.columns)
  id.columns.sans.div <- setdiff(id.columns, c('DIV','file.name'))
  
  # Set output folder
  if(!dir.exists(file.path(project.output.dir,"output"))) dir.create(file.path(project.output.dir,"output"))
  
  # Get the parameter data files
  div_data_files <- unlist(lapply(div_data_file_folders, list.files, pattern = '.csv$', full.names = TRUE, recursive = FALSE))
  
  
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
    div_data[, assay := 'nfa']
    div_data <- add_wllq_by_well(div_data, wllq.tb.by.well.file, num_rows_per_plate, num_columns_per_plate, all_DIVs = use_DIVs)
  }
  
  
  # Add points at DIV 2 -----------------------------------------------------
  
  # For consistency with old methods in Shafer lab, 
  # add a point of 0 at DIV == 2 for every parameter value
  if (add_DIV2_values_of_0) {
    div2_fill_0_table <- unique(div_data, by = id.columns.sans.div)
    div2_fill_0_table[, `:=`(DIV = 2,
                             value_by_DIV = 0)]
    div2_fill_0_table[, file.name := 'not recorded - dummy data at DIV 2']
    div2_fill_0_table[, wllq_notes_by_well := 'not recorded - dummy data at DIV 2']
    div_data <- rbind(div_data, div2_fill_0_table)
    rm(div2_fill_0_table)
    
    use_DIVs <- union(2, use_DIVs)
  }  
  
  
  # Interpolate values for non-standard DIV recordings ---------------------------------
  
  # NOTE - this step serves to interpolate non-standard DIV, NOTE for estimating missing DIV
  # For example, recordings usually made on DIV 5, 7, 9, and 12
  # but due to e.g., extreme weather on DIV 9, recordings were made on DIV 5, 7, 10, and 12 for a given plate
  # In this case, we want to interpolate what the value would have been on DIV 9
  # using the data from DIV 10 and 12
  
  if (interpolate_diff_divs) {
    
    # Cycle through the plates that have non-standard DIV one by one
    div_data[, date_plate := paste0(date,'_',Plate.SN)]
    update_plates <- div_data[!DIV %in% use_DIVs, unique(date_plate)]
    
    for (date_platei in update_plates) {
      cat(date_platei,"\n")
      
      # Identify DIV to add for affected plate
      add.DIVs <- div_data[date_plate == date_platei, setdiff(use_DIVs, DIV)]
      
      # For each DIV to add, run the linearInterpolate function
      for (add.DIV in add.DIVs) {
        add.dat <- linear_interpolate_DIV(div_data[date_plate == date_platei], DIV.to.interpolate = add.DIV)
        div_data <- rbind(div_data, add.dat)
      }
      
      # Remove the non-standard DIVs
      div_data <- div_data[DIV %in% use_DIVs]
    }
  } else {
    cat('The following non-standard DIVs are present in the div_data:\n')
    div_data[!DIV %in% use_DIVs, .N, by = .(DIV)]
    warning(paste0("\nData from some plates have non-standard DIVs."))
  }
  
  # Missing DIV check ---------------------------------------------------
  
  div_data[, DIV := as.numeric(DIV)]
  div_data[, divs_with_wllq_1 := paste0(sort(unique(DIV[wllq_by_well == 1])),collapse = ","), by = c(id.columns.sans.div)]
  use_DIVs_str <- paste0(sort(unique(use_DIVs)),collapse = ",")
  if(nrow(div_data[divs_with_wllq_1 != use_DIVs_str]) > 0) {
    cat('Some plates/wells do not have data (with wllq == 1) for all of the expected DIV:\n')
    div_data[divs_with_wllq_1 != use_DIVs_str, .(expected_DIVs = use_DIVs_str, .N), by = .(divs_with_wllq_1)]
    cat('Wllq will be set to 0 for wells/plates with missing DIVs\n')
    div_data[divs_with_wllq_1 != use_DIVs_str, `:=`(wllq_by_well = 0,
                                                    wllq_notes_by_well = 'some DIV are missing',
                                                    wllq_ref_by_well = 'wllq set to 0 in prepare_parameter_values()')]
    
  }
  
  
  # Old estimate missing DIV prep stuff:
  
  # check if ea plate has a recording for each of use_DIVs ------------------------
  # original method, assuming every plate had same missing DIVs
  # date_plates <- unique(div_data$date_plate)
  # for (date_platei in date_plates) {
  #   plate_divs <- div_data[date_plate == date_platei, sort(unique(DIV))]
  #   missing_divs <- setdiff(use_DIVs, plate_divs)
  #   if (length(missing_divs) > 0) {
  #     cat(paste0("There is no data for ",sub("_"," ",date_platei), " DIV ",paste0(missing_divs,collapse=","),"\n"))
  #     
  #     if (length(missing_divs) > 1) {
  #       warning(paste0("No data will be used from ",date_platei))
  #       div_data <- div_data[date_plate != date_platei]
  #     }
  #     else {
  #       cat("Values will be estimated from corresponding wells in other plates in same culture.\n")
  #       # Generate values for missing DIV by the median of other plates from this DIV
  #       div_data <- estimate_missing_DIV(dat = div_data, date_platei, missing_divs)
  #     }
  #   }
  # }
  # div_data <- div_data[, date_plate := NULL]
  cat("\nChecking that every plate has a recording on DIV",use_DIVs,"...\n")
  div_data[, well_id := paste(date, Plate.SN, well, sep = "_")]
  wells_missing_div <- div_data[, .(DIV_flag = ifelse(length(setdiff(use_DIVs, unique(DIV)))>0, 1, 0),
                                    missing_DIV = list(setdiff(use_DIVs, unique(DIV)))), by = c("date","Plate.SN","date_plate","well","well_id")][DIV_flag == 1]
  if (nrow(wells_missing_div) == 0) {
    cat("Every well has data from DIVs",use_DIVs,"\n")
  } else {
    check_plates <- wells_missing_div[, unique(date_plate)]
    for(date_platei in check_plates) {
      missing_divs <- wells_missing_div[date_plate == date_platei, unique(unlist(missing_DIV))]
      cat(date_platei,'is missing some recordings.\n')
      
      # loop through each missing_div on this plate (will check if multiple DIV estimated afterwards)
      # That way I still have estimated values, even if wllq set to 0
      for (add.DIV in missing_divs) {
        # Generate values for missing DIV by the median of other plates from this DIV
        div_data <- estimate_missing_DIV(dat = div_data, date_platei, add.DIV)
      }
    }
    
    # If more than 1 DIV value had to be estimated for a given well, set wllq==0
    cat("\nSetting wllq to 0 for plates with estimated values for more than 1 DIV...\n")
    check_multiple_missing <- function(wllq_notes_vector) {
      if (sum(grepl("estimated as median from corresponding wells",wllq_notes_vector)) > 1)
        paste0(wllq_notes_vector, "Multiple recordings missing; ")
      else
        # no changes needed for wllqnotes
        wllq_notes_vector
    }
    div_data[, wllq_notes_by_well := lapply(.SD, check_multiple_missing), .SDcols = "wllq_notes_by_well", by = "well_id"]
    div_data[grepl("Multiple recordings missing",wllq_notes_by_well), wllq_by_well := 0] # could I merge this step with above?
    cat(div_data[grepl("Multiple recordings missing",wllq_notes_by_well), length(unique(Plate.SN))],'plates affected.\n')
  }
  div_data <- div_data[, c("date_plate","well_id") := NULL]

  
  # Clean up -----------------------------------------------------------
  
  # rename a few items for consistency with previous code
  div_data[parameter == 'Mutual.Information', parameter := 'mi']
  
  # Confirm that have the same number of parameter values for every well and DIV
  # (i.e., confirm that have both MI and prepared data files for each plate)
  params.per.well.div <- div_data[, .(num_params = .N), by = c(id.columns)]
  max.num.params.per.well.div <- params.per.well.div[, max(num_params)]
  rows.with.missing.params <- params.per.well.div[num_params < max.num.params.per.well.div]
  if(nrow(rows.with.missing.params) > 0) {
    browser()
    cat('\nThe follow wells/DIVs appear to be missing some parameter values:\n\n')
    print(rows.with.missing.params)
    stop('Some parameter values may be missing for some wells/DIVs')
  }  
  
  # Convert the dose column to character
  # This is to ensure that no minor discrepancies in the decimal values 
  # from different input data files interfere with the groupings 
  # when group using "by" to get the AUC below
  # (Note that TCPL uses 3 sig figs for the concentration. 5 are saved here)
  div_data[, c(dose_column_name) := as.character(signif(get(dose_column_name), 5)) ]
  
  
  # Save file with all data by DIV ------------------------------------------
  
  # Save a snapshot of the combined prepared data, 
  # with the added/interpolated rows where DIV where missing
  div_data_wide <- data.table::dcast(div_data, ... ~ parameter, value.var = 'value_by_DIV')
  write.csv(div_data_wide, file = file.path(project.output.dir,"output",paste0(project_name,"_parameters_by_DIV.csv")), row.names = FALSE)
  rm(div_data_wide)
}