## MEA neuronal ontogeny AUC calculations
## Adapted from Chris Frank - February-March 2016
## Adapated from Mahmoud 2019 
## Last edited by Amy Carpenter July 2023

parameter_values_to_AUC <- function(project.output.dir, 
                                    project_name, 
                                    div_data_file_folders,
                                    id.columns = c("date","Plate.SN","DIV","well","trt","dose","units","file.name"),
                                    dose_column_name = NULL,
                                    use_DIVs = c(5,7,9,12),
                                    interpolate_diff_divs = TRUE,
                                    add_DIV2_values_of_0 = TRUE,
                                    wllq.tb.by.well.file = NULL,
                                    num_rows_per_plate = 6,
                                    num_columns_per_plate = 8) {
  

  
  
  # Fill all NA values with zero --------------------------------------------
  
  # Replace all NAs with zeros for AUC calculations 
  # **This may be undesirable for MEA parameters that are derived from other parameters.
  # But has been the standard practice for the past several years
  div_data[is.na(value_by_DIV), value_by_DIV := 0]
  
  
  # Calculate the AUC -------------------------------------------------------
  
  # Sort the rows by DIV
  div_data[, DIV := as.numeric(DIV)] # make sure DIV is numeric, not character
  div_data <- div_data[order(DIV)]
  
  # Calculate the trapezoidal area under the curve
  auc_data <- div_data[, .(value_AUC = trapz(x = DIV,
                                             y = value_by_DIV)),
                       by = c(id.columns.sans.div,'parameter')]
  auc_data[, value_AUC := round(value_AUC, digits = 6)] # for consistency with older functions
  auc_data[, parameter := paste0(parameter,'_auc')]
  
  # Make auc_data wide and save 
  auc_data_wide <- data.table::dcast(auc_data, ... ~ parameter, value.var = 'value_AUC')
  write.csv(auc_data_wide,
            file.path(project.output.dir, "output", paste0(project_name,"_AUC.csv")),
            row.names = FALSE)
  
  cat("\n", paste0(project_name,"_AUC.csv") ," is ready\n",sep = "")
  
}

