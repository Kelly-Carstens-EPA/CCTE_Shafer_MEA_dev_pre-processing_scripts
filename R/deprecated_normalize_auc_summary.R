## function to normalize by plate for each ontogeny parameter
# Adapated from AUC_analysis_Revised_pipeline
# Last edited Feb 26, 2021 by Amy Carpenter
# This funciton is deprecated - we now use the ToxCast Pipeline R package to normalize the data

deprecated_normalize_auc_summary <- function(summary_table, direction=set_direction, plates_with_bad_controls = c()) {
  
  summary_table <- as.data.frame(summary_table)
  
  # standardize dose column, to prevent issues when select controls where dose == 0
  original_dose_class <- class(summary_table$dose)
  summary_table$dose <- as.numeric(summary_table$dose)
  summary_table$dose <- sprintf("%.5f", summary_table$dose)
  
  # Split data set by plate (date as well because plateIDs get reused)
  per_plate_split <- split(summary_table, interaction(summary_table[,"date"], summary_table[,"Plate.SN"], drop=TRUE))
  
  # Normalize each plate by percent of control median and combine into single data frame
  norm_plates <- list() # initialize list
  
  # Loop through each plate
  for (i in per_plate_split) {
    
    if (!(i[,"Plate.SN"][[1]] %in% plates_with_bad_controls)) { # Exclude any individual plates where untreated wells are off 
      
      # Loop through each ontogeny parameter
      for (j in grep('_auc',names(summary_table), val = T)) {
        
        cntrls <- (subset(i, dose == sprintf("%.5f", 0))[,j]) # Get vector of control values for that DIV on that plate
        i[,j] <- (i[,j] / median(na.omit(cntrls)))*100 # Divide all values on that plate by controls median
      }
      
    } else { # For plates singled-out above, use culture median instead of plate median
      
      for (j in grep('_auc',names(summary_table), val = T)) {
        
        cntrls <- (subset(summary_table, dose == sprintf("%.5f", 0) & date==i[,"date"][[1]])[,j]) # identify same-date control wells
        i[,j] <- (i[,j] / median(na.omit(cntrls)))*100  # Divide by same-culture control wells median     
      }
    }
    
    norm_plates[[length(norm_plates)+1]] <- i  # Add to growing list of normalized plates
  }
  
  rm(i,j) 
  norm_plates <- do.call(rbind, norm_plates) #Re-form one table of values
  
  # if set_direction is up, take positive difference as negative
  if (direction=="up") {
    norm_plates[,grep('_auc',names(norm_plates), val = T)] <- 100-(norm_plates[,grep('_auc',names(norm_plates), val = T)]-100)
  }
  
  # return 'dose' column to original class
  norm_plates$dose <- do.call(paste0('as.',original_dose_class), args = list(norm_plates$dose))
  
  #print(norm_plates) #checkpoint for normalization
  norm_plates
}
