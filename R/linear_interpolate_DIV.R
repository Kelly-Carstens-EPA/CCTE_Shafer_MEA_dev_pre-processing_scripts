# functions to handle non-standard DIV recordings
# before calculating AUC values
# last updated July 2023 by Amy Carpenter

# linear interpolation function
# notes:
# - this function finds the approximate value between 2 other values
# - NA values are set to 0 in order to find the slope

linear_interpolate_DIV <- function(dati, DIV.to.interpolate) {
  
  # find the DIVs from original recordings that are above and below DIV.to.interpolate
  dati[, DIV := as.numeric(DIV)]
  all.divs.recorded <- dati[!grepl("interpolation",file.name),sort(unique(DIV))]
  if (DIV.to.interpolate < min(all.divs.recorded) | DIV.to.interpolate > max(all.divs.recorded))
    stop(paste0('DIV to interpolate (',DIV.to.interpolate,') is outside the range of DIVs recorded (',
                paste0(all.divs.recorded,collapse = ","),').\nCannot interpolate this DIV.'))
  DIV.lower <- max(all.divs.recorded[all.divs.recorded < DIV.to.interpolate])
  DIV.upper <- min(all.divs.recorded[all.divs.recorded > DIV.to.interpolate])
  cat("\tInterpolating values for DIV",DIV.to.interpolate, "from DIV",DIV.lower,"and",DIV.upper,"...\n")
  
  # For the purposes of interpolation, set NAs to 0 (debatable choice... but it's what we do right now)
  dati[is.na(value_by_DIV), value_by_DIV := 0]
  
  # Create table with values for lower and upper DIV
  id.columns <- setdiff(names(dati), c('DIV','file.name','value_by_DIV'))
  id.columns <- id.columns[!grepl('wllq',id.columns)]
  DIV.lower.dat <- dati[DIV == DIV.lower]
  DIV.upper.dat <- dati[DIV == DIV.upper]
  new.dat <- merge(DIV.lower.dat, DIV.upper.dat, by = id.columns,
                   suffixes = c('.lower','.upper'), all = T)
  stopifnot(nrow(new.dat[is.na(DIV.lower) | is.na(DIV.upper)]) == 0)
  
  # calculate the slope between DIV.lower to DIV.upper, add estimated value
  new.dat[, slope := (value_by_DIV.upper - value_by_DIV.lower)/(DIV.upper - DIV.lower)]
  new.dat[, value_by_DIV := value_by_DIV.lower + slope*(DIV.to.interpolate - DIV.lower)]
  new.dat[, DIV := DIV.to.interpolate]
  
  # Update wllq_by_well
  new.dat[, wllq_by_well := pmin(wllq_by_well.lower, wllq_by_well.upper, na.rm = T)]
  new.dat[, wllq_notes_by_well := paste0(paste0("Linear interpolation from DIV",DIV.lower," to DIV",DIV.upper, "; "),
                                         ifelse(!is.na(wllq_notes_by_well.lower),wllq_notes_by_well.lower,''),
                                         ifelse(!is.na(wllq_notes_by_well.lower) & !is.na(wllq_notes_by_well.upper),'; ',''),
                                         ifelse(!is.na(wllq_notes_by_well.upper),wllq_notes_by_well.upper,''))]
  new.dat[, wllq_ref_by_well := paste0(ifelse(!is.na(wllq_ref_by_well.lower),wllq_ref_by_well.lower,''),
                                       ifelse(!is.na(wllq_ref_by_well.lower) & !is.na(wllq_ref_by_well.upper),'; ',''),
                                       ifelse(!is.na(wllq_ref_by_well.upper),wllq_ref_by_well.upper,''))]
  
  new.dat[, .N, by = .(wllq_by_well, wllq_notes_by_well, wllq_ref_by_well)]
  
  # New file.name
  new.dat$file.name <- paste0("linear interpolation from DIV",DIV.lower," to DIV",DIV.upper)
  new.dat <- new.dat[, .SD, .SDcols = names(dati)]
  
  return(new.dat)
}
