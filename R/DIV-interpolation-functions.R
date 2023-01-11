# functions to handle non-standard DIV recordings
# before calculating AUC values

# basic function to move problem DIV to normal DIV
forceStandardDIV <- function(fdat, nonstandard_divs, target_divs) {
  for (i in 1:length(nonstandard_divs)) {
    fdat[DIV == nonstandard_divs[i], DIV := target_divs[i]]
  }
  return(fdat)
}

# linear interpolation function
# notes:
# - this function finds the approximate value between 2 other values
# - NA values are set to 0 in order to find the slope
linearInterpolateDIV <- function(dat2, new.DIV, remove.DIV = NULL) {
  
  # add DIV 2 = 0 data
  dat_div2 <- dat2[DIV == unique(dat2$DIV)[1]]
  dat_div2[, DIV := 2]
  id.cols <- intersect(names(dat2), c("date","Plate.SN","well","trt","dose","units","DIV","file.name","date_plate","wllq_by_well","wllq_notes_by_well"))
  endpoints_cols <- setdiff(names(dat_div2), id.cols)
  dat_div2[, c(endpoints_cols) := list(rep(0, .N))]
  dat_div2[, file.name := "DIV2_0_values"]
  dat_div2[, `:=`(wllq_by_well = 1, wllq_notes_by_well = "")]
  dat2 <- rbind(dat2, dat_div2)
  
  # find the DIVs from original recordings (or DIV 2) that are above and below new.DIV
  all.divs <- sort(c(new.DIV, dat2[!grepl("linear_interpolation",file.name), unique(DIV)]))
  lower.DIV <- all.divs[which(all.divs == new.DIV) - 1]
  upper.DIV <- all.divs[which(all.divs == new.DIV) + 1]
  cat("Interpolating values for DIV",new.DIV, "from DIV",lower.DIV,"and",upper.DIV,"\n")
  
  # create long table with upper and lower DIVs
  calc.dat2 <- dat2[DIV %in% c(lower.DIV, upper.DIV)]
  # calc.dat2[, `:=`(nAE = as.numeric(nAE), nABE = as.numeric(nABE), ns.n = as.numeric(ns.n))]
  calc.dat2.m <- melt(calc.dat2, id.vars = id.cols,
                      variable.name = "endpoint", value.name = "endpoint_value", na.rm=F, variable.factor = F)
  calc.dat2.m[, file.name := NULL]
  
  # for the purposes of interpolation, set NA to 0 (consider for all endpoints though...)
  calc.dat2.m[is.na(endpoint_value), endpoint_value := 0]
  
  # make wide data table
  calc.dat2.w <- dcast(calc.dat2.m, ... ~ DIV, value.var=list("endpoint_value","wllq_by_well","wllq_notes_by_well"), 
                       fun.aggregate = list(unique, min, function(x) paste(unique(wllq_notes_by_well),sep=""))) # possibly could just use "unique" for all, since there should only be 1 value
  calc.dat2.w$wllq_by_well <- do.call(pmin, calc.dat2.w[, .SD, .SDcols = grep("wllq_min",names(calc.dat2.w), val = T)]) # each col of the DT will be read as a separate argument to pmin
  # calc.dat2.w$wllq_notes_by_well <- do.call(paste0, args = list(unique(calc.dat2.w[, .SD, .SDcols = grep("wllq_notes_by_well_function",names(calc.dat2.w), val = T)]), collapse=""))
  setnames(calc.dat2.w, old = paste0("wllq_notes_by_well_function_",c(lower.DIV, upper.DIV)), new = paste0("wllq_notes_by_well_",c("lower.DIV", "upper.DIV")))
  calc.dat2.w[, wllq_notes_by_well := ifelse(wllq_notes_by_well_lower.DIV == wllq_notes_by_well_upper.DIV, 
                                     wllq_notes_by_well_lower.DIV, 
                                     paste0(wllq_notes_by_well_lower.DIV, wllq_notes_by_well_upper.DIV))] # trying to get teh unique wllq_notes_by_well, there is probably a better way to do this...
  setnames(calc.dat2.w, old = paste0("endpoint_value_unique_",c(lower.DIV,upper.DIV)), new = c("lower_col","upper_col"))
  calc.dat2.w <- calc.dat2.w[, .SD, .SDcols = names(calc.dat2.w)[!grepl(paste0("(",lower.DIV,")|(",upper.DIV,")|(lower\\.DIV)|(upper\\.DIV)"),names(calc.dat2.w))]] # removing unneeded columns
  
  # earlier method, before had wllq
  # calc.dat2.w <- dcast(calc.dat2.m, ... ~ DIV, value.var="endpoint_value")
  # setnames(calc.dat2.w, old = as.character(lower.DIV), new = "lower_col")
  # setnames(calc.dat2.w, old = as.character(upper.DIV), new = "upper_col")
  
  # calculate the slope between lower.DIV to upper.DIV, add estimated value
  calc.dat2.w[, slope := (upper_col - lower_col)/(upper.DIV - lower.DIV)]
  calc.dat2.w[, new_col := lower_col + slope*(new.DIV - lower.DIV)]
  
  # make endpoints columns again, filling with new_col values
  add_dat <- dcast(calc.dat2.w, date + Plate.SN + well + trt + dose + units + wllq_by_well + wllq_notes_by_well ~ endpoint, value.var = "new_col")
  
  # add new DIV column and file.name
  add_dat$DIV <- new.DIV
  add_dat$file.name <- paste0("linear_interpolation_from_DIV",lower.DIV,"_to_DIV",upper.DIV)
  add_dat[, wllq_notes_by_well := paste0(wllq_notes_by_well, "Linear interpolation from DIV",lower.DIV," to DIV",upper.DIV, "; ")]
  add_dat[, date_plate := paste(date, Plate.SN, sep = "_")]
  
  # add new data
  dat2 <- rbind(dat2, add_dat)
  
  # quick visual verification of the result
  plot(meanfiringrate ~ DIV, dat2[dose == 0], pch = "", ylab = "Mean Firing Rate (Hz) in control wells",
       main = paste0("Verification of Linear Interpolation of DIV",new.DIV,
       " from DIV",lower.DIV," and DIV",upper.DIV,"\nin control wells on ",unique(dat2$Plate.SN)),
       ylim = c(0, max(dat2[dose == 0, meanfiringrate])*1.1))
  for (welli in dat2[dose == 0, unique(well)]) {
    points(meanfiringrate ~ DIV, dat2[well == welli][order(DIV)], type = "l")
    points(meanfiringrate ~ DIV, dat2[well == welli][order(DIV)], pch = 19, col = "black")
    points(meanfiringrate ~ DIV, dat2[well == welli & grepl("linear_interpolation",file.name)][order(DIV)], pch = 19, col = "blue")
    points(meanfiringrate ~ DIV, dat2[well == welli & !DIV %in% c(2,5,7,9,12)][order(DIV)], pch = 19, col = "gray70")
    points(meanfiringrate ~ DIV, dat2[well == welli & wllq_by_well == 0][order(DIV)], pch = 19, col = rgb(1,0,0,alpha=0.7), cex=1.5)
    
    # # now I'm getting extra... trying to make an AUC plot with only the interpolated values
    # auc_pts <- dat2[well == welli & grepl("linear_interpolation",file.name)][order(DIV), .(DIV, meanfiringrate)]
    # polygon(x = c(2,))
  }
  legend(x = "topleft", legend = c("interpolated DIV","non-standard DIV to remove","wllq_by_well==0"), col = c("blue","gray70",rgb(1,0,0,alpha=0.7)), 
         pch = rep(19,length(col)), bg = "transparent")
  
  # remove added DIV 2 data
  dat2 <- dat2[DIV != 2]
  
  # if specified, remove the nonstandard DIV
  if (!is.null(remove.DIV)) {
    dat2 <- dat2[DIV != remove.DIV]
  }
  
  return(dat2)
}

#--------------- other things to try:
# - pick an arbitrary curve, fit it to data that is missing DIVs, then get approx DIV values

# - use existing data to determine the best curve to use for each endoint
## e.g., if hill usually fits best, then forcing a hill model fit will produce that shape
## hmm, could I even anchor/restrict the range of the x-value of the inflection points/critical points?

# - much later - see if better to calculate AUC from curve, instead of trapezoidal
