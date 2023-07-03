dataset_checks <- function(dat) {
  
  # for transitioning to assigning spids later
  remove_spid_col <- FALSE
  if (!is.element("spid",names(dat))) {
    remove_spid_col <- TRUE
    dat[, spid := treatment]
  }
  
  # this section is to confirm that the data has been processed correctly
  cat("\nFinal Checks\n")
  cat("Number of cultures dates:",dat[, length(unique(sub("_.+$","",apid)))])
  cat("\nRange of culture dates:", dat[, range(sub("_.+$","",apid))] )
  cat("\nNumber of plates tested:",dat[, length(unique(apid))])
  cat("\nNumber of compounds tested:",dat[wllt == "t", length(unique(spid))])
  cat("\nNA or infinite rvals for AUC or cytotoxicity endpoints:\n")
  print(dat[(is.na(rval) | is.infinite(rval)) & !grepl('DIV',acsn), .N, by = "wllq"])
  cat("\nWllq counts for all data points:\n")
  print(dat[, .N, by = "wllq"]) # note if wllq is NA anywhere
  cat("Number of unique acsn's present:",length(unique(dat$acsn))," (expect 2 cyto + 17 AUC + 4*17 DIV = 87)\n")
  check.points <- dcast(dat[, .N, by = c("acsn","apid")], apid ~ acsn, value.var = "N", fill = 0)
  setnames(check.points, old = names(check.points), new = sub("CCTE_Shafer_MEA_dev_","",names(check.points)))
  cat(paste0("The following plates don't have the expected number of points (48):\n"))
  standard_cols <- setdiff(names(check.points), c("apid"))
  pts_flag <- FALSE
  for (apidi in unique(check.points$apid)) {
    if (check.points[apid == apidi, any(.SD != 48), .SDcols = c(standard_cols)]) {
      pts_flag <- TRUE
      MEA_pts <- check.points[apid == apidi, .(sort(unique(.SD))), .SDcols = setdiff(standard_cols, c("AB","LDH"))]
      check.points[apid == apidi, MEA_pts := paste0(sort(unique(unlist(MEA_pts))),collapse=",")]
      print(check.points[apid == apidi, .SD, .SDcols = intersect(c("apid","AB","LDH","MEA_pts"),names(check.points))])
    }
  }
  if(!pts_flag) {
    cat("(all plates have the expected number of points for each assay component)\n")
  }
  
  # check number of controls per plate
  cat("\nApid/acsn pairs without 6 control wells:\n")
  print(dat[wllt == "n", .N, by = c("acsn","apid")][N != 6])
  
  # range by acsn
  cat("\nRange of rval's by acsn:\n")
  print(dat[wllq == 1, .(min = format(min(rval,na.rm=T),digits=2,scientific = F), 
                   median = format(median(rval,na.rm=T),digits=2,scientific = F),
                   max = format(max(rval,na.rm=T),digits=2,scientific = F),
                   num_NA = sum(is.na(rval))), by = "acsn"][order(acsn)])
  
  # PLOTS to visually confirm results
  
  # view all by plate
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "t" & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12"], vertical = T, pch = 1, method = "jitter", las = 2, cex.axis = 0.75,
             col = "cornflowerblue", main = paste0(project_name," NFA Mean Firing Rate DIV12 by Plate"))
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "n" & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean_DIV12"], vertical = T, pch = 19, method = "jitter", las = 2, cex.axis = 0.75,
             add = T)
  legend(x = "topright", legend = c("control","all treated"), col = c("black","cornflowerblue"), pch = c(19,1), bg = "transparent")
  
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "t" & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean"], vertical = T, pch = 1, method = "jitter", las = 2, cex.axis = 0.75,
             col = "cornflowerblue", main = paste0(project_name," NFA Mean Firing Rate AUC by Plate"))
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "n" & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean"], vertical = T, pch = 19, method = "jitter", las = 2, cex.axis = 0.75,
             add = T)
  legend(x = "topright", legend = c("control","all treated"), col = c("black","cornflowerblue"), pch = c(19,1), bg = "transparent")
  
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "t" & acsn == "CCTE_Shafer_MEA_dev_active_electrodes_number_DIV12"], vertical = T, pch = 1, method = "jitter", las = 2, cex.axis = 0.75,
             col = "cornflowerblue", main = paste0(project_name," NFA # Active Electrodes DIV12 by Plate"))
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "n" & acsn == "CCTE_Shafer_MEA_dev_active_electrodes_number_DIV12"], vertical = T, pch = 19, method = "jitter", las = 2, cex.axis = 0.75,
             add = T)
  legend(x = "topright", legend = c("control","all treated"), col = c("black","cornflowerblue"), pch = c(19,1), bg = "transparent")
  
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "t" & acsn == "CCTE_Shafer_MEA_dev_active_electrodes_number"], vertical = T, pch = 1, method = "jitter", las = 2, cex.axis = 0.75,
             col = "cornflowerblue", main = paste0(project_name," NFA # Active Electrodes AUC by Plate"))
  stripchart(rval ~ sub("_","\n",apid), dat[wllq == 1 & wllt == "n" & acsn == "CCTE_Shafer_MEA_dev_active_electrodes_number"], vertical = T, pch = 19, method = "jitter", las = 2, cex.axis = 0.75,
             add = T)
  legend(x = "topright", legend = c("control","all treated"), col = c("black","cornflowerblue"), pch = c(19,1), bg = "transparent")
  
  # define 'plotdat' - of the AUC MFR, with specialized conc group labels
  plotdat <- dat[acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean"]
  plotdat[, conc_grp := ifelse(wllt == "n",paste0(treatment,"\n",conc),signif(conc,1))]
  conc_grps <- unique(plotdat$conc_grp)
  plotdat$conc_grp <- factor(plotdat$conc_grp, levels = c(grep("\n",conc_grps,val = T),sort(unique(as.numeric(conc_grps[!grepl("\n",conc_grps)])))), ordered = T)
  
  # view all compounds together by dose
  stripchart(rval ~ conc_grp, plotdat[wllq == 1], vertical = T, pch = 1, method = "jitter", las = 2,
             main = paste0("Mean Firing Rate AUC by dose for all compounds in ",project_name), ylab = "CCTE_Shafer_MEA_dev_firing_rate_mean (AUC)", xlab = "conc")
  if (plotdat[, any(wllq==0)])
    stripchart(rval ~ conc_grp, plotdat[wllq == 0], vertical = T, pch = 1, method = "jitter",
               add = T, col = "red")
  legend(x = "topright", legend = c("wllq==1","wllq==0"), col = c("black","red"), pch = c(1,1), bg = "transparent")
  
  # find a compound that is likely to be a positive and plot dose response
  dat[, max_conc_by_spid := as.character(max(conc)), by = .(spid)]
  plot_spid <- dat[as.character(conc) == max_conc_by_spid & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean", .(med_rval = median(rval)), by = "spid"][med_rval == min(med_rval), spid[1]]
  plot_plates <- dat[spid == plot_spid, unique(apid)]
  stripchart(rval ~ conc_grp, plotdat[apid %in% plot_plates & (spid == plot_spid | wllt == "n") & wllq == 1], vertical = T, pch = 19, las = 2,
             col = rgb(0.1,0.1,0.1,0.5),
             ylim = range(dat[wllq == 1 & acsn == "CCTE_Shafer_MEA_dev_firing_rate_mean",rval]), ylab = "CCTE_Shafer_MEA_dev_firing_rate_mean (AUC)",
             xlab = "conc", main = paste0("Example Down Response:\n",dat[spid == plot_spid,unique(treatment)]," Mean Firing Rate AUC Dose Response"))
  if (plotdat[apid %in% plot_plates & (spid == plot_spid | wllt == "n"), any(wllq==0)])
    stripchart(rval ~ conc_grp, plotdat[apid %in% plot_plates & (spid == plot_spid | wllt == "n") & wllq == 0], vertical = T, pch = 19, las = 2,
               add = TRUE, col = rgb(0.9,0,0,0.5))
  legend(x = "topright", legend = c("wllq==1","wllq==0"), col = c(rgb(0.1,0.1,0.1,0.5),rgb(0.9,0,0,0.5)), pch = c(19,19), bg = "transparent")
  
  # Cytotox
  plotdat <- dat[grepl("(LDH)|(AB)",acsn)]
  plotdat[, conc_grp := ifelse(wllt == "n",paste0(treatment,"\n",conc),signif(conc,1))]
  conc_grps <- unique(plotdat$conc_grp)
  plotdat$conc_grp <- factor(plotdat$conc_grp, levels = c(grep("\n",conc_grps,val = T),sort(unique(as.numeric(conc_grps[!grepl("\n",conc_grps)])))), ordered = T)
  stripchart(rval ~ conc_grp, plotdat[wllq == 1 & grepl("AB",acsn)], las = 2,
             vertical = TRUE, pch = 1, method = "jitter", xlab = "conc", main = paste0("AB Blank-Corrected Values for ",project_name,"\nwhere wllq == 1"))
  if (nrow(plotdat[wllq == 1 & grepl("LDH",acsn)]) > 0) {
    stripchart(rval ~ conc_grp, plotdat[wllq == 1 & grepl("LDH",acsn)], las = 2,
               vertical = TRUE, pch = 1, method = "jitter", xlab = "conc", main = paste0("LDH Blank-Corrected Values for ",project_name,"\nwhere wllq == 1"))
  }
  
  if(remove_spid_col) dat[, spid := NULL]
  dat[, max_conc_by_spid := NULL]
  return(0)
 
}


# function designed for wide-format input data
dataset_checks_wide <- function(dat, normalized = FALSE, direction = '') {
  arg_table_name <- as.character(substitute(dat)) # in example, I saw deparse(substitute(dat)) as well
  dat <- as.data.table(dat)
  cat("\nFinal Checks for",arg_table_name,"\n")
  cat("Number of cultures dates:",dat[, length(unique(date))])
  cat("\nRange of culture dates:", dat[, range(date)] )
  cat("\nNumber of plates tested:",dat[, length(unique(Plate.SN))])
  cat("\nNumber of compounds tested:",dat[dose != 0, length(unique(treatment))])
  cat("\nNumber of NA values where wllq==1:\n")
  print(dat[wllq == 1, lapply(.SD, function(coli) sum(is.na(coli))), .SDcols = grep('_auc',names(dat))])
  
  cat("\nWllq breakdown for all points:\n")
  print(dat[, .N, by = "wllq"]) # note if wllq is NA anywhere
  
  cat("Plates that don't have exactly 48 wells:\n")
  print(dat[, .N, by = .(Plate.SN, date)][N != 48])
  
  # check number of controls per plate
  cat("\nPlates that don't have exactly 6 control wells:\n")
  print(dat[dose == 0 & wllq == 1, .N, by = c('date','Plate.SN')][N != 6])
  
  # PLOTS to visually confirm results
  if (normalized) {
    # adding section page to distinguish the normalized values
    plot.new()
    text(0.5, 0.5, labels = paste0(project_name, '\nNormalized ', direction, ' AUC Visualizations'), cex = 2)
  }
  
  # view all by plate
  dat[, apid := paste0(date,"_",Plate.SN)]
  stripchart(meanfiringrate_auc ~ sub("_","\n",apid), dat[wllq == 1 & dose != 0], vertical = T, pch = 1, method = "jitter", las = 2, cex.axis = 0.75,
             col = "cornflowerblue", main = paste0(project_name,if(normalized) paste0(" Normalized ",direction)," Mean Firing Rate AUC by Plate"))
  stripchart(meanfiringrate_auc ~ sub("_","\n",apid), dat[wllq == 1 & dose == 0], vertical = T, pch = 19, method = "jitter", las = 2, cex.axis = 0.75,
             add = T)
  legend(x = "topright", legend = c("control","all treated"), col = c("black","cornflowerblue"), pch = c(19,1), bg = "transparent")
  
  stripchart(nAE_auc ~ sub("_","\n",apid), dat[wllq == 1 & dose != 0], vertical = T, pch = 1, method = "jitter", las = 2, cex.axis = 0.75,
             col = "cornflowerblue", main = paste0(project_name,if(normalized) paste0(" Normalized ",direction)," # Active Electrodes AUC by Plate"))
  stripchart(nAE_auc ~ sub("_","\n",apid), dat[wllq == 1 & dose == 0], vertical = T, pch = 19, method = "jitter", las = 2, cex.axis = 0.75,
             add = T)
  legend(x = "topright", legend = c("control","all treated"), col = c("black","cornflowerblue"), pch = c(19,1), bg = "transparent")
  
  # define 'dat' - of the AUC MFR, with specialized conc group labels
  dat[, conc_grp := signif(dose,1)]
  conc_grps <- unique(dat$conc_grp)
  dat$conc_grp <- factor(dat$conc_grp, levels = sort(unique(dat$conc_grp)), ordered = T)
  
  # view all compounds together by dose
  stripchart(meanfiringrate_auc ~ conc_grp, dat[wllq == 1], vertical = T, pch = 1, method = "jitter", las = 2,
             main = paste0(if(normalized) paste0(" Normalized ",direction," "),"Mean Firing Rate AUC by dose for all compounds in ",project_name), ylab = paste0(if(normalized) paste0(" Normalized ",direction), 'Mean Firing Rate AUC'), xlab = "conc")
  if (dat[, any(wllq==0)])
    stripchart(meanfiringrate_auc ~ conc_grp, dat[wllq == 0], vertical = T, pch = 1, method = "jitter",
               add = T, col = "red")
  legend(x = "topright", legend = c("wllq==1","wllq==0"), col = c("black","red"), pch = c(1,1), bg = "transparent")
  
  # find a compound that is likely to be a positive and plot dose response
  if (grepl('[Uu]p',direction)) {
    plot_trt <- dat[dose == max(dose), .(med_rval = median(meanfiringrate_auc)), by = "treatment"][med_rval == max(med_rval), treatment[1]]
  }
  else {
    plot_trt <- dat[dose == max(dose), .(med_rval = median(meanfiringrate_auc)), by = "treatment"][med_rval == min(med_rval), treatment[1]]
  }
  plot_plates <- dat[treatment == plot_trt, unique(.SD), .SDcols = c('Plate.SN','date')]
  setkey(dat, Plate.SN, date)
  dat[J(plot_plates)]
  stripchart(meanfiringrate_auc ~ conc_grp, dat[J(plot_plates)][wllq == 1 & (treatment == plot_trt | dose == 0)], vertical = T, pch = 19, las = 2,
             col = rgb(0.1,0.1,0.1,0.5),
             ylim = range(dat[wllq == 1, meanfiringrate_auc]), ylab = paste0(if(normalized) paste0(" Normalized ",direction), 'Mean Firing Rate AUC'),
             xlab = "conc", main = paste0("trt: ",dat[treatment == plot_trt,unique(treatment)],if(normalized) paste0(" Normalized ",direction)," Mean Firing Rate AUC Dose Response"))
  if (dat[J(plot_plates)][(treatment == plot_trt | dose == 0), any(wllq == 0)])
    stripchart(rval ~ conc_grp, dat[J(plot_plates)][wllq == 0 & (treatment == plot_trt | dose == 0)], vertical = T, pch = 19, las = 2,
               add = TRUE, col = rgb(0.9,0,0,0.5))
  legend(x = "topright", legend = c("wllq==1","wllq==0"), col = c(rgb(0.1,0.1,0.1,0.5),rgb(0.9,0,0,0.5)), pch = c(19,19), bg = "transparent")
  
}
