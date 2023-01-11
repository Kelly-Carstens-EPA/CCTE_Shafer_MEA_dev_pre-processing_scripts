confirm_concs <- function(dat, spidmap, con, expected_target_concs = c(0.03,0.1,0.3,1,3,10,30), update_concs_without_prompt = FALSE) {
  
  require(RMySQL)
  
  cat("Checking conc's:\n")
  
  # check if any conc's are NA
  if (dat[wllt == 't', any(is.na(conc))]) {
    stop(paste("\nThe following treatments have conc NA:",paste0(dat[is.na(conc), unique(treatment)],collapse=",")))
  }
  
  # get the stck from invitrodb - sometimes spidmap has duplicate stock_conc entries
  query_term <- paste0("SELECT * FROM sample WHERE spid IN('",paste(spidmap[!is.na(spid),unique(spid)],collapse="','",sep=""),"');")
  sample_info <- dbGetQuery(con, query_term)
  dbDisconnect(con)
  spidmap <- merge(spidmap, sample_info, by = "spid", all.x = T)
  if (nrow(spidmap[is.na(stkc)]) > 0) {
    cat(paste0("The following spids were not found in invitrodb: ",paste0(spidmap[is.na(stkc),unique(spid)],collapse=", "),
                   "\nWill use 'stock_conc' provided in spidmap file instead."))
    spidmap[is.na(stkc), stkc := stock_conc]
  }
  # there can be multiple stock_concs listed for a given spid in file, so eliminating duplicates here and using invitrodb stkc
  spidmap <- spidmap[, unique(.SD), .SDcols = c("spid","treatment","stkc","stkc_unit","expected_stock_conc")] 
  # I am trusting that there is only 1 stkc for each spid lsited in invitrodb, so I won't check for that
  
  # compare the concentrations
  cat("\nAll compounds are assumed to have conc's",expected_target_concs,"\n(You can change this by updating the argument 'expected_target_concs' of the function confirm_concs()).\n")
  compare_concs <- merge(spidmap[, .(stkc, expected_stock_conc, spidmap_guess_concs = paste0(signif(stkc/expected_stock_conc*expected_target_concs,3),collapse=",")), by = "spid"],
                         dat[wllt == "t", .(source_concs = paste0(sort(unique(signif(conc,3))),collapse=","), num_concs = length(unique(conc))), by = c("spid","treatment")], 
                         by = "spid", all.y = TRUE)
  cat('\nFYI, the following stock conc\'s pulled from invitrodb do not match the expected stock conc:\n')
  print(compare_concs[signif(stkc,4) != signif(expected_stock_conc,4), .(spid, treatment, stkc, expected_stock_conc, spidmap_guess_concs, source_concs, num_concs)])
  
  if(nrow(compare_concs[source_concs != spidmap_guess_concs | is.na(spidmap_guess_concs)]) > 0) {
    cat("The concentrations for the following compounds might need to be corrected:\n")
    # removing this feature for now -> user can do manually if needed, and check that it is correct
    # compare_concs$probably_partially_conc_corrected <- sapply(strsplit(compare_concs$source_concs,split=","), function(x) length(x) > length(expected_target_concs))
    print(compare_concs[source_concs != spidmap_guess_concs | is.na(spidmap_guess_concs)][order(num_concs)])
    
    if (update_concs_without_prompt) response <- "y"
    else response <- readline(prompt = "Update conc's where source_concs != spidmap_guess_concs with conc := signif(stkc/expected_stock_conc*source_concs, 3)? (y/n): ")
    
    if (response %in% c("y","Y","yes","Yes")) {
      
      compare_concs[source_concs != spidmap_guess_concs | is.na(spidmap_guess_concs), need_to_update_concs := TRUE]
      
      dat[, conc_org := conc]
      dat <- merge(dat, compare_concs[, .(spid, need_to_update_concs, stkc, expected_stock_conc)], by = "spid", all.x = TRUE)
      
      # for spid's with 'probably_partially_conc_corrected', standardize the conc's first:
      # cat("Standardizing concs where 'probably_partially_conc_corrected'==TRUE...\n")
      # dat[probably_partially_conc_corrected == TRUE, conc := signif(conc, digits = 1)]
      # dat[, conc_standardized := conc]
  
      # now correct the conc's
      cat("Correcting conc's...\n")
      dat[need_to_update_concs == TRUE, conc := signif(stkc/expected_stock_conc*conc, 3)]
      update_summary <- dat[need_to_update_concs == TRUE, .(stkc, concs_in_source_dat = paste0(unique(conc_org),collapse=", ")),
                            by = c("spid","treatment","conc","stkc","expected_stock_conc")][order(spid,conc), .(treatment, spid, stkc, expected_stock_conc, 
                                                                                                                concs_in_source_dat, conc_updated = format(conc,digits=4,scientific=F))]
      if (update_concs_without_prompt) {
        cat("conc's that changed:\n")
        print(update_summary[signif(as.numeric(concs_in_source_dat),3) != signif(as.numeric(conc_updated),3)]) # display conc's that were actually updated
      }
      else {
        # assign("update_summary",update_summary, envir = .GlobalEnv)
        cat("View the table 'update_summary' to confirm that the concentration-corrections are correct.\n")
        cat("If it looks correct, enter c to continue. Else Q to quit and fix.\n")
        browser()
        # response <- readline(prompt = "Does conc correction look correct for each compound and dose? (y/n): ")
        # if (!(response %in% c("y","Y","yes","Yes"))) browser()
      }
      
      dat[, c("conc_org","stkc","need_to_update_concs","expected_stock_conc") := list(NULL)] # remove added columns
    }
    else {
      stop("Update conc's, then re-run\n")
    }
  }
  else {
    cat("\nAll compounds have the expected concentration-corrected values\n")
  }
  return(dat)
}
