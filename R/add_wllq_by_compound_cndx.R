add_wllq_by_compound_cndx <- function() {
  # 
  # 
  # # Add wllq notes by compound + cndx ------------------------------------
  # 
  # 
  # wllq.tb.by.trt.file <- paste0('projects/well_quality_tables/',project_name,'_well_quality_assignments_by_compound_cndx_culture_date.csv')
  # if (file.exists(wllq.tb.by.trt.file)) {
  #   
  #   wllq.tb.by.trt <- as.data.table(read.csv(wllq.tb.by.trt.file, colClasses = 'character'))
  #   wllq.tb.by.trt[culture_date == '', culture_date := NA_character_]
  #   
  #   # Rename compound name column to match how it appears in dat
  #   setnames(wllq.tb.by.trt, old = 'compound_name', new = 'Plate.Position') # because that is the compound name column in this particular data set
  #   
  #   # Expand where cndx is "all"
  #   # (assuming same cndx's for all compounds)
  #   setnames(wllq.tb.by.trt, old = 'cndx', new = 'cndx_org')
  #   indicies.all.cndx <- which(wllq.tb.by.trt$cndx_org == 'all')
  #   full.cndx.series <- data.table('cndx' = sort(unique(dat$cndx)))
  #   wllq.tb.full.cndx <- rbindlist(lapply(indicies.all.cndx, function(i) cbind(wllq.tb.by.trt[i], full.cndx.series)))
  #   wllq.tb.by.trt[, cndx := cndx_org]
  #   wllq.tb.by.trt <- rbind(wllq.tb.by.trt[cndx_org != 'all'],
  #                           wllq.tb.full.cndx)
  #   wllq.tb.by.trt[, cndx_org := NULL]
  #   wllq.tb.by.trt[, cndx := as.numeric(cndx)]
  #   
  #   # Expand rows for each "affected assay"
  #   affected_assays_list <- unique(unlist(stringi::stri_split(wllq.tb.by.trt$affected_assays, regex = '[,;][ ]*')))
  #   wllq.tb.by.trt <- rbindlist(lapply(affected_assays_list,
  #                                      function(assay_typei) cbind(wllq.tb.by.trt[grepl(assay_typei,affected_assays), 
  #                                                                                 .SD, .SDcols = setdiff(names(wllq.tb.by.trt),'affected_assays')],
  #                                                                  'assay' = assay_typei)))
  #   
  #   # Get the minimum wllq and collapse across notes for every date, treatment, and cndx, and assay
  #   wllq.tb.by.trt <- wllq.tb.by.trt[, .(wllq_by_trt = min(wllq),
  #                                        wllq_notes_by_trt = paste0(unique(wllq_notes), collapse = "; "),
  #                                        wllq_ref_by_trt = paste0(unique(wllq_ref), collapse = "; ")),
  #                                    by = c(names(wllq.tb.by.trt)[!grepl('wllq',names(wllq.tb.by.trt))])]
  #   stopifnot(nrow(wllq.tb.by.trt[, .N, by = .(culture_date, Plate.Position, cndx, assay)][N > 1]) == 0) # empty, good
  #   
  #   # Prepare to merge with dat
  #   dat[, dat.id := 1]
  #   wllq.tb.by.trt[, wllq.by.trt.id := .I]
  #   
  #   # Merge with dat where culture_date is not defined (where wllq note applies to all cultures for given treatment and cndx)
  #   dat <- merge(dat, wllq.tb.by.trt[is.na(culture_date) | culture_date %in% c('',' '), .SD, .SDcols = setdiff(names(wllq.tb.by.trt),'culture_date')],
  #                by = c('assay','Plate.Position','cndx'), all = T)
  #   # Any rows from wllq.tb.by.trt that didn't match a row in dat (based on columns in "by")?
  #   stopifnot(nrow(dat[is.na(dat.id)]) == 0) # should be empty
  #   
  #   # Add flag if a wllq note was applied to multiple cultures
  #   dat <- dat[order(culture_date)]
  #   check.tb <- dat[!is.na(wllq.by.trt.id), .(num_cultures_tested = length(unique(culture_date)), 
  #                                             culture_dates = paste0(unique(culture_date),collapse = ", "),
  #                                             conc_list = paste0(unique(conc),collapse = ', ')), 
  #                   by = .(Plate.Position, cndx, wllq_notes_by_trt_short = paste0(stri_sub(wllq_notes_by_trt, from  = 1, to = 20),'...'))][num_cultures_tested > 1][order(Plate.Position, cndx)]
  #   if (nrow(check.tb) > 0) {
  #     cat('Note that multiple cultures and (possibly multiple concentrations) are affected by\nthe wllq notes by compound-cndx for the following:\n\n')
  #     print(check.tb)
  #     cat('\nIf these wllq_notes should not apply to every culture and/or conc shown above,\nthen specify the culture_date in the well_quality_assignments_by_compound_cndx_culture_date.csv table.\n\n')
  #   }
  #   
  #   # Merge with dat where culture_date is defined
  #   dat <- merge(dat, wllq.tb.by.trt[!(is.na(culture_date) | culture_date %in% c('',' '))],
  #                by = c('culture_date','assay','Plate.Position','cndx'), all = T)
  #   # Any rows from wllq.tb.by.trt that didn't match a row in dat (based on columns in "by")?
  #   stopifnot(nrow(dat[is.na(dat.id)]) == 0) # should be empty
  #   
  #   # Merge wllq_ from merge 1 and merge 2, make sure merge all wllq notes and get min wllq for each well
  #   # Combine the wllq_ columns from 2 merges, get the minimum wllq
  #   dat[, `:=`(wllq_by_trt = pmin(wllq_by_trt.x, wllq_by_trt.y, na.rm = T),
  #              wllq_notes_by_trt = paste0(ifelse(!is.na(wllq_notes_by_trt.x),wllq_notes_by_trt.x,''),
  #                                         ifelse(!is.na(wllq_notes_by_trt.x) & !is.na(wllq_notes_by_trt.y),'; ',''),
  #                                         ifelse(!is.na(wllq_notes_by_trt.y),wllq_notes_by_trt.y,'')),
  #              wllq_ref_by_trt = paste0(ifelse(!is.na(wllq_ref_by_trt.x),wllq_ref_by_trt.x,''),
  #                                       ifelse(!is.na(wllq_ref_by_trt.x) & !is.na(wllq_ref_by_trt.y),'; ',''),
  #                                       ifelse(!is.na(wllq_ref_by_trt.y),wllq_ref_by_trt.y,''))
  #   )]
  #   # Make blank values NA to standardize
  #   dat[wllq_notes_by_trt == '', wllq_notes_by_trt := NA_character_]
  #   dat[wllq_ref_by_trt == '', wllq_ref_by_trt := NA_character_]
  #   # confirm no rows duplicated
  #   stopifnot(nrow(dat[, .N, by = .(group_integer, culture_date, plate, rowi, coli, acnm)][N > 1]) == 0) # empty, good
  #   
  #   # Give some output to visually confirm what happened
  #   cat('wllq updates by treatment_cndx.\nConfirm that wllqs apply to anticipated treatment-conc-test group/date combinations:\n\n')
  #   print(dat[!is.na(wllq_by_trt), .(num_pts = .N), by = .(Plate.Position, group_integer, culture_date, cndx, conc, wllq_by_trt, wllq_notes_by_trt_short = stri_sub(wllq_notes_by_trt, from = 1, to = 20), wllq_ref_by_trt)][order(Plate.Position, group_integer, culture_date, as.numeric(cndx), conc)])
  #   
  #   # Remove columns no longer needed
  #   dat[, paste0(c('wllq_by_trt','wllq_notes_by_trt','wllq_ref_by_trt','wllq.by.trt.id'),rep(c('.x','.y'), each = 4)) := NULL]
  #   dat[, c('dat.id') := NULL]
  #   
  # } else {
  #   cat(paste0('projects/',project_name,'_well_quality_assignments_by_culture_treatment.csv'),' does not exist. wllq_by_trt will be set to NA.\n')
  #   dat[, c('wllq_by_trt') := NA_integer_]
  #   dat[, c('wllq_notes_by_trt','wllq_ref_by_trt') := NA_character_] 
  # }
  # 
  
}