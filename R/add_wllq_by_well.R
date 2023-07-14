

add_wllq_by_well <- function(dat, wllq.tb.by.well.file, num_rows_per_plate, num_columns_per_plate, all_DIVs = c(5,7,9,12)) {
  
  # read wllq table
  wllq.tb.by.well <- as.data.table(read.csv(wllq.tb.by.well.file, colClasses = 'character'))
  
  # Check if 'assay' column is present, to match with affected_assays in wllq table
  stopifnot('assay' %in% names(dat))
  
  # Expand rows for each "affected assay"
  # only include 'assays' that are present in dat
  dat[, assay := tolower(assay)]
  wllq.tb.by.well[, affected_assays := tolower(affected_assays)]
  affected_assays_list <- unique(unlist(stringi::stri_split(wllq.tb.by.well$affected_assays, regex = '[,;][ ]*')))
  affected_assays_list <- intersect(affected_assays_list, unique(dat$assay))
  if (length(affected_assays_list) == 0) {
    cat('no notes in ',basename(wllq.tb.by.well.file),
        ' apply to the assays in dat (',
        paste0(sort(unique(dat$assay)),collapse = ","),
        '). Defaulting to wllq_by_well = 1\n')
    dat[, `:=`(wllq_by_well = 1,
               wllq_notes_by_well = NA_character_,
               wllq_ref_by_well = NA_character_)]
    return(dat)
  }
  # expand affected_assays into table with 1 assay per row
  wllq.tb.by.well <- rbindlist(lapply(affected_assays_list,
                                      function(assay_typei) cbind(wllq.tb.by.well[grepl(assay_typei,affected_assays)],
                                                                  'assay' = assay_typei)))
  wllq.tb.by.well[, affected_assays := NULL]
  
  # Expand rows where well == 'all' (for a given plate)
  indicies.all.wells <- which(wllq.tb.by.well$well == 'all')
  all.wells.tb <- data.table('rowi' = rep(1:num_rows_per_plate, each = num_columns_per_plate),
                             'coli' = rep(1:num_columns_per_plate, times = num_rows_per_plate))
  wllq.tb.full.plates <- rbindlist(lapply(indicies.all.wells, 
                                          function(i) cbind(wllq.tb.by.well[i], all.wells.tb)))
  wllq.tb.by.well <- rbind(wllq.tb.by.well[well != 'all'],
                           wllq.tb.full.plates, fill = T)
  
  # For the wells specifically defined in wllq.tb, get the rowi and coli
  wllq.tb.by.well[is.na(rowi), rowi := match(stri_extract_first(toupper(well), regex = '[A-Z]{1,}'), LETTERS)]
  wllq.tb.by.well[is.na(coli), coli := as.numeric(stri_extract_first(well, regex = '[0-9]{1,}'))]
  wllq.tb.by.well[, well := NULL]
  
  
  # Expand rows where DIV == 'all' (if applicable)
  if ('DIV' %in% names(wllq.tb.by.well)) {
    setnames(wllq.tb.by.well, old = 'DIV', new = 'DIV_org')
    indicies.all.DIVs <- which(wllq.tb.by.well$DIV_org == 'all')
    all.DIVs.tb <- data.table('DIV' = all_DIVs)
    wllq.tb.full.DIVs <- rbindlist(lapply(indicies.all.DIVs, 
                                          function(i) cbind(wllq.tb.by.well[i], all.DIVs.tb)))
    wllq.tb.by.well <- rbind(wllq.tb.by.well[DIV_org != 'all'],
                             wllq.tb.full.DIVs, fill = T)    
  }
  wllq.tb.by.well[is.na(DIV), DIV := as.numeric(DIV_org)]
  stopifnot(nrow(wllq.tb.by.well[is.na(DIV)]) == 0)
  wllq.tb.by.well[, DIV_org := NULL]
  
  
  # If there are multiple rows in wllq.tb that apply to the same well/assay/DIV,
  # Get the minimum wllq and collapse across all notes
  id.columns <- names(wllq.tb.by.well)[!grepl('wllq',names(wllq.tb.by.well))]
  wllq.tb.by.well <- wllq.tb.by.well[, .(wllq_by_well = min(wllq),
                                         wllq_notes_by_well = paste0(unique(wllq_notes), collapse = "; "),
                                         wllq_ref_by_well = paste0(unique(wllq_ref), collapse = "; ")),
                                     by = c(id.columns)]
  stopifnot(nrow(wllq.tb.by.well[, .N, by = c(id.columns)][N > 1]) == 0) # checking 1 row per id.columns combination
  
  # If not present, add rowi and coli to dat
  if (!'coli' %in% names(dat))
    dat[, `:=`(rowi = match(stri_extract_first(toupper(well), regex = '[A-Z]{1,}'), LETTERS),
               coli = as.numeric(stri_extract_first(well, regex = '[0-9]{1,}')))]
  
  # Make all id.columns the same type in both tables, to ensure that they will merge
  dat[, c(id.columns) := lapply(.SD, as.character), .SDcols = id.columns]
  wllq.tb.by.well[, c(id.columns) := lapply(.SD, as.character), .SDcols = id.columns]
  
  # Merge wllq.tb with dat
  dat[, dat.id := 1]
  wllq.tb.by.well[, wllq.tb.by.well.id := 1]
  dat <- merge(dat, wllq.tb.by.well,
               by = id.columns, all = T)
  
  # Check for any rows from wllq.tb.by.well that didn't match a row in dat
  # (indicates possible type-o in wllq.tb)
  if(nrow(dat[is.na(dat.id)]) > 0) {
    cat('\nSome rows in ',basename(wllq.tb.by.well.file), ' did not match any rows in dat:\n')
    print(dat[is.na(dat.id)])
    stop(paste0('Some rows in ',basename(wllq.tb.by.well.file), ' did not match any rows in dat'))
  }
  dat[, c('dat.id','wllq.tb.by.well.id') := NULL]
  
  # Default to wllq_by_well == 1 were not specified in wllq.tb.by.well
  dat[is.na(wllq_by_well), wllq_by_well := 1]
  
  dat
}



# Add wllq notes by compound + cndx ------------------------------------


wllq.tb.by.trt.file <- paste0('projects/well_quality_tables/',project_name,'_well_quality_assignments_by_compound_cndx_culture_date.csv')
if (file.exists(wllq.tb.by.trt.file)) {
  
  wllq.tb.by.trt <- as.data.table(read.csv(wllq.tb.by.trt.file, colClasses = 'character'))
  wllq.tb.by.trt[culture_date == '', culture_date := NA_character_]
  
  # Rename compound name column to match how it appears in dat
  setnames(wllq.tb.by.trt, old = 'compound_name', new = 'Plate.Position') # because that is the compound name column in this particular data set
  
  # Expand where cndx is "all"
  # (assuming same cndx's for all compounds)
  setnames(wllq.tb.by.trt, old = 'cndx', new = 'cndx_org')
  indicies.all.cndx <- which(wllq.tb.by.trt$cndx_org == 'all')
  full.cndx.series <- data.table('cndx' = sort(unique(dat$cndx)))
  wllq.tb.full.cndx <- rbindlist(lapply(indicies.all.cndx, function(i) cbind(wllq.tb.by.trt[i], full.cndx.series)))
  wllq.tb.by.trt[, cndx := cndx_org]
  wllq.tb.by.trt <- rbind(wllq.tb.by.trt[cndx_org != 'all'],
                          wllq.tb.full.cndx)
  wllq.tb.by.trt[, cndx_org := NULL]
  wllq.tb.by.trt[, cndx := as.numeric(cndx)]
  
  # Expand rows for each "affected assay"
  affected_assays_list <- unique(unlist(stringi::stri_split(wllq.tb.by.trt$affected_assays, regex = '[,;][ ]*')))
  wllq.tb.by.trt <- rbindlist(lapply(affected_assays_list,
                                     function(assay_typei) cbind(wllq.tb.by.trt[grepl(assay_typei,affected_assays), 
                                                                                .SD, .SDcols = setdiff(names(wllq.tb.by.trt),'affected_assays')],
                                                                 'assay' = assay_typei)))
  
  # Get the minimum wllq and collapse across notes for every date, treatment, and cndx, and assay
  wllq.tb.by.trt <- wllq.tb.by.trt[, .(wllq_by_trt = min(wllq),
                                       wllq_notes_by_trt = paste0(unique(wllq_notes), collapse = "; "),
                                       wllq_ref_by_trt = paste0(unique(wllq_ref), collapse = "; ")),
                                   by = c(names(wllq.tb.by.trt)[!grepl('wllq',names(wllq.tb.by.trt))])]
  stopifnot(nrow(wllq.tb.by.trt[, .N, by = .(culture_date, Plate.Position, cndx, assay)][N > 1]) == 0) # empty, good
  
  # Prepare to merge with dat
  dat[, dat.id := 1]
  wllq.tb.by.trt[, wllq.by.trt.id := .I]
  
  # Merge with dat where culture_date is not defined (where wllq note applies to all cultures for given treatment and cndx)
  dat <- merge(dat, wllq.tb.by.trt[is.na(culture_date) | culture_date %in% c('',' '), .SD, .SDcols = setdiff(names(wllq.tb.by.trt),'culture_date')],
               by = c('assay','Plate.Position','cndx'), all = T)
  # Any rows from wllq.tb.by.trt that didn't match a row in dat (based on columns in "by")?
  stopifnot(nrow(dat[is.na(dat.id)]) == 0) # should be empty
  
  # Add flag if a wllq note was applied to multiple cultures
  dat <- dat[order(culture_date)]
  check.tb <- dat[!is.na(wllq.by.trt.id), .(num_cultures_tested = length(unique(culture_date)), 
                                            culture_dates = paste0(unique(culture_date),collapse = ", "),
                                            conc_list = paste0(unique(conc),collapse = ', ')), 
                  by = .(Plate.Position, cndx, wllq_notes_by_trt_short = paste0(stri_sub(wllq_notes_by_trt, from  = 1, to = 20),'...'))][num_cultures_tested > 1][order(Plate.Position, cndx)]
  if (nrow(check.tb) > 0) {
    cat('Note that multiple cultures and (possibly multiple concentrations) are affected by\nthe wllq notes by compound-cndx for the following:\n\n')
    print(check.tb)
    cat('\nIf these wllq_notes should not apply to every culture and/or conc shown above,\nthen specify the culture_date in the well_quality_assignments_by_compound_cndx_culture_date.csv table.\n\n')
  }
  
  # Merge with dat where culture_date is defined
  dat <- merge(dat, wllq.tb.by.trt[!(is.na(culture_date) | culture_date %in% c('',' '))],
               by = c('culture_date','assay','Plate.Position','cndx'), all = T)
  # Any rows from wllq.tb.by.trt that didn't match a row in dat (based on columns in "by")?
  stopifnot(nrow(dat[is.na(dat.id)]) == 0) # should be empty
  
  # Merge wllq_ from merge 1 and merge 2, make sure merge all wllq notes and get min wllq for each well
  # Combine the wllq_ columns from 2 merges, get the minimum wllq
  dat[, `:=`(wllq_by_trt = pmin(wllq_by_trt.x, wllq_by_trt.y, na.rm = T),
             wllq_notes_by_trt = paste0(ifelse(!is.na(wllq_notes_by_trt.x),wllq_notes_by_trt.x,''),
                                        ifelse(!is.na(wllq_notes_by_trt.x) & !is.na(wllq_notes_by_trt.y),'; ',''),
                                        ifelse(!is.na(wllq_notes_by_trt.y),wllq_notes_by_trt.y,'')),
             wllq_ref_by_trt = paste0(ifelse(!is.na(wllq_ref_by_trt.x),wllq_ref_by_trt.x,''),
                                      ifelse(!is.na(wllq_ref_by_trt.x) & !is.na(wllq_ref_by_trt.y),'; ',''),
                                      ifelse(!is.na(wllq_ref_by_trt.y),wllq_ref_by_trt.y,''))
  )]
  # Make blank values NA to standardize
  dat[wllq_notes_by_trt == '', wllq_notes_by_trt := NA_character_]
  dat[wllq_ref_by_trt == '', wllq_ref_by_trt := NA_character_]
  # confirm no rows duplicated
  stopifnot(nrow(dat[, .N, by = .(group_integer, culture_date, plate, rowi, coli, acnm)][N > 1]) == 0) # empty, good
  
  # Give some output to visually confirm what happened
  cat('wllq updates by treatment_cndx.\nConfirm that wllqs apply to anticipated treatment-conc-test group/date combinations:\n\n')
  print(dat[!is.na(wllq_by_trt), .(num_pts = .N), by = .(Plate.Position, group_integer, culture_date, cndx, conc, wllq_by_trt, wllq_notes_by_trt_short = stri_sub(wllq_notes_by_trt, from = 1, to = 20), wllq_ref_by_trt)][order(Plate.Position, group_integer, culture_date, as.numeric(cndx), conc)])
  
  # Remove columns no longer needed
  dat[, paste0(c('wllq_by_trt','wllq_notes_by_trt','wllq_ref_by_trt','wllq.by.trt.id'),rep(c('.x','.y'), each = 4)) := NULL]
  dat[, c('dat.id') := NULL]
  
} else {
  cat(paste0('projects/',project_name,'_well_quality_assignments_by_culture_treatment.csv'),' does not exist. wllq_by_trt will be set to NA.\n')
  dat[, c('wllq_by_trt') := NA_integer_]
  dat[, c('wllq_notes_by_trt','wllq_ref_by_trt') := NA_character_] 
}

