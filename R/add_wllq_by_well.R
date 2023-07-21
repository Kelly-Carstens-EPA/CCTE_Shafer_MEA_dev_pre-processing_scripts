

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
  wllq.tb.by.well <- rbind(wllq.tb.by.well[!well %in% 'all'],
                           wllq.tb.full.plates, fill = T)
  wllq.tb.by.well[!well %in% 'all', rowi := match(stri_extract_first(toupper(well), regex = '[A-Z]{1,}'), LETTERS)]
  wllq.tb.by.well[!well %in% 'all', coli := as.numeric(stri_extract_first(well, regex = '[0-9]{1,}'))]
  wllq.tb.by.well[, well := NULL]
  
  
  # Expand rows where DIV == 'all' (if applicable)
  if ('DIV' %in% names(wllq.tb.by.well)) {
    setnames(wllq.tb.by.well, old = 'DIV', new = 'DIV_org')
    indicies.all.DIVs <- which(wllq.tb.by.well$DIV_org == 'all')
    all.DIVs.tb <- data.table('DIV' = all_DIVs)
    wllq.tb.full.DIVs <- rbindlist(lapply(indicies.all.DIVs, 
                                          function(i) cbind(wllq.tb.by.well[i], all.DIVs.tb)))
    wllq.tb.by.well <- rbind(wllq.tb.by.well[!DIV_org %in% 'all'],
                             wllq.tb.full.DIVs, fill = T)    
    wllq.tb.by.well[!DIV_org %in% 'all', DIV := as.numeric(DIV_org)]
    stopifnot(nrow(wllq.tb.by.well[assay %in% 'nfa' & is.na(DIV)]) == 0)
    wllq.tb.by.well[, DIV_org := NULL]
  }

  
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
    warning(paste0('Some rows in ',basename(wllq.tb.by.well.file), ' did not match any rows in dat'))
  }
  dat[, c('dat.id','wllq.tb.by.well.id') := NULL]
  
  # Default to wllq_by_well == 1 were not specified in wllq.tb.by.well
  dat[is.na(wllq_by_well), wllq_by_well := 1]
  
  dat
}