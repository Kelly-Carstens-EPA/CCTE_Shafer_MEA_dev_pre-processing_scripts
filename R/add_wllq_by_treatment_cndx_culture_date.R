add_wllq_by_treatment_cndx_culture_date <- function(dat, 
                                                    wllq.tb.by.trt.file) {
  
  # read wllq table
  wllq.tb.by.trt <- as.data.table(read.csv(wllq.tb.by.trt.file, colClasses = 'character'))
  
  # Check if 'assay' column is present, to match with affected_assays in wllq table
  stopifnot('assay' %in% names(dat))
  
  # Expand rows for each "affected assay"
  # only include 'assays' that are present in dat
  dat[, assay := toupper(assay)]
  wllq.tb.by.trt[, affected_assays := toupper(affected_assays)]
  affected_assays_list <- unique(unlist(stringi::stri_split(wllq.tb.by.trt$affected_assays, regex = '[,;][ ]*')))
  affected_assays_list <- intersect(affected_assays_list, unique(dat$assay))
  if (length(affected_assays_list) == 0) {
    cat('no notes in ',basename(wllq.tb.by.trt.file),
        ' apply to the assays in dat (',
        paste0(sort(unique(dat$assay)),collapse = ","),
        '). Defaulting to wllq_by_trt = 1\n')
    dat[, `:=`(wllq_by_trt = 1,
               wllq_notes_by_trt = NA_character_,
               wllq_ref_by_trt = NA_character_)]
    return(dat)
  }
  # expand affected_assays into table with 1 assay per row
  wllq.tb.by.trt <- rbindlist(lapply(affected_assays_list,
                                     function(assay_typei) cbind(wllq.tb.by.trt[grepl(assay_typei,affected_assays)],
                                                                 'assay' = assay_typei)))
  wllq.tb.by.trt[, affected_assays := NULL]
  
  
  # If there are multiple rows in wllq.tb that apply to the same well/assay/DIV,
  # Get the minimum wllq and collapse across all notes
  id.columns <- names(wllq.tb.by.trt)[!grepl('wllq',names(wllq.tb.by.trt))]
  wllq.tb.by.trt <- wllq.tb.by.trt[, .(wllq_by_trt = min(wllq_by_trt),
                                       wllq_notes_by_trt = paste0(unique(wllq_notes_by_trt), collapse = "; "),
                                       wllq_ref_by_trt = paste0(unique(wllq_ref_by_trt), collapse = "; ")),
                                   by = c(id.columns)]
  stopifnot(nrow(wllq.tb.by.trt[, .N, by = c(id.columns)][N > 1]) == 0) # checking 1 row per id.columns combination
  
  # Make all id.columns the same type in both tables, to ensure that they will merge
  dat[, c(id.columns) := lapply(.SD, as.character), .SDcols = id.columns]
  wllq.tb.by.trt[, c(id.columns) := lapply(.SD, as.character), .SDcols = id.columns]  
  
  # Merge wllq.tb with dat
  dat[, dat.id := 1]
  wllq.tb.by.trt[, wllq.tb.by.trt.id := 1]
  dat <- merge(dat, wllq.tb.by.trt,
               by = id.columns, all = T)
  
  # Check for any rows from wllq.tb.by.trt that didn't match a row in dat
  # (indicates possible type-o in wllq.tb)
  if(nrow(dat[is.na(dat.id)]) > 0) {
    cat('\nSome rows in ',basename(wllq.tb.by.trt.file), ' did not match any rows in dat:\n')
    print(dat[is.na(dat.id), .N, by = c(id.columns)])
    warning(paste0('Some rows in ',basename(wllq.tb.by.trt.file), ' did not match any rows in dat. These will be ignored.'))
    dat <- dat[!is.na(dat.id)]
  }
  dat[, c('dat.id','wllq.tb.by.trt.id') := NULL]
  
  # Default to wllq_by_trt == 1 where not specified in wllq.tb.by.trt
  dat[is.na(wllq_by_trt), wllq_by_trt := 1]
  
  dat
}