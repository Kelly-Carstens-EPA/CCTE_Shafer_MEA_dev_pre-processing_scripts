# this script will run the first main steps of the data preparation
# it will check if a step has already been completed before doing it
# it will select files, create h5files, create prepared data, AUC table, and collect cytotox data

check_existing <- function(path, pattern, pause_between_steps) {
  # save_objects <- c("check_existing","main.output.dir",
  #                   "project_name","pause_between_steps","save_notes_graphs",
  #                   "default_ControlTreatmentName",
  #                   "different_vehicleControlCompounds",
  #                   "different_vehicleControls",
  #                   "spidmap_file",
  #                   "spid_sheet",
  #                   "scripts.dir",
  #                   "root_output_dir")
  # rm(list = setdiff(ls(parent.frame()), save_objects), envir = parent.frame()) # free up memory...
  if (pause_between_steps) {
    next_step <- readline(prompt = "Continue? (y/n): ")
    if (next_step == "n") {
      stop("User elected to stop.")
    }
  }
  resp <- "r" # default response if there are no existing files for this step
  num_files <- length(list.files(path, pattern, recursive = F)) # check if any of these files already exist
  if (num_files > 0) {
    cat(num_files, "files already exist.")
    if (num_files == 1) cat(" (",list.files(path, pattern, recursive = F),")")
    if (pause_between_steps) 
      repeat {
        resp <- readline(prompt = "Do you want to Continue with only the current files, Remake all files, Append to existing files, or Quit? (c/r/a/q): ")
        if (resp %in% c("c","r","a","q")) break
      }
    else 
      resp <- "c"
  }
  if (resp == "q") {
    stop("No error - user elected to stop.")
  }
  assign("append",switch(resp,
                         "a" = TRUE,
                         "r" = FALSE,
                         "c" = FALSE),
         envir = parent.frame())
  return(resp)
}

setwd(scripts.dir)
main.output.dir <- file.path(root_output_dir, project_name)
if (!dir.exists(main.output.dir)) dir.create(main.output.dir)

# h5_conversion.R
cat("\n- Create h5 files:\n")
resp <- check_existing(path = file.path(main.output.dir,"h5files"), pattern = "\\.h5", pause_between_steps)
if (resp %in% c("r","a")) {
  source('spike_list_functions.R')
  remake_all <- !append
  source('h5_conversion.R')
  cat("h5files are ready in folder",file.path(main.output.dir,"h5files"),"\n")
  rm(list = setdiff(ls(), keep_items))
}

# create_ont_csv
cat("\n- Calculate the components:\n")
resp <- check_existing(path = file.path(main.output.dir,"prepared_data"), pattern = "\\.csv", pause_between_steps)
if (resp %in% c("r","a")) {
  source('create_ont_csv.R')
  source('create_burst_ont_Data.R')
  source('local.corr.all.ont.ae.filter.R')
  create_ont_csv(basepath = main.output.dir, get_h5Files_under_basepath = TRUE, remake_all = !append)
  rm(list = setdiff(ls(), keep_items))
}

# normalized mutual information calculation
cat("\n- Calculate the Mutual Information:\n")
resp <- check_existing(path = file.path(main.output.dir,"All_MI"), pattern = "\\.csv", pause_between_steps)
if (resp %in% c("r","a")) {
  source('spikeLoadRoutines.R')
  source('nmi2_final.R')
  source('nmi_wrapper.R')
  source('MI_script_all.R')
  run_mi_functions(basepath = main.output.dir, get_h5Files_under_basepath = TRUE, remake_all = !append)
  rm(list = setdiff(ls(), keep_items))
}

# burst parameter to AUC
cat("\n- Check over component values by DIV, calculate AUC:\n")
resp <- check_existing(path = file.path(main.output.dir, "output"), pattern = "_AUC", pause_between_steps)
if (resp %in% c("r","a")) {
  # append has no effect here. This fun is relatively fast, so will remake all regardless
  source('DIV-interpolation-functions.R')
  source('estimate_missing_DIV.R')
  source('burst_parameter_to_AUC.R')
  rm(list = setdiff(ls(), keep_items))
}

# cytotox prep
cat("\n- Extract the cytotoxicity data from Calculations files:\n")
resp <- check_existing(path = file.path(main.output.dir, "output"), pattern = "_cytotox", pause_between_steps)
if (resp %in% c("r","a")) {
  source('cytotox_prep06.R')
  run_cytotox_functions(basepath = main.output.dir, get_files_from_log = TRUE, filename = paste0(project_name,"_cytotox.csv"), 
                        append = append)
  rm(list = setdiff(ls(), keep_items))
}

cat("\n'source_steps.R' is complete.\n")

