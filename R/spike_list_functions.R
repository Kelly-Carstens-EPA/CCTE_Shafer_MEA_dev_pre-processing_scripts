# spike_list_functions.R
# Diana Hall
# May 5, 2015
# purpose: functions to run the pipeline with spike list files
# one electrode and one timestampe per row
# Adapted by Amy Carpenter


spkList2list <-function (file, recording_duration_sec = 900.00) {
  
  data.raw<-read.csv(file,header=F,colClasses=c("character", "NULL", "character","character","character")) # make electrode column char, not factor
  
  # remove the rows after the "Well Information" tag
  # adapted from "CG Additions 1/13/201"
  if (any(grepl("Well Information",data.raw$V1))) {
    well_information_row <- which(data.raw$V1 == "Well Information")
    data.raw <- data.raw[1:(well_information_row - 1),]
  }
  
  # remove first column and get the header
  # (using this instead of header=T, because sometimes the colnames are not in the first row in the spike list file)
  data.raw <- data.raw[, c(2,3,4)]
  header_row <- which(data.raw[,1] == "Time (s)")
  colnames(data.raw) <- data.raw[header_row,]
  data.raw <- data.raw[-c(header_row),]
  
  # remove any remaning empty rows, then convert to time and amplitude to numeric
  data.raw <- data.raw[data.raw$Electrode != "",] # works even if no rows in Electrode are == ""
  data.raw[,1]<-as.numeric(as.character(data.raw[,1]))
  data.raw[,3]<-as.numeric(as.character(data.raw[,3]))
  
  #remove NA
  ind.want<-which(!is.na(data.raw[,1]) )
  if (length( ind.want )>0){
    data.raw2<-data.frame(
      elect<-data.raw[ ind.want ,"Electrode"],
      timestamps<-data.raw[ ind.want ,"Time (s)"],
      stringsAsFactors = FALSE
    )
    rm(data.raw)
    data.raw2 <- data.raw2[order(data.raw2$timestamps),]
    
    # if the total time from first spike to last spike is more than 3 minutes short of recording_duration_sec seconds, flag it
    last_time <- tail(data.raw2$timestamps, n=1)
    first_time <- head(data.raw2$timestamps, n=1)
    if (last_time - first_time < (recording_duration_sec - 3*60)) {
      stop(paste0("\n",file," only goes from ",first_time," to ",last_time," seconds\n(over 3 minutes short of recording_duration_sec)"))
    }
    
    # remove any points more than recording_duration_sec seconds after the first recorded spike
    data.raw2 <- data.raw2[data.raw2$timestamps < first_time + recording_duration_sec,]
    
    # order data frame by electrode
    data.raw2<-data.raw2[order(data.raw2$elect), ]
    
    spikes<-split(data.raw2$timestamps, data.raw2$elect)
    rm(data.raw2)
  } else {
    spikes<-NULL
  }
  spikes
}


# function to convert spike list to h5 file
axion.spkList.to.h5<-function(spkListFile, plate.chem.info, h5.dir, recording_duration_sec = 900.00){
  
  # Create h5file name
  h5file_basename <- sub("_spike_list.csv$",".h5",basename(spkListFile))
  h5file <- file.path(h5.dir, h5file_basename)
  
  # Read the spike list file with the function spkList2list
  spikes.sep <- lapply(spkListFile, spkList2list,
                       recording_duration_sec = recording_duration_sec)
  
  # Display summary table of the spikes (optional)
  summary.table <- t(sapply(spikes.sep, sjemea::axion.spikesum2) )
  rownames(summary.table) <- gsub("_spike_list.csv", "", basename(spkListFile))
  ma <- do.call("rbind", lapply(spikes.sep, sjemea::axion.spikestodf))
  #s2 is a list with all the channels and spikes under each channel
  s2 <- split(ma$time, ma$elec)
  numelec <- length(s2)
  total.spikes <- sum(sapply(s2, length))
  time.ranges <- sapply(s2, range)
  time.min <- min(time.ranges[1, ])
  time.max <- max(time.ranges[2, ])
  #printf formats the text and variables to output ready formats
  # cat(printf("Total number of spikes: %d\n", total.spikes))
  # cat(printf("Unique number of electrodes: %d\n", numelec))
  # cat(printf("Time range [%.3f %.3f] (seconds)\n", time.min, 
  #            time.max))
  print(summary.table)
  
  # Save the spikes information and meta data to the h5file
  map.to.h5.dh(s2, plate.chem.info, h5file)
  
  # Add the summary.table to the h5file
  d <- as.data.frame(summary.table)
  d2 <- data.frame(file = rownames(summary.table), d, 
                   stringsAsFactors = FALSE)
  h5write(d2, path.expand(h5file), "summary.table")

  h5file
  
}










