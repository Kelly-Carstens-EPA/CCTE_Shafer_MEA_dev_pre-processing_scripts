# spike_list_functions.R
# Diana Hall
# May 5, 2015
# purpose: functions to run the pipeline with spike list files
# one electrode and one timestampe per row



spkList2list <-function (file) {
    
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
      
      # if the total time from first spike to last spike is more than 3 minutes short of 900.00 seconds, flag it
      last_time <- tail(data.raw2$timestamps, n=1)
      first_time <- head(data.raw2$timestamps, n=1)
      if (last_time - first_time < (900.00 - 3*60)) {
        cat(paste0("\n",file," only goes from ",first_time," to ",last_time," seconds\n"))
        cat("Continue with this spike list file anyways? (Only do this if you know why the recording is significantly less than 900sec\n")
        resp <- readline(prompt = "(y/n): ")
        if (!(resp %in% c("y","Y"))) {
          stop("Update spike list file selection")
        }
      }
      
      # remove any points more than 900.00 seconds after the first recorded spike
      data.raw2 <- data.raw2[data.raw2$timestamps < first_time + 900.00,]
      
      # order data frame by electrode
      data.raw2<-data.raw2[order(data.raw2$elect), ]
      
      spikes<-split(data.raw2$timestamps, data.raw2$elect)
      rm(data.raw2)
    } else {
      spikes<-NULL
    }
    spikes
  }



axion.spkList.to.h5<-function(key, spkListFile, chem.info, debug=T){
  #function to convert spike list to h5 file
  
  #remove '_spike_list' from file name
  key <- sub("_spike_list.*$","",key)
  
  # # correct for issues with ( ) in file names
  # h5file <- gsub("\\(|\\)", "_", sprintf("%s/%s.h5", h5.dir, key))
  # if ( substring(basename(h5file), nchar(basename(h5file))-3, nchar(basename(h5file))-3)=="_" ){
  #   h5file<-paste( dirname(h5file) ,"/" ,
  #     unlist( strsplit(basename(h5file), split="_.h5") ) ,
  #                 ".h5", sep="")
  # }
  
  # what if we just allow ()? let's give it a try
  h5file <- sprintf("%s/%s.h5", h5.dir, key)

  #f is a list of all files
  f <- spkListFile

  #get spikes
  spikes.sep <- lapply(f, spkList2list)
  
  short.filenames <- gsub("_spike_list.csv", "", basename(f))
  summary.table <- t(sapply(spikes.sep, sjemea::axion.spikesum2) )
  rownames(summary.table) <- short.filenames
  ma <- do.call("rbind", lapply(spikes.sep, sjemea::axion.spikestodf))
  #s2 is a list with all the channels and spikes under each channel
  s2 <- split(ma$time, ma$elec)
  numelec <- length(s2)
  total.spikes <- sum(sapply(s2, length))
  time.ranges <- sapply(s2, range)
  time.min <- min(time.ranges[1, ])
  time.max <- max(time.ranges[2, ])
  #printf formats the text and variables to output ready formats
  #cat contatenates files and then prints them
  cat(printf("Total number of spikes: %d\n", total.spikes))
  cat(printf("Unique number of electrodes: %d\n", numelec))
  cat(printf("Time range [%.3f %.3f] (seconds)\n", time.min, 
             time.max))
  print(summary.table)
  
  map.to.h5.dh(s2, chem.info, h5file)
  
  if (debug) {
    d <- as.data.frame(summary.table)
    d2 <- data.frame(file = rownames(summary.table), d, 
                     stringsAsFactors = FALSE)
    h5write(d2, path.expand(h5file), "summary.table")
  }
  
  
  h5file
  
  
}










