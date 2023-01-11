# This script was copied from NFA Spike List to mc0 R Scripts Bit Bucket repository, branch MI-related-updates on 04/17/2020 9:57am
# this script was edited in order to get the plate and data from the file name, instead of the h5file, as is done in the current master branch

# This script assembles the spike data from the h5 files for the Mutual Information
# This script is configured to extract the DIV values from the 4th tag in the h5 filename
# div<-strsplit(file.name,split = "[_\\(\\)]")[[1]][4] (see below)
# Change that line if needed

require('rhdf5')
require('gtools')

load.spikedata_final <- function(pseudoSamplingRate,fileName){
  # Loads spike data for a particular recording from h5 file.
  # Designed for [6,8] plate array with 4x4 MEAs in each well.
  # Returns a [6,8] data.frame of sparse [16,N] matrices. 
  
  x<-list()
  meta<-list()
  for (j in 1:length(fileName)){
    
    test_h5<-try(h5read(fileName[[j]], '/Well'), silent=T) # checking capitalization/version of objects in fileName[[j]]
    if (inherits(test_h5,"try-error")){
      
      spikes <- h5read(fileName[[j]], '/spikes')
      well <- h5read(fileName[[j]], '/well')
      names <- h5read(fileName[[j]], '/names')
      sCount <- h5read(fileName[[j]], '/sCount')
      treatment<-h5read(fileName[[j]], '/treatment')
      dose<-h5read(fileName[[j]], '/dose')
      units<-h5read(fileName[[j]], '/units')
      # added the following 2 lines Mar 26, 2020 to work with changes in chem.info.3 in h5 file - 04/17/2020: removed for this script
      # plate<-h5read(fileName[[j]], '/Plate.SN')
      # date<-h5read(fileName[[j]], '/Experiment.Date')
      # will get DIV from filename
      #div<-h5read(fileName[[j]], '/DIV')
    }  else{
      spikes <- h5read(fileName[[j]], '/spikes')
      well <- h5read(fileName[[j]], '/Well')
      names <- h5read(fileName[[j]], '/names')
      sCount <- h5read(fileName[[j]], '/sCount')
      treatment<-h5read(fileName[[j]], '/Treatment')
      dose<-h5read(fileName[[j]], '/Dose')
      units<-h5read(fileName[[j]], '/Units')
      plate<-h5read(fileName[[j]], '/Plate.SN')
      date<-h5read(fileName[[j]], '/Experiment.Date')
    }
    file.name<-basename(fileName[[j]])
    
    t <- seq(min(spikes[]),max(spikes[]),by=1/pseudoSamplingRate)
    t0 <- min(spikes[])
    
    spikeStreamArray <- as.list(matrix(nrow = 6, ncol = 8))
    dim(spikeStreamArray) <- c(6,8)
    # Initialize spiking array of zeros for each well
    for (rr in 1:6) {
      for (cc in 1:8) {
        spikeStreamArray[[rr,cc]] <- matrix(0,nrow = 16, ncol = (length(t)))
      }
    }
    
    # Populate each electrode in each well.
    # Electrode conversion chart:
    # #| 1  2  3  4
    # --------------
    # 1| 1  2  3  4
    # 2| 5  6  7  8
    # 3| 9  10 11 12
    # 4| 13 14 15 16
    
    recallIndex = 1
    
    for (ii in 1:length(sCount[])) {
      rr <- substr(names[ii],1,1) # the spiking electrode names
      # Well row number:
      # Convert character to ascii code value.. then normalize:
      # [A,B,C,D,E,F] -> [65,66,67,68,69,70] -> '' - 64 ->[1,2,3,4,5,6]
      rr <- as.integer(asc(rr))-64L
      # Well column number:
      cc <- strtoi(substr(names[ii],2,2))
      # Electrode Number (see conversion chart above)
      ee <- strtoi(substr(names[ii],4,4))-1
      ee <- ee*4 + strtoi(substr(names[ii],5,5))
      
      # Now that we know where we are, we assign spikes at each time
      # that they occur
      #
      # Retrieve the spiking times:
      spikeTimes <- spikes[recallIndex:(recallIndex+sCount[ii]-1)]
      # Convert to indices:
      spikeTimes <- floor((spikeTimes-t0)*pseudoSamplingRate)+1 # convert time in s to num of 'samples'/bins after start
      # Assign binary 1 values to array where spikes occur:
      spikeStreamArray[[rr,cc]][ee,spikeTimes] <- 1
      
      # Update the recallIndex value so we crawl ahead through spikes
      recallIndex <- recallIndex + sCount[ii]
    }
    x[[j]]<-spikeStreamArray
    
    wells_vector<-c()
    trt_vector<-c()
    dose_vector<-c()
    unit_vector<-c()
    plate_vector<-c()
    date_vector<-c()
    div_vector<-c()
    file.name_vector<-c()
    
    if (inherits(test_h5,"try-error")){
      
      
      for (i in 1:dim(well)){
        wells_vector[i]<-well[i]
        trt_vector[i]<-treatment[i]
        dose_vector[i]<-dose[i]
        unit_vector[i]<-units[i]
        # plate_vector[i]<-plate[i]
        # date_vector[i]<-date[i]
        #div_vector[i]<-div[i]
      }
      # add 04/17/2020 to collect date and plate from h5file name instead of from h5 body
      date<-strsplit(file.name,split = "_")[[1]][2]
      plate<-strsplit(file.name,split = "_")[[1]][3]
      div<-strsplit(file.name,split = "[_\\(\\)]")[[1]][4]
      
      div_vector<-rep(div,length(wells_vector))
      plate_vector<-rep(plate,length(wells_vector))
      date_vector<-rep(date,length(wells_vector))
      
    } else{
      
      for (i in 1:dim(well)){
        wells_vector[i]<-well[i]
        trt_vector[i]<-treatment[i]
        dose_vector[i]<-dose[i]
        unit_vector[i]<-units[i]
        plate_vector[i]<-plate[i]
        date_vector[i]<-date[i]
        #div_vector[i]<-div[i]
      }
      div<-substr(fileName[[j]],nchar(fileName[[j]])-11,nchar(fileName[[j]])-10)
      #plate<-substr(fileName[[j]],nchar(fileName[[j]])-22,nchar(fileName[[j]])-13)
      #date<-substr(fileName[[j]],nchar(fileName[[j]])-31,nchar(fileName[[j]])-23)
      
      div_vector<-rep(div,length(wells_vector))
    }
    file.name_vector<-rep(file.name, length(wells_vector))
    
    meta[[j]]<-cbind(date_vector,plate_vector,div_vector,wells_vector,trt_vector,dose_vector,unit_vector,file.name_vector)
    
  }
  
  # add the meta data to the end of x
  L <- length(x)
  for (i in 1:L) {
    x[[L+i]] <- meta[[i]]
  }
  
  # Return the collection of recordings across wells
  # for the given plate.
  return(x)
}