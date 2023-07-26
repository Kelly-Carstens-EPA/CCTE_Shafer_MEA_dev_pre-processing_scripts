#' Create burst ontogeny data
#' 
#' Function to calculate 16 parameter values from MEA recordings. 
#'
#' @param h5Files Input files (spike list files converted to h5files, as created by the function h5_conversion())
#' @param csv.filename.AEfilt Name of output file with parameter values calculated based on active-electrodes filter only. Only created if AEfile = T. (see AEfile below)
#' @param csv.filename.ABEfilt Name of output file with parameter values calculated based on active-elecrodes and actively bursting electrodes filter for appropriate parameters
#' @param save.rdata Binary whether to save the object that contains the parameter values and additional information as an RData file
#' @param add.silent.wells Binary of whether to add data rows to the output table for the wells that did not have any activity. For these wells, 
#' parameters values in the list 'may.be.zero' are set to 0 and parameter values in the list 'may.not.be.zero' are set to NA.
#' @param AEfile Binary whether to save a file for which only an "active" electrodes filter has been applied to the appropriate parameters.
#' Note that this function will always save a file in which the appropriate parameters are filtered to "active" electrodes AND the following 4 parameters are filtered
#' to "actively bursting" electrodes: mean.isis.abe, per.spikes.in.burst.abe, mean.dur.abe, and mean.IBIs.abe. 
#' @param remake_all Binary whether to calculate the parameter values for all input h5files (TRUE), or only for the h5files that do not already have a corresponding .csv file in the 
#'
#' @details 
#' This function was originally created by Diana Hall as part of the meadq package. The function here has been edited to reference 
#' local.corr.all.ont.ae.filter (a modified version of a script not in meadq package)
#' for the calculation of the correlation coefficient 12/13/2019.
#'
#' @return
#' @export
#'
#' @examples
create_burst_ont_Data <-
  function(h5Files,  csv.filename.AEfilt, csv.filename.ABEfilt, save.rdata=F, add.silent.wells=T , AEfile = F, remake_all = TRUE){
    write.header <- T
    
    # get a list of the date-plate combinations selected
    date_plates <- unique( lapply(h5Files, function(x) strsplit(basename(x), split = "_")[[1]][2:3]) )
    num.plates = length( date_plates )
    
    for (cur.plate in 1:num.plates){
      
      write.header = T
      #####  Read in and make data
      s=list()
      date <- date_plates[[cur.plate]][1]
      plate <- date_plates[[cur.plate]][2]
      h5Files.wanted = h5Files[grep(x=h5Files,pattern= paste0(date,"_",plate,"_"))]
      
      #order according to DIV
      temp.l<-sapply(sapply(strsplit(basename(h5Files.wanted),split="_"), 
                            function(x) x[4]) , nchar)
      temp.div= substring( 
        sapply( strsplit(basename(h5Files.wanted),split="[_\\(\\)]"), function(x) x[4]),
        rep(4,length(h5Files.wanted)),temp.l)
      
      h5Files.wanted = h5Files.wanted[order(as.numeric(temp.div)) ]
      
      # if this data has already been processed and  remake_all is not TRUE, then skip this file
      check_file <- paste( paste( csv.filename.ABEfilt, 
                    strsplit(unique(basename(h5Files.wanted)),split="_")[[1]][2] ,
                    strsplit(unique(basename(h5Files.wanted)),split="_")[[1]][3],sep="_"), ".csv", sep="" ) 
      if(!remake_all && file.exists(check_file)) {
        next
      }
      
      for (cur.file in 1:length( h5Files.wanted ) ){
        
        s[[cur.file]]=h5.read.spikes( h5Files.wanted[cur.file] )
        # get meta-data
        g=c(); #remove any information stored in g
        g<-h5read(path.expand( h5Files.wanted[cur.file] ), name = "/")
        if ( is.element( "treatment", tolower( names(g) ) ) ){
          names(g)<-tolower(names(g))
          names(g)[names(g)=="scount"]<-"sCount"
          s[[cur.file]]$cw = substring(s[[cur.file]]$channel,1,2 )
          s[[cur.file]]$treatment = g$treatment         
          names(s[[cur.file]]$treatment) = g$well
          s[[cur.file]]$dose = g$dose
          names(s[[cur.file]]$dose) = g$well
          s[[cur.file]]$sCount = g$sCount
          names(s[[cur.file]]$sCount) = g$names
          s[[cur.file]]$units = g$units
          names(s[[cur.file]]$units) = g$well # appears to list all the wells, regardless of whether or not they spiked
          s[[cur.file]]$well = g$well 
          s[[cur.file]]$summary.table = g$summary.table
          s[[cur.file]]$array = g$array     
          s[[cur.file]]$DIV = strsplit( basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4]
        } else{
          print ('Error: meta-data not available in .h5 file');
          
        }
        
        # calculating parameters related to bursting
        allb=list()
        allb <- lapply(s[[cur.file]]$spikes, mi.find.bursts)
        s[[cur.file]]$allb <- allb       
        s[[cur.file]]$bs<-calc.burst.summary( s[[cur.file]]  )
        
        h=calc.all.isi(s[[cur.file]], allb)
        s[[cur.file]]$isi = h
        
        
        nspikes.old <- calculate.network.spikes(s[[cur.file]])
        
        s[[cur.file]]$ns.all <- nspikes.old$ns.all
        
        # this script get the parameter values of nsN and ns.T from each well individually, so those values could be changed individually
        nspikes <- summary.network.spikes.dh(s[[cur.file]],nspikes.old, ns.E = 2, sur)
        
        s[[cur.file]]$ns.all <- nspikes$ns.all
        
        # local correlation for one electrode with all other electrodes
        #  need to know active electrodes to filter calculation of r
        ae.index.l<- split((s[[cur.file]]$meanfiringrate*60>=5), f=s[[cur.file]]$cw )
        ae.index.v<- unlist( split((s[[cur.file]]$meanfiringrate*60>=5), f=s[[cur.file]]$cw ) )
        
        local.cor=c()
        local.cor<-local.corr.all.ont.ae.filter(s, t.p=cur.file, ae.index.v)
        # local.cor<-local.corr.all.ont(s, t.p=cur.file) # the original r calculation
        s[[cur.file]]$local.cor<-local.cor
        
        # # cv.time and network - these endpoints are not correct
        ## so we are no longer including these at this time
        # temp.list=c()
        # temp.list<- cv.timestamps(s[[cur.file]]$spikes, s, cur.file)
        # s[[cur.file]]$cv.time = temp.list[[1]]
        # names(s[[cur.file]]$cv.time) = unique(s[[cur.file]]$cw)
        # s[[cur.file]]$cv.network = temp.list[[2]]
        # names(s[[cur.file]]$cv.network) = unique(s[[cur.file]]$cw)
        
        
        
        
        ####  AE Filtered data frame
        ##+++++++++++++++++  make data frame        
        ae.index.l<- split((s[[cur.file]]$meanfiringrate*60>=5), f=s[[cur.file]]$cw )
        ae.index.v<- unlist( split((s[[cur.file]]$meanfiringrate*60>=5), f=s[[cur.file]]$cw ) )
        # we need the wells that have at least one AE
        well.names<-unique( subset(s[[cur.file]]$cw, ae.index.v) )
        num.wells = length(unique( subset(s[[cur.file]]$cw, ae.index.v) ) )
        well.indices = which(is.element(s[[cur.file]]$well, well.names ))
        
        if( grepl( strsplit(basename(s[[cur.file]]$file), split="_")[[1]][4],
                   pattern="DIV" ) ){
          DIV <- rep( substring(strsplit(basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4],
                                4, nchar(strsplit(basename(s[[cur.file]]$file), 
                                                  split="[_\\(\\)]")[[1]][4])  ), num.wells)           
        } else {
          DIV<-rep(strsplit(basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4],
                   num.wells)
        } 
        
        df=c() #erase prior data
        if ( num.wells>0 ){
          
          df<-data.frame(
            date = rep( unlist(strsplit(basename(s[[cur.file]]$file), split="_"))[2]  , 
                        num.wells) , 
            Plate.SN = rep( strsplit(basename(s[[cur.file]]$file), split="_")[[1]][3], 
                            num.wells) ,
            DIV = as.numeric( DIV) ,
            well = well.names ,
            trt = s[[cur.file]]$treatment[ well.names ]  , 
            dose = s[[cur.file]]$dose[ well.names ]  ,
            units = s[[cur.file]]$units[ well.names ]  ,
            meanfiringrate=sapply(by(subset(s[[cur.file]]$meanfiringrate, ae.index.v),
                                     subset(s[[cur.file]]$cw, ae.index.v), mean, na.rm=T), 
                                  function(x) x[1]) ,
            
            burst.per.min = sapply( 
              by( subset(s[[cur.file]]$bs$bursts.per.min, ae.index.v), 
                  subset(s[[cur.file]]$cw, ae.index.v), mean, na.rm=T),
              function(x) x[1]) ,
            
            mean.isis = sapply(by(subset( s[[cur.file]]$bs$mean.isi, ae.index.v),
                                  subset(s[[cur.file]]$cw, ae.index.v) , mean, na.rm=T), 
                               function(x) x[1]) ,
            
            per.spikes.in.burst =sapply(by(subset(s[[cur.file]]$bs$per.spikes.in.burst, 
                                                  ae.index.v),
                                           subset(s[[cur.file]]$cw, ae.index.v) , 
                                           mean, na.rm=T), function(x) x[1]) ,
            
            mean.dur =sapply(by(subset( s[[cur.file]]$bs$mean.dur, ae.index.v),
                                subset( s[[cur.file]]$cw, ae.index.v) , mean, na.rm=T), 
                             function(x) x[1]) ,
            
            mean.IBIs =sapply(by(subset( s[[cur.file]]$bs$mean.IBIs, ae.index.v) ,
                                 subset( s[[cur.file]]$cw, ae.index.v), mean, na.rm=T), 
                              function(x) x[1]) ,
            
            nAE = unlist( lapply( 
              by( s[[cur.file]]$meanfiringrate*60, s[[cur.file]]$cw, 
                  function(x) x>=5), sum ) )[well.names] ,
            
            nABE= unlist( lapply( 
              by( s[[cur.file]]$bs$bursts.per.min, s[[cur.file]]$cw, 
                  function(x) x>=0.5), sum ) )[well.names] ,
            
            ns.n = sapply(s[[cur.file]]$ns.all, 
                          function(x) x$brief['n'])[paste(well.names, "n", sep=".")] , 
            ns.peak.m = sapply(s[[cur.file]]$ns.all , 
                               function(x) x$brief['peak.m'])[ 
                                 paste(well.names, "peak.m", sep=".") ] ,
            ns.durn.m = sapply(s[[cur.file]]$ns.all, 
                               function(x) x$brief['durn.m'])[ 
                                 paste(well.names, "durn.m", sep=".") ] ,
            ns.percent.of.spikes.in.ns = sapply(s[[cur.file]]$ns.all, 
                                                function(x) x$brief['percent.of.spikes.in.ns'])[ paste(well.names, "percent.of.spikes.in.ns", sep=".") ] ,
            ns.mean.insis = sapply(s[[cur.file]]$ns.all, 
                                   function(x) x$brief['mean.insis'])[ paste(well.names, "mean.insis", sep=".") ] ,
            ns.durn.sd = sapply(s[[cur.file]]$ns.all, 
                                function(x) x$brief['durn.sd'])[ paste(well.names, "durn.sd", sep=".") ] ,
            ns.mean.spikes.in.ns = sapply(s[[cur.file]]$ns.all, 
                                          function(x) x$brief['mean.spikes.in.ns'])[ paste(well.names, "mean.spikes.in.ns", sep=".") ] ,
            
            r=sapply(s[[cur.file]]$local.cor, mean, na.rm=T)[well.names] ,
            
            file.name = rep( basename(s[[cur.file]]$file), num.wells )  
          ) #end data frame
          
          # write data to .csv file, 
          # wait to write for burst analysis
          
          if (AEfile) {
            write.table(  df, file = paste( paste( csv.filename.AEfilt,
                                                   strsplit(basename(s[[cur.file]]$file),split="_")[[1]][2],
                                                   strsplit(basename(s[[cur.file]]$file),split="_")[[1]][3],sep="_"), ".csv", sep="" ),
                          sep=",", append = !write.header, col.names=write.header, row.names=F )
          }
          
          
          #end of if FALSE
        } # end of if num.wells>0
        
        
        
        
        
        
        
        #####  ABE Filtered Data
        
        # +++++++++++++make data frame
        ae.index.l<- split((s[[cur.file]]$meanfiringrate*60>=5), f=s[[cur.file]]$cw )
        ae.index.v<- unlist( split((s[[cur.file]]$meanfiringrate*60>=5), f=s[[cur.file]]$cw ) )
        
        abe.index<- unlist( split((s[[cur.file]]$bs$bursts.per.min>=0.5), f=s[[cur.file]]$cw ) )
        # Restrict all ABE to also be AEs (there actually was an instance when this was not true!)
        abe.index = abe.index & ae.index.v
        
        # we need the wells that we have at least one AE
        well.names<-unique( subset(s[[cur.file]]$cw, ae.index.v) )
        num.wells = length(unique( subset(s[[cur.file]]$cw, ae.index.v) ) )
        well.indices = which(is.element(s[[cur.file]]$well, well.names ))
        
        
        if ( sum(abe.index)>0 ){
          
          # we need the wells that have at least one ABE
          well.names.abe<-unique( subset(s[[cur.file]]$cw, abe.index) )
          num.wells.abe = length(unique( subset(s[[cur.file]]$cw, abe.index) ) )
          well.indices.abe = which(is.element(s[[cur.file]]$well, well.names.abe )) 
          
          #make data frame
          if( grepl( strsplit(basename(s[[cur.file]]$file), split="_")[[1]][4],
                     pattern="[DIV]" ) ){
            DIV <- rep( substring(strsplit(basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4],
                                  4, nchar(strsplit(basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4])  ), 
                        num.wells)           
          } else {
            DIV<-rep(unlist( strsplit(basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4] ),
                     num.wells)
          } 
          
          
          mean.isis.abe = sapply(by(subset( s[[cur.file]]$bs$mean.isi, abe.index),
                                    subset(s[[cur.file]]$cw, abe.index) , mean, na.rm=T), 
                                 function(x) x[1]) 
          
          per.spikes.in.burst.abe =sapply(by(subset(s[[cur.file]]$bs$per.spikes.in.burst, abe.index),
                                             subset(s[[cur.file]]$cw, abe.index) , mean, na.rm=T),
                                          function(x) x[1]) 
          
          mean.dur.abe =sapply(by(subset( s[[cur.file]]$bs$mean.dur, abe.index),
                                  subset( s[[cur.file]]$cw, abe.index) , mean, na.rm=T), 
                               function(x) x[1]) 
          
          mean.IBIs.abe =sapply(by(subset( s[[cur.file]]$bs$mean.IBIs, abe.index) ,
                                   subset( s[[cur.file]]$cw, abe.index), mean, na.rm=T), 
                                function(x) x[1]) 
          
          if ( num.wells>num.wells.abe ){
            
            mean.isis<-rep(NA, length(well.names))
            names(mean.isis)<-well.names
            mean.isis[well.names.abe]<-mean.isis.abe
            
            per.spikes.in.burst<-rep(NA, length(well.names))
            names( per.spikes.in.burst )<-well.names
            per.spikes.in.burst[well.names.abe]<-per.spikes.in.burst.abe
            
            mean.dur<-rep(NA, length(well.names))
            names(mean.dur)<-well.names
            mean.dur[well.names.abe]<-mean.dur.abe
            
            mean.IBIs<-rep(NA, length(well.names))
            names(mean.IBIs)<-well.names
            mean.IBIs[well.names.abe]<-mean.IBIs.abe
          } else {
            mean.isis <- mean.isis.abe  
            
            per.spikes.in.burst <-per.spikes.in.burst.abe 
            
            mean.dur<-mean.dur.abe 
            
            mean.IBIs<-mean.IBIs.abe  
            
          }
          
          
          df2=c()
          df2<-data.frame(
            date = df$date ,
            Plate.SN = df$Plate.SN ,
            DIV = as.numeric( df$DIV ) ,
            well = df$well ,
            trt = df$trt  ,
            dose = df$dose  ,
            units = df$units  ,
            meanfiringrate = df$meanfiringrate ,
            
            burst.per.min = df$burst.per.min ,
            
            mean.isis = mean.isis ,
            
            per.spikes.in.burst = per.spikes.in.burst  ,
            
            mean.dur = mean.dur ,
            
            mean.IBIs = mean.IBIs  ,
            
            nAE = df$nAE ,
            
            nABE = df$nABE ,
            
            ns.n = df$ns.n , #nspikes$ns.all$A8$brief
            ns.peak.m = df$ns.peak.m ,
            ns.durn.m = df$ns.durn.m   ,
            ns.percent.of.spikes.in.ns = df$ns.percent.of.spikes.in.ns ,
            ns.mean.insis = df$ns.mean.insis ,
            ns.durn.sd = df$ns.durn.sd ,
            ns.mean.spikes.in.ns = df$ns.mean.spikes.in.ns ,
            
            
            r= df$r ,
            
            file.name = df$file.name
          )
          
          may.be.zero<-c("nAE", "nABE", "meanfiringrate", "burst.per.min",
                         "per.spikes.in.burst", 
                         "ns.percent.of.spikes.in.ns","ns.n", 
                         "per.spikes.in.burst", "r")
          
          # add in well level number for wells that don't appear in calculation
          if (add.silent.wells){ 
            silent.wells<-setdiff( s[[cur.file]]$well, df2$well )
            
            # get ID info from filename, not df2 
            # Because df2 might not have been created if num.wells==0
            date_value <- strsplit(basename(s[[cur.file]]$file), split="_")[[1]][2]
            Plate.SN_value <- strsplit(basename(s[[cur.file]]$file), split="_")[[1]][3]
            DIV_value <- as.numeric(strsplit(basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4])
            
            date.sw<-rep(date_value, length(silent.wells) )
            Plate.SN.sw<-rep(Plate.SN_value, length(silent.wells) )
            DIV.sw<-rep(DIV_value, length(silent.wells) )
            
            trt.sw<- s[[1]]$treatment[silent.wells]
            dose.sw<-s[[1]]$dose[silent.wells]
            units.sw<-s[[1]]$units[silent.wells]
            file.name.sw<-rep(basename(s[[cur.file]]$file), length(silent.wells) )
            
            for (i in 1:length(may.be.zero)){
              assign(paste(may.be.zero[i],"sw",sep="."), 
                     rep(0, length(silent.wells) ) ) 
            }
            
            # don't have names(df2), if df2 has not been created
            may.not.be.zero <- c("mean.isis", "mean.dur", "mean.IBIs",
                                 "ns.peak.m", "ns.durn.m", "ns.mean.insis",
                                 "ns.durn.sd", "ns.mean.spikes.in.ns")
            for (i in 1:length(may.not.be.zero) ){
              assign(paste(may.not.be.zero[i],"sw",sep="."), 
                     rep(NA, length(silent.wells) ) ) 
            }
            
            df2.sw<-
              data.frame(date=date.sw, Plate.SN=Plate.SN.sw, 
                         DIV=DIV.sw, well= silent.wells, trt = trt.sw, 
                         dose = dose.sw,
                         units= units.sw, meanfiringrate = meanfiringrate.sw, 
                         burst.per.min=burst.per.min.sw, mean.isis=mean.isis.sw,
                         per.spikes.in.burst=per.spikes.in.burst.sw, 
                         mean.dur=mean.dur.sw, mean.IBIs=mean.IBIs.sw, nAE=nAE.sw,
                         nABE=nABE.sw, ns.n=ns.n.sw, ns.peak.m=ns.peak.m.sw, 
                         ns.durn.m=ns.durn.m.sw, 
                         ns.percent.of.spikes.in.ns=ns.percent.of.spikes.in.ns.sw,
                         ns.mean.insis=ns.mean.insis.sw, 
                         ns.durn.sd=ns.durn.sd.sw, 
                         ns.mean.spikes.in.ns=ns.mean.spikes.in.ns.sw, 
                         r=r.sw, 
                         file.name=file.name.sw )
            
            df2<-rbind(df2, df2.sw)
            
          } # end of if add silent wells
          
          #change NA to 0 for endpoints that may be zero
          for(cur.var in may.be.zero ){
            df2[ is.na(df2[,cur.var]) , cur.var]<-0
          }
          
          
          # write data to .csv file
          write.table(  df2, 
                        file= paste0( paste( csv.filename.ABEfilt, 
                                             strsplit(basename(s[[cur.file]]$file),split="_")[[1]][2],
                                             strsplit(basename(s[[cur.file]]$file),split="_")[[1]][3],sep="_"), ".csv") ,
                        sep=",", append = !write.header, col.names=write.header, row.names=F )
          
          
          
        } else { 
          # (sum(abe.index) <= 0)
          df2<-df
          df2$mean.isis<-rep(NA, length(df2$well))
          df2$mean.IBIs<-rep(NA, length(df2$well))
          df2$per.spikes.in.burst<-rep(NA, length(df2$well))
          df2$mean.dur<-rep(NA, length(df2$well))
          
          
          
          
          may.be.zero<-c("nAE", "nABE", "meanfiringrate", "burst.per.min",
                         "per.spikes.in.burst", 
                         "ns.percent.of.spikes.in.ns","ns.n", 
                         "per.spikes.in.burst", "r")
          
          # add in well level number for wells that don't appear in calculation
          if (add.silent.wells){           
            silent.wells<-setdiff( s[[cur.file]]$well, df2$well )
            
            # get ID info from filename, not df2 
            # Because df2 might not have been created if num.wells==0
            date_value <- strsplit(basename(s[[cur.file]]$file), split="_")[[1]][2]
            Plate.SN_value <- strsplit(basename(s[[cur.file]]$file), split="_")[[1]][3]
            DIV_value <- strsplit(basename(s[[cur.file]]$file), split="[_\\(\\)]")[[1]][4]
            
            date.sw<-rep(date_value, length(silent.wells) )
            Plate.SN.sw<-rep(Plate.SN_value, length(silent.wells) )
            DIV.sw<-rep(DIV_value, length(silent.wells) )
            
            trt.sw<- s[[1]]$treatment[silent.wells]
            dose.sw<-s[[1]]$dose[silent.wells]
            units.sw<-s[[1]]$units[silent.wells]
            file.name.sw<-rep(basename(s[[cur.file]]$file), length(silent.wells) )
            for (i in 1:length(may.be.zero)){
              assign(paste(may.be.zero[i],"sw",sep="."), 
                     rep(0, length(silent.wells) ) ) 
            }
            # don't have names(df2), if df2 has not been created
            may.not.be.zero <- c("mean.isis", "mean.dur", "mean.IBIs",
                                 "ns.peak.m", "ns.durn.m", "ns.mean.insis",
                                 "ns.durn.sd", "ns.mean.spikes.in.ns")
            
            for (i in 1:length(may.not.be.zero) ){
              assign(paste(may.not.be.zero[i],"sw",sep="."), 
                     rep(NA, length(silent.wells) ) ) 
            }
            
            df2.sw<-
              data.frame(date=date.sw, Plate.SN=Plate.SN.sw, 
                         DIV=DIV.sw, well= silent.wells, trt = trt.sw, 
                         dose = dose.sw,
                         units= units.sw, meanfiringrate = meanfiringrate.sw, 
                         burst.per.min=burst.per.min.sw, mean.isis=mean.isis.sw,
                         per.spikes.in.burst=per.spikes.in.burst.sw, 
                         mean.dur=mean.dur.sw, mean.IBIs=mean.IBIs.sw, nAE=nAE.sw,
                         nABE=nABE.sw, ns.n=ns.n.sw, ns.peak.m=ns.peak.m.sw, 
                         ns.durn.m=ns.durn.m.sw, 
                         ns.percent.of.spikes.in.ns=ns.percent.of.spikes.in.ns.sw,
                         ns.mean.insis=ns.mean.insis.sw, 
                         ns.durn.sd=ns.durn.sd.sw, 
                         ns.mean.spikes.in.ns=ns.mean.spikes.in.ns.sw, 
                         r=r.sw, 
                         file.name=file.name.sw )
            
            df2<-rbind(df2, df2.sw)
            
          } # end of if add silent wells
          
          #change NA to 0 for endpoints that may be zero
          for(cur.var in may.be.zero ){
            df2[ is.na(df2[,cur.var]) , cur.var]<-0
          }
          
          
          
          # write data to .csv file
          write.table(  df2, 
                        file= paste0( paste( csv.filename.ABEfilt, 
                                             strsplit(basename(s[[cur.file]]$file),split="_")[[1]][2],
                                             strsplit(basename(s[[cur.file]]$file),split="_")[[1]][3],sep="_"), ".csv") ,
                        sep=",", append = !write.header, col.names=write.header, row.names=F )
          
          
        }#end of if (sum(ABEind>0))
        
        
        
        write.header=F
      } #end loop through h5Files in current plate
      
      
      if (save.rdata ){
        save(s, file = paste(h5.dir ,  paste("/s",
                                             strsplit(basename(s[[cur.file]]$file),split="_")[[1]][2] ,
                                             strsplit(basename(s[[cur.file]]$file),split="_")[[1]][3] , 
                                             paste('thru',s[[cur.file]]$DIV, ".Rdata", sep="") ,
                                             sep="_" ), sep="")    ) 
      }
      
      print(paste('-----------------       Done with file # ', cur.file, sep= "") )
      
    } # end of loop through plates
    
    
  }
