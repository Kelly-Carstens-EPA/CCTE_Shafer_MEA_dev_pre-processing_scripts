# This script was copied from NFA Spike List to mc0 R Scripts Bit Bucket repository, branch MI-related-updates on 04/17/2020 9:57am
# This script prepares the data, then calls the nmi function to calcuation the mutual information

require('pracma')
library('compiler')

nmi_wrapper<-function(mutdata){
  
  apply_fun<-function(c){
    mutinfo<-sapply(c,nmi)
  }
  
  L<-length(mutdata)
  
  # mutdata is called 'x' in spikeLoadRoutines
  # the first half of x contains the recording data
  # the second half of x contains the meta data for each recording
  
  # collect the meta data
  meta <- mutdata[(L/2 + 1) : L]
  mutdata2 <- mutdata[1 : (L/2)]
  rm(mutdata)
  
  #### I take the transpose here so that mutual info values match up with the correct wells (I think, need to have Ken double check this)
  # Amy, March 2020: The previous comment was added by someone else. We have used the code as it is, so I believe it is correct
  mutdata3<-list()
  for (i in 1:length(mutdata2)){       
    mutdata3[[i]]<-t(mutdata2[[i]])
  }
  rm(mutdata2)
  
  #mutdata3<-sapply(mutdata2,t) doesn't do what I would expect, which is why I used the for loop above instead 
  
  testing<-sapply(mutdata3,apply_fun)
  
  # transform the data in "testing" into a table with the following column names
  heading<-c("date","Plate.SN","DIV","well","trt","dose","units","file.name", "Mutual Information")
  output <- list()
  for (i in 1:ncol(testing)) {
    
    reci <- data.frame(meta[[i]],testing[,i])
    colnames(reci) <- heading
    output <- rbind(output, reci)
    
  }
  
  return(output)
}