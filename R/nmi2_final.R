require('pracma')
library('compiler')

nmi<-function(x){  
  
  #n=size(x,1)
  #N=size(x,2)
  x <- as.matrix(x)
  n <- dim(x)[1]
  N <- dim(x)[2]
  
  x[x!=0]<-1
  
  P= diag(x%*%t(x))/N;  
  
  H = -P*log2(P) - (1-P)*log2(1-P);
  
  entropyCutOff = -.999*log2(.999) - (1-.999)*log2(1-.999);
  
  x<-as.matrix(x[!is.na(H),])
  H<-H[!is.na(H)]
  x=x[H > entropyCutOff,]
  H=H[H > entropyCutOff]
  nn <- dim(x)[1]
  # nn=size(x,1)
  #if (nn<=1){
  #Info=0
  #return(Info)
#}
  if (is.null(nn)){           ### added this line because some null values for nn were stopping the function from running
    Info=0
    return(Info)
  }else if (nn<=1){
    Info=0
    return(Info)
  }
#X<-c()  
#for (i in 1:nrow(t(x))){  
# X[i]=strtoi(gsub(" ","",gsub(",","",toString(t(x)[i,]))),base=2)+1
#}
 bin2Dec <- function(m) {sum(2^(which(rev(m) == 1)-1))}
 X <- apply(x,2,bin2Dec)+1

 

 


 PP=accumarray(X,rep(1,length(X)))
 PP<-PP[PP!=0]
 PP<-PP/N
 HH = sum(-PP*log2(PP))


#H[is.nan(H)]<-0        
#HH[is.nan(HH)]<-0

#H[is.na(H)]<-0        
#HH[is.na(HH)]<-0      
 
 
 Info=sum(H)-HH
 Info=(1/(nn-1))*Info
 return(Info)
}

nmi <-cmpfun(nmi)