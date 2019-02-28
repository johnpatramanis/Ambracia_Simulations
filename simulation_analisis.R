plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))
class(density)


b=rnorm(100,100,2)
a=density(b,n=11)$y

setwd("C:/Users/John/Desktop/Ambracia Sims/parameters")

for (i in 56:60){

############################################################
#PARAMETERS input  

ParametersFile <- paste ("PARAMETERS_",i, sep = "", collapse = NULL)
con  <- file(ParametersFile, open = "r")

ParametersList <- list()
j=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  
  ParametersList[j] <- strsplit(line,'\t')
  j <- j + 1
}

close(con)  
print(ParametersList) 

  
  
  
  
  
  
############################################################  
#COMUStats input

ComusFile <- paste ("COMUSTATS_",i, sep = "", collapse = NULL)
con  <- file(ComusFile, open = "r")

ComusList <- list()
j=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  
  ComusList[j] <- strsplit(line,'\t')
  j <- j + 1
}

close(con)  

comusdataframe=data.frame(matrix(as.numeric(unlist(ComusList[-1])),nrow=length(ComusList)-1,byrow=T))

j=1
newcomusdistributions=list()
for (k in 1:ncol(comusdataframe)){
  

newcomusdistributions[j]=list(density(comusdataframe[,k],n=11)$x)
j=j+1
  
  
}










##########################################################
#f3 input

F3file <- paste ("f3FINAL_",i,".txt", sep= "", collapse = NULL)
con  <- file(F3file, open = "r")

F3list <- list()
j=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  
  F3list[j] <- strsplit(line,'\t')
  j <- j + 1
}

close(con) 

f3dataframe=data.frame(matrix(as.numeric(unlist(F3list)),nrow=length(F3list),byrow=T))
f3distribution=f3dataframe[,1]
newf3distribution=density(f3distribution,n=11)$x








#########################################################
#PCA clustering

PCAfile <- paste ("PCA_CLUSTERING_",i, sep = "", collapse = NULL)
con  <- file(PCAfile, open = "r")

PCAlist <- list()
j=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  
  PCAlist[j] <- strsplit(line,'\t')
  j <- j + 1
}

close(con) 


#write.table(t(as.data.frame(ParametersList[1])),file="FOR_ABC", quote=F,sep="\t",row.names=F,col.names=F)


cat(paste(ParametersList[1],ParametersList[2],collapse=''),file='FOR_ABC',append=TRUE,sep="\n")


}

