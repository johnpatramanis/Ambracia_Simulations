#plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))
#class(density)




setwd("C:/Users/John/Desktop/Ambracia Sims/parameters")


#cat(paste('who_col','N_locals','N_metropolis','N_initial_colony','N_final_colony','r_locals','r_metropolis','r_colony','m_locals_metropolis0','m_locals_colony','m_metropolis_locals','m_metropolis_colony','m_colony_locals','m_colony_metropolis',collapse=' \t'),file='FOR_ABC',append=TRUE,sep="\t")

#cat(paste(paste('f3',seq(1,11),sep='_'),collapse=' \t'),file='FOR_ABC',append=TRUE,sep="\t")




for (i in 0:2499){
print(i)
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

  
migrations = list(unlist(stringr::str_extract_all(ParametersList[[2]],'([0-9]+.[0-9]+)'),recursive = TRUE))
  
  
  
  
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
  
mycomusdataframe=na.omit(comusdataframe[,k])
  
newcomusdistributions[j]=list(density(mycomusdataframe,n=11)$x)
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

if (length(F3list)==0){
newf3distribution=list(rep('NA',11))
flag=1
   
}

if (length(F3list)>1){
 
f3dataframe=data.frame(matrix(as.numeric(unlist(F3list)),nrow=length(F3list),byrow=T))
f3distribution=f3dataframe[,1]
newf3distribution=density(f3distribution,n=11)$x
flag=0
}

if (length(F3list)==1){
  newf3distribution=list(rep('NA',11))
  flag=1
}









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
PCAlist=list(as.numeric(unlist(PCAlist)))
close(con) 



for (w in 1:length(ParametersList[[1]])){

cat(paste(ParametersList[[1]][w],'\t'),file='FOR_ABC',append=TRUE,sep='\t')
}
for (w in 1:length(migrations[[1]])){
  #for (c in 1:length(ParametersList[[2]][w])){
    
    cat(paste(migrations[[1]][w],'\t'),file='FOR_ABC',append=TRUE,sep='\t')
    
  #}
}
for (w in 1:length(PCAlist[[1]])){
  
  cat(paste(PCAlist[[1]][w],'\t'),file='FOR_ABC',append=TRUE,sep='\t')
}

if (flag==0){
for (w in 1:length(newf3distribution)){
  
  cat(paste(newf3distribution[w],'\t'),file='FOR_ABC',append=TRUE,sep='\t')
}
}
if (flag==1){
for (w in 1:length(newf3distribution)){
  
cat(paste(newf3distribution[[w]],'\t'),file='FOR_ABC',append=TRUE,sep='\t')  
  
}  
  
  
  
}





for (c in 1:length(newcomusdistributions)){
  
  for (w in 1:length(newcomusdistributions[[c]])){
    
  cat(paste(newcomusdistributions[[c]][w],'\t'),file='FOR_ABC',append=TRUE,sep='\t')  
  }
  
  
}



cat('',file='FOR_ABC',append=TRUE,sep='\n')


}

system('grep -vwE "NA" FOR_ABC > FOR_ABC_CLEAN')

a <- read.table("FOR_ABC_CLEAN", h=F)

dim(a)# check dims

params <- a[-1,1:14]   #leave first row out for testing, take the rest
stats <- a[-1,-(1:14)] # << ,<<


dim(stats) # check dims to make sure


headers=read.table("LABELS_FOR_ABC", h=F,sep='\t')
head(heades)#check headers

names(params) <- headers[1,1:14]
names(stats) <- mynames.stats[1,-(1:14)]

library(abc)



