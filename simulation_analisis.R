plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))
class(density)


b=rnorm(100,100,2)
a=density(b,n=11)$y

for (i in 56:99){
############################################################
#PARAMETERS input  
  
ParametersFile <- paste ("PARAMETERS_",i, sep = "", collapse = NULL)
con  <- file(ComusFile, open = "r")

ParametersList <- list()
j=0
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) 
{
  
  ParametersList[j] <- strsplit(oneLine,'\t')
  j <- j + 1

  
    
  
} 
close(con)  
  
  
  
  
  
  
  
  
############################################################  
#COMUStats input
  
ComusFile <- paste ("COMUSTATS_",i, sep = "", collapse = NULL)
con  <- file(ComusFile, open = "r")

ComusList <- list()
j=0
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) 
  {
  
  ComusList[j] <- strsplit(oneLine,'\t')
  j <- j + 1
  
  
  } 
close(con)




##########################################################
#f3 input

F3file <- paste ("f3FINAL_",i,".txt", sep= "", collapse = NULL)
con  <- file(F3file, open = "r")

F3list <- list()
j=0
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) 
{
  
  F3list[j] <- strsplit(oneLine,'\t')
  j <- j + 1
  
  
} 
close(con)









#########################################################
#PCA clustering

PCAfile <- paste ("PCA_CLUSTERING_",i, sep = "", collapse = NULL)
con  <- file(PCAfile, open = "r")

PCAlist <- list()
j=0
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) 
{
  
  PCAlist[j] <- strsplit(oneLine,'\t')
  j <- j + 1
  
  
} 
close(con)






}

