plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))
class(density)


b=rnorm(100,100,2)
a=density(b,n=11)$y
a
myFiles <- list.files(pattern="COMUSTATS_*")            
myfile1=read.csv(file='COMUSTATS_6',sep='\t',header=TRUE, row.names=NULL)
head(myfile1)





inputFile <- "COMUSTATS_6"
con  <- file(inputFile, open = "r")

dataList <- list()
ecdfList <- list()

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  
j=strsplit(oneLine,'\t')
print(j)
  
} 

