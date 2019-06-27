
setwd("/home/kluser2/datasets/ambracia_sims/all_free_small2/")


system('grep -vwE "NA" FOR_ABC > FOR_ABC_CLEAN')
library(abc)
library(Metrics)


a <- read.table("FOR_ABC_CLEAN", h=F)
dim(a)
LOGIT_MATRIX=matrix(nrow=13,ncol=2)
for (Z in 1:13){

LOGIT_MATRIX[Z,1]=max(a[,Z+1])
LOGIT_MATRIX[Z,2]=min(a[,Z+1])
    
}


predicted=vector()
actual=vector()
mean_diff=vector()


for (j in 1:dim(a)[1]){
  
  



# check dims

params <- a[-j,2:14]   #leave one out
stats <- a[-j,-(1:14)] # << ,<<
test <- a[j,-(1:14)]
test_params <-a[j,2:14]

dim(params)
dim(stats) # check dims to make sure






  
myabc <- abc(target=test, param=params, sumstat=stats, tol=0.1, method="ridge", hcorr=TRUE)

sum=summary(myabc)


for (w in 1:13){

PredictionsFile <- paste ("PREDICTIONS_MEAN_",w, sep = "", collapse = NULL)

print(j)
predicted[j]=as.numeric(sum[4,w])
actual[j]=as.numeric(test_params[w])

cat(paste(actual[j],'\t'),file=PredictionsFile,append=TRUE,sep='\t')
cat(paste(predicted[j],'\t'),file=PredictionsFile,append=TRUE,sep='\t')
cat('',file=PredictionsFile,append=TRUE,sep='\n')

}

}

print(mae(actual,predicted))


