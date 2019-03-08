#plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))
#class(density)




setwd("C:/Users/John/Desktop/Ambracia Sims/")


system('grep -vwE "NA" FOR_ABC > FOR_ABC_CLEAN')
library(abc)
library(Metrics)


a <- read.table("FOR_ABC_CLEAN", h=F)
dim(a)
predicted=vector()
actual=vector()
mean_diff=vector()

#PLOTTING
#plot(params[,7],stats[,7])
#abline(lm(params[,7]~stats[,7]), col="red")
#lines(lowess(params[,7],stats[,7]), col="blue")
#scatter_plot <- ggplot(foot_height, aes(foot, height))
#scatter_plot + geom_point() + labs(x = "foot length (cm)", y = "height (cm)") + geom_smooth(method="lm")









for (j in 1:dim(a)[1]){
  
  



# check dims

params <- a[-j,2:14]   #leave one out
stats <- a[-j,-(1:14)] # << ,<<
test <- a[j,-(1:14)]
test_params <-a[j,2:14]

dim(params)
dim(stats) # check dims to make sure


#headers=read.table("LABELS_FOR_ABC", h=F,sep='\t')
#head(headers)#check headers

#names(params) <- headers[1,1:14]
#names(stats) <- mynames.stats[1,-(1:14)]



mysample <- 1:248

  
myabc <- abc(target=test[,mysample], param=params, sumstat=stats[,mysample], tol=0.1, method="ridge", hcorr=TRUE)

sum=summary(myabc)
predicted[j]=as.numeric(sum[4,4])
actual[j]=as.numeric(test_params[4])

}

print(mae(actual,predicted))

