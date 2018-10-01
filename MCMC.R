
# Normal likelihood and standardised variables (independent & dependent) Script1
source('MCMC Functions/Jags-Ymet-XmetMulti-Mrobust(Norm-Stand).R')

# Normal likelihood and non-standardised binary variables (dependent and Count variable still standardised) Script2
source('MCMC Functions/Jags-Ymet-XmetMulti-Mrobust(Norm-NONStand).R')

# Normal likelihood and all variables non-standardised Script3
source('MCMC Functions/Jags-Ymet-XmetMulti-Mrobust(Norm-RAW).R')

# JAGS Utilities
source('MCMC Functions/DBDA2E-utilities.R')

library(rjags)


# clean dataframe
bf = readRDS("BlackFriday_Clean.rds")

# all variables minus Purch_Avg
xName = names(bf)[-3]

# all variables minus Purch_Avg, years in city, marital status
xName2 = names(bf)[-c(1,3,11:14)]

# get train and test data (for xName)
set.seed(123)
random_rows = sample(nrow(bf),10)
test_bf = as.matrix(bf[random_rows,]) # keep only independent variables used in xName2
train_bf = bf[-random_rows,]

# get train and test data (for xName2)
set.seed(123)
random_rows = sample(nrow(bf),10)
test_bf = as.matrix(bf[random_rows,-c(1,11:14)]) # keep only independent variables used in xName2
train_bf = bf[-random_rows,]

reg_chains = genMCMC(train_bf, xName = xName2, yName = "Purch_Avg",
                     numSavedSteps = 100000, xPred = test_bf[,-2], thinSteps = 15)



parameterNames = varnames(reg_chains) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=reg_chains , parName=parName , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( reg_chains )
show(summaryInfo)

# Display posterior information:
plotMCMC( reg_chains , data=bf , xName=xName2 , yName="Purch_Avg" , 
          pairsPlot=TRUE , showCurve=FALSE )
#------------------------------------------------------------------------------- 

#saveRDS(reg_chains,file = "Chain Objects/Script2_100000steps_15thin_5000adapt_5000burn_xName_TEST-TRAIN.rds")



# extracting prediction info and finding MSE
model_predictions = summaryInfo[c("pred[1]","pred[2]","pred[3]","pred[4]",
                                  "pred[5]","pred[6]","pred[7]","pred[8]",
                                  "pred[9]","pred[10]"), c("Mode","HDIlow","HDIhigh")]

model_predictions = as.data.frame(cbind(test_bf[,2], model_predictions))
colnames(model_predictions)[1:2] = c("Actual","Prediction") 

model_predictions$Error = model_predictions$Actual - model_predictions$Prediction

MSE = mean(model_predictions$Error^2)
