source('Jags-Ymet-XmetMulti-Mrobust(Norm).R')
source('DBDA2E-utilities.R')


xPred = matrix(c(0,15,1,1,0,0,0,1,0,1,1,0,0), nrow = 1, ncol = 13, byrow = TRUE,
               dimnames = list(c("pred1"),
                               c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
                                 "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
                                 "beta[11]","beta[12]","beta[13]")))


bf = as.data.frame(bf[,-"Purch_Sum"])
xName = names(bf)[-3]



reg_chains = genMCMC(bf, xName = xName, yName = "Purch_Avg",
                     numSavedSteps = 1000, xPred = xPred)



parameterNames = varnames(reg_chains) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=reg_chains , parName=parName , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( reg_chains )
show(summaryInfo)

# Display posterior information:
plotMCMC( reg_chains , data=bf , xName=xName , yName="Purch_Avg" , 
          pairsPlot=TRUE , showCurve=FALSE )
#------------------------------------------------------------------------------- 

#saveRDS(reg_chains,file = "Output Chains/Part B/reg_chains4000000(2).rds")
