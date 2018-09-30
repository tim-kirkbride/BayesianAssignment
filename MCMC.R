
# Normal likelihood and standardised variables (independent & dependent) Script1
source('MCMC Functions/Jags-Ymet-XmetMulti-Mrobust(Norm-Stand).R')

# Normal likelihood and non-standardised binary variables (dependent and Count variable still standardised) Script2
source('MCMC Functions/Jags-Ymet-XmetMulti-Mrobust(Norm-NONStand).R')

# Normal likelihood and all variables non-standardised Script3
source('MCMC Functions/Jags-Ymet-XmetMulti-Mrobust(Norm-RAW).R')

# JAGS Utilities
source('MCMC Functions/DBDA2E-utilities.R')

library(rjags)

# made up prediction (all independent variables)
xPred = matrix(c(0,15,1,1,0,0,0,1,0,1,1,0,0,
                 0,130,0,1,0,1,0,1,0,0,1,0,1), nrow = 2, ncol = 13, byrow = TRUE,
               dimnames = list(c("pred1","pred2"),
                               c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
                                 "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
                                 "beta[11]","beta[12]","beta[13]")))

# made up prediction (minus years in city, marital status)
xPred2 = matrix(c(15,1,1,0,0,0,1,0,
                  130,0,1,0,1,0,1,0), nrow = 2, ncol = 8, byrow = TRUE,
               dimnames = list(c("pred1","pred2"),
                               c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
                                 "beta[6]","beta[7]","beta[8]")))


bf = readRDS("BlackFriday_Clean.rds")

# all variables minus Purch_Avg
xName = names(bf)[-3]

# all variables minus Purch_Avg, years in city, marital status
xName2 = names(bf)[-c(1,3,11:14)]

reg_chains = genMCMC(bf, xName = xName2, yName = "Purch_Avg",
                     numSavedSteps = 100000, xPred = xPred2, thinSteps = 15)



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

#saveRDS(reg_chains,file = "Chain Objects/Script2_100000steps_15thin_5000adapt_5000burn_xName2.rds")
