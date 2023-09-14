# Script for ensemble Environmental Niche Modeling for red spruce Picea rubens Sarg. in Eastern North America
# described in:
# Lachmuth et al. (2023). Novel genomic offset metrics integrate local adaptation into habitat suitability 
# forecasts and inform assisted migration. Ecological Monographs

# written by S. Lachmuth at the Appalachian Lab, Frostburg, MD, USA 2020-2023
#
# Code is provided as is, without support 


# Requires additional code "varImportance.R" to be sourced below
# "varImportance.R" is only available on request 
# Not publicly provided as not originally written by Lachmuth et al.

# Load data ---------------------------------------------------------------
# Selected climate variables:
selVarsSDM<- c("CMD","DD_0","DD18","eFFP","EXT","MAR","MSP","PAS","PET","RH","TD") 

# Read red spruce occurence and climate data
dat_allClim<-read.csv(file = "clim_occurence_redSpruce.csv")
dat<-rubensDat_ecoReg[,c(names(rubensDat_ecoReg)[1:5],selVarsSDM)]
summary(dat)

# Reat spatial blocking folds
spatBlockFolds<-dget(paste0(datapath,"spatBlockFolds_2.5.Robj"))



# Fit and split-validate HSM --------------------------------------------

# Variables for storing results
nCV <- 5 # Number of cross-validations
nRow <- nrow(dat)
nMod<-5 # Number of models in ensemble
nVar<-length(selVarsSDM) # Number of climate variables
nStat<-3 # Number of evaluation statistics: AUC, TSS, Kappa

# Make table for parallelization
tab <- expand.grid(1:5,c("GLM","GAM","GBM","MARS","RF"))



## Start loop --------------------------------------------------------------
# Set up cluster
require(doParallel)
cl <- makeCluster(30)
registerDoParallel(cl)

# Run foreach loop
fitMod<-foreach(k = 1:nrow(tab), .packages=c("MASS","gam","gbm","earth","randomForest","biomod2","Hmisc","utils")) %dopar%{
  source("varImportance.R") #adjust path # Code can be provided on request.
  set.seed(42)

  cvNum <- tab[k,1]
  mod <- as.character(tab[k, 2])

  # Make list to store the results
  lll <- list(cv=cvNum, mod=mod, Best_models=NULL,Test_results=rep(0,3), Var_imp=rep(0,nVar), Pred_results=NULL)

  # Separate the original data in calibration and evaluation subsets spatBlockFolds
  fold<-spatBlockFolds[[cvNum]]
  calib<-dat[fold[[1]],]
  eval<-dat[fold[[2]],]



### GLM ---------------------------------------------------------------------
if(tab[k,]$Var2=="GLM"){
  require(MASS)
  require(biomod2)
  require(utils)
  require(Hmisc)
  
  glmStart <- glm(PresAbs~. , data=calib[,c("PresAbs",selVarsSDM)],
                  family=binomial, 
                  trace = FALSE, 
                  control=glm.control(maxit=100,epsilon = 1e-08))
  glm.formula <-makeFormula("PresAbs",calib[,selVarsSDM],"quadratic",interaction.level=0) 
  glmBest <- stepAIC(glmStart,
                     scope = list(lower = glmStart, upper = glm.formula),
                     data = calib, direction = "both", trace = FALSE, k = 2, control=glm.control(maxit=100,epsilon = 1e-08))
  lll$Best_models <- glmBest
  
  
  lll$Var_imp <- var_importance(data = calib[,6:ncol(calib)], model = glmBest, iterations_num = 10)[,2]
  
  
  # prediction on the evaluation data and evaluation using AUC, TSS and KAPPA 
  Pred_test <-  predict(glmBest, eval, type="response") 
  lll$Test_results[1] <- somers2(Pred_test,eval$PresAbs)["C"]
  Test_TSS<-Find.Optim.Stat(Stat = "TSS", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[2] <- Test_TSS[1,1]
  Test_KAPPA<-Find.Optim.Stat(Stat = "KAPPA", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[3] <- Test_KAPPA[1,1]
  
  # prediction on the total dataset
  lll$Pred_results <- predict(glmBest, dat, type="response")
  test<-data.frame(test="test")
  #write.csv(test,file="filename.csv") #adjust filename
}


### GAM ---------------------------------------------------------------------
if(tab[k,]$Var2=="GAM"){
  require(MASS)
  require(gam)
  require(utils)
  require(Hmisc)
  require(biomod2)
  require(dplyr)
  
  dir.create(paste0(datapath,"gamData_",cvNum,"/"),showWarnings = F)
  setwd(paste0(datapath,"gamData_",cvNum,"/"))
  dput(calib[,c("PresAbs",selVarsSDM)],"scopeStart.txt")
  gamStart <- gam(PresAbs ~ 1 ,  data=dget("scopeStart.txt"),
                  family=binomial) 
  # For gam warning see here (harmless):
  #https://stackoverflow.com/questions/57664927/warning-in-gam-with-release-of-r-3-6-1
  gamBest <-step.Gam(gamStart,
                     biomod2:::.scope(dget("scopeStart.txt"),"s",c(2,3)), 
                     trace=T, 
                     direction = "both")
  
  lll$Best_models <- gamBest
  setwd(wd)
  
  lll$Var_imp <- var_importance(data = calib[,6:ncol(calib)], model = gamBest, iterations_num = 10)[,2]
  
  # prediction on the evaluation data and evaluation using the AUC and TSS approach
  Pred_test <-  predict(gamBest, eval, type="response")
  lll$Test_results[1] <- somers2(Pred_test,eval$PresAbs)["C"]
  Test_TSS<-Find.Optim.Stat(Stat = "TSS", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[2] <- Test_TSS[1,1]
  Test_KAPPA<-Find.Optim.Stat(Stat = "KAPPA", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[3] <- Test_KAPPA[1,1]
  
  # prediction on the total dataset
  lll$Pred_results <- predict(gamBest, dat, type="response")
  test<-data.frame(test="test")
  #write.csv(test,file="filename.csv") #adjust filename
}


### GBM ---------------------------------------------------------------------
if(tab[k,]$Var2=="GBM"){
  require(MASS)
  require(gbm)
  require(biomod2)
  require(utils)
  require(Hmisc)
  gbmBest<- gbm(PresAbs~., 
                data=calib[,c("PresAbs",selVarsSDM)], 
                distribution  ="bernoulli",
                interaction.depth=(nVar+1), 
                shrinkage = 0.01, 
                n.trees=10000,
                cv.folds=10, 
                bag.fraction = 0.4)  #new values May 21
  lll$Best_models <- gbmBest
  lll$Var_imp <- var_importance(data = calib[,6:ncol(calib)], model = gbmBest, iterations_num = 10)[,2]
  
  
  # prediction on the evaluation data and evaluation using the AUC and TSS approach
  Pred_test <-  predict(gbmBest, eval, type="response") 
  lll$Test_results[1] <- somers2(Pred_test,eval$PresAbs)["C"]
  Test_TSS<-Find.Optim.Stat(Stat = "TSS", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[2] <- Test_TSS[1,1]
  Test_KAPPA<-Find.Optim.Stat(Stat = "KAPPA", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[3] <- Test_KAPPA[1,1]
  
  
  # prediction on the total dataset
  lll$Pred_results <- predict(gbmBest, dat, type="response")
  test<-data.frame(test="test")
  #write.csv(test,file="filename.csv") #adjust filename
}



### MARS --------------------------------------------------------------------
if(tab[k,]$Var2=="MARS"){
  require(MASS)
  require(earth)
  require(biomod2) 
  require(utils)
  require(Hmisc)
  # default: nfold = 0, no additonal cross-validation
  marsBest <- earth(PresAbs ~ ., 
                    data = calib[,c("PresAbs",selVarsSDM)],
                    degree = 1, 
                    pmethod = 'exhaustive',
                    nfold = 10,
                    glm = list(family = binomial))
  lll$Best_models <- marsBest
  lll$Var_imp <- var_importance(data = calib[,6:ncol(calib)], model = marsBest, iterations_num = 10)[,2]
  
  # prediction on the evaluation data and evaluation using the AUC and TSS approach
  Pred_test <-  predict(marsBest, eval, type="response")
  lll$Test_results[1] <- somers2(Pred_test, eval$PresAbs)["C"]
  Test_TSS<-Find.Optim.Stat(Stat = "TSS", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[2] <- Test_TSS[1,1]
  Test_KAPPA<-Find.Optim.Stat(Stat = "KAPPA", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[3] <- Test_KAPPA[1,1]
  
  
  # prediction on the total dataset
  lll$Pred_results <- predict(marsBest, dat, type="response")
  test<-data.frame(test="test")
  #write.csv(test,file="filename.csv") #adjust filename
}


### Random forest -----------------------------------------------------------
if(tab[k,]$Var2=="RF"){
  require(randomForest)
  require(biomod2) 
  require(utils)
  require(Hmisc)
  
  RFBest = randomForest(x = calib[,selVarsSDM], 
                        y = as.factor(calib$PresAbs), 
                        ntree = 5000, 
                        keep.forest=T, importance = TRUE)
  
  lll$Best_models <-RFBest
  lll$Var_imp <- var_importance_RF(data = calib[,6:ncol(calib)], model = RFBest, iterations_num = 10)[,2]
  
  # prediction on the evaluation data and evaluation using the AUC and TSS approach
  Pred_test <-  predict(RFBest, eval, type="prob")[,2]
  
  lll$Test_results[1] <- somers2(Pred_test,eval$PresAbs)["C"]  
  Test_TSS<-Find.Optim.Stat(Stat = "TSS", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[2] <- Test_TSS[1,1]
  Test_KAPPA<-Find.Optim.Stat(Stat = "KAPPA", Fit = Pred_test*1000, Obs = eval$PresAbs )
  lll$Test_results[3] <- Test_KAPPA[1,1]
  
  # prediction on the total dataset
  lll$Pred_results <- predict(RFBest, dat, type="prob")[,2]
  test<-data.frame(test="test")
  #write.csv(test,file="filename.csv") #adjust filename
}

return(lll)

}  

stopCluster(cl)


# Save
#save(fitMod,file="fitMod.RData"))

