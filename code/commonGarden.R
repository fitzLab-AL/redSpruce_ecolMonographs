# Script for common-garden analyses to evaluate the predictive performance of scaled genomic offsets as
# described in:
# Lachmuth et al. (2023). Novel genomic offset metrics integrate local adaptation into habitat suitability 
# forecasts and inform assisted migration. Ecological Monographs

# written by S. Lachmuth at the Appalachian Lab, Frostburg, MD, USA 2020-2023
#
# Code is provided as is, without support 


# Read garden data --------------------------------------------------------
gardDat<-read.csv(file = "commonGarden.csv")


# Population and garden level models ------------------------------------------
summary(gardDat)
predictor<-c("rawClimDist","rawClimDist_sigma","Offset","Offset_sigma")
garden<-c("MD","NC","VT")


# Make array to store models
arr_models <- array(list(type = any),c(length(predictor),length(garden)), dimnames=list(predictor, garden)) 


# Prepare results table
results_garden<-data.frame(Predictor=character(0),Garden=character(0),
                               intercept_pop=numeric(0),slope_pop=numeric(0),R2_pop=numeric(0),R2adj_pop=numeric(0),
                               F_pop=numeric(0), numdf_pop=numeric(0), dendf_pop=numeric(0),
                               p_pop=numeric(0), AIC_pop=numeric(0))



for(i in 1:length(predictor)){
  for(j in 1:length(garden)){
  
      if(predictor[i]=="rawClimDist"){
        LM<-lm(Growth~rawClimDist, data=gardDat[gardDat$Garden==garden[j],])
      }   
      if(predictor[i]=="rawClimDist_sigma"){
        LM<-lm(Growth~rawClimDist_sigma, data=gardDat[gardDat$Garden==garden[j],])
      }
      if(predictor[i]=="Offset"){
        LM<-lm(Growth~Offset, data=gardDat[gardDat$Garden==garden[j],])
      }   
      if(predictor[i]=="Offset_sigma"){
        LM<-lm(Growth~Offset_sigma, data=gardDat[gardDat$Garden==garden[j],])
      }
      
      intercept_pop<-summary(LM)$coefficients[1,1]
      slope_pop<-summary(LM)$coefficients[2,1]
      R2_pop<-summary(LM)$r.squared
      R2adj_pop<-summary(LM)$adj.r.squared
      AIC_pop<-AIC(LM)
      F_pop<-round(summary(LM)$fstatistic[1],3)
      numdf_pop<-summary(LM)$fstatistic[2]
      dendf_pop<-summary(LM)$fstatistic[3]
      p_pop<-dropterm(LM,test="F")$'Pr(F)'[2]
      
      
      # Save models
      arr_models[[predictor[i], garden[j]]]<-LM
      
      # Save stats
      results_garden<-rbind(pop_results_garden,data.frame(Predictor=predictor[i],Garden=garden[j],
                                                              intercept_pop=intercept_pop,slope_pop=slope_pop,R2_pop=R2_pop,R2adj_pop=R2adj_pop,
                                                              F_pop=F_pop, numdf_pop=numdf_pop, dendf_pop=dendf_pop,p_pop=p_pop, AIC_pop=AIC_pop))
      
    }
  }



