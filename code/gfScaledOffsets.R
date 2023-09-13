# Scripts to model genomic variation, calculate and scale genomic offsets and calculate Donor and Recipient Importance
# described in Lachmuth et al. (2023) Ecological Monographs

# written by S. Lachmuth at the Appalachian Lab, Frostburg, MD, USA 2020-2023
#
# Code is provided as is, without support 



# Settings ----------------------------------------------------------------
nCores<-50
nBreak<-100



# Load required R packages ------------------------------------------------
library(adehabitatLT)
library(extendedForest) 
library(fields)
library(gradientForest)
library(reshape)
library(parallel)
library(doParallel)


# GRADIENT FOREST MODELING -----------------------------------------------------

# Read climatic data for sampled red spruce populations
climGF<-read.csv(file = paste0(datapath,"/forGitHub_MS1/clim_sprucePops.csv"))


# Read allele frequencies of 240 candidate SNPs for sampled red spruce populations
snpScores<-read.csv(file = paste0(datapath,"/forGitHub_MS1/snpScores_sprucePops.csv")) # pops are in same order as in climate data file


# Model of allele frequency turnover along climatic gradients
maxLevel <- log2(0.368*nrow(climGF)/2) #account for correlations, see ?gradientForest 

# Fit gf models for candidate SNPs 
gfMod <- gradientForest(cbind(climGF, snpScores), predictor.vars=colnames(climGF),
                           response.vars=colnames(snpScores), ntree=5000, 
                           maxLevel=maxLevel, trace=T, corr.threshold=0.5)




# CONTEMPORARY SPATIAL OFFSETS --------------------------------------------
# Read climate data for all geographic grid cells (2.5 arcmin) with red spruce presence records

clim_sprucePres_xy<-fread(paste0(datapath,"forGitHub_MS1/clim_sprucePres.csv"), stringsAsFactors = T, header = T)
# just climate data:
clim_sprucePres<-clim_sprucePres_xy[,-(1:3)]


# transform climate data using gradient forest model predictions
transClim<-predict(gfMod,clim_sprucePres)


# Calculate spatial offsets
# run in parallel
breakIt <- split(1:nrow(transClim), cut(1:nrow(transClim), nCores, labels = FALSE))


cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
spatOffset <- foreach(i = 1:length(breakIt), .packages=c("fields")) %dopar%{
  
  # Calculate pairwise offset matrix
  mat<-rdist(transClim[breakIt[[i]],],transClim)
  return(mat)
}
stopCluster(cl)

spatOffset_list<- do.call(rbind, spatOffset)
spatOffset_df<-as.data.frame(spatOffset_list)
colnames(spatOffset_df)<-rownames(transClim)


# melt into long format
require(reshape)
spatOffset_mat<-as.matrix(spatOffset_df)
spatOffset_lon <- melt(spatOffset_mat)[melt(lower.tri(spatOffset_mat))$value,]


# Get empirical cumulative density function
ecdfSpatOffset<-ecdf(spatOffset_lon$value)


# SPATIO-TEMPORAL OFFSETS -------------------------------------------------
# Spatio-temporal offsets describe the G - E disruption for a (donor) population being 
# transferred from one grid cell under current climate to any other grid cell (recipient) under future climate.
# Here, we use all grid cells with extant spruce populations as donors and all grid cells in eco-regions
# with extant spruce populations as well as adjacent eco-regions as recipients.


# Read current climate:
clim_sprucePres_xy<-fread(paste0(datapath,"forGitHub_MS1/clim_sprucePres.csv"), stringsAsFactors = T, header = T)
# just climate data:
clim_sprucePres<-clim_sprucePres_xy[,-(1:3)]


# Read future climate:
futClim_studyArea_xy<-fread(paste0(datapath,"/forGitHub_MS1/futClim_studyArea.csv"), stringsAsFactors = T, header = T)
head(futClim_studyArea_xy)
# just climate data:
futClim_studyArea<-futClim_studyArea_xy[,-(1:3)]


## GF transform climate data
# Current transformed climate (coordinates required for mapping - not part of this script)
currClim_trans_xy <- data.frame(clim_sprucePres_xy[,c("x","y","cellindex")], predict(gfMod,clim_sprucePres))
currClim_trans_xy <- currClim_trans_xy[order(currClim_trans_xy$cellindex),]

# Future transformed climate (coordinates required for mapping - not part of this script)
futClim_trans_xy <- data.frame(futClim_studyArea_xy[,c("x","y","cellindex")], predict(gfMod,futClim_studyArea))
futClim_trans_xy <- futClim_trans_xy[order(futClim_trans_xy$cellindex),]


# Remove xy coordinates and cell index
currClim_trans<-currClim_trans_xy[,-(1:3)]
futClim_trans <-futClim_trans_xy[,-(1:3)]


### Calculate raw offsets between all donor (current) and recipient (future) cells ------------------
breakItAll <- split(1:nrow(futClim_trans), cut(1:nrow(futClim_trans), nBreak, labels = FALSE))

# Run in parallel:
cl <- makeCluster(nCores)
registerDoParallel(cl)

allOffset <- foreach(i = 1:length(breakItAll), .packages=c("fields")) %dopar%{
  
  # Calculate pairwise offset matrix
  mat_allOffset<-rdist(futClim_trans[breakItAll[[i]],],currClim_trans)
  return(mat_allOffset)
  
}
stopCluster(cl)

# Call and format data
allOffset_list<- do.call(rbind, allOffset)
allOffset_df<-as.data.frame(allOffset_list)
# Name columns and rows according to cell indices (important for mapping later on (not part of this script))
colnames(allOffset_df)<-futClim_trans_xy$cellindex 
rownames(allOffset_df)<-futClim_trans_xy$cellindex 



# clear memory
rm(spatOffset)
rm(spatOffset_list)
rm(allOffset)
rm(allOffset_list)



# SCALE SPATIO-TEMPORAL OFFSETS -------------------------------------

## Re-express raw offsets as probability based on spatial offset ecdf  -------------

# Prep. for loop
allOffset_prob<-allOffset_df
rm(allOffset_df)

breakIt_prob<- split(1:nrow(allOffset_prob), cut(1:nrow(allOffset_prob), nBreak, labels = FALSE))

# Run in parallel
cl <- makeCluster(nCores)
registerDoParallel(cl)

allOffset_prob_loop <- foreach(i = 1:length(breakIt_prob), .packages=c("stats")) %dopar%{
  allOffset_prob_sub<-allOffset_prob[breakIt_prob[[i]],]
  for(j in 1:ncol(allOffset_prob_sub)) {
    allOffset_prob_sub[ , j] <- ecdfSpatOffset(allOffset_prob_sub[ , j]) 
  } 
  return(allOffset_prob_sub)
}
stopCluster(cl)

# Call and format data
allOffset_prob_list<-do.call(rbind, allOffset_prob_loop)
allOffset_prob_df<-as.data.frame(allOffset_prob_list)

rm(allOffset_prob)
rm(allOffset_prob_list)


## Re-express raw offsets as quantiles (sigmas) of chi distribution with 1 degree of freedom (half-normal distribution)  -------------
# Prep. for loop
allOffset_sigma<-allOffset_prob_df
rm(allOffset_prob_df)

breakIt_sigma <- split(1:nrow(allOffset_sigma), cut(1:nrow(allOffset_sigma), nBreak, labels = FALSE))

# Run in parallel
cl <- makeCluster(nCores)
registerDoParallel(cl)

allOffset_sigma_loop <- foreach(i = 1:length(breakIt_qchi), .packages=c("adehabitatLT")) %dopar%{
  allOffset_sigma_sub<-allOffset_sigma[breakIt_qchi[[i]],]
  for(j in 1:ncol(allOffset_sigma_sub)) {       
    allOffset_sigma_sub[ , j] <- adehabitatLT::qchi(allOffset_sigma_sub[ , j],1) 
  } 
  return(allOffset_sigma_sub)
}
stopCluster(cl)


# Call and format data
allOffset_sigma_list<-do.call(rbind, allOffset_sigma_loop)
allOffset_sigma_df<-as.data.frame(allOffset_sigma_list)



# DONOR & RECIPIENT IMPORTANCE --------------------------------------------

## Apply sigma threshold ---------------------------------------------------
# Set "not to exceed" offset threshold (here: 1 sigma) and calculate transferability matrix
sigmaTH<-1


# Prep. for loop
allOffset_sigmaTH<- as.data.frame(matrix(NA,ncol=ncol(allOffset_sigma_df),nrow=nrow(allOffset_sigma_df))) 
names(allOffset_sigmaTH)<-names(allOffset_sigma_df)

breakIt_sigma <- split(1:nrow(allOffset_sigma), cut(1:nrow(allOffset_sigma), nBreak, labels = FALSE))

# Run in parallel
cl <- makeCluster(nCores)
registerDoParallel(cl)

allOffset_sigmaTH_loop <- foreach(i = 1:length(breakIt_sigma)) %dopar%{
  allOffset_sigmaTH_sub<-allOffset_sigma[breakIt_sigma[[i]],]
  for  (j in 1:ncol(allOffset_sigmaTH_sub)){
    allOffset_sigmaTH_sub[,j]<-ifelse(allOffset_sigmaTH_sub[,j]<=sigmaTH,1,0)
  }
  return(allOffset_sigmaTH_sub)
}
stopCluster(cl)

# Call and format data
allOffset_sigmaTH_list<- do.call(rbind, allOffset_sigmaTH_loop)
allOffset_sigmaTH_df<-as.data.frame(allOffset_sigmaTH_list)
row.names(allOffset_sigmaTH_df)<-row.names(allOffset_sigma_df)




## Calculate Donor Importance ----------------------------------------------
# Here: recipient area = entire study area
# For smaller donor or recipient areas subset allOffset_sigmaTH_df accordingly (based on cell indices)

# Make Donor Importance dataframe
# Cell indices and xy coordinates are required for mapping (not part of this script)
donImp_xy<-data.frame(clim_sprucePres_xy[,c("x","y","cellindex")])

# Calculate donor importance
donImp_xy$donImp<-colSums(allOffset_sigmaTH_df,na.rm = T) 
donImp_xy$percDonImp<-donImp_xy$donImp/nrow(allOffset_sigmaTH_df)*100 # as percentage



## Calculate Recipient Importance ------------------------------------------
# Here: recipient area = entire study area
# For smaller donor or recipient areas subset allOffset_sigmaTH_df accordingly (based on cell indices)

# Make Recipient Importance dataframe
# Cell indices and xy coordinates are required for mapping
recImp_xy<-data.frame(futClim_studyArea_xy[,c("x","y","cellindex")]) 

# Calculate Recipient Importance
recImp_xy$recImp<-rowSums(allOffset_sigmaTH_df,na.rm = T)
recImp_xy$percRecImp<-recImp_xy$recImp/ncol(allOffset_sigmaTH_df)*100 # as percentage
