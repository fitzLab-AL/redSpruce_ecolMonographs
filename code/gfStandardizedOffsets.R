# Scripts to model genomic turnover, calculate and standardize genomic offsets and calculate Donor and Recipient Importance
# described in Lachmuth et al. (resubmitted) Ecological Monographs

# written by S Lachmuth at the Appalachian Lab, Frostburg, MD, USA 2020-2023
#
# Code is provided as is, without support 



# Settings ----------------------------------------------------------------
nCores<-50



# Load required R packages ------------------------------------------------
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

# New error message: .Error in rbind(structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#                              number of columns of matrices must match (see arg 51)



# CONTEMPORARY SPATIAL OFFSETS --------------------------------------------
# Read climate data for all geographic grid cells (2.5 arcmin) with red spruce presence records

# TO-DO: provide just selected climate data!!!

clim_sprucePres_xy<-fread(paste0(datapath,"forGitHub_MS1/clim_sprucePres.csv"), stringsAsFactors = T, header = T)
# just climate data:
clim_sprucePres<-clim_sprucePres_xy[,-(1:3)]


# transform climate data with gradient forest model
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

spatOffsetGF<- do.call(rbind, spatOffset)
spatOffset_df<-as.data.frame(spatOffsetGF)
colnames(spatOffset_df)<-rownames(transClim)
spatOffset_df[1:10,1:10]


# melt into long format
require(reshape)
spatOffset_mat<-as.matrix(spatOffset_df)
spatOffset_lon <- melt(spatOffset_mat)[melt(lower.tri(spatOffset_mat))$value,]


# Get empirical cumulative density function
ecdfSpatOffset<-ecdf(spatOffset_lon$value)


# SPATIO-TEMPORAL OFFSETS -------------------------------------------------



# STANDARDIZE SPATIO-TEMPORAL OFFSETS -------------------------------------


## Re-express raw offsets as quantiles of spatial offset ecdf  -------------


## Re-express raw offsets as quantiles (sigmas) of chi distribution with 1 degree of freedom (half-normal distribution)  -------------



# DONOR & RECIPIENT IMPORTANCE --------------------------------------------

## Apply sigma threshold ---------------------------------------------------
# Set "not to exceed" offset threshold (here: 1 sigma) and calculate transferability matrix


## Calculate donor importance ----------------------------------------------


## Calculate recipient importance ------------------------------------------


