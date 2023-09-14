# redSpruce_ecolMonographs
R code and formated data sets from: Susanne Lachmuth, Thibaut Capblancq, Anoob Prakash, Stephen R. Keller, Matthew C. Fitzpatrick. Novel genomic offset metrics integrate local adaptation into habitat suitability forecasts and inform assisted migration. Ecological Monographs (accepted).

***

## DATA
For the availibilty of raw data see "Open Research Statement" in the manuscript.

## List of files in the /data folder:

- **commonGarden.csv:** *Red spruce height growth population means from three common gardens as well as (raw and scaled) genomic offsets and climate transfer distances*
- **clim_occurence_redSpruce.csv:** *Red spruce presence-absence and climate data used for Environmental Niche Modeling*
- **clim_sprucePops.csv:** *Climate data of our sampled red spruce populations. See "Methods" section of the manuscript and Appendices S1 and S2 for details on the data sources.*
- **clim_sprucePres.csv:** *Climate data for extant red spruce populations. See "Methods" section of the manuscript and Appendices S1 for details on the data sources.*
- **futClim_studyArea.csv:** *End of 21st century climate data (under SSP5-85) for our study area. See "Methods" section of the manuscript for details on the original data source.*
- **snpScores_sprucePops.csv:** *Spatial blocking folds for split-validation in Environmental Niche Modeling.*
- **spatBlockFolds_2.5.Robj:** *SNP data at candidate loci for our sampled populations (see "Methods" section of the Manuscripts for details on the selection of candidate loci).*


***

## R SCRIPTS
We provide code for fitting an ensemble of Envionmental Niche Models for red spruce in eastern North America, Gradient Forest modeling of genomic variation, calculating and scaling genomic offsets, and evaluating their predictive performance in three common gardens. Calculation of Donor and Recipient Importance is included in the 'gfScaledOffsets' script.

***

- **commonGarden.R:** *Evaluation of the predictive performance of raw and scaled genomic offsets as well as climate transfer distances in three common gardens*
- **ENM.R:** * Fitting an ensemble of Envionmental Niche Models for red spruce in eastern North America*
- **gfScaledOffsets.R:** * Gradient Forest modeling of genomic variation, scaling of genomic offsets, and calculation of Donor and Recipient Importance*

***
