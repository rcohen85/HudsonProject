setwd("C:/Users/rec297/Documents/GitHub/HudsonProject")
library(soundecology)
library(seewave)
library(parallel)
library(tuneR)
library(stringr)
library(lubridate)
source("ACI_RC.R")
source("multiple_sounds_RC.R")

# Use soundecology package to calculate Bioacoustic Index (Boelman et al. 2007), 
# Acoustic Diversity Index (Villanueva-Rivera et al. 2011), Acoustic Evenness Index
# (Villanueva-Rivera et al. 2011), and Acoustic Complexity Index (Pieretti et al. 2011)
# modified to debug the acoustic_complexity function


dep = c('Test')
inDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/HawthorneValley/FEP_MVB_2022'
outDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/HawthorneValley/FEP_MVB_2022/AcousticIndices'

for (i in 1:length(dep)){

# multiple_sounds_RC(directory=paste(inDir,'/',dep[i],sep=""),flac=TRUE,
#                 resultfile = paste(outDir,'/',dep[i],'_BI_5_10K.csv',sep=""),
#                 soundindex = "bioacoustic_index",
#                 min_freq = 5000, max_freq = 10000,fft_w=1024,no_cores=4)
#   
# multiple_sounds_RC(directory=paste(inDir,'/',dep[i],sep=""),flac=TRUE,
#                 resultfile = paste(outDir,'/',dep[i],'_BI_2_5K.csv',sep=""),
#                 soundindex = "bioacoustic_index",
#                 min_freq = 2000, max_freq = 5000,fft_w=1024,no_cores=4)

# 
# multiple_sounds_RC(directory=paste(inDir,'/',dep,sep=""),
#                 resultfile = paste(outDir,'/',dep,'ADI.csv',sep=""),
#                 soundindex = "acoustic_diversity", max_freq = 4000,no_cores=8)
# 
# multiple_sounds_RC(directory=paste(inDir,'/',dep,sep=""),
#                 resultfile = paste(outDir,'/',dep,'AEI.csv',sep=""),
#                 soundindex = "acoustic_evenness", max_freq = 4000,no_cores=8)
# 
  browser()
multiple_sounds_RC(directory=paste(inDir,'/',dep[i],sep=""),flac=TRUE,
                resultfile = paste(outDir,'/',dep[i],'_ACI_1s.csv',sep=""),
                soundindex = "acoustic_complexity", min_freq = 2000, max_freq = 15000,no_cores=1,dur=60,j=1,fft_w=1024)

# 
# multiple_sounds_RC(directory=inDir,
#                    resultfile = paste(outDir,'/',dep,'H_30s_Fish.csv',sep=""),
#                    soundindex = "H",no_cores=1,dur=30)
# multiple_sounds_RC(directory=inDir,
#                    resultfile = paste(outDir,'/',dep,'SH_30s_NoFish.csv',sep=""),
#                    soundindex = "SH",no_cores=1,dur=30)
# multiple_sounds_RC(directory=inDir,
#                    resultfile = paste(outDir,'/',dep,'TH_30s_Fish.csv',sep=""),
#                    soundindex = "TH",msmooth=c(10,50),no_cores=1,dur=30)

}
