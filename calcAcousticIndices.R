library(soundecology)
library(seewave)
library(parallel)
library(tuneR)
source("ACI_RC.R")
source("multiple_sounds_RC.R")

# Use soundecology package to calculate Bioacoustic Index (Boelman et al. 2007), 
# Acoustic Diversity Index (Villanueva-Rivera et al. 2011), Acoustic Evenness Index
# (Villanueva-Rivera et al. 2011), and Acoustic Complexity Index (Pieretti et al. 2011)
# modified to debug the acoustic_complexity function

dep = 'BC04'
inDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/AcousticData_WAVE/BlackCreek/Swift/WAVE'
outDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/SoundscapeMetrics/BlackCreek'

# multiple_sounds_RC(directory=paste(inDir,'/',dep,sep=""),
#                 resultfile = paste(outDir,'/',dep,'BI.csv',sep=""),
#                 soundindex = "bioacoustic_index", max_freq = 4000,no_cores=8)
# 
multiple_sounds_RC(directory=paste(inDir,'/',dep,sep=""),
                resultfile = paste(outDir,'/',dep,'ADI.csv',sep=""),
                soundindex = "acoustic_diversity", max_freq = 4000,no_cores=8)

multiple_sounds_RC(directory=paste(inDir,'/',dep,sep=""),
                resultfile = paste(outDir,'/',dep,'AEI.csv',sep=""),
                soundindex = "acoustic_evenness", max_freq = 4000,no_cores=8)

# multiple_sounds_RC(directory=paste(inDir,'/',dep,sep=""), 
#                 resultfile = paste(outDir,'/',dep,'ACI.csv',sep=""),
#                 soundindex = "acoustic_complexity", max_freq = 4000,no_cores=8)

