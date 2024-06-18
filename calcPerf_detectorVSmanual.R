library(stringr)
library(dplyr)
library(lubridate)
library(tuneR)
library(ggplot2)
library(caret)
library(pracma)

# Assumes annotation files & detection files each correspond to a single test audio file

# DCLDE, original data
# anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/DCLDE_Data/Test/SelectionTables/Ch01_wTimestamps'
# audioDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/DCLDE_Data/Test/Channel_1_wTimestamps'
anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/DCLDE_Data/Test/SelectionTables/Ch01'
audioDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/DCLDE_Data/Test/Channel_1'

# DCLDE, freq shifted
# anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/DCLDE_Data/Test/SelectionTables/Ch01_wTimestamps_Shifted'
# audioDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/DCLDE_Data/Test/Channel_1_wTimestamps_Shifted'

# MA, original data
# anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MA_Data/Test/SelectionTables/ExhaustivelyAnnotated'
# audioDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MA_Data/Test/SoundFiles/ExhaustivelyAnnotated'

# MA, freq shifted
# anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MA_Data/Test/ShiftedSelectionTables/ExhaustivelyAnnotated'
# audioDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MA_Data/Test/ShiftedSoundFiles/ExhaustivelyAnnotated'

# MD, original data
# anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MD_Data/Test/IndividualSelTables'
# audioDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MD_Data/Test/SoundFiles'

# MD, freq shifted
# anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MD_Data/Test/ShiftedSelTables'
# audioDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/MD_Data/Test/ShiftedSoundFiles'

# Koogu
# detDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/Koogu/Testing/20240530-132959/Detections_Thr025'
# saveDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/Koogu/Testing/20240530-132959'

# BirdNET
# detDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/BirdNET/Testing/V2.4/20240530_3/Detections'
# saveDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/BirdNET/Testing/V2.4/20240530_3'

detDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/Shyam_NARWDet_Test/Detections'
saveDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/Shyam_NARWDet_Test'

binSize = 3 # duration of detector bins (s); if -1, will estimate detector bin size from detector output
minO = 0.4 # minimum temp1oral overlap (% of annotation duration) to count a detection window as containing a manually annotated call
chans = 1 # either specify channel(s) of interest, or leave as an empty list to analyze all available channels
SNRcol = c()#SNR.NIST.Quick..dB.'

##### 
annotation_files <- dir(anotDir,pattern='.txt')
detection_files <- list.files(detDir, pattern = ".txt")
audio_files <- dir(audioDir)

# # Get durations of files in list - ONLY WORKS FOR WAVE/MP3 FILES
durs = numeric()
for (i in 1:numel(audio_files)){
  durs = rbind(durs, (readWave(paste(audioDir,'/',audio_files[i],sep=""),header=TRUE)$samples/readWave(paste(audioDir,'/',audio_files[i],sep=""),header=TRUE)$sample.rate))}

allTimeBins = numeric()
allData = numeric()

for (i in 1:length(audio_files)) { # step through each test audio file and evaluate TP/FP/TN/FN
  
  if (str_detect(audio_files[i],'.wav')){
    ext = '.wav'
  } else if (str_detect(audio_files[i],'.flac')){
    ext = '.flac'
  } else if (str_detect(audio_files[i],'.aif')){
    ext = '.aif'
  }
  
  # Find corresponding ground truth annotations and detector output selection table
  # thisTimeStamp = str_extract(annotation_files[i],'\\d{8}_\\d{6}')
  matchingDets = str_which(detection_files,str_remove(audio_files[i],ext))
  matchingAnnots = str_which(annotation_files,str_remove(audio_files[i],ext))
  
  if (!isempty(matchingAnnots)){
    file_annotation <- paste(anotDir, annotation_files[matchingAnnots], sep = "/")
    annotations <-   read.table(file_annotation, header = TRUE, sep = "\t")
    
    if (dim(annotations)[1]>0){
      # Remove annotation type 'waveform'
      if ("View" %in% colnames(annotations) & "Spectrogram 1" %in% annotations$View & "Waveform 1" %in% annotations$View){
        annotations <- subset(annotations, annotations$View == "Spectrogram 1")
      }
      if (!("Delta.Time..s." %in% colnames(annotations))){
        annotations$Delta.Time..s. = annotations$End.Time..s. - annotations$Begin.Time..s.
      }
      
      # Sort into ascending start time order, reset row numbers
      annotations = annotations[order(annotations$Begin.Time..s.),]
      rownames(annotations) = 1:nrow(annotations)
    }
  } else {
    annotations = data.frame()
  }
  
  if (!isempty(matchingDets)){
    # Load detections
    file_detec <-  paste(detDir,detection_files[matchingDets],sep = "/")
    detections <- read.table(file_detec, header = TRUE, sep = "\t")
    
    if (dim(detections)[1]>0){
      # Remove detection type 'waveform'
      if (("View" %in% colnames(detections) & "Spectrogram 1" %in% detections$View & "Waveform 1" %in% detections$View)){
        detections <- subset(detections, detections$View == "Spectrogram 1")
      } # add Delta Time column if not already present
      if (!("Delta.Time..s." %in% colnames(detections))){
        detections$Delta.Time..s. = detections$End.Time..s. - detections$Begin.Time..s.
      }
      if (!("Tags" %in% colnames(detections)) & 'Species.Code' %in% colnames(detections)){
        detections$Tags = detections$Species.Code
      }
      if (!("Tags" %in% colnames(detections)) & 'Label' %in% colnames(detections)){
        detections$Tags = detections$Label
      }
      if (!("Score" %in% colnames(detections))){
        detections$Score = detections$Confidence
      }
      noCallInd = str_which(detections$Tags,'nocall')
      if (!isempty(noCallInd)){
      detections = detections[-noCallInd,]}
      
      if (nrow(detections)>0){
      # Sort into ascending start time order, reset row numbers
      detections = detections[order(detections$Begin.Time..s.),]
      rownames(detections) = 1:nrow(detections)
      } else {
        detections = data.frame()
      }
    }
  } else {
    detections = data.frame()
  }
  
  
  # Establish time bins consistent with how the detector saw the data 
  if (binSize>0){
    timeBins = seq(0,durs[i],by=binSize)
  }else{
    timeBins = seq(0,durs[i],by=round(min(detections$End.Time..s.-detections$Begin.Time..s.),digits=2))}
  allTimeBins = c(allTimeBins,timeBins)
  
  # Trim any annotations that exceed the last time bin 
  # (there may be a tiny bit of data at the end of the file which is not seen by the detector because it's not long enough to be a full spectrogram)
  tooLong = which(annotations$End.Time..s.>timeBins[length(timeBins)])
  annotations$End.Time..s[tooLong] = timeBins[length(timeBins)]
  
  # Determine full set of channels any annotations or detections exist in
  if (isempty(chans)){
  allChans = sort(unique(c(annotations$Channel,detections$Channel)))
  } else { allChans = chans}
  
  temp2 = numeric()
  
  for (j in 1:length(allChans)){
    
    temp1 = data.frame(Time=timeBins)
    temp1$File=audio_files[i]
    temp1$Channel=j
    temp1$SNR=NA
    temp1$ALabel=NA
    temp1$DLabel=NA
    temp1$DScore=NA
    
    anInd = which(annotations$Channel==allChans[j])
    detInd = which(detections$Channel==allChans[j])
    
    # for each time bin, note existing annotation and/or detection labels, mark if TP/FP/TN/FN
    for (k in 1:(length(timeBins)-1)){
      
      # find any annotations which overlap with this bin
      whichAnInds = which(annotations$Begin.Time..s.[anInd]<timeBins[k+1] & annotations$End.Time..s[anInd]>timeBins[k]) 
      Alabels = character()
      SNR = numeric()
      if (length(whichAnInds)>0){
        for (l in 1:length(whichAnInds)){
          # determine if overlap is sufficient
          # overlap = as.numeric(difftime(min(timeBins[k+1],annotations$End.Time..s.[anInd[whichAnInds[l]]]),max(annotations$Begin.Time..s.[anInd[whichAnInds[l]]],timeBins[k])))/annotations$Delta.Time..s.[anInd[whichAnInds[l]]]
          overlap = (min(timeBins[k+1],annotations$End.Time..s.[anInd[whichAnInds[l]]]) - max(annotations$Begin.Time..s.[anInd[whichAnInds[l]]],timeBins[k]))/annotations$Delta.Time..s.[anInd[whichAnInds[l]]]
          if (overlap >=minO){
            Alabels = c(Alabels,annotations$Tags[anInd[whichAnInds[l]]])
            SNR = c(SNR,annotations[anInd[whichAnInds[l]],SNRcol])
          }
        }
        names(Alabels) = seq_along(Alabels)
        names(SNR) = seq_along(SNR)
        if (length(Alabels)>0){
          temp1[k,'ALabel'][[1]] = list(Alabels)
          temp1[k,'SNR'][[1]] = list(SNR)
        }
      }
      
      # find any detections which overlap with this bin
      whichDetInds = which(detections$Begin.Time..s.[detInd]<timeBins[k+1] & detections$End.Time..s.[detInd]>timeBins[k]) 
      Dlabels = character()
      Dscores = numeric()
      if (length(whichDetInds)>0){
        for (l in 1:length(whichDetInds)){
          Dlabels = c(Dlabels,detections$Tags[detInd[whichDetInds[l]]])
          Dscores = c(Dscores,detections$Score[detInd[whichDetInds[l]]])
        }
        names(Dlabels) = seq_along(Dlabels)
        names(Dscores) = seq_along(Dscores)
        
        temp1[k,'DLabel'][[1]]= list(Dlabels)
        temp1[k,'DScore'][[1]] = list(Dscores)
      }
      
    }
    
    if (j==1){
      temp2 = temp1
    }else{
      temp2 = rbind(temp2,temp1)
    } 
  }
  
  if (i==1){
    allData = temp2
  }else{
    allData = rbind(allData,temp2)
  }
  
}

## Tally TP/FP/TN/FN across all classes and compute performance metrics
thresh = c(0.25,0.5,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.98,0.99)
metMat = matrix(nrow=length(thresh),ncol=11)
colnames(metMat) = c("Thresh",'nCalls','nTP','nFP','nTN','nFN','A','P','R','F1','FPR')

allALabels = stack(setNames(allData$ALabel,seq_along(allData$ALabel)))[2:1]
names(allALabels) = c('Index','Label')
if (!isempty(SNRcol)){
allALabels[,3] = unlist(allData$SNR)
names(allALabels) = c('Index','Label','SNR')}
allALabels = allALabels[!is.na(allALabels$Label),]
rownames(allALabels) = seq(length=dim(allALabels)[1])
allALabels$Index = as.numeric(allALabels$Index)

allDLabels = stack(setNames(allData$DLabel,seq_along(allData$DLabel)))[2:1]
allDLabels[,3] = unlist(allData$DScore)
names(allDLabels) = c('Index','Label','Score')
allDLabels = allDLabels[!is.na(allDLabels$Label),]
rownames(allDLabels) = seq(length=nrow(allDLabels))
allDLabels$Index = as.numeric(allDLabels$Index)

metMat[,2] = length(allALabels$Index)
allLabels = unlist(unique(c(allALabels$Label,allDLabels$Label)))

# confMat = matrix(nrow=length(allLabels),ncol=length(allLabels))
# rownames(confMat) = paste('true_',allLabels,sep="")
# colnames(confMat) = paste('pred_',allLabels,sep="")

for (j in 1:length(thresh)){
  
  goodInds = which(allDLabels$Score>=thresh[j]) # indices in allDLabels where detection scores exceeded confidence threshold
  binInds = allDLabels$Index[goodInds] # indices in allData where detection scores exceeded confidence threshold
  badInds = setdiff(1:length(allData$Time),binInds) # indices in allData with no detection score exceeding confidence threshold
  
  # Count up the easy FPs & TNs
  FP = length(intersect(binInds,which(is.na(allData$ALabel)))) # how many detections above the threshold occur in bins w NO annotation labels?
  TN = length(intersect(badInds,which(is.na(allData$ALabel)))) # how many times with NO detection above the threshold also have NO annotation labels?
  
  # Tally up instances of matched and mismatched labels in each bin with a detection and/or annotation
  TP = 0
  FN = 0
  for (i in 1:length(allLabels)){
    thisLabDet = allDLabels$Index[goodInds[which(allDLabels$Label[goodInds]==allLabels[i])]]  # indices in allData with this detection label and scores exceeding confidence threshold
    thisLabAn = allALabels$Index[allALabels$Label==allLabels[i]] # indices in allData with this annotation label
    otherLabAn = setdiff(allALabels$Index[allALabels$Label!=allLabels[i]],thisLabAn) # indices in allData with only some other annotation label
    TP = TP + length(which(thisLabDet %in% thisLabAn)) # count bins where detection label and annotation label were the same
    FN = FN + length(setdiff(thisLabAn,thisLabDet)) # count bins where there was an annotation with this label, but no detection with this label
    FP = FP + length(which(thisLabDet %in% otherLabAn)) # add to FPs instances with this detection label, but a different annotation label
    
    # confMat[i,i] = length(which(thisLabDet %in% thisLabAn))
    # if (length(allLabels)>=2){
    #   otherLabs = setdiff(seq_along(allLabels),i)
    #   for (k in otherLabs){
    #     otherLabAn = setdiff(allALabels$Index[allALabels$Label==allLabels[otherLabs]],thisLabAn) # indices in allData with only this other annotation label
    #   }
    # }
  }
  
  metMat[j,3] = TP
  metMat[j,4] = FP
  metMat[j,5] = TN
  metMat[j,6] = FN
  metMat[j,7] = round((TP+TN)/(TP+TN+FP+FN),3) #Accuracy
  metMat[j,8] = round(TP/(TP+FP),3) # Precision
  metMat[j,9] = round(TP/(TP+FN),3) # Recall
  metMat[j,10] = round((2*TP)/((2*TP)+FP+FN),3) # F1 Score
  metMat[j,11] = round(FP/(FP+TN),3) # FPR
  
}

metMat = as.data.frame(metMat)
metMat$Thresh = thresh
write.table(metMat,paste(saveDir,'/PerformanceMetrics.txt',sep=""))
cat(paste('Accuracy: ',as.character(min(metMat$A,na.rm=TRUE)*100),'-',as.character(max(metMat$A,na.rm=TRUE)*100),
          '\nPrecision: ',as.character(min(metMat$P,na.rm=TRUE)*100),'-',as.character(max(metMat$P,na.rm=TRUE)*100),
          '\nRecall: ',as.character(min(metMat$R,na.rm=TRUE)*100),'-',as.character(max(metMat$R,na.rm=TRUE)*100),
          '\nF1: ',as.character(min(metMat$F1,na.rm=TRUE)*100),'-',as.character(max(metMat$F1,na.rm=TRUE)*100),
          '\nFPR: ',as.character(min(metMat$FPR,na.rm=TRUE)*100),'-',as.character(max(metMat$FPR,na.rm=TRUE)*100),sep=""))

# Compute precision/recall as a function of SNR
if (!isempty(SNRcol)){
SNRvals = round(linspace(min(allALabels$SNR),max(allALabels$SNR),n=20),digits=2)
SNReval = matrix(nrow=length(SNRvals),ncol=8)
colnames(SNReval) = c('SNRthresh','N','TP','FP','TN','FN','P','R')

for (j in 1:length(SNRvals)){
  
  goodInds = which(allALabels$SNR>=SNRvals[j]) # indices in allALabels where annotation SNR exceededed threshold
  binInds = allALabels$Index[goodInds] # indices in allData where annotation SNR exceededed threshold
  badInds = setdiff(1:length(allData$Time),binInds) # indices in allData with no annotation SNR exceededing threshold
  
  # Count up the easy FNs & TNs
  FN = length(intersect(binInds,which(is.na(allData$DLabel)))) # how many annotations above the threshold occur in bins w NO detection labels?
  TN = length(intersect(badInds,which(is.na(allData$ALabel)))) # how many times with NO annotation above the threshold also have NO detection labels?
  
  # Tally up instances of matched and mismatched labels in each bin with a detection and/or annotation
  TP = 0
  FP = 0
  for (i in 1:length(allLabels)){
    thisLabDet = allDLabels$Index[allDLabels$Label==allLabels[i]]  # indices in allData with this detection label 
    thisLabAn = allALabels$Index[goodInds[which(allALabels$Label[goodInds]==allLabels[i])]] # indices in allData with this annotation label and SNR exceeding threshold
    otherLabAn = union(badInds,setdiff(allALabels$Index[allALabels$Label!=allLabels[i]],thisLabAn)) # indices in allData with only some other annotation label or SNR below threshold
    TP = TP + length(which(thisLabDet %in% thisLabAn)) # count bins where detection label and annotation label were the same
    FN = FN + length(setdiff(thisLabAn,thisLabDet)) # count bins where there was an annotation with this label and exceeding SNR thresh, but no detection with this label
    FP = FP + length(which(thisLabDet %in% otherLabAn)) # add to FPs instances with this detection label, but a different annotation label
    
    SNReval[j,2] = length(goodInds)
    SNReval[j,3] = TP
    SNReval[j,4] = FP
    SNReval[j,5] = TN
    SNReval[j,6] = FN
    SNReval[j,7] = round(TP/(TP+FP),3) # Precision
    SNReval[j,8] = round(TP/(TP+FN),3) # Recall
  }
  
}

SNReval = as.data.frame(SNReval)
SNReval$SNRthresh = SNRvals
write.table(SNReval,paste(saveDir,'/SNReval.txt',sep=""))
}

## Plot performance curves
# PR curve vs confidence score
ggplot(metMat,aes(label=Thresh))+
        geom_point(aes(x=R,y=P))+
        geom_path(aes(x=R,y=P))+
        geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
        coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
        labs(title=paste('PR Curve, Min Overlap = ',minO*100,'%',sep=""),
             x='Recall',
             y='Precision')
ggsave(filename=paste(saveDir,'/PR_conf.png',sep=""))

if (!isempty(SNRcol)){
# PR curve vs SNR
ggplot(SNReval,aes(label=SNRthresh))+
  geom_point(aes(x=R,y=P))+
  geom_path(aes(x=R,y=P))+
  geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title=paste('PR Curve, Min Overlap = ',minO*100,'%',sep=""),
       x='Recall',
       y='Precision')
ggsave(filename=paste(saveDir,'/PR_SNR.png',sep=""))
}

# ROC curve vs conf
ggplot(metMat,aes(label=Thresh))+
        geom_point(aes(x=FPR,y=R))+
        geom_path(aes(x=FPR,y=R))+
        geom_text(aes(x=FPR,y=R),hjust = 0, nudge_x = 0.0005)+
        coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
        labs(title=paste('ROC Curve, Min Overlap = ',minO*100,'%',sep=""),
             x='FPR',
             y='Recall')
ggsave(filename=paste(saveDir,'/ROC.png',sep=""))


# Plot precision vs thresh
ggplot(metMat)+
        geom_point(aes(x=Thresh,y=P))+
        geom_path(aes(x=Thresh,y=P))+
        coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
        labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
             x='Threshold',
             y='Precision')
ggsave(filename=paste(saveDir,'/PvThresh.png',sep=""))

# Plot recall vs thresh
ggplot(metMat)+
        geom_point(aes(x=Thresh,y=R))+
        geom_path(aes(x=Thresh,y=R))+
        coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
        labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
             x='Threshold',
             y='Recall')
ggsave(filename=paste(saveDir,'/RvThresh.png',sep=""))

# # Assumes 1 annotation file & 1 detection file for each test audio file
# # Expects time stamps in annotation & detection file names
# 
# anotDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/SelectionTables/20240730'
# detDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/Koogu/Testing/20240321-140450'
# listDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/Sounds'
# minO = 0.6 # minimum temp1oral overlap (% of annotation duration) to count a detection window as containing a manually annotated call
# # durs = 15 # duration of sound files (minutes)
# 
# ##### 
# annotation_files <- dir(anotDir,pattern='.txt')
# detection_files <- list.files(detDir, pattern = ".txt")
# list_file = list.files(listDir,pattern='_forwardSlash.txt')
# 
# # Find corresponding list file
# file_list <-  paste(listDir,list_file,sep = "/")
# audioPaths <- unlist(read.table(file_list, header = FALSE, sep = "\t"))
# 
# # Get start times of test audio file set
# fileStarts = parse_date_time(str_extract(audioPaths,"\\d{8}_\\d{6}"),'Ymd_HMS')
# 
# # # Get durations of files in list - ONLY WORKS FOR WAVE/MP3 FILES
# durs = numeric()
# for (i in 1:numel(audioPaths)){
#   durs = rbind(durs, (readWave(audioPaths[i],header=TRUE)$samples/readWave(audioPaths[i],header=TRUE)$sample.rate))}
# 
# allTimeBins = numeric()
# allData = numeric()
# 
# for (i in 1:length(annotation_files)) { # step through each test audio file and evaluate TP/FP/TN/FN
#   
#   file_annotation <- paste(anotDir, annotation_files[i], sep = "/")
#   annotations <-   read.table(file_annotation, header = TRUE, sep = "\t")
#   
#   # Remove annotation type 'waveform'
#   if ("View" %in% colnames(annotations) & "Spectrogram 1" %in% annotations$View & "Waveform 1" %in% annotations$View){
#   annotations <- subset(annotations, annotations$View == "Spectrogram 1")
#   }
#   
#   # Get file start timestamps to convert annotation begin times to absolute timestamps
#   # A_fileStarts = parse_date_time(str_extract(annotations$Begin.File,"\\d{8}_\\d{6}"),"Ymd_HMS")
#   A_fileStarts = parse_date_time(str_extract(file_annotation,"\\d{8}_\\d{6}"),"Ymd_HMS")
#   # annotations$StartTime = A_fileStarts + dseconds(annotations$File.Offset..s.)
#   # annotations$EndTime = annotations$StartTime + dseconds(annotations$Delta.Time..s.)
#   annotations$StartTime = A_fileStarts + dseconds(annotations$Begin.Time..s.)
#   annotations$EndTime = annotations$StartTime + dseconds(annotations$End.Time..s.-annotations$Begin.Time..s.)
#   
#   # Find corresponding detector output selection table and test audio file
#   thisTimeStamp = str_extract(annotation_files[i],'\\d{8}_\\d{6}')
#   matchingDets = str_which(detection_files,thisTimeStamp)
#   matchingAudio = str_which(audioPaths,str_extract(file_annotation,"\\d{8}_\\d{6}"))
#   
#   # Sort into ascending start time order, reset row numbers
#   annotations = annotations[order(annotations$StartTime),]
#   rownames(annotations) = 1:nrow(annotations)
#   
#   # Load detections
#   file_detec <-  paste(detDir,detection_files[matchingDets],sep = "/")
#   detections <- read.table(file_detec, header = TRUE, sep = "\t")
#   
#   # Remove detection type 'waveform'
#   if (("View" %in% colnames(detections) & "Spectrogram 1" %in% detections$View & "Waveform 1" %in% detections$View)){
#     detections <- subset(detections, detections$View == "Spectrogram 1")
#   }
#   
#   # Get file start timestamps to convert detection begin times to absolute timestamps
#   D_fileStarts = parse_date_time(str_extract(file_detec,"\\d{8}_\\d{6}"),"Ymd_HMS")
#   # detections$StartTime = D_fileStarts + dseconds(detections$File.Offset..s.)
#   # detections$EndTime = detections$StartTime + dseconds(detections$Delta.Time..s.)
#   detections$StartTime = D_fileStarts + dseconds(detections$Begin.Time..s.)
#   detections$EndTime = detections$StartTime + dseconds(detections$End.Time..s.-detections$Begin.Time..s.)
#   
#   # Sort into ascending start time order, reset row numbers
#   detections = detections[order(detections$StartTime),]
#   rownames(detections) = 1:nrow(detections)
#   
#   # Establish time bins consistent with how the detector saw the data 
#   seqFun <- Vectorize(function(x,y) seq.POSIXt(x,y,by=min(detections$End.Time..s.-detections$Begin.Time..s.)))
#   timeBins = as_datetime(c(seqFun(fileStarts[matchingAudio],(fileStarts[matchingAudio]+dseconds(durs[i])))))
#   allTimeBins = c(allTimeBins,timeBins)
#   
#   # Trim any annotations that exceed the last time bin 
#   # (there may be a tiny bit of data at the end of the file which is not seen by the detector because it's not long enough to be a full spectrogram)
#   tooLong = which(annotations$EndTime>timeBins[length(timeBins)])
#   annotations$EndTime[tooLong] = timeBins[length(timeBins)]
#   
#   # Determine full set of channels any annotations or detections exist in
#   allChans = sort(unique(c(annotations$Channel,detections$Channel)))
#   
#   temp2 = numeric()
#   
#   for (j in 1:length(allChans)){
#     
#     temp1 = data.frame(Time=timeBins)
#     temp1$Channel=j
#     temp1$ALabel=NA
#     temp1$DLabel=NA
#     temp1$DScore=NA
#     
#     anInd = which(annotations$Channel==allChans[j])
#     detInd = which(detections$Channel==allChans[j])
#     
#     # for each time bin, note existing annotation and/or detection labels, mark if TP/FP/TN/FN
#     for (k in 1:(length(timeBins)-1)){
#       
#       # find any annotations which overlap with this bin
#       whichAnInds = which(annotations$StartTime[anInd]<timeBins[k+1] & annotations$EndTime[anInd]>timeBins[k]) 
#       Alabels = character()
#       if (length(whichAnInds)>0){
#         for (l in 1:length(whichAnInds)){
#           # determine if overlap is sufficient
#           overlap = as.numeric(difftime(min(timeBins[k+1],annotations$EndTime[anInd[whichAnInds[l]]]),max(annotations$StartTime[anInd[whichAnInds[l]]],timeBins[k])))/annotations$Delta.Time..s.[anInd[whichAnInds[l]]]
#           if (overlap >=minO){
#             Alabels = c(Alabels,annotations$Tags[anInd[whichAnInds[l]]])
#           }
#         }
#         names(Alabels) = seq_along(Alabels)
#         if (length(Alabels)>0){
#           temp1[k,'ALabel'][[1]] = (Alabels)
#         }
#       }
#       
#       # find any detections which overlap with this bin
#       whichDetInds = which(detections$StartTime[detInd]<timeBins[k+1] & detections$EndTime[detInd]>timeBins[k]) 
#       Dlabels = character()
#       Dscores = numeric()
#       if (length(whichDetInds)>0){
#         for (l in 1:length(whichDetInds)){
#           Dlabels = c(Dlabels,detections$Tags[detInd[whichDetInds[l]]])
#           Dscores = c(Dscores,detections$Score[detInd[whichDetInds[l]]])
#         }
#         names(Dlabels) = seq_along(Dlabels)
#         names(Dscores) = seq_along(Dscores)
#         
#         temp1[k,'DLabel'][[1]]= list(Dlabels)
#         temp1[k,'DScore'][[1]] = list(Dscores)
#       }
#       
#     }
#     
#     if (j==1){
#       temp2 = temp1
#     }else{
#       temp2 = rbind(temp2,temp1)
#     } 
#   }
#   
#   if (i==1){
#     allData = temp2
#   }else{
#     allData = rbind(allData,temp2)
#   }
# 
# }
# 
# ## Calculate performance metrics, confusion matrix
# thresh = c(0.25,0.5,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.98,0.99)
# metMat = matrix(nrow=length(thresh),ncol=11)
# colnames(metMat) = c("Thresh",'nCalls','nTP','nFP','nTN','nFN','A','P','R','F1','FPR')
# 
# allALabels = stack(setNames(allData$ALabel,seq_along(allData$ALabel)))[2:1]
# names(allALabels) = c('Index','Label')
# allALabels = allALabels[!is.na(allALabels$Label),]
# rownames(allALabels) = seq(length=nrow(allALabels))
# allALabels$Index = as.numeric(allALabels$Index)
# allDLabels = stack(setNames(allData$DLabel,seq_along(allData$DLabel)))[2:1]
# allDLabels[,3] = unlist(allData$DScore)
# names(allDLabels) = c('Index','Label','Score')
# allDLabels = allDLabels[!is.na(allDLabels$Label),]
# rownames(allDLabels) = seq(length=nrow(allDLabels))
# allDLabels$Index = as.numeric(allDLabels$Index)
# 
# metMat[,2] = length(allALabels$Index)
# allLabels = unlist(unique(c(allALabels$Label,allDLabels$Label)))
# 
# # confMat = matrix(nrow=length(allLabels),ncol=length(allLabels))
# # rownames(confMat) = paste('true_',allLabels,sep="")
# # colnames(confMat) = paste('pred_',allLabels,sep="")
# 
# for (j in 1:length(thresh)){
#   
#   goodInds = which(allDLabels$Score>=thresh[j]) # indices in allDLabels where detection scores exceeded confidence threshold
#   binInds = allDLabels$Index[goodInds] # indices in allData where detection scores exceeded confidence threshold
#   badInds = setdiff(1:length(allData$Time),binInds) # indices in allData with no detection score exceeding confidence threshold
# 
#   # Count up the easy FPs & TNs
#   FP = length(intersect(binInds,which(is.na(allData$ALabel)))) # how many detections occur in bins w no annotation labels?
#   TN = length(intersect(which(is.na(allData$ALabel)),badInds)) # how many times with NO detection also have no annotation labels?
#   
#   # Tally up instances of matched and mismatched labels in each bin with a detection and/or annotation
#   TP = 0
#   FN = 0
#   for (i in 1:length(allLabels)){
#     thisLabDet = allDLabels$Index[goodInds[which(allDLabels$Label[goodInds]==allLabels[i])]]  # indices in allData with this detection label and scores exceeding confidence threshold
#     thisLabAn = allALabels$Index[allALabels$Label==allLabels[i]] # indices in allData with this annotation label
#     otherLabAn = setdiff(allALabels$Index[allALabels$Label!=allLabels[i]],thisLabAn) # indices in allData with only some other annotation label
#     TP = TP + length(which(thisLabDet %in% thisLabAn)) # count bins where detection label and annotation label were the same
#     FN = FN + length(setdiff(thisLabAn,thisLabDet)) # count bins where there was an annotation with this label, but no detection with this label
#     FP = FP + length(which(thisLabDet %in% otherLabAn)) # add to FPs instances with this detection label, but a different annotation label
#     
#     # confMat[i,i] = length(which(thisLabDet %in% thisLabAn))
#     # if (length(allLabels)>=2){
#     #   otherLabs = setdiff(seq_along(allLabels),i)
#     #   for (k in otherLabs){
#     #     otherLabAn = setdiff(allALabels$Index[allALabels$Label==allLabels[otherLabs]],thisLabAn) # indices in allData with only this other annotation label
#     #   }
#     # }
#   }
#   
#   metMat[j,3] = TP
#   metMat[j,4] = FP
#   metMat[j,5] = TN
#   metMat[j,6] = FN
#   metMat[j,7] = round((TP+TN)/(TP+TN+FP+FN),3) #Accuracy
#   metMat[j,8] = round(TP/(TP+FP),3) # Precision
#   metMat[j,9] = round(TP/(TP+FN),3) # Recall
#   metMat[j,10] = round((2*TP)/((2*TP)+FP+FN),3) # F1 Score
#   metMat[j,11] = round(FP/(FP+TN),3) # FPR
#   
# }
# 
# metMat = as.data.frame(metMat)
# metMat$Thresh = thresh
# cat(paste('Accuracy: ',as.character(min(metMat$A)*100),'-',as.character(max(metMat$A)*100),
#           '\nPrecision: ',as.character(min(metMat$P)*100),'-',as.character(max(metMat$P)*100),
#           '\nRecall: ',as.character(min(metMat$R)*100),'-',as.character(max(metMat$R)*100),
#           '\nF1: ',as.character(min(metMat$F1)*100),'-',as.character(max(metMat$F1)*100),
#           '\nFPR: ',as.character(min(metMat$FPR)*100),'-',as.character(max(metMat$FPR)*100),sep=""))
# 
# 
# labInd = which(sapply(allData$DLabel,FUN=function(x) 'Clack Thunder' %in% x))
# 
# 
# ## Plot Precision-Recall and ROC curves
# # PR curve
# print(ggplot(metMat,aes(label=Thresh))+
#         geom_point(aes(x=R,y=P))+
#         geom_path(aes(x=R,y=P))+
#         geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('PR Curve, Min Overlap = ',minO*100,'%',sep=""),
#              x='Recall',
#              y='Precision'))
# 
# # ROC curve
# print(ggplot(metMat,aes(label=Thresh))+
#         geom_point(aes(x=FPR,y=R))+
#         geom_path(aes(x=FPR,y=R))+
#         geom_text(aes(x=FPR,y=R),hjust = 0, nudge_x = 0.0005)+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('ROC Curve, Min Overlap = ',minO*100,'%',sep=""),
#              x='FPR',
#              y='Recall'))
# 
# # Plot precision vs thresh
# print(ggplot(metMat)+
#         geom_point(aes(x=Thresh,y=P))+
#         geom_path(aes(x=Thresh,y=P))+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
#              x='Threshold',
#              y='Precision'))
# 
# # Plot recall vs thresh
# print(ggplot(metMat)+
#         geom_point(aes(x=Thresh,y=R))+
#         geom_path(aes(x=Thresh,y=R))+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
#              x='Threshold',
#              y='Recall'))


########
# anotDir = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/PrevAnalysis_ManualDetections'
# detDir = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/KooguNARW_DetectorOutput'
# listDir = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/prevAnalysis_ListFiles'
# # minCount = 0.1 # minimum duration (s) of annotation window which must overlap w detector window to count as a bin containing a call
# minO = 0.8 # minimum temp1oral overlap (% of annotation duration) to count a detection window as containing a manually annotated call
# durs = 15 # duration of sound files (minutes)
# 
# #####
# annotation_files <- dir(anotDir,pattern='.txt')
# detection_files <- list.files(detDir, pattern = "allChans.txt")
# list_files = list.files(listDir,pattern='.txt')
# perfDF = matrix(nrow=length(annotation_files),ncol=7)
# dep = list()
# 
# for (i in 1:length(annotation_files)) {
# 
#   TP = 0
#   FN = 0
#   FP = 0
#   TN = 0
# 
#   file_annotation <- paste(anotDir, annotation_files[i], sep = "/")
#   d = str_sub(annotation_files[i],1,4)
#   dep = rbind(dep,d)
#   annotations <-   read.table(file_annotation, header = TRUE, sep = "\t")
# 
#   # Remove annotation type 'waveform'
#   annotations <- subset(annotations, annotations$View == "Spectrogram 1")
# 
#   # Get file start timestamps to convert annotation begin times to absolute timestamps
#   fileStarts = parse_date_time(str_extract(annotations$Begin.File,"\\d{8}_\\d{6}"),"Ymd_HMS")
#   annotations$StartTime = fileStarts + dseconds(annotations$File.Offset..s.)
#   annotations$EndTime = annotations$StartTime + dseconds(annotations$Delta.Time..s.)
# 
#   # Sort into ascending start time order, reset row numbers
#   annotations = annotations[order(annotations$StartTime),]
#   rownames(annotations) = 1:nrow(annotations)
# 
#   # Find corresponding detector output selection table
#   # matchingFile = str_which(detection_files,annotation_files[i])
#   matchingFile = str_which(detection_files,d)
#   file_detec <-  paste(detDir,detection_files[matchingFile],sep = "/")
#   detections <- read.table(file_detec, header = TRUE, sep = "\t")
# 
#   # Remove detection type 'waveform'
#   detections <- subset(detections, detections$View == "Spectrogram 1")
# 
#   # Get file start timestamps to convert detection begin times to absolute timestamps
#   fileStarts = parse_date_time(str_extract(detections$Begin.File,"\\d{8}_\\d{6}"),"Ymd_HMS")
#   detections$StartTime = fileStarts + dseconds(detections$File.Offset..s.)
#   detections$EndTime = detections$StartTime + dseconds(detections$Delta.Time..s.)
# 
#   # Sort into ascending start time order, reset row numbers
#   detections = detections[order(detections$StartTime),]
#   rownames(detections) = 1:nrow(detections)
# 
#   # Find corresponding list file
#   matchingFile = str_which(list_files,d)
#   file_list <-  paste(listDir,list_files[matchingFile],sep = "/")
#   fileList <- read.table(file_list, header = FALSE, sep = "\t")
# 
#   # Get start times of file set
#   # dataStart = parse_date_time(str_extract(fileList[1,],"\\d{8}_\\d{6}"),'Ymd_HMS')
# 
#   # # Establish time bins consistent with how the detector saw the data
#   # timeBins = seq.POSIXt(from=dataStart,to=max(detections$EndTime,annotations$EndTime),by=detections$Delta.Time..s.[1])
# 
#   # Get start times of file set
#   fileStarts = parse_date_time(str_extract(fileList[,1],"\\d{8}_\\d{6}"),'Ymd_HMS')
# 
#   # # Get durations of files in list - ONLY WORKS FOR WAVE/MP3 FILES
#   # durs = soundgen::analyze(x=fileList)$summary$duration
# 
#   # Establish time bins consistent with how the detector saw the data
#   seqFun <- Vectorize(function(x,y) seq.POSIXt(x,y,by=detections$Delta.Time..s.[1]))
#   timeBins = as_datetime(c(seqFun(fileStarts,(fileStarts+dseconds((durs*60)-1)))))
# 
#   # # Chop annotations into short time bins to see which overlap with detections
#   # seqFun <- Vectorize(function(x,y) seq(x,y,by=0.1))
#   # smallSteps = seqFun(annotations$StartTime,annotations$EndTime)
# 
#   # Determine full set of channels any annotations or detections exist in
#   allChans = sort(unique(c(annotations$Channel,detections$Channel)))
#   # allChans = c(1,2,7,10)
# 
#   for (j in 1:length(allChans)){
# 
#     temp11 = data.frame()
#     temp12 = data.frame()
# 
#     # Find annotations in this channel
#     anInd = which(annotations$Channel==allChans[j])
#     # if (length(anInd)>0){
#     #   whichBinsA = histc(unlist(smallSteps[anInd]),as.numeric(timeBins))
#     #   binswAnnots = timeBins[which(whichBinsA$cnt>=(minCount/0.1))]
#     # } else {
#     #   whichBinsA = numeric()
#     #   binswAnnots = numeric()
#     # }
# 
#     if (length(anInd)>0){
#       whichBinsA = numeric()
#       for (k in 1:length(anInd)){
#         whichStart = which(timeBins<=annotations$StartTime[anInd[k]]) # which bin contains annotation start time
#         whichStart = whichStart[length(whichStart)]
#         whichEnd = which(timeBins>=annotations$EndTime[anInd[k]]) # which bin contains annotation end time
#         whichEnd = whichEnd[1]
# 
#         if (whichEnd-whichStart==1){ # if annotation falls within a single bin
#           whichBinsA = c(whichBinsA,whichStart)
#         } else if (whichEnd-whichStart>1){ # if annotation spans more than one bin, only count bins that have enough of the annotation in them
#           if (timeBins[whichStart+1]-annotations$StartTime[anInd[k]] >= minO*annotations$Delta.Time..s.[anInd[k]]){
#             whichBinsA = c(whichBinsA,whichStart)
#           }
#           if (whichEnd-whichStart>2){
#             whichBinsA = c(whichBinsA,seq(whichStart+1,whichEnd-2,by=1))
#           }
#           if (annotations$EndTime[anInd[k]]-timeBins[whichEnd-1] >= minO*annotations$Delta.Time..s.[anInd[k]]){
#             whichBinsA = c(whichBinsA,whichEnd-1)
#           }
#         }
#       }
#       binswAnnots = timeBins[whichBinsA]
#     } else {
#       whichBinsA = numeric()
#       binswAnnots = numeric()
#     }
# 
#     # In case of multiple annotations occurring in the same bin, remove duplicate bin counts
#     binswAnnots = unique(binswAnnots)
# 
#     # Determine which bins have detections in this channel
#     detInd = which(detections$Channel==allChans[j])
#     # whichBinsD = histc(as.numeric(detections$StartTime[detInd]),as.numeric(timeBins))
#     # binswDets = unique(timeBins[whichBinsD$bin])
#     whichBinsD = timeBins %in% detections$StartTime[detInd]
#     binswDets = timeBins[whichBinsD]
# 
#     sharedCalls = intersect(binswAnnots,binswDets)
#     missedCalls = setdiff(binswAnnots,as_datetime(sharedCalls))
#     badDets = setdiff(binswDets,sharedCalls)
#     trueAbs = setdiff(timeBins,c(sharedCalls,missedCalls,badDets))
# 
#     TP = TP + length(sharedCalls)
#     FN = FN + length(missedCalls)
#     FP = FP + length(badDets)
#     TN = TN + length(trueAbs)
# 
#     TPind = which(detections$StartTime[detInd] %in% sharedCalls)
#     FPind = which(detections$StartTime[detInd] %in% badDets)
# 
#     if (length(detInd)>0){ # start collecting information for this channel
#       temp11 = data.frame(TimeBin = binswDets) # time bin
#       temp11$Channel = allChans[j]             # channel
#       temp11$True = 0                          # Was there a call truly present?
#       temp11$True[TPind] = 1
#       temp11$Detected = 1                      # Was a call detected?
#       temp11$Conf = detections$Score[detInd]   # confidence scores of detections
#     }
#     if (length(missedCalls)>0){
#       temp12 = data.frame(TimeBin = as_datetime(missedCalls))
#       temp12$Channel = allChans[j]
#       temp12$True = 1
#       temp12$Detected = 0
#       temp12$Conf = NA
#     }
# 
#     if (j==1){
#       eval(parse(text=paste(d,' = rbind(temp11,temp12)',sep="")))
#     } else {
#       eval(parse(text=paste(d,' = do.call(\"rbind\",list(',d,',temp11,temp12))',sep="")))
#     }
# 
#   }
# 
#   perfDF[i,2] = TP
#   perfDF[i,3] = FP
#   perfDF[i,4] = FN
#   perfDF[i,5] = TN
# 
# }
# 
# eval(parse(text=paste(d,'= arrange(',d,',TimeBin,Channel)',sep="")))
# 
# perfDF[,6] = round(perfDF[,2]/(perfDF[,2]+perfDF[,3]),2)
# perfDF[,7] = round(perfDF[,2]/(perfDF[,2]+perfDF[,4]),2)
# perfDF = as.data.frame(perfDF)
# colnames(perfDF) = c('Dep','TP','FP','FN','TN','Precision','Recall')
# perfDF$Dep = dep
# 
# 
# ## Make Precision-Recall & ROC curves for each deployment
# 
# thresh = c(0.25,0.5,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.98,0.99)
# # thresh = seq(0,1,by=0.1)
# deps = list(MD01,MD02,MD03,MD04)
# depNames = c('MD01','MD02','MD03','MD04')
# 
# for (i in 1:length(deps)){
# 
#   metMat = matrix(nrow=length(thresh),ncol=3)
#   colnames(metMat) = c('P','R','FPR')
# 
#   for (j in 1:length(thresh)){
# 
#     goodInds = which(deps[[i]]$Conf>=thresh[j])
#     badInds = setdiff(1:length(deps[[i]]$TimeBin),goodInds)
#     TP = length(which(deps[[i]]$True[goodInds]==1 & deps[[i]]$Detected[goodInds]==1))
#     FP = length(which(deps[[i]]$True[goodInds]==0 & deps[[i]]$Detected[goodInds]==1))
#     FN = length(which(deps[[i]]$True[badInds]==1))
#     TN = perfDF$TN[i] + length(which(deps[[i]]$Detected[badInds]==1))
# 
#     metMat[j,1] = TP/(TP+FP) # Precision
#     metMat[j,2] = TP/(TP+FN) # Recall
#     metMat[j,3] = FP/(FP+TN) # FPR
# 
#   }
# 
#   metMat = as.data.frame(metMat)
#   metMat$Thresh = thresh
# 
#   # PR curve
#   print(ggplot(metMat,aes(label=Thresh))+
#           geom_point(aes(x=R,y=P))+
#           geom_path(aes(x=R,y=P))+
#           geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
#           coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#           labs(title=paste(depNames[i],' PR Curve, Min Overlap = ',minO*100,'%',sep=""),
#                x='Recall',
#                y='Precision'))
# 
#   # ROC curve
#   print(ggplot(metMat,aes(label=Thresh))+
#           geom_point(aes(x=FPR,y=R))+
#           geom_path(aes(x=FPR,y=R))+
#           geom_text(aes(x=FPR,y=R),hjust = 0, nudge_x = 0.0005)+
#           coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#           labs(title=paste(depNames[i],' ROC Curve, Min Overlap = ',minO*100,'%',sep=""),
#                x='FPR',
#                y='Recall'))
# 
#   # Plot precision vs thresh
#   print(ggplot(metMat)+
#           geom_point(aes(x=Thresh,y=P))+
#           geom_path(aes(x=Thresh,y=P))+
#           coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#           labs(title=paste(depNames[i],' , Min Overlap = ',minO*100,'%',sep=""),
#                x='Threshold',
#                y='Precision'))
# 
#   # Plot recall vs thresh
#   print(ggplot(metMat)+
#           geom_point(aes(x=Thresh,y=R))+
#           geom_path(aes(x=Thresh,y=R))+
#           coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#           labs(title=paste(depNames[i],' , Min Overlap = ',minO*100,'%',sep=""),
#                x='Threshold',
#                y='Recall'))
# }
# 
# 
# 
# ## Make a single Precision-Recall curve for all test data
# 
# allDat = rbind(MD01,MD02,MD03,MD04)
# thresh = c(0.25,0.5,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.98,0.99)
# metMat = matrix(nrow=length(thresh),ncol=3)
# colnames(metMat) = c('P','R','FPR')
# 
# for (j in 1:length(thresh)){
# 
#   goodInds = which(allDat$Conf>=thresh[j])
#   badInds = setdiff(1:length(allDat$TimeBin),goodInds)
#   TP = length(which(allDat$True[goodInds]==1 & allDat$Detected[goodInds]==1))
#   FP = length(which(allDat$True[goodInds]==0 & allDat$Detected[goodInds]==1))
#   FN = length(which(allDat$True[badInds]==1))
#   TN = perfDF$TN[i] + length(which(allDat$Detected[badInds]==1))
# 
#   metMat[j,1] = TP/(TP+FP) # Precision
#   metMat[j,2] = TP/(TP+FN) # Recall
#   metMat[j,3] = FP/(FP+TN) # FPR
# 
# }
# 
# metMat = as.data.frame(metMat)
# metMat$Thresh = thresh
# 
# # PR curve
# print(ggplot(metMat,aes(label=Thresh))+
#         geom_point(aes(x=R,y=P))+
#         geom_path(aes(x=R,y=P))+
#         geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('PR Curve, Min Overlap = ',minO*100,'%',sep=""),
#              x='Recall',
#              y='Precision'))
# 
# # ROC curve
# print(ggplot(metMat,aes(label=Thresh))+
#         geom_point(aes(x=FPR,y=R))+
#         geom_path(aes(x=FPR,y=R))+
#         geom_text(aes(x=FPR,y=R),hjust = 0, nudge_x = 0.0005)+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('ROC Curve, Min Overlap = ',minO*100,'%',sep=""),
#              x='FPR',
#              y='Recall'))
# 
# # Plot precision vs thresh
# print(ggplot(metMat)+
#         geom_point(aes(x=Thresh,y=P))+
#         geom_path(aes(x=Thresh,y=P))+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
#              x='Threshold',
#              y='Precision'))
# 
# # Plot recall vs thresh
# print(ggplot(metMat)+
#         geom_point(aes(x=Thresh,y=R))+
#         geom_path(aes(x=Thresh,y=R))+
#         coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
#         labs(title=paste('Min Overlap = ',minO*100,'%',sep=""),
#              x='Threshold',
#              y='Recall'))
