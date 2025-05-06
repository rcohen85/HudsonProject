# Script to compare ground truth selection tables against detector output selection tables and
# compute various detector performance evaluation metrics
# Note: Assumes ground truth selection files & detection selection files each correspond to a single test audio file

# Settings -------------------
# Path to ground truth selection tables
anotDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/AtlStur/ThunderClassifier_Training/KooguModel/SecondIteration/TestData/SelectionTables'
# name of column containing ground truth labels
labCol = 'Tags' 
# Path to audio test files
audioDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/AtlStur/ThunderClassifier_Training/KooguModel/SecondIteration/TestData/Audio'
# Audio file extension
fileExt = ".wav"
# Path to detector output selection tables
detDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/AtlStur/ThunderClassifier_Training/KooguModel/SecondIteration/ModelEval/20241211-142047/SelTabs'
# name of column containing detection labels
detCol = "Tags" 
# Path to save performance metrics & plots
saveDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/AtlStur/ThunderClassifier_Training/KooguModel/SecondIteration/ModelEval/20241211-142047'

# duration of detector bins (s); if -1, will estimate detector bin size from detector output
binSize = 2
# minimum temp1oral overlap (% of annotation duration) to count a detection window as containing a manually annotated call
minO = 0.5
# either specify channel(s) of interest, or leave as an empty list to analyze all available channels
chans = 1 

# Calculations ------------------------- 
library(stringr)
library(tuneR)
library(ggplot2)
library(pracma)

annotation_files <- dir(anotDir,pattern='.txt')
detection_files <- list.files(detDir, pattern = ".txt")
audio_files <- dir(audioDir,pattern=fileExt,recursive=TRUE)

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
    annotations <-   read.table(file_annotation, header = TRUE, sep = "\t",check.names = FALSE)
    
    if (dim(annotations)[1]>0){
      # Remove annotation type 'waveform'
      if ("View" %in% colnames(annotations) & "Spectrogram 1" %in% annotations$View & "Waveform 1" %in% annotations$View){
        annotations <- subset(annotations, annotations$View == "Spectrogram 1")
      }
      if (!("Delta Time (s)" %in% colnames(annotations))){
        annotations$"Delta Time (s)" = annotations$"End Time (s)" - annotations$"Begin Time (s)"
      }
      
      # Sort into ascending start time order, reset row numbers
      annotations = annotations[order(annotations$"Begin Time (s)"),]
      rownames(annotations) = 1:nrow(annotations)
    }
  } else {
    annotations = data.frame()
  }
  
  if (!isempty(matchingDets)){
    # Load detections
    file_detec <-  paste(detDir,detection_files[matchingDets],sep = "/")
    detections <- read.table(file_detec, header = TRUE, sep = "\t",check.names = FALSE)
    
    if (dim(detections)[1]>0){
      # Remove detection type 'waveform'
      if (("View" %in% colnames(detections) & "Spectrogram 1" %in% detections$View & "Waveform 1" %in% detections$View)){
        detections <- subset(detections, detections$View == "Spectrogram 1")
      } # add Delta Time column if not already present
      if (!("Delta Time (s)" %in% colnames(detections))){
        detections$"Delta Time (s)" = detections$"End Time (s)" - detections$"Begin Time (s)"
      }
      if (!("Score" %in% colnames(detections))){
        detections$Score = detections$Confidence
      }
      noCallInd = str_which(detections[[detCol]],'nocall')
      if (!isempty(noCallInd)){
      detections = detections[-noCallInd,]}
      
      if (nrow(detections)>0){
      # Sort into ascending start time order, reset row numbers
      detections = detections[order(detections$"Begin Time (s)"),]
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
    timeBins = seq(0,durs[i],by=round(min(detections$"End Time (s)"-detections$"Begin Time (s)"),digits=2))}
  allTimeBins = c(allTimeBins,timeBins)
  
  # Trim any annotations that exceed the last time bin 
  # (there may be a tiny bit of data at the end of the file which is not seen by the detector because it's not long enough to be a full spectrogram)
  tooLong = which(annotations$"End Time (s)">timeBins[length(timeBins)])
  annotations$"End Time (s)"[tooLong] = timeBins[length(timeBins)]
  
  # Determine full set of channels any annotations or detections exist in
  if (isempty(chans)){
  allChans = sort(unique(c(annotations$Channel,detections$Channel)))
  } else { allChans = chans}
  
  temp2 = numeric()
  
  for (j in 1:length(allChans)){
    
    temp1 = data.frame(Time=timeBins)
    temp1$File=audio_files[i]
    temp1$Channel=j
    temp1$ALabel=NA
    temp1$DLabel=NA
    temp1$DScore=NA
    
    anInd = which(annotations$Channel==allChans[j])
    detInd = which(detections$Channel==allChans[j])
    
    # for each time bin, note existing annotation and/or detection labels, mark if TP/FP/TN/FN
    for (k in 1:(length(timeBins)-1)){
      
      # find any annotations which overlap with this bin
      whichAnInds = which(annotations$"Begin Time (s)"[anInd]<timeBins[k+1] & annotations$"End Time (s)"[anInd]>timeBins[k]) 
      Alabels = character()
      if (length(whichAnInds)>0){
        for (l in 1:length(whichAnInds)){
          # determine if overlap is sufficient
          overlap = (min(timeBins[k+1],annotations$"End Time (s)"[anInd[whichAnInds[l]]]) - max(annotations$"Begin Time (s)"[anInd[whichAnInds[l]]],timeBins[k]))/annotations$"Delta Time (s)"[anInd[whichAnInds[l]]]
          if (overlap >=minO){
            Alabels = c(Alabels,annotations[anInd[whichAnInds[l]],labCol])
          }
        }
        names(Alabels) = seq_along(Alabels)
        if (length(Alabels)>0){
          temp1[k,'ALabel'][[1]] = list(Alabels)
        }
      }
      
      # find any detections which overlap with this bin
      whichDetInds = which(detections$"Begin Time (s)"[detInd]<timeBins[k+1] & detections$"End Time (s)"[detInd]>timeBins[k]) 
      Dlabels = character()
      Dscores = numeric()
      if (length(whichDetInds)>0){
        for (l in 1:length(whichDetInds)){
          Dlabels = c(Dlabels,detections[detInd[whichDetInds[l]],detCol])
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

# Tally TP/FP/TN/FN across all classes and compute performance metrics
thresh = c(0.1,0.15,0.25,0.5,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.98,0.99)
metMat = matrix(nrow=length(thresh),ncol=11)
colnames(metMat) = c("Thresh",'nCalls','nTP','nFP','nTN','nFN','A','P','R','F1','FPR')

allALabels = stack(setNames(allData$ALabel,seq_along(allData$ALabel)))[2:1]
names(allALabels) = c('Index','Label')

allALabels = allALabels[!is.na(allALabels$Label),]
rownames(allALabels) = seq(length=dim(allALabels)[1])
allALabels$Index = as.numeric(allALabels$Index)

allDLabels = stack(setNames(allData$DLabel,seq_along(allData$DLabel)))[2:1]
allDLabels[,3] = unlist(allData$DScore)
names(allDLabels) = c('Index','Label','Score')
allDLabels = allDLabels[!is.na(allDLabels$Label),]
rownames(allDLabels) = seq(length=nrow(allDLabels))
allDLabels$Index = as.numeric(allDLabels$Index)

allLabels = unlist(unique(c(allALabels$Label,allDLabels$Label)))

# Tally up instances of matched and mismatched labels in each bin with a detection and/or annotation

for (i in 1:length(allLabels)){
  
  metMat[,2] = length(which(allALabels$Label==allLabels[i]))
  
  for (j in 1:length(thresh)){
    
    goodInds = which(allDLabels$Score>=thresh[j]) # indices in allDLabels where detection scores exceeded confidence threshold
    binInds = allDLabels$Index[goodInds] # indices in allData where detection scores exceeded confidence threshold
    badInds = setdiff(1:length(allData$Time),binInds) # indices in allData with no detection score exceeding confidence threshold
    
    thisDetLab = allDLabels$Index[goodInds[which(allDLabels$Label[goodInds]==allLabels[i])]]  # indices in allData with this detection label and scores exceeding confidence threshold
    noDetLab = setdiff(seq(1,nrow(allData)),thisDetLab) # indices in allData with NO detection of this label & exceeding confidence thresh
    thisAnLab = allALabels$Index[allALabels$Label==allLabels[i]] # indices in allData with this annotation label
  
    notThisAnLab = setdiff(seq(1,nrow(allData)),thisAnLab) # indices in allData WITHOUT this annotation label
    TP = length(which(thisDetLab %in% thisAnLab)) # count bins where detection label and annotation label were the same
    FN = length(setdiff(thisAnLab,thisDetLab)) # count bins where there was an annotation with this label, but no detection with this label
  
    FP = length(intersect(thisDetLab,notThisAnLab)) # count bins with this detection label above threshold, but no/other annotation label
    TN = length(intersect(noDetLab,notThisAnLab)) # count bins with NO detection above threshold and also NO annotation for this label
    
    
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
  write.table(metMat,paste(saveDir,'/PerformanceMetrics_',str_remove(allLabels[i],' '),'.txt',sep=""))
  cat(paste('Label: ',allLabels[i],
            '\nAccuracy: ',as.character(min(metMat$A,na.rm=TRUE)*100),'-',as.character(max(metMat$A,na.rm=TRUE)*100),
            '\nPrecision: ',as.character(min(metMat$P,na.rm=TRUE)*100),'-',as.character(max(metMat$P,na.rm=TRUE)*100),
            '\nRecall: ',as.character(min(metMat$R,na.rm=TRUE)*100),'-',as.character(max(metMat$R,na.rm=TRUE)*100),
            '\nF1: ',as.character(min(metMat$F1,na.rm=TRUE)*100),'-',as.character(max(metMat$F1,na.rm=TRUE)*100),
            '\nFPR: ',as.character(min(metMat$FPR,na.rm=TRUE)*100),'-',as.character(max(metMat$FPR,na.rm=TRUE)*100),'\n',sep=""))
  
}


# Plot performance curves ---------------------
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

