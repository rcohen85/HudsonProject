library(stringr)
# library(pracma)
library(dplyr)
library(lubridate)
# library(soundgen)
library(ggplot2)

anotDir = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/PrevAnalysis_ManualDetections'
detDir = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/KooguNARW_DetectorOutput'
listDir = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/prevAnalysis_ListFiles'
# minCount = 0.1 # minimum duration (s) of annotation window which must overlap w detector window to count as a bin containing a call
minO = 0.1 # minimum temporal overlap (% of annotation duration) to count a detection window as containing a manually annotated call
durs = 15 # duration of sound files (minutes)

##### 
annotation_files <- dir(anotDir,pattern='.txt')
detection_files <- list.files(detDir, pattern = "allChans.txt")
list_files = list.files(listDir,pattern='.txt')
perfDF = matrix(nrow=length(annotation_files),ncol=7)
dep = list()

for (i in 1:length(annotation_files)) {
  
  TP = 0
  FN = 0
  FP = 0
  TN = 0
  
  file_annotation <- paste(anotDir, annotation_files[i], sep = "/")
  d = str_sub(annotation_files[i],1,4)
  dep = rbind(dep,d)
  annotations <-   read.table(file_annotation, header = TRUE, sep = "\t")
  
  # Remove annotation type 'waveform'
  annotations <- subset(annotations, annotations$View == "Spectrogram 1")
  
  # Get file start timestamps to convert annotation begin times to absolute timestamps
  fileStarts = parse_date_time(str_extract(annotations$Begin.File,"\\d{8}_\\d{6}"),"Ymd_HMS")
  annotations$StartTime = fileStarts + dseconds(annotations$File.Offset..s.)
  annotations$EndTime = annotations$StartTime + dseconds(annotations$Delta.Time..s.)
  
  # Sort into ascending start time order, reset row numbers
  annotations = annotations[order(annotations$StartTime),]
  rownames(annotations) = 1:nrow(annotations)
  
  # Find corresponding detector output selection table
  matchingFile = str_which(detection_files,d)
  file_detec <-  paste(detDir,detection_files[matchingFile],sep = "/")
  detections <- read.table(file_detec, header = TRUE, sep = "\t")
  
  # Remove detection type 'waveform'
  detections <- subset(detections, detections$View == "Spectrogram 1")
  
  # Get file start timestamps to convert detection begin times to absolute timestamps
  fileStarts = parse_date_time(str_extract(detections$Begin.File,"\\d{8}_\\d{6}"),"Ymd_HMS")
  detections$StartTime = fileStarts + dseconds(detections$File.Offset..s.)
  detections$EndTime = detections$StartTime + dseconds(detections$Delta.Time..s.)
  
  # Sort into ascending start time order, reset row numbers
  detections = detections[order(detections$StartTime),]
  rownames(detections) = 1:nrow(detections)
  
  # Find corresponding list file
  matchingFile = str_which(list_files,d)
  file_list <-  paste(listDir,list_files[matchingFile],sep = "/")
  fileList <- read.table(file_list, header = FALSE, sep = "\t")
  
  # Get start times of file set
  # dataStart = parse_date_time(str_extract(fileList[1,],"\\d{8}_\\d{6}"),'Ymd_HMS')

  # # Establish time bins consistent with how the detector saw the data
  # timeBins = seq.POSIXt(from=dataStart,to=max(detections$EndTime,annotations$EndTime),by=detections$Delta.Time..s.[1])

  # Get start times of file set
  fileStarts = parse_date_time(str_extract(fileList[,1],"\\d{8}_\\d{6}"),'Ymd_HMS')
  
  # # Get durations of files in list - ONLY WORKS FOR WAVE/MP3 FILES
  # durs = soundgen::analyze(x=list(fileList))$summary$duration
  
  # Establish time bins consistent with how the detector saw the data (assumed to be bin)
  seqFun <- Vectorize(function(x,y) seq.POSIXt(x,y,by=detections$Delta.Time..s.[1]))
  timeBins = as_datetime(c(seqFun(fileStarts,(fileStarts+dseconds((durs*60)-1)))))
  
  # # Chop annotations into short time bins to see which overlap with detections
  # seqFun <- Vectorize(function(x,y) seq(x,y,by=0.1))
  # smallSteps = seqFun(annotations$StartTime,annotations$EndTime)
  
  # Determine full set of channels any annotations or detections exist in
  allChans = sort(unique(c(annotations$Channel,detections$Channel)))
  # allChans = c(1,2,7,10)
  
  for (j in 1:length(allChans)){
    
    temp1 = data.frame()
    temp2 = data.frame()
    
    # Find annotations in this channel
    anInd = which(annotations$Channel==allChans[j])
    # if (length(anInd)>0){
    #   whichBinsA = histc(unlist(smallSteps[anInd]),as.numeric(timeBins))
    #   binswAnnots = timeBins[which(whichBinsA$cnt>=(minCount/0.1))]
    # } else {
    #   whichBinsA = numeric()
    #   binswAnnots = numeric()
    # }
    
    if (length(anInd)>0){
      whichBinsA = numeric()
      for (k in 1:length(anInd)){
        whichStart = which(timeBins<=annotations$StartTime[anInd[k]]) # which bin contains annotation start time
        whichStart = whichStart[length(whichStart)]
        whichEnd = which(timeBins>=annotations$EndTime[anInd[k]]) # which bin contains annotation end time
        whichEnd = whichEnd[1]
        
        if (whichEnd-whichStart==1){ # if annotation falls within a single bin
          whichBinsA = c(whichBinsA,whichStart)
        } else if (whichEnd-whichStart>1){ # if annotation spans more than one bin, only count bins that have enough of the annotation in them
          if (timeBins[whichStart+1]-annotations$StartTime[anInd[k]] >= minO*annotations$Delta.Time..s.[anInd[k]]){
            whichBinsA = c(whichBinsA,whichStart)
          }
          if (whichEnd-whichStart>2){
            whichBinsA = c(whichBinsA,seq(whichStart+1,whichEnd-2,by=1))
          }
          if (annotations$EndTime[anInd[k]]-timeBins[whichEnd-1] >= minO*annotations$Delta.Time..s.[anInd[k]]){
            whichBinsA = c(whichBinsA,whichEnd-1)
          }
        }
      }
      binswAnnots = timeBins[whichBinsA]
    } else {
      whichBinsA = numeric()
      binswAnnots = numeric()
    }
    
    # In case of multiple annotations occurring in the same bin, remove duplicate bin counts
    binswAnnots = unique(binswAnnots)
    
    # Determine which bins have detections in this channel
    # detInd = which(detections$Channel==allChans[j])
    # whichBinsD = histc(as.numeric(detections$StartTime[detInd]),as.numeric(timeBins))
    # binswDets = unique(timeBins[whichBinsD$bin])
    whichBinsD = timeBins %in% detections$StartTime[detInd]
    binswDets = timeBins[whichBinsD]
    
    sharedCalls = intersect(binswAnnots,binswDets)
    missedCalls = setdiff(binswAnnots,as_datetime(sharedCalls))
    badDets = setdiff(binswDets,sharedCalls)
    trueAbs = setdiff(timeBins,c(sharedCalls,missedCalls,badDets))
    
    TP = TP + length(sharedCalls)
    FN = FN + length(missedCalls)
    FP = FP + length(badDets)
    TN = TN + length(trueAbs)
    
    TPind = which(detections$StartTime[detInd] %in% sharedCalls)
    FPind = which(detections$StartTime[detInd] %in% badDets)
    
    if (length(detInd)>0){ # start collecting information for this channel
      temp1 = data.frame(TimeBin = binswDets) # time bin
      temp1$Channel = allChans[j]             # channel
      temp1$True = 0                          # Was there a call truly present?
      temp1$True[TPind] = 1
      temp1$Detected = 1                      # Was a call detected?
      temp1$Conf = detections$Score[detInd]   # confidence scores of detections
    }
    if (length(missedCalls)>0){
      temp2 = data.frame(TimeBin = as_datetime(missedCalls))
      temp2$Channel = allChans[j]
      temp2$True = 1
      temp2$Detected = 0
      temp2$Conf = NA
    }
    
    if (j==1){
      eval(parse(text=paste(d,' = rbind(temp1,temp2)',sep="")))
    } else {
      eval(parse(text=paste(d,' = do.call(\"rbind\",list(',d,',temp1,temp2))',sep="")))
    }
    
  }
  
  perfDF[i,2] = TP
  perfDF[i,3] = FP
  perfDF[i,4] = FN
  perfDF[i,5] = TN
  
}

eval(parse(text=paste(d,'= arrange(',d,',TimeBin,Channel)',sep="")))

perfDF[,6] = round(perfDF[,2]/(perfDF[,2]+perfDF[,3]),2)
perfDF[,7] = round(perfDF[,2]/(perfDF[,2]+perfDF[,4]),2)
perfDF = as.data.frame(perfDF)
colnames(perfDF) = c('Dep','TP','FP','FN','TN','Precision','Recall')
perfDF$Dep = dep


## Make Precision-Recall & ROC curves

thresh = c(0.25,0.5,0.65,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.98,0.99)
# thresh = seq(0,1,by=0.1)
deps = list(MD01,MD02,MD03,MD04)
depNames = c('MD01','MD02','MD03','MD04')

for (i in 1:length(deps)){
  
  metMat = matrix(nrow=length(thresh),ncol=3)
  colnames(metMat) = c('P','R','FPR')
  
  for (j in 1:length(thresh)){
    
    goodInds = which(deps[[i]]$Conf>=thresh[j])
    badInds = setdiff(1:length(deps[[i]]$TimeBin),goodInds)
    TP = length(which(deps[[i]]$True[goodInds]==1 & deps[[i]]$Detected[goodInds]==1))
    FP = length(which(deps[[i]]$True[goodInds]==0 & deps[[i]]$Detected[goodInds]==1))
    FN = length(which(deps[[i]]$True[badInds]==1))
    TN = perfDF$TN[i] + length(which(deps[[i]]$Detected[badInds]==1))
    
    metMat[j,1] = TP/(TP+FP) # Precision
    metMat[j,2] = TP/(TP+FN) # Recall
    metMat[j,3] = FP/(FP+TN) # FPR
    
  }
  
  metMat = as.data.frame(metMat)
  metMat$Thresh = thresh
  
  # PR curve
  print(ggplot(metMat,aes(label=Thresh))+
    geom_point(aes(x=R,y=P))+
    geom_path(aes(x=R,y=P))+
    geom_text(aes(x=R,y=P),hjust = 0, nudge_x = 0.0005)+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste(depNames[i],' PR Curve, Min Overlap = ',minO*100,'%',sep=""),
         x='Recall',
         y='Precision'))
  
  # ROC curve
  print(ggplot(metMat,aes(label=Thresh))+
    geom_point(aes(x=FPR,y=R))+
    geom_path(aes(x=FPR,y=R))+
    geom_text(aes(x=FPR,y=R),hjust = 0, nudge_x = 0.0005)+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste(depNames[i],' ROC Curve, Min Overlap = ',minO*100,'%',sep=""),
         x='FPR',
         y='Recall'))
  
  # Plot precision vs thresh
  print(ggplot(metMat)+
    geom_point(aes(x=Thresh,y=P))+
    geom_path(aes(x=Thresh,y=P))+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste(depNames[i],' , Min Overlap = ',minO*100,'%',sep=""),
         x='Threshold',
         y='Precision'))
  
  # Plot recall vs thresh
  print(ggplot(metMat)+
    geom_point(aes(x=Thresh,y=R))+
    geom_path(aes(x=Thresh,y=R))+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
    labs(title=paste(depNames[i],' , Min Overlap = ',minO*100,'%',sep=""),
         x='Threshold',
         y='Recall'))
}


