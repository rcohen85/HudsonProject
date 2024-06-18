library(stringr)

inPath = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/NewAnnotationTables/MD04_ground_truth.txt'
outDir = 'U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/NewClassifierTrainTest/TrainAnnotations'


bigTable = read.table(inPath,sep="\t",check.names=FALSE,header=TRUE)
bigTable$`Begin Time (s)` = bigTable$`File Offset (s)`
bigTable$`End Time (s)` = bigTable$`Begin Time (s)` + bigTable$`Delta Time (s)`
bigTable$Tags = 'NARW'
trainAudio = unique(bigTable$`Begin File`)

for (i in 1:length(trainAudio)){
  
  fileInd = str_which(bigTable$`Begin File`,trainAudio[i])
  newTable = bigTable[fileInd,]
  saveName = str_replace(trainAudio[i],'.aif','.txt')
  write.table(newTable,paste(outDir,'/',saveName,sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  
}