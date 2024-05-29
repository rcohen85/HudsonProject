library(pracma)
library(stringr)
library(ggplot2)

inDir = '/Users/rec297/Documents/CCB/DCLDE/NoiseProp'
numThreshVals = 13
fileList = list.files(path=inDir,pattern='.txt')
noiseProp = c('4x','3x','Balanced','Balanced','4x','3x')
# SNR = c('All','All','15','All','15','15')
# Fs = c('Original','Shifted','Original','Shifted','Original','Shifted')
# tempRes = c('Original, 30s Win','Original','Shifted','Original, 30s Win',
#             'Original','Shifted','Original, 30s Win','Original','Shifted')
# dataset = c('MD','CC','MD','DCLDE','DCLDE','CC') # SNR comp
dataset = c('MD','MD','MD','DCLDE','DCLDE','DCLDE') # noise prop
# dataset = c('MD','CC','CC','MD','DCLDE','DCLDE') # Fs shift
# dataset = c('MD','MD','CC','CC','CC','MD','DCLDE','DCLDE','DCLDE') # temp res

PRvals = matrix(data=NA,nrow=numThreshVals*numel(fileList),ncol=5)
j=1
for (i in 1:numel(fileList)){
  tb = read.table(paste(inDir,'/',fileList[i],sep=""))
  PRvals[j:(j+numThreshVals-1),1] = tb$P
  PRvals[j:(j+numThreshVals-1),2] = tb$R
  PRvals[j:(j+numThreshVals-1),3] = str_replace(fileList[i],'.txt','')
  # PRvals[j:(j+numThreshVals-1),4] = rep(SNR[i],numThreshVals)
  PRvals[j:(j+numThreshVals-1),4] = rep(noiseProp[i],numThreshVals)
  # PRvals[j:(j+numThreshVals-1),4] = rep(Fs[i],numThreshVals)
  # PRvals[j:(j+numThreshVals-1),4] = rep(tempRes[i],numThreshVals)
  PRvals[j:(j+numThreshVals-1),5] = rep(dataset[i],numThreshVals)
  j=j+numThreshVals
}

PRvals = as.data.frame(PRvals)
colnames(PRvals) = c('P','R','Model','noiseProp','Dataset')
PRvals$P = as.numeric(PRvals$P)
PRvals$R = as.numeric(PRvals$R)
PRvals$Model = as.factor(PRvals$Model)
# PRvals$SNR = as.factor(PRvals$SNR)
PRvals$noiseProp = as.factor(PRvals$noiseProp)
# PRvals$Fs = as.factor(PRvals$Fs)
# PRvals$tempRes = as.factor(PRvals$tempRes)
PRvals$Dataset = as.factor(PRvals$Dataset)

print(ggplot(PRvals)+
        geom_point(aes(x=R,y=P,shape=Dataset))+
        geom_path(aes(x=R,y=P,group=Model,colour=noiseProp))+
        coord_cartesian(xlim=c(0,1),ylim=c(0,0.35))+
        labs(x='Recall',
             y='Precision'))
