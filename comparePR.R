library(pracma)
library(stringr)
library(ggplot2)

inDir = 'P:/users/cohen_rebecca_rec297/CCB/DCLDE2024/ComparisonFigures/TrainSetSize'
numThreshVals = 13
fileList = list.files(path=inDir,pattern='.txt')

# noiseProp = c('1:1','4:1','3:1','1:1','1:1','3:1','2:1','2:1','2:1','3:1','4:1')
trainSize = c('800','400','16k','400','400')
# SNR = c('0','0','15','15','15','0')
# Fs = c('Original','Sped Up 10x','Original','Sped Up 10x','Original','Sped Up 10x')
# tempRes = c('Original', '30s Win Eval','Sped Up 10x','Original', '30s Win Eval',
#             'Sped Up 10x','Original', '30s Win Eval','Sped Up 10x')
# mod = c('SNR18, 10x','10x','SNR15, 10x','None','None','None')
# mod = c('BirdNET','BirdNET','BirdNET','Shiu et al.','Shiu et al.','Shiu et al.')

# dataset = c('MD','CC','DCLDE') # baseline
# dataset = c('MD','MD','MD','CC','DCLDE','CC','CC','DCLDE','MD','DCLDE','DCLDE') # noise prop
dataset = c('MD','DCLDE','CC','MD','CC') # train set size
# dataset = c('MD','CC','MD','DCLDE','CC','DCLDE') # SNR comp
# dataset = c('MD','CC','CC','MD','DCLDE','DCLDE') # Fs shift
# dataset = c('MD','MD','CC','CC','CC','MD','DCLDE','DCLDE','DCLDE') # temp res
# dataset = c('CC','MD','DCLDE','DCLDE','CC','MD')

PRvals = matrix(data=NA,nrow=numThreshVals*numel(fileList),ncol=5)
ROCvals = matrix(data=NA,nrow=numThreshVals*numel(fileList),ncol=5)
j=1
for (i in 1:numel(fileList)){
  tb = read.table(paste(inDir,'/',fileList[i],sep=""))
  PRvals[j:(j+numThreshVals-1),1] = tb$P
  PRvals[j:(j+numThreshVals-1),2] = tb$R
  PRvals[j:(j+numThreshVals-1),3] = str_replace(fileList[i],'.txt','')
  # PRvals[j:(j+numThreshVals-1),4] = rep(noiseProp[i],numThreshVals)
  PRvals[j:(j+numThreshVals-1),4] = rep(trainSize[i],numThreshVals)
  # PRvals[j:(j+numThreshVals-1),4] = rep(SNR[i],numThreshVals)
  # PRvals[j:(j+numThreshVals-1),4] = rep(Fs[i],numThreshVals)
  # PRvals[j:(j+numThreshVals-1),4] = rep(tempRes[i],numThreshVals)
  # PRvals[j:(j+numThreshVals-1),4] = rep(mod[i],numThreshVals)
  PRvals[j:(j+numThreshVals-1),5] = rep(dataset[i],numThreshVals)
  
  ROCvals[j:(j+numThreshVals-1),1] = tb$R
  ROCvals[j:(j+numThreshVals-1),2] = tb$FPR
  ROCvals[j:(j+numThreshVals-1),3] = str_replace(fileList[i],'.txt','')
  # ROCvals[j:(j+numThreshVals-1),4] = rep(noiseProp[i],numThreshVals)
  ROCvals[j:(j+numThreshVals-1),4] = rep(trainSize[i],numThreshVals)
  # ROCvals[j:(j+numThreshVals-1),4] = rep(SNR[i],numThreshVals)
  # ROCvals[j:(j+numThreshVals-1),4] = rep(Fs[i],numThreshVals)
  # ROCvals[j:(j+numThreshVals-1),4] = rep(tempRes[i],numThreshVals)
  # ROCvals[j:(j+numThreshVals-1),4] = rep(mod[i],numThreshVals)
  ROCvals[j:(j+numThreshVals-1),5] = rep(dataset[i],numThreshVals)
  
  j=j+numThreshVals
}

PRvals = as.data.frame(PRvals)
colnames(PRvals) = c('P','R','Model','trainSize','Dataset')
PRvals$P = as.numeric(PRvals$P)
PRvals$R = as.numeric(PRvals$R)
PRvals$Model = as.factor(PRvals$Model)
# PRvals$noiseProp = as.factor(PRvals$noiseProp)
PRvals$trainSize = as.factor(PRvals$trainSize)
# PRvals$SNR = as.factor(PRvals$SNR)
# PRvals$Fs = as.factor(PRvals$Fs)
# PRvals$tempRes = as.factor(PRvals$tempRes)
# PRvals$mod = as.factor(PRvals$mod)
PRvals$Dataset = as.factor(PRvals$Dataset)

ROCvals = as.data.frame(ROCvals)
colnames(ROCvals) = c('R','FPR','Model','trainSize','Dataset')
ROCvals$R = as.numeric(ROCvals$R)
ROCvals$FPR = as.numeric(ROCvals$FPR)
ROCvals$Model = as.factor(ROCvals$Model)
# ROCvals$noiseProp = as.factor(ROCvals$noiseProp)
ROCvals$trainSize = as.factor(ROCvals$trainSize)
# ROCvals$SNR = as.factor(ROCvals$SNR)
# ROCvals$Fs = as.factor(ROCvals$Fs)
# ROCvals$tempRes = as.factor(ROCvals$tempRes)
# ROCvals$mod = as.factor(ROCvals$mod)
ROCvals$Dataset = as.factor(ROCvals$Dataset)

print(ggplot(PRvals)+
        geom_path(aes(x=R,y=P,group=Model,colour=trainSize),linewidth=1.05)+
        # scale_color_manual(values=c('#00BA38','#F8766D','#619CFF'))+
        geom_point(aes(x=R,y=P,shape=Dataset),size=1.4)+
        coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
        labs(x='Recall',
             y='Precision')+
        guides(color=guide_legend(title='Tuning Set Size')))

print(ggplot(ROCvals)+
        geom_path(aes(x=FPR,y=R,group=Model,colour=trainSize),linewidth=1.05)+
        geom_point(aes(x=FPR,y=R,shape=Dataset),size=1.4)+
        coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
        labs(x='False Positive Rate',
             y='Recall')+
        guides(color=guide_legend(title='Tuning Set Size')))
