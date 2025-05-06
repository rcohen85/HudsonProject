library(gtools)
library(ggplot2)
library(stringr)


inDir = 'C:/Users/rec297/Box/Hawthorne Valley/Eval Selection Tables'
pTP = 0.88
selTabs = list.files(inDir,pattern='*.txt',recursive=TRUE)

for (i in 1:length(selTabs)){
  
  evalTab = read.table(paste(inDir,selTabs[i],sep='/'), header = TRUE, sep = "\t",stringsAsFactors = TRUE,comment.char="",fill=TRUE)
  tRow = which(evalTab$Eval==1)
  fRow = which(evalTab$Eval==0)
  # evalTab = evalTab[c(tRow,fRow),c('Begin.Date.Time','Begin.File','Eval','Score','Label')]
  evalTab = evalTab[c(tRow,fRow),c('Eval','Score','Label')]
  
  if (i==1){
    bigTab = evalTab
  } else {
    bigTab = rbind(bigTab,evalTab)
  }
  
  
  # TPind = which(evalTab$TP.FP == 'TP')
  # FPind = which(evalTab$TP.FP == 'FP')
  # noEval = which(evalTab$TP.FP == '')
  # 
  # evalTab$TP.FP[TPind] = 1
  # evalTab$TP.FP[FPind] = 0
  # evalTab = evalTab[-noEval,]
  # evalTab$TP.FP = as.numeric(evalTab$TP.FP)
  # evalTab$Score[evalTab$Score==1] = 0.995
  # evalTab$Score = logit(evalTab$Score,min=0,max=1)
  
}

# rowdate = as.data.frame(paste(as.character(bigTab$Begin.Date.Time),as.character(bigTab$Begin.File)))
# uniqueSels = unique(rowdate)
# if (nrow(uniqueSels)==nrow(rowdate)){
  row.names(bigTab) = seq(1,nrow(bigTab))
  bigTab$Eval = droplevels(as.factor(bigTab$Eval))
  bigTab$Eval = as.numeric(bigTab$Eval)
  bigTab$Eval = bigTab$Eval - 1
  bigTab$LogitScore = log(bigTab$Score/(1-bigTab$Score))
  bigTab$Label = as.character(bigTab$Label)
  species = unique(bigTab$Label)
  
  for (i in 1:length(species)){
    thisSpec = str_which(bigTab$Label,species[i])
    
    mod = glm(Eval~Score,data=bigTab[thisSpec,],family="binomial")
    
    
    predData = data.frame(Score=seq(0,1,0.001))
    # predData = data.frame(LogitScore=seq(-3,7,0.1))
    predData$Preds = predict(mod,predData,type='r')
    thresh = (log(pTP/(1-pTP))-mod$coefficients[1])/mod$coefficients[2]
    
    plot(Eval~Score,bigTab[thisSpec,],
         main=paste(species[i],', N Eval = ',length(thisSpec),', ',pTP*100,'% TP Score Threshold = ',round(thresh,2),sep=""),
         ylab='pr(BirdNET prediction is correct)',xlab='Confidence Score',
         xlim=range(predData$Score),pch=16,cex=1.5,col=rgb(0,0,0,.2))
    lines(predData$Preds~predData$Score,lwd=2,col=rgb(0,0.75,1,0.5))
    abline(v=thresh,col='red',lwd=2)
    
    # print(ggplot(bigTab[thisSpec,], aes(x=Score, y=Eval)) + geom_point() +
    #   stat_smooth(method="glm", color="blue",se=FALSE,
    #               method.args = list(family=binomial))+
    #     geom_vline(xintercept=thresh,color="red",linewidth=1)+
    #   labs(title=paste(species[i],', N Eval = ',length(thisSpec),', 90% TP Score Threshold = ',round(thresh,2),sep="")))
    
    
  }
# } else {
#   print('Duplicate evaluations present')
# }
