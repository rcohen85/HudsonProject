library(mgcv)
library(gtools)
library(ggplot2)

selTab = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/AtlStur/ThunderClassifier_Output/AS04/AS04_Model_20240124_5.txt'

evalTab = read.table(selTab, header = TRUE, sep = "\t")

TPind = which(evalTab$TP.FP == 'TP')
FPind = which(evalTab$TP.FP == 'FP')
noEval = which(evalTab$TP.FP == '')

evalTab$TP.FP[TPind] = 1
evalTab$TP.FP[FPind] = 0
evalTab = evalTab[-noEval,]
evalTab$TP.FP = as.numeric(evalTab$TP.FP)
# evalTab$Score[evalTab$Score==1] = 0.995
# evalTab$Score = logit(evalTab$Score,min=0,max=1)

mod = gam(TP.FP~s(Score,bs='ts',k=4),data=evalTab,family=binomial())
# predDat = data.frame(Score=logit(seq(0.005,0.995,by=0.005),min=0,max=1))
predDat = data.frame(Score=seq(0.005,0.995,by=0.005))
predDat$FP = predict.gam(mod,predDat,type="response")

ggplot(predDat,aes(x=Score,y=FP))+geom_path()
