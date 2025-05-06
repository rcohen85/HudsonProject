library(lubridate)
library(stringr)
library(ggplot2)
library(dplyr)

dateRange = c(as_datetime('2023/07/08'),as_datetime('2023/07/25'))

# Water Chestnut sensor
inDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/ChestnutSAVSoundscapes/Dissolved_Oxygen/DO.Temp.2023.waterchestnut'
fileList = dir(inDir,pattern='*.txt')

for (i in 1:length(fileList)){
DO = read.table(paste(inDir,fileList[i],sep="/"),
                skip=2,header=TRUE,fill=TRUE,sep=',',check.names=FALSE)
DO[,'Time (sec)'] = as_datetime(DO[,'Time (sec)'])
if (i==1){
  trapaData = DO
} else {
  trapaData = rbind(trapaData,DO)
}
}

colnames(trapaData) = c('Datetime','BatteryVolts','Temp','DO_mgl','Q')

# dr = which(trapaData$Datetime>=dateRange[1] & trapaData$Datetime<=dateRange[2])
# trapaData = trapaData[dr,]

# SAV sensor

inFile = "W:/projects/2022_NOAA-NERRS_HudsonNY_144488/ChestnutSAVSoundscapes/Dissolved_Oxygen/DO.Temp.SAV River location/226037.csv"

SAVdata = read.table(inFile,header=TRUE,fill=TRUE,sep=',',check.names=FALSE)
SAVdata$DateTimeStamp = parse_date_time(SAVdata$DateTimeStamp,'mdy HM')
NPdata = str_which(SAVdata$StationCode,'hudnpwq')
SAVdata = SAVdata[NPdata,]

plotDF = data.frame(trapaData$Datetime,trapaData$DO_mgl,rep('Chestnut',nrow(trapaData)))
colnames(plotDF) = c('Date','DO','Site')
plotDF2 = data.frame(SAVdata$DateTimeStamp,SAVdata$DO_mgl,rep('SAV1',nrow(SAVdata)))
colnames(plotDF2) = c('Date','DO','Site')
plotDF = rbind(plotDF,plotDF2)

# Plot
minDO = min(plotDF$DO[plotDF$Date>=dateRange[1] & plotDF$Date<=dateRange[2]])
maxDO = max(plotDF$DO[plotDF$Date>=dateRange[1] & plotDF$Date<=dateRange[2]])

ggplot(plotDF,aes(x=Date,y=DO,colour=Site))+
  geom_point()+
  scale_color_manual(labels=c('WC1','SAV1'),values=c('blue2','chartreuse3'))+
  coord_cartesian(xlim=c(dateRange[1],dateRange[2]),
                  ylim=c(minDO,maxDO*1.1))+
  xlab("")+
  ylab("Dissolved Oxygen (mg/l)")+ 
  theme(axis.title.y=element_text(size=16),
        axis.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position = c(0.9, 0.85))+
  guides(color = guide_legend(reverse=TRUE))


ggplot(plotDF,aes(x=Site,y=DO))+
  geom_boxplot(notch=TRUE)+
  xlab("")+
  ylab("Dissolved Oxygen (mg/l)")


plotDF %>% group_by(Site) %>% summarize(med=median(DO,na.rm=TRUE))

