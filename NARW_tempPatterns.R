library(stringr)
library(lubridate)
library(pracma)
library(ggplot2)
library(ggbreak)
library(suncalc)

MD01 = read.table('U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/NewAnnotationTables/Combined_ground_truth.txt', sep="\t",header=TRUE,check.names=FALSE)
stDate = '2014-11-23'
endDate = '2016-07-14'

fileStarts = parse_date_time(str_extract(MD01$`Begin File`,"\\d{8}_\\d{6}"),'Ymd_HMS')
callStarts = fileStarts + dseconds(MD01$`File Offset (s)`)
callStarts_local = with_tz(callStarts,tzone="America/New_York")
timeBins = as_datetime(seq.Date(from=as.Date(stDate),to=as.Date(endDate),by=7))
whichBins = histc(as.numeric(callStarts_local),as.numeric(timeBins))

plotDF = data.frame(Bins=timeBins,Counts=whichBins$cnt)

# Plot daily call counts
print(ggplot(plotDF,aes(x=Bins,y=Counts))+
  geom_col(fill=4,alpha=0.6)+
  scale_x_continuous(breaks=c(as_datetime(c('2014-11-15','2014-12-15','2015-01-15','2015-02-15',
                                            '2015-03-15','2015-04-15','2016-02-15','2016-03-15',
                                            '2016-04-15','2016-05-15','2016-06-15','2016-07-15'))),
                     limits=c(as_datetime(c('2014-11-03','2016-07-15'))),
                     name="")+
  labs(title="Manually-Annotated NARW Upcalls")+
    # theme(axis.text.x=element_text(angle=45))+
    scale_x_break(c(as_datetime('2015-04-19'),as_datetime('2016-02-01')),scales='free',space=0.2)+
    theme_light())

# Convert call start times to normalized time of day
dayData = getSunlightTimes(date=seq.Date(as.Date(stDate)-1,as.Date(endDate)+1,by=1),
                           lat=38.303000,lon=-74.505000,
                           keep=c("nauticalDawn","sunrise","sunset","nauticalDusk"),
                           tz="America/New_York")

sunrise = as.numeric(dayData$sunrise)
sunset = as.numeric(dayData$sunset)
sunrise_next = as.numeric(c(dayData$sunrise[2:nrow(dayData)],dayData$sunrise[nrow(dayData)]))

dayDurs = (sunset - sunrise)/60
nightDurs = (sunrise_next - sunset)/60

normTime = data.frame(ToD = rep(NA,length(callStarts_local)))
for (i in 1:length(callStarts_local)){
  thisCall = as.numeric(callStarts_local[i])
  if (sum(thisCall >= sunrise & thisCall < sunset)){ # call falls in the daytime
    dayInd = which(thisCall >= sunrise & thisCall < sunset)
    normTime$ToD[i] = -((sunset[dayInd]-thisCall)/60)/dayDurs[dayInd]
  } else if (sum(thisCall >= sunset & thisCall < sunrise_next)){ # call falls in the nighttime
    nightInd = which(thisCall >= sunset & thisCall < sunrise_next)
    normTime$ToD[i] = ((thisCall-sunset[nightInd])/60)/nightDurs[nightInd]
  }
}
# callHours = data.frame(Hour=hour(callStarts_local))

# Plot hourly call counts
# ggplot(callHours)+geom_histogram(aes(x=Hour),bins=24)
print(ggplot(normTime)+
  geom_histogram(aes(x=ToD),bins=24,color='deepskyblue3',fill=4,alpha=0.6)+
  scale_x_continuous(breaks=c(-1,0,1),label=c("Sunrise","Sunset","Sunrise"))+
  labs(y="Counts",x=NULL,title="Normalized Time of Day")+
  theme_light())


