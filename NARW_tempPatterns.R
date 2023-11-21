library(stringr)
library(lubridate)
library(pracma)
library(ggplot2)

MD01 = read.table('U:/projects/2013_UnivMD_Maryland_71485/KooguNARWDetEval/NewAnnotationTables/MD01_ground_truth.txt', sep="\t",header=TRUE,check.names=FALSE)
fileStarts = parse_date_time(str_extract(MD01$`Begin File`,"\\d{8}_\\d{6}"),'Ymd_HMS')
callStarts = fileStarts + dseconds(MD01$`File Offset (s)`)
timeBins = as_datetime(seq.Date(from=as.Date('2014-11-23'),to=as.Date('2015-04-30'),by=1))
whichBins = histc(as.numeric(callStarts),as.numeric(timeBins))

plotDF = data.frame(Bins=timeBins,Counts=whichBins$cnt)

# Plot daily call counts
ggplot(plotDF)+geom_col(aes(x=Bins,y=Counts))

hourBins = seq(1,24,1)
callHours = data.frame(Hour=hour(callStarts))

# Plot hourly call counts
ggplot(callHours)+geom_histogram(aes(x=Hour),binwidth=1)
