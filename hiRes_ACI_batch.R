setwd("C:/Users/rec297/Documents/GitHub/HudsonProject")
library(soundecology)
library(seewave)
library(parallel)
library(tuneR)
library(stringr)
library(lubridate)
source("ACI_RC.R")

dep = c('Swamp_Forest')
inDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/HawthorneValley/FEP_MVB_2022'
outDir = 'W:/projects/2022_NOAA-NERRS_HudsonNY_144488/HawthorneValley/FEP_MVB_2022/AcousticIndices'
soundindex = "acoustic_complexity"
min_freq = 2000
max_freq = 15000
no_cores=1
dur=60
j=1
fft_w=1024
flac=TRUE
fileheader <- c("FILENAME,SAMPLINGRATE,BIT,DURATION,START_S,CHANNELS,INDEX,FFT_W,MIN_FREQ,MAX_FREQ,J,LEFT_CHANNEL,RIGHT_CHANNEL")

#Open flac files
get_wav <- function(directory, flacsoundfile){
  #is it windows?
  if (.Platform$OS.type == "windows"){
    from_file = paste(directory, "\\", flacsoundfile, sep = "")
  }else{
    from_file = paste(directory, "/", flacsoundfile, sep = "")
  }
  
  wav_file = paste(strtrim(flacsoundfile, nchar(flacsoundfile) - 5), "wav", sep = ".")
  
  file.copy(from = from_file, to = flacsoundfile)
  wav2flac(flacsoundfile, reverse = TRUE, overwrite = TRUE)
  # file.remove(flacsoundfile)
  if (file.exists(flacsoundfile)) file.remove(flacsoundfile)
  if (file.exists(wav_file)){
    return(wav_file)
  }else{
    return(NA)
  }
}

for (i in 1:length(dep)){
  
  #Start timer
  time0 <- proc.time()
  
  directory=paste(inDir,'/',dep[i],sep="")
  resultfile = paste(outDir,'/',dep[i],'_ACI_',dur,'s.csv',sep="")
  cat(fileheader, file = resultfile, append = FALSE)
  
  #Are the files flac
  if (flac == TRUE){
    wav_files <- dir(path = directory, pattern = "flac$", ignore.case = TRUE)
    if (length(wav_files) == 0) {
      stop(paste("Could not find any .flac files in the specified directory:\n    ", directory))
    }
  }else{
    wav_files <- dir(path = directory, pattern = "wav$", ignore.case = TRUE)
    if (length(wav_files) == 0) {
      stop(paste("Could not find any .wav files in the specified directory:\n    ", directory))
    }
  }
  
  
  for (l in 1:length(wav_files)){
    
    if (flac == TRUE){
      soundfile_path <- get_wav(directory, wav_files[l])
    }else{
      if (.Platform$OS.type == "windows"){
        soundfile_path = paste(directory, "\\", wav_files[l], sep="")
      }else{
        soundfile_path = paste(directory, "/", wav_files[l], sep="")
      }
    }
    
    nSamp = length(readWave(soundfile_path)@left)
    Fs = readWave(soundfile_path)@samp.rate
    if (readWave(soundfile_path)@stereo){
      no_channels = 2
    } else {no_channels = 1}
    stepSize = dur*Fs
    chunks = floor(nSamp/stepSize)
    from = 1
    units = "samples"
    fileStart = as_datetime(str_extract(wav_files[l],'\\d{8}_\\d{6}'))
    return_list = list()
    durs = list()
    timeStamps = list()
    start_s = list()
    
    for(k in 1:floor(chunks)){
      start_s = c(start_s,((k-1)*dur)+1)
      to = from+stepSize-1
      this_soundfile <- readWave(soundfile_path, from = from, to = to, units = units)
      durs = rbind(durs,round(length(this_soundfile@left)/this_soundfile@samp.rate, 2))
      timeStamps = rbind(timeStamps,as.data.frame(fileStart+((k-1)*dur)+(seq(1,floor(dur/j)))))
      temp_list <- ACI_RC(this_soundfile, min_freq = min_freq, max_freq = max_freq, j = j, fft_w = fft_w)
      from = from + stepSize
      
      if(k==1){
        return_list = temp_list
      }else{
        return_list = mapply(rbind,return_list,temp_list,SIMPLIFY=FALSE)
      }
    }
    
    if(k==floor(chunks) & (nSamp%%stepSize)>0){
      start_s = c(start_s,((k)*dur)+1)
      this_soundfile <- readWave(soundfile_path, from = from, to = nSamp, units = units)
      durs = rbind(durs,round(length(this_soundfile@left)/this_soundfile@samp.rate, 2))
      timeStamps = rbind(timeStamps,as.data.frame(fileStart+((k-1)*dur)+floor(dur/j)+1))
      temp_list <- ACI_RC(this_soundfile, ...)
      return_list = mapply(rbind,return_list,temp_list,SIMPLIFY=FALSE)
      file.remove(soundfile_path)
    }
    
    if (flac == TRUE){
      if (file.exists(soundfile_path)) file.remove(soundfile_path)
    }
    
    outList = list(this_res = paste("\n", wav_files[l], ",", this_soundfile@samp.rate, ",", this_soundfile@bit, ",", durs, 
                                    ",",start_s,',', no_channels, ",", soundindex, ",", fft_w, ",", min_freq, ",", max_freq, ",",
                                    j, ",", return_list$AciTotAll_left, ",", return_list$AciTotAll_right, sep=""),
                   ACImat = return_list$aci_left_matrix,
                   timeStamps = timeStamps)
    
    rm(this_soundfile)
    cat('Done with file ',l,' of ',length(wav_files),'\n')
    
    cat(outList$this_res, file = resultfile, append = TRUE)
    if (l==1){
      allACImat = outList$ACImat
      TSDF = outList$timeStamps
    } else{
      allACImat = rbind(allACImat,outList$ACImat)
    TSDF = rbind(TSDF,outList$timeStamps)}
  }
  
  
  cat('Saving acoustic complexity matrix\n',dim(allACImat),'\n')
  saveName = paste(outDir,'/',dep[i],'_ACI_',j,'s',sep="")
  # saveName = str_sub(resultfile,1,str_locate(resultfile,'.csv')[1]-1)
  binWidth = (max_freq-min_freq)/dim(allACImat)[2]
  binFreqs = (seq(0,dim(allACImat)[2]-1)*binWidth)+min_freq
  colnames(TSDF) = c('TimeStamps')
  save(allACImat,TSDF,binFreqs,max_freq,min_freq,file=paste(saveName,'.Rdata',sep=""))
  
  #Stop timer
  time1 <- proc.time() - time0
  cat(paste(" The analysis of ", length(wav_files), " files took ", round(time1["elapsed"], 2), " seconds\n\n", sep = ""))
  
}