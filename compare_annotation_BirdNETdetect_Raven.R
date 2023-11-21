# install this package
#install.packages("data.table")   
library("data.table")
library(stringr)

# This function uses raven selection tables
# BirdNET detection table need to have the
# annotation filename included as part of it's name, for example
# annotation file : "Rausu03_44100Hz_20220320_220617.txt"
# BirdNET detection file: "Rausu03_44100Hz_20220320_220617.BirdNET.selection.table_model3.txt"
#

# Run from here
comp_AnDetRav <- function(direct_anot, direct_detect) {

annotation_files <- dir(direct_anot,pattern='.txt')
  
detection_files <- list.files(direct_detect, pattern = "txt")

res_fin <- data.frame()
sum_fin <- data.frame()

for (txt_file in 1:length(annotation_files)) {

file_annotation <- paste(direct_anot, annotation_files[txt_file], sep = "/")
dep = str_sub(annotation_files[txt_file],1,4)
annotations <-   read.table(file_annotation, header = TRUE, sep = "\t")

# Remove annotation type 'waveform'
annotations <- subset(annotations, annotations$View == "Spectrogram 1")

# Find corresponding detector output selection table
matchingFile = str_which(detection_files,dep)
file_detec <-  paste(direct_detect,detection_files[matchingFile],sep = "/")
detections <- read.table(file_detec, header = TRUE, sep = "\t")

# Remove detection type 'waveform'
detections <- subset(detections, detections$View == "Spectrogram 1")

col_a_begin <- colnames(annotations)[grep("Begin.T", colnames(annotations))]
col_a_end <- colnames(annotations)[grep("End.T", colnames(annotations))]

col_d_begin <- colnames(detections)[grep("Begin.T", colnames(detections))]
col_d_end <- colnames(detections)[grep("End", colnames(detections))]

key_ann <- data.table(annotations)
setkeyv(key_ann, c(col_a_begin, col_a_end))


res_loop <- foverlaps(data.table(detections), 
          key_ann, type="any") ## return overlap indices

res_loop$anot_file <- annotation_files[txt_file]

res_fin <- rbind(res_fin, res_loop)


sum_loop <- data.frame(annotation_file = annotation_files[txt_file],
  n_tot_calls = length(unique(key_ann$Selection)),
  n_calls_detected = length(unique(res_loop$Selection)),
  n_calls_missed = length(which(!unique(annotations$Selection) %in% unique(res_loop$Selection))),
  n_detections =  nrow(detections),
  n_TP = length(which(is.na(res_loop$Selection) == FALSE)),
  n_FP = length(which(is.na(res_loop$Selection) == TRUE)))
  
sum_fin <- rbind(sum_fin, sum_loop)

}

list_res <- list()
list_res[[1]] <- res_fin
list_res[[2]] <- sum_fin

return(list_res)
}

# To here

# Now let's use the function

direct_anot <- "test/direct_anot/" # The directory path containing the annotation files
direct_detect <- "test/direct_detections//" # The directory path containing BirdNET detections

# Run the function and save the results in an object
res <- comp_AnDetRav(direct_anot = direct_anot, direct_detect = direct_detect)

# Final total table
res[[1]]

# Final summaries
res[[2]]
  