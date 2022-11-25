library(data.table)
#Find all of the fimo files
FIMO_files = list.files('./')
FIMO_files = FIMO_files[grep('_fimo', FIMO_files)]
#Subset out the top 400k plus strand FIMO locations
for (z in 1:length(FIMO_files)) {
  FIMO = fread(paste(FIMO_files[z], sep = ''),
               sep = 'auto', skip = 0, header = TRUE, fill = FALSE)
  FIMO = FIMO[FIMO$strand=='+']
  FIMO = FIMO[-c(grep('_', FIMO$sequence_name))]
  if (length(grep('chrM', FIMO$sequence_name)) > 0) {
  FIMO = FIMO[-c(grep('chrM', FIMO$sequence_name))]
  }
  print(unique(FIMO$sequence_name))
  print(nrow(FIMO))
  FIMO = FIMO[order(-score)][1:400000]
  write.table(FIMO, file = paste(substr(FIMO_files[z], 1, nchar(FIMO_files[z])-9),
                                 '_plus_400k_fimo.txt', sep = ''),
              quote=FALSE, sep = '\t', row.names=FALSE, na = "")
}
#Subset out the top 400k minus strand FIMO locations
for (z in 1:length(FIMO_files)) {
  FIMO = fread(paste(FIMO_files[z], sep = ''),
               sep = 'auto', skip = 0, header = TRUE, fill = FALSE)
  FIMO = FIMO[FIMO$strand=='-']
  FIMO = FIMO[-c(grep('_', FIMO$sequence_name))]
  if (length(grep('chrM', FIMO$sequence_name)) > 0) {
    FIMO = FIMO[-c(grep('chrM', FIMO$sequence_name))]
  }
  print(unique(FIMO$sequence_name))
  print(nrow(FIMO))
  FIMO = FIMO[order(-score)][1:400000]
  write.table(FIMO, file = paste(substr(FIMO_files[z], 1, nchar(FIMO_files[z])-9),
                                 '_minus_400k_fimo.txt', sep = ''),
              quote=FALSE, sep = '\t', row.names=FALSE, na = "")
}
