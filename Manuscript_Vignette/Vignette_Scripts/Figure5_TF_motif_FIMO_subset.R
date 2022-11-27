library(data.table)
FIMO_files = list.files('./')
FIMO_files = FIMO_files[grep('rm_FIMO.txt', FIMO_files)]
#Make 400k plus/minus strand reads
for (z in 1:length(FIMO_files)) {
  #Load in FIMO file    
  FIMO = fread(paste('./', FIMO_files[z], sep = ''),
               sep = 'auto', skip = 0, header = TRUE, fill = FALSE)
  #Remove alternative and mitochondrial chromosomes
  FIMO = FIMO[-c(grep('_', FIMO$sequence_name))]
  FIMO = FIMO[-c(grep('chrM', FIMO$sequence_name))]
  print(unique(FIMO$sequence_name))
  
  FIMO_minus = FIMO[FIMO$strand=='-']
  FIMO_minus = FIMO_minus[order(-score)][1:400000]
  write.table(FIMO_minus, file = paste(substr(FIMO_files[z], 1, nchar(FIMO_files[z])-9),
                                       '_minus_400k_fimo.txt', sep = ''), 
              quote=FALSE, sep = '\t', row.names=FALSE, na = "")
  
  FIMO_plus = FIMO[FIMO$strand=='+']
  FIMO_plus = FIMO_plus[order(-score)][1:400000]
  write.table(FIMO_plus, file = paste(substr(FIMO_files[z], 1, nchar(FIMO_files[z])-9),
                                      '_plus_400k_fimo.txt', sep = ''), 
              quote=FALSE, sep = '\t', row.names=FALSE, na = "")
}