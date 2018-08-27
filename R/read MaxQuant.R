# read MaxQuant

read.MaxQuant <- function(file){
  read.table(file, header = T, sep = "\t", fill = T)
}
