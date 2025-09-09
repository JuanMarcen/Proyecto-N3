rm(list = ls())
library(lubridate)
library(dplyr)

# Creation of stations data from csv
data_dir <- 'C:/Users/jumar/OneDrive/Escritorio/N3/datos confe'

files <- list.files(data_dir, full.names = TRUE)

estaciones <- c('EM71', 'P084', 'P018', 'E085', 'R026',
                'EM08', 'R037', 'R062',
                'EM70', 'P021', 'P087', 'P088',
                'A042', 'P023', 'P024', 'R036')

stations <- data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(stations) <- c('STAID', 'NAME', 'LAT', 'LON')
stations$STAID <- estaciones

for (i in 1:length(estaciones)){
  station <- estaciones[i]
  
  file <- files[grepl(station, files) & grepl("PQUIN", files)]
  
  borrar <- read.csv(file, skip = 13, nrows = 2, header = FALSE, sep = ';')
  
  lon <- as.numeric(gsub(',', '.', gsub('.*:\t', '', borrar[1, ])))
  
  lat <- as.numeric(gsub(',', '.', gsub('.*:\t', '', borrar[2, ])))
  
  borrar <- read.csv(file, skip = 1, nrows = 1, header = FALSE, sep = ';')
  borrar[, 1] <- iconv(borrar[, 1], from = "latin1", to = "UTF-8") #accents
  
  name <- gsub('.*:\t', '', borrar[, 1])
  
  stations[i, 'LAT'] <- lat
  stations[i, 'LON'] <- lon
  stations[i, 'NAME'] <- name
  
  rm(list = c('borrar', 'file', 'lat', 'lon', 'name', 'station', 'i'))
  
}

saveRDS(stations, 'stations.rds')
