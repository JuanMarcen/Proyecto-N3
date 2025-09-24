rm(list = ls())
library(lubridate)
library(dplyr)

# Creation of stations data from csv
data_dir <- 'C:/Users/jumar/OneDrive/Escritorio/N3/datos confe'

files <- list.files(data_dir, full.names = TRUE)

estaciones <- c('EM71', 'P084', 'P018', 'E085', 'R026',
                'EM08', 'R037', 'R062',
                'EM70', 'P021', 'P087', 'P088',
                'A042', 'P023', 'P024', 'R036',
                'A266', 'A126', 'A058', 'A287', 'A288', 'EM12', 'EM07', 'R040', 'E283', 'R038', 'P019', 'P020')

stations <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 7))
colnames(stations) <- c('STAID', 'NAME', 'LAT', 'LON', 'X', 'Y', 'Z')
stations$STAID <- estaciones

for (i in 1:length(estaciones)){
  station <- estaciones[i]
  
  file <- files[grepl(station, files) & grepl("PQUIN", files)]
  
  cat(estaciones[i], ':\t', file, '\n')
  
  borrar <- read.csv(file, skip = 13, nrows = 2, header = FALSE, sep = ';')
  
  lon <- as.numeric(gsub(',', '.', gsub('.*:\t', '', borrar[1, 1])))
  
  lat <- as.numeric(gsub(',', '.', gsub('.*:\t', '', borrar[2, 1])))
  
  borrar <- read.csv(file, skip = 3, nrows = 2, header = FALSE, sep = ';')
  
  X <- as.numeric(gsub(',', '.', gsub('.*:\t', '', borrar[1, 1])))
  Y <- as.numeric(gsub(',', '.', gsub('.*:\t', '', borrar[2, 1])))
    
  borrar <- read.csv(file, skip = 1, nrows = 1, header = FALSE, sep = ';')
  borrar[, 1] <- iconv(borrar[, 1], from = "latin1", to = "UTF-8") #accents
  
  name <- gsub('.*:\t', '', borrar[, 1])
  
  stations[i, 'LAT'] <- lat
  stations[i, 'LON'] <- lon
  stations[i, 'X'] <- X
  stations[i, 'Y'] <- Y
  stations[i, 'NAME'] <- name
  
  rm(list = c('borrar', 'file', 'lat', 'lon', 'name', 'station', 'i', 'X', 'Y'))
  
}

#AÑADIR ALTURA DE LAS RESTANTES 12
stations$Z <- c(655, 855, 1394.5, 292, 761, 
               779, 1412.8, 735,
               1011, 1180, 1280, 1095,
               872.4, 1173.4, 1008.2, 1229)

stations$color <- c(rep('blue', times = 5), rep('red', times = 3),
                    rep('forestgreen', times = 4), rep('red', times = 4),
                    rep('forestgreen', times = 12))

# añadido de Ateca
stations <- rbind(
  stations,
  data.frame(
    STAID = 'A126', 
    NAME = 'ATECA', 
    LAT = 41.322944, 
    LON = -1.800639, 
    X = 600374.3, 
    Y = 4575302.1, 
    Z = 600, 
    color = 'black')
)

saveRDS(stations, 'stations.rds')

stations.text <- stations[, c('STAID', 'NAME', 'X', 'Y', 'Z')]
write.table(stations.text, "stations.txt", sep = "\t", row.names = FALSE)
write.csv(stations.text, "stations.csv", row.names = FALSE)
