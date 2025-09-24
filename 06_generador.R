rm(list = ls())


MHO <- readRDS('MHO.rds')
MHQ <- readRDS('MHQ.rds')
MDO <- readRDS('MDO.rds')
MDQ <- readRDS('MDQ.rds')

stations <- readRDS('stations.rds')
estaciones <- stations$STAID[-nrow(stations)]
#fast check
for (station in estaciones){
  m <- MDQ[[station]]$mho
  X <- MDQ[[station]]$X
  
  cat('Estacion', station, '\n')
  cat('fv: ', length(m$mu.fv), '\n')
  cat('X: ', dim(X)[1], '\n')
  cat('Bien: ', length(m$mu.fv)==dim(X)[1], '\n')
}
