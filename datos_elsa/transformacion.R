rm(list=ls())

#transformación de datos a como deseo 

#geopotenciales
lon <- readRDS("datos_elsa/lon.rds")
lat <- readRDS("datos_elsa/lat.rds")
g_700_12pm_60_23 <- readRDS("datos_elsa/g_700_12pm_60_23.rds")
g_500_12pm_60_23 <- readRDS("datos_elsa/g_500_12pm_60_23.rds")
g_300_12pm_60_23 <- readRDS("datos_elsa/g_300_12pm_60_23.rds")
Date <- readRDS("datos_elsa/Date.rds")

library(dplyr)
library(lubridate)
Date1<-as.Date(Date,format = "%m/%d/%Y")

#quitar 29 febrero
ind <- which(format(Date1, '%m-%d') != '02-29')

Date <- Date1[ind]
lat <- lat[ind]
lon <- lon[ind]
g300 <- g_300_12pm_60_23[ind]
g500 <- g_500_12pm_60_23[ind]
g700 <- g_700_12pm_60_23[ind]

#coordenadas unicas 11x16
coord <- data.frame(lat=lat,lon=lon)
coord <- unique(coord)
coord$station <- paste(coord$lat, coord$lon, sep = "_")

head(coord)

coord<-coord[order(-coord$lat),] #de izq a derecha
head(coord)
orden_estaciones <- unique(coord$station)
head(orden_estaciones)

#g300
g300 <- data.frame(Date = Date, lat = lat, lon = lon, g300 = g300)
g300$station <- paste0(round(g300$lat),"N.",abs(round(g300$lon)),ifelse(g300$lon>-0.5,"E","W"))
g300$lat<-NULL
g300$lon<-NULL

# Convertir el data frame de largo a ancho
g300 <- reshape(g300, 
                    idvar = "Date", 
                    timevar = "station", 
                    direction = "wide")

#colnames(g300)[-1] <- sub("^g300\\.", "", colnames(g300)[-1])
head(g300)
dim(g300)

#g500
g500 <- data.frame(Date = Date, lat = lat, lon = lon, g500 = g500)
g500$station <- paste0(round(g500$lat),"N.",abs(round(g500$lon)),ifelse(g500$lon>-0.5,"E","W"))
g500$lat<-NULL
g500$lon<-NULL

# Convertir el data frame de largo a ancho
g500 <- reshape(g500, 
                    idvar = "Date", 
                    timevar = "station", 
                    direction = "wide")


#g700
g700 <- data.frame(Date = Date, lat = lat, lon = lon, g700 = g700)
g700$station <- paste0(round(g700$lat),"N.",abs(round(g700$lon)),ifelse(g700$lon>-0.5,"E","W"))
g700$lat<-NULL
g700$lon<-NULL

# Convertir el data frame de largo a ancho
g700 <- reshape(g700, 
                    idvar = "Date", 
                    timevar = "station", 
                    direction = "wide")


saveRDS(coord, "datos_elsa/coord_g.rds")
saveRDS(g300, "datos_elsa/g300.rds")
saveRDS(g500, "datos_elsa/g500.rds")
saveRDS(g700, "datos_elsa/g700.rds")

# temperaturas
rm(list = ls())
lon <- readRDS("datos_elsa/lon_t.rds")
lat <- readRDS("datos_elsa/lat_t.rds")
t_300_12pm_60_23 <- readRDS("datos_elsa/t_300_12pm_60_23.rds")
t_500_12pm_60_23 <- readRDS("datos_elsa/t_500_12pm_60_23.rds")
t_700_12pm_60_23 <- readRDS("datos_elsa/t_700_12pm_60_23.rds")
Date <- readRDS("datos_elsa/Date_t.rds")

Date1<-as.Date(Date,format = "%m/%d/%Y")

#quitar 29 febrero
ind <- which(format(Date1, '%m-%d') != '02-29')

Date <- Date1[ind]
lat <- lat[ind]
lon <- lon[ind]
t300 <- t_300_12pm_60_23[ind]
t500 <- t_500_12pm_60_23[ind]
t700 <- t_700_12pm_60_23[ind]

#coordenadas unicas 11x16
coord <- data.frame(lat=lat,lon=lon)
coord <- unique(coord)
coord$station <- paste(coord$lat, coord$lon, sep = "_")

head(coord)

coord<-coord[order(-coord$lat),] #de izq a derecha
head(coord)
orden_estaciones <- unique(coord$station)
head(orden_estaciones)

#t300
t300 <- data.frame(Date = Date, lat = lat, lon = lon, t300 = t300)
t300$station <- paste0(round(t300$lat),"N.",abs(round(t300$lon)),ifelse(t300$lon>-0.5,"E","W"))
t300$lat<-NULL
t300$lon<-NULL

# Convertir el data frame de largo a ancho
t300 <- reshape(t300,
                idvar = "Date",
                timevar = "station",
                direction = "wide")



#t500
t500 <- data.frame(Date = Date, lat = lat, lon = lon, t500 = t500)
t500$station <- paste0(round(t500$lat),"N.",abs(round(t500$lon)),ifelse(t500$lon>-0.5,"E","W"))
t500$lat<-NULL
t500$lon<-NULL

# Convertir el data frame de largo a ancho
t500 <- reshape(t500, 
                idvar = "Date", 
                timevar = "station", 
                direction = "wide")



#t700
t700 <- data.frame(Date = Date, lat = lat, lon = lon, t700 = t700)
t700$station <- paste0(round(t700$lat),"N.",abs(round(t700$lon)),ifelse(t700$lon>-0.5,"E","W"))
t700$lat<-NULL
t700$lon<-NULL

# Convertir el data frame de largo a ancho
t700 <- reshape(t700, 
                idvar = "Date", 
                timevar = "station", 
                direction = "wide")



saveRDS(coord, "datos_elsa/coord_t.rds")
saveRDS(t300, "datos_elsa/t300.rds")
saveRDS(t500, "datos_elsa/t500.rds")
saveRDS(t700, "datos_elsa/t700.rds")
