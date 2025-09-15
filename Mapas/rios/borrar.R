#to plot the map, centered in aragon
aux.xlim <- c(600000, 640000)
aux.ylim <- c(4480000, 4800000)
#to plot the map, showing the Mediterranea sea
aux.xlim <- c(600000, 840000)
aux.ylim <- c(4480000, 4800000)

library(sp)
install.packages("Mapas/rios/rgdal_1.6-7.tar.gz", repos = NULL, type = "source")
library(rgdal)
library(sp)
library(sf)
plot(provPoly[ c(1:111)[!is.element(1:111, c(28, 29, 32, 42, 58, 73, 77, 80, 81, 82, 83))] , ], 
     axes=F, cex.axis=0.5, xlim=aux.xlim , ylim=aux.ylim)
#ordenes generales, con los rios
image(elevGrid, col= escala.colores  , breaks=brks, add=T)
#Para superponer la cuenca del Ebro, traida de GEOPORTAL ebro
basura <- readOGR(dsn='Mapas/rios', layer='Limite_Cuenca_Ebro')
plot(basura, add=T)
basura1 <- readOGR(dsn='.', layer='Red_hidrografica')
plot(basura1, add=T, col=4)

# versiones más nuevas
library(sf)
library(terra)

# Provincias
provPoly <- st_read("Mapas/rios/provpoly.shp")

# Límite cuenca
basura <- st_read("Mapas/rios/Limite_Cuenca_Ebro.shp")

# Red hidrográfica
basura1 <- st_read("Mapas/rios/Red_hidrografica.shp")

# Raster
elevGrid <- rast("Mapas/rios/mde5000.asc")

# Pintar
plot(elevGrid, xlim=aux.xlim, ylim=aux.ylim, col=terrain.colors(20))
plot(st_geometry(provPoly), add=TRUE)
plot(st_geometry(basura), add=TRUE)
plot(st_geometry(basura1), add=TRUE, col="blue")



¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡
#Pintar la zona

map('world', xlim=c(-5,4), ylim=c(39.5,44))
abline(v=c(-2,-1,0,0.5))
abline(h=c(40:44))
#Superponer la cuenca

perfil.CHE <- readOGR(dsn='.', layer='Limite_Cuenca_Ebro')
perfil.CHE.longitud.latitud <- spTransform(perfil.CHE, CRS("+proj=longlat"))
plot(perfil.CHE.longitud.latitud , add=T)

points(existencias.Tx$long, existencias.Tx$lat, cex=0.4)

text(existencias.Tx$long, existencias.Tx$lat, labels= existencias.Tx$Indicativo, cex=0.4)

#los observatorios usados para el ajuste del modelo
points(df.Aragon.summary.Tmax[, c(3,2)], pch=16, col=3)