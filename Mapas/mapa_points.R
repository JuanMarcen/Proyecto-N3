rm(list = setdiff(ls(), 'stations'))

stations <- readRDS('stations.rds')

#load packages
library(dplyr)
library(sf)
library(sp)
library(rnaturalearth)
library(ggplot2)
library(viridis)
library(mapSpain)
library(tidyverse)
library(rnaturalearthdata)

# calcular background, limits
limits <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = data.frame(X = c(-2.5, -1), Y = c(40.3, 42)),
      #coords = data.frame(X = c(-2, -1), Y = c(40.3, 42)),
      #coords = data.frame(X = c(-3.2, 0.2), Y = c(39, 44)),
      data = data.frame(X = c(-2.5, -1), Y = c(40.3, 42)),
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

# background
world_map <- ne_countries(scale = "large", returnclass = 'sf')
european_union <- c("Algeria", "Andorra", "France", "Gibraltar", "Morocco", "Portugal", "Spain")
european_union_map <- 
  world_map %>% 
  filter(name %in% european_union)
background <- st_transform(european_union_map, 2062)

# stations
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")], 
      data = stations[c("STAID", "NAME", "LON", "LAT", "X", "Y", 'Z', 'color')],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

# GEOPOTENTIAL POINTS
cpoints <- data.frame(
  LON = c(-1, -1, -1, -1, 0, 0, 0, 0, -2, -2, -2, -2, -3, -3, -3, -3),
  LAT = rep(40:43, times = 4)
)

cpoints <- data.frame(
  LON = c(-1, -1, -2, -2),
  LAT = c(41, 42, 41, 42)
)

cpoints <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = cpoints[c("LON", "LAT")], 
      data = cpoints[c("LON", "LAT")],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)


#AEMET
library(readxl)
AEMET <- read_excel("C:/Users/jumar/OneDrive/Escritorio/N3/datos AEMET/estaciones-precipitación-horaria.xlsx")
AEMET <- AEMET %>%
  filter(
    NOM_PROV %in% c('SORIA', 'ZARAGOZA', 'TERUEL', 'NAVARRA') &
      AÑO_ÚLTIMO_DATO_PRECIPITACIÓN_HORARIA >= 2023 &
      COORDENADA_X >= min(stations$X) &
      COORDENADA_X <= max(stations$X) &
      COORDENADA_Y >= min(stations$Y) &
      COORDENADA_Y <= max(stations$Y) &
      FUNCIONA == 'Sí'& 
      NÚMERO_MESES_PRECIPITACIÓN_HORARIA >= 0.9 * (AÑO_ÚLTIMO_DATO_PRECIPITACIÓN_HORARIA - AÑO_PRIMER_DATO_PRECIPITACIÓN_HORARIA + 1) * 12
  )

AEMET_sf <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = AEMET[c("COORDENADA_X", "COORDENADA_Y")], 
      data = AEMET,
      proj4string = CRS("EPSG:32630") # UTM original
    ),
    'sf'
  ),
  2062  # Destino
)

# cleaning of workspace
rm(list = setdiff(ls(), c("background", "limits", "stations", "cpoints", "AEMET_sf")))

# Physical map of desired zone
hypsobath <- esp_get_hypsobath()
hypsobath <- st_transform(hypsobath, 2062)

bath_tints <- colorRampPalette( 
  rev( c( "#D8F2FE", "#C6ECFF", "#B9E3FF",
          "#ACDBFB", "#A1D2F7", "#96C9F0",
          "#8DC1EA", "#84B9E3", "#79B2DE",
          "#71ABD8" ) ) )

hyps_tints <- colorRampPalette( 
  rev( c( "#F5F4F2", "#E0DED8", "#CAC3B8", "#BAAE9A",
          "#AC9A7C", "#AA8753", "#B9985A", "#C3A76B",
          "#CAB982", "#D3CA9D", "#DED6A3", "#E8E1B6",
          "#EFEBC0", "#E1E4B5", "#D1D7AB", "#BDCC96",
          "#A8C68F", "#94BF8B", "#ACD0A5") ) )

levels <- sort(unique(hypsobath$val_inf))
br_bath <- length(levels[levels < 0])
br_terrain <- length(levels) - br_bath
pal <- c(bath_tints((br_bath)), hyps_tints((br_terrain) ))

hypsobath$val_inf[hypsobath$val_inf < 0] <- -1
hypsobath$val_inf <- as.factor(hypsobath$val_inf)
levels(hypsobath$val_inf)[1] <- "< 0"

rios <- st_read("Mapas/rios/Red_hidrografica.shp")
rios <- st_transform(rios, 2062)

save(stations,
     background,
     limits,
     hypsobath,
     pal,
     rios,
     file = 'Mapas/data_mapas.RData')



map_zone <- ggplot(hypsobath) +
  geom_sf(aes(fill = val_inf), color = NA) +
  geom_sf(data = rios, color = "#40B6ED", size = 0.5) +
  coord_sf(xlim = st_coordinates(limits)[,1], 
           ylim = st_coordinates(limits)[,2]) + 
  scale_fill_manual(name = "Elevación", values = pal[c(7, 8:17)],
                    breaks = levels(hypsobath$val_inf),
                    guide = guide_legend(reverse = TRUE)) +
  xlab("Longitud") + ylab("Latitud") +
  geom_point(aes(x = X, y = Y, color = stations$color), 
             data = data.frame(st_coordinates(stations)),
             size = 3) +
  # ESQUINAS GEOPOTENCIAL
  geom_point(aes(x = X, y = Y),
             data = data.frame(st_coordinates(cpoints)),
             col = 'black',
             size = 3,
             shape = 15) +
  # AEMET
  geom_point(aes(x = X, y = Y),
             data = data.frame(st_coordinates(AEMET_sf)),
             col = 'purple',
             size = 3,
             shape = 17) +
  #labels stations
  ggrepel::geom_label_repel(aes(x = X, y = Y, label = stations$STAID, color = stations$color),
                            size = 2.5, #original size = 3.5
                            position = 'identity', label.size = 0.025,
                            max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
                            data = data.frame(st_coordinates(stations)),
                            seed = 23) +
  #labels aemet
  ggrepel::geom_label_repel(aes(x = X, y = Y, label = AEMET_sf$INDICATIVO, color = 'purple'),
                            size = 2.5, #original size = 3.5
                            position = 'identity', label.size = 0.025,
                            max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
                            data = data.frame(st_coordinates(AEMET_sf)),
                            seed = 23) +
  scale_color_identity() + 
  ggtitle(label = 'Ateca y alrededores') 
  
map_zone


# map_zone_con_rios <- map_zone +
#   geom_sf(data = rios, color = "blue", size = 0.5)  # color azul para ríos
# 
# map_zone_con_rios
# mapa España con puntos (más adelante)
# cleaning of workspace
rm(list = setdiff(ls(), c("background", "limits", "stations")))
spain_points <- function(values, stations, title, title_legend){
  
  coords <- data.frame(
    Longitude = stations$LON, Latitude = stations$LAT)
  
  data <- cbind(coords, values)
  
  data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)
  
  # Extraer puntos que superan el umbral
  highlight_df <- subset(data, values > 0.05)
  
  # Convertimos coordenadas para que estén en el mismo sistema que el mapa (EPSG:2062)
  highlight_sf <- st_as_sf(highlight_df, coords = c("Longitude", "Latitude"), crs = 4326)
  highlight_sf <- st_transform(highlight_sf, 2062)
  highlight_coords <- st_coordinates(highlight_sf)
  highlight_df_proj <- cbind(highlight_df, X = highlight_coords[,1], Y = highlight_coords[,2])
  
  # Capa base
  p <- ggplot(data = background) +
    geom_sf(fill = "antiquewhite") +
    
    # Puntos con colores por valores
    geom_sf(data = data_sf, aes(color = values), size = 2) +
    
    # Borde rojo para valores > 0.05
    geom_point(data = highlight_df_proj,
               aes(x = X, y = Y),
               shape = 21,
               color = "red",    # color del borde
               fill = NA,        # sin relleno
               size = 3,
               stroke = 1.2) +
    
    xlab("Longitude (º)") + 
    ylab("Latitude (º)") + 
    ggtitle(title) +
    coord_sf(xlim = st_coordinates(limits)[,1], ylim = st_coordinates(limits)[,2]) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6, angle = 90),
          axis.title = element_text(size = 10, face = "bold"))  +
    
    # scale_color_gradientn(
    #   colors = c("#2c0072", "#ff6ec7", "#ffe680", "#fdae61", "#f0f0f0"),
    #   values = scales::rescale(c(0, 0.05, 0.10, 0.20, 0.30)),
    #   name = title_legend,
    #   limits = c(0, 0.30),
    #   breaks = seq(0, 0.30, 0.05)
    # ) +
    
    scale_color_viridis_c(
      option = "plasma",  # también puedes probar "plasma", "cividis", "viridis"
      name = title_legend,
      limits = c(0, 1),  # porque quieres de 0 a 1
      breaks = seq(0, 1, 0.1)
    ) +
    
    theme(
      legend.position = "bottom",
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 8)
    ) +
    
    guides(color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 15,
      barheight = 0.8
    ))
  
  return(p)
}

spain_points(values = rep(1, times = length(stations)),
             stations, 
             title = 'Ateca y alrededores',
             title_legend = NULL)


spain_points_names <- function(stations, title){
  
  # Pasar estaciones a sf
  data_sf <- st_as_sf(stations, coords = c("LON", "LAT"), crs = 4326)
  
  # Plot
  p <- ggplot(data = background) +
    geom_sf(fill = "antiquewhite") +
    
    # Estaciones en punto negro
    geom_sf(data = data_sf, color = "black", size = 2) +
    
    # Nombres de estaciones
    geom_sf_text(
      data = data_sf,
      aes(label = STAID),
      size = 2, nudge_y = 0.5
    ) +
    
    coord_sf(xlim = st_coordinates(limits)[,1], ylim = st_coordinates(limits)[,2]) +
    
    xlab("Longitude (º)") + 
    ylab("Latitude (º)") + 
    ggtitle(title) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6, angle = 90),
          axis.title = element_text(size = 10, face = "bold"))
  
  return(p)
}

spain_points_names(stations, title = 'Ateca y alrededores')








