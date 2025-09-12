rm(list = ls())

load('data.RData')
stations <- readRDS('stations.rds')

library(dplyr)
#----Existence number of data months (hourly and 15')----
par(mfrow = c(4,4))
aux.ind <- df_hours$t + (df_hours$mes - 0.5)/12

aux.x <- tapply(aux.ind, aux.ind, mean )
for(i in estaciones){
  aux.y <- tapply(!is.na(df_hours[[paste0(i, '.p')]]), aux.ind, sum)
  plot(aux.x, aux.y, main = i,
       xlab = 'Year',
       ylab = 'No NA')
  
}

aux.ind <- df_minutes$t + (df_minutes$mes - 0.5)/12

aux.x <- tapply(aux.ind, aux.ind, mean )
for(i in estaciones){
  aux.y <- tapply(!is.na(df_minutes[[paste0(i, '.p')]]), aux.ind, sum)
  plot(aux.x, aux.y, main = i,
       xlab = 'Year',
       ylab = 'No NA')
  
}

dev.off()

#----Descriptive analysis of 0's and Non-0's----
# HOURLY AND 15' (FIRST HOURLY)
# I need info for each station

# objective. 1 dataframe for each station
df.desc.mes <- function(station, data){
  
  #info per months (due to possible stationality)
  df <- as.data.frame(matrix(NA, ncol = 1, nrow = 12))
  colnames(df) <- 'Mes'
  df$Mes <- sort(unique(data$mes))
  
  df_station <- data[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p'))]
  
  # things to compute
  # 1. rel. freq of 0's for each month (not taking into account NA's)
  cant0 <- tapply(df_station[[paste0(station, '.p')]] == 0, df_station$mes, sum, na.rm = T)
  total <- tapply(!is.na(df_station[, paste0(station, '.p')]), df_station$mes, sum)
  freq0 <- as.data.frame(cant0 / total * 100)
  
  df$f.rel.0 <- round(freq0[, 1], 3)
  
  # 2. For values greater than 0
  df_no0 <- df_station %>%
    filter(.data[[paste0(station, ".p")]] > 0)
  
  # 2.1 media
  media <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0$mes, mean, na.rm = T))
  df$media <- round(media[, 1], 3)
  
  # 2.2 mediana
  mediana <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0$mes, median, na.rm = T))
  df$mediana <- round(mediana[, 1], 3)
  
  # 2.3 cuantiles 0.90, 0.99 y 0.95
  quant.geq.0 <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0$mes, quantile, probs = c(0.9),na.rm = T))
  df$q0.90 <- round(quant.geq.0[, 1], 3)
  
  quant.geq.0 <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0$mes, quantile, probs = c(0.95),na.rm = T))
  df$q0.95 <- round(quant.geq.0[, 1], 3)
  
  quant.geq.0 <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0$mes, quantile, probs = c(0.99),na.rm = T))
  df$q0.99 <- round(quant.geq.0[, 1], 3)
  
  # 2.4 máximo
  max <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0$mes, max, na.rm = T))
  df$max <- round(max[, 1], 3)
  
  return(df)
}

for (station in estaciones){
  assign(paste0('df.h.desc.', station), df.desc.mes(station, df_hours))
}   

# in case we want the info for data every 15'
for (station in estaciones){
  assign(paste0('df.min.desc.', station), df.desc.mes(station, df_minutes))
}   

# table writing for Word
library(writexl)
write_xlsx(df.h.desc.R036, "borrar.xlsx")

# 3 Aproxximate distribution (all and per month) for non 0's
dens.no.0 <- function(station, data){
  
  df_station <- data[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p'))]
  
  df_no0 <- df_station %>%
    filter(.data[[paste0(station, ".p")]] > 0)
  
  dens <- density(df_no0[, paste0(station, ".p")], 
                  from = 0,
                  to   = max(df_no0[, paste0(station, ".p")]))
  plot(dens, lwd = 2, main = paste("Densidades de", station), xlab = "Valor")
  
  col_mes <- rainbow(12)
  for (i in 1:12){
    ind <- which(df_no0$mes == i)
    dens <- density(df_no0[ind, paste0(station, ".p")], 
                    from = 0,
                    to   = max(df_no0[ind, paste0(station, ".p")]))
    lines(dens, col = col_mes[i])
    
  }
  
  legend("topright", legend = c("Total", paste("Mes", 1:12)),
         col = c("black", col_mes), lwd = 2)
}

for (station in estaciones){
  dens.no.0(station, df_hours)
}                        
                        
# interesting graphs
# all together graphs of different columns of the previous
# frequency of 0
library(Polychrome)
cols <- createPalette(16, c("#FF0000", "#0000FF", "#00FF00", "#FFFF00",
                            '#FF00FF', '#FFA500', '#228B22', '#A020F0',
                            '#ADD8E6', '#EE82EE'))


graph.col.df <- function(stations, name.df, col, ylab, cols){
  
  dfs <- mget(paste0(name.df, estaciones), envir = .GlobalEnv)
  valores <- unlist(lapply(dfs, function(x) x[[col]]))
  min <- min(valores)
  max <- max(valores)
  
  station <- stations[1]
  aux.df <- get(paste0(name.df, station))
  plot(aux.df$Mes, aux.df[[col]], type = 'b', pch = 19,
       ylim = c(min, max),
       xlab = 'Mes', ylab = ylab,
       col = cols[1])
  for (i in 2:length(stations)){
    station <- stations[i]
    aux.df <- get(paste0(name.df, station))
    lines(aux.df$Mes, aux.df[[col]], type = 'b', pch = 19,
          col = cols[i])
  }
  
  legend("topleft", legend = stations,
         col = cols, lwd = 2, ncol = 2)
  
}

graph.col.df(estaciones, 'df.h.desc.', 'media', 'Media mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'mediana', 'Mediana mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'q0.90', 'Cuantil 0.90 mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'q0.95', 'Cuantil 0.95 mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'q0.99', 'Cuantil 0.99 mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'max', 'Máximo mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'f.rel.0', 'Freq.rel.0 mensual', cols)

graph.col.df(estaciones, 'df.min.desc.', 'media', 'Media mensual', cols)
graph.col.df(estaciones, 'df.min.desc.', 'mediana', 'Mediana mensual', cols)
graph.col.df(estaciones, 'df.min.desc.', 'q0.90', 'Cuantil 0.90 mensual', cols)
graph.col.df(estaciones, 'df.min.desc.', 'q0.95', 'Cuantil 0.95 mensual', cols)
graph.col.df(estaciones, 'df.min.desc.', 'q0.99', 'Cuantil 0.99 mensual', cols)
graph.col.df(estaciones, 'df.min.desc.', 'max', 'Máximo mensual', cols)
graph.col.df(estaciones, 'df.min.desc.', 'f.rel.0', 'Freq.rel.0 mensual', cols)

# comparision of maximum between hours and 15'
par(mfrow = c(4,4))
for (i in 1:length(estaciones)){
  station <- estaciones[i]
  df.h <- get(paste0('df.h.desc.', station))
  df.min <- get(paste0('df.min.desc.', station))
  
  plot(df.h$Mes, df.h$max, type = 'b', pch = 19,
       ylim = c(0, 60),
       xlab = 'Mes', ylab = 'Máximo mensual', main = station,
       col = cols[i])
  lines(df.min$Mes, df.min$max, type = 'b', pch = 19)
  
  legend('topleft', legend = c('h', 'min'), 
         col = c(cols[i], 'black'), lwd = 2, 
         cex = 0.8,
         pt.cex = 0.8,
         pch = 19)
}

# Possible clusterings
load('Mapas/data_mapas.RData')

df.stations.col <- function(stations, name.df, col){
  station <- stations[1]
  df.aux <- get(paste0(name.df, station))
  df.aux <- df.aux[, c('Mes', col)]
  colnames(df.aux)[2] <- paste0(colnames(df.aux)[2], '.', station)
  for (i in 2:length(stations)){
    df.aux.2 <- get(paste0(name.df, stations[i]))
    df.aux <- cbind(df.aux, df.aux.2[[col]])
    colnames(df.aux)[i+1] <- paste0(col, stations[i])
  }
  
  return(df.aux)
}

df.f.rel.0 <- df.stations.col(estaciones, 'df.h.desc.', 'f.rel.0')
df.media <- df.stations.col(estaciones, 'df.h.desc.', 'media')
df.mediana <- df.stations.col(estaciones, 'df.h.desc.', 'mediana')
df.q0.90 <- df.stations.col(estaciones, 'df.h.desc.', 'q0.90')
df.q0.95 <- df.stations.col(estaciones, 'df.h.desc.', 'q0.95')
df.q0.99 <- df.stations.col(estaciones, 'df.h.desc.', 'q0.99')
df.max <- df.stations.col(estaciones, 'df.h.desc.', 'max')

dendogram <- function(df, method.dist, method.clust){
  
  dist.aux <- dist(t(df[, 2:ncol(df)]), method = method.dist)
  hc <- hclust(dist.aux, method = method.clust)
  plot(hc)
  
  return(hc)
}

hc.f.rel.0 <- dendogram(df.f.rel.0, 'euclidean', 'single')
clust.f.rel.0 <- cutree(hc.f.rel.0, k = 2)
hc.media <- dendogram(df.media, 'euclidean', 'single')
clust.media <- cutree(hc.media, k = 2)
hc.mediana <- dendogram(df.mediana, 'euclidean', 'single')
clust.mediana <- cutree(hc.mediana, k = 2)
hc.q0.90 <- dendogram(df.q0.90, 'euclidean', 'single')
clust.q0.90 <- cutree(hc.q0.90, k = 2)
hc.q0.95 <- dendogram(df.q0.95, 'euclidean', 'single')
clust.q0.95 <- cutree(hc.q0.95, k = 2)
hc.q0.99 <- dendogram(df.q0.99, 'euclidean', 'single')
clust.q0.99 <- cutree(hc.q0.99, k = 2)
hc.max <- dendogram(df.max, 'euclidean', 'single')
clust.max <- cutree(hc.max, k = 2)

stations$grupo.f.rel.0 <- factor(c(clust.f.rel.0, 0))

# dibujar en mapa esto
library(ggplot2)
map_zone <- ggplot(hypsobath) +
  geom_sf(aes(fill = val_inf), color = NA) +
  coord_sf(xlim = st_coordinates(limits)[,1], 
           ylim = st_coordinates(limits)[,2]) + 
  scale_fill_manual(name = "Elevación", values = pal[c(7, 8:17)],
                    breaks = levels(hypsobath$val_inf),
                    guide = guide_legend(reverse = TRUE)) +
  xlab("Longitud") + ylab("Latitud") +
  geom_point(aes(x = X, y = Y, color = stations$grupo.f.rel.0), 
             data = data.frame(st_coordinates(stations))) +
  ggrepel::geom_label_repel(aes(x = X, y = Y, label = stations$STAID, 
                                color = stations$grupo.f.rel.0), 
                            size = 3.5,
                            position = 'identity', label.size = 0.025,
                            max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
                            data = data.frame(st_coordinates(stations)),
                            seed = 23) +
  scale_color_manual(values = c('0' = 'black', '1' = 'red', 
                                '2' = 'blue', '3' = 'forestgreen'),
                     name = 'Clusters') + 
  ggtitle(label = 'Grupos según frecuencia relativa 0')

map_zone



# probatinas
require(factoextra)
fviz_dist(dist.f.rel.0, show_labels = TRUE)

library(corrplot)
corrplot(cor(df.f.rel.0[, c(2:ncol(df.f.rel.0))], method = 'pearson'),
         method = 'shade',
         order = 'hclust',
         hclust.method = 'centroid')

#----NA's analysis----
# Absolute frequency of missing data per month and year 
# wih red line at level 0.25 --> 25% missing data
aux.ind <- df_minutes$t + (df_minutes$mes - 0.5)/12
aux.x <- tapply(aux.ind, aux.ind, mean )

par(mfrow = c(4,4))
for (station in estaciones){
  cont.na <- tapply(is.na(df_minutes[, paste0(station, '.p')]), aux.ind, sum)
  total <- tapply(df_minutes[, paste0(station, '.p')], aux.ind, length)
  aux.y <- cont.na / total 
  
  plot(aux.x, aux.y, pch = 19, 
       xlab = 'Year', ylab = 'Freq. abs. NA', main = station)
  abline(h = 0.25, col = 'red')
}



#----eventos extremos, para ver cuadno son las fechas --> dependencias espaciales----
basura1 <- apply(df_hours[, paste0(estaciones, '.p')] > 0, 1, mean, na.rm = T)
hist(basura)
basura2 <- apply(df_hours[, paste0(estaciones, '.p')], 1, mean, na.rm = T)
summary(basura[basura > 0])
hist(basura[basura > 0])
basura3 <- apply(df_hours[, paste0(estaciones, '.p')], 1, max, na.rm = T)
summary(basura3[basura3 > 0])
hist(basura3[basura3 > 0])