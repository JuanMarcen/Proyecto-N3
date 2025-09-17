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
df.desc.mes <- function(station, data){ #añadir data de rachas para menor tiempo de computacion 
  
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
  
  # r y lambda de una gamma 
  var <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0$mes, var, na.rm = T))
  df$var <- var[, 1]
  
  df$r <- media[, 1]^2 / var[, 1] 
  df$lambda <- df$r / media[, 1]
  
  df$coef.var <- 1 / sqrt(df$r)
  
  # RACHAS (3 h --> 3 columnas mas)
  # RACHAS (6H --> 6 COLUMNAS MAS)
  
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
graph.col.df(estaciones, 'df.h.desc.', 'r', 'Parámetro de forma mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'coef.var', 'Coeficiente de variación mensual', cols)

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

hc.f.rel.0 <- dendogram(df.f.rel.0, 'euclidean', 'ward.D2')
clust.f.rel.0 <- cutree(hc.f.rel.0, k = 2)
hc.media <- dendogram(df.media, 'euclidean', 'ward.D2')
clust.media <- cutree(hc.media, k = 2)
hc.mediana <- dendogram(df.mediana, 'euclidean', 'ward.D2')
clust.mediana <- cutree(hc.mediana, k = 3)
hc.q0.90 <- dendogram(df.q0.90, 'euclidean', 'ward.D2')
clust.q0.90 <- cutree(hc.q0.90, k = 2)
hc.q0.95 <- dendogram(df.q0.95, 'euclidean', 'ward.D2')
clust.q0.95 <- cutree(hc.q0.95, k = 3)
hc.q0.99 <- dendogram(df.q0.99, 'euclidean', 'ward.D2')
clust.q0.99 <- cutree(hc.q0.99, k = 2)
hc.max <- dendogram(df.max, 'euclidean', 'ward.D2')
clust.max <- cutree(hc.max, k = 2)

stations$grupo.f.rel.0 <- factor(c(clust.f.rel.0, 0))
stations$grupo.media <- factor(c(clust.media, 0))
stations$grupo.mediana <- factor(c(clust.mediana, 0))
stations$grupo.q0.90 <- factor(c(clust.q0.90, 0))
stations$grupo.q0.95 <- factor(c(clust.q0.95, 0))
stations$grupo.q0.99 <- factor(c(clust.q0.99, 0))
stations$grupo.max <- factor(c(clust.max, 0))

# dibujar en mapa esto
library(ggplot2)
library(sf)
library(sp)
map_hc <- function(col, title){
  map_zone <- ggplot(hypsobath) +
    geom_sf(aes(fill = val_inf), color = NA) +
    geom_sf(data = rios, color = "#40B6ED", size = 0.5) +
    coord_sf(xlim = st_coordinates(limits)[,1], 
             ylim = st_coordinates(limits)[,2]) + 
    scale_fill_manual(name = "Elevación", values = pal[c(7, 8:17)],
                      breaks = levels(hypsobath$val_inf),
                      guide = guide_legend(reverse = TRUE)) +
    xlab("Longitud") + ylab("Latitud") +
    geom_point(aes(x = X, y = Y, color = stations[[paste0('grupo.', col)]]), 
               data = data.frame(st_coordinates(stations))) +
    ggrepel::geom_label_repel(aes(x = X, y = Y, label = stations$STAID, 
                                  color = stations[[paste0('grupo.', col)]]), 
                              size = 3.5,
                              position = 'identity', label.size = 0.025,
                              max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
                              data = data.frame(st_coordinates(stations)),
                              seed = 23) +
    scale_color_manual(values = c('0' = 'black', '1' = 'red', 
                                  '2' = 'blue', '3' = 'forestgreen'),
                       name = 'Clusters') + 
    ggtitle(label = paste('Grupos según', title))
  
  map_zone
}

map_hc('f.rel.0', 'frecuencia relativa 0')
map_hc('media', 'media')
map_hc('mediana', 'mediana')
map_hc('q0.90', 'q0.90')
map_hc('q0.95', 'q0.95')
map_hc('q0.99', 'q0.99')
map_hc('max', 'max')

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
df.eventos <- data.frame(matrix(NA, nrow = 1, ncol = 10))
colnames(df.eventos) <- c('t', 'mes', 'dia.mes', 'h', 
                          'NAs', 'sitios.lluvia', 'frac.lluvia', 
                          'mean', 'max', 'station.max')

# para cada estudio seleccionamos 5 eventos + 2horas antes y dos horas depsues
# --> 25 filas
aux.df <- data.frame(matrix(NA, nrow = 25, ncol = 10))
colnames(aux.df) <- c('t', 'mes', 'dia.mes', 'h', 
                          'NAs', 'sitios.lluvia', 'frac.lluvia', 
                          'mean', 'max', 'station.max')

medias.lluvia.todos <- apply(df_hours[, paste0(estaciones, '.p')] > 0, 1, mean, na.rm = T)
summary(medias.lluvia.todos)

# medias
medias.total <- apply(df_hours[, paste0(estaciones, '.p')], 1, mean, na.rm = T)
summary(medias.total[medias.total > 0])
ev.medias <- sort(medias.total[medias.total> 0], decreasing = T)[1:5]
aux.ind <- which(medias.total %in% ev.medias)
aux.ind <- unique(unlist(lapply(aux.ind, function(i) {
  seq(max(1, i-2), min(nrow(df_hours), i+2))
})))

media.total.ev <- df_hours[aux.ind, ]
aux.df[, c('t', 'mes', 'dia.mes', 'h')] <- media.total.ev[, c('t', 'mes', 'dia.mes', 'h')]
aux.df$NAs <- apply(media.total.ev, 1, function(x) sum(is.na(x)))
aux.df$sitios.lluvia <- apply(media.total.ev[, paste0(estaciones, '.p')] > 0,
                              1, sum, na.rm = T)
aux.df$frac.lluvia <- aux.df$sitios.lluvia / (length(estaciones) - aux.df$NAs)
aux.df$mean <- apply(media.total.ev[, paste0(estaciones, '.p')],
                     1, mean, na.rm = T)
aux.df$max <- apply(media.total.ev[, paste0(estaciones, '.p')],
                    1, max, na.rm = T)
aux.df$station.max <- sapply(1:length(aux.df$max), function(i) {
  names(media.total.ev)[which(media.total.ev[i, ] == aux.df$max[i])[1]]  # [1] para tomar la primera coincidencia si hay varias
})

  
df.eventos <- rbind(df.eventos, aux.df)
df.eventos <- df.eventos[-1, ]

#maximos 
max.total <- apply(df_hours[, paste0(estaciones, '.p')], 1, max, na.rm = T)
summary(max.total[max.total > 0])
quantile(max.total[max.total > 0], probs = 0.99)
max.total.max <- df_hours %>%
  filter(max.total == max(max.total, na.rm = TRUE))
aux.df[, c('t', 'mes', 'dia.mes', 'h')] <- max.total.max[, c('t', 'mes', 'dia.mes', 'h')]
aux.df$NAs <- sum(is.na(max.total.max))
aux.df$sitios.lluvia <- sum(max.total.max[, paste0(estaciones, '.p')] > 0, na.rm = T)
aux.df$frac.lluvia <- aux.df$sitios.lluvia / (length(estaciones) - aux.df$NAs)
aux.df$mean <- mean(unlist(max.total.max[, paste0(estaciones, '.p')]), na.rm = T)
aux.df$max <- max(unlist(max.total.max[, paste0(estaciones, '.p')]), na.rm = T)
aux.df$station.max <- names(max.total.max[, paste0(estaciones, '.p')])[which.max(unlist(max.total.max[, paste0(estaciones, '.p')]))] 

df.eventos <- rbind(df.eventos, aux.df)

#pintar en mapa
media.total.ev$A126.p <- rep(0, length.out = dim(media.total.ev)[1])

df.mapa <-  data.frame(
  st_coordinates(stations),
  STAID = stations$STAID
)
df.mapa <- cbind(df.mapa, t(media.total.ev[, c(paste0(estaciones, '.p'), 'A126.p')]))

nombres <- c('ev.media-2h', 'ev.media-1h', 'ev.media', 'ev.media+1h', 'ev.media+2h')
nombres.final <- c()
for (i in 1:5){
  aux <- paste0(i, '.', nombres)
  nombres.final <- c(nombres.final, aux)
}

colnames(df.mapa)[4:ncol(df.mapa)] <- nombres.final

mapa.ind <- function(data, event.no, event.type , lag = ''){
  col <- paste0(event.no, '.ev.', event.type, lag)
  event <- paste0(event.no, '.ev.', event.type)
  
  if (event.type == 'media'){
    fecha <- media.total.ev[(event.no * 5 - 1), c('dia.mes', 'mes', 't')]
  }
  
  map_zone <- ggplot(hypsobath) +
    geom_sf(aes(fill = val_inf), color = NA) +
    geom_sf(data = rios, color = "#40B6ED", size = 0.5) +
    coord_sf(xlim = st_coordinates(limits)[,1], 
             ylim = st_coordinates(limits)[,2]) + 
    scale_fill_manual(name = "Elevación", values = pal[c(7, 8:17)],
                      breaks = levels(hypsobath$val_inf),
                      guide = 'none') +
    xlab("Longitud") + ylab("Latitud") +
    
    # no NA
    geom_point(aes(x = X, y = Y, 
                   size = data[!is.na(data[[col]]),col], 
                   color = stations$color[!is.na(data[[col]])]), 
               data = data[!is.na(data[[col]]), ]) +
    scale_size_continuous(name = "Lluvia", 
                          limits = range(data[[event]], na.rm = TRUE)) +
    # NA
    geom_point(aes(x = X, y = Y, 
                   color = stations$color[is.na(data[[col]])]), 
               data = data[is.na(data[[col]]), ], 
               shape = 4, size = 2, stroke = 1.5) +
    
    ggrepel::geom_label_repel(aes(x = X, y = Y, 
                                  label = stations$STAID, 
                                  color = stations$color), 
                              size = 3.5,
                              position = 'identity', label.size = 0.025,
                              max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
                              data = data,
                              seed = 23) +
    
    scale_color_identity() +
    ggtitle(label = paste(paste(fecha, collapse = '/'), lag))
  
  return(map_zone)
}

mapa.ind(df.mapa, 5, 'media')


mapa.ev <- function(data, event){
  
  m1 <- mapa.ind(data, event, '-2h')
  m2 <- mapa.ind(data, event, '-1h')
  m3 <- mapa.ind(data, event, '')
  m4 <- mapa.ind(data, event, '+1h')
  m5 <- mapa.ind(data, event, '+2h')
  
  mapa <- ggpubr:: ggarrange(m1, m2, m3, m4, m5, 
                             ncol = 5,
                             common.legend = T,
                             legend = 'bottom')
  return(mapa)
}

mapa.ev(df.mapa, '5.ev.media')

map_zone <- ggplot(hypsobath) +
  geom_sf(aes(fill = val_inf), color = NA) +
  geom_sf(data = rios, color = "#40B6ED", size = 0.5) +
  coord_sf(xlim = st_coordinates(limits)[,1], 
           ylim = st_coordinates(limits)[,2]) + 
  scale_fill_manual(name = "Elevación", values = pal[c(7, 8:17)],
                    breaks = levels(hypsobath$val_inf),
                    guide = guide_legend(reverse = TRUE)) +
  xlab("Longitud") + ylab("Latitud") +
  
  geom_point(aes(x = X, y = Y, 
                 size = media.val, 
                 color = stations$color[!is.na(df.mapa$media.val)]), 
             data = df.mapa[!is.na(df.mapa$media.val), ]) +
  
  geom_point(aes(x = X, y = Y, 
                 color = stations$color[is.na(df.mapa$media.val)]), 
             data = df.mapa[is.na(df.mapa$media.val), ], 
             shape = 4, size = 2, stroke = 1.5) +
  
  ggrepel::geom_label_repel(aes(x = X, y = Y, 
                                label = stations$STAID, 
                                color = stations$color), 
                            size = 3.5,
                            position = 'identity', label.size = 0.025,
                            max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
                            data = df.mapa,
                            seed = 23) +
  
  scale_color_identity() +
  ggtitle(label = paste('Evento extremo'))

map_zone


#----Rachas----
# METER EN LA FUNCION DE SCAR DATOS 
station <- estaciones[1]

valor <- df_hours$EM71.p

library(zoo)
racha <- rollapply(valor, width = 3, FUN = sum, align = 'left')

aux <- df_hours

aux <- aux %>%
  mutate(
    across(all_of(paste0(estaciones, ".p")), 
           ~ rollapply(.x, width = 3, FUN = sum, align = "left", fill = NA),
           .names = "racha3h_{.col}")
  ) # align left, me dice que empieza en la hora actual y suma las dos siguientes

aux.ind <- rep(c(1:3), times = dim(df_hours)[1]/3)
aux$ind <- aux.ind

aux.marcador <- is.element(aux$ind, 3)
unique(aux[aux.marcador, 'h'])
