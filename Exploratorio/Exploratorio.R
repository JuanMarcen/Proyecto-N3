rm(list = ls())

load('data.RData')
stations <- readRDS('stations.rds')

library(dplyr)

#----Existence number of data months (hourly and 15')----
par(mfrow = c(7,4))
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
library(MASS)
df.desc <- function(station, data, tipo){ #añadir data de rachas para menor tiempo de computacion 
  
  if (tipo == 'mes'){
    n <- 12
  }else if (tipo == 'h'){
    n <- 24
  }
  
  #info per months (due to possible stationality)
  df <- as.data.frame(matrix(NA, ncol = 1, nrow = n))
  colnames(df) <- tipo
  df[[tipo]] <- sort(unique(data[[tipo]]))
  
  df_station <- data[, c('t', 'l', 'mes', 'dia.mes', 'h', paste0(station, '.p'))]
  
  # things to compute
  # 1. rel. freq of 0's for each month (not taking into account NA's)
  cant0 <- tapply(df_station[[paste0(station, '.p')]] == 0, df_station[[tipo]], sum, na.rm = T)
  total <- tapply(!is.na(df_station[, paste0(station, '.p')]), df_station[[tipo]], sum)
  freq0 <- as.data.frame(cant0 / total * 100)
  
  df$f.rel.0 <- round(freq0[, 1], 3)
  
  # 2. For values greater than 0
  df_no0 <- df_station %>%
    filter(.data[[paste0(station, ".p")]] > 0)
  
  # 2.1 media
  media <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], mean, na.rm = T))
  df$media <- round(media[, 1], 3)
  
  # 2.2 mediana
  mediana <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], median, na.rm = T))
  df$mediana <- round(mediana[, 1], 3)
  
  # 2.3 cuantiles 0.90, 0.99 y 0.95
  quant.geq.0 <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], quantile, probs = c(0.9),na.rm = T))
  df$q0.90 <- round(quant.geq.0[, 1], 3)
  
  quant.geq.0 <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], quantile, probs = c(0.95),na.rm = T))
  df$q0.95 <- round(quant.geq.0[, 1], 3)
  
  quant.geq.0 <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], quantile, probs = c(0.99),na.rm = T))
  df$q0.99 <- round(quant.geq.0[, 1], 3)
  
  # 2.4 máximo
  max <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], max, na.rm = T))
  df$max <- round(max[, 1], 3)
  
  # r y lambda de una gamma 
  var <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], var, na.rm = T))
  df$var <- var[, 1]
  
  fit <- tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], function(x) fitdistr(x, 'gamma'))
  df$shape.mle <- sapply(fit, function(x) x$estimate['shape'])
  df$rate.mle <- sapply(fit, function(x) x$estimate['rate'])
  
  #aux1 <- log(media[, 1])
  #aux2 <- data.frame(tapply(df_no0[, paste0(station, '.p')], df_no0[[tipo]], function(x) mean(log(x), na.rm = T)))
  #s <- aux1 - aux2[, 1]
  
  df$r <- media[, 1]^2 / var[, 1] 
  #df[['r.2']] <- (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
  df$lambda <- df$r / media[, 1]
  #df$lambda.2 <- df$r.2 / media[, 1]
  
  df$coef.var <- 1 / sqrt(df$r)
  #df$coef.var.2 <- 1 / sqrt(df$r.2)
  df$coef.var.mle <- 1 / sqrt(df$shape.mle)
  
  # RACHAS (3 h --> 3 columnas mas)
  # RACHAS (6H --> 6 COLUMNAS MAS)
  
  return(df)
}

for (station in estaciones){
  assign(paste0('df.h.desc.', station), df.desc(station, df_hours, tipo = 'mes'))
}   

for (station in estaciones){
  assign(paste0('df.h.desc.h.', station), df.desc(station, df_hours, tipo = 'h'))
}   

# for (station in estaciones){
#   assign(paste0('df.h.desc.l.', station), df.desc(station, df_hours, tipo = 'l'))
# }   

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
  
  legend("topright", legend = c("Total", paste("mes", 1:12)),
         col = c("black", col_mes), lwd = 2)
}

for (station in estaciones){
  dens.no.0(station, df_hours)
}                        
                        
# interesting graphs
# all together graphs of different columns of the previous
# frequency of 0
library(Polychrome)
cols <- createPalette(28, c("#FF0000", "#0000FF", "#00FF00", "#FFFF00",
                            '#FF00FF', '#FFA500', '#228B22', '#A020F0',
                            '#ADD8E6', '#EE82EE'))


graph.col.df <- function(stations, name.df, col, ylab, cols, tipo = 'mes'){
  
  dfs <- mget(paste0(name.df, estaciones), envir = .GlobalEnv)
  valores <- unlist(lapply(dfs, function(x) x[[col]]))
  min <- min(valores)
  max <- max(valores)
  
  station <- stations[1]
  aux.df <- get(paste0(name.df, station))
  plot(aux.df[[tipo]], aux.df[[col]], type = 'b', pch = 19,
       ylim = c(min, max),
       xlab = tipo, ylab = ylab,
       col = cols[1])
  for (i in 2:length(stations)){
    station <- stations[i]
    aux.df <- get(paste0(name.df, station))
    lines(aux.df[[tipo]], aux.df[[col]], type = 'b', pch = 19,
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
graph.col.df(estaciones, 'df.h.desc.', 'shape.mle', 'Parámetro de forma mensual', cols)
graph.col.df(estaciones, 'df.h.desc.', 'coef.var.mle', 'Coeficiente de variación mensual', cols)

graph.col.df(estaciones, 'df.h.desc.h.', 'coef.var.mle', 'Coeficiente de variación horario', cols, tipo = 'h')

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
  
  plot(df.h$mes, df.h$max, type = 'b', pch = 19,
       ylim = c(0, 60),
       xlab = 'mes', ylab = 'Máximo mensual', main = station,
       col = cols[i])
  lines(df.min$mes, df.min$max, type = 'b', pch = 19)
  
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
  df.aux <- df.aux[, c('mes', col)]
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

# dibujar en mapa esto (DE ESTO NO HAY NADA EN EL EXPLORATORIO)
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

#----Coef.var según altitud (equivalente a decir por estación)----


#----NA's analysis----
# Absolute frequency of missing data per month and year 
# wih red line at level 0.25 --> 25% missing data
aux.ind <- df_minutes$t + (df_minutes$mes - 0.5)/12
aux.x <- tapply(aux.ind, aux.ind, mean )

par(mfrow = c(7,4))
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
medias.total <- apply(df_hours[, paste0(estaciones, '.p')], 1, mean, na.rm = F) 
#na.rm = F obliga a mediciones en todas estaciones
summary(medias.total[medias.total > 0])
ev.medias <- sort(medias.total[medias.total> 0], decreasing = T)[1:6] #6 porque dos son consecutivas
aux.ind <- which(medias.total %in% ev.medias)
aux.ind <- unique(unlist(lapply(aux.ind, function(i) {
  seq(max(1, i-2), min(nrow(df_hours), i+2))
})))
aux.ind <- aux.ind[-11]

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
max.total <- apply(df_hours[, paste0(estaciones, '.p')], 1, max, na.rm = F)
summary(max.total[max.total > 0])
ev.max <- sort(max.total[max.total> 0], decreasing = T)[1:5]
aux.ind <- which(max.total %in% ev.max)
aux.ind <- unique(unlist(lapply(aux.ind, function(i) {
  seq(max(1, i-2), min(nrow(df_hours), i+2))
})))

max.total.ev <- df_hours[aux.ind, ]
aux.df[, c('t', 'mes', 'dia.mes', 'h')] <- max.total.ev[, c('t', 'mes', 'dia.mes', 'h')]
aux.df$NAs <- apply(max.total.ev, 1, function(x) sum(is.na(x)))
aux.df$sitios.lluvia <- apply(max.total.ev[, paste0(estaciones, '.p')] > 0,
                              1, sum, na.rm = T)
aux.df$frac.lluvia <- aux.df$sitios.lluvia / (length(estaciones) - aux.df$NAs)
aux.df$mean <- apply(max.total.ev[, paste0(estaciones, '.p')],
                     1, mean, na.rm = T)
aux.df$max <- apply(max.total.ev[, paste0(estaciones, '.p')],
                    1, max, na.rm = T)
aux.df$station.max <- sapply(1:length(aux.df$max), function(i) {
  names(max.total.ev)[which(max.total.ev[i, ] == aux.df$max[i])[1]]  # [1] para tomar la primera coincidencia si hay varias
})


df.eventos <- rbind(df.eventos, aux.df)

#pintar en mapa

df.mapa <-  data.frame(
  st_coordinates(stations),
  STAID = stations$STAID
)
df.mapa <- cbind(df.mapa, t(media.total.ev[, c(paste0(estaciones, '.p'))]),
                 t(max.total.ev[, c(paste0(estaciones, '.p'))]))

nombres <- c('ev.media-2h', 'ev.media-1h', 'ev.media', 'ev.media+1h', 'ev.media+2h')
nombres.final <- c()
for (i in 1:5){
  aux <- paste0(i, '.', nombres)
  nombres.final <- c(nombres.final, aux)
}
colnames(df.mapa)[4:(4 + 25 - 1)] <- nombres.final

nombres <- c('ev.max-2h', 'ev.max-1h', 'ev.max', 'ev.max+1h', 'ev.max+2h')
nombres.final <- c()
for (i in 1:5){
  aux <- paste0(i, '.', nombres)
  nombres.final <- c(nombres.final, aux)
}
colnames(df.mapa)[29:ncol(df.mapa)] <- nombres.final

mapa.ind <- function(data, event.no, event.type , lag = ''){
  col <- paste0(event.no, '.ev.', event.type, lag)
  event <- paste0(event.no, '.ev.', event.type)
  
  if (event.type == 'media'){
    fecha <- media.total.ev[(event.no * 5 - 1), c('dia.mes', 'mes', 't')]
    h <- media.total.ev[(event.no * 5 - 1), c('h')]
  }
  
  if (event.type == 'max'){
    fecha <- max.total.ev[(event.no * 5 - 1), c('dia.mes', 'mes', 't')]
    h <- max.total.ev[(event.no * 5 - 1), c('h')]
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
    # geom_point(aes(x = X, y = Y, 
    #                color = stations$color[is.na(data[[col]])]), 
    #            data = data[is.na(data[[col]]), ], 
    #            shape = 4, size = 2, stroke = 1.5) +
    
    ggrepel::geom_label_repel(aes(x = X, y = Y, 
                                  label = data[[col]], 
                                  color = stations$color), 
                              size = 3.5,
                              position = 'identity', label.size = 0.025,
                              max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
                              data = data,
                              seed = 23) +
    
    scale_color_identity() +
    ggtitle(label = paste0(paste(fecha, collapse = '/'), '-', h, 'h'))
  
  return(map_zone)
}

mapa.ind(df.mapa, 5, 'max')


mapa.ev <- function(data, event.no, event.type){
  
  m1 <- mapa.ind(data, event.no, event.type, '-2h')
  m2 <- mapa.ind(data, event.no, event.type, '-1h')
  m3 <- mapa.ind(data, event.no, event.type, '')
  m4 <- mapa.ind(data, event.no, event.type, '+1h')
  m5 <- mapa.ind(data, event.no, event.type, '+2h')
  
  mapa <- ggpubr:: ggarrange(m1, m2, m3, m4, m5, 
                             ncol = 5,
                             common.legend = T,
                             legend = 'bottom')
  return(mapa)
}

ggsave(
  filename = "Mapas/ev.extremos/1.ev.media.png", 
  plot = mapa.ev(df.mapa, 1, 'media'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/2.ev.media.png", 
  plot = mapa.ev(df.mapa, 2, 'media'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/3.ev.media.png", 
  plot = mapa.ev(df.mapa, 3, 'media'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/4.ev.media.png", 
  plot = mapa.ev(df.mapa, 4, 'media'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/5.ev.media.png", 
  plot = mapa.ev(df.mapa, 5, 'media'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/1.ev.max.png", 
  plot = mapa.ev(df.mapa, 1, 'max'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/2.ev.max.png", 
  plot = mapa.ev(df.mapa, 2, 'max'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/3.ev.max.png", 
  plot = mapa.ev(df.mapa, 3, 'max'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/4.ev.max.png", 
  plot = mapa.ev(df.mapa, 4, 'max'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)

ggsave(
  filename = "Mapas/ev.extremos/5.ev.max.png", 
  plot = mapa.ev(df.mapa, 5, 'max'), 
  width = 12*1.5,
  height = 6,     
  dpi = 300       
)


# map_zone <- ggplot(hypsobath) +
#   geom_sf(aes(fill = val_inf), color = NA) +
#   geom_sf(data = rios, color = "#40B6ED", size = 0.5) +
#   coord_sf(xlim = st_coordinates(limits)[,1], 
#            ylim = st_coordinates(limits)[,2]) + 
#   scale_fill_manual(name = "Elevación", values = pal[c(7, 8:17)],
#                     breaks = levels(hypsobath$val_inf),
#                     guide = guide_legend(reverse = TRUE)) +
#   xlab("Longitud") + ylab("Latitud") +
#   
#   geom_point(aes(x = X, y = Y, 
#                  size = media.val, 
#                  color = stations$color[!is.na(df.mapa$media.val)]), 
#              data = df.mapa[!is.na(df.mapa$media.val), ]) +
#   
#   geom_point(aes(x = X, y = Y, 
#                  color = stations$color[is.na(df.mapa$media.val)]), 
#              data = df.mapa[is.na(df.mapa$media.val), ], 
#              shape = 4, size = 2, stroke = 1.5) +
#   
#   ggrepel::geom_label_repel(aes(x = X, y = Y, 
#                                 label = stations$STAID, 
#                                 color = stations$color), 
#                             size = 3.5,
#                             position = 'identity', label.size = 0.025,
#                             max.time = 0.5, max.iter = 1000000, max.overlaps = 100,
#                             data = df.mapa,
#                             seed = 23) +
#   
#   scale_color_identity() +
#   ggtitle(label = paste('Evento extremo'))
# 
# map_zone


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

#----HERRAMIENTA ANALISIS CRITICO----
# análisis de días críticos
# filtar al periodo comun de las horas
library(lubridate)
n <- length(estaciones)
aux <- apply(df_hours[, paste0(estaciones, '.p')], 1, function(row) all(!is.na(row)))
first <- which(aux)[1]
last <- tail(which(aux), 1)
date.first <- as.Date(paste0(df_hours[first, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
date.last <- as.Date(paste0(df_hours[last, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
per.comun.h <- c(date.first, date.last)
rm(list = setdiff(ls(), c('estaciones', 'df_hours', 'df_minutes', 'stations', 'per.comun.h')))

df_hours$date <- as.Date(paste(df_hours$t, df_hours$mes, df_hours$dia.mes, sep = "-"), format = '%Y-%m-%d')
df_hours <- df_hours %>% filter(date >= per.comun.h[1] & date <= per.comun.h[2])

df_minutes$date <- as.Date(paste(df_minutes$t, df_minutes$mes, df_minutes$dia.mes, sep = "-"), format = '%Y-%m-%d')
df_minutes <- df_minutes %>% filter(date >= per.comun.h[1] & date <= per.comun.h[2])

# calculo de máximo horario regional
cols <- paste0(estaciones, ".p")
max.reg <- apply(df_hours[, paste0(estaciones, '.p')], 1, max, na.rm = T)
df_hours$max.reg <- max.reg
# guardado de nombres en que estacion
df_hours$estacion.max <- apply(df_hours[, cols], 1, function(x) {
  if (all(is.na(x))) {
    return(NA)  # si toda la fila es NA
  } else {
    idx <- which.max(x)      # posición del máximo
    return(estaciones[idx])  # nombre de la estación
  }
})
#filtrado para maximos regionales > 0
df_hours.max.gr.0 <- df_hours[df_hours$max.reg > 0, ]
thresh <- quantile(df_hours.max.gr.0[['max.reg']], prob = 0.99)
#lo que pasan el umbral
df_hours.above.thresh <- df_hours.max.gr.0[df_hours.max.gr.0[['max.reg']] > thresh, ]

#ver para esos días los datos 15'
#fijarse en date y h
df_minutes.above.thresh <- df_minutes %>%
  semi_join(df_hours.above.thresh, by = c("date", "h"))

#funcion o bucle que me devuelva para cada hora si TRUE o no segun la proporcion
ind.crit <- c()
for (i in 1:nrow(df_hours.above.thresh)){
  value.h <- df_hours.above.thresh$max.reg[i]
  col <- paste0(df_hours.above.thresh$estacion.max[i], '.p')
  date <- df_hours.above.thresh$date[i]
  h <- df_hours.above.thresh$h[i]
  
  values.min <- df_minutes.above.thresh %>%
    filter(date == !!date, h == !!h) %>%
    pull(!!sym(col))
  
  if(sum(values.min/value.h > 0.85) >= 1){
    ind.crit[i] <- TRUE
  }else{
    ind.crit[i] <- FALSE
  }
}

crit.h <- df_hours.above.thresh[ind.crit, ]
crit.min <- df_minutes.above.thresh %>%
  semi_join(crit.h, by = c("date", "h"))

#ahora para cada día crítico pintar en mapa junto a dos horas anteriores y posteriores
i <- 1
date <- crit.h$date[i]
h <- crit.h$h[i]
#match con el df_original
i.og <- which(df_hours$date == date & df_hours$h == h)


#maps of this days
load('Mapas/data_mapas.RData')
library(dplyr)
library(ggplot2)
library(st)
library(sf)
library(future.apply)
plan(multisession)
map.base <- ggplot(hypsobath) +
  geom_sf(aes(fill = val_inf), color = NA) +
  geom_sf(data = rios, color = "#40B6ED", size = 0.5) +
  coord_sf(xlim = st_coordinates(limits)[,1],
           ylim = st_coordinates(limits)[,2]) +
  scale_fill_manual(name = "Elevation", values = pal[c(7, 8:17)],
                    breaks = levels(hypsobath$val_inf),
                    guide = 'none') +
  xlab("Longitude") + ylab("Latitude")

mapa.ind <- function(stations, data.h, crit.hour, crit.hour.center,
                     map.base,
                     calculate.val.min = FALSE,
                     data.min = NULL){
  
  data.map <- data.frame(
    st_coordinates(stations),
    STAID = stations$STAID,
    p.h = t(data.h[crit.hour, paste0(estaciones, '.p')]),
    p.h.center = t(data.h[crit.hour.center, paste0(estaciones, '.p')])
  )
  
  colnames(data.map)[(ncol(data.map) - 1): ncol(data.map)] <- c('p.h', 'p.h.center')
  
  # ---- crear columna categórica con las ETIQUETAS (no nombres de color) ----
  data.map <- data.map %>%
    mutate(cat_ph = case_when(
      p.h >= 0   & p.h < 1    ~ "Low (0–1)",
      p.h >= 1   & p.h <= 7.5 ~ "Medium (1–7.5)",
      p.h > 7.5               ~ "High (>7.5)",
      TRUE                    ~ NA_character_
    ))
  
  # forzar orden de niveles de la leyenda (aquí pides High, Medium, Low)
  data.map$cat_ph <- factor(data.map$cat_ph,
                            levels = c("High (>7.5)", "Medium (1–7.5)", "Low (0–1)"))
  
  # 🔧 truco: añadir filas "fantasma" para mantener leyenda completa
  niveles_faltantes <- setdiff(levels(data.map$cat_ph), unique(data.map$cat_ph))
  if (length(niveles_faltantes) > 0) {
    data.fantasma <- data.frame(
      X = NA, Y = NA, p.h = NA, p.h.center = NA, STAID = NA,
      cat_ph = factor(niveles_faltantes, levels = levels(data.map$cat_ph))
    )
    data.map <- rbind(data.map, data.fantasma)
  }
  
  fecha <- data.h$date[crit.hour]
  h <- data.h$h[crit.hour]
  station <- data.map$STAID[which.max(data.map$p.h)]
  
  aux.df <- data.h[c(crit.hour.center -2, crit.hour.center - 1,
                     crit.hour.center,
                     crit.hour.center + 1, crit.hour.center + 2), ]
  range <- range(aux.df[, paste0(stations$STAID, '.p')], na.rm = T)
  
  if (calculate.val.min == TRUE){
    values.min <- data.min[which(data.min$date == fecha &
                                   data.min$h == h), paste0(station, '.p')]
  }else{
    values.min <- NULL
  }
  
  if (!is.null(values.min)){
    title <- paste0(fecha, '-', h, 'h', '\n', station, ". 15' data: ", 
                    paste0(values.min, collapse = '-'))
  }else{
    title <- paste0(fecha, '-', h, 'h')
  }
  
  map_zone <- map.base +
    
    # usar cat_ph EN aes() directamente
    geom_point(
      aes(x = X, y = Y,
          size = p.h,
          color = cat_ph),
      data = data.map
    ) +
    scale_size_continuous(
      name = "Rain (center)",
      range = c(2, 8),
      limits = range,
      guide = "none"
    ) +
    
    # asignar colores a las ETIQUETAS y forzar orden en la leyenda
    scale_color_manual(
      name = "Hourly rain",
      values = c(
        "High (>7.5)" = "red",
        "Medium (1–7.5)" = "blue",
        "Low (0–1)" = "forestgreen"
      ),
      breaks = c("High (>7.5)", "Medium (1–7.5)", "Low (0–1)"),
      drop = FALSE
    ) +
    
    ggrepel::geom_label_repel(
      aes(x = X, y = Y, 
          label = p.h, 
          color = cat_ph),   # usar misma variable para color de etiquetas
      size = 3.5,
      position = 'identity', label.size = 0.025,
      max.time = 0.3, max.iter = 20000, max.overlaps = 100,
      data = data.map,
      seed = 23
    ) +
    
    ggtitle(label = title)
  
  return(map_zone)
}


mapa.ev <- function(stations, data.h, crit.hour, data.min, map.base){
  
  m1 <- mapa.ind(stations, data.h = data.h, crit.hour = crit.hour - 2,
                 crit.hour.center = crit.hour, map.base = map.base)
  m2 <- mapa.ind(stations, data.h = data.h, crit.hour = crit.hour - 1,
                 crit.hour.center = crit.hour, map.base = map.base)
  m3 <- mapa.ind(stations, data.h = data.h, crit.hour = crit.hour,
                 crit.hour.center = crit.hour,
                 calculate.val.min = T, data.min = data.min, map.base = map.base)
  m4 <- mapa.ind(stations, data.h = data.h, crit.hour = crit.hour + 1,
                 crit.hour.center = crit.hour, map.base = map.base)
  m5 <- mapa.ind(stations, data.h = data.h, crit.hour = crit.hour + 2,
                 crit.hour.center = crit.hour, map.base = map.base)
  
  mapa <- ggpubr:: ggarrange(m1, m2, m3, m4, m5, 
                             ncol = 5,
                             common.legend = T,
                             legend = 'bottom')
  return(mapa)
}


crit.hours.map <- function(stations, period, data.h, data.min, quantile, prop,
                           map.base){
  
  stations.name <- stations$STAID
  
  data.h$date <- as.Date(paste(data.h$t, data.h$mes, data.h$dia.mes, sep = "-"), 
                         format = '%Y-%m-%d')
  data.h <- data.h %>% filter(date >= period[1] & date <= period[2])
  
  data.min$date <- as.Date(paste(data.min$t, data.min$mes, data.min$dia.mes, sep = "-"), 
                           format = '%Y-%m-%d')
  data.min <- data.min %>% filter(date >= period[1] & date <= period[2])
  
  # calculo de máximo horario regional
  cols <- paste0(stations.name, ".p")
  max.reg <- apply(data.h[, paste0(stations.name, '.p')], 1, max, na.rm = T)
  data.h$max.reg <- max.reg
  # guardado de nombres en que estacion
  data.h$estacion.max <- apply(data.h[, cols], 1, function(x) {
    if (all(is.na(x))) {
      return(NA)  # si toda la fila es NA
    } else {
      idx <- which.max(x)      # posición del máximo
      return(stations.name[idx])  # nombre de la estación
    }
  })
  #filtrado para maximos regionales > 0
  data.h.max.gr.0 <- data.h[data.h$max.reg > 0, ]
  thresh <- quantile(data.h.max.gr.0[['max.reg']], prob = quantile)
  #lo que pasan el umbral
  data.h.above.thresh <- data.h.max.gr.0[data.h.max.gr.0[['max.reg']] > thresh, ]
  
  #ver para esos días los datos 15'
  #fijarse en date y h
  data.min.above.thresh <- data.min %>%
    semi_join(data.h.above.thresh, by = c("date", "h"))
  
  #funcion o bucle que me devuelva para cada hora si TRUE o no segun la proporcion
  ind.crit <- c()
  for (i in 1:nrow(data.h.above.thresh)){
    value.h <- data.h.above.thresh$max.reg[i]
    col <- paste0(data.h.above.thresh$estacion.max[i], '.p')
    date <- data.h.above.thresh$date[i]
    h <- data.h.above.thresh$h[i]
    
    values.min <- data.min.above.thresh %>%
      filter(date == !!date, h == !!h) %>%
      pull(!!sym(col))
    
    if(sum(values.min/value.h > prop) >= 1){
      ind.crit[i] <- TRUE
    }else{
      ind.crit[i] <- FALSE
    }
  }
  
  crit.h <- data.h.above.thresh[ind.crit, ]
  crit.min <- data.min.above.thresh %>%
    semi_join(crit.h, by = c("date", "h"))
  
  cat('We find ', nrow(crit.h), ' possible critical hours')
  #ahora para cada día crítico pintar en mapa junto a dos horas anteriores y posteriores
  res_list <- future_lapply(1:nrow(crit.h), function(i) {
    date <- crit.h$date[i]
    h <- crit.h$h[i]
    i.og <- which(data.h$date == date & data.h$h == h)
    
    
    
    plot_i <- mapa.ev(stations, 
                      data.h = data.h, crit.hour = i.og, data.min = data.min,
                      map.base = map.base)
    
    ggsave(
      filename = paste0("Mapas/crit.h/mapa_", i, "_", date, "_", h, "h.png"),
      plot = plot_i,
      width = 16*1.75, height = 4*1.75, dpi = 150
    )
    
    return(NULL)
  })
  
  
  return(crit.h)
}


basura <- crit.hours.map(stations, per.comun.h, df_hours, df_minutes, 0.99, 0.95,
                         map.base = map.base)

plan(sequential)

#medida de similitud --> var regional + sd REGIONAL + COEF VAR (sd / mean)
# usar el periodo comun de las horas
period <- per.comun.h
# calculo de medidas
data <- df_hours
#filtar per.comun
data$date <- as.Date(paste(data$t, data$mes, data$dia.mes, sep = "-"), format = '%Y-%m-%d')
data <- data %>% filter(date >= period[1] & date <= period[2])

data$var <- apply(data[, paste0(estaciones, '.p')], 1, var, na.rm = T) 
#var = 0 --> misma lluvia en todas estaciones
data$CV <- apply(data[, paste0(estaciones, '.p')], 1, 
                 FUN = function(x) sd(x, na.rm=T) / mean(x, na.rm = T))
#CV = NaN --> mean = 0 y sd = 0

#filtrado de var > 0 y ver su distribucion
data.var.gr.0 <- data[data$var > 0 & !is.na(data$var), ]
data.var.gr.0[which.max(data.var.gr.0[, 'var']), c('date', 'h', paste0(estaciones, '.p'))]


# ver la varianza en las horas problematicas (obtenidas anteriormente)
data.var.crit.h <- data.var.gr.0 %>%
  semi_join(basura, by = c('date', 'h'))
boxplot(data.var.gr.0[, 'var'])
points(rep(1, times = nrow(basura)), data.var.crit.h[['var']], pch = 19, col = 'red')
abline(h = quantile(data.var.gr.0[['var']], 
                    probs = seq(0.9, 1, by = 0.01)), 
       col = 'blue')
plot(density(data.var.gr.0[, 'var']))
abline(v = data.var.crit.h[['var']], col = 'blue')

# var y cv
plot(data.var.gr.0[['CV']], data.var.gr.0[['var']])
points(data.var.crit.h[['CV']], data.var.crit.h[['var']], pch = 19, col = 'red')

# CV
boxplot(data.var.gr.0[, 'CV'])
points(rep(1, times = nrow(basura)), data.var.crit.h[['CV']], pch = 19, col = 'red')
