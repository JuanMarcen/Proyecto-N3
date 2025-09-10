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

#----Descriptive analysis of 0's and Non-0's----
# HOURLY AND 15' (FIRST HOURLY)
# I need info for each station
# example (1 station)
# objective. 1 dataframe for each station
sta <- estaciones[1]

#info per months (due to possible stationality)
del <- as.data.frame(matrix(NA, ncol = 1, nrow = 12))
colnames(del) <- 'Mes'
del$Mes <- sort(unique(df_hours$mes))

df_station <- df_hours[, c('t', 'l', 'mes', 'dia.mes', paste0(sta, '.p'))]

# things to compute
# 1. rel. freq of 0's for each month (not taking into account NA's)
cant0 <- tapply(df_station[[paste0(sta, '.p')]] == 0, df_station$mes, sum, na.rm = T)
total <- tapply(!is.na(df_station[, paste0(sta, '.p')]), df_station$mes, sum)
freq0 <- as.data.frame(cant0 / total * 100)

del$freq.rel.0 <- freq0[, 1]

# 2. For values greater than 0
# 2.1 Aproxximate distribution (all and per month)
#they should be really asymmetric
df_no0 <- df_station %>%
  filter(.data[[paste0(sta, ".p")]] > 0)


dens <- density(df_no0[, paste0(sta, ".p")], 
                from = min(df_no0[, paste0(sta, ".p")]),
                to   = max(df_no0[, paste0(sta, ".p")]))
plot(dens, lwd = 2, main = paste("Densidades de", sta), xlab = "Valor")

col_mes <- rainbow(12)
for (i in 1:12){
  ind <- which(df_no0$mes == i)
  dens <- density(df_no0[ind, paste0(sta, ".p")], 
                  from = min(df_no0[ind, paste0(sta, ".p")]),
                  to   = max(df_no0[ind, paste0(sta, ".p")]))
  lines(dens, col = col_mes[i])
  
}

legend("topright", legend = c("Total", paste("Mes", 1:12)),
       col = c("black", col_mes), lwd = 2)


# 2.2 media
mean.geq.0 <- data.frame(tapply(df_no0[, paste0(sta, '.p')], df_no0$mes, mean, na.rm = T))
del$mean.geq.0 <- mean.geq.0[, 1]

# 2.3 mediana
med.geq.0 <- data.frame(tapply(df_no0[, paste0(sta, '.p')], df_no0$mes, median, na.rm = T))
del$med.geq.0 <- med.geq.0[, 1]

# 2.4 cuantiles 0.90, 0.99 y 0.95
quant.geq.0 <- data.frame(tapply(df_no0[, paste0(sta, '.p')], df_no0$mes, quantile, probs = c(0.9),na.rm = T))
del$q0.90.geq.0 <- quant.geq.0[, 1]

quant.geq.0 <- data.frame(tapply(df_no0[, paste0(sta, '.p')], df_no0$mes, quantile, probs = c(0.95),na.rm = T))
del$q0.95.geq.0 <- quant.geq.0[, 1]

quant.geq.0 <- data.frame(tapply(df_no0[, paste0(sta, '.p')], df_no0$mes, quantile, probs = c(0.99),na.rm = T))
del$q0.99.geq.0 <- quant.geq.0[, 1]

# 2.5 máximo
max.geq.0 <- data.frame(tapply(df_no0[, paste0(sta, '.p')], df_no0$mes, max, na.rm = T))
del$max.geq.0 <- max.geq.0[, 1]                        
                        
                        

# eventos extremos, para ver cuadno son las fechas --> dependencias espaciales
basura1 <- apply(df_hours[, paste0(estaciones, '.p')] > 0, 1, mean, na.rm = T)
hist(basura)
basura2 <- apply(df_hours[, paste0(estaciones, '.p')], 1, mean, na.rm = T)
summary(basura[basura > 0])
hist(basura[basura > 0])
basura3 <- apply(df_hours[, paste0(estaciones, '.p')], 1, max, na.rm = T)
summary(basura3[basura3 > 0])
hist(basura3[basura3 > 0])