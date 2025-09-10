rm(list = ls())

load('data.RData')

# Existence number of data months (hourly and 15')
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

basura1 <- apply(df_hours[, paste0(estaciones, '.p')] > 0, 1, mean, na.rm = T)
hist(basura)
basura2 <- apply(df_hours[, paste0(estaciones, '.p')], 1, mean, na.rm = T)
summary(basura[basura > 0])
hist(basura[basura > 0])
basura3 <- apply(df_hours[, paste0(estaciones, '.p')], 1, max, na.rm = T)
summary(basura3[basura3 > 0])
hist(basura3[basura3 > 0])