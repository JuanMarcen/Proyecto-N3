rm(list = ls())
estaciones.nuevo <- readRDS("C:/Users/jumar/Downloads/estaciones.nuevo.rds")
names(estaciones.nuevo)

summary(estaciones.nuevo$CETINA)
aux.df <- estaciones.nuevo$CETINA
summary(aux.df[is.element(aux.df$AÑO, 2000:2025) & aux.df$Precipitacion > 10, 'Precipitacion']/10)

load('data.RData')

aux.df.2 <- df_days[, c('t', 'l', 'mes', 'dia.mes', 'A126.p')]
summary(aux.df.2[is.element(aux.df.2$t, 2000:2025) & aux.df.2$A126.p > 1, 'A126.p'])

tapply(aux.df$Precipitacion, aux.df$AÑO, function(x) sum(is.na(x)))


#---
aux.df <- estaciones.nuevo$DAROCA
summary(aux.df[is.element(aux.df$AÑO, 2000:2025) & aux.df$Precipitacion > -1, 'Precipitacion']/10)
quantile(aux.df[is.element(aux.df$AÑO, 2000:2025) & aux.df$Precipitacion > -1, 'Precipitacion']/10, 
         probs = c(0.90, 0.95, 0.99))

load('data.RData')

aux.df.2 <- df_days[, c('t', 'l', 'mes', 'dia.mes', 'A042.p')]
summary(aux.df.2[is.element(aux.df.2$t, 2000:2025) & aux.df.2$A042.p > -1, 'A042.p'])
quantile(aux.df.2[is.element(aux.df.2$t, 2000:2025) & aux.df.2$A042.p > -1, 'A042.p'],
         probs = c(0.90, 0.95, 0.99), na.rm = T)

qqplot(aux.df[is.element(aux.df$AÑO, 2000:2025) & aux.df$Precipitacion > 10, 'Precipitacion']/10,
       aux.df.2[is.element(aux.df.2$t, 2000:2025) & aux.df.2$A042.p > 1, 'A042.p'])
abline(c(0,1), col = 'red')
