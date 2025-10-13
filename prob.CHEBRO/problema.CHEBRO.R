# comparación lluvias AEMET CHEBRO
rm(list = ls())
estaciones.nuevo <- readRDS("C:/Users/jumar/Downloads/estaciones.nuevo.rds")
load('data.RData')

# P021, P023 -- TORNOS
summary(estaciones.nuevo$TORNOS)
aux.df <- estaciones.nuevo$TORNOS
tapply(aux.df$Precipitacion, aux.df$AÑO, function(x) sum(is.na(x)))
summary(aux.df[is.element(aux.df$AÑO, 2000:2014) 
               & aux.df$Precipitacion > 10
               , 'Precipitacion']/10)

aux.df.2 <- df_days[, c('t', 'l', 'mes', 'dia.mes', 'P021.p')]
summary(aux.df.2[is.element(aux.df.2$t, 2000:2014) 
                & aux.df.2$P021.p > 1
                 , 'P021.p'])

aux.df.3 <- df_days[, c('t', 'l', 'mes', 'dia.mes', 'P023.p')]
summary(aux.df.3[is.element(aux.df.3$t, 2000:2014) 
                 & aux.df.3$P023.p > 1
                 , 'P023.p'])

quantile(aux.df[is.element(aux.df$AÑO, 2000:2014) & aux.df$Precipitacion > -1, 'Precipitacion']/10, 
         probs = c(0.90, 0.95, 0.99), na.rm = T)
quantile(aux.df.2[is.element(aux.df.2$t, 2000:2014) & aux.df.2$P021.p > -1, 'P021.p'], 
         probs = c(0.90, 0.95, 0.99), na.rm = T)
quantile(aux.df.3[is.element(aux.df.3$t, 2000:2014) & aux.df.3$P023.p > -1, 'P023.p'], 
         probs = c(0.90, 0.95, 0.99), na.rm = T)

qqplot(aux.df[is.element(aux.df$AÑO, 2000:2014) & aux.df$Precipitacion > 10, 'Precipitacion']/10,
       aux.df.2[is.element(aux.df.2$t, 2000:2014) & aux.df.2$P021.p > 1, 'P021.p'],
       xlab = 'TORNOS', ylab = 'P021', main = 'QQ plot TORNOS - P021')
abline(c(0,1), col = 'red')

qqplot(aux.df[is.element(aux.df$AÑO, 2000:2014) & aux.df$Precipitacion > 10, 'Precipitacion']/10,
       aux.df.3[is.element(aux.df.3$t, 2000:2014) & aux.df.3$P023.p > 1, 'P023.p'],
       xlab = 'TORNOS', ylab = 'P023', main = 'QQ plot TORNOS - P023')
abline(c(0,1), col = 'red')

# A042 -- DAROCA
summary(estaciones.nuevo$DAROCA)
aux.df <- estaciones.nuevo$DAROCA
tapply(aux.df$Precipitacion, aux.df$AÑO, function(x) sum(is.na(x)))
summary(aux.df[is.element(aux.df$AÑO, 2000:2024) 
               & aux.df$Precipitacion > 10
               , 'Precipitacion']/10)

aux.df.2 <- df_days[, c('t', 'l', 'mes', 'dia.mes', 'A042.p')]
summary(aux.df.2[is.element(aux.df.2$t, 2000:2024) 
                 & aux.df.2$A042.p > 1
                 , 'A042.p'])

quantile(aux.df[is.element(aux.df$AÑO, 2000:2024) & aux.df$Precipitacion > -1, 'Precipitacion']/10, 
         probs = c(0.90, 0.95, 0.99), na.rm = T)
quantile(aux.df.2[is.element(aux.df.2$t, 2000:2024) & aux.df.2$A042.p > -1, 'A042.p'], 
         probs = c(0.90, 0.95, 0.99), na.rm = T)

qqplot(aux.df[is.element(aux.df$AÑO, 2000:2024) & aux.df$Precipitacion > 10, 'Precipitacion']/10,
       aux.df.2[is.element(aux.df.2$t, 2000:2024) & aux.df.2$A042.p > 1, 'A042.p'],
       xlab = 'DAROCA', ylab = 'A042', main = 'QQ plot DAROCA - A042')
abline(c(0,1), col = 'red')

# A126 -- CETINA

summary(estaciones.nuevo$CETINA)
aux.df <- estaciones.nuevo$CETINA
tapply(aux.df$Precipitacion, aux.df$AÑO, function(x) sum(is.na(x)))
summary(aux.df[is.element(aux.df$AÑO, 2001:2014) 
               & aux.df$Precipitacion > 10
               , 'Precipitacion']/10)

aux.df.2 <- df_days[, c('t', 'l', 'mes', 'dia.mes', 'A126.p')]
summary(aux.df.2[is.element(aux.df.2$t, 2001:2014) 
                 & aux.df.2$A126.p > 1
                 , 'A126.p'])

quantile(aux.df[is.element(aux.df$AÑO, 2001:2014) & aux.df$Precipitacion > -1, 'Precipitacion']/10, 
         probs = c(0.90, 0.95, 0.99), na.rm = T)
quantile(aux.df.2[is.element(aux.df.2$t, 2001:2014) & aux.df.2$A126.p > -1, 'A126.p'], 
         probs = c(0.90, 0.95, 0.99), na.rm = T)

qqplot(aux.df[is.element(aux.df$AÑO, 2001:2014) & aux.df$Precipitacion > 10, 'Precipitacion']/10,
       aux.df.2[is.element(aux.df.2$t, 2001:2014) & aux.df.2$A126.p > 1, 'A126.p'],
       xlab = 'CETINA', ylab = 'A126', main = 'QQ plot CETINA - A126')
abline(c(0,1), col = 'red')
