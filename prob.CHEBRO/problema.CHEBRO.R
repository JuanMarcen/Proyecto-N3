# comparación lluvias AEMET CHEBRO
rm(list = ls())
estaciones.nuevo <- readRDS("C:/Users/jumar/Downloads/estaciones.nuevo.rds")
load('data.RData')

library(dplyr)
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


#days with most rain
aux.ind <- order(aux.df$Precipitacion[is.element(aux.df$AÑO, 2000:2014)], decreasing = T)[1:10]
aux.ind <- (unlist(lapply(aux.ind, function(i) {
  seq(max(1, i-1), min(nrow(df_hours), i+1))
})))
prec.df <- aux.df[is.element(aux.df$AÑO, 2000:2014), ][aux.ind, ]

aux.df.2 <- df_days_2[, c('t_meteo', 'l_meteo', 'P021.p')]
aux.df.3 <- df_days_2[, c('t_meteo', 'l_meteo', 'P023.p')]
prec.df <- left_join(
  prec.df, aux.df.2[, -4],
  by = c('DIA' = 'l_meteo', 'AÑO' = 't_meteo')
)
prec.df <- left_join(
  prec.df, aux.df.3[, -4],
  by = c('DIA' = 'l_meteo', 'AÑO' = 't_meteo')
)
prec.df$Precipitacion <- prec.df$Precipitacion/10

library(ggplot2)
library(tidyr)  # para pivot_longer

# Si no tienes el índice aún
prec.df$index <- 1:nrow(prec.df)

# Reformamos los datos al formato largo
prec_long <- prec.df %>%
  pivot_longer(
    cols = c(Precipitacion, P021.p, P023.p),
    names_to = "Variable",
    values_to = "Valor"
  )

# Gráfico agrupado
ggplot(prec_long, aes(x = factor(index), y = Valor, fill = Variable)) +
  geom_col(position = "dodge") +  # barras lado a lado
  geom_vline(xintercept = seq(3.5, length(unique(prec_long$index)), by = 3),
             linetype = "dashed", color = "gray40") +
  labs(x = "Event",
       y = "Daily precipitation",
       title = "Daily precipitation comparison of extreme event",
       fill = "Variable") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(breaks = seq(2, nrow(prec.df), by = 3), labels = prec.df[seq(2, nrow(prec.df), by = 3), 'FECHA'])


# summaries en días más tochos
summary.min <- function(prec.df, station.chebro, event){
  par(mfrow = c(1,2))
  aux.ind <- seq(2, 30, by = 3)
  aux.df.4 <- df_minutes %>%
    filter(
      t_meteo == prec.df$AÑO[aux.ind[event]] & 
        l_meteo == prec.df$DIA[aux.ind[event]]
    )
  cat('\nSummary datos minutales ', station.chebro, ' fecha: ', as.character(prec.df$FECHA[event]), '\n')
  print(summary(aux.df.4[[station.chebro]]))
  boxplot(aux.df.4[[station.chebro]], 
          main = paste(station.chebro, '-', prec.df$FECHA[event]),
          ylab = station.chebro)
  plot(aux.df.4[[station.chebro]], xaxt = 'n', type = 'b', pch = 19,
       xlab = 'h', ylab = station.chebro, 
       main = paste(station.chebro, '-', prec.df$FECHA[event]),)
  axis(side = 1, at = seq(1, nrow(aux.df.4), by = 4), labels = c(8:23, 0:7))
}

for (i in 1:10){
  summary.min(prec.df, 'P023.p', i)
}

for (i in 1:10){
  summary.min(prec.df, 'P021.p', i)
}






#----A042 -- DAROCA----
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

#days with most rain
aux.ind <- order(aux.df$Precipitacion[is.element(aux.df$AÑO, 2000:2024)], decreasing = T)[1:10]
aux.ind <- (unlist(lapply(aux.ind, function(i) {
  seq(max(1, i-1), min(nrow(df_hours), i+1))
})))
prec.df <- aux.df[is.element(aux.df$AÑO, 2000:2024), ][aux.ind, ]

aux.df.2 <- df_days_2[, c('t_meteo', 'l_meteo', 'A042.p')]


prec.df <- left_join(
  prec.df, aux.df.2[, -4],
  by = c('DIA' = 'l_meteo', 'AÑO' = 't_meteo')
)

prec.df$Precipitacion <- prec.df$Precipitacion/10

library(ggplot2)
library(tidyr)  # para pivot_longer

# Si no tienes el índice aún
prec.df$index <- 1:nrow(prec.df)

# Reformamos los datos al formato largo
prec_long <- prec.df %>%
  pivot_longer(
    cols = c(Precipitacion, A042.p),
    names_to = "Variable",
    values_to = "Valor"
  )

# Gráfico agrupado
ggplot(prec_long, aes(x = factor(index), y = Valor, fill = Variable)) +
  geom_col(position = "dodge") +  # barras lado a lado
  geom_vline(xintercept = seq(3.5, length(unique(prec_long$index)), by = 3),
             linetype = "dashed", color = "gray40") +
  labs(x = "Event",
       y = "Daily precipitation",
       title = "Daily precipitation comparison of extreme event",
       fill = "Variable") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(breaks = seq(2, nrow(prec.df), by = 3), labels = prec.df[seq(2, nrow(prec.df), by = 3), 'FECHA'])


for (i in 1:10){
  summary.min(prec.df, 'A042.p', i)
}

#----A126 -- CETINA----

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


aux.ind <- order(aux.df$Precipitacion[is.element(aux.df$AÑO, 2001:2014)], decreasing = T)[1:10]
aux.ind <- (unlist(lapply(aux.ind, function(i) {
  seq(max(1, i-1), min(nrow(df_hours), i+1))
})))
prec.df <- aux.df[is.element(aux.df$AÑO, 2001:2014), ][aux.ind, ]

aux.df.2 <- df_days_2[, c('t_meteo', 'l_meteo', 'A126.p')]

prec.df <- left_join(
  prec.df, aux.df.2[, -4],
  by = c('DIA' = 'l_meteo', 'AÑO' = 't_meteo')
)

prec.df$Precipitacion <- prec.df$Precipitacion/10

library(ggplot2)
library(tidyr)  # para pivot_longer

# Si no tienes el índice aún
prec.df$index <- 1:nrow(prec.df)

# Reformamos los datos al formato largo
prec_long <- prec.df %>%
  pivot_longer(
    cols = c(Precipitacion, A126.p),
    names_to = "Variable",
    values_to = "Valor"
  )

# Gráfico agrupado
ggplot(prec_long, aes(x = factor(index), y = Valor, fill = Variable)) +
  geom_col(position = "dodge") +  # barras lado a lado
  geom_vline(xintercept = seq(3.5, length(unique(prec_long$index)), by = 3),
             linetype = "dashed", color = "gray40") +
  labs(x = "Event",
       y = "Daily precipitation",
       title = "Daily precipitation comparison of extreme event",
       fill = "Variable") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(breaks = seq(2, nrow(prec.df), by = 3), 
                   labels = prec.df[seq(2, nrow(prec.df), by = 3), 'FECHA']) +
  scale_fill_manual(values = c("skyblue", "orange"), 
                    labels = c('A126','CETINA'))
  

for (i in 1:10){
  summary.min(prec.df, 'A126.p', i)
}




####----
aux.df <- estaciones.nuevo$CETINA
summary(df_minutes$EM71.p)
summary(df_minutes$P020.p)
summary(df_minutes$P018.p)
apply(df_minutes[, 7:ncol(df_minutes)], 2, max, na.rm = T)
