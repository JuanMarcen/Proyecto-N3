rm(list = ls())


MHO <- readRDS('MHO.rds') #muy pesado (demasiado)
MHQ <- readRDS('MHQ.rds')
MDO <- readRDS('MDO.rds')
MDQ <- readRDS('MDQ.rds')

library(lubridate)
library(dplyr)
library(pROC)
library(gamlss)

stations <- readRDS('stations.rds')
estaciones <- stations$STAID
#fast check
for (station in estaciones){
  m <- MDO[[station]]$M5
  X <- MDO[[station]]$X
  
  cat('Estacion', station, '\n')
  cat('fv: ', length(m$fitted.values), '\n')
  cat('X: ', dim(X)[1], '\n')
  cat('Bien: ', length(m$fitted.values)==dim(X)[1], '\n')
}

# all in check
#----DAILY GENERATOR IN P019----
# this station looks good in the controls carried out
# is the one with the most data available (since 1997)
station <- 'P019'
mdo <- MDO[[station]]$M5
mdq <- MDQ[[station]]$M6

Xo <- MDO[[station]]$X
Xo$date <- as.Date(paste(Xo$t, Xo$mes, Xo$dia.mes, sep = "-"), 
                   format = '%Y-%m-%d')

Xq <- MDQ[[station]]$X
Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"), 
                   format = '%Y-%m-%d')

df.gen <- Xo[, c('t', 'mes', 'dia.mes')] %>%
  inner_join(Xq, by = c('t', 'mes', 'dia.mes')) %>%
  select(t, mes, dia.mes) %>%
  mutate(
    date = as.Date(paste(t, mes, dia.mes, sep = "-"),
                   format = '%Y-%m-%d')
  ) %>%
  relocate(date, .before = t)

ind1 <- which(Xo$date %in% df.gen$date)
ind2 <- which(Xq$date %in% df.gen$date)

df.gen <- df.gen %>%
   mutate(
    prob.p = mdo$fitted.values[ind1], #change this for predict of a new df?
    mu.fv = mdq$mu.fv[ind2],
    sigma.fv = mdq$sigma.fv[ind2],
    shape.fv = 1 / mdq$sigma.fv[ind2]^2,
    rate.fv = 1 / mdq$sigma.fv[ind2]^2 / mdq$mu.fv[ind2] #shape/mu
  )

roc <- roc(Xo$Y, mdo$fitted.values)
c <- as.numeric(coords(roc, "best", ret = "threshold"))

df.gen <- df.gen %>%
  rowwise() %>%
  mutate(
    p.day = ifelse(prob.p > c, rgamma(1, shape = shape.fv, rate = rate.fv) , 0)
  ) %>%
  ungroup()

plot(density(df.gen$p.day), col = 'red', lwd = 2)
lines(density(Xq$P019.p[ind2]), col = 'blue', lwd = 2)  
