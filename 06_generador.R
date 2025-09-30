rm(list = ls())


MHO <- readRDS('MHO.rds') #muy pesado (demasiado) #por que tengo donde más datos
MHQ <- readRDS('MHQ.rds')
MDO <- readRDS('MDO.rds')
MDQ <- readRDS('MDQ.rds')

library(lubridate)
library(dplyr)
library(pROC)
library(gamlss)
library(tidyr)

stations <- readRDS('stations.rds')
estaciones <- stations$STAID

global_df <- readRDS('global_df.rds')
global_df$t <- year(global_df$date)
global_df <- global_df[, -which(colnames(global_df) %in% c('zg300.', 'zg500.', 'zg700.', 
                                                           'zt300.', 'zt500.', 'zt700.'))]

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
daily.generator <- function(station, 
                      data.mo, mo, 
                      data.mq, mq, 
                      data.df,
                      n.sim){
  
  mo <- data.mo[[station]][[mo]]
  mq <- data.mq[[station]][[mq]]
  
  Xo <- data.mo[[station]]$X
  Xo$date <- as.Date(paste(Xo$t, Xo$mes, Xo$dia.mes, sep = "-"), 
                     format = '%Y-%m-%d')
  
  Xq <- data.mq[[station]]$X
  Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"), 
                     format = '%Y-%m-%d')
  
  # in order to predict correctly (USE THIS ONE TO PREDICT)
  mq2 <- gamlss(mq$mu.formula, 
                sigma.fo = mq$sigma.formula, 
                family = GA, 
                data = Xq,trace = F)
  
  
  #dataframe with climatic data
  X_data <- data.df %>%
    filter(STAID == station) %>%
    select(date) %>%
    inner_join(Xo, by = 'date') %>%
    as.data.frame()
  X_data <- X_data[, names(Xq)] # predict in gamlss, requires the same order in the columns
  
  
  #new values of mu and sigma 
  mu.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'mu', type = 'response')
  )
  
  sigma.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'sigma', type = 'response')
  )
  
  # df.gen <- Xo[, c('t', 'mes', 'dia.mes')] %>%
  #   inner_join(Xq, by = c('t', 'mes', 'dia.mes')) %>%
  #   select(t, mes, dia.mes) %>%
  #   mutate(
  #     date = as.Date(paste(t, mes, dia.mes, sep = "-"),
  #                    format = '%Y-%m-%d')
  #   ) %>%
  #   relocate(date, .before = t)
  
  
  # ind1 <- which(Xo$date %in% df.gen$date)
  # ind2 <- which(df.gen$date %in% Xq$date)
  
  df.gen <- Xo %>% select(date)
  
  df.gen <- df.gen %>%
    mutate(
      prob.p = mo$fitted.values, #change this for predict of a new df?
      mu.fv = mu.new,
      sigma.fv = sigma.new,
      shape.fv = 1 / sigma.new^2,
      rate.fv = 1 / sigma.new^2 / mu.new, #shape/mu
      p.day.obs = X_data[[paste0(station, '.p')]]
    )
  #rain = 1 or no rain = 0
  # df.gen <- df.gen %>%
  #   rowwise() %>%
  #   mutate(
  #     rain = ifelse(runif(1, min = 0, max = 1) <= prob.p, 1, 0)
  #   ) %>%
  #   ungroup()
  # 
  # df.gen <- df.gen %>%
  #   mutate(
  #     
  #   ) %>%
  #   rowwise() %>%
  #   mutate(
  #     p.day.sim = ifelse(rain == 1, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
  #   ) %>%
  #   ungroup()
  #   
  
  
  # df.gen <- df.gen %>%
  #   mutate(
  #     p.day.obs = Xq[ind2, paste0(station, '.p')]
  #   ) %>%
  #   rowwise() %>%
  #   mutate(
  #     p.day.sim = list(
  #       if (prob.p > c) rgamma(n.sim, shape = shape.fv, rate = rate.fv) else rep(0, n.sim)
  #     )
  #   ) %>%
  #   ungroup()
  # 
  # df.gen <- df.gen %>%
  #   unnest_wider(p.day.sim, names_sep = ".")
  
  # Matriz de Bernoulli (lluvia o no)
  n_days <- nrow(df.gen)
  rain_matrix <- matrix(runif(n_days * n.sim), nrow = n_days, ncol = n.sim) <= df.gen$prob.p
  
  # Matriz de gamma (simulaciones de lluvia)
  gamma_matrix <- matrix(rgamma(n_days * n.sim, shape = rep(df.gen$shape.fv, n.sim),
                                rate = rep(df.gen$rate.fv, n.sim)), 
                         nrow = n_days, ncol = n.sim)
  
  # Asignamos 0 donde no llueve
  gamma_matrix[!rain_matrix] <- 0
  
  # Convertimos la matriz a columnas en df.gen
  sim_df <- as.data.frame(gamma_matrix)
  names(sim_df) <- paste0("p.day.sim.", 1:n.sim)
  
  df.gen <- cbind(df.gen, sim_df)
  
  return(df.gen)
}
  


basura <- daily.generator(station = estaciones[4], 
                          data.mo = MDO, mo = 'M5', 
                          data.mq = MDQ, mq = 'M6', 
                          data.df = global_df, 
                          n.sim = 100)
plot(density(basura[['p.day.obs']]), col = 'blue', lwd = 2,
     xlim = c(0, max(basura[,paste0('p.day.sim.', 1:100)])))
for(i in 1:100){
  lines(
    density(basura[[paste0('p.day.sim.', i)]]), 
    col = 'red', lwd = 2
    )
}
lines(density(basura[['p.day.obs']]), col = 'blue', lwd = 2)

max <- apply(basura[, c('p.day.obs', paste0('p.day.sim.', 1:100))], 2, max)
boxplot(max)
points(max[1], pch = 19, col = 'red')

hourly.generator <- function(station, 
                            data.mo, mo, 
                            data.mq, mq, 
                            n.sim){
  
  mo <- data.mo[[station]][[mo]]
  mq <- data.mq[[station]][[mq]]
  
  Xo <- data.mo[[station]]$X
  Xo$datetime <- ISOdatetime(year = Xo$t,
                             month = Xo$mes,
                             day = Xo$dia.mes,
                             hour = Xo$h,
                             min = 0,
                             sec = 0,
                             tz = "UTC")
  
  Xq <- data.mq[[station]]$X
  Xq$datetime <- ISOdatetime(year = Xq$t,
                             month = Xq$mes,
                             day = Xq$dia.mes,
                             hour = Xq$h,
                             min = 0,
                             sec = 0,
                             tz = "UTC")
  
  df.gen <- Xo[, c('t', 'mes', 'dia.mes', 'h')] %>%
    inner_join(Xq, by = c('t', 'mes', 'dia.mes', 'h')) %>%
    select(t, mes, dia.mes, h) %>%
    mutate(
      datetime = ISOdatetime(year = t,
                             month = mes,
                             day = dia.mes,
                             hour = h,
                             min = 0,
                             sec = 0,
                             tz = "UTC")
    ) %>%
    select(datetime, t, mes, dia.mes, h)
  
  ind1 <- which(Xo$datetime %in% df.gen$datetime)
  ind2 <- which(Xq$datetime %in% df.gen$datetime)
  
  df.gen <- df.gen %>%
    mutate(
      prob.p = mo$fitted.values[ind1], #change this for predict of a new df?
      mu.fv = mq$mu.fv[ind2],
      sigma.fv = mq$sigma.fv[ind2],
      shape.fv = 1 / mq$sigma.fv[ind2]^2,
      rate.fv = 1 / mq$sigma.fv[ind2]^2 / mq$mu.fv[ind2] #shape/mu
    )
  
  roc <- suppressMessages(
    roc(Xo$Y, mo$fitted.values)
  )
  c <- as.numeric(coords(roc, "best", ret = "threshold"))
  
  df.gen <- df.gen %>%
    mutate(
      p.h.obs = Xq[ind2, paste0(station, '.p')]
    ) %>%
    rowwise() %>%
    mutate(
      p.h.sim = list(
        if (prob.p > c) rgamma(n.sim, shape = shape.fv, rate = rate.fv) else rep(0, n.sim)
      )
    ) %>%
    ungroup()
  
  df.gen <- df.gen %>%
    unnest_wider(p.h.sim, names_sep = ".")
  
  gc()
  
  return(df.gen)
}

basura2 <- hourly.generator(estaciones[25],
                            MHO, 'M5',
                            MHQ, 'M6', 
                            n.sim = 10)
plot(density(basura2[['p.h.obs']]), col = 'blue', lwd = 2)
for(i in 1:10){
  lines(
    density(basura2[[paste0('p.h.sim.', i)]]), 
    col = 'red', lwd = 2
  )
}
lines(density(basura2[['p.h.obs']]), col = 'blue', lwd = 2)

