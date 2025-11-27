# COMPARISON OF MODELS USING DIFFERENT METRICS: SCC, RMSE, RB
rm(list = setdiff(ls(), c('MHQ', 'MDQ', 'common.models.final',
                          'per.comun.day', 'per.comun.h', 'estaciones')))
# load data and functions
load('data.RData')
rm(list = setdiff(ls(), c('estaciones')))
MDQ <- qread('MDQ.qs') # models and data for each model
MHQ <- qread('MHQ.qs')
source('methods.R')

X.MHQ <- qread('X.MHQ.qs')
X.MDQ <- qread('X.MDQ.qs')

library(lubridate)
library(dplyr)
# period of reference (common period of all stations)
ref.period <- as.Date(c('12/05/2011', '30/11/2023'), format = '%d/%m/%Y')



# Anual mean precipitation in JJA
value.obs.sim.station <- function(station, data, ref.period, models.list, model, 
                                 type, quantile = NULL,
                                 n.sim = 100, months = NULL,
                                 adjusted.ref.period = FALSE){
  # observed values
  X <- data[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = "%Y-%m-%d")
  x.obs <- X[, c('date', paste0(station, '.p'))]
  
  if (!is.null(months)){
    ind.ref <- which(x.obs$date >= ref.period[1] & x.obs$date <= ref.period[2]
                     & month(x.obs$date) %in% months)
  }else{
    ind.ref <- which(x.obs$date >= ref.period[1] & x.obs$date <= ref.period[2])
  }
  
  
  x.obs <- x.obs[ind.ref, ]
  
  
  
  #quantiles, mean, max, ... anything that has an R built in function
  x.obs.final <- tryCatch({
    
    fun <- match.fun(type)  # error
    
    if (type == "quantile") {
      x.obs %>% 
        group_by(year(date)) %>%
        summarise(
          mean.p = fun(.data[[paste0(station, '.p')]], na.rm = TRUE, probs = quantile)
        )
    } else {
      x.obs %>% 
        group_by(year(date)) %>%
        summarise(
          mean.p = fun(.data[[paste0(station, '.p')]], na.rm = TRUE)
        )
    }
    
  }, error = function(e) {
    # Si hay error, ejecuta esta alternativa:
    if (type == 'freq'){
      x.obs %>%
        group_by(year(date)) %>%
        summarise(
          freq = sum(.data[[paste0(station, '.p')]] > 0.1) / n()
        ) %>%
        ungroup()
    } else if (type == 'intensity'){
      x.obs %>%
        group_by(year(date)) %>%
        summarise(
          intensity = mean(.data[[paste0(station, '.p')]][.data[[paste0(station, '.p')]] > 0.1], na.rm = TRUE)
        ) %>%
        ungroup()
    } else {
      stop("'type' no válido. Usa 'una función de R-base entre '', 'freq' o 'intensity'.")
    }
    
  })
  

  
  
  x.obs.final <- mean(x.obs.final[[2]])
  #simulation values
  mod <- models.list[[station]][[model]] #desired model
  
  if (adjusted.ref.period == TRUE){
    ind.ref <- 1:length(mod$sigma.fv)
  }else{
    ind.ref <- ind.ref
  }
  
  vector.sim <- c()
  for (i in 1:n.sim){
    shape <- 1 / mod$sigma.fv[ind.ref]^2
    rate <- shape / mod$mu.fv[ind.ref]
    y.sim <- rgamma(length(ind.ref), shape = shape, rate = rate)
    
    y.sim.df <- data.frame(
      date = x.obs$date,
      y.sim = y.sim
    )
    
    y.sim.final <- tryCatch({
      
      fun <- match.fun(type)  # error
      
      if (type == "quantile") {
        y.sim.df %>% 
          group_by(year(date)) %>%
          summarise(
            mean.p = fun(y.sim, na.rm = TRUE, probs = quantile)
          )
      } else {
        y.sim.df %>% 
          group_by(year(date)) %>%
          summarise(
            mean.p = fun(y.sim, na.rm = TRUE)
          )
      }
      
    }, error = function(e) {
      # Si hay error, ejecuta esta alternativa:
      if (type == 'freq'){
        y.sim.df %>%
          group_by(year(date)) %>%
          summarise(
            freq = sum(y.sim > 0.1) / n()
          ) %>%
          ungroup()
      } else if (type == 'intensity'){
        y.sim.df %>%
          group_by(year(date)) %>%
          summarise(
            intensity = mean(y.sim[y.sim > 0.1], na.rm = TRUE)
          ) %>%
          ungroup()
      } else {
        stop("'type' no válido. Usa 'una función de R-base entre '', 'freq' o 'intensity'.")
      }
      
    })
    
    y.sim.final <- mean(y.sim.final[[2]])
    vector.sim[i] <- y.sim.final
  }
  
  
  return(c(x.obs.final, vector.sim))
}


value.obs.sim.full <- function(estaciones, data, ref.period, models.list, model, 
                               type, quantile = NULL,
                              n.sim = 100, months = NULL,
                              adjusted.ref.period = FALSE){
  N <- length(estaciones)
  df <- data.frame(matrix(NA, ncol = n.sim + 2, nrow = length(estaciones)))
  colnames(df) <- c('station', 'x.obs', paste0('y.sim.', 1:n.sim))
  df$station <- estaciones
  
  names.fill <- c('x.obs', paste0('y.sim.', 1:n.sim))
  
  for (i in 1:length(estaciones)){
    #cat('Cálculo para la estación: ', estaciones[i], '\n')
    value.station <- value.obs.sim.station(estaciones[i],
                                           data = data,
                                           ref.period = ref.period,
                                           models.list = models.list,
                                           model = model,
                                           n.sim = n.sim,
                                           months = months,
                                           type = type, 
                                           quantile = quantile,
                                           adjusted.ref.period = adjusted.ref.period)
    df[i, names.fill] <- value.station
  }
  
  return(df)
}

mean.metric <- function(metric, data){
  
  fun <- match.fun(metric)

  values <- c()
  
  x.obs <- data[[2]] # observed values will be always in the second column
  
  for (i in 3:ncol(data)){ #simulation values rest of columns 3 to the
    y.sim <- data[[i]]
    values[i-2] <- fun(x.obs, y.sim)
  }
  
  final.value <- mean(values)
  
  return(final.value)
}

#MHQ
df.amounts.M5 <- value.obs.sim.full(estaciones[-4], ref.period, MHQ, 'M4', n.sim = 100, type = 'mean')
df.amounts.M6 <- value.obs.sim.full(estaciones[-4], ref.period, MHQ, 'M6', n.sim = 100, type = 'mean')

mean.metric('SCC', df.amounts.M6)  
mean.metric('SCC', df.amounts.M5)
mean.metric('RMSE', df.amounts.M6)
mean.metric('RMSE', df.amounts.M5)
mean.metric('RB', df.amounts.M6)
mean.metric('RB', df.amounts.M5)

#MDQ
value.obs.sim.full(estaciones, X.MDQ, per.comun.day, MDQ, model = 'M8', type = 'mean')
day.amounts.M5 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M5', n.sim = 100, type = 'mean')
day.amounts.M6 <- value.obs.sim.full(estaciones, X.MDQ, per.comun.day, MDQ, 'M6', n.sim = 100, type = 'mean')

day.q0.90.M5 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M9', n.sim = 100, type = 'quantile', quantile = 0.90)
day.q0.90.M6 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M6', n.sim = 100, type = 'quantile', quantile = 0.90)

day.q0.95.M5 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M5', n.sim = 100, type = 'quantile', quantile = 0.95)
day.q0.95.M6 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M6', n.sim = 100, type = 'quantile', quantile = 0.95)

day.q0.99.M5 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M5', n.sim = 100, type = 'quantile', quantile = 0.99)
day.q0.99.M6 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M6', n.sim = 100, type = 'quantile', quantile = 0.99)

day.freq.M5 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M5', n.sim = 100, type = 'freq')
day.freq.M6 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M6', n.sim = 100, type = 'freq')

day.int.M5 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M5', n.sim = 100, type = 'intensity')
day.int.M6 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M6', n.sim = 100, type = 'intensity')

mean.metric('SCC', day.amounts.M6)  
mean.metric('SCC', day.amounts.M5)
mean.metric('RMSE', day.amounts.M6)
mean.metric('RMSE', day.amounts.M5)
mean.metric('RB', day.amounts.M6)
mean.metric('RB', day.amounts.M5)

model.comp.df <- function(estaciones, ref.period, data, models.list, m1, m2, n.sim = 100){
  amounts.m1 <- value.obs.sim.full(estaciones, data, ref.period,
                                   models.list, m1, type = 'mean', 
                                   adjusted.ref.period = F)
  amounts.m2 <- value.obs.sim.full(estaciones, data, ref.period,
                                   models.list, m2, type = 'mean', 
                                   adjusted.ref.period = F)
  
  q0.90.m1 <- value.obs.sim.full(estaciones, data, ref.period,
                                 models.list, m1, type = 'quantile', 
                                 adjusted.ref.period = F,
                                 quantile = 0.90)
  q0.90.m2 <- value.obs.sim.full(estaciones, data, ref.period,
                                 models.list, m2, type = 'quantile', 
                                 adjusted.ref.period = F,
                                 quantile = 0.90)

  q0.95.m1 <- value.obs.sim.full(estaciones, data, ref.period,
                                 models.list, m1, type = 'quantile', 
                                 adjusted.ref.period = F,
                                 quantile = 0.95)
  q0.95.m2 <- value.obs.sim.full(estaciones, data, ref.period,
                                 models.list, m2, type = 'quantile', 
                                 adjusted.ref.period = F,
                                 quantile = 0.95)

  q0.99.m1 <-value.obs.sim.full(estaciones, data, ref.period,
                                models.list, m1, type = 'quantile', 
                                adjusted.ref.period = F,
                                quantile = 0.99)
  q0.99.m2 <- value.obs.sim.full(estaciones, data, ref.period,
                                 models.list, m2, type = 'quantile', 
                                 adjusted.ref.period = F,
                                 quantile = 0.99)

  freq.m1 <- value.obs.sim.full(estaciones, data, ref.period,
                                models.list, m1, type = 'freq', 
                                adjusted.ref.period = F)
  freq.m2 <-value.obs.sim.full(estaciones, data, ref.period,
                               models.list, m2, type = 'freq', 
                               adjusted.ref.period = F)

  int.m1 <- value.obs.sim.full(estaciones, data, ref.period,
                               models.list, m1, type = 'intensity', 
                               adjusted.ref.period = F)
  int.m2 <- value.obs.sim.full(estaciones, data, ref.period,
                               models.list, m2, type = 'intensity', 
                               adjusted.ref.period = F)

  df <- data.frame(
    type = c('Amount', 'Frequency', 'Intensity', 'q0.90', 'q0.95', 'q0.99')
  )

  metrics <- c('SCC', 'RMSE', 'RB')
  for (i in metrics){
    df[[paste0(i, '.', m1)]] <- c(mean.metric(i, amounts.m1),
                                  mean.metric(i, freq.m1),
                                  mean.metric(i, int.m1),
                                  mean.metric(i, q0.90.m1),
                                  mean.metric(i, q0.95.m1),
                                  mean.metric(i, q0.99.m1))
    df[[paste0(i, '.', m2)]] <- c(mean.metric(i, amounts.m2),
                                  mean.metric(i, freq.m2),
                                  mean.metric(i, int.m2),
                                  mean.metric(i, q0.90.m2),
                                  mean.metric(i, q0.95.m2),
                                  mean.metric(i, q0.99.m2))
  }


  df <- df[, c(1, 2, 4, 6, 3, 5, 7)]
  return(df)
}

df.day.M6.M9 <- model.comp.df(estaciones, per.comun.day, X.MDQ, MDQ, 'M6', 'M9', n.sim = 100)
df.hour.M7.M8 <- model.comp.df(estaciones[-4], per.comun.h, X.MHQ, MHQ, 'M7', 'M8', n.sim = 100)

library(writexl)
write_xlsx(df.day.M6.M9, "df.day.M6.M9.xlsx")
write_xlsx(df.hour.M7.M8, "df.hour.M7.M8.xlsx")

#results for common models
hour.amount.common.MHQ <- value.obs.sim.full(estaciones, X.MHQ, per.comun.h,
                                            common.models.final, 'MHQ', type = 'mean',
                                            adjusted.ref.period = TRUE)
hour.q0.90.common.MHQ <- value.obs.sim.full(estaciones, X.MHQ, per.comun.h,
                                            common.models.final, 'MHQ', type = 'quantile', quantile = 0.90,
                                            adjusted.ref.period = TRUE)
hour.q0.95.common.MHQ <- value.obs.sim.full(estaciones, X.MHQ, per.comun.h,
                                            common.models.final, 'MHQ', type = 'quantile', quantile = 0.95,
                                            adjusted.ref.period = TRUE)
hour.q0.99.common.MHQ <- value.obs.sim.full(estaciones, X.MHQ, per.comun.h,
                                            common.models.final, 'MHQ', type = 'quantile', quantile = 0.99,
                                            adjusted.ref.period = TRUE)
hour.freq.common.MHQ <- value.obs.sim.full(estaciones, X.MHQ, per.comun.h,
                                            common.models.final, 'MHQ', type = 'freq',
                                            adjusted.ref.period = TRUE)
hour.int.common.MHQ <- value.obs.sim.full(estaciones, X.MHQ, per.comun.h,
                                          common.models.final, 'MHQ', type = 'intensity',
                                          adjusted.ref.period = TRUE)


round(mean.metric('SCC', hour.amount.common.MHQ), 5)
round(mean.metric('RMSE', hour.amount.common.MHQ), 5)
round(mean.metric('RB', hour.amount.common.MHQ), 5)



