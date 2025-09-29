# COMPARISON OF MODELS USING DIFFERENT METRICS: SCC, RMSE, RB
rm(list = ls())

# load data and functions
load('data.RData')
rm(list = setdiff(ls(), c('estaciones')))
MDQ <- readRDS('MDQ.rds') # models and data for each model
MHQ <- readRDS('MHQ.rds')
source('methods.R')

library(lubridate)
library(dplyr)
# period of reference (common period of all stations)
ref.period <- as.Date(c('12/05/2011', '30/11/2023'), format = '%d/%m/%Y')



# Anual mean precipitation in JJA
n <- year(ref.period[2]) - year(ref.period[1]) + 1
mean.obs.sim.station <- function(station, ref.period, data, model, 
                                 n.sim = 100, months = NULL){
  # observed values
  X <- data[[station]]$X
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = "%Y-%m-%d")
  x.obs <- X[, c('date', paste0(station, '.p'))]
  
  if (!is.null(months)){
    ind.ref <- which(x.obs$date >= ref.period[1] & x.obs$date <= ref.period[2]
                     & month(x.obs$date) %in% months)
  }else{
    ind.ref <- which(x.obs$date >= ref.period[1] & x.obs$date <= ref.period[2])
  }
  
  
  x.obs <- x.obs[ind.ref, ]
  
  x.obs.final <- x.obs %>% 
    group_by(year(date)) %>%
    summarise(mean.p = quantile(.data[[paste0(station, '.p')]], na.rm = TRUE, probs = 0.90))
  
  x.obs.final <- mean(x.obs.final[['mean.p']])
  #simulation values
  mod <- data[[station]][[model]] #desired model
  
  vector.sim <- c()
  for (i in 1:n.sim){
    shape <- 1 / mod$sigma.fv[ind.ref]^2
    rate <- shape / mod$mu.fv[ind.ref]
    y.sim <- rgamma(length(ind.ref), shape = shape, rate = rate)
    
    y.sim <- data.frame(
      date = x.obs$date,
      y.sim = y.sim
    )
    
    y.sim.final <- y.sim %>% 
      group_by(year(date)) %>%
      summarise(mean.p = mean(y.sim, na.rm = TRUE))
    
    y.sim.final <- mean(y.sim.final[['mean.p']])
    vector.sim[i] <- y.sim.final
  }
  
  
  return(c(x.obs.final, vector.sim))
}

mean.obs.sim.full <- function(estaciones, ref.period, data, model, 
                              n.sim = 100, months = NULL){
  N <- length(estaciones)
  df <- data.frame(matrix(NA, ncol = n.sim + 2, nrow = length(estaciones)))
  colnames(df) <- c('station', 'x.obs', paste0('y.sim.', 1:n.sim))
  df$station <- estaciones
  
  names.fill <- c('x.obs', paste0('y.sim.', 1:n.sim))
  
  for (i in 1:length(estaciones)){
    mean.station <- mean.obs.sim.station(estaciones[i],
                                         ref.period = ref.period,
                                         data = data,
                                         model = model,
                                         n.sim = n.sim,
                                         months = months)
    df[i, names.fill] <- mean.station
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
df.amounts.M5 <- mean.obs.sim.full(estaciones[-4], ref.period, MHQ, 'M4', n.sim = 100)
df.amounts.M6 <- mean.obs.sim.full(estaciones[-4], ref.period, MHQ, 'M6', n.sim = 100)

mean.metric('SCC', df.amounts.M6)  
mean.metric('SCC', df.amounts.M5)
mean.metric('RMSE', df.amounts.M6)
mean.metric('RMSE', df.amounts.M5)
mean.metric('RB', df.amounts.M6)
mean.metric('RB', df.amounts.M5)

#MDQ
df.amounts.M5 <- mean.obs.sim.full(estaciones, ref.period, MDQ, 'M5', n.sim = 100)
df.amounts.M6 <- mean.obs.sim.full(estaciones, ref.period, MDQ, 'M6', n.sim = 100)

mean.metric('SCC', df.amounts.M6)  
mean.metric('SCC', df.amounts.M5)
mean.metric('RMSE', df.amounts.M6)
mean.metric('RMSE', df.amounts.M5)
mean.metric('RB', df.amounts.M6)
mean.metric('RB', df.amounts.M5)
