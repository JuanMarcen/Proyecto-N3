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
value.obs.sim.station <- function(station, ref.period, data, model, 
                                 type, quantile = NULL,
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
  mod <- data[[station]][[model]] #desired model
  
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

a <- value.obs.sim.station(station, ref.period, MDQ, 'M6', type = 'freq', n.sim = 10)

value.obs.sim.full <- function(estaciones, ref.period, data, model, 
                               type, quantile = NULL,
                              n.sim = 100, months = NULL){
  N <- length(estaciones)
  df <- data.frame(matrix(NA, ncol = n.sim + 2, nrow = length(estaciones)))
  colnames(df) <- c('station', 'x.obs', paste0('y.sim.', 1:n.sim))
  df$station <- estaciones
  
  names.fill <- c('x.obs', paste0('y.sim.', 1:n.sim))
  
  for (i in 1:length(estaciones)){
    cat('Cálculo para la estación: ', estaciones[i], '\n')
    value.station <- value.obs.sim.station(estaciones[i],
                                         ref.period = ref.period,
                                         data = data,
                                         model = model,
                                         n.sim = n.sim,
                                         months = months,
                                         type = type, 
                                         quantile = quantile)
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
df.amounts.M5 <- value.obs.sim.full(estaciones[-4], ref.period, MHQ, 'M4', n.sim = 100)
df.amounts.M6 <- value.obs.sim.full(estaciones[-4], ref.period, MHQ, 'M6', n.sim = 100)

mean.metric('SCC', df.amounts.M6)  
mean.metric('SCC', df.amounts.M5)
mean.metric('RMSE', df.amounts.M6)
mean.metric('RMSE', df.amounts.M5)
mean.metric('RB', df.amounts.M6)
mean.metric('RB', df.amounts.M5)

#MDQ
df.amounts.M5 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M5', n.sim = 100, type = 'freq')
df.amounts.M6 <- value.obs.sim.full(estaciones, ref.period, MDQ, 'M6', n.sim = 100, type = 'freq')

mean.metric('SCC', df.amounts.M6)  
mean.metric('SCC', df.amounts.M5)
mean.metric('RMSE', df.amounts.M6)
mean.metric('RMSE', df.amounts.M5)
mean.metric('RB', df.amounts.M6)
mean.metric('RB', df.amounts.M5)

