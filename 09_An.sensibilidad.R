rm(list = setdiff(ls(), c('common.models.final', 'MDQ', 'MHQ', 'X.MDQ',
                          'X.MHQ', 'estaciones', 'per.comun.day',
                          'per.comun.h')))

# Sensitivity analysis
library(qs)
X.MDQ <- qread('X.MDQ.qs')
X.MHQ <- qread('X.MHQ.qs')

#----Phase 1:Simulation based on observed data ----
library(dplyr)
library(lubridate)
library(Polychrome)
cols <- createPalette(28, c("#FF0000", "#0000FF", "#00FF00", "#FFFF00",
                            '#FF00FF', '#FFA500', '#228B22', '#A020F0',
                            '#ADD8E6', '#EE82EE'))

df.q <- function(stations, data, model.list, model,
                 period, 
                 quantiles, month = NULL){
  
  df <- data.frame(
    matrix(NA, ncol = length(quantiles) + 1, nrow = length(stations))
  )
  names.q <- paste0('q', quantiles)
  names(df) <- c('station', names.q)
  df$station <- stations
  
  for (station in stations){
    X <- data[[station]]
    X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = '-'), 
                      format = '%Y-%m-%d')
    X <- X %>% 
      filter(date >= period[1] & date <= period[2])
    
    x.obs <- X[[paste0(station, '.p')]]
    
    if (!is.null(mes)){
      ind <- which(X$mes %in% mes)
    }else{
      ind <- 1:dim(X)[1]
    }
    
    q.obs <- quantile(x.obs[ind], probs = quantiles)
    names(q.obs) <- names.q
    
    df[which(stations == station), names.q] <- q.obs
  }
  
  return(df)
}

df.q.MDQ.JJA <- df.q(estaciones, X.MDQ, common.models.final, 'MDQ', 
                     per.comun.day, c(0.05, 0.50, 0.90, 0.95, 0.99),
                     month = c(6, 7, 8))

boxplot.q.sim <- function(station, data, model.list, model, period,
                          quantiles, n.sim, month = NULL, plot = FALSE){
  
  mod <- model.list[[station]][[model]]
  X <- data[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = '-'), 
                    format = '%Y-%m-%d')
  X <- X %>% 
    filter(date >= period[1] & date <= period[2])
  
  x.obs <- X[[paste0(station, '.p')]]
  
  if (!is.null(month)){
    ind <- which(X$mes %in% mes)
  }
  
  mu <- mod$mu.fv[ind]
  shape <- 1 / mod$sigma.fv[ind]^2
  rate <- shape / mu
  
  names.q <- paste0('q', quantiles)
  
  q.obs <- quantile(x.obs[ind], probs = quantiles)
  names(q.obs) <- names.q
  
  # simulation based on observed data (model fitted values of gamma distr.)
  q.sim <- data.frame(matrix(NA, ncol = length(q.obs)))
  for (i in 1:n.sim){
    u <- rgamma(length(x.obs[ind]), shape = shape, rate = rate)
    q <- quantile(u, probs = quantiles)
    names(q) <- names.q
    q.sim <- rbind(q.sim, q)
  }
  
  q.sim[1, ] <- q.obs
  
  if (plot == TRUE){
  boxplot(q.sim[-1, ],
          at = q.obs,               
          names = paste0(quantiles, '.obs'),  
          xlim = c(q.obs[1]-0.5, q.obs[5]+0.5),
          ylim = c(0, max(q.sim)),
          col = "lightblue",
          main = paste("Boxplots alineados con q.obs", station),
          ylab = "Valores simulados",
          xlab = "Cuantiles observados")
  lines(q.obs, q.obs, col = "red", pch = 19, cex = 1.3, type = 'b')
  }
  
  q.obs.matrix <- matrix(q.obs, nrow = n.sim, ncol = length(quantiles), byrow = T)
  prob.gr.obs <- q.sim[-1,] > q.obs.matrix
  prob.gr.obs <- colSums(prob.gr.obs) / n.sim
  
  return(prob.gr.obs)
}

set.seed(05052002)
boxplot.q.sim(estaciones[1], X.MDQ, common.models.final, 'MDQ', 
              per.comun.day, 
              quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
              n.sim = 1000, month = c(6,7,8), plot = TRUE)



df.q.sim.gr.obs <- function(stations, data, model.list, model, period,
                            quantiles, n.sim, month = NULL, plot = FALSE){
  df <- data.frame(
    matrix(NA, ncol = length(quantiles) + 1, nrow = length(stations))
  )
  names.q <- paste0('q', quantiles)
  names(df) <- c('station', names.q)
  df$station <- stations
  
  for (station in stations){
    df[which(stations == station), -1] <- boxplot.q.sim(station = station, 
                                                        data = data, 
                                                        model.list = model.list, 
                                                        model = model, 
                                                        period = period, 
                                                        quantiles = quantiles,
                                                        n.sim = n.sim,
                                                        month = month)
                                                    
  }
  
  return(df)
}
set.seed(05052002)
df.q.sim.gr.obs.MDQ.DJF <- df.q.sim.gr.obs(estaciones, 
                                           X.MDQ, common.models.final, 'MDQ', 
                                           per.comun.day, 
                                           quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
                                           n.sim = 100, month = c(12,1,2))

boxplot(df.q.sim.gr.obs.MDQ.JJA[, -1])
plot(1:5, df.q.sim.gr.obs.MDQ.JJA[1, -1], type = 'b', pch = 19,
     xaxt = 'n', xlab = 'Quantiles', ylab = 'P(q.sim > q.obs)',
     ylim = c(0,1), col = cols[1])
axis(1, at = 1:5, labels = paste0('q', c(0.05, 0.50, 0.90, 0.95, 0.99)))
for(i in 1:length(estaciones)){
  lines(1:5, df.q.sim.gr.obs.MDQ.JJA[i, -1], pch = 19, type = 'b', col = cols[i])
}
legend("topleft", legend = estaciones,
       col = cols, lwd = 2, ncol = 2)
