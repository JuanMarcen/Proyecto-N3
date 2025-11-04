rm(list = setdiff(ls(), c('common.models.final', 'MDQ', 'MHQ', 'X.MDQ',
                          'X.MHQ', 'estaciones', 'per.comun.day',
                          'per.comun.h')))

# Sensitivity analysis
library(qs)
X.MDQ <- qread('X.MDQ.qs')
X.MHQ <- qread('X.MHQ.qs')
#X.MDO <- qread('X.MDO.qs')
X.MHO <- qread('X.MHO.qs')

load('data.RData')

#----GENERATORS----
library(gamlss)
sanitize_formula <- function(f) {
  ftxt <- paste(deparse(f), collapse = " ")
  ftxt <- gsub("\\b[[:alnum:]_.]+\\$", "", ftxt)
  ftxt <- gsub("poly\\(\\s*([^,\\)]+)\\s*,\\s*([0-9]+)\\s*\\)",
               "poly(\\1, \\2, raw = TRUE)", ftxt, perl = TRUE)
  as.formula(ftxt, env = environment())
}
generator.qty.obs <- function(L, shape, rate, n.sim){
  y.sim <- data.frame(matrix(NA, nrow = L, ncol = n.sim))
  for (i in 1:n.sim){
    u <- rgamma(L, shape = shape, rate = rate)
    y.sim[, i] <- u
  }
  names(y.sim) <- paste0('y.sim.', 1:n.sim)
  
  return(y.sim)
}


#
generator.qty.oc <- function(station, data.mq, data.mo, 
                             model.list, model.q, model.o, 
                             period, n.sim, ocurrence = TRUE,
                             years = NULL){
  
  Xo <- data.mo[[station]]
  Xo$date <- as.Date(paste(Xo$t, Xo$mes, Xo$dia.mes, sep = "-"), 
                     format = '%Y-%m-%d')
  
  Xo <- Xo %>% filter(date >= period[1] & date <= period[2])
  
  Xq <- data.mq[[station]]      
  
  Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"), 
                     format = '%Y-%m-%d')
  Xq <- Xq %>% filter(date >= period[1] & date <= period[2])
  Xq <- Xq[, -ncol(Xq)]
  
  mq <- model.list[[station]][[model.q]]
  
  aux.mu.formula <- mq$mu.formula
  aux.mu.formula <- as.formula(
    paste(as.character(aux.mu.formula[2]), '~', 
          paste(labels(terms(aux.mu.formula)), collapse = '+'))
  )
  
  aux.sigma.formula <- mq$sigma.formula
  aux.sigma.formula <- as.formula(
    paste('~', as.character(aux.sigma.formula[2]))
  )
  
  
  mu.form <- sanitize_formula(aux.mu.formula)
  sigma.form <- sanitize_formula(aux.sigma.formula)
  
  mo <- model.list[[station]][[model.o]]
  
  #re ajuste modelo cantidad
  mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                family = GA, data = Xq, trace = FALSE)
  
  # por hecho que hay lag (sino mirar programa 06_generador.R)
  # nuevos mu y sigma necesarios para la primera fecha
  # filtrar por años si el valor es no nulo
  
  if (!is.null(years)){
    ind.year <- which(Xo$t %in% years)
  }else{
    ind.year <- nrow(Xo)
  }
  aux.Xo <- Xo[ind.year, names(Xq)]
  
  cat('Dia 1\n')
  cat('mu..')
  mu.new <- suppressWarnings(
    predict(mq2, newdata = aux.Xo[1, , drop = FALSE], what = 'mu', type = 'response',
            data = Xq)
  )
  cat('sigma..\n')
  sigma.new <- suppressWarnings(
    predict(mq2, newdata = aux.Xo[1, , drop = FALSE], what = 'sigma', type = 'response',
            data = Xq)
  )
  
  df.gen <- aux.Xo %>%
    select(t, l, dia.mes, mes)
  
  df.gen <- df.gen %>%
    mutate(
      prob.p = mo$fitted.values[ind.year], #fitted values of model MDO(phase2) or dynamically(phase3)
      mu.fv = NA,
      sigma.fv = NA,
      shape.fv = NA,
      rate.fv = NA,
      p.day.obs = aux.Xo[[paste0(station, '.p')]]
    )
  
  
  
  # parámetros para el día 1 (estos son vectores escalares porque representan "la estructura")
  df.gen[1, c(5:(ncol(df.gen) - 1))] <- c(
    prob.p = mo$fitted.values[1],
    mu.fv = mu.new,
    sigma.fv = sigma.new,
    shape.fv = 1 / sigma.new^2,
    rate.fv = 1 / sigma.new^2 / mu.new
  )
  
  # Inicializar matriz de simulaciones: filas = días, columnas = simulaciones
  n_days <- nrow(df.gen)
  sim_matrix <- matrix(NA_real_, nrow = n_days, ncol = n.sim)
  colnames(sim_matrix) <- paste0("p.day.sim.", 1:n.sim)
  
  # Simulación día 1 (para todas las n.sim)
  rain_vec <- runif(n.sim) <= df.gen$prob.p[1]
  gamma_vec <- rgamma(n.sim, shape = rep(df.gen$shape.fv[1], n.sim),
                      rate = rep(df.gen$rate.fv[1], n.sim))
  gamma_vec[!rain_vec] <- 0
  sim_matrix[1, ] <- gamma_vec
  
  # Índices de simulaciones activas (las que aún no han sido descartadas)
  active_sims <- seq_len(n.sim)
  
  # Umbral para considerar un valor "absurdo" — aquí uso el percentil 0.999 histórico
  umbral <- max(Xo[[paste0(station, '.p')]])
  # Si prefieres un valor fijo, reemplaza la línea anterior, p.ej. umbral <- 2000
  
  # Loop día a día
  for (i in 2:nrow(df.gen)) {
    cat('Dia ', i, ' \n')
    
    if (length(active_sims) == 0) {
      cat("No quedan simulaciones activas en el día ", i, ". Terminando.\n")
      break
    }
    
    # Preparamos aux con una fila repetida por cada simulación activa
    aux <- aux.Xo[i, , drop = FALSE]
    aux <- aux[rep(1, length(active_sims)), , drop = FALSE]
    
    # Lag = valor simulado el día anterior para cada simulación activa
    lag_values <- sim_matrix[i - 1, active_sims]
    
    # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
    bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > 500))
    if (length(bad_idx_rel) > 0) {
      # Mapear a índices absolutos en el conjunto original de simulaciones
      bad_idx_abs <- active_sims[bad_idx_rel]
      cat("  ⚠️ Detectadas simulaciones inválidas en día ", i - 1, ": ",
          paste0(bad_idx_abs, collapse = ", "), " -> descartando\n")
      # Dejamos NA explícito (ya lo es) en la fila pasada; pero por claridad:
      sim_matrix[i - 1, bad_idx_abs] <- NA_real_
      # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
      active_sims <- setdiff(active_sims, bad_idx_abs)
    }
    
    # Si tras eliminar no queda ninguna simulación válida, salimos
    if (length(active_sims) == 0) {
      cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
      break
    }
    
    # Reconstruir aux con la nueva longitud de simulaciones activas
    aux <- Xo[i, , drop = FALSE]
    aux <- aux[rep(1, length(active_sims)), , drop = FALSE]
    aux[[paste0(station, '.p.lag')]] <- as.numeric(sim_matrix[i - 1, active_sims])
    
    aux <- aux[, names(Xq)]
    aux2 <- aux
    
    # obtención de mu y sigma solo para simulaciones activas
    cat('mu..')
    aux2$mu.new <- suppressWarnings(
      predict(mq2, newdata = aux, what = 'mu', type = 'response',
              data = Xq)
    )
    cat('sigma..')
    aux2$sigma.new <- suppressWarnings(
      predict(mq2, newdata = aux, what = 'sigma', type = 'response',
              data = Xq)
    )
    
    # obtención de probabilidad de lluvia solo para simulaciones activas
    # cat('prob..\n')
    if (ocurrence == TRUE){
      aux2$prob.dia <- predict(mo, newdata = aux, type = 'response')
    }else{
      aux2$prob.dia <- df.gen$prob.p[i]
    }
    # Simular solo para las filas activas
    aux2 <- aux2 %>%
      select(prob.dia, mu.new, sigma.new) %>%
      rowwise() %>%
      mutate(
        shape.fv = 1 / sigma.new^2,
        rate.fv = 1 / sigma.new^2 / mu.new
      ) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(
        sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
      ) %>%
      ungroup()
    
    # Guardar las simulaciones activas en la matriz en la posición correspondiente
    sim_matrix[i, active_sims] <- as.numeric(aux2$sim)
    
    # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
  } # end for days
  
  # Convertir sim_matrix en data.frame con nombres adecuados y añadir a df.gen
  sim_df <- as.data.frame(sim_matrix)
  names(sim_df) <- paste0("p.day.sim.", 1:n.sim)
  df.gen <- cbind(df.gen, sim_df)
  
  #assign(paste0('df.gen.day.', station), df.gen)
  
  # rm(list = c('df.gen', 'sim_df', 'sim_matrix', 
  #             'rain_vec', 'gamma_vec', 'n_days', 'sigma.new',
  #             'mu.new', 'X_data', 'Xo', 'Xq', 'mq2',
  #             'mu.form', 'sigma.form', 'mo', 'mq'))
  
  cat(length(active_sims), ' han sobrevivido a la simulación')
  
  return(df.gen[, paste0("p.day.sim.", 1:n.sim)])
}


basura <- generator.qty.oc(estaciones[1], X.MHQ, X.MHO, common.models.final,
                            'MHQ', 'MHO', per.comun.h, n.sim = 100, ocurrence = T,
                            years = 2015)

basura.sim <- basura[[3]]
rain.yes <- apply(basura.sim, 2, function(x) sum(x > 0)/nrow(basura.sim))
station <- estaciones[1]
boxplot(rain.yes)

points(1, sum(common.models.final[[station]][['MDO']]$y > 0)/length(common.models.final[[station]][['MDO']]$y)
       , col = 'red',
       pch = 19)

#----SENSITIVITY ANALYSIS----
library(dplyr)
library(lubridate)
library(Polychrome)
cols <- createPalette(28, c("#FF0000", "#0000FF", "#00FF00", "#FFFF00",
                            '#FF00FF', '#FFA500', '#228B22', '#A020F0',
                            '#ADD8E6', '#EE82EE'))
#particion según el lag
partition.lag <- function(df){
  aux.df <- data.frame(matrix(NA, ncol = ncol(df), nrow = nrow(df)))
  colnames(aux.df) <- colnames(df)
  for (j in 1:ncol(df)){
    x <- df[, j]
    aux.x <- x
    for (i in 2:length(x)){
      if (x[i-1] == 0){
        aux.x[i] <- FALSE
      }else{
        aux.x[i] <- T
      }
    }
    aux.df[, j] <- as.logical(aux.x)
  }
  
  return(aux.df)
}
#indices que se corresponden con el cuarto de dias con mayor vaLor de la variable influyente
top.quarter.max.z.value <- function(mod, indx){
  z.values <- summary(mod)$coefficients[, 'z value']
  
  aux <- grep('z', names(z.values))
  max.z.value <- z.values[aux][which.max(abs(z.values[aux]))]
  var.name <- gsub("poly\\(([^,]+),.*", "\\1", names(max.z.value))
  
  data <- mod$data
  top.quarter <- which(data[indx, var.name] >= quantile(data[indx, var.name], 0.75))
  
  return(list(var.name = var.name, top.quarter = top.quarter))
}

# dataframe of observed quantiles
df.q <- function(stations, data, model.list, model,
                 period, 
                 quantiles, month = NULL,
                 boxplot = FALSE){
  
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
    
    if (!is.null(month)){
      ind <- which(X$mes %in% month)
    }else{
      ind <- 1:dim(X)[1]
    }
    
    q.obs <- quantile(x.obs[ind], probs = quantiles)
    names(q.obs) <- names.q
    
    df[which(stations == station), names.q] <- q.obs
  }
  
  if (boxplot == TRUE){
    boxplot(df[, -1], ylab = 'Precipitation',
            main = paste0('Observed quantiles months ', paste(month, collapse = '-')))
  }
  
  return(df)
}

par(mfrow = c(2, 2))
df.q(estaciones, X.MDQ, common.models.final, 'MDQ', 
     per.comun.day, c(0.05, 0.50, 0.90, 0.95, 0.99),
     month = c(6, 7, 8), boxplot = T)
df.q(estaciones, X.MDQ, common.models.final, 'MDQ', 
     per.comun.day, c(0.05, 0.50, 0.90, 0.95, 0.99),
     month = c(9, 10, 11), boxplot = T)
df.q(estaciones, X.MDQ, common.models.final, 'MDQ', 
     per.comun.day, c(0.05, 0.50, 0.90, 0.95, 0.99),
     month = c(12, 1, 2), boxplot = T)
df.q(estaciones, X.MDQ, common.models.final, 'MDQ', 
     per.comun.day, c(0.05, 0.50, 0.90, 0.95, 0.99),
     month = c(3, 4, 5), boxplot = T)


boxplot.q.sim <- function(station, data, model.list, model, period,
                          quantiles, n.sim, month = NULL, years = NULL, 
                          seasons = NULL, 
                          generator,
                          plot = FALSE, 
                          extreme.lag = FALSE, n.days = NULL,
                          partition.lag = FALSE,
                          top.quarter.z.value = FALSE){
  
  mod <- model.list[[station]][[model]]
  X <- data[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = '-'), 
                    format = '%Y-%m-%d')
  X <- X %>% 
    filter(date >= period[1] & date <= period[2])
  
  if (!is.null(years)){
    X <- X %>%
      filter(t %in% years)
  }
  
  x.obs <- X[[paste0(station, '.p')]]
  
  if (!is.null(seasons)){
    ind.list <- list(
      ind = 1:dim(X)[1],
      ind.jja = which(X$mes %in% c(6, 7, 8)),
      ind.son = which(X$mes %in% c(9, 10, 11)),
      ind.djf = which(X$mes %in% c(12, 1, 2)),
      ind.mam = which(X$mes %in% c(3, 4, 5))
    )
    
  }
  else if (!is.null(month)){
    ind.list <- list(ind = which(X$mes %in% month))
  }else{
    ind.list <- list(ind = 1:dim(X)[1])
  }
  
  ind <- ind.list[['ind']]
  mu <- mod$mu.fv[ind]
  shape <- 1 / mod$sigma.fv[ind]^2
  rate <- shape / mu
  
  names.q <- paste0('q', quantiles)
  
  # SIMULATIONS
  if(generator == 1){
    # como es ind, si season no es nulo, me hace todo
    y.sim <- generator.qty.obs(length(x.obs[ind]), shape, rate, n.sim = n.sim)
  }
  
  if (generator == 2){
    y.sim <- generator.qty.oc(station = station,
                              data.mo = X.MDO,
                              data.mq = X.MDQ,
                              model.list = model.list,
                              model.q = 'MDQ',
                              model.o = 'MDO',
                              period = period,
                              n.sim = n.sim,
                              ocurrence = FALSE)
    aux.marcador <- which(apply(y.sim, 2, function(x) sum(is.na(x))) > 0)
    if (length(aux.marcador) >= 1){
      y.sim <- y.sim[, -aux.marcador]
    }  
  }
  
  if (generator == 3){
    y.sim <- generator.qty.oc(station = station,
                              data.mo = X.MHO,
                              data.mq = X.MHQ,
                              model.list = model.list,
                              model.q = 'MHQ',
                              model.o = 'MHO',
                              period = period,
                              n.sim = n.sim, 
                              years = years)
    aux.marcador <- which(apply(y.sim, 2, function(x) sum(is.na(x))) > 0)
    if (length(aux.marcador) >= 1){
      y.sim <- y.sim[, -aux.marcador]
    }
    
  }
  
  mdo <- model.list[[station]][['MDO']]
  
  aux <- c()
  aux.lag.T <- c()
  aux.lag.F <- c()
  aux.top.quarter <- c()
  
  aux.equal <- c()
  aux.lag.T.equal <- c()
  aux.lag.F.equal <- c()
  aux.top.quarter.equal <- c()
  
  season.names <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
  
  i <- 1 #for season names in plot quantiles
  par(mfrow = c(3,2))
  for (indx in ind.list){
    q.obs <- quantile(x.obs[indx], probs = quantiles)
    names(q.obs) <- names.q
    #quantiles of simulation
    q.sim <- t(apply(y.sim[indx, ], 2, quantile, probs = quantiles, na.rm = T))

    # if (i == 1){
    #   par(mfrow = c(1,1))
    # }
    # if (i == 2){
    #   par(mfrow =c(2,2))
    # }
    
    if (plot == TRUE){
      if (!is.null(seasons)){
        title <- paste0("Simulated quantiles vs obs. quantiles ", 
                        station, '-', model, 
                        ' ', season.names[i])
      }
      else if (!is.null(month)){
        title <- paste0("Simulated quantiles vs obs. quantiles ", 
                         station, '-', model, 
                         ' (months ', paste0(month, collapse = '-'), ')')
      }
      else{
        title <- paste0("Simulated quantiles vs obs. quantiles ", 
                        station, '-', model)
      }
      boxplot(q.sim,
              at = q.obs,               
              names = paste0(quantiles, '.obs'),  
              xlim = c(q.obs[1]-0.5, q.obs[5]+0.5),
              ylim = c(0, max(q.sim)),
              col = "lightblue",
              main = title,
              ylab = "Simulated values",
              xlab = "Observed quantiles")
      lines(q.obs, q.obs, col = "red", pch = 19, cex = 1.3, type = 'b')
      i <- i + 1
    }
    
    q.obs.matrix <- matrix(q.obs, nrow = nrow(q.sim), ncol = length(quantiles), byrow = T)
    prob.gr.obs <- q.sim > q.obs.matrix
    prob.gr.obs <- colSums(prob.gr.obs) / (nrow(q.sim))
    
    aux <- c(aux, prob.gr.obs) #here it ends
    
    #equal vales
    prob.equal.obs <- q.sim == q.obs.matrix
    prob.equal.obs <- colSums(prob.equal.obs) / (nrow(q.sim))
    
    aux.equal <- c(aux.equal, prob.equal.obs)

    
    #extreme lag day
    if (extreme.lag == TRUE){
      aux.X <- X[indx, ]
      aux.ind <- order(aux.X[[paste0(station, '.p.lag')]], decreasing = T)[1:n.days]
      
      df2 <- data.frame(
        date = aux.X[aux.ind, 'date'],
        x.obs = x.obs[indx][aux.ind],
        x.lag.obs = aux.X[aux.ind , paste0(station, '.p.lag')]
      )
      
      df2 <- cbind(df2, y.sim[indx, ][aux.ind, ])
      # for (i in 1:n.sim){
      #   u <- rgamma(length(x.obs[ind][aux.ind]), shape = shape, rate = rate)
      #   y.sim <- cbind(y.sim, u)
      # }
      # colnames(y.sim)[4:ncol(y.sim)] <- paste0('y.sim.', 1:n.sim)
      
      if (plot == TRUE){
        if (!is.null(seasons)){
          title <- paste0("Extreme lag days ", 
                          station, '-', model, 
                          ' ', season.names[i - 1])
        }else if (!is.null(month)){
          title <- paste0("Extreme lag days ", 
                          station, '-', model, 
                          ' (months ', paste0(month, collapse = '-'), ')')
        }else{
          title <- paste0("Extreme lag days ", 
                          station, '-', model)
        }
        boxplot(t(df2[, -c(1, 2, 3)]), 
                names = df2$date, ylab = 'Precipitation',
                xlab = 'Date', 
                main = title)
        points(1: nrow(df2), df2$x.obs, 
               col = 'red', pch = 19)
      }
    }
    
    #partition of lags
    if (partition.lag == TRUE){
      aux.y.sim <- y.sim[indx, ] #simulacion en indices
      y.sim.part.lag <- partition.lag(aux.y.sim) #partition en esa simulacion
      
      # T and F partitions of simulations
      y.sim.lag.T <- aux.y.sim
      y.sim.lag.T[!y.sim.part.lag] <- NA
      
      y.sim.lag.F <- aux.y.sim
      y.sim.lag.F[y.sim.part.lag == T] <- NA
      
      #simulation quantiles
      q.sim.T <- t(apply(y.sim.lag.T, 2, quantile, probs = quantiles, na.rm = T))
      q.sim.F <- t(apply(y.sim.lag.F, 2, quantile, probs = quantiles, na.rm = T))
      
      #observed quantiles of partitions
      x.obs.part <- partition.lag(data.frame(x.obs[indx]))
      
      x.obs.lag.T <- x.obs[indx]
      x.obs.lag.T[!x.obs.part] <- NA
      
      x.obs.lag.F <- x.obs[indx]
      x.obs.lag.F[x.obs.part == T] <- NA
      
      q.obs.T <- quantile(x.obs.lag.T, probs = quantiles, na.rm = T)
      q.obs.F <- quantile(x.obs.lag.F, probs = quantiles, na.rm = T)
      
      if (plot == TRUE){
        if (!is.null(seasons)){
          title <- paste0("Simulated quantiles vs obs. quantiles ", 
                          station, '-', model, 
                          ' ', season.names[i - 1])
        }
        else if (!is.null(month)){
          title <- paste0("Simulated quantiles vs obs. quantiles ", 
                          station, '-', model, 
                          ' (months ', paste0(month, collapse = '-'), ')')
        }
        else{
          title <- paste0("Simulated quantiles vs obs. quantiles ", 
                          station, '-', model)
        }
        
        boxplot(q.sim.T,
                at = q.obs.T,               
                names = paste0(quantiles, '.obs'),  
                xlim = c(q.obs.T[1]-0.5, q.obs.T[5]+0.5),
                ylim = c(0, max(q.sim.T)),
                col = "lightblue",
                main = paste(title, '(lag == TRUE)'),
                ylab = "Simulated values",
                xlab = "Observed quantiles")
        lines(q.obs.T, q.obs.T, col = "red", pch = 19, cex = 1.3, type = 'b')
        
        boxplot(q.sim.F,
                at = q.obs.F,               
                names = paste0(quantiles, '.obs'),  
                xlim = c(q.obs.F[1]-0.5, q.obs.F[5]+0.5),
                ylim = c(0, max(q.sim.F)),
                col = "lightblue",
                main = paste(title, '(lag == FALSE)'),
                ylab = "Simulated values",
                xlab = "Observed quantiles")
        lines(q.obs.F, q.obs.F, col = "red", pch = 19, cex = 1.3, type = 'b')
      }
      
      #df3
      q.obs.matrix.lag.T <- matrix(q.obs.T, nrow = nrow(q.sim.T), ncol = length(quantiles), byrow = T)
      prob.gr.obs.lag.T <- q.sim.T > q.obs.matrix.lag.T
      prob.gr.obs.lag.T <- colSums(prob.gr.obs.lag.T) / (nrow(q.sim.T))
      
      aux.lag.T <- c(aux.lag.T, prob.gr.obs.lag.T)
      
      #df3.equal
      prob.equal.obs.lag.T <- q.sim.T == q.obs.matrix.lag.T
      prob.equal.obs.lag.T <- colSums(prob.equal.obs.lag.T) / (nrow(q.sim.T))
      
      aux.lag.T.equal <- c(aux.lag.T.equal, prob.equal.obs.lag.T)
      
      #df4
      q.obs.matrix.lag.F <- matrix(q.obs.F, nrow = nrow(q.sim.F), ncol = length(quantiles), byrow = T)
      prob.gr.obs.lag.F <- q.sim.F > q.obs.matrix.lag.F
      prob.gr.obs.lag.F <- colSums(prob.gr.obs.lag.F) / (nrow(q.sim.F))
      
      aux.lag.F <- c(aux.lag.F, prob.gr.obs.lag.F)
      
      #df4.equal
      prob.equal.obs.lag.F <- q.sim.F == q.obs.matrix.lag.F
      prob.equal.obs.lag.F <- colSums(prob.equal.obs.lag.F) / (nrow(q.sim.F))
      
      aux.lag.F.equal <- c(aux.lag.F.equal, prob.equal.obs.lag.F)
    }
    
    if (top.quarter.z.value == TRUE){
      aux.y.sim <- y.sim[indx, ]#ismulacion en indices
      # indices de este subcojunto
      var.name <- top.quarter.max.z.value(mdo, indx)[[1]]
      top.quarter <- top.quarter.max.z.value(mdo, indx)[[2]] 
      
      #guardado de simulaciones del top quarter
      y.sim.top.quarter <- aux.y.sim[top.quarter, ]
      q.sim.top.quarter <- t(apply(y.sim.top.quarter, 2, 
                                   quantile, probs = quantiles, na.rm = T))
      x.obs.top.quarter <- x.obs[indx][top.quarter]
      q.obs.top.quarter <- quantile(x.obs.top.quarter, probs = quantiles, na.rm = T)
      
      if (plot == TRUE){
        if (!is.null(seasons)){
          title <- paste0("Simulated quantiles vs obs. quantiles ", 
                          station, '-', model, 
                          ' ', season.names[i - 1])
        }
        else if (!is.null(month)){
          title <- paste0("Simulated quantiles vs obs. quantiles ", 
                          station, '-', model, 
                          ' (months ', paste0(month, collapse = '-'), ')')
        }
        else{
          title <- paste0("Simulated quantiles vs obs. quantiles ", 
                          station, '-', model)
        }
        boxplot(q.sim.top.quarter,
                at = q.obs.top.quarter,               
                names = paste0(quantiles, '.obs'),  
                xlim = c(q.obs.top.quarter[1]-0.5, q.obs.top.quarter[5]+0.5),
                ylim = c(0, max(q.sim.top.quarter)),
                col = "lightblue",
                main = paste(title, '- Top quarter', var.name),
                ylab = "Simulated values",
                xlab = "Observed quantiles")
        lines(q.obs.top.quarter, q.obs.top.quarter, col = "red", pch = 19, cex = 1.3, type = 'b')
        
        
      }
      
      #df 5
      q.obs.matrix.top.quarter <- matrix(q.obs.top.quarter, nrow = nrow(q.sim.top.quarter),
                                         ncol = length(quantiles), byrow = T)
      prob.gr.obs.top.quarter <- q.sim.top.quarter > q.obs.matrix.top.quarter
      prob.gr.obs.top.quarter <- colSums(prob.gr.obs.top.quarter) / (nrow(q.sim.top.quarter))
      
      aux.top.quarter <- c(aux.top.quarter, prob.gr.obs.top.quarter)

      #df5.equal
      prob.equal.obs.top.quarter <- q.sim.top.quarter == q.obs.matrix.top.quarter
      prob.equal.obs.top.quarter <- colSums(prob.equal.obs.top.quarter) / (nrow(q.sim.top.quarter))
      
      aux.top.quarter.equal <- c(aux.top.quarter.equal, prob.equal.obs.top.quarter)
      
      if (extreme.lag == TRUE){
        aux.X <- X[indx, ][top.quarter, ]
        aux.ind <- order(aux.X[[paste0(station, '.p.lag')]], decreasing = T)[1:n.days]
        
        df2.top.quarter <- data.frame(
          date = aux.X[aux.ind, 'date'],
          x.obs = x.obs[indx][top.quarter][aux.ind],
          x.lag.obs = aux.X[aux.ind , paste0(station, '.p.lag')]
        )
        
        df2.top.quarter <- cbind(df2.top.quarter, y.sim[indx, ][top.quarter, ][aux.ind, ])
        
        if (plot == TRUE){
          if (!is.null(seasons)){
            title <- paste0("Extreme lag days ", 
                            station, '-', model, 
                            ' ', season.names[i - 1])
          }else if (!is.null(month)){
            title <- paste0("Extreme lag days ", 
                            station, '-', model, 
                            ' (months ', paste0(month, collapse = '-'), ')')
          }else{
            title <- paste0("Extreme lag days ", 
                            station, '-', model)
          }
          boxplot(t(df2.top.quarter[, -c(1, 2, 3)]), 
                  names = df2.top.quarter$date, ylab = 'Precipitation',
                  xlab = 'Date', 
                  main = paste(title, 'Top quarter', var.name))
          points(1: nrow(df2.top.quarter), df2.top.quarter$x.obs, 
                 col = 'red', pch = 19)
        }
      }
    }
  }
  
  df <- data.frame(matrix(aux, ncol = length(quantiles), byrow = T))
  df.equal <- data.frame(matrix(aux.equal, ncol = length(quantiles), byrow = T))
  if(!is.null(seasons)){
    rownames(df) <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
  }
  colnames(df) <- names.q
  
  if (partition.lag == TRUE){
    df3 <- data.frame(matrix(aux.lag.T, ncol = length(quantiles), byrow = T))
    df4 <- data.frame(matrix(aux.lag.F, ncol = length(quantiles), byrow = T))
    
    df3.equal <- data.frame(matrix(aux.lag.T.equal, ncol = length(quantiles), byrow = T))
    df4.equal <- data.frame(matrix(aux.lag.F.equal, ncol = length(quantiles), byrow = T))
    if(!is.null(seasons)){
      rownames(df3) <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
      rownames(df4) <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
      rownames(df3.equal) <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
      rownames(df4.equal) <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
    }
    colnames(df3) <- names.q
    colnames(df4) <- names.q
    colnames(df3.equal) <- names.q
    colnames(df4.equal) <- names.q
  }
  
  if (top.quarter.z.value == TRUE){
    df5 <- data.frame(matrix(aux.top.quarter, ncol = length(quantiles), byrow = T))
    df5.equal <- data.frame(matrix(aux.top.quarter.equal, ncol = length(quantiles), byrow = T))
    if(!is.null(seasons)){
      rownames(df5) <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
      rownames(df5.equal) <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
    }
    colnames(df5) <- names.q
    colnames(df5.equal) <- names.q
    
  }
  
  if(extreme.lag == T & partition.lag == FALSE){
    df <- list(
      prob.gr.q.obs = df,
      extreme.lag = df2,
      y.sim = y.sim
    )
  }
  else if (extreme.lag == F & partition.lag == T){
    df <- list(
      prob.gr.q.obs = df,
      prob.equal.q.obs = df.equal,
      y.sim = y.sim,
      prob.gr.q.obs.lag.T = df3,
      prob.gr.q.obs.lag.F = df4,
      prob.equal.q.obs.lag.T = df3.equal,
      prob.equal.q.obs.lag.F = df4.equal)
  }
  else if (extreme.lag == T & partition.lag == T){
    df <- list(
      prob.gr.q.obs = df,
      prob.equal.q.obs = df.equal,
      y.sim = y.sim,
      extreme.lag = df2,
      prob.gr.q.obs.lag.T = df3,
      prob.gr.q.obs.lag.F = df4,
      prob.equal.q.obs.lag.T = df3.equal,
      prob.equal.q.obs.lag.F = df4.equal)
  }
  else{
    df <- list(
      prob.gr.q.obs = df,
      prob.equal.q.obs = df.equal,
      y.sim = y.sim
    )
  }
  
  if (top.quarter.z.value == TRUE){
    df[['prob.gr.obs.top.quarter']] <- df5
    df[['prob.equal.obs.top.quarter']] <- df5.equal
    if(extreme.lag == TRUE){
      df[['extreme.lag.top.quarter']] <- df2.top.quarter
    }
  }
  
  return(df)
}

set.seed(05052002)
# par(mfrow= c(2,2))
basura <- boxplot.q.sim(estaciones[estaciones == 'A126'], X.MHO, common.models.final, 'MHQ', 
              per.comun.h, 
              quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
              n.sim = 100, month = NULL, years = 2015,
              seasons = TRUE,
              generator = 3,
              plot = TRUE, 
              extreme.lag = TRUE,
              n.days = 8, 
              partition.lag = TRUE,
              top.quarter.z.value = F)

boxplot.q.sim(estaciones[1], X.MDQ, common.models.final, 'MDQ', 
              per.comun.day, 
              quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
              n.sim = 100, month = c(9,10, 11), plot = TRUE)
boxplot.q.sim(estaciones[1], X.MDQ, common.models.final, 'MDQ', 
              per.comun.day, 
              quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
              n.sim = 100, month = c(12,1,2), plot = TRUE)
boxplot.q.sim(estaciones[1], X.MDQ, common.models.final, 'MDQ', 
              per.comun.day, 
              quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
              n.sim = 100, month = c(3,4,5), plot = TRUE)



df.q.sim.gr.obs <- function(stations, data, model.list, model, period,
                            quantiles, n.sim, month = NULL, 
                            seasons = NULL, generator, 
                            plot = FALSE,
                            cols = NULL){
  
  if (!is.null(seasons)){
    n <- length(quantiles) * 5
    names.q <- paste0('q', quantiles)
    names.q <- as.vector(outer(names.q, 
                    c('ALL', 'JJA', 'SON', 'DJF', 'MAM'), paste, sep = "."))
  }else{
    n <- length(quantiles)
    names.q <- paste0('q', quantiles)
  }
  df <- data.frame(
    matrix(NA, ncol = n + 1, nrow = length(stations))
  )
  
  names(df) <- c('station', names.q)
  df$station <- stations
  
  for (station in stations){
     aux.df <- boxplot.q.sim(station = station, 
                            data = data, 
                            model.list = model.list, 
                            model = model, 
                            period = period, 
                            quantiles = quantiles,
                            n.sim = n.sim,
                            month = month,
                            seasons = seasons,
                            generator = generator)
  
  
     v <- as.vector(t(aux.df[[1]]))
     names(v) <- as.vector(outer(colnames(aux.df[[1]]), 
                                   rownames(aux.df[[1]]), paste, sep = "."))
     df[which(stations == station), -1] <- v
                                                    
  }
  
  
  
  if (plot == TRUE){
    if(is.null(seasons)){
      boxplot(df[, -1], xlab = 'Quantiles', ylab = 'P(q.sim > q.obs)',
              main = paste('P(q.sim > q.obs) months',
                           paste0(month, collapse = '-'), '-', model))
      # plot(1:5, df[1, -1], type = 'b', pch = 19,
      #      xaxt = 'n', xlab = 'Quantiles', ylab = 'P(q.sim > q.obs)',
      #      ylim = c(0,1), col = cols[1])
      #axis(1, at = 1:5, labels = paste0('q', c(0.05, 0.50, 0.90, 0.95, 0.99)))
      for(i in 1:length(stations)){
        lines(1:5, df[i, -1], pch = 19, type = 'b', col = cols[i])
      }
      abline(h = 0.25, lty = 2, lwd = 2)
      abline(h = 0.75, lty = 2, lwd = 2)
      legend("topleft", legend = stations,
             col = cols, lwd = 2, ncol = 2, cex = 0.5)
    }else{
      par(mfrow = c(1,1))
      boxplot(df[, -1], xlab = 'Quantiles', ylab = 'P(q.sim > q.obs)',
              main = paste('P(q.sim > q.obs) by seasons', '-', model),
              cex.axis = 0.7)
      # plot(1:5, df[1, -1], type = 'b', pch = 19,
      #      xaxt = 'n', xlab = 'Quantiles', ylab = 'P(q.sim > q.obs)',
      #      ylim = c(0,1), col = cols[1])
      #axis(1, at = 1:5, labels = paste0('q', c(0.05, 0.50, 0.90, 0.95, 0.99)))
      # for(i in 1:length(stations)){
      #   lines(1:25, df[i, -1], pch = 19, type = 'b', col = cols[i])
      # }
      abline(h = 0.25, lty = 2, lwd = 2)
      abline(h = 0.75, lty = 2, lwd = 2)
      # legend("topleft", legend = stations,
      #        col = cols, lwd = 2, ncol = 2, cex = 0.5)
    }
    
    
  }
  
  return(df)
}

set.seed(05052002)
#the only one that correpsonds to the boxplot plots is the first (due to the seed)
df.q.sim.gr.obs.gen.3.ALL <- df.q.sim.gr.obs(estaciones, 
                                             X.MDO, common.models.final, 'MDQ', 
                                             per.comun.day, 
                                             quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
                                             n.sim = 100, month = NULL,
                                             seasons = T,
                                             generator = 3,
                                             plot = T, cols = cols
                                             )


df.q.sim.gr.obs(estaciones, 
                X.MDQ, common.models.final, 'MDQ', 
                per.comun.day, 
                quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
                n.sim = 100, month = c(9, 10, 11),
                plot = T, cols = cols)
df.q.sim.gr.obs(estaciones, 
                X.MDQ, common.models.final, 'MDQ', 
                per.comun.day, 
                quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
                n.sim = 100, month = c(12, 1, 2),
                plot = T, cols = cols)
df.q.sim.gr.obs(estaciones, 
                X.MDQ, common.models.final, 'MDQ', 
                per.comun.day, 
                quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
                n.sim = 100, month = c(3, 4, 5),
                plot = T, cols = cols)



extreme.lag.sim <- function(station, data, model.list, model, period,
                            n.days, n.sim, month = NULL, plot = FALSE){
  mod <- model.list[[station]][[model]]
  X <- data[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = '-'), 
                    format = '%Y-%m-%d')
  X <- X %>% 
    filter(date >= period[1] & date <= period[2])
  
  x.obs <- X[[paste0(station, '.p')]]
  
  if (!is.null(month)){
    ind <- which(X$mes %in% month)
  }else{
    ind <- 1:dim(X)[1]
  }
  
  # observed data with greatest lags
  aux.X <- X[ind, ]
  aux.ind <- order(aux.X[[paste0(station, '.p.lag')]], decreasing = T)[1:n.days]
  
  # simulate only in those days
  mu <- mod$mu.fv[ind][aux.ind]
  shape <- 1 / mod$sigma.fv[ind][aux.ind]^2
  rate <- shape / mu
  
  # simulation based on observed data (model fitted values of gamma distr.)
  y.sim <- data.frame(
    date = aux.X[aux.ind, 'date'],
    x.obs = x.obs[ind][aux.ind],
    x.lag.obs = aux.X[aux.ind , paste0(station, '.p.lag')]
    )
  
  y.sim <- cbind(y.sim, generator.qty.obs(length(x.obs[ind][aux.ind]), 
                                          shape, rate, n.sim = n.sim))
  # for (i in 1:n.sim){
  #   u <- rgamma(length(x.obs[ind][aux.ind]), shape = shape, rate = rate)
  #   y.sim <- cbind(y.sim, u)
  # }
  # colnames(y.sim)[4:ncol(y.sim)] <- paste0('y.sim.', 1:n.sim)
  
  if (plot == TRUE){
    boxplot(t(y.sim[, -c(1, 2, 3)]), 
            names = y.sim$date, ylab = 'Precipitation',
            xlab = 'Date', 
            main = paste0('Extreme lag days ', station, ' (months', 
                         paste0(month, collapse = '-'), ') - ', model))
    points(1: nrow(y.sim), y.sim$x.obs, 
           col = 'red', pch = 19)
  }
  
  return(y.sim)
}

set.seed(05052002)
extreme.lag.sim(estaciones[1], X.MDQ, common.models.final, 
               'MDQ', per.comun.day, n.days = 8,
               n.sim = 100, 
               month = c(6, 7, 8), 
               plot = T)
extreme.lag.sim(estaciones[1], X.MDQ, common.models.final, 
                'MDQ', per.comun.day, n.days = 8,
                n.sim = 100, 
                month = c(9, 10, 11), 
                plot = T)
extreme.lag.sim(estaciones[1], X.MDQ, common.models.final, 
                'MDQ', per.comun.day, n.days = 8,
                n.sim = 100, 
                month = c(12, 1, 2), 
                plot = T)
extreme.lag.sim(estaciones[1], X.MDQ, common.models.final, 
                'MDQ', per.comun.day, n.days = 8,
                n.sim = 100, 
                month = c(3, 4, 5), 
                plot = T)





# EXTRA PARTITIONS----

# basura <- generator.qty.oc(station = estaciones[1],
#                            data.mq = X.MDQ,
#                            data.mo = X.MDO,
#                            model.list = common.models.final,
#                            model.q = 'MDQ',
#                            model.o = 'MDO',
#                            period = per.comun.day,
#                            n.sim = 100,
#                            ocurrence = TRUE)

basura.sim <- basura[[3]]

#partition of simulation with lag 0 or >0
partition.lag <- function(df){
  aux.df <- data.frame(matrix(NA, ncol = ncol(df), nrow = nrow(df)))
  colnames(aux.df) <- colnames(df)
  for (j in 1:ncol(df)){
    x <- df[, j]
    aux.x <- x
    for (i in 2:length(x)){
      if (x[i-1] == 0){
        aux.x[i] <- FALSE
      }else{
        aux.x[i] <- T
      }
    }
    aux.df[, j] <- as.logical(aux.x)
  }
  
  return(aux.df)
}


basura.sim.part.lag <- partition.lag(basura.sim)

basura.sim.lag.T <- basura.sim
basura.sim.lag.T[!basura.sim.part.lag] <- NA

basura.sim.lag.F <- basura.sim
basura.sim.lag.F[basura.sim.part.lag == T] <- NA
#substract first day of partitions
basura.sim.lag.F <- basura.sim.lag.F[-1, ]
basura.sim.lag.T <- basura.sim.lag.T[-1, ]

#quantiles of partitions
quantiles <- c(0.05, 0.5, 0.90, 0.95, 0.99)

q.sim.T <- t(apply(basura.sim.lag.T, 2, quantile, probs = quantiles, na.rm = T))
q.sim.F <- t(apply(basura.sim.lag.F, 2, quantile, probs = quantiles, na.rm = T))

# observed quantiles (in function up)
data <- X.MDO
X <- data[[station]]
period <- per.comun.day
X <- data[[station]]
X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = '-'),                                      format = '%Y-%m-%d')
X <- X %>%     filter(date >= period[1] & date <= period[2])
x.obs <- X[[paste0(station, '.p')]]

#partition of observed quantiles
x.obs.part <- partition.lag(data.frame(x.obs))

x.obs.lag.T <- x.obs
x.obs.lag.T[!x.obs.part] <- NA

x.obs.lag.F <- x.obs
x.obs.lag.F[x.obs.part == T] <- NA

q.obs.T <- quantile(x.obs.lag.T, probs = quantiles, na.rm = T)
q.obs.F <- quantile(x.obs.lag.F, probs = quantiles, na.rm = T)

# if plot = T
boxplot(q.sim.F,
        at = q.obs.F,               
        names = paste0(quantiles, '.obs'),  
        xlim = c(q.obs.F[1]-0.5, q.obs.F[5]+0.5),
        ylim = c(0, max(q.sim.F)),
        col = "lightblue",
        main = 'Hola',
        ylab = "Simulated values",
        xlab = "Observed quantiles")
lines(q.obs.F, q.obs.F, col = "red", pch = 19, cex = 1.3, type = 'b')


# studying the Z values of the models of ocurrence 
#(which seems to be causing the problems)
# becuase of the graphs of extreme values
top.quarter.max.z.value <- function(mod, indx){
  z.values <- summary(mod)$coefficients[, 'z value']
  
  aux <- grep('z', names(z.values))
  max.z.value <- z.values[aux][which.max(abs(z.values[aux]))]
  var.name <- gsub("poly\\(([^,]+),.*", "\\1", names(max.z.value))
  
  data <- mod$data
  top.quarter <- which(data[indx, var.name] >= quantile(data[indx, var.name], 0.75))
  
  return(top.quarter)
}


top.quarter.2 <- top.quarter.max.z.value(m, indx = ind.list$ind.djf)

y.sim <- basura[[2]]
y.sim <- y.sim[ind.list$ind.djf, ]
ind <- top.quarter.2[[2]]
y.sim.top.quarter <- y.sim[ind, ]
q.sim.top.quarter <- t(apply(y.sim.top.quarter, 2, 
                             quantile, probs = quantiles, na.rm = T))
x.obs.top.quarter <- x.obs[ind.list$ind.djf][ind]
q.obs.top.quarter <- quantile(x.obs.top.quarter, probs = quantiles, na.rm = T)

boxplot(q.sim.top.quarter,
        at = q.obs.top.quarter,               
        names = paste0(quantiles, '.obs'),  
        xlim = c(q.obs.top.quarter[1]-0.5, q.obs.top.quarter[5]+0.5),
        ylim = c(0, max(q.sim.top.quarter)),
        col = "lightblue",
        main = 'Hola',
        ylab = "Simulated values",
        xlab = "Observed quantiles")
lines(q.obs.top.quarter, q.obs.top.quarter, col = "red", pch = 19, cex = 1.3, type = 'b')


q.obs.matrix.top.quarter <- matrix(q.obs.top.quarter, nrow = nrow(q.sim.top.quarter), 
                                   ncol = length(quantiles), byrow = T)
prob.gr.obs.top.quarter <- q.sim.top.quarter > q.obs.matrix.top.quarter
prob.gr.obs.top.quarter <- colSums(prob.gr.obs.top.quarter) / (nrow(q.sim.top.quarter))

# iguales ?
prob.gr.obs.top.quarter.2 <- q.sim.top.quarter == q.obs.matrix.top.quarter
prob.gr.obs.top.quarter.2 <- colSums(prob.gr.obs.top.quarter.2) / (nrow(q.sim.top.quarter))


## EXTRA PARA EL HORARIO (AUX.XO YA EJECUTADO DENTRO DEL GENERADOR)
quantile(aux.Xo[aux.Xo$mes %in% c(9,10,11), 'EM71.p'], probs = 0.99)


#----GENERADOR HORARIO DEF----
# AÑADIDO LA CORRECCIÓN DIARIA
# TODOS DÍAS LLUEVE
# DÍAS DONDE NO LLUEVE NO LLUEVE
#CONSTRUCCION DATA.FRAME DATOS correctos!
library(dplyr)
cs <- function(t,harmonics=1, total) {
  # if(min(t) <0 | max(t) > 1){ stop(" t must be in [0,1] range")}
  if(min(harmonics) <1){stop("harmonics > = 1")}
  ret <- numeric(0)
  for ( i in harmonics) {
    ret <- cbind( ret, cos(2*pi*i*t/total), sin(2*pi*i*t/total))
  }
  if (missing(harmonics)) cnames <- c('c','s')
  else {
    cnames <- paste( c("c","s"), rep(harmonics, each = 2),sep=".")
  }
  colnames(ret) <- cnames
  ret
}

global_df <- readRDS('global_df.rds')
global_df <- global_df[, -which(colnames(global_df) %in% c('zg300.', 'zg500.', 'zg700.', 
                                                           'zt300.', 'zt500.', 'zt700.'))]

standardization <- function(x, sum.to, correction = FALSE){
  if (correction == FALSE){
    factor <- sum.to / sum(x)
    x <- x * factor
  }else{
    x[x <= 0.2 & x > 0] <- 0.2
    remainder <- sum.to - sum(x[x == 0.2])
    factor <- remainder/sum(x[x > 0.2])
    x[x > 0.2] <- x[x > 0.2] * factor
  }
  
  return(x)
}

RAIN.GENERATOR.og <- function(station, data.h, data.day, period,
                             models.list, data.mq, mo, mq, 
                             ocurrence = TRUE, years = NULL, 
                             n.sim = 100, type){
  
  if (type == 'hour'){
    #calculation of harmonics
    l <- 1:365
    h <- 0:23
    harm_l <- data.frame(l = l,
                         cs(l, harmonics = 1:4, 365))
    colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
    harm_h <- data.frame(h = h,
                         cs(0:23, harmonics = 1:4, 24))
    colnames(harm_h)[2:ncol(harm_h)] <- paste0(colnames(harm_h)[2:ncol(harm_h)], '.h')
    
    #construction of dataframe where we compute the simulations
    data.h <- data.h %>%
      left_join(harm_l, by = 'l') %>%
      left_join(harm_h, by = 'h')
    
    station.p <- paste0(station, '.p')
    
    data.h <- data.h[, c('t', 'l', 'mes', 'dia.mes', 'h', paste0(station, '.p'),
                         colnames(harm_l)[2:ncol(harm_l)],
                         colnames(harm_h)[2:ncol(harm_h)])]
    
    data.day <- data.day[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p') )]
    
    colnames(data.day)[colnames(data.day) == station.p] <- paste0(station.p, '.day')
    
    data.h <- data.h %>%
      left_join(data.day, by = c('t', 'l', 'mes', 'dia.mes'))
    
    data.h <- data.h %>% 
      mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]]))
    
    data.h <- data.h %>%
      mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0))%>% 
      relocate(Y, .after = !!station.p)
    
    data.h$date <- as.Date(paste(data.h$t, data.h$mes, data.h$dia.mes, sep = "-"), 
                           format = '%Y-%m-%d')
    
    data.h <- data.h %>%
      mutate(
        gr.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] >= 2),
        less.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] < 2)
      )
    
    data.h <- data.h %>% filter(date >= period[1] & date <= period[2])
    
    Xq <- data.mq[[station]]
    
    Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"),
                       format = '%Y-%m-%d')
    Xq <- Xq %>% filter(date >= period[1] & date <= period[2])
    Xq <- Xq[, -ncol(Xq)]
    
    # quantity model
    mq <- models.list[[station]][[mq]]
    
    # formulas rewritten for new dataframe
    aux.mu.formula <- mq$mu.formula
    aux.mu.formula <- as.formula(
      paste(as.character(aux.mu.formula[2]), '~', 
            paste(labels(terms(aux.mu.formula)), collapse = '+'))
    )
    
    aux.sigma.formula <- mq$sigma.formula
    aux.sigma.formula <- as.formula(
      paste('~', as.character(aux.sigma.formula[2]))
    )
    
    
    mu.form <- sanitize_formula(aux.mu.formula)
    sigma.form <- sanitize_formula(aux.sigma.formula)
    
    #ocurrence model
    mo <- models.list[[station]][[mo]]
    
    #re ajuste modelo cantidad (debido a las nuevas formulas)
    mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                  family = GA, data = Xq, trace = FALSE)
    
    # filtrar por años si el valor es no nulo
    if (!is.null(years)){
      ind.year <- which(data.h$t %in% years)
    }else{
      ind.year <- nrow(data.h)
    }
    aux.data.h <- data.h[ind.year, ]
    
    #number of days with rain
    n.days.rain <- length(
      unique(
        aux.data.h[aux.data.h[[paste0(station.p, '.day')]] > 0, 'date']
      )
    )
    
    aux.data.h <- data.h[ind.year, ] # for correct computation of prediction
    
    
    # INICIO DE LA SIMULACION
    # MATRIZ DE SIMULACIONES
    n_days <- nrow(aux.data.h)
    #n.sim <- 100
    sim_matrix <- matrix(NA_real_, nrow = n_days, ncol = n.sim)
    colnames(sim_matrix) <- paste0("p.day.sim.", 1:n.sim)
    
    # relleno de 0's para días sin lluvia observada
    ind.no.rain <- is.element(aux.data.h[[paste0(station.p, '.day')]], 0)
    sim_matrix[ind.no.rain, ] <- 0
    
    # Umbral para considerar un valor "absurdo" — aquí uso el percentil 0.999 histórico
    umbral <- max(aux.data.h[[paste0(station, '.p')]], na.rm = T)
    
    # hourly simulation for days with rain
    for (day in 1:n.days.rain){
      # la matriz se va a ir acutalizando. 
      # primera hora del primer dia con lluvia observada
      aux <- apply(sim_matrix, 1, function(row) all(is.na(row)))
      first.hour <- which(aux)[1]
      first.day <- first.hour:(first.hour + 23)
      
      #luego cunado se vuelva a llamar se corresponderá con el segundo dia
      
      # Índices de simulaciones activas (las que aún no han sido descartadas)
      active_sims <- seq_len(n.sim)
      
      # Si prefieres un valor fijo, reemplaza la línea anterior, p.ej. umbral <- 2000
      
      # Loop hora a hora dentro del dia 
      # hace un loop para cada día
      cat('Día', day, '\n')
      # dentro de este bucle debería tener cuidado con que alguno explote
      # de momento confiamos en que no explote ningun dia
      for (i in first.hour:(first.hour + 23)) {
        cat('Hora ', i - first.hour, ': ')
        
        if (length(active_sims) == 0) {
          cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
          break
        }
        
        # Preparamos aux con una fila repetida por cada simulación activa
        aux <- aux.data.h[i, , drop = FALSE]
        aux <- aux[rep(1, n.sim), , drop = FALSE]
        
        # Lag = valor simulado el día anterior para cada simulación activa
        lag_values <- sim_matrix[i - 1, 1:n.sim]
        
        # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
        bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
        # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
        if (length(bad_idx_rel) > 0) {
          # Mapear a índices absolutos en el conjunto original de simulaciones
          bad_idx_abs <- active_sims[bad_idx_rel]
          cat("  ⚠️ Detectadas simulaciones inválidas en hora ", i - first.hour - 1, ": ",
              paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral\n")
          # lag.values pasa a ser el umbral en aquellos que lo supera 
          #(aunque el df orignial aparezca el valor simulado)
          lag_values[bad_idx_abs] <- umbral
          # sim_matrix[i - 1, bad_idx_abs] <- NA_real_
          # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
          active_sims <- setdiff(active_sims, bad_idx_abs)
        }
        
        # Si tras eliminar no queda ninguna simulación válida, salimos
        if (length(active_sims) == 0) {
          cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
          break
        }
        
        # Reconstruir aux con la nueva longitud de simulaciones activas
        # en el caso de poner el umbral en el lag, simepre voy a tener el numero de simulaciones
        aux <- aux.data.h[i, , drop = FALSE]
        aux <- aux[rep(1, n.sim), , drop = FALSE]
        aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
        
        aux.prob <- aux
        aux <- aux[, names(Xq)]
        aux2 <- aux
        
        # obtención de mu y sigma solo para simulaciones activas
        cat('mu..')
        aux2$mu.new <- suppressWarnings(
          predict(mq2, newdata = aux, what = 'mu', type = 'response',
                  data = Xq)
        )
        cat('sigma..\n')
        aux2$sigma.new <- suppressWarnings(
          predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                  data = Xq)
        )
        
        # obtención de probabilidad de lluvia solo para simulaciones activas
        # cat('prob..\n')
        if (ocurrence == TRUE){
          aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
        }else{
          aux2$prob.dia <- df.gen$prob.p[i]
        }
        # Simular solo para las filas activas
        aux2 <- aux2 %>%
          select(prob.dia, mu.new, sigma.new) %>%
          rowwise() %>%
          mutate(
            shape.fv = 1 / sigma.new^2,
            rate.fv = 1 / sigma.new^2 / mu.new
          ) %>%
          ungroup() %>%
          rowwise() %>%
          mutate(
            sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
          ) %>%
          ungroup()
        
        # Guardar las simulaciones activas en la matriz en la posición correspondiente
        sim_matrix[i, ] <- as.numeric(aux2$sim)
        
        # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
      } # end for days
      
      
      # CORRECION 1
      # en todas simulaciones debe haber una hora que llueva
      # simulaciones en las que ese día no ha llovido
      sim.no.rain <- which(apply(sim_matrix[first.day, ], 2, function(x) sum(x > 0)) == 0)
      active_sims.repeat <- sim.no.rain
      #repetir simulaciones hasta que sim.no.rain = vacío
      while(length(sim.no.rain) != 0){
        for (i in (first.hour):(first.hour + 23)) {
          cat('Hora ', i - first.hour, ': ')
          
          if (length(active_sims.repeat) == 0) {
            cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
            break
          }
          
          # Preparamos aux con una fila repetida por cada simulación activa
          aux <- aux.data.h[i, , drop = FALSE]
          aux <- aux[rep(1, length(active_sims.repeat)), , drop = FALSE]
          
          # Lag = valor simulado el día anterior para cada simulación activa
          lag_values <- sim_matrix[i - 1, active_sims.repeat]
          
          # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
          bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
          # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
          if (length(bad_idx_rel) > 0) {
            # Mapear a índices absolutos en el conjunto original de simulaciones
            bad_idx_abs <- active_sims.repeat[bad_idx_rel]
            cat("  ⚠️ Detectadas simulaciones inválidas en hora ", 
                i - first.hour - 1, ": ",
                paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral \n")
            # Dejamos NA explícito (ya lo es) en la fila pasada; pero por claridad:
            #sim_matrix[i - 1, bad_idx_abs] <- NA_real_
            lag_values[bad_idx_rel] <- umbral
            # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
            # active_sims.repeat <- setdiff(active_sims.repeat, bad_idx_abs)
          }
          
          # Si tras eliminar no queda ninguna simulación válida, salimos
          if (length(active_sims.repeat) == 0) {
            cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
            break
          }
          
          # Reconstruir aux con la nueva longitud de simulaciones activas
          aux <- aux.data.h[i, , drop = FALSE]
          aux <- aux[rep(1, length(sim.no.rain)), , drop = FALSE]
          aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
          
          aux.prob <- aux
          aux <- aux[, names(Xq)]
          aux2 <- aux
          
          # obtención de mu y sigma solo para simulaciones activas
          cat('mu..')
          aux2$mu.new <- suppressWarnings(
            predict(mq2, newdata = aux, what = 'mu', type = 'response',
                    data = Xq)
          )
          cat('sigma..\n')
          aux2$sigma.new <- suppressWarnings(
            predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                    data = Xq)
          )
          
          # obtención de probabilidad de lluvia solo para simulaciones activas
          # cat('prob..\n')
          if (ocurrence == TRUE){
            aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
          }else{
            aux2$prob.dia <- df.gen$prob.p[i]
          }
          # Simular solo para las filas activas
          aux2 <- aux2 %>%
            select(prob.dia, mu.new, sigma.new) %>%
            rowwise() %>%
            mutate(
              shape.fv = 1 / sigma.new^2,
              rate.fv = 1 / sigma.new^2 / mu.new
            ) %>%
            ungroup() %>%
            rowwise() %>%
            mutate(
              sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
            ) %>%
            ungroup()
          
          # Guardar las simulaciones activas en la matriz en la posición correspondiente
          sim_matrix[i, sim.no.rain] <- as.numeric(aux2$sim)
          
          # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
        }
        
        sim.no.rain <- which(apply(sim_matrix[first.day, ], 2, function(x) sum(x > 0)) == 0)
        active_sims.repeat <- sim.no.rain
      }
      
      # CORRECION 2
      #normalizar con resepcto a lo que ha llovido
      rain.day <- aux.data.h[first.hour, paste0(station.p, '.day')]
      # factors <- rain.day / apply(sim_matrix[first.day, ], 2, sum)
      # 
      # sim_matrix[first.day, ] <- sweep(sim_matrix[first.day, ], 2, factors, `*`)
      
      sim_matrix[first.day, ] <- apply(sim_matrix[first.day, ], 2, 
                                       FUN = standardization, sum.to = rain.day,
                                       correction = TRUE)
      
    }
    
    # Convertir sim_matrix en data.frame con nombres adecuados y añadir a df.gen
    sim_df <- as.data.frame(sim_matrix)
    names(sim_df) <- paste0("y.sim.", 1:n.sim)
    
    cat(length(active_sims), ' han sobrevivido a la simulación')
  }else if(type =='day'){
    #data.h not needed
    #corrections not needed
    #calculation of harmonics
    l <- 1:365
    harm_l <- data.frame(l = l,
                         cs(l, harmonics = 1:4, 365))
    colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
    
    #construction of dataframe where we compute the simulations
    data.day <- data.day %>%
      left_join(harm_l, by = 'l') 
    
    station.p <- paste0(station, '.p')
    
    data.day <- data.day[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p'),
                         colnames(harm_l)[2:ncol(harm_l)])]
    
    data.day <- data.day %>% 
      mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]]))
    
    data.day <- data.day %>%
      mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0))%>% 
      relocate(Y, .after = !!station.p)
    
    data.day$date <- as.Date(paste(data.day$t, data.day$mes, data.day$dia.mes, sep = "-"), 
                           format = '%Y-%m-%d')
    
    # data.h <- data.h %>%
    #   mutate(
    #     gr.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] >= 2),
    #     less.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] < 2)
    #   )
    # 
    data.day <- data.day %>% filter(date >= period[1] & date <= period[2])
    
    aux.global_df <- global_df[global_df$STAID == station, ]
    
    data.day <- data.day %>%
      left_join(aux.global_df[, c(2, grep('z', names(aux.global_df)))], 
                by = 'date', keep = FALSE)
    
    Xq <- data.mq[[station]]
    
    Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"),
                       format = '%Y-%m-%d')
    Xq <- Xq %>% filter(date >= period[1] & date <= period[2])
    Xq <- Xq[, -ncol(Xq)]
    
    # quantity model
    mq <- models.list[[station]][[mq]]
    
    # formulas rewritten for new dataframe
    aux.mu.formula <- mq$mu.formula
    aux.mu.formula <- as.formula(
      paste(as.character(aux.mu.formula[2]), '~', 
            paste(labels(terms(aux.mu.formula)), collapse = '+'))
    )
    
    aux.sigma.formula <- mq$sigma.formula
    aux.sigma.formula <- as.formula(
      paste('~', as.character(aux.sigma.formula[2]))
    )
    
    
    mu.form <- sanitize_formula(aux.mu.formula)
    sigma.form <- sanitize_formula(aux.sigma.formula)
    
    #ocurrence model
    mo <- models.list[[station]][[mo]]
    
    #re ajuste modelo cantidad (debido a las nuevas formulas)
    mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                  family = GA, data = Xq, trace = FALSE)
    
    # filtrar por años si el valor es no nulo
    if (!is.null(years)){
      ind.year <- which(data.day$t %in% years)
    }else{
      ind.year <- nrow(data.day)
    }
    aux.data.day <- data.day[ind.year, ]
    
    # INICIO DE LA SIMULACION
    # MATRIZ DE SIMULACIONES
    n_days <- nrow(aux.data.day)
    sim_matrix <- matrix(NA_real_, nrow = n_days, ncol = n.sim)
    colnames(sim_matrix) <- paste0("p.day.sim.", 1:n.sim)
    
    # relleno de 0's para días sin lluvia observada
    ind.no.rain <- is.element(aux.data.day[[paste0(station.p)]], 0)
    sim_matrix[ind.no.rain, ] <- 0
    # dias que llueve
    n.days.rain <- sum(!ind.no.rain)
    
    # Umbral para considerar un valor "absurdo" — aquí uso el percentil 0.999 histórico
    umbral <- max(aux.data.day[[paste0(station, '.p')]], na.rm = T)
    
    # hourly simulation for days with rain
    for (i in 1:n.days.rain){
      # la matriz se va a ir acutalizando. 
      # primera hora del primer dia con lluvia observada
      aux <- apply(sim_matrix, 1, function(row) all(is.na(row)))
      first.day <- which(aux)[1]
      
      #luego cunado se vuelva a llamar se corresponderá con el segundo dia
      
      # Índices de simulaciones activas (las que aún no han sido descartadas)
      active_sims <- seq_len(n.sim)
      
      # Si prefieres un valor fijo, reemplaza la línea anterior, p.ej. umbral <- 2000
      
      # hace un loop para cada día
      cat('Día', i)
      
      
      if (length(active_sims) == 0) {
        cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
        break
      }
      
      # Preparamos aux con una fila repetida por cada simulación activa
      aux <- aux.data.day[first.day, , drop = FALSE]
      aux <- aux[rep(1, n.sim), , drop = FALSE]
      
      # Lag = valor simulado el día anterior para cada simulación activa
      lag_values <- sim_matrix[first.day - 1, 1:n.sim]
      
      # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
      bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
      # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
      if (length(bad_idx_rel) > 0) {
        # Mapear a índices absolutos en el conjunto original de simulaciones
        bad_idx_abs <- active_sims[bad_idx_rel]
        cat("  ⚠️ Detectadas simulaciones inválidas en hora ", i - first.hour - 1, ": ",
            paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral\n")
        # lag.values pasa a ser el umbral en aquellos que lo supera 
        #(aunque el df orignial aparezca el valor simulado)
        lag_values[bad_idx_abs] <- umbral
        # sim_matrix[i - 1, bad_idx_abs] <- NA_real_
        # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
        active_sims <- setdiff(active_sims, bad_idx_abs)
      }
      
      # Si tras eliminar no queda ninguna simulación válida, salimos
      if (length(active_sims) == 0) {
        cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
        break
      }
      
      # Reconstruir aux con la nueva longitud de simulaciones activas
      # en el caso de poner el umbral en el lag, simepre voy a tener el numero de simulaciones
      aux <- aux.data.day[first.day, , drop = FALSE]
      aux <- aux[rep(1, n.sim), , drop = FALSE]
      aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
      
      aux.prob <- aux
      aux <- aux[, names(aux) %in% names(Xq)]
      aux2 <- aux
      
      # obtención de mu y sigma solo para simulaciones activas
      cat('mu..')
      aux2$mu.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'mu', type = 'response',
                data = Xq)
      )
      cat('sigma..\n')
      aux2$sigma.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                data = Xq)
      )
      
      # obtención de probabilidad de lluvia solo para simulaciones activas
      # cat('prob..\n')
      if (ocurrence == TRUE){
        aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
      }else{
        aux2$prob.dia <- df.gen$prob.p[i]
      }
      # Simular solo para las filas activas
      aux2 <- aux2 %>%
        select(prob.dia, mu.new, sigma.new) %>%
        rowwise() %>%
        mutate(
          shape.fv = 1 / sigma.new^2,
          rate.fv = 1 / sigma.new^2 / mu.new
        ) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
        ) %>%
        ungroup()
      
      # Guardar las simulaciones activas en la matriz en la posición correspondiente
      sim_matrix[first.day, ] <- as.numeric(aux2$sim)
      
      # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
    } # end for days
    
    # Convertir sim_matrix en data.frame con nombres adecuados y añadir a df.gen
    sim_df <- as.data.frame(sim_matrix)
    names(sim_df) <- paste0("y.sim.", 1:n.sim)
    
    cat(length(active_sims), ' han sobrevivido a la simulación')
    
  }else{
    stop('Type invalid. Use hour or day.')
  }
  
  
  return(sim_df)
}
RAIN.GENERATOR.exp.1 <- function(station, data.h, data.day, period,
                              models.list, data.mq, mo, mq, 
                              ocurrence = TRUE, years = NULL, 
                              n.sim = 100, type){
  
  if (type == 'hour'){
    #calculation of harmonics
    l <- 1:365
    h <- 0:23
    harm_l <- data.frame(l = l,
                         cs(l, harmonics = 1:4, 365))
    colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
    harm_h <- data.frame(h = h,
                         cs(0:23, harmonics = 1:4, 24))
    colnames(harm_h)[2:ncol(harm_h)] <- paste0(colnames(harm_h)[2:ncol(harm_h)], '.h')
    
    #construction of dataframe where we compute the simulations
    data.h <- data.h %>%
      left_join(harm_l, by = 'l') %>%
      left_join(harm_h, by = 'h')
    
    station.p <- paste0(station, '.p')
    
    data.h <- data.h[, c('t', 'l', 'mes', 'dia.mes', 'h', paste0(station, '.p'),
                         colnames(harm_l)[2:ncol(harm_l)],
                         colnames(harm_h)[2:ncol(harm_h)])]
    
    data.day <- data.day[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p') )]
    
    colnames(data.day)[colnames(data.day) == station.p] <- paste0(station.p, '.day')
    
    data.h <- data.h %>%
      left_join(data.day, by = c('t', 'l', 'mes', 'dia.mes'))
    
    data.h <- data.h %>% 
      mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]]))
    
    data.h <- data.h %>%
      mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0))%>% 
      relocate(Y, .after = !!station.p)
    
    data.h$date <- as.Date(paste(data.h$t, data.h$mes, data.h$dia.mes, sep = "-"), 
                           format = '%Y-%m-%d')
    
    data.h <- data.h %>%
      mutate(
        gr.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] >= 2),
        less.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] < 2)
      )
    
    data.h <- data.h %>% filter(date >= period[1] & date <= period[2])
    
    Xq <- data.mq[[station]]
    
    Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"),
                       format = '%Y-%m-%d')
    Xq <- Xq %>% filter(date >= period[1] & date <= period[2])
    Xq <- Xq[, -ncol(Xq)]
    
    # quantity model
    mq <- models.list[[station]][[mq]]
    
    # formulas rewritten for new dataframe
    aux.mu.formula <- mq$mu.formula
    aux.mu.formula <- as.formula(
      paste(as.character(aux.mu.formula[2]), '~', 
            paste(labels(terms(aux.mu.formula)), collapse = '+'))
    )
    
    aux.sigma.formula <- mq$sigma.formula
    aux.sigma.formula <- as.formula(
      paste('~', as.character(aux.sigma.formula[2]))
    )
    
    
    mu.form <- sanitize_formula(aux.mu.formula)
    sigma.form <- sanitize_formula(aux.sigma.formula)
    
    #ocurrence model
    mo <- models.list[[station]][[mo]]
    
    #re ajuste modelo cantidad (debido a las nuevas formulas)
    mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                  family = GA, data = Xq, trace = FALSE)
    
    # filtrar por años si el valor es no nulo
    if (!is.null(years)){
      ind.year <- which(data.h$t %in% years)
    }else{
      ind.year <- nrow(data.h)
    }
    aux.data.h <- data.h[ind.year, ]
    
    #number of days with rain
    n.days.rain <- length(
      unique(
        aux.data.h[aux.data.h[[paste0(station.p, '.day')]] > 0, 'date']
      )
    )
    
    aux.data.h <- data.h[ind.year, ] # for correct computation of prediction
    
    
    # INICIO DE LA SIMULACION
    # MATRIZ DE SIMULACIONES
    n_days <- nrow(aux.data.h)
    #n.sim <- 100
    sim_matrix <- matrix(NA_real_, nrow = n_days, ncol = n.sim)
    colnames(sim_matrix) <- paste0("p.day.sim.", 1:n.sim)
    
    # relleno de 0's para días sin lluvia observada
    ind.no.rain <- is.element(aux.data.h[[paste0(station.p, '.day')]], 0)
    sim_matrix[ind.no.rain, ] <- 0
    
    # Umbral para considerar un valor "absurdo" — aquí uso el percentil 0.999 histórico
    umbral <- max(aux.data.h[[paste0(station, '.p')]], na.rm = T) * 2
    
    # hourly simulation for days with rain
    for (day in 1:n.days.rain){
      # la matriz se va a ir acutalizando. 
      # primera hora del primer dia con lluvia observada
      aux <- apply(sim_matrix, 1, function(row) all(is.na(row)))
      first.hour <- which(aux)[1]
      first.day <- first.hour:(first.hour + 23)
      
      #luego cunado se vuelva a llamar se corresponderá con el segundo dia
      
      # Índices de simulaciones activas (las que aún no han sido descartadas)
      active_sims <- seq_len(n.sim)
      
      # Si prefieres un valor fijo, reemplaza la línea anterior, p.ej. umbral <- 2000
      
      # Loop hora a hora dentro del dia 
      # hace un loop para cada día
      cat('Día', day, '\n')
      # dentro de este bucle debería tener cuidado con que alguno explote
      # de momento confiamos en que no explote ningun dia
      for (i in first.hour:(first.hour + 23)) {
        cat('Hora ', i - first.hour, ': ')
        
        if (length(active_sims) == 0) {
          cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
          break
        }
        
        # Preparamos aux con una fila repetida por cada simulación activa
        aux <- aux.data.h[i, , drop = FALSE]
        aux <- aux[rep(1, n.sim), , drop = FALSE]
        
        # Lag = valor simulado el día anterior para cada simulación activa
        lag_values <- sim_matrix[i - 1, 1:n.sim]
        
        # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
        bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
        # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
        if (length(bad_idx_rel) > 0) {
          # Mapear a índices absolutos en el conjunto original de simulaciones
          bad_idx_abs <- active_sims[bad_idx_rel]
          cat("  ⚠️ Detectadas simulaciones inválidas en hora ", i - first.hour - 1, ": ",
              paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral\n")
          # lag.values pasa a ser el umbral en aquellos que lo supera 
          #(aunque el df orignial aparezca el valor simulado)
          lag_values[bad_idx_abs] <- umbral
          # sim_matrix[i - 1, bad_idx_abs] <- NA_real_
          # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
          active_sims <- setdiff(active_sims, bad_idx_abs)
        }
        
        # Si tras eliminar no queda ninguna simulación válida, salimos
        if (length(active_sims) == 0) {
          cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
          break
        }
        
        # Reconstruir aux con la nueva longitud de simulaciones activas
        # en el caso de poner el umbral en el lag, simepre voy a tener el numero de simulaciones
        aux <- aux.data.h[i, , drop = FALSE]
        aux <- aux[rep(1, n.sim), , drop = FALSE]
        aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
        
        aux.prob <- aux
        aux <- aux[, names(Xq)]
        aux2 <- aux
        
        # obtención de mu y sigma solo para simulaciones activas
        cat('mu..')
        aux2$mu.new <- suppressWarnings(
          predict(mq2, newdata = aux, what = 'mu', type = 'response',
                  data = Xq)
        )
        cat('sigma..\n')
        aux2$sigma.new <- suppressWarnings(
          predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                  data = Xq)
        )
        
        # obtención de probabilidad de lluvia solo para simulaciones activas
        # cat('prob..\n')
        if (ocurrence == TRUE){
          aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
        }else{
          aux2$prob.dia <- df.gen$prob.p[i]
        }
        # Simular solo para las filas activas
        aux2 <- aux2 %>%
          select(prob.dia, mu.new, sigma.new) %>%
          rowwise() %>%
          mutate(
            shape.fv = 1 / sigma.new^2,
            rate.fv = 1 / sigma.new^2 / mu.new
          ) %>%
          ungroup() %>%
          rowwise() %>%
          mutate(
            sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
          ) %>%
          ungroup()
        
        # Guardar las simulaciones activas en la matriz en la posición correspondiente
        sim_matrix[i, ] <- as.numeric(aux2$sim)
        
        # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
      } # end for days
      
      
      # CORRECION 1
      # en todas simulaciones debe haber una hora que llueva
      # simulaciones en las que ese día no ha llovido
      sim.no.rain <- which(apply(sim_matrix[first.day, ], 2, function(x) sum(x > 0)) == 0)
      active_sims.repeat <- sim.no.rain
      #repetir simulaciones hasta que sim.no.rain = vacío
      while(length(sim.no.rain) != 0){
        for (i in (first.hour):(first.hour + 23)) {
          cat('Hora ', i - first.hour, ': ')
          
          if (length(active_sims.repeat) == 0) {
            cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
            break
          }
          
          # Preparamos aux con una fila repetida por cada simulación activa
          aux <- aux.data.h[i, , drop = FALSE]
          aux <- aux[rep(1, length(active_sims.repeat)), , drop = FALSE]
          
          # Lag = valor simulado el día anterior para cada simulación activa
          lag_values <- sim_matrix[i - 1, active_sims.repeat]
          
          # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
          bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
          # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
          if (length(bad_idx_rel) > 0) {
            # Mapear a índices absolutos en el conjunto original de simulaciones
            bad_idx_abs <- active_sims.repeat[bad_idx_rel]
            cat("  ⚠️ Detectadas simulaciones inválidas en hora ", 
                i - first.hour - 1, ": ",
                paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral \n")
            # Dejamos NA explícito (ya lo es) en la fila pasada; pero por claridad:
            #sim_matrix[i - 1, bad_idx_abs] <- NA_real_
            lag_values[bad_idx_rel] <- umbral
            # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
            active_sims.repeat <- setdiff(active_sims.repeat, bad_idx_abs)
          }
          
          # Si tras eliminar no queda ninguna simulación válida, salimos
          if (length(active_sims.repeat) == 0) {
            cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
            break
          }
          
          # Reconstruir aux con la nueva longitud de simulaciones activas
          aux <- aux.data.h[i, , drop = FALSE]
          aux <- aux[rep(1, length(sim.no.rain)), , drop = FALSE]
          aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
          
          aux.prob <- aux
          aux <- aux[, names(Xq)]
          aux2 <- aux
          
          # obtención de mu y sigma solo para simulaciones activas
          cat('mu..')
          aux2$mu.new <- suppressWarnings(
            predict(mq2, newdata = aux, what = 'mu', type = 'response',
                    data = Xq)
          )
          cat('sigma..\n')
          aux2$sigma.new <- suppressWarnings(
            predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                    data = Xq)
          )
          
          # obtención de probabilidad de lluvia solo para simulaciones activas
          # cat('prob..\n')
          if (ocurrence == TRUE){
            aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
          }else{
            aux2$prob.dia <- df.gen$prob.p[i]
          }
          # Simular solo para las filas activas
          aux2 <- aux2 %>%
            select(prob.dia, mu.new, sigma.new) %>%
            rowwise() %>%
            mutate(
              shape.fv = 1 / sigma.new^2,
              rate.fv = 1 / sigma.new^2 / mu.new
            ) %>%
            ungroup() %>%
            rowwise() %>%
            mutate(
              sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
            ) %>%
            ungroup()
          
          # Guardar las simulaciones activas en la matriz en la posición correspondiente
          sim_matrix[i, sim.no.rain] <- as.numeric(aux2$sim)
          
          # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
        }
        
        sim.no.rain <- which(apply(sim_matrix[first.day, ], 2, function(x) sum(x > 0)) == 0)
        active_sims.repeat <- sim.no.rain
      }
      
      # CORRECION 2
      #normalizar con resepcto a lo que ha llovido
      rain.day <- aux.data.h[first.hour, paste0(station.p, '.day')]
      factors <- rain.day / apply(sim_matrix[first.day, ], 2, sum)
      
      sim_matrix[first.day, ] <- sweep(sim_matrix[first.day, ], 2, factors, `*`)
      
    }
    
    # Convertir sim_matrix en data.frame con nombres adecuados y añadir a df.gen
    sim_df <- as.data.frame(sim_matrix)
    names(sim_df) <- paste0("y.sim.", 1:n.sim)
    
    cat(length(active_sims), ' han sobrevivido a la simulación')
  }else if(type =='day'){
    #data.h not needed
    #corrections not needed
    #calculation of harmonics
    l <- 1:365
    harm_l <- data.frame(l = l,
                         cs(l, harmonics = 1:4, 365))
    colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
    
    #construction of dataframe where we compute the simulations
    data.day <- data.day %>%
      left_join(harm_l, by = 'l') 
    
    station.p <- paste0(station, '.p')
    
    data.day <- data.day[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p'),
                             colnames(harm_l)[2:ncol(harm_l)])]
    
    data.day <- data.day %>% 
      mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]]))
    
    data.day <- data.day %>%
      mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0))%>% 
      relocate(Y, .after = !!station.p)
    
    data.day$date <- as.Date(paste(data.day$t, data.day$mes, data.day$dia.mes, sep = "-"), 
                             format = '%Y-%m-%d')
    
    # data.h <- data.h %>%
    #   mutate(
    #     gr.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] >= 2),
    #     less.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] < 2)
    #   )
    # 
    data.day <- data.day %>% filter(date >= period[1] & date <= period[2])
    
    aux.global_df <- global_df[global_df$STAID == station, ]
    
    data.day <- data.day %>%
      left_join(aux.global_df[, c(2, grep('z', names(aux.global_df)))], 
                by = 'date', keep = FALSE)
    
    Xq <- data.mq[[station]]
    
    Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"),
                       format = '%Y-%m-%d')
    Xq <- Xq %>% filter(date >= period[1] & date <= period[2])
    Xq <- Xq[, -ncol(Xq)]
    
    # quantity model
    mq <- models.list[[station]][[mq]]
    
    # formulas rewritten for new dataframe
    aux.mu.formula <- mq$mu.formula
    aux.mu.formula <- as.formula(
      paste(as.character(aux.mu.formula[2]), '~', 
            paste(labels(terms(aux.mu.formula)), collapse = '+'))
    )
    
    aux.sigma.formula <- mq$sigma.formula
    aux.sigma.formula <- as.formula(
      paste('~', as.character(aux.sigma.formula[2]))
    )
    
    
    mu.form <- sanitize_formula(aux.mu.formula)
    sigma.form <- sanitize_formula(aux.sigma.formula)
    
    #ocurrence model
    mo <- models.list[[station]][[mo]]
    
    #re ajuste modelo cantidad (debido a las nuevas formulas)
    mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                  family = GA, data = Xq, trace = FALSE)
    
    # filtrar por años si el valor es no nulo
    if (!is.null(years)){
      ind.year <- which(data.day$t %in% years)
    }else{
      ind.year <- nrow(data.day)
    }
    aux.data.day <- data.day[ind.year, ]
    
    # INICIO DE LA SIMULACION
    # MATRIZ DE SIMULACIONES
    n_days <- nrow(aux.data.day)
    sim_matrix <- matrix(NA_real_, nrow = n_days, ncol = n.sim)
    colnames(sim_matrix) <- paste0("p.day.sim.", 1:n.sim)
    
    # relleno de 0's para días sin lluvia observada
    ind.no.rain <- is.element(aux.data.day[[paste0(station.p)]], 0)
    sim_matrix[ind.no.rain, ] <- 0
    # dias que llueve
    n.days.rain <- sum(!ind.no.rain)
    
    # Umbral para considerar un valor "absurdo" — aquí uso el percentil 0.999 histórico
    umbral <- max(aux.data.day[[paste0(station, '.p')]], na.rm = T)
    
    # hourly simulation for days with rain
    for (i in 1:n.days.rain){
      # la matriz se va a ir acutalizando. 
      # primera hora del primer dia con lluvia observada
      aux <- apply(sim_matrix, 1, function(row) all(is.na(row)))
      first.day <- which(aux)[1]
      
      #luego cunado se vuelva a llamar se corresponderá con el segundo dia
      
      # Índices de simulaciones activas (las que aún no han sido descartadas)
      active_sims <- seq_len(n.sim)
      
      # Si prefieres un valor fijo, reemplaza la línea anterior, p.ej. umbral <- 2000
      
      # hace un loop para cada día
      cat('Día', i)
      
      
      if (length(active_sims) == 0) {
        cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
        break
      }
      
      # Preparamos aux con una fila repetida por cada simulación activa
      aux <- aux.data.day[first.day, , drop = FALSE]
      aux <- aux[rep(1, n.sim), , drop = FALSE]
      
      # Lag = valor simulado el día anterior para cada simulación activa
      lag_values <- sim_matrix[first.day - 1, 1:n.sim]
      
      # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
      bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
      # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
      if (length(bad_idx_rel) > 0) {
        # Mapear a índices absolutos en el conjunto original de simulaciones
        bad_idx_abs <- active_sims[bad_idx_rel]
        cat("  ⚠️ Detectadas simulaciones inválidas en hora ", i - first.hour - 1, ": ",
            paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral\n")
        # lag.values pasa a ser el umbral en aquellos que lo supera 
        #(aunque el df orignial aparezca el valor simulado)
        lag_values[bad_idx_abs] <- umbral
        # sim_matrix[i - 1, bad_idx_abs] <- NA_real_
        # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
        active_sims <- setdiff(active_sims, bad_idx_abs)
      }
      
      # Si tras eliminar no queda ninguna simulación válida, salimos
      if (length(active_sims) == 0) {
        cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
        break
      }
      
      # Reconstruir aux con la nueva longitud de simulaciones activas
      # en el caso de poner el umbral en el lag, simepre voy a tener el numero de simulaciones
      aux <- aux.data.day[first.day, , drop = FALSE]
      aux <- aux[rep(1, n.sim), , drop = FALSE]
      aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
      
      aux.prob <- aux
      aux <- aux[, names(aux) %in% names(Xq)]
      aux2 <- aux
      
      # obtención de mu y sigma solo para simulaciones activas
      cat('mu..')
      aux2$mu.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'mu', type = 'response',
                data = Xq)
      )
      cat('sigma..\n')
      aux2$sigma.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                data = Xq)
      )
      
      # obtención de probabilidad de lluvia solo para simulaciones activas
      # cat('prob..\n')
      if (ocurrence == TRUE){
        aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
      }else{
        aux2$prob.dia <- df.gen$prob.p[i]
      }
      # Simular solo para las filas activas
      aux2 <- aux2 %>%
        select(prob.dia, mu.new, sigma.new) %>%
        rowwise() %>%
        mutate(
          shape.fv = 1 / sigma.new^2,
          rate.fv = 1 / sigma.new^2 / mu.new
        ) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
        ) %>%
        ungroup()
      
      # Guardar las simulaciones activas en la matriz en la posición correspondiente
      sim_matrix[first.day, ] <- as.numeric(aux2$sim)
      
      # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
    } # end for days
    
    # Convertir sim_matrix en data.frame con nombres adecuados y añadir a df.gen
    sim_df <- as.data.frame(sim_matrix)
    names(sim_df) <- paste0("y.sim.", 1:n.sim)
    
    cat(length(active_sims), ' han sobrevivido a la simulación')
    
  }else{
    stop('Type invalid. Use hour or day.')
  }
  
  
  return(sim_df)
}
RAIN.GENERATOR.exp.2 <- function(station, data.h, data.day, period,
                              models.list.mo, models.list.mq, data.mq, mo, mq, 
                              ocurrence = TRUE, years = NULL, 
                              n.sim = 100, type){
  
  if (type == 'hour'){
    #calculation of harmonics
    l <- 1:365
    h <- 0:23
    harm_l <- data.frame(l = l,
                         cs(l, harmonics = 1:4, 365))
    colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
    harm_h <- data.frame(h = h,
                         cs(0:23, harmonics = 1:4, 24))
    colnames(harm_h)[2:ncol(harm_h)] <- paste0(colnames(harm_h)[2:ncol(harm_h)], '.h')
    
    #construction of dataframe where we compute the simulations
    data.h <- data.h %>%
      left_join(harm_l, by = 'l') %>%
      left_join(harm_h, by = 'h')
    
    station.p <- paste0(station, '.p')
    
    data.h <- data.h[, c('t', 'l', 'mes', 'dia.mes', 'h', paste0(station, '.p'),
                         colnames(harm_l)[2:ncol(harm_l)],
                         colnames(harm_h)[2:ncol(harm_h)])]
    
    data.day <- data.day[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p') )]
    
    colnames(data.day)[colnames(data.day) == station.p] <- paste0(station.p, '.day')
    
    data.h <- data.h %>%
      left_join(data.day, by = c('t', 'l', 'mes', 'dia.mes'))
    
    data.h <- data.h %>% 
      mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]]))
    
    data.h <- data.h %>%
      mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0))%>% 
      relocate(Y, .after = !!station.p)
    
    data.h$date <- as.Date(paste(data.h$t, data.h$mes, data.h$dia.mes, sep = "-"), 
                           format = '%Y-%m-%d')
    
    data.h <- data.h %>%
      mutate(
        gr.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] >= 2),
        less.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] < 2)
      )
    
    data.h <- data.h %>% filter(date >= period[1] & date <= period[2])
    
    Xq <- data.mq[[station]]
    
    Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"),
                       format = '%Y-%m-%d')
    Xq <- Xq %>% filter(date >= period[1] & date <= period[2])
    Xq <- Xq[, -ncol(Xq)]
    
    # quantity model
    mq <- models.list.mq[[station]][[mq]]
    
    # formulas rewritten for new dataframe
    aux.mu.formula <- mq$mu.formula
    aux.mu.formula <- as.formula(
      paste(as.character(aux.mu.formula[2]), '~', 
            paste(labels(terms(aux.mu.formula)), collapse = '+'))
    )
    
    aux.sigma.formula <- mq$sigma.formula
    aux.sigma.formula <- as.formula(
      paste('~', as.character(aux.sigma.formula[2]))
    )
    
    
    mu.form <- sanitize_formula(aux.mu.formula)
    sigma.form <- sanitize_formula(aux.sigma.formula)
    
    #ocurrence model
    mo <- models.list.mo[[station]][[mo]]
    
    #re ajuste modelo cantidad (debido a las nuevas formulas)
    mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                  family = GA, data = Xq, trace = FALSE)
    
    # filtrar por años si el valor es no nulo
    if (!is.null(years)){
      ind.year <- which(data.h$t %in% years)
    }else{
      ind.year <- nrow(data.h)
    }
    aux.data.h <- data.h[ind.year, ]
    
    #number of days with rain
    n.days.rain <- length(
      unique(
        aux.data.h[aux.data.h[[paste0(station.p, '.day')]] > 0, 'date']
      )
    )
    
    aux.data.h <- data.h[ind.year, ] # for correct computation of prediction
    
    
    # INICIO DE LA SIMULACION
    # MATRIZ DE SIMULACIONES
    n_days <- nrow(aux.data.h)
    #n.sim <- 100
    sim_matrix <- matrix(NA_real_, nrow = n_days, ncol = n.sim)
    colnames(sim_matrix) <- paste0("p.day.sim.", 1:n.sim)
    
    # relleno de 0's para días sin lluvia observada
    ind.no.rain <- is.element(aux.data.h[[paste0(station.p, '.day')]], 0)
    sim_matrix[ind.no.rain, ] <- 0
    
    # Umbral para considerar un valor "absurdo" — aquí uso el percentil 0.999 histórico
    umbral <- max(aux.data.h[[paste0(station, '.p')]], na.rm = T)
    
    # hourly simulation for days with rain
    for (day in 1:n.days.rain){
      # la matriz se va a ir acutalizando. 
      # primera hora del primer dia con lluvia observada
      aux <- apply(sim_matrix, 1, function(row) all(is.na(row)))
      first.hour <- which(aux)[1]
      first.day <- first.hour:(first.hour + 23)
      
      #luego cunado se vuelva a llamar se corresponderá con el segundo dia
      
      # Índices de simulaciones activas (las que aún no han sido descartadas)
      active_sims <- seq_len(n.sim)
      
      # Si prefieres un valor fijo, reemplaza la línea anterior, p.ej. umbral <- 2000
      
      # Loop hora a hora dentro del dia 
      # hace un loop para cada día
      cat('Día', day, '\n')
      # dentro de este bucle debería tener cuidado con que alguno explote
      # de momento confiamos en que no explote ningun dia
      for (i in first.hour:(first.hour + 23)) {
        cat('Hora ', i - first.hour, ': ')
        
        if (length(active_sims) == 0) {
          cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
          break
        }
        
        # Preparamos aux con una fila repetida por cada simulación activa
        aux <- aux.data.h[i, , drop = FALSE]
        aux <- aux[rep(1, n.sim), , drop = FALSE]
        
        # Lag = valor simulado el día anterior para cada simulación activa
        lag_values <- sim_matrix[i - 1, 1:n.sim]
        
        # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
        bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
        # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
        if (length(bad_idx_rel) > 0) {
          # Mapear a índices absolutos en el conjunto original de simulaciones
          bad_idx_abs <- active_sims[bad_idx_rel]
          cat("  ⚠️ Detectadas simulaciones inválidas en hora ", i - first.hour - 1, ": ",
              paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral\n")
          # lag.values pasa a ser el umbral en aquellos que lo supera 
          #(aunque el df orignial aparezca el valor simulado)
          lag_values[bad_idx_abs] <- umbral
          # sim_matrix[i - 1, bad_idx_abs] <- NA_real_
          # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
          active_sims <- setdiff(active_sims, bad_idx_abs)
        }
        
        # Si tras eliminar no queda ninguna simulación válida, salimos
        if (length(active_sims) == 0) {
          cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
          break
        }
        
        # Reconstruir aux con la nueva longitud de simulaciones activas
        # en el caso de poner el umbral en el lag, simepre voy a tener el numero de simulaciones
        aux <- aux.data.h[i, , drop = FALSE]
        aux <- aux[rep(1, n.sim), , drop = FALSE]
        aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
        
        aux.prob <- aux
        aux <- aux[, names(Xq)]
        aux2 <- aux
        
        # obtención de mu y sigma solo para simulaciones activas
        cat('mu..')
        aux2$mu.new <- suppressWarnings(
          predict(mq2, newdata = aux, what = 'mu', type = 'response',
                  data = Xq)
        )
        cat('sigma..\n')
        aux2$sigma.new <- suppressWarnings(
          predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                  data = Xq)
        )
        
        # obtención de probabilidad de lluvia solo para simulaciones activas
        # cat('prob..\n')
        if (ocurrence == TRUE){
          aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
        }else{
          aux2$prob.dia <- df.gen$prob.p[i]
        }
        # Simular solo para las filas activas
        aux2 <- aux2 %>%
          select(prob.dia, mu.new, sigma.new) %>%
          rowwise() %>%
          mutate(
            shape.fv = 1 / sigma.new^2,
            rate.fv = 1 / sigma.new^2 / mu.new
          ) %>%
          ungroup() %>%
          rowwise() %>%
          mutate(
            sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
          ) %>%
          ungroup()
        
        # Guardar las simulaciones activas en la matriz en la posición correspondiente
        sim_matrix[i, ] <- as.numeric(aux2$sim)
        
        # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
      } # end for days
      
      
      # CORRECION 1
      # en todas simulaciones debe haber una hora que llueva
      # simulaciones en las que ese día no ha llovido
      sim.no.rain <- which(apply(sim_matrix[first.day, ], 2, function(x) sum(x > 0)) == 0)
      active_sims.repeat <- sim.no.rain
      #repetir simulaciones hasta que sim.no.rain = vacío
      while(length(sim.no.rain) != 0){
        for (i in (first.hour):(first.hour + 23)) {
          cat('Hora ', i - first.hour, ': ')
          
          if (length(active_sims.repeat) == 0) {
            cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
            break
          }
          
          # Preparamos aux con una fila repetida por cada simulación activa
          aux <- aux.data.h[i, , drop = FALSE]
          aux <- aux[rep(1, length(active_sims.repeat)), , drop = FALSE]
          
          # Lag = valor simulado el día anterior para cada simulación activa
          lag_values <- sim_matrix[i - 1, active_sims.repeat]
          
          # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
          bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
          # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
          if (length(bad_idx_rel) > 0) {
            # Mapear a índices absolutos en el conjunto original de simulaciones
            bad_idx_abs <- active_sims.repeat[bad_idx_rel]
            cat("  ⚠️ Detectadas simulaciones inválidas en hora ", 
                i - first.hour - 1, ": ",
                paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral \n")
            # Dejamos NA explícito (ya lo es) en la fila pasada; pero por claridad:
            #sim_matrix[i - 1, bad_idx_abs] <- NA_real_
            lag_values[bad_idx_rel] <- umbral
            # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
            active_sims.repeat <- setdiff(active_sims.repeat, bad_idx_abs)
          }
          
          # Si tras eliminar no queda ninguna simulación válida, salimos
          if (length(active_sims.repeat) == 0) {
            cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
            break
          }
          
          # Reconstruir aux con la nueva longitud de simulaciones activas
          aux <- aux.data.h[i, , drop = FALSE]
          aux <- aux[rep(1, length(sim.no.rain)), , drop = FALSE]
          aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
          
          aux.prob <- aux
          aux <- aux[, names(Xq)]
          aux2 <- aux
          
          # obtención de mu y sigma solo para simulaciones activas
          cat('mu..')
          aux2$mu.new <- suppressWarnings(
            predict(mq2, newdata = aux, what = 'mu', type = 'response',
                    data = Xq)
          )
          cat('sigma..\n')
          aux2$sigma.new <- suppressWarnings(
            predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                    data = Xq)
          )
          
          # obtención de probabilidad de lluvia solo para simulaciones activas
          # cat('prob..\n')
          if (ocurrence == TRUE){
            aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
          }else{
            aux2$prob.dia <- df.gen$prob.p[i]
          }
          # Simular solo para las filas activas
          aux2 <- aux2 %>%
            select(prob.dia, mu.new, sigma.new) %>%
            rowwise() %>%
            mutate(
              shape.fv = 1 / sigma.new^2,
              rate.fv = 1 / sigma.new^2 / mu.new
            ) %>%
            ungroup() %>%
            rowwise() %>%
            mutate(
              sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
            ) %>%
            ungroup()
          
          # Guardar las simulaciones activas en la matriz en la posición correspondiente
          sim_matrix[i, sim.no.rain] <- as.numeric(aux2$sim)
          
          # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
        }
        
        sim.no.rain <- which(apply(sim_matrix[first.day, ], 2, function(x) sum(x > 0)) == 0)
        active_sims.repeat <- sim.no.rain
      }
      
      # CORRECION 2
      #normalizar con resepcto a lo que ha llovido
      rain.day <- aux.data.h[first.hour, paste0(station.p, '.day')]
      factors <- rain.day / apply(sim_matrix[first.day, ], 2, sum)
      
      sim_matrix[first.day, ] <- sweep(sim_matrix[first.day, ], 2, factors, `*`)
      
    }
    
    # Convertir sim_matrix en data.frame con nombres adecuados y añadir a df.gen
    sim_df <- as.data.frame(sim_matrix)
    names(sim_df) <- paste0("y.sim.", 1:n.sim)
    
    cat(length(active_sims), ' han sobrevivido a la simulación')
  }else if(type =='day'){
    #data.h not needed
    #corrections not needed
    #calculation of harmonics
    l <- 1:365
    harm_l <- data.frame(l = l,
                         cs(l, harmonics = 1:4, 365))
    colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
    
    #construction of dataframe where we compute the simulations
    data.day <- data.day %>%
      left_join(harm_l, by = 'l') 
    
    station.p <- paste0(station, '.p')
    
    data.day <- data.day[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p'),
                             colnames(harm_l)[2:ncol(harm_l)])]
    
    data.day <- data.day %>% 
      mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]]))
    
    data.day <- data.day %>%
      mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0))%>% 
      relocate(Y, .after = !!station.p)
    
    data.day$date <- as.Date(paste(data.day$t, data.day$mes, data.day$dia.mes, sep = "-"), 
                             format = '%Y-%m-%d')
    
    # data.h <- data.h %>%
    #   mutate(
    #     gr.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] >= 2),
    #     less.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] < 2)
    #   )
    # 
    data.day <- data.day %>% filter(date >= period[1] & date <= period[2])
    
    aux.global_df <- global_df[global_df$STAID == station, ]
    
    data.day <- data.day %>%
      left_join(aux.global_df[, c(2, grep('z', names(aux.global_df)))], 
                by = 'date', keep = FALSE)
    
    Xq <- data.mq[[station]]
    
    Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"),
                       format = '%Y-%m-%d')
    Xq <- Xq %>% filter(date >= period[1] & date <= period[2])
    Xq <- Xq[, -ncol(Xq)]
    
    # quantity model
    mq <- models.list[[station]][[mq]]
    
    # formulas rewritten for new dataframe
    aux.mu.formula <- mq$mu.formula
    aux.mu.formula <- as.formula(
      paste(as.character(aux.mu.formula[2]), '~', 
            paste(labels(terms(aux.mu.formula)), collapse = '+'))
    )
    
    aux.sigma.formula <- mq$sigma.formula
    aux.sigma.formula <- as.formula(
      paste('~', as.character(aux.sigma.formula[2]))
    )
    
    
    mu.form <- sanitize_formula(aux.mu.formula)
    sigma.form <- sanitize_formula(aux.sigma.formula)
    
    #ocurrence model
    mo <- models.list[[station]][[mo]]
    
    #re ajuste modelo cantidad (debido a las nuevas formulas)
    mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                  family = GA, data = Xq, trace = FALSE)
    
    # filtrar por años si el valor es no nulo
    if (!is.null(years)){
      ind.year <- which(data.day$t %in% years)
    }else{
      ind.year <- nrow(data.day)
    }
    aux.data.day <- data.day[ind.year, ]
    
    # INICIO DE LA SIMULACION
    # MATRIZ DE SIMULACIONES
    n_days <- nrow(aux.data.day)
    sim_matrix <- matrix(NA_real_, nrow = n_days, ncol = n.sim)
    colnames(sim_matrix) <- paste0("p.day.sim.", 1:n.sim)
    
    # relleno de 0's para días sin lluvia observada
    ind.no.rain <- is.element(aux.data.day[[paste0(station.p)]], 0)
    sim_matrix[ind.no.rain, ] <- 0
    # dias que llueve
    n.days.rain <- sum(!ind.no.rain)
    
    # Umbral para considerar un valor "absurdo" — aquí uso el percentil 0.999 histórico
    umbral <- max(aux.data.day[[paste0(station, '.p')]], na.rm = T)
    
    # hourly simulation for days with rain
    for (i in 1:n.days.rain){
      # la matriz se va a ir acutalizando. 
      # primera hora del primer dia con lluvia observada
      aux <- apply(sim_matrix, 1, function(row) all(is.na(row)))
      first.day <- which(aux)[1]
      
      #luego cunado se vuelva a llamar se corresponderá con el segundo dia
      
      # Índices de simulaciones activas (las que aún no han sido descartadas)
      active_sims <- seq_len(n.sim)
      
      # Si prefieres un valor fijo, reemplaza la línea anterior, p.ej. umbral <- 2000
      
      # hace un loop para cada día
      cat('Día', i)
      
      
      if (length(active_sims) == 0) {
        cat("No quedan simulaciones activas en la hora ", i, ". Terminando.\n")
        break
      }
      
      # Preparamos aux con una fila repetida por cada simulación activa
      aux <- aux.data.day[first.day, , drop = FALSE]
      aux <- aux[rep(1, n.sim), , drop = FALSE]
      
      # Lag = valor simulado el día anterior para cada simulación activa
      lag_values <- sim_matrix[first.day - 1, 1:n.sim]
      
      # Detectar cuáles son inválidos ahora (NA, Inf, o > umbral)
      bad_idx_rel <- which(!is.finite(lag_values) | is.na(lag_values) | (lag_values > umbral))
      # CORREGIR. SI PASA UMBRAL--> LAG PARA LA SIGUIENTE SERÁ EL MAXIMO OBSERVADO
      if (length(bad_idx_rel) > 0) {
        # Mapear a índices absolutos en el conjunto original de simulaciones
        bad_idx_abs <- active_sims[bad_idx_rel]
        cat("  ⚠️ Detectadas simulaciones inválidas en hora ", i - first.hour - 1, ": ",
            paste0(bad_idx_abs, collapse = ", "), " -> lag = umbral\n")
        # lag.values pasa a ser el umbral en aquellos que lo supera 
        #(aunque el df orignial aparezca el valor simulado)
        lag_values[bad_idx_abs] <- umbral
        # sim_matrix[i - 1, bad_idx_abs] <- NA_real_
        # Eliminamos esas simulaciones de active_sims para que no se usen como lag en adelante
        active_sims <- setdiff(active_sims, bad_idx_abs)
      }
      
      # Si tras eliminar no queda ninguna simulación válida, salimos
      if (length(active_sims) == 0) {
        cat("No quedan simulaciones activas tras depuración en día ", i, ". Terminando.\n")
        break
      }
      
      # Reconstruir aux con la nueva longitud de simulaciones activas
      # en el caso de poner el umbral en el lag, simepre voy a tener el numero de simulaciones
      aux <- aux.data.day[first.day, , drop = FALSE]
      aux <- aux[rep(1, n.sim), , drop = FALSE]
      aux[[paste0(station, '.p.lag')]] <- as.numeric(lag_values)
      
      aux.prob <- aux
      aux <- aux[, names(aux) %in% names(Xq)]
      aux2 <- aux
      
      # obtención de mu y sigma solo para simulaciones activas
      cat('mu..')
      aux2$mu.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'mu', type = 'response',
                data = Xq)
      )
      cat('sigma..\n')
      aux2$sigma.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'sigma', type = 'response',
                data = Xq)
      )
      
      # obtención de probabilidad de lluvia solo para simulaciones activas
      # cat('prob..\n')
      if (ocurrence == TRUE){
        aux2$prob.dia <- predict(mo, newdata = aux.prob, type = 'response')
      }else{
        aux2$prob.dia <- df.gen$prob.p[i]
      }
      # Simular solo para las filas activas
      aux2 <- aux2 %>%
        select(prob.dia, mu.new, sigma.new) %>%
        rowwise() %>%
        mutate(
          shape.fv = 1 / sigma.new^2,
          rate.fv = 1 / sigma.new^2 / mu.new
        ) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
        ) %>%
        ungroup()
      
      # Guardar las simulaciones activas en la matriz en la posición correspondiente
      sim_matrix[first.day, ] <- as.numeric(aux2$sim)
      
      # Nota: las columnas correspondientes a simulaciones descartadas permanecen en NA
    } # end for days
    
    # Convertir sim_matrix en data.frame con nombres adecuados y añadir a df.gen
    sim_df <- as.data.frame(sim_matrix)
    names(sim_df) <- paste0("y.sim.", 1:n.sim)
    
    cat(length(active_sims), ' han sobrevivido a la simulación')
    
  }else{
    stop('Type invalid. Use hour or day.')
  }
  
  
  return(sim_df)
}
#---ANÁLISIS OCURRENCIA GENERADORES----
set.seed(05052002)
y.sim.og.cor <- RAIN.GENERATOR.og('A126', df_hours, df_days, per.comun.h, 
                          common.models.final, X.MHQ, 'MHO', 'MHQ',
                          ocurrence = TRUE, years = 2015,
                          type = 'hour')

set.seed(05052002)
y.sim.exp.1 <- RAIN.GENERATOR.exp.1('A126', df_hours, df_days, per.comun.h, 
                              common.models.final, X.MHQ, 'MHO', 'MHQ',
                              ocurrence = TRUE, years = 2015,
                              type = 'hour')
set.seed(05052002)
y.sim.exp.2 <- RAIN.GENERATOR.exp.2('A126', df_hours, df_days, per.comun.h, 
                                    models.list.mo = common.models.final, 
                                    models.list.mq = MHQ, X.MHQ, 'MHO', 'M8',
                                    ocurrence = TRUE, years = 2015,
                                    type = 'hour')

# analysis of effects for another MHQ
m <- common.models.final[['A126']][['MHQ']]
summary(m)
library(gam)
X <- X.MHQ[['A126']]

X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"),
                   format = '%Y-%m-%d')
X <- X %>% filter(date >= per.comun.h[1] & date <= per.comun.h[2])

plot(gam(formula = A126.p ~ s(A126.p.lag), data = X, family = Gamma(link = 'log')))
abline(v = 8)

#new model
common.models.final[['A126']][['MHQ.2']] <- update(m, formula = .~. - poly(A126.p.lag, 2) +
                   I(A126.p.lag >= 8):I(A126.p.lag) +
                   I(A126.p.lag < 8):I(log(pmax(A126.p.lag, 1e-06))), 
                 what = 'mu')

set.seed(05052002)
y.sim.exp.3.cor <- RAIN.GENERATOR.og('A126', df_hours, df_days, per.comun.h, 
                              common.models.final, X.MHQ, 'MHO', 'MHQ.2',
                              ocurrence = TRUE, years = 2015,
                              type = 'hour')



hourly.simulations <- list(
  y.sim.og = y.sim.og,
  y.sim.og.cor = y.sim.og.cor,
  y.sim.exp.1 = y.sim.exp.1,
  y.sim.exp.2 = y.sim.exp.2,
  y.sim.exp.3 = y.sim.exp.3,
  y.sim.exp.3.cor = y.sim.exp.3.cor
)
qsave(hourly.simulations, 'hourly.simulations.qs')
#--
set.seed(05052002)
basura <- RAIN.GENERATOR('A126', df_hours, df_days, per.comun.day, 
                           common.models.final, X.MDQ, 'MDO', 'MDQ',
                           ocurrence = TRUE, years = 2015,
                           type = 'day')

ocurrence.analysis <- function(type, station, data.h, data.day, 
                               period, years, y.sim){
  if (type == 'hour'){
    station.p <- paste0(station, '.p')
    X <- data.h[, c('t', 'l', 'mes', 'dia.mes', 'h', paste0(station, '.p'))]
    
    aux.X <- data.day[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p') )]
    
    colnames(aux.X)[colnames(aux.X) == station.p] <- paste0(station.p, '.day')
    
    X <- X %>%
      left_join(aux.X, by = c('t', 'l', 'mes', 'dia.mes'))
    
    X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                      format = '%Y-%m-%d')
    
    X <- X %>% filter(date >= period[1] & date <= period[2] &
                        t %in% years)
    ind.rain <- !is.element(X[[paste0(station.p, '.day')]], 0)
    
    n <- nrow(y.sim[ind.rain, ])
    boxplot(apply(y.sim[ind.rain, ], 2, function(x) sum(x == 0)/n))
    points(1, sum(X[ind.rain, station.p] == 0)/n, pch = 19, col = 'red')
    
    y.sim.rain <- y.sim[ind.rain, ]
    y.sim.rain <- cbind(X$date[ind.rain], y.sim.rain)
    print(
      head(cbind(y.sim.rain[y.sim.rain$y.sim.1 > 0, c(1,2)], X[ind.rain, paste0(station.p, '.day')][y.sim.rain$y.sim.1 > 0]))
    )
    #llueve todos días simulados? En teoría si
    print(sum(X[ind.rain, paste0(station.p, '.day')] >0))
    
    colnames(y.sim.rain)[1] <- 'date'
    n.days <- length(unique(y.sim.rain$date))
    
    y.sim.rain <- y.sim.rain %>%
      group_by(date) %>%
      summarise(across(1:(ncol(y.sim) - 1), ~ any(.x > 0))) %>%
      ungroup()
    
    sim.rain <- apply(y.sim.rain[, 2:ncol(y.sim.rain)], 2, sum)
    plot(sim.rain, ylim =c(50, n.days))
    abline(h = n.days, col = 'red')
    
  }else if (type == 'day'){
    station.p <- paste0(station, '.p')
    X <- data.h[, c('t', 'l', 'mes', 'dia.mes', paste0(station, '.p'))]
    
    X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                      format = '%Y-%m-%d')
    
    X <- X %>% filter(date >= period[1] & date <= period[2] &
                        t %in% years)
    
    n <- nrow(y.sim)
    boxplot(apply(y.sim, 2, function(x) sum(x == 0)/n))
    points(1, sum(X[, station.p] == 0)/n, pch = 19, col = 'red')
    
  }else{
    stop('Type invalid. Use hour or day.')
  }
}

ocurrence.analysis(type = 'hour', station = 'A126', data.h = df_hours,
                   data.day = df_days, period = per.comun.h, years = 2015,
                   y.sim = y.sim.exp.3.cor)

ocurrence.analysis(type = 'day', station = 'A126', data.h = df_hours,
                   data.day = df_days, period = per.comun.day, years = 2015,
                   y.sim = basura)


#extr
m <- common.models.final[['A126']][['MHO']]
coef <- m$coefficients [c(12,13)]
plot(gam(formula = m$formula, data = m$data))
data.aux <- m$data
aux.x <- 0:50/25
plot(aux.x, coef[1]*log(aux.x) + coef[2] * aux.x, type = 'l')


m <- common.models.final[['A126']][['MHQ']]
summary(m)

summary(m$sigma.fv)

aux.x <- 0:23
plot(aux.x, -0.11564 * sin(2 * pi * aux.x / 24) + 0.06046 * cos(2 * pi * aux.x / 24), type = 'l')
aux.d <- 1:365
plot(aux.d, -0.04777 * sin(2 * pi * aux.d / 365) - 0.21749 * cos(2 * pi * aux.d / 365), type = 'l')
aux.x <- 1:50
plot(aux.x, 0.46089* log(aux.x), type = 'l')


plot(gam(formula = m$mu.formula, data = X))

#----ANÁLISIS GENERADOR HORARIO----
station <- 'A126'

boxplot.q.sim <- function(station, data.mo, data.h, period, 
                          years = NULL, months = NULL, seasons = NULL,
                          y.sim, 
                          quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
                          plot = TRUE){
  #datos donde se ajusta el model, i.e, dias de lluvia
  X <- data.mo[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = '-'), 
                    format = '%Y-%m-%d')
  X <- X %>% 
    filter(date >= period[1] & date <= period[2])
  
  data.h$date <- as.Date(paste(data.h$t, data.h$mes, data.h$dia.mes, sep = "-"), 
                         format = '%Y-%m-%d')
  data.h <- data.h %>% filter(date >= period[1] & date <= period[2])
  
  if (!is.null(years)){
    X <- X %>%
      filter(t %in% years)
    data.h <- data.h %>%
      filter(t %in% years)
  }
  
  x.obs <- X[[paste0(station, '.p')]]
  
  #obtencion del subset de simulaciones con lluvias
  aux.ind <- which(data.h$date %in% unique(X$date))
  y.sim.rain <- y.sim[aux.ind, ]
  # ahora X y y.sim.rain tiene misma longitud
  
  if (!is.null(seasons)){
    ind.list <- list(
      ind = 1:dim(X)[1],
      ind.jja = which(X$mes %in% c(6, 7, 8)),
      ind.son = which(X$mes %in% c(9, 10, 11)),
      ind.djf = which(X$mes %in% c(12, 1, 2)),
      ind.mam = which(X$mes %in% c(3, 4, 5))
    )
    
  }else if (!is.null(month)){
    ind.list <- list(ind = which(X$mes %in% month))
  }else{
    ind.list <- list(ind = 1:dim(X)[1])
  }
  
  names.q <- paste0('q.', quantiles)
  
  season.names <- c('ALL', 'JJA', 'SON', 'DJF', 'MAM')
  
  i <- 1
  for (indx in ind.list){
    x.obs.indx <- x.obs[indx]
    q.obs <- quantile(x.obs.indx[x.obs.indx > 0], probs = quantiles)
    names(q.obs) <- names.q
    #quantiles of simulation
    y.sim.indx <- y.sim.rain[indx, ]
    #q.sim <- t(apply(y.sim[indx, ], 2, quantile, probs = quantiles, na.rm = T))
    q.sim <- t(apply(y.sim.indx, 2, function(col) {
      pos <- col[col > 0]
      if (length(pos) > 0) {
        quantile(pos, probs = quantiles, na.rm = TRUE)
      } else {
        rep(NA, length(quantiles))  # si no hay positivos, devolvemos NAs
      }
    }))
    if (plot == TRUE){
      if (!is.null(seasons)){
        title <- paste0("Simulated quantiles vs obs. quantiles ", 
                        station, '-', 'MHQ', 
                        ' ', season.names[i])
      }
      else if (!is.null(month)){
        title <- paste0("Simulated quantiles vs obs. quantiles ", 
                        station, '-', 'MHQ', 
                        ' (months ', paste0(month, collapse = '-'), ')')
      }
      else{
        title <- paste0("Simulated quantiles vs obs. quantiles ", 
                        station, '-', 'MHQ')
      }
      boxplot(q.sim,
              at = q.obs,               
              names = paste0(quantiles, '.obs'),  
              xlim = c(q.obs[1]-0.5, q.obs[5]+0.5),
              ylim = c(0, max(q.sim)),
              col = "lightblue",
              main = title,
              ylab = "Simulated values",
              xlab = "Observed quantiles")
      lines(q.obs, q.obs, col = "red", pch = 19, cex = 1.3, type = 'b')
      i <- i + 1
    }
  }
}

boxplot.q.sim('A126', X.MHO, df_hours, per.comun.h,
              years = 2015, seasons = T, y.sim = y.sim.exp.3.cor, plot = T)


#comparación de máximos y quantiles
library(ggplot2)
library(dplyr)
library(tidyr)
X <- X.MHO[['A126']]
X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"),
                  format = '%Y-%m-%d')
X <- X %>% filter(date >= per.comun.h[1] & date <= per.comun.h[2] 
                  & t == 2015)

library(ggplot2)
library(dplyr)
library(tidyr)

num.comparison <- function(fun.to.do, sims.list, X, var = "A126.p", quantile = NULL) {
  # Determinar si se trata de un cuantil
  is.quantile <- (fun.to.do == "quantile")
  fun <- match.fun(ifelse(is.quantile, "quantile", fun.to.do))
  
  # Filtrar solo valores positivos en las observaciones
  x_pos <- X[[var]][X[[var]] > 0]
  
  # Calcular valor observado (solo positivos)
  if (is.quantile) {
    if (is.null(quantile)) stop("Debes especificar el argumento 'quantile' cuando fun.to.do = 'quantile'")
    obs <- fun(x_pos, probs = quantile, na.rm = TRUE)
  } else {
    obs <- fun(x_pos, na.rm = TRUE)
  }
  
  # Calcular estadístico para cada conjunto de simulaciones (solo valores positivos)
  df_sim <- lapply(names(sims.list), function(name) {
    sims <- sims.list[[name]]
    sim_vals <- apply(
      sims, 2,
      function(col) {
        col_pos <- col[col > 0]
        if (is.quantile) {
          fun(col_pos, probs = quantile, na.rm = TRUE)
        } else {
          fun(col_pos, na.rm = TRUE)
        }
      }
    )
    data.frame(Simulation = 1:length(sim_vals),
               Value = sim_vals,
               Set = name)
  }) %>% bind_rows()
  
  # Crear gráfico con ggplot2
  ggplot(df_sim, aes(x = Simulation, y = Value, color = Set)) +
    geom_line(alpha = 0.8) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = obs, color = "black", linewidth = 0.8) +
    labs(
      title = ifelse(is.quantile,
                     paste0("Quantile (", quantile, ")"),
                     fun.to.do),
      x = "Simulations",
      y = "Value",
      color = "Dataset"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
}


num.comparison('max', 
               list(y.sim.og = y.sim.og,
                    y.sim.og.cor = y.sim.og.cor,
                    y.sim.exp.3 = y.sim.exp.3,
                    y.sim.exp.3.cor = y.sim.exp.3.cor), 
               X, var = 'A126.p')
num.comparison('quantile', 
               list(y.sim.og = y.sim.og,
                    y.sim.og.cor = y.sim.og.cor,
                    y.sim.exp.3 = y.sim.exp.3,
                    y.sim.exp.3.cor = y.sim.exp.3.cor), 
               X, var = 'A126.p',
               quantile = 0.9)

basura <- t(apply(y.sim.og.cor, 2, function(col) {
  pos <- col[col > 0]
  if (length(pos) > 0) {
    quantile(pos, probs = 0.9, na.rm = TRUE)
  } else {
    rep(NA, length(1))  # si no hay positivos, devolvemos NAs
  }
}))
boxplot(t(basura))
points(1, quantile(X[['A126.p']][X[['A126.p']] > 0], probs = 0.90), col = 'red', pch = 19)

# Cantidades más aplatanadas?
X <- X.MHO[[station]]
X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = '-'), 
                  format = '%Y-%m-%d')
X <- X %>% 
  filter(date >= per.comun.h[1] & date <= per.comun.h[2])

df_hours$date <- as.Date(paste(df_hours$t, df_hours$mes, df_hours$dia.mes, sep = "-"), 
                       format = '%Y-%m-%d')
data.h <- df_hours %>% filter(date >= per.comun.h[1] & date <= per.comun.h[2])

if (!is.null(years)){
  X <- X %>%
    filter(t %in% years)
  data.h <- data.h %>%
    filter(t %in% years)
}

x.obs <- X[[paste0(station, '.p')]]

#obtencion del subset de simulaciones con lluvias
aux.ind <- which(data.h$date %in% unique(X$date))
y.sim.rain <- y.sim[aux.ind, ]


plot(1:length(x.obs), x.obs, col = 'black')
