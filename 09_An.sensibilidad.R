rm(list = setdiff(ls(), c('common.models.final', 'MDQ', 'MHQ', 'X.MDQ',
                          'X.MHQ', 'estaciones', 'per.comun.day',
                          'per.comun.h')))

# Sensitivity analysis
library(qs)
X.MDQ <- qread('X.MDQ.qs')
X.MHQ <- qread('X.MHQ.qs')
X.MDO <- qread('X.MDO.qs')
X.MHO <- qread('X.MHO.qs')

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
                             period, n.sim, ocurrence = TRUE){
  
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
  aux.Xo <- Xo[, names(Xq)]
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
      prob.p = mo$fitted.values, #fitted values of model MDO(phase2) or dynamically(phase3)
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
                          quantiles, n.sim, month = NULL, 
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
                              data.mo = X.MDO,
                              data.mq = X.MDQ,
                              model.list = model.list,
                              model.q = 'MDQ',
                              model.o = 'MDO',
                              period = period,
                              n.sim = n.sim)
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
basura <- boxplot.q.sim(estaciones[1], X.MDO, common.models.final, 'MDQ', 
              per.comun.day, 
              quantiles = c(0.05, 0.50, 0.90, 0.95, 0.99),
              n.sim = 100, month = NULL,
              seasons = TRUE,
              generator = 3,
              plot = TRUE, 
              extreme.lag = TRUE,
              n.days = 8, 
              partition.lag = TRUE,
              top.quarter.z.value = TRUE)

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
