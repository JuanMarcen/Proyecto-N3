rm(list = ls())

# MHO <- readRDS('MHO.rds') #muy pesado (demasiado) #por que tengo donde más datos
# MHQ <- readRDS('MHQ.rds')
# MDO <- readRDS('MDO.rds')
# MDQ <- readRDS('MDQ.rds')

library(qs)
library(lubridate)
library(dplyr)
library(pROC)
library(gamlss)
# qsave(MHO, "MHO.qs")
# qsave(MHQ, "MHQ.qs")
# qsave(MDO, "MDO.qs")
# qsave(MDQ, "MDQ.qs")
MHO <- qread('MHO.qs')
MHQ <- qread('MHQ.qs')
MDO <- qread('MDO.qs')
MDQ <- qread('MDQ.qs')



stations <- readRDS('stations.rds')
estaciones <- stations$STAID

#only for daily
global_df <- readRDS('global_df.rds')
global_df$t <- year(global_df$date)
global_df <- global_df[, -which(colnames(global_df) %in% c('zg300.', 'zg500.', 'zg700.', 
                                                           'zt300.', 'zt500.', 'zt700.'))]

#fast check
# for (station in estaciones){
#   m <- MDO[[station]]$M5
#   X <- MDO[[station]]$X
#   
#   cat('Estacion', station, '\n')
#   cat('fv: ', length(m$fitted.values), '\n')
#   cat('X: ', dim(X)[1], '\n')
#   cat('Bien: ', length(m$fitted.values)==dim(X)[1], '\n')
# }

# all in check
#----DAILY GENERATOR IN P019----
# this station looks good in the controls carried out
# is the one with the most data available (since 1997)

sanitize_formula <- function(f) {
  ftxt <- paste(deparse(f), collapse = " ")
  ftxt <- gsub("\\b[[:alnum:]_.]+\\$", "", ftxt)
  ftxt <- gsub("poly\\(\\s*([^,\\)]+)\\s*,\\s*([0-9]+)\\s*\\)",
               "poly(\\1, \\2, raw = TRUE)", ftxt, perl = TRUE)
  as.formula(ftxt, env = environment())
}

#as function it doesnt work, so we do a foor loop deleting data
aux.mo <- 'M5'
aux.mq <- 'M6'
n.sim <- 100
set.seed(05052002)
for (station in estaciones[1]){
  mo <- MDO[[station]][[aux.mo]]
  mq <- MDQ[[station]][[aux.mq]]
  
  Xo <- MDO[[station]]$X
  Xo$date <- as.Date(paste(Xo$t, Xo$mes, Xo$dia.mes, sep = "-"), 
                     format = '%Y-%m-%d')
  
  Xq <- MDQ[[station]]$X        
  
  # Crear copias locales de las fórmulas
  mu.form <- sanitize_formula(mq$mu.formula)
  sigma.form <- sanitize_formula(mq$sigma.formula)
  
  print(mu.form)
  print(sigma.form)
  
  cat('1..')
  # Ajustar el modelo
  mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                family = GA, data = Xq, trace = FALSE)
  
  #dataframe with climatic data
  # en el fondo no neceasrio porque depende del dia anterior no?
  # solo necesario si no está la variable lag en ninguno modelos
  # igual mas cómodo generalizar
  if (length(grep(paste0(station,'.p.lag'), labels(terms(mu.form)))) != 1 & 
      length(grep(paste0(station,'.p.lag'), labels(terms(formula(mo))))) != 1){
    X_data <- global_df %>%
      filter(STAID == station) %>%
      select(date) %>%
      inner_join(Xo, by = 'date') %>%
      as.data.frame()
    X_data <- X_data[, names(Xq)] # predict in gamlss, requires the same order in the columns
    
    cat('1..')
    mu.new <- suppressWarnings(
      predict(mq2, newdata = X_data, what = 'mu', type = 'response')
    )
    cat('2..')
    sigma.new <- suppressWarnings(
      predict(mq2, newdata = X_data, what = 'sigma', type = 'response')
    )
    cat('3..')
    
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
    
    assign(paste0('df.gen.day.', station), df.gen)
    
    rm(list = c('df.gen', 'sim_df', 'gamma_matrix',
                'rain_matrix', 'n_days', 'sigma.new',
                'mu.new', 'X_data', 'Xo', 'Xq', 'mq2',
                'mu.form', 'sigma.form', 'mo', 'mq'))
  }else{
    #simular día a día
    X_data <- global_df %>%
      filter(STAID == station) %>%
      select(date) %>%
      inner_join(Xo, by = 'date') %>%
      as.data.frame()
    X_data <- X_data[, names(Xq)] # predict in gamlss, requires the same order in the columns
    
    #nueva mu y sigma necesarios para la primera fecha
    cat('Dia 1\n')
    cat('mu..')
    mu.new <- suppressWarnings(
      predict(mq2, newdata = X_data[1, ], what = 'mu', type = 'response')
    )
    cat('sigma..\n')
    sigma.new <- suppressWarnings(
      predict(mq2, newdata = X_data[1, ], what = 'sigma', type = 'response')
    )
    
    df.gen <- Xo %>% select(date)
    
    df.gen <- df.gen %>%
      mutate(
        prob.p = NA, #change this for predict of a new df?
        mu.fv = NA,
        sigma.fv = NA,
        shape.fv = NA,
        rate.fv = NA, #shape/mu
        p.day.obs = X_data[[paste0(station, '.p')]]
      )
    
    #va a ser na en el resto porque tengo 100 valores disitntos para cada 
    df.gen[1, c(2:(ncol(df.gen) - 1))] <- c(
      prob.p = mo$fitted.values[1], #change this for predict of a new df?
      mu.fv = mu.new,
      sigma.fv = sigma.new,
      shape.fv = 1 / sigma.new^2,
      rate.fv = 1 / sigma.new^2 / mu.new
    )
    
    #simulación dia 1
    rain_matrix <- matrix(runif(1 * n.sim), nrow = 1, ncol = n.sim) <= df.gen$prob.p[1]
    
    # Matriz de gamma (simulaciones de lluvia)
    gamma_matrix <- matrix(rgamma(1 * n.sim, shape = rep(df.gen$shape.fv[1], n.sim),
                                  rate = rep(df.gen$rate.fv[1], n.sim)), 
                           nrow = 1, ncol = n.sim)
    
    # Asignamos 0 donde no llueve
    gamma_matrix[!rain_matrix] <- 0
    
    # Convertimos la matriz a columnas en df.gen
    sim_df <- as.data.frame(gamma_matrix)
    names(sim_df) <- paste0("p.day.sim.", 1:n.sim)
    
    sim_df <- sim_df[rep(1, nrow(df.gen)), , drop = F]  # repite la fila
    sim_df[-1, ] <- NA                    # ponemos NA en todas menos la primera
    
    # Unir columnas
    df.gen <- cbind(df.gen, sim_df)
    
    for (i in 2:nrow(df.gen)){
      cat('Dia ', i, ' \n')
      # prueba para el día 2 (generalizado)
      # crear un nuevo data frame de 100 (simulaciones) x columnas X_data
      aux <- X_data[i, , drop = FALSE]
      aux <-aux[rep(1, n.sim), , drop = FALSE]
      aux[[paste0(station, '.p.lag')]] <- as.numeric(df.gen[i - 1, paste0("p.day.sim.", 1:n.sim)])
      
      aux <- aux[, names(Xq)]
      aux2 <- aux
      #mirar tema NA's
      #obtención de mu y sigma para cada uno de estos
      cat('mu..')
      aux2$mu.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'mu', type = 'response')
      )
      cat('sigma..')
      aux2$sigma.new <- suppressWarnings(
        predict(mq2, newdata = aux, what = 'sigma', type = 'response')
      )
      
      #obtencion de probabilidad de lluvia
      cat('prob..\n')
      aux2$prob.dia <- predict(mo, newdata = aux, type = 'response')
      
      aux2 <- aux2 %>% 
        select(prob.dia, mu.new, sigma.new) %>%
        rowwise() %>%
        mutate(
          shape.fv = 1 / sigma.new^2,
          rate.fv = 1 /sigma.new^2 / mu.new
        ) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          sim = ifelse(runif(1) <= prob.dia, rgamma(1, shape = shape.fv, rate = rate.fv), 0)
        ) %>% 
        ungroup()
      
      df.gen[i, paste0("p.day.sim.", 1:n.sim)] <- as.numeric(aux2$sim)
    }
    

  }

  
}


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
  
  # Crear copias locales de las fórmulas
  mu.form <- sanitize_formula(mq$mu.formula)
  sigma.form <- sanitize_formula(mq$sigma.formula)
  
  print(mu.form)
  print(sigma.form)
  
  cat('1..')
  # Ajustar el modelo
  mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
                family = GA, data = Xq, trace = FALSE)
  
  #dataframe with climatic data
  X_data <- data.df %>%
    filter(STAID == station) %>%
    select(date) %>%
    inner_join(Xo, by = 'date') %>%
    as.data.frame()
  X_data <- X_data[, names(Xq)] # predict in gamlss, requires the same order in the columns
 
  cat('1..')
  mu.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'mu', type = 'response')
  )
  cat('2..')
  sigma.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'sigma', type = 'response')
  )
  cat('3..')

  
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


# NU FUNCIONA, MIRAR DEL TEMA DE ENTORNOS NS 
# cunaod lo hago dentro con una estacion funciona, cuando cambio de esacion ya no 
# basura <- daily.generator(station = estaciones[2], 
#                           data.mo = MDO, mo = 'M5', 
#                           data.mq = MDQ, mq = 'M6', 
#                           data.df = global_df, 
#                           n.sim = 100)

basura <- get(paste0('df.gen.day.',station))
plot(density(basura[['p.day.obs']]), col = 'blue', lwd = 2,
     xlim = c(0, max(basura[,paste0('p.day.sim.', 1:n.sim)])))
for(i in 1:n.sim){
  lines(
    density(basura[[paste0('p.day.sim.', i)]]), 
    col = 'red', lwd = 2
    )
}
lines(density(basura[['p.day.obs']]), col = 'blue', lwd = 2)

max <- apply(basura[, c('p.day.obs', paste0('p.day.sim.', 1:n.sim))], 2, max)
boxplot(max)
points(max[1], pch = 19, col = 'red')
abline(h = quantile(max, probs = 0.90), lty = 2, col = 'blue')
abline(h = quantile(max, probs = 0.95), lty = 2, col = 'blue')
abline(h = quantile(max, probs = 0.05), lty = 2, col = 'blue')
abline(h = quantile(max, probs = 0.1), lty = 2, col = 'blue')

#----HOURLY GENERATOR----
hourly.generator <- function(station, 
                            data.mo, mo, 
                            data.mq, mq, 
                            n.sim){
  cat('1..')
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
  cat('2..')
  Xq <- data.mq[[station]]$X
  cat('3..')
  Xq$datetime <- ISOdatetime(year = Xq$t,
                             month = Xq$mes,
                             day = Xq$dia.mes,
                             hour = Xq$h,
                             min = 0,
                             sec = 0,
                             tz = "UTC")
  
  # in order to predict correctly (USE THIS ONE TO PREDICT)
  environment(mq$mu.formula) <- environment()
  environment(mq$sigma.formula) <- environment()
  mq2 <- gamlss(mq$mu.formula, 
                sigma.fo = mq$sigma.formula, 
                family = GA, 
                data = Xq, trace = F)
  cat('4..')
  
  #sum(abs(mq$mu.fv - mq2$mu.fv) < 1e-6)
  
  # X_data: smae columns as Xq
  
  X_data <- Xo[, names(Xq)]
  cat('5..')
  
  #new values of mu and sigma 
  mu.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'mu', type = 'response')
  )
  
  cat('6..')
  sigma.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'sigma', type = 'response')
  )
  
  df.gen <- Xo %>%
    select(datetime)
    
  
  df.gen <- df.gen %>%
    mutate(
      prob.p = mo$fitted.values, #change this for predict of a new df?
      mu.fv = mu.new,
      sigma.fv = sigma.new,
      shape.fv = 1 / sigma.new^2,
      rate.fv = 1 / sigma.new^2 / mu.new, #shape/mu
      p.h.obs = Xo[[paste0(station, '.p')]]
    )
  
  # n.simulations
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
  names(sim_df) <- paste0("p.h.sim.", 1:n.sim)
  
  df.gen <- cbind(df.gen, sim_df)
  
  return(df.gen)
}


aux.mo <- 'M5'
aux.mq <- 'M6'
n.sim <- 100
for (station in estaciones[1]){
  cat('1..')
  mo <- MHO[[station]][[aux.mo]]
  mq <- MHQ[[station]][[aux.mq]]
  
  Xo <- MHO[[station]]$X
  Xo$datetime <- ISOdatetime(year = Xo$t,
                             month = Xo$mes,
                             day = Xo$dia.mes,
                             hour = Xo$h,
                             min = 0,
                             sec = 0,
                             tz = "UTC")
  cat('2..')
  Xq <- MHQ[[station]]$X
  cat('3..')
  Xq$datetime <- ISOdatetime(year = Xq$t,
                             month = Xq$mes,
                             day = Xq$dia.mes,
                             hour = Xq$h,
                             min = 0,
                             sec = 0,
                             tz = "UTC")
  
  # in order to predict correctly (USE THIS ONE TO PREDICT)
  mu.form <- sanitize_formula(mq$mu.formula)
  sigma.form <- sanitize_formula(mq$sigma.formula)
  
  print(mu.form)
  print(sigma.form)
  
  mq2 <- gamlss(mu.form, 
                sigma.fo = sigma.form, 
                family = GA, 
                data = Xq, trace = F)
  cat('4..')
  
  #sum(abs(mq$mu.fv - mq2$mu.fv) < 1e-6)
  
  # X_data: smae columns as Xq
  
  X_data <- Xo[, names(Xq)]
  cat('5..')
  
  #new values of mu and sigma 
  mu.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'mu', type = 'response')
  )
  
  cat('6..')
  sigma.new <- suppressWarnings(
    predict(mq2, newdata = X_data, what = 'sigma', type = 'response')
  )
  
  df.gen <- Xo %>%
    select(datetime)
  
  
  df.gen <- df.gen %>%
    mutate(
      prob.p = mo$fitted.values, #change this for predict of a new df?
      mu.fv = mu.new,
      sigma.fv = sigma.new,
      shape.fv = 1 / sigma.new^2,
      rate.fv = 1 / sigma.new^2 / mu.new, #shape/mu
      p.h.obs = Xo[[paste0(station, '.p')]]
    )
  
  # n.simulations
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
  names(sim_df) <- paste0("p.h.sim.", 1:n.sim)
  
  df.gen <- cbind(df.gen, sim_df)
  
  assign(paste0('df.gen.h.', station), df.gen)
  
  rm(list = c('df.gen', 'sim_df', 'gamma_matrix',
              'rain_matrix', 'n_days', 'sigma.new',
              'mu.new', 'X_data', 'Xo', 'Xq', 'mq2',
              'mu.form', 'sigma.form', 'mo', 'mq'))
}


basura2 <- hourly.generator(station = estaciones[1],
                            data.mo = MHO, mo = 'M5',
                            data.mq = MHQ, mq = 'M6', 
                            n.sim = 100)

basura2 <- get(paste0('df.gen.h.',station))
plot(density(basura2[['p.h.obs']]), col = 'blue', lwd = 2, 
     xlim = c(0, max(basura2[,paste0('p.h.sim.', 1:n.sim)])))
for(i in 1:n.sim){
  lines(
    density(basura2[[paste0('p.h.sim.', i)]]), 
    col = 'red', lwd = 2
  )
}
lines(density(basura2[['p.h.obs']]), col = 'blue', lwd = 2)
max <- apply(basura2[, c('p.h.obs', paste0('p.h.sim.', 1:n.sim))], 2, max)
boxplot(max)
points(max[1], pch = 19, col = 'red')
abline(h = quantile(max, probs = 0.90), lty = 2, col = 'blue')
abline(h = quantile(max, probs = 0.95), lty = 2, col = 'blue')
abline(h = quantile(max, probs = 0.05), lty = 2, col = 'blue')
abline(h = quantile(max, probs = 0.1), lty = 2, col = 'blue')

#HACE FATAL LOS MÁXIMOS... PORQUE NO HEMOS PUESTO LA RESTRICCION DE QUE LA SUMA HORARIA DEBER IGUAL A LA DIARIA