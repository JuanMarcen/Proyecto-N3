rm(list = ls())

library(qs)
MDO <- qread('MDO.qs')
MDQ <- qread('MDQ.qs')
MHO <- qread('MHO.qs')
MHQ <- qread('MHQ.qs')
stations <- readRDS('stations.rds')
estaciones <- stations$STAID
# selección (manual) del modelo que quiero para cada uno
# criterio que sigo

# para el caso de los de ocurrencia: AUC bueno, num.covariables pequeño
# para los de cantidad: Num.covariables pequeño ?
library(pROC)
library(ggplot2)

# MODELOS DE OCURRENCIA
plot.auc.n.vars <- function(mo, estaciones, data){
  df <- data.frame(
    station = estaciones
  )
  
  for (i in 1:length(estaciones)){
    station <- estaciones[i]
    mod <- data[[station]][[mo]]
    X <- data[[station]]$X
    
    roc <- suppressMessages(
      roc(X$Y, predict(mod, type = 'response'))
    )
    
    df[i, 'AUC'] <- auc(roc)
    
    df[i, 'n.vars'] <- length(data[[station]][[paste0('vars.', mo)]]) - 1
    
    
  }
  
  g <- ggplot(aes(x = AUC, y = n.vars), data = df) + 
    geom_point(size = 3) + 
    geom_text(aes(label = 1:nrow(df)), vjust = -0.5, size = 5) +
    scale_y_continuous(breaks = seq(min(df$n.vars), max(df$n.vars), by = 1))+ 
    labs(x = "Área bajo la curva (AUC)", y = "Número de variables")
  
  print(g)
  
  return(df)
}

mo <- 'M9'
df.mho <- plot.auc.n.vars('M5', estaciones, MHO)
# MHO escogido == 19 --> df.comp
# MHO escogido == 1 --> df.comp.2
df.mdo <- plot.auc.n.vars('M10', estaciones, MDO)
# MDO escogido == 11

# MODELOS DE CANTIDAD
n.vars.mdq <- c()
n.vars.mhq <- c()
mq <- 'M8'
for (station in estaciones){
  #n.vars.mdq <- c(n.vars.mdq, length(MDQ[[station]][[paste0('vars.', mq)]]) -1)
  n.vars.mhq <- c(n.vars.mhq, length(MHQ[[station]][[paste0('vars.', mq)]]) -1)
}
which.min(n.vars.mhq)
# MHQ escogido == 15
which.min(n.vars.mdq)
#MDQ escogido == 8

#----Periodo comun (distinto en horas y días)----
# sabemos que los diarios como mucho a 2023
load('data.RData')
global_df <- readRDS('global_df.rds')

library(lubridate)
n <- length(estaciones)
aux <- apply(df_hours[, paste0(estaciones, '.p')], 1, function(row) all(!is.na(row)))
first <- which(aux)[1]
last <- tail(which(aux), 1)
date.first <- as.Date(paste0(df_hours[first, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
date.last <- as.Date(paste0(df_hours[last, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
per.comun.h <- c(date.first, date.last)


aux <- apply(df_days[, paste0(estaciones, '.p')], 1, function(row) all(!is.na(row)))
first <- which(aux)[1]
date.first <- as.Date(paste0(df_days[first, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
aux2 <- apply(global_df[, (8:ncol(global_df))], 1, function(row) all(!is.na(row)))
last <- tail(which(aux2), 1)
date.last <- global_df$date[last]

per.comun.day <- c(date.first, date.last)

#----Ajuste modelos en periodo comun----
# ajuste del modelo escogido
# ajuste según selección en el periodo común
# comparar entre ellos

# ajuste de los modelos comunes en todo el periodo comun
library(dplyr)
library(gamlss)
#library(logistf)
X.MDO <- qread('X.MDO.qs')
mod.comun <- function(station, per.comun, data, mod, chosen, 
                      type = 'log', subtype = NULL, log_file = "errores_modelo.txt",
                      data.used = NULL) {
  
  data_name <- deparse(substitute(data))
  
  # Inicializamos o añadimos cabecera
  header <- paste0("=== Registro de ADVERTENCIAS en mod.comun ", data_name, " ===\n")
  if(!file.exists(log_file)) {
    writeLines(header, log_file)
  } else {
    cat(header, file = log_file, append = TRUE)
  }
  
  # auxiliar para registrar advertencias
  log_problem <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- paste0("[", timestamp, "] Estación: ", station, " - ", msg, "\n")
    cat(line, file = log_file, append = TRUE)
  }
  
  
  # PREPROCESADO
  if (!is.null(data.used)){
    X <- data.used[[station]]
  }else{
    X <- data[[station]]$X
  }
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = '%Y-%m-%d')
  X <- X %>% filter(date >= per.comun[1] & date <= per.comun[2])
  
  mod.per.comun <- NULL
  
  # --- Ajuste del modelo ---
  if (type == 'log') {
    vars <- labels(terms(data[[estaciones[chosen]]][[mod]]$formula))
    vars <- gsub(estaciones[chosen], station, vars)
    formula <- as.formula(paste('Y ~', paste(vars, collapse = '+')))
    # formula modification for no divergence
    # formula <- update(
    #   formula,
    #   as.formula(
    #     paste0(". ~ . - poly(", station, ".p.lag, 3) + ", station, ".p.lag")
    #   )
    # )
    #formula modification for no divergence
    # formula <- update(
    #   formula,
    #   as.formula(
    #     paste0(". ~ . - s.2.h - c.2.h")
    #   )
    # )
    
    mod.per.comun <- tryCatch({
      withCallingHandlers(
        glm(formula, family = binomial(link = "logit"), data = X,
            control = glm.control(epsilon = 1e-8, maxit = 100, trace = FALSE)),
        warning = function(w) {
          # registramos la advertencia
          log_problem(paste("ADVERTENCIA (GLM):", conditionMessage(w)))
          # intentamos mufflear, pero no fallaremos si no existe el restart
          tryCatch(invokeRestart("muffleWarning"), error = function(e) NULL)
        }
      )
    },
    error = function(e) {
      # Ignoramos errores completamente (no se imprimen ni registran)
      return(NULL)
    })
    
  } else if (type == 'gamma') {
    vars <- labels(terms(data[[estaciones[chosen]]][[mod]]$mu.formula))
    vars <- gsub(estaciones[chosen], station, vars)
    formula <- as.formula(paste(paste0(station, '.p'), '~', paste(vars, collapse = '+')))
    
    mod.per.comun <- tryCatch({
      withCallingHandlers(
        {
          if (subtype == 'day') {
            gamlss(formula, 
                   sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)), 
                   family = GA, data = X, trace = FALSE)
          } else if (subtype == 'hour') {
            gamlss(formula, 
                   sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                     I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
                   family = GA, data = X, trace = FALSE)
          } else {
            stop("Subtype not valid. Use 'day' or 'hour'.")
          }
        },
        warning = function(w) {
          log_problem(paste("ADVERTENCIA (GAMLSS):", conditionMessage(w)))
          tryCatch(invokeRestart("muffleWarning"), error = function(e) NULL)
        }
      )
    },
    error = function(e) {
      return(NULL)
    })
    
  } else {
    stop("Type not valid. Use 'log' or 'gamma'.")
  }
  
  return(mod.per.comun)
}




common.models <- list()
for(station in estaciones){
  cat('Estación ', station, '\n')
  #common.models[[station]][['MHO.pc']] <- mod.comun(station, per.comun.h, MHO, 'M5', 19, type = 'log')
  #common.models[[station]][['MHO.pc.2']] <- mod.comun(station, per.comun.h, MHO, 'M5', 1, type = 'log')
  common.models[[station]][['MHO.pc.3']] <- mod.comun(station, per.comun.h, MHO, 'M9', 19, type = 'log')
  # common.models[[station]][['MDO.pc']] <- mod.comun(station, per.comun.day, MDO, 'M5', 11, type = 'log')
  # common.models[[station]][['MDO.pc.2']] <- mod.comun(station, per.comun.day, MDO, 'M7', 11, type = 'log')
  common.models[[station]][['MDO.pc.3']] <- mod.comun(station, per.comun.day, MDO, 'M8', 11, type = 'log')
  common.models[[station]][['MDO.pc.4']] <- mod.comun(station, per.comun.day, MDO, 'M10', 19, type = 'log',
                                                      data.used = X.MDO)
  # common.models[[station]][['MHQ.pc']] <- mod.comun(station, per.comun.h, MHQ, 'M6', 15, type = 'gamma', subtype = 'hour')
  common.models[[station]][['MHQ.pc.2']] <- mod.comun(station, per.comun.h, MHQ, 'M8', 6, type = 'gamma', subtype = 'hour')
  # common.models[[station]][['MHQ.pc.3']] <- mod.comun(station, per.comun.h, MHQ, 'M8', 8, type = 'gamma', subtype = 'hour')
  # common.models[[station]][['MHQ.pc.4']] <- mod.comun(station, per.comun.h, MHQ, 'M8', 11, type = 'gamma', subtype = 'hour')
  common.models[[station]][['MDQ.pc']] <- mod.comun(station, per.comun.day, MDQ, 'M6', 8, type = 'gamma', subtype = 'day')
}

# DIVERGENT STATIONS
mho.div.stations <- c('E085', 'R026', 'R062', 'P088', 'A042', 'P023', 'P024', 'A266', 'A058', 'EM12', 'P019')
mho.2.div.stations <- c('E085', 'R026', 'R062', 'P088', 'P023', 'P024', 'A266', 'A126', 'A058', 'EM12', 'P019')
mdo.div.stations <- c('R037')
#trasnformadno a linela el lag no hay errores en MDO!!

vars.div <- function(stations, models, model.type){
  vars.div.list <- list()
  for (station in stations){
    data.aux <- models[[station]][[model.type]]$data
    formula.aux <- models[[station]][[model.type]]$formula
    mod.aux <- glm(formula.aux, data = data.aux, family = binomial(logit))
    mod.aux.2 <- step(mod.aux)
    
    vars.div.list[[station]] <- setdiff(labels(terms(mod.aux$formula)), labels(terms(mod.aux.2$formula)))
  }
  
  return(vars.div.list)
}
vars.div.mho <- vars.div(mho.div.stations, common.models, 'MHO.pc')
vars.div.mho.2 <- vars.div(mho.2.div.stations, common.models, 'MHO.pc.2')
vars.div.mdo <- vars.div(mdo.div.stations, common.models, 'MDO.pc')

# qsave(vars.div.mho, 'vars.div.mho.qs')
# qsave(vars.div.mho.2, 'vars.div.mho.2.qs')
# qsave(vars.div.mdo, 'vars.div.mdo.qs')

vars.div.mho <- qread('vars.div.mho.qs')
vars.div.mho.2 <- qread('vars.div.mho.2.qs')
vars.div.mdo <- qread('vars.div.mdo.qs')

station <- 'R037' #P021 #R062
station.p <- paste0(station, '.p')
data.aux <- common.models[[station]][['MDO.pc']]$data
max(data.aux[[station.p]])
summary(data.aux[[paste0(station.p, '.day')]])

formula.aux <- common.models[[station]][['MDO.pc']]$formula
print(formula.aux)
mod.aux <- glm(formula.aux, data = data.aux, family = binomial(link = 'logit'))
mod.aux
library(gam)
basura <- gam(formula = Y ~ s(R037.p.lag) , data = data.aux, family = binomial)
basura1 <- plot(basura$data$R037.p.lag, basura$fitted.values)
plot(basura)

library(logistf)
mod.aux.firth <- logistf(formula.aux, data.aux)
mod.aux.firth
-2 * mod.aux.firth$loglik['full'] + 2 * length(mod.aux.firth$coefficients)
# mod.aux.2 <- update(mod.aux, formula = as.formula(
#   paste('.~. - ', paste0('poly(', station, '.p.lag, 3)'), '+', paste0(station, '.p.lag'),
#         '-', paste0('poly(', station, '.p.day, 3)'), '+', paste0(station, '.p.day'))
# ))

aux <- data.frame(dfbetas(mod.aux))
thresh <- 2 / sqrt(nrow(data.aux))
vars <- names(aux)[-1]
par(mfrow = c(2,4))
plot(aux[, vars[1]], main = vars[1], ylab = '')
abline(h = thresh, col = 'red')
abline(h = -thresh, col = 'red')
for (i in 2:length(vars)){
  plot(aux[, vars[i]], ylab = '', main = vars[i])
  abline(h = thresh, col = 'red')
  abline(h = -thresh, col = 'red')
}

data.aux[which.min(aux$poly.R037.p.lag..3.1), paste0(station.p, '.lag')]

#effects in each station
library(gam)
for(station in estaciones){
  station.p <- paste0(station, '.p')
  X <- X.MHO[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = '%Y-%m-%d')
  data.aux <- X %>% filter(date >= per.comun.h[1] & date <= per.comun.h[2])
  
  basura <- gam(formula = as.formula(paste0(station.p, ' ~ s(', station.p, '.lag)')) , 
                data = as.data.frame(data.aux), family = Gamma(link = 'log'))
  basura1 <- plot(basura$data[[paste0(station.p, '.lag')]], basura$fitted.values)
}


#----ajuste de modelos por selección AIC en el periodo común----
# recuperar funciones de selección
# selección ocurrencia
harmonics.l <- list(
  h1 = c('s.1.l', 'c.1.l'),
  h2 = c('s.2.l', 'c.2.l'),
  h3 = c('s.3.l', 'c.3.l'),
  h4 = c('s.4.l', 'c.4.l')
)
harmonics.h <- list(
  h1 = c('s.1.h', 'c.1.h'),
  h2 = c('s.2.h', 'c.2.h'),
  h3 = c('s.3.h', 'c.3.h'),
  h4 = c('s.4.h', 'c.4.h')
)
step_rlog <- function(initial_model,
                      data,
                      vars,
                      harmonics.l,
                      harmonics.h, step = 'AIC'){ 
  
  if (step == 'AIC'){
    k <- 2
  }else if (step == 'BIC'){
    n <- dim(data)[1]
    k <- log(n)
  }
  
  mod.aux <- initial_model
  cat('Initial model ', deparse(formula(mod.aux)), '\n')
  cat(step, ': ', AIC(mod.aux, k = k), '\n')
  
  scope.aux <- update(formula(mod.aux), paste('. ~ . +', paste(vars, collapse = '+')))
  
  
  mod.aux <- step(mod.aux, scope = scope.aux, direction = 'both', trace = FALSE, k = k)
  
  cat('Model after ', step, 'step algorithm for variables: ', deparse(formula(mod.aux)), '\n')
  cat(step, ': ', AIC(mod.aux, k = k), '\n')
  
  for (h in harmonics.l){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- glm(formula.aux, data = data, family = binomial(logit))
    
    if (AIC(mod.temp, k = k) < AIC(mod.aux, k = k)){
      cat('Added: ',paste(h, collapse = '+'), '\n')
      cat(step, ': ', AIC(mod.temp, k = k), '\n')
      
      mod.aux <- mod.temp
    }else{
      break
    }
  }
  
  if (!is.null(harmonics.h)){
    for (h in harmonics.h){
      
      formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
      
      mod.temp <- glm(formula.aux, data = data, family = binomial(logit))
      
      if (AIC(mod.temp, k = k) < AIC(mod.aux, k = k)){
        cat('Added: ',paste(h, collapse = '+'), '\n')
        cat(step, ': ', AIC(mod.temp, k = k), '\n')
        
        mod.aux <- mod.temp
      }else{
        break
      }
    }
  }
  
  
  cat('\nFinal model: ', deparse(formula(mod.aux)), '\n')
  cat(step, ': ', AIC(mod.aux, k = k), '\n\n')
  
  return(mod.aux)
}

step_glm <- function(initial_model,
                     data,
                     vars,
                     harmonics.l, 
                     harmonics.h ){ 
  
  
  # if(type == 'day'){
  #   sigma.formula <- as.formula('~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365))')
  # }else if(type == 'hour'){
  #   sigma.formula <- as.formula('~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
  #                               I(sin(2*pi*h/24)) + I(cos(2*pi*h/24))')
  # }else{
  #   stop("type not valis. Use 'day' or 'hour'.")
  # }
  
  # print(sigma.formula)
  
  mod.aux <- initial_model
  cat('Initial model ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n')
  
  scope.aux <- update(formula(mod.aux), paste('. ~ . +', paste(vars, collapse = '+')))
  
  mod.aux <- step(mod.aux, scope = scope.aux, direction = 'both', trace = F)
  
  cat('Model after ', 'AIC ', 'step algorithm for variables: ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n')
  
  
  for (h in harmonics.l){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- gamlss(formula.aux, 
                       sigma.fo = mod_null$sigma.formula, 
                       family = GA, 
                       data = data,
                       trace = F)
    
    if (AIC(mod.temp) < AIC(mod.aux)){
      cat('Added: ',paste(h, collapse = '+'), '\n')
      cat('AIC: ', AIC(mod.temp), '\n')
      
      mod.aux <- mod.temp
    }else{
      break
    }
  }
  
  if (!is.null(harmonics.h)){
    for (h in harmonics.h){
      
      formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
      
      mod.temp <- gamlss(formula.aux, 
                         sigma.fo = mod_null$sigma.formula, 
                         family = GA, 
                         data = data,
                         trace = F)
      
      if (AIC(mod.temp) < AIC(mod.aux)){
        cat('Added: ',paste(h, collapse = '+'), '\n')
        cat('AIC: ', AIC(mod.temp), '\n')
        
        mod.aux <- mod.temp
      }else{
        break
      }
    }
  }
  
  
  
  cat('\nFinal model: ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n\n')
  
  return(mod.aux)
}

# seleccion en todas estaciones
# degrees of polynomials
degrees.p.lag <- c(3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(degrees.p.lag) <- estaciones
degrees.p.day <- c(3, 1, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(degrees.p.day) <- estaciones

#grados cantidad
degrees.p.lag2 <- c(3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 2, 3,
                    3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3)
names(degrees.p.lag2) <- estaciones
degrees.p.day2 <- c(1, 3, 3, 3, 3, 3, 3, 3, 1, 3, 1, 1, 1, 1, 1, 3,
                    3, 2, 2, 3, 1, 1, 3, 1, 1, 2, 1, 1)

names(degrees.p.day2) <- estaciones

deg_list <- readRDS('deg_list.rds')
#degrees lag of daily models
deg_lag <- c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
             3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(deg_lag) <- estaciones

for (station in estaciones){
  cat('Ajuste modelos estación ', station, '\n')
  #ocurrencia horaria
  cat('\nModelo de ocurrencia horario\n')
  X <- MHO[[station]]$X
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                    format = '%Y-%m-%d')
  
  ind <- which(X$date >= per.comun.h[1] & X$date <= per.comun.h[2])
  
  deg.lag <- degrees.p.lag[station]
  deg.day <- degrees.p.day[station]
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = MHO[[station]]$X[ind, ])
  
  common.models[[station]][['MHO.pc.sel']] <- step_rlog(mod_null, 
                                              data = MHO[[station]]$X[ind, ], 
                                              vars = c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                                       paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                                       paste0(station,'.p.lag:',station,'.p.day')), 
                                              harmonics.l, harmonics.h)
  
  #ocurrencia diaria
  cat('\nModelo ocurrencia diario \n')
  X <- MDO[[station]]$X
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                    format = '%Y-%m-%d')
  ind <- which(X$date >= per.comun.day[1] & X$date <= per.comun.day[2])
  
  #deg.list <- deg_list[[station]]
  deg.lag <- deg_lag[station]
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = MDO[[station]]$X[ind, ])
  
  # aux <- paste0('poly(', names(deg.list)[1], ', ', deg.list[1], ')')
  # for (i in 2:length(deg.list)){
  #   aux <- c(aux,  paste0('poly(', names(deg.list)[i], ', ', deg.list[i], ')'))
  # }
  
  # list[[station]][['MDO.pc.sel']] <- step_rlog(mod_null,
  #                                                       data = X,
  #                                                       vars = c(aux,
  #                                                                paste0('poly(', station, '.p.lag, ', deg.lag, ')')),
  #                                                       harmonics.l, harmonics.h)
  common.models[[station]][['MDO.pc.sel']] <- step_rlog(mod_null,
                                               data = MDO[[station]]$X[ind, ],
                                               vars = c(colnames(X)[16:(ncol(X) - 1)],
                                                        paste0('poly(', station, '.p.lag, ', deg.lag, ')')),
                                               harmonics.l, harmonics.h = NULL)
  
  
  # cantidad horario
  cat('\nModelo de cantidad horario\n')
  X <- MHQ[[station]]$X
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                    format = '%Y-%m-%d')
  
  ind <- which(X$date >= per.comun.h[1] & X$date <= per.comun.h[2])
  
  deg.lag <- degrees.p.lag2[station]
  deg.day <- degrees.p.day2[station]
  
  mod_null <- gamlss(as.formula(paste(paste0(station, '.p'), '~ 1')), 
                     sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                       I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
                     family = GA, 
                     data = MHQ[[station]]$X[ind, ],
                     trace = F)
  
  common.models[[station]][['MHQ.pc.sel']] <- step_glm(mod_null, 
                                            data = MHQ[[station]]$X[ind, ], 
                                            vars = c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                                     paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                                     paste0(station,'.p.lag:',station,'.p.day')), 
                                            harmonics.l, harmonics.h)
  
  #cantidad diario
  cat('\nModelo de cantidad diario\n')
  X <- MDQ[[station]]$X
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                    format = '%Y-%m-%d')
  
  ind <- which(X$date >= per.comun.day[1] & X$date <= per.comun.day[2])
  
  mod_null <- gamlss(as.formula(paste(paste0(station, '.p'), '~ 1')), 
                     sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)), 
                     family = GA, 
                     data = MDQ[[station]]$X[ind, ],
                     trace = F)
  
  common.models[[station]][['MDQ.pc.sel']] <- step_glm(mod_null, 
                                                       data = MDQ[[station]]$X[ind, ], 
                                                       vars = c(colnames(X)[16:(ncol(X) - 1)],
                                                               paste0('poly(', station, '.p.lag, ', deg.lag, ')')), 
                                                       harmonics.l, harmonics.h = NULL)
}

# qsave(common.models, 'common.models.qs', preset = 'archive', compress_level = 22)
# 
# common.models <- qread('common.models.qs')




#----comparacion modelos----
library(pROC)
library(ggplot2)

comp.models <- function(estaciones, mod.type, models.list, period, common, sel, AUC = NULL){
  df <- data.frame(
    station = estaciones
  )
  
  for (i in 1:length(estaciones)){
    station <- estaciones[i]
    
    X <- mod.type[[station]]$X
    X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                      format = '%Y-%m-%d')
    
    ind <- which(X$date >= period[1] & X$date <= period[2])
    
    mod.common <- models.list[[station]][[common]]
    mod.sel <- models.list[[station]][[sel]]
    
    if (!is.null(AUC)){
      roc.mod.common <- suppressMessages(
        roc(X$Y[ind], predict(mod.common, type = 'response'))
        )
      roc.mod.sel <- suppressMessages(
        roc(X$Y[ind], predict(mod.sel, type = 'response'))
      )
      
      df[i, 'AUC.mod.common'] <- auc(roc.mod.common)
      df[i, 'AUC.mod.sel'] <- auc(roc.mod.sel)
  
    }
    
    df[i, 'AIC.mod.common'] <- AIC(mod.common)
    df[i, 'AIC.mod.sel'] <- AIC(mod.sel)
  }
  
  
  return(df)
  
}


df.MHO <- comp.models(estaciones, MHO, common.models, per.comun.h, 'MHO.pc', 'MHO.pc.sel', AUC = TRUE)
df.MDO <- comp.models(estaciones, MDO, common.models, per.comun.day, 'MDO.pc', 'MDO.pc.sel', AUC = TRUE)
df.MHQ <- comp.models(estaciones, MHQ, common.models, per.comun.h, 'MHQ.pc', 'MHQ.pc.sel')
df.MDQ <- comp.models(estaciones, MDQ, common.models, per.comun.day, 'MDQ.pc', 'MDQ.pc.sel')

#guardado de cosas 
df.comp <- list(
  MHO = df.MHO,
  MDO = df.MDO,
  MHQ = df.MHQ,
  MDQ = df.MDQ 
)

qsave(df.comp, 'df.comp.qs')
df.comp <- qread('df.comp.qs')

#construcción de otra comparación (SOLO MHO) y juntar al df.comp grande
comp.models.2 <- function(estaciones, mod.type, models.list, period, common, AUC = NULL){
  df <- data.frame(
    station = estaciones
  )
  
  for (i in 1:length(estaciones)){
    station <- estaciones[i]
    
    X <- mod.type[[station]]$X
    X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                      format = '%Y-%m-%d')
    
    ind <- which(X$date >= period[1] & X$date <= period[2])
    
    mod.common <- models.list[[station]][[common]]
    
    if (!is.null(AUC)){
      roc.mod.common <- suppressMessages(
        roc(X$Y[ind], predict(mod.common, type = 'response'))
      )
      
      df[i, 'AUC.mod.common.3'] <- auc(roc.mod.common)
      
    }
    
    df[i, 'AIC.mod.common.3'] <- AIC(mod.common)
  }
  
  
  return(df)
  
}

df.MHQ.2 <- comp.models.2(estaciones, MHQ, common.models, per.comun.h, 'MHQ.pc.4', AUC = NULL)
df.comp$MHQ <- cbind(df.comp$MHQ, df.MHQ.2[, -1])
df.comp$MHQ <- df.comp$MHQ[, c(1,2,3,4,6,5)]
colnames(df.comp$MHQ)[5] <- 'AIC.mod.common.4'

df.MDO.4 <- comp.models.2(estaciones, MDO, common.models, per.comun.day, 'MDO.pc.4', AUC = T)
df.comp$MDO <- cbind(df.comp$MDO, df.MDO.4[, -1])
df.comp$MDO <- df.comp$MDO[, c(1,2,3,4, 10,5,6,7,8,11,9)] 

qsave(df.comp, 'df.comp.qs')
df.comp <- qread('df.comp.qs')

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggnewscale)

heat.map <- function(df, n.metrics = 1, title){
  scale_row_values <- function(values){
    rng <- range(values, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(values))) # valor medio si todo igual
    scales::rescale(values, from = rng)
  }
  if (n.metrics == 2){
    df_long <- df %>%
      pivot_longer(cols = -station, names_to = "variable", values_to = "value") %>%
      mutate(Tipo = ifelse(grepl("AUC", variable), "AUC", "AIC")) %>%
      left_join(df %>% select(station, AIC.mod.sel), by = "station") %>%
      mutate(value_rel = ifelse(Tipo == "AIC", value / AIC.mod.sel, value)) %>%
      group_by(station, Tipo) %>%
      mutate(value_row_rel = scale_row_values(value_rel)) %>%
      ungroup() %>%
      mutate(
        station = factor(station, levels = unique(station)),
        variable = factor(variable, levels = unique(variable))
      )
    
    ggplot() +
      # AUC
      geom_tile(data = df_long %>% filter(Tipo == "AUC"),
                aes(x = variable, y = station, fill = value), color = "white") +
      geom_text(data = df_long %>% filter(Tipo == 'AUC'),
                aes(x = variable, y = station, label = round(value, 3))) +
      scale_fill_gradientn(
        colors = c("blue", "lightblue", "orange", "red"),
        limits = c(min(df_long[which(df_long$Tipo == 'AUC'), 'value']),max(df_long[which(df_long$Tipo == 'AUC'), 'value'])),
        name = "AUC",
      ) +
      new_scale_fill() +  # reinicia escala
      
      # AIC
      geom_tile(data = df_long %>% filter(Tipo == "AIC"),
                aes(x = variable, y = station, fill = value_row_rel), color = "white") +
      geom_text(data = df_long %>% filter(Tipo == 'AIC'),
                aes(x = variable, y = station, label = round(value_rel, 4))) +
      scale_fill_gradientn(
        colors = c("#277DF5", "white", "red"),
        name = "Row relative\nAIC",
        breaks = c(0, 0.5, 1),
        labels = c("min", "", "max")
      ) +
      
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(label = title)
  }
  else if (n.metrics == 1){
    df_long <- df %>%
      pivot_longer(cols = -station, names_to = "variable", values_to = "value") %>%
      left_join(df %>% select(station, AIC.mod.sel), by = "station") %>%
      mutate(value_rel = value / AIC.mod.sel) %>%
      group_by(station) %>%
      mutate(value_row_rel = scale_row_values(value_rel)) %>%
      ungroup() %>%
      mutate(
        station = factor(station, levels = unique(station)),
        variable = factor(variable, levels = unique(variable))
      )
    
    ggplot(df_long, aes(x = variable, y = station, fill = value_row_rel)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value_rel, 4))) +
      scale_fill_gradientn(
        colors = c("#277DF5", "white", "red"),
        name = "Row relative\nAIC",
        breaks = c(0, 0.5, 1),
        labels = c("min", "", "max")
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(label = title)
  }
  else{
    stop("'n.metrics not valid. Use 1 or 2")
  }
}

heat.map(df.comp$MHO, n.metrics = 2, title = 'MHO')
heat.map(df.comp$MDO, n.metrics = 2, title = 'MDO')
heat.map(df.comp$MHQ, n.metrics = 1, title = 'MHQ')
heat.map(df.comp$MDQ, n.metrics = 1, title = 'MDQ')

# model MHO not really good for any of the 2 chosen ones
# the rest dont look bad

#----Threshold of MHO models (heat maps of AIC and compare to selected one)----
library(pROC)
df.thresh <- function(stations, common.models, model.type, 
                      thresholds){
  df <- data.frame(matrix(NA, ncol = length(thresholds) * 2 + 1, nrow = length(estaciones)))
  colnames(df) <- c('station', paste0('AIC.', thresholds), paste0('AUC.', thresholds))
  df[, 'station'] <- stations
  
  for (j in 1:length(stations)){
    station <- stations[j]
    cat('Station: ', station, '\n')
    station.p <- paste0(station, '.p')
    data.aux <- common.models[[station]][[model.type]]$data
    formula.aux <- common.models[[station]][[model.type]]$formula
    mod.aux <- glm(formula.aux, data = data.aux, family = binomial(link = 'logit'))
    
    row.aic <- c()
    row.auc <- c()
    for(i in 1:length(thresholds)){
      #AIC
      aux <- update(mod.aux, formula = as.formula(paste0('.~. -poly(', station.p, '.lag, 3) + 
                  I(', station.p, '.lag <', thresholds[i], '):I(log(pmax(', station.p, '.lag, 1e-6))) + 
                  I(', station.p, '.lag >=',  thresholds[i], '):I(', station.p, '.lag)')))
      row.aic <- c(row.aic, aux$aic)
      
      #AUC
      roc <- suppressMessages(
        roc(data.aux$Y, predict(aux, type = 'response'))
      )
      row.auc <- c(row.auc, auc(roc))
    }
    
    df[j, paste0('AIC.', thresholds)] <- row.aic
    df[j, paste0('AUC.', thresholds)] <- row.auc
    
  }
  
  return(df)
}

df.thresh.MHO.pc.3 <- df.thresh(estaciones, common.models, 'MHO.pc.3', 
                                seq(0.5, 5, by = 0.5))
df.thresh.MHO.pc.3$AIC.mod.sel <- df.comp$MHO$AIC.mod.sel
df.thresh.MHO.pc.3$AUC.mod.sel <- df.comp$MHO$AUC.mod.sel
heat.map(df.thresh.MHO.pc.3, n.metrics = 2, 'MHO thresholds intervals')

df.thresh.MHO.pc.3.log <- df.thresh(estaciones, common.models, 'MHO.pc.3', 
                                seq(0.5, 5, by = 0.5))
df.thresh.MHO.pc.3.log$AIC.mod.sel <- df.comp$MHO$AIC.mod.sel
df.thresh.MHO.pc.3.log$AUC.mod.sel <- df.comp$MHO$AUC.mod.sel
heat.map(df.thresh.MHO.pc.3.log, n.metrics = 2, 'MHO thresholds intervals')

df.thresh.MDO.pc.log <- df.thresh(estaciones, common.models, 'MDO.pc', 
                                    seq(7, 13, by = 1))
df.thresh.MDO.pc.log$AIC.mod.common.2 <- df.comp$MDO$AIC.mod.common.2
df.thresh.MDO.pc.log$AUC.mod.common.2 <- df.comp$MDO$AUC.mod.common.2
df.thresh.MDO.pc.log$AIC.mod.common.2.1 <- df.comp$MDO$AIC.mod.common.2.1
df.thresh.MDO.pc.log$AUC.mod.common.2.1 <- df.comp$MDO$AUC.mod.common.2.1
df.thresh.MDO.pc.log$AIC.mod.sel <- df.comp$MDO$AIC.mod.sel
df.thresh.MDO.pc.log$AUC.mod.sel <- df.comp$MDO$AUC.mod.sel
heat.map(df.thresh.MDO.pc.log, n.metrics = 2, 'MDO thresholds intervals')




#comparison of AIC and AUC
mapa_calor <- function(df, tipo = c("AUC", "AIC")) {
  tipo <- match.arg(tipo)  # Forzar que sea uno de los dos
  
  # Pasar a formato largo
  df_long <- melt(df, id.vars = "station")
  
  if (tipo == "AUC") {
    # Escalado global entre 0.5 y 1
    ggplot(df_long, aes(x = variable, y = station, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), size = 5) +
      scale_fill_gradientn(
        colors = c("blue", "lightblue", "orange", "red"),
        limits = c(min(df[, -1]),max(df[, -1])),
        name = "AUC"
      ) +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    
  } else if (tipo == "AIC") {
    # 1️⃣ Relativizar cada fila con respecto a la primera columna
    df_rel <- df
    
    df_long <- melt(df_rel, id.vars = "station")
    
    # 2️⃣ Calcular valores normalizados por fila solo para la escala de color
    df_long <- df_long %>%
      group_by(station) %>%
      mutate(color_val = (value - min(value)) / (max(value) - min(value))) %>%
      ungroup()
    
    # 3️⃣ Graficar usando color_val para el gradiente y value para el texto
    ggplot(df_long, aes(x = variable, y = station, fill = color_val)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 1)), size = 5) +
      scale_fill_gradientn(
        colors = c("#277DF5", "white", "red"),
        name = "Row relative\nAIC",
        breaks = c(0, 0.5, 1),
        labels = c("min", "", "max")
      ) +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
  }
}

thresholds <- seq(0.5, 5, by = 0.5)
AIC.df <- cbind(df.thresh.MHO.pc.3$station,
                df.thresh.MHO.pc.3[, head(grep('AIC', colnames(df.thresh.MHO.pc.3)), -1)],
                df.thresh.MHO.pc.3.log[, head(grep('AIC', colnames(df.thresh.MHO.pc.3.log)), -1)])
colnames(AIC.df) <- c('station', paste0('AIC.lin.', thresholds), paste0('AIC.log.', thresholds))
mapa_calor(AIC.df, 'AIC')
mapa_calor(AIC.df[, c('station', paste0('AIC.log.', thresholds))], 'AIC')

AUC.df <- cbind(df.thresh.MHO.pc.3$station,
                df.thresh.MHO.pc.3[, head(grep('AUC', colnames(df.thresh.MHO.pc.3)), -1)],
                df.thresh.MHO.pc.3.log[, head(grep('AUC', colnames(df.thresh.MHO.pc.3.log)), -1)])
colnames(AUC.df) <- c('station', paste0('AUC.lin.', thresholds), paste0('AUC.log.', thresholds))
mapa_calor(AUC.df, 'AUC')
mapa_calor(AUC.df[, c('station', paste0('AUC.log.', thresholds))], 'AUC')
mapa_calor(df.thresh.MDO.pc.log[, c(1, head(grep('AIC', colnames(df.thresh.MDO.pc.log)), -1))], 'AIC')

#----final common models----
common.models.final <- list()
for (station in estaciones){
  cat('Station: ', station, '\n')
  station.p <- paste0(station, '.p')
  
  #MHO
  data.aux <- common.models[[station]][['MHO.pc.3']]$data
  formula.aux <- common.models[[station]][['MHO.pc.3']]$formula
  mod.aux <- glm(formula.aux, data = data.aux, family = binomial(link = 'logit'))
  
  common.models.final[[station]][['MHO']] <- update(mod.aux, 
              formula = as.formula(paste0('.~. -poly(', station.p, '.lag, 2) + 
              I(', station.p, '.lag <', 2, '):I(log(pmax(', station.p, '.lag, 1e-6))) + 
              I(', station.p, '.lag >=',  2, '):I(', station.p, '.lag)')))
  
  #MDO
  common.models.final[[station]][['MDO']] <- common.models[[station]][['MDO.pc.3']]
  common.models.final[[station]][['MDO.2']] <- common.models[[station]][['MDO.pc.4']]
  data.aux <- common.models[[station]][['MDO.pc.3']]$data
  formula.aux <- common.models[[station]][['MDO.pc.3']]$formula
  mod.aux <- glm(formula.aux, data = data.aux, family = binomial(link = 'logit'))
  common.models.final[[station]][['MDO.3']] <- update(
    mod.aux, 
    formula = as.formula(paste0('.~. -poly(', station.p, '.lag, 2) +
                                I(', station.p, '.lag <', 9, '):I(log(pmax(', station.p, '.lag, 1e-6))) + 
                                I(', station.p, '.lag >=',  9, '):I(', station.p, '.lag)')))
  
  #MHQ
  common.models.final[[station]][['MHQ']] <- common.models[[station]][['MHQ.pc.2']]
  
  #MDQ
  common.models.final[[station]][['MDQ']] <- common.models[[station]][['MDQ.pc']]
 
  
}


#----CORRECCIÓN MHO----
X.MHO <- qread('X.MHO.qs')
for (station in estaciones){
  X <- X.MHO[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = '%Y-%m-%d')
  X <- X %>% filter(date >= per.comun.h[1] & date <= per.comun.h[2])
  
  X <- X %>%
    mutate(
      gr.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] >= 2),
      less.thresh = as.numeric(.data[[paste0(station, '.p.lag')]] < 2)
    )
  
  mod <- common.models.final[[station]]$MHO
  
  common.models.final[[station]]$MHO <- update(
    mod, formula = as.formula(
      paste0(
        '.~. - I(', station, '.p.lag < 2):I(log(pmax(', station, 
        '.p.lag, 1e-06))) - I(', station, '.p.lag >= 2):I(',
        station, '.p.lag) + less.thresh:I(log(pmax(', station, 
        '.p.lag, 1e-06))) + gr.thresh:I(', station, '.p.lag)'
      )
      ),
    data = X
  )
}

#----analisis de comunalidades----
# analyisis of CI of coefficients of models (ONLY COMMON MODELS)
data.common.models <- list()
for (station in estaciones){
  cat('Estación ', station, '\n')
  #data.common.models[[station]][['MHO.pc.vars']] <- common.models[[station]][['MHO.pc']]$coefficients
  # data.common.models[[station]][['MHO.pc.IC']] <- confint.default(common.models[[station]][['MHO.pc']])
  # data.common.models[[station]][['MHO.pc.2.vars']] <- common.models[[station]][['MHO.pc.2']]$coefficients
  # data.common.models[[station]][['MHO.pc.2.IC']] <- confint.default(common.models[[station]][['MHO.pc.2']])
  data.common.models[[station]][['MHO.pc.3.vars']] <- common.models[[station]][['MHO.pc.3']]$coefficients
  data.common.models[[station]][['MHO.pc.3.IC']] <- confint.default(common.models[[station]][['MHO.pc.3']])
  # data.common.models[[station]][['MHO.pc.sel.vars']] <- common.models[[station]][['MHO.pc.sel']]$coefficients
  #data.common.models[[station]][['MHO.pc.sel.IC']] <- confint.default(common.models[[station]][['MHO.pc.sel']])
  
  # data.common.models[[station]][['MDO.pc.vars']] <- common.models[[station]][['MDO.pc']]$coefficients
  # data.common.models[[station]][['MDO.pc.IC']] <- confint.default(common.models[[station]][['MDO.pc']])
  data.common.models[[station]][['MDO.pc.2.vars']] <- common.models[[station]][['MDO.pc.2']]$coefficients
  data.common.models[[station]][['MDO.pc.2.IC']] <- confint.default(common.models[[station]][['MDO.pc.2']])
  data.common.models[[station]][['MDO.pc.3.vars']] <- common.models[[station]][['MDO.pc.3']]$coefficients
  data.common.models[[station]][['MDO.pc.3.IC']] <- confint.default(common.models[[station]][['MDO.pc.3']])
  # data.common.models[[station]][['MDO.pc.sel.vars']] <- common.models[[station]][['MDO.pc.sel']]$coefficients
  #data.common.models[[station]][['MDO.pc.sel.IC']] <- confint.default(common.models[[station]][['MDO.pc.sel']])
  
  # data.common.models[[station]][['MDQ.pc.mu.vars']] <- common.models[[station]][['MDQ.pc']]$mu.coefficients
  # data.common.models[[station]][['MDQ.pc.sigma.vars']] <- common.models[[station]][['MDQ.pc']]$sigma.coefficients
  # data.common.models[[station]][['MDQ.pc.IC']] <- confint(common.models[[station]][['MDQ.pc']])
  # data.common.models[[station]][['MDQ.pc.sel.mu.vars']] <- common.models[[station]][['MDQ.pc.sel']]$mu.coefficients
  # data.common.models[[station]][['MDQ.pc.sel.sigma.vars']] <- common.models[[station]][['MDQ.pc.sel']]$sigma.coefficients
  #data.common.models[[station]][['MDQ.pc.sel.IC']] <- confint(common.models[[station]][['MDQ.pc.sel']])
  
  # data.common.models[[station]][['MHQ.pc.mu.vars']] <- common.models[[station]][['MHQ.pc']]$mu.coefficients
  # data.common.models[[station]][['MHQ.pc.sigma.vars']] <- common.models[[station]][['MHQ.pc']]$sigma.coefficients
  # data.common.models[[station]][['MHQ.pc.IC']] <- confint(common.models[[station]][['MHQ.pc']])
  # data.common.models[[station]][['MHQ.pc.sel.mu.vars']] <- common.models[[station]][['MHQ.pc.sel']]$mu.coefficients
  # data.common.models[[station]][['MHQ.pc.sel.sigma.vars']] <- common.models[[station]][['MHQ.pc.sel']]$sigma.coefficients
  #data.common.models[[station]][['MHQ.pc.sel.IC']] <- confint(common.models[[station]][['MHQ.pc.sel']])
}

qsave(data.common.models, 'data.common.models.qs')
data.common.models <- qread('data.common.models.qs')

plots.CI <- function(data, df.comp, model.type, 
                     metric.submodel, stations, 
                     filter = FALSE, quantity = FALSE){
  if (quantity == FALSE){
    if(filter == TRUE){
      out <- which.max(df.comp[[model.type]][[metric.submodel]])
      
      stations.names <- stations$STAID
      
      aux <- stations.names[-out]
      df <- data.frame(
        station = factor(aux, levels = aux),
        Z = stations[which(stations$STAID %in% aux), 'Z']
      )
      
      df$station <- factor(df$station, levels = df$station[order(df$Z)])
    }else{
      stations.names <- stations$STAID
      
      aux <- stations.names
      df <- data.frame(
        station = factor(aux, levels = aux),
        Z = stations[which(stations$STAID %in% aux), 'Z']
      )
      
      df$station <- factor(df$station, levels = df$station[order(df$Z)])
    }
    
    #number of vars
    L <- length(data[[aux[1]]][[paste0(model.type, '.pc.vars')]])
    
    for (j in 1:L){
      for (i in 1:length(aux)){
        station <- aux[i]
        
        vars <- names(data[[station]][[paste0(model.type, '.pc.vars')]])
        var <- vars[j]
        
        df[i, 'value'] <- data[[station]][[paste0(model.type, '.pc.vars')]][var]
        
        df[i, c("low", "high")] <- data[[station]][[paste0(model.type, '.pc.IC')]][var, ]
        
      }
      
      #print(dim(df))
      # cat(var, '\n')
      max <- max(df[['low']])
      min <- min(df[['high']])
      # cat(min, '-', max, '\n')
      #cat(sum(df[['low']] <= max) == nrow(df) & sum(df[['high']] >= min) == nrow(df), '\n')
      
      if (grepl('.lag', var) == T){
        var.name <- gsub(aux[length(aux)], 'station', var)
      }
      else if (grepl('.day', var) == T){
        var.name <- gsub(aux[length(aux)], 'station', var)
      }
      else{
        var.name <- var
      }
      
      g <- ggplot(df, aes(x = station, y = value))
      
      if (sum(df[['low']] <= max) == nrow(df) & sum(df[['high']] >= min) == nrow(df)){
        g <- g + geom_rect(aes(xmin = -Inf, xmax = Inf, 
                               ymin = min(c(min, max)) - 0.01, 
                               ymax = max(c(min, max)) + 0.01),
                       fill = "grey", alpha = 0.2)
      }
      
      g <- g + 
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.5) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
        # geom_hline(yintercept = min, linetype = 'dashed', color = 'blue') +
        # geom_hline(yintercept = max, linetype = 'dashed', color = 'red') +
        theme_minimal(base_size = 14) +
        labs(
          x = "Station",
          y = "Estimated coefficient",
          title = paste("CI of MDO coefficient:", var.name)
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      
      #sum(df[['low']] < max) == length(aux) & sum(df[['high']] > min) == length(aux)
      print(g)
    }
  }
  else{
    if(filter == TRUE){
      out <- which.max(df.comp[[model.type]][[metric.submodel]])
      
      stations.names <- stations$STAID
      
      aux <- stations.names[-out]
      df <- data.frame(
        station = factor(aux, levels = aux),
        Z = stations[which(stations$STAID %in% aux), 'Z']
      )
      
      df$station <- factor(df$station, levels = df$station[order(df$Z)])
    }
    else{
      stations.names <- stations$STAID
      
      aux <- stations.names
      df <- data.frame(
        station = factor(aux, levels = aux),
        Z = stations[which(stations$STAID %in% aux), 'Z']
      )
      
      df$station <- factor(df$station, levels = df$station[order(df$Z)])
    }
    
    #number of mu vars
    L <- length(data[[aux[1]]][[paste0(model.type, '.pc.mu.vars')]])
    
    for (j in 1:L){
      for (i in 1:length(aux)){
        station <- aux[i]
        
        vars <- names(data[[station]][[paste0(model.type, '.pc.mu.vars')]])
        var <- vars[j]
        
        df[i, 'value'] <- data[[station]][[paste0(model.type, '.pc.mu.vars')]][var]
        
        df[i, c("low", "high")] <- data[[station]][[paste0(model.type, '.pc.IC')]][paste0('mu.', var), ]
        
        
      }
      
      max <- max(df[['low']])
      min <- min(df[['high']])
      
      if (grepl('.lag', var) == T){
        var.name <- gsub(aux[length(aux)], 'station', var)
      }
      else if (grepl('.day', var) == T){
        var.name <- gsub(aux[length(aux)], 'station', var)
      }
      else{
        var.name <- var
      }
      
      g <- ggplot(df, aes(x = station, y = value))
      
      if (sum(df[['low']] <= max) == nrow(df) & sum(df[['high']] >= min) == nrow(df)){
        g <- g + geom_rect(aes(xmin = -Inf, xmax = Inf, 
                               ymin = min(c(min, max)) - 0.01, 
                               ymax = max(c(min, max)) + 0.01),
                           fill = "grey", alpha = 0.2)
      }
      
      g <- g +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.5) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
        theme_minimal(base_size = 14) +
        labs(
          x = "Station",
          y = "Estimated coefficient",
          title = paste("CI of", model.type, "mu.coefficient:", var.name)
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(g)
    }
    
    #length of sigma vars
    L <- length(data[[aux[1]]][[paste0(model.type, '.pc.sigma.vars')]])
    
    for (j in 1:L){
      for (i in 1:length(aux)){
        station <- aux[i]
        
        vars <- names(data[[station]][[paste0(model.type, '.pc.sigma.vars')]])
        var <- vars[j]
        
        df[i, 'value'] <- data[[station]][[paste0(model.type, '.pc.sigma.vars')]][var]
        
        df[i, c("low", "high")] <- data[[station]][[paste0(model.type, '.pc.IC')]][paste0('sigma.', var), ]
        
      }
      
      max <- max(df[['low']])
      min <- min(df[['high']])
      
      if (grepl('.lag', var) == T){
        var.name <- gsub(aux[length(aux)], 'station', var)
      }
      else if (grepl('.day', var) == T){
        var.name <- gsub(aux[length(aux)], 'station', var)
      }
      else{
        var.name <- var
      }
      
      g <- ggplot(df, aes(x = station, y = value))
      
      if (sum(df[['low']] <= max) == nrow(df) & sum(df[['high']] >= min) == nrow(df)){
        g <- g + geom_rect(aes(xmin = -Inf, xmax = Inf, 
                               ymin = min(c(min, max)) - 0.01, 
                               ymax = max(c(min, max)) + 0.01),
                           fill = "grey", alpha = 0.2)
      }
      
      g <- g +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.5) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
        theme_minimal(base_size = 14) +
        labs(
          x = "Station",
          y = "Estimated coefficient",
          title = paste("CI of", model.type, "sigma.coefficient:", var.name)
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(g)
    }
  }
  
}

plots.CI(data.common.models, df.comp, 'MDO', 'AIC.mod.common', stations, filter = T)
plots.CI(data.common.models, df.comp, 'MHO', 'AIC.mod.common', stations, filter = T)
plots.CI(data.common.models, df.comp, 'MHQ', 'AIC.mod.common', stations, quantity = T)
plots.CI(data.common.models, df.comp, 'MDQ', 'AIC.mod.common', stations, quantity = T)

#MDO
# var1
aux <- estaciones[-which(estaciones == 'R037')]
df <- data.frame(
  station = factor(aux, levels = aux),
  Z = stations[which(stations$STAID %in% aux), 'Z']
)

df$station <- factor(df$station, levels = df$station[order(df$Z)])

#quito la mala

for (i in 1:length(aux)){
  station <- aux[i]
  
  vars <- names(data.common.models[[station]][['MDO.pc.vars']])
  var <- vars[19]
  
  df[i, 'value'] <- data.common.models[[station]][['MDO.pc.vars']][var]
  
  df[i, c("low", "high")] <- data.common.models[[station]][['MDO.pc.IC']][var, ]
  
}

ggplot(df, aes(x = station, y = value)) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.5) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme_minimal(base_size = 14) +
  labs(
    x = "Estación",
    y = "Coeficiente estimado",
    title = paste("Intervalos de confianza del coeficiente MDO:", var)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#----AUC separados----

df <- data.frame(
  station = estaciones
)
library(pROC)
for (i in 1:length(estaciones)){
  station <- estaciones[i]
  
  X <- X.MDO[[station]]
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                    format = '%Y-%m-%d')
  
  ind <- which(X$date >= per.comun.day[1] & X$date <= per.comun.day[2])
  
  mod.common <- common.models.final[[station]][['MDO']]
  
  roc.mod.common <- suppressMessages(
    roc(X$Y[ind], predict(mod.common, type = 'response'))
  )
  
  df[i, 'AUC.mod.common'] <- auc(roc.mod.common)
    
}
summary(df$AUC.mod.common)  
  