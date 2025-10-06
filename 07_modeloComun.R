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

mo <- 'M5'
df.mho <- plot.auc.n.vars('M5', estaciones, MHO)
# MHO escogido == 19
df.mdo <- plot.auc.n.vars('M5', estaciones, MDO)
# MDO escogido == 11

# MODELOS DE CANTIDAD
n.vars.mdq <- c()
n.vars.mhq <- c()
mq <- 'M6'
for (station in estaciones){
  n.vars.mdq <- c(n.vars.mdq, length(MDQ[[station]][[paste0('vars.', mq)]]) -1)
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
mod.comun <- function(station, per.comun, data, mod, chosen, type = 'log', subtype = NULL){
  X <- data[[station]]$X
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), 
                    format = '%Y-%m-%d')
  X <- X %>%
    filter(date >= per.comun[1] & date <= per.comun[2])
  
  if (type == 'log'){
    vars <- labels(terms(data[[estaciones[chosen]]][[mod]]$formula))
    vars <- gsub(estaciones[chosen], station, vars)
    formula <- as.formula(paste('Y ~', 
                                paste(vars, collapse = '+')))
    
    
    #cat('Fórmula del modelo común: \n', deparse(formula), '\n')
    
    mod.per.comun <- glm(formula, family = binomial(logit), data = X)
  }else if (type == 'gamma'){
    
    vars <- labels(terms(data[[estaciones[chosen]]][[mod]]$mu.formula))
    vars <- gsub(estaciones[chosen], station, vars)
    formula <- as.formula(paste(paste0(station, '.p'), '~', 
                                paste(vars, collapse = '+')))
    
    if (subtype == 'day'){
      
      mod.per.comun <- gamlss(formula, 
                              sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)), 
                              family = GA, 
                              data = X,
                              trace = F)
      
    }else if(subtype == 'hour'){
      mod.per.comun <- gamlss(formula, 
                              sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                                I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
                              family = GA, 
                              data = X,
                              trace = F)
    }else{
      stop("Subtype not valid. Use 'day' or 'hour'.")
    }
    
  }else{
    stop("Type not valid. Use 'log' or 'gamma'.")
  }
  
  
  return(mod.per.comun)
}

common.models <- list()
for(station in estaciones){
  cat('Estación ', station, '\n')
  common.models[[station]][['MHO.pc']] <- mod.comun(station, per.comun.h, MHO, 'M5', 19, type = 'log')
  common.models[[station]][['MDO.pc']] <- mod.comun(station, per.comun.day, MDO, 'M5', 11, type = 'log')
  common.models[[station]][['MHQ.pc']] <- mod.comun(station, per.comun.h, MHQ, 'M6', 15, type = 'gamma', subtype = 'hour')
  common.models[[station]][['MDQ.pc']] <- mod.comun(station, per.comun.day, MDQ, 'M6', 8, type = 'gamma', subtype = 'day')
}

# ajuste de modelos por selección AIC en el periodo común
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

qsave(common.models, 'common.models.qs', preset = 'archive', compress_level = 22)

common.models <- qread('common.models.qs')


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

#----guardado de coss solo necesarias----
df.comp <- list(
  MHO = df.MHO,
  MDO = df.MDO,
  MHQ = df.MHQ,
  MDQ = df.MDQ 
)

qsave(df.comp, 'df.comp.qs')

data.common.models <- list()
for (station in estaciones){
  cat('Estación ', station, '\n')
  data.common.models[[station]][['MHO.pc.vars']] <- common.models[[station]][['MHO.pc']]$coefficients
  data.common.models[[station]][['MHO.pc.IC']] <- confint.default(common.models[[station]][['MHO.pc']])
  data.common.models[[station]][['MHO.pc.sel.vars']] <- common.models[[station]][['MHO.pc.sel']]$coefficients
  data.common.models[[station]][['MHO.pc.sel.IC']] <- confint.default(common.models[[station]][['MHO.pc.sel']])
  
  data.common.models[[station]][['MDO.pc.vars']] <- common.models[[station]][['MDO.pc']]$coefficients
  data.common.models[[station]][['MDO.pc.IC']] <- confint.default(common.models[[station]][['MDO.pc']])
  data.common.models[[station]][['MDO.pc.sel.vars']] <- common.models[[station]][['MDO.pc.sel']]$coefficients
  data.common.models[[station]][['MDO.pc.sel.IC']] <- confint.default(common.models[[station]][['MDO.pc.sel']])
  
  data.common.models[[station]][['MDQ.pc.mu.vars']] <- common.models[[station]][['MDQ.pc']]$mu.coefficients
  data.common.models[[station]][['MDQ.pc.sigma.vars']] <- common.models[[station]][['MDQ.pc']]$sigma.coefficients
  data.common.models[[station]][['MDQ.pc.IC']] <- confint(common.models[[station]][['MDQ.pc']])
  data.common.models[[station]][['MDQ.pc.sel.mu.vars']] <- common.models[[station]][['MDQ.pc.sel']]$mu.coefficients
  data.common.models[[station]][['MDQ.pc.sel.sigma.vars']] <- common.models[[station]][['MDQ.pc.sel']]$sigma.coefficients
  data.common.models[[station]][['MDQ.pc.sel.IC']] <- confint(common.models[[station]][['MDQ.pc.sel']])
  
  data.common.models[[station]][['MHQ.pc.mu.vars']] <- common.models[[station]][['MHQ.pc']]$mu.coefficients
  data.common.models[[station]][['MHQ.pc.sigma.vars']] <- common.models[[station]][['MHQ.pc']]$sigma.coefficients
  data.common.models[[station]][['MHQ.pc.IC']] <- confint(common.models[[station]][['MHQ.pc']])
  data.common.models[[station]][['MHQ.pc.sel.mu.vars']] <- common.models[[station]][['MHQ.pc.sel']]$mu.coefficients
  data.common.models[[station]][['MHQ.pc.sel.sigma.vars']] <- common.models[[station]][['MHQ.pc.sel']]$sigma.coefficients
  data.common.models[[station]][['MHQ.pc.sel.IC']] <- confint(common.models[[station]][['MHQ.pc.sel']])
}


qsave(data.common.models, 'data.common.models.qs')
