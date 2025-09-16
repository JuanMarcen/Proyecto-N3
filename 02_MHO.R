# model hourly ocurrence
# substract days in which it doesn't rain
rm(list = ls())

load('data.RData')


# harmonics
cs <- function(t,harmonics=1) {
  # if(min(t) <0 | max(t) > 1){ stop(" t must be in [0,1] range")}
  if(min(harmonics) <1){stop("harmonics > = 1")}
  ret <- numeric(0)
  for ( i in harmonics) {
    ret <- cbind( ret, cos(2*pi*i*t/365), sin(2*pi*i*t/365))
  }
  if (missing(harmonics)) cnames <- c('c','s')
  else {
    cnames <- paste( c("c","s"), rep(harmonics, each = 2),sep=".")
  }
  colnames(ret) <- cnames
  ret
}

l <- 1:365
h <- 0:23
harm_l <- data.frame(l = l,
                     cs(l, harmonics = 1:4))
colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
harm_h <- data.frame(h = h,
                        cs(0:23, harmonics = 1:4))
colnames(harm_h)[2:ncol(harm_h)] <- paste0(colnames(harm_h)[2:ncol(harm_h)], '.h')


library(dplyr)
df_hours <- df_hours %>%
  left_join(harm_l, by = 'l') %>%
  left_join(harm_h, by = 'h')

# creation of final X matrix -- Generalize
# for all stations
# Saving design matrix (X) and model 
mho_list <- list()
X_list <- list()
for (station in estaciones){
  station.p <- paste0(station, '.p')
  X <- df_hours[, c('t', 'l', 'mes', 'dia.mes', 'h', station.p,
                            colnames(harm_l)[2:ncol(harm_l)],
                            colnames(harm_h)[2:ncol(harm_h)])]
  
  # whole day
  p_day <- df_days[,  c('t', 'l', 'mes', 'dia.mes', station.p)]
  colnames(p_day)[colnames(p_day) == station.p] <- paste0(station.p, '.day')
  
  X <- X %>%
    left_join(p_day, by = c('t', 'l', 'mes', 'dia.mes'))
  
  # lags
  X <- X %>% 
    mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]]))
  
  #eliminate days with no rain
  X_final <- X[-which(X[[paste0(station.p,'.day')]] == 0), ]
  
  # hourly indicator of rain
  X_final <- X_final %>%
    mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0)) %>% 
    relocate(Y, .after = !!station.p) %>%
    as.data.frame() %>% na.omit()
  
  #na omit necesario?
  
  # MHO
  X_list[[station]] <- X_final
  formula <- as.formula(
    paste('Y ~', paste(colnames(X_list[[station]])[8:ncol(X_list[[station]])], collapse = '+'))
  )
  
  mho_list[[station]] <- glm(formula = formula, family = binomial(logit), data = X_list[[station]])
  
  
  # assign(paste0('X.', station), X_final)
  # assign(paste0('mho.', station), mho)
  
  # rm(list = c('X', 'X_final', 'mho', 'formula', 'p_day', 'station.p', 'station'))
}

# estudio de covariables en los modelos
## ejemplo para EM71
# mod <- get(paste0('mho.', estaciones[1]))
# X <- get(paste0('X.', estaciones[1]))
mod <- mho_list[[estaciones[15]]]
X <- X_list[[estaciones[15]]]
mod_rebuilt <- glm(formula = formula(mod), family = binomial(link = "logit"), data = X)
step(mod_rebuilt, direction = "backward")

#----actualización de modelos----
#borrador
basura <- update(mod, data = X, formula = .~. + poly(P024.p.lag, 3))
basura <- update(mho.EM71, data = X.EM71, formula = .~. + 
                   poly(EM71.p.day, 2) + poly(EM71.p.lag, 2))

library(gam)
aux.marcador <- is.element(X$mes, c(6,7,8)) #meses de mayor lluvia al parecer
mod1 <- gam(formula = as.formula('Y ~ s(EM71.p.day)'), data = X[aux.marcador, ])
summary(mod1)
plot(gam(formula = as.formula('Y ~ s(EM71.p.day)'), data = X[aux.marcador, ]))
plot(gam(formula = as.formula('Y ~ s(EM71.p.lag, df = 3)'), data = X[aux.marcador, ]))

#decisión del grado de polinomio (a mano)
for (station in estaciones){
  X <- X_list[[station]]
  aux.marcador <- is.element(X$mes, c(6,7,8)) #meses de mayor lluvia al parecer
  #variable lag
  mod.lag <- gam(formula = as.formula(paste0('Y ~ s(', station, '.p.lag)')), data = X[aux.marcador, ])
  plot(mod.lag, main = paste(station))
  #variable del dia
  mod.day <- gam(formula = as.formula(paste0('Y ~ s(', station, '.p.day)')), data = X[aux.marcador, ])
  plot(mod.day, main = paste(station))
}

degrees.p.lag <- c(3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(degrees.p.lag) <- estaciones
degrees.p.day <- c(3, 1, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(degrees.p.day) <- estaciones

# selección variables (automático) según AIC
#mod_null <- glm(Y ~ 1, family = binomial(logit), data = X)
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
                      harmonics.h){ 
  # argumentos:
  # initial_model: modelo inicial (modelo nulo)
  # data: datos a utilizar (contiene explicada y explicativas)
  # vars: covariables no armónicas que se desean introducir en el modelo (char)
  # incluidas siempre? pues al modelo initial 
  # harmonics: lista de armónicos
  
  mod.aux <- initial_model
  cat('Initial model ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n')
  
  #  asumimos que las covariables siempre son escogidas
  # if (!is.null(vars)) {
  #   for (var in vars){
  #     formula.aux <- update(formula(mod.aux), paste(". ~ . +", var))
  #     mod.temp <- glm(formula.aux, data = data, family = binomial(link = "logit"))
  #     
  #     if (AIC(mod.temp) < AIC(mod.aux)){
  #       cat('Added: ', var, '\n')
  #       cat('AIC: ', AIC(mod.temp), '\n')
  #       mod.aux <- mod.temp
  #     }
  #   }
  #   
  # }
  
  
  for (h in harmonics.l){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- glm(formula.aux, data = data, family = binomial(logit))
    
    if (AIC(mod.temp) < AIC(mod.aux)){
      cat('Added: ',paste(h, collapse = '+'), '\n')
      cat('AIC: ', AIC(mod.temp), '\n')
      
      mod.aux <- mod.temp
    }else{
      break
    }
  }
  
  for (h in harmonics.h){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- glm(formula.aux, data = data, family = binomial(logit))
    
    if (AIC(mod.temp) < AIC(mod.aux)){
      cat('Added: ',paste(h, collapse = '+'), '\n')
      cat('AIC: ', AIC(mod.temp), '\n')
      
      mod.aux <- mod.temp
    }else{
      break
    }
  }
  
  cat('\nFinal model: ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n\n')
  
  return(mod.aux)
}

#actualización de modelos
for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  deg.lag <- degrees.p.lag[station]
  deg.day <- degrees.p.day[station]
  
  formula_null <- as.formula(paste0('Y ~ ', paste0('poly(', station, '.p.day, ',deg.day, ')'), '+',
                                    paste0('poly(', station, '.p.lag, ',deg.lag, ')')))
  
  mod_null <- glm(formula_null, family = binomial(logit), data = X_list[[station]])
  
  mho_list[[station]] <- step_rlog(mod_null, 
                                   data = X_list[[station]], 
                                   vars = c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                            paste0('poly(', station, '.p.lag, ',deg.lag, ')')), 
                                   harmonics.l, harmonics.h)
}

# for (station in estaciones){
#   print(summary(mho_list[[station]]))
# }

# guardado de una lista definitiva por estacion
# estacion --> modelo, df, variables escogidas

MHO <- list()
for (station in estaciones){
  MHO[[station]][['mho']] <- mho_list[[station]]
  MHO[[station]][['vars']] <- mho_list[[station]]$coefficients
  MHO[[station]][['X']] <- X_list[[station]]
} 

#----Estudio comunalidades----
# stand by 

#----evaluation #put in a Rmarkdow----
library(pROC)

for (station in estaciones){
  mho <- get(paste0('mho.', station))
  X <- get(paste0('X.', station))
  
  summary(mho)
  
  boxplot(mho$fitted.values ~ X$Y,
          xlab = 'Y real',
          ylab = 'P(Y = 1| X)',
          main = station)
  
  pred <- predict(mho, type = 'response')
  roc <- roc(X$Y, pred)
  plot(roc, col="blue", main = paste('Curva ROC', station), 
       print.thres = T,
       print.auc = T)
    
  thres <- coords(roc, 'best')[1]
  
  pred_class <- ifelse(pred > thres, 1, 0)
  table <- table(Predicho = cut(pred,c(0,thres,1)), Real = X$Y)
  
  print(table)
  
  print(round(100*prop.table(table, 2)), 2)
  
  rm(list = c('X', 'mho', 'pred', 'rho', 'station', 'thres', 'pred_class', 'table'))
}


