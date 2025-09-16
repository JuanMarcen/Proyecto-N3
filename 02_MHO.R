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
mod <- mho_list[[estaciones[1]]]
X <- X_list[[estaciones[1]]]
mod_rebuilt <- glm(formula = formula(mod), family = binomial(link = "logit"), data = X)
step(mod_rebuilt, direction = "backward")

# actualización de modelos
basura <- update(mod, data = X, formula = .~. + poly(EM71.p.day, 3))
basura <- update(mho.EM71, data = X.EM71, formula = .~. + 
                   poly(EM71.p.day, 2) + poly(EM71.p.lag, 2))

library(gam)
aux.marcador <- is.element(X$mes, c(6,7,8)) #meses de mayor lluvia al parecer
plot(gam(formula = as.formula('Y ~ s(EM71.p.day)'), data = X[aux.marcador, ]))
plot(gam(formula = as.formula('Y ~ s(EM71.p.lag)'), data = X[aux.marcador, ]))

# selección variables (automático) según BIC
mod_null <- glm(Y ~ 1, family = binomial(logit), data = X)
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
  # harmonics: lista de armónicos
  
  
  
  mod.aux <- initial_model
  cat('Initial model ', deparse(formula(mod.aux)), '\n')
  cat('BIC: ', BIC(mod.aux), '\n')
  
  if (!is.null(vars)) {
    for (var in vars){
      formula.aux <- update(formula(mod.aux), paste(". ~ . +", var))
      mod.temp <- glm(formula.aux, data = data, family = binomial(link = "logit"))
      
      if (BIC(mod.temp) < BIC(mod.aux)){
        cat('Added: ', var, '\n')
        cat('BIC: ', BIC(mod.temp), '\n')
        mod.aux <- mod.temp
      }
    }
    
  }
  
  
  for (h in harmonics.l){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- glm(formula.aux, data = data, family = binomial(logit))
    
    if (BIC(mod.temp) < BIC(mod.aux)){
      cat('Added: ',paste(h, collapse = '+'), '\n')
      cat('BIC: ', BIC(mod.temp), '\n')
      
      mod.aux <- mod.temp
    }else{
      break
    }
  }
  
  for (h in harmonics.h){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- glm(formula.aux, data = data, family = binomial(logit))
    
    if (BIC(mod.temp) < BIC(mod.aux)){
      cat('Added: ',paste(h, collapse = '+'), '\n')
      cat('BIC: ', BIC(mod.temp), '\n')
      
      mod.aux <- mod.temp
    }else{
      break
    }
  }
  
  cat('\nFinal model: ', deparse(formula(mod.aux)), '\n')
  cat('BIC: ', BIC(mod.aux), '\n')
  
  return(mod.aux)
}

basura <- step_rlog(mod_null, data = X, c('EM71.p.day', 'EM71.p.lag'), harmonics.l, harmonics.h)

step(mod, direction = 'backward')

# evaluation #put in a Rmarkdow
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


