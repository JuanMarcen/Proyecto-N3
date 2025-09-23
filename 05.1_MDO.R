# MODELS OF DAIRLY DATA --> OCURRENCE

rm(list = ls())

load('data.RData')

library(lubridate)
# In the covariates matrix we add the climate variables
global_df <- readRDS('global_df.rds')
global_df$t <- year(global_df$date)
global_df <- global_df[, -which(colnames(global_df) %in% c('zg300.', 'zg500.', 'zg700.', 
                                                           'zt300.', 'zt500.', 'zt700.'))]

#example 1 station
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

l <- 1:365
harm_l <- data.frame(l = l,
                     cs(l, harmonics = 1:4, 365))
colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')

library(dplyr)
df_days <- df_days %>%
  left_join(harm_l, by = 'l') 

#----Data building----
X_list <- list()
for (station in estaciones){
  station.p <- paste0(station, '.p')
  
  X <- df_days[, c('t', 'l', 'mes', 'dia.mes', station.p,
                   colnames(harm_l)[2:ncol(harm_l)])]
  
  # lags
  X <- X %>% 
    mutate(!!paste0(station.p, '.lag') := lag(.data[[station.p]])) %>%
    na.omit()
  
  X <- X %>%
    left_join(global_df[global_df$STAID == station, c(3,4,8:ncol(global_df))],
              by =c('t','l'))
  
  X_list[[station]] <- X %>% 
    mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0)) %>%
    relocate(Y, .after = !!station.p) %>%
    as.data.frame() %>%
    na.omit()
}

#----Variable selection----
harmonics.l <- list(
  h1 = c('s.1.l', 'c.1.l'),
  h2 = c('s.2.l', 'c.2.l'),
  h3 = c('s.3.l', 'c.3.l'),
  h4 = c('s.4.l', 'c.4.l')
)

step_rlog <- function(initial_model,
                      data,
                      vars,
                      harmonics.l){ 
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
  if (!is.null(vars)) {
    for (var in vars){
      formula.aux <- update(formula(mod.aux), paste(". ~ . +", var))
      mod.temp <- glm(formula.aux, data = data, family = binomial(link = "logit"))
      
      if (AIC(mod.temp) < AIC(mod.aux)){
        cat('Added: ', var, '\n')
        cat('AIC: ', AIC(mod.temp), '\n')
        mod.aux <- mod.temp
      }
    }
    
  }
  
  
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
  
  
  cat('\nFinal model: ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n\n')
  
  return(mod.aux)
}

mdo_list <- list()
for (station in estaciones){
  cat('Estación: ',station, '\n\n')
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = X_list[[station]])
  
  mdo_list[[station]] <- step_rlog(mod_null,
                                   X_list[[station]],
                                   colnames(X_list[[station]])[15:ncol(X)],
                                   harmonics.l)
  
  
}

#guardado en otro data frame
MDO <- list()
for (station in estaciones){
  MDO[[station]][['mdo']] <- mdo_list[[station]]
  MDO[[station]][['vars']] <- mdo_list[[station]]$coefficients
  MDO[[station]][['X']] <- X_list[[station]]
} 

saveRDS(MDO, 'MDO.rds')

library(gam)
plot(gam(formula = as.formula(paste('Y ~ s(', colnames(X_list[[station]])[36], ')')), data = X))

#----Comunalidades---
stations <- readRDS('stations.rds')
stations.no.ateca <- stations[-nrow(stations), ]

vars <- list()
for (station in estaciones){
  vars[['todas']][[station]] <- names(MDO[[station]]$vars)
}

for (grupo in unique(stations.no.ateca$color)){
  
  aux <- stations[stations$color == grupo, 'STAID']
  
  for (station in aux){
    vars[[grupo]][[station]] <- names(MDO[[station]]$vars)
  }
  
  rm('aux')
}

#todas
common <- Reduce(intersect, vars[['todas']])
common.medio.ebro <- Reduce(intersect, vars[['blue']])
common.alto.jalon <- Reduce(intersect, vars[['forestgreen']])
common.bajo.jalon <- Reduce(intersect, vars[['red']])

#df de valores de los coeficientes
df.common <- function(common.vars, group = NULL){
  df <- data.frame(matrix(NA, nrow = length(common.vars)))
  colnames(df) <- 'var'
  df$var <- common.vars
  
  if(!is.null(group)){
    for (station in stations$STAID[stations$color == group]){
      df[[station]] <- MDO[[station]]$vars[common.vars]
    }
  }else{
    for (station in estaciones){
      df[[station]] <- MDO[[station]]$vars[common.vars]
    }
  }
  
  return(df)
}

df.common.todas <- df.common(common, group = NULL)
df.common.medio.ebro <- df.common(common.medio.ebro, 'blue')
df.common.alto.jalon <- df.common(common.alto.jalon, 'forestgreen')
df.common.bajo.jalon <- df.common(common.bajo.jalon, 'red')

#not very useful
# for (var in common){
#   aux <- t(df.common[df.common$var == var, -1])
#   plot(aux, type = 'b', 
#        xaxt = 'n', 
#        xlab = '', 
#        ylab = 'Value', 
#        main = paste('Valor coeficiente', var))
#   axis(1, at = 1:length(estaciones), labels = estaciones, las = 2)
# }
# rownames(aux)

#distance matrix
require(factoextra)
dist.plots <- function(df, common.vars){
  for (var in common.vars){
    dist <- dist(t(df[df[['var']] == var, -1]), method = 'euclidean')
    g <- fviz_dist(dist, show_labels = TRUE, order = TRUE) +
      ggtitle(paste('Diferencia del valor de', var))
    print(g)
  }
}

dist.plots(df.common.todas, common)
dist.plots(df.common.medio.ebro, common.medio.ebro)
dist.plots(df.common.alto.jalon, common.alto.jalon)
dist.plots(df.common.bajo.jalon, common.bajo.jalon)

#numerical summaries
num.sum <- function(df, common.vars){
  for (var in common.vars){
    cat('Resumen númerico ', var, ':\n', 
        t(summary(t(df[df[['var']] == var, -1]))), '\n\n')
  }
  
}

num.sum(df.common.todas, common)
num.sum(df.common.medio.ebro, common.medio.ebro)
num.sum(df.common.alto.jalon, common.alto.jalon)
num.sum(df.common.bajo.jalon, common.bajo.jalon)


# clustering entre coeficientes, sin tener en cuanta laz zonas predefinidas


#----Evalutation of goodness of fit----
library(pROC)

for (station in estaciones){
  mdo <- mdo_list[[station]]
  X <- X_list[[station]]
  
  summary(mdo)
  
  boxplot(mdo$fitted.values ~ X$Y,
          xlab = 'Y real',
          ylab = 'P(Y = 1| X)',
          main = station)
  
  pred <- predict(mdo, type = 'response')
  roc <- roc(X$Y, pred)
  plot(roc, col="blue", main = paste('Curva ROC', station), 
       print.thres = T,
       print.auc = T)
  
  thres <- coords(roc, 'best')[1]
  
  pred_class <- ifelse(pred > thres, 1, 0)
  table <- table(Predicho = cut(pred,c(0,thres,1)), Real = X$Y)
  
  print(table)
  
  print(round(100*prop.table(table, 2)), 2)
}

