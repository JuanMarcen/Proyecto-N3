# MODELS OF DAIRLY DATA --> OCURRENCE

rm(list = ls())

load('data.RData')

library(lubridate)
library(dplyr)
# In the covariates matrix we add the climate variables
global_df <- readRDS('global_df.rds')
global_df$t <- year(global_df$date)
global_df <- global_df[, -which(colnames(global_df) %in% c('zg300.', 'zg500.', 'zg700.', 
                                                           'zt300.', 'zt500.', 'zt700.'))]
global_df <- global_df %>%
  mutate(across(8:ncol(.), ~ lag(.), .names = "{.col}.lag"))

# cálculo de periodo comun
n <- length(estaciones)
aux <- apply(df_days[, paste0(estaciones, '.p')], 1, function(row) all(!is.na(row)))

first <- which(aux)[1]
last <- tail(which(aux), 1)
date.first <- as.Date(paste0(df_days[first, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
date.last <- as.Date(paste0(df_days[last, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
#datos ERA5 hasta 2023 creo
aux2 <- apply(global_df[, (8:ncol(global_df))], 1, function(row) all(!is.na(row)))
last <- tail(which(aux2), 1)
date.last2 <- global_df$date[last]


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
M1_list <- list()
M2_list <- list()
M3_list <- list()
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
  
  formula_M1 <- as.formula(paste('Y ~', paste(c(colnames(harm_l)[2:ncol(harm_l)]),
                                              collapse = '+')))
  
  ind.all <-  grep("z", names(X_list[[station]]))
  ind.lag <- grep("^z.*\\.lag$", names(X_list[[station]]))
  ind.no.lag <- setdiff(ind.all, ind.lag)
  p.lag <- which(names(X_list[[station]]) == paste0(station.p, '.lag'))
  formula_M2 <- as.formula(paste('Y ~', paste(
    colnames(X_list[[station]])[c(p.lag, ind.no.lag)], collapse = '+'
  )))

  
  formula_M3 <- as.formula(
    paste('Y ~', paste(
      colnames(X_list[[station]])[c(7:15, ind.no.lag)], collapse = '+')
  ))

  cat('Ajuste modelo M1: ', deparse(formula_M1), '\n')
  M1_list[[station]] <- glm(formula = formula_M1, family = binomial(logit), data = X_list[[station]])

  cat('Ajuste modelo M2: ', deparse(formula_M2), '\n')
  M2_list[[station]] <- glm(formula = formula_M2, family = binomial(logit), data = X_list[[station]])

  cat('Ajuste modelo M3: ', deparse(formula_M3), '\n')
  M3_list[[station]] <- glm(formula = formula_M3, family = binomial(logit), data = X_list[[station]])
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
  # CAMBIAR A QUE HAGA STEP 
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
  scope.aux <- update(formula(mod.aux), paste('. ~ . +', paste(vars, collapse = '+')))
  
  mod.aux <- step(mod.aux, scope = scope.aux, direction = 'both', trace = FALSE)
  
  cat('Model after ', 'AIC ', 'step algorithm for variables: ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n')
  
  
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


# YA LO HEMOS HECHO Y GUARDADO CON LO ANTIGUO, NO ES NECESARIO CAMBIAR NADA,
# SI QUISIERAMOS VOLVER A LANZAR HASTA EL MODELO m8, HABRÍA QUE FILTRAR SIN LOS LAGS
M4_list <- list()
for (station in estaciones[1]){
  cat('Estación: ',station, '\n\n')
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = X_list[[station]])
  
  M4_list[[station]] <- step_rlog(mod_null,
                                   X_list[[station]],
                                   colnames(X_list[[station]])[15:ncol(X)],
                                   harmonics.l)
  
  
}

#guardado en otro data frame
MDO <- list()
for (station in estaciones){
  MDO[[station]][['M1']] <- M1_list[[station]]
  MDO[[station]][['vars.M1']] <- M1_list[[station]]$coefficients
  MDO[[station]][['M2']] <- M2_list[[station]]
  MDO[[station]][['vars.M2']] <- M2_list[[station]]$coefficients
  MDO[[station]][['M3']] <- M3_list[[station]]
  MDO[[station]][['vars.M3']] <- M3_list[[station]]$coefficients
  MDO[[station]][['M3']] <- M3_list[[station]]
  MDO[[station]][['vars.M3']] <- M3_list[[station]]$coefficients
  MDO[[station]][['M4']] <- M4_list[[station]]
  MDO[[station]][['vars.M4']] <- M4_list[[station]]$coefficients
  MDO[[station]][['X']] <- X_list[[station]]
} 
rm(list = c('M1_list', 'M2_list', 'M3_list', 'M4_list', 'X_list'))
saveRDS(MDO, 'MDO.rds')
MDO <- readRDS('MDO.rds')



# M5 Y M6 ....
library(gam)
deg_list <- list()
deg_list.new <- list() #0 es igual a log
for (station in estaciones){
  vars <- names(MDO[[station]]$vars.M4)
  vars_era5 <- vars[grepl('z', vars)]
  deg_list[[station]] <- rep(0, times = length(vars_era5))
  names(deg_list[[station]]) <- vars_era5
  deg_list.new[[station]] <- rep(0, times = length(vars_era5))
  names(deg_list.new[[station]]) <- vars_era5
  # for (var in vars_era5){
  #   plot(gam(formula = as.formula(paste('Y ~ s(', var, ')')), data = MDO[[station]]$X),
  #        main = paste(station, var, sep = '-'))
  # }
}


for (station in estaciones) {
  for (var in rev(names(deg_list[[station]]))){
    plot(gam(formula = as.formula(paste('Y ~ s(', var, ')')), data = MDO[[station]]$X),
         main = paste(station, var, sep = '-'))
  }
  cat("Introduce los valores del vector de grados de la estación", station, 
      '\n', "(separados por espacios, y pulsa Enter al final):\n")
  vec <- scan(what = numeric(), quiet = TRUE)
  aux <- names(deg_list[[station]])
  deg_list.new[[station]] <- vec
  names(deg_list.new[[station]]) <- aux
  dev.off()
}

saveRDS(deg_list, 'deg_list.rds')
deg_list <- readRDS('deg_list.rds')

for (station in estaciones){
  cat(paste0(station, '.p.lag') %in% names(MDO[[station]]$vars.M4), '\n')
}

for (station in rev(estaciones)){
  plot(gam(formula = as.formula(paste0('Y ~ s(', station, '.p.lag)')), data = MDO[[station]]$X,
           family = binomial), main = station)
}

deg_lag <- c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
             3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(deg_lag) <- estaciones

M5_list <- list()
for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  deg.list <- deg_list[[station]]
  deg.lag <- deg_lag[station]
  #formula_null <- as.formula(paste0('Y ~ ', paste0('poly(', station, '.p.day, ',deg.day, ')'), '+',
  #                                  paste0('poly(', station, '.p.lag, ',deg.lag, ')')))
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = MDO[[station]]$X)
  
  aux <- paste0('poly(', names(deg.list)[1], ', ', deg.list[1], ')')
  for (i in 2:length(deg.list)){
    aux <- c(aux,  paste0('poly(', names(deg.list)[i], ', ', deg.list[i], ')'))
  }
  
  M5_list[[station]] <- step_rlog(mod_null, 
                                  data = MDO[[station]]$X, 
                                  vars = c(aux, 
                                           paste0('poly(', station, '.p.lag, ', deg.lag, ')')), 
                                  harmonics.l)
  
}

for (station in estaciones){
  MDO[[station]][['M5']] <- M5_list[[station]]
  MDO[[station]][['vars.M5']] <- M5_list[[station]]$coefficients
}

saveRDS(MDO, 'MDO.rds')

# M6 (interacción armonicos y lag) + step
M6_list <- list()
for (station in estaciones){
  M5 <- MDO[[station]]$M5
  vars <- c(labels(terms(M5$formula)), 
            paste0(station, '.p.lag:', c('c.1.l', 's.1.l')))
  vars <- setdiff(vars, unlist(harmonics.l))
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = MDO[[station]]$X)
  
  M6_list[[station]] <- step_rlog(mod_null, 
                                  data = MDO[[station]]$X, 
                                  vars = vars, 
                                  harmonics.l)
}

for (station in estaciones){
  MDO[[station]][['M6']] <- M6_list[[station]]
  MDO[[station]][['vars.M6']] <- M6_list[[station]]$coefficients
}
rm('M6_list')
saveRDS(MDO, 'MDO.rds')


# M7: lag lineal
M7_list <- list()
for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  deg.list <- deg_list[[station]]
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = MDO[[station]]$X)
  
  aux <- paste0('poly(', names(deg.list)[1], ', ', deg.list[1], ')')
  for (i in 2:length(deg.list)){
    aux <- c(aux,  paste0('poly(', names(deg.list)[i], ', ', deg.list[i], ')'))
  }
  
  M7_list[[station]] <- step_rlog(mod_null, 
                                  data = MDO[[station]]$X, 
                                  vars = c(aux, 
                                           paste0(station, '.p.lag')), 
                                  harmonics.l)
}

for (station in estaciones){
  MDO[[station]][['M7']] <- M7_list[[station]]
  MDO[[station]][['vars.M7']] <- M7_list[[station]]$coefficients
}
rm('M7_list')
qsave(MDO, 'MDO.qs')

# M8: lag cuadratico
M8_list <- list()
for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  deg.list <- deg_list[[station]]
  deg.lag <- 2
  #formula_null <- as.formula(paste0('Y ~ ', paste0('poly(', station, '.p.day, ',deg.day, ')'), '+',
  #                                  paste0('poly(', station, '.p.lag, ',deg.lag, ')')))
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = MDO[[station]]$X)
  
  aux <- paste0('poly(', names(deg.list)[1], ', ', deg.list[1], ')')
  for (i in 2:length(deg.list)){
    aux <- c(aux,  paste0('poly(', names(deg.list)[i], ', ', deg.list[i], ')'))
  }
  
  M8_list[[station]] <- step_rlog(mod_null, 
                                  data = MDO[[station]]$X, 
                                  vars = c(aux, 
                                           paste0('poly(', station, '.p.lag, ', deg.lag, ')')), 
                                  harmonics.l)
  
}

for (station in estaciones){
  MDO[[station]][['M8']] <- M8_list[[station]]
  MDO[[station]][['vars.M8']] <- M8_list[[station]]$coefficients
}
rm('M8_list')
qsave(MDO, 'MDO.qs')


# M9: SELECCIÓN CON LOS LAGS
X.MDO <- qread('X.MDO.qs')
M9_list <- list()
for (station in estaciones){
  cat('Estación: ',station, '\n\n')
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = X.MDO[[station]])
  
  M9_list[[station]] <- step_rlog(mod_null,
                                  X_list[[station]],
                                  colnames(X_list[[station]])[15:ncol(X)],
                                  harmonics.l)
  
  
}

for (station in estaciones){
  MDO[[station]][['M9']] <- M9_list[[station]]
  MDO[[station]][['vars.M9']] <- M9_list[[station]]$coefficients
}
rm('M9_list')
qsave(MDO, 'MDO.qs')




#M10: ESTUDIAR TRANSFORMACIONES
# PARA EL LAG ERA DEGREE 3, PERO EN COMMON MODELS LO HEMOS CAMBIADO AL 2
library(gam)
deg_list.lag <- list()
#0 es igual a log
for (station in estaciones[20]){
  vars <- names(MDO[[station]]$vars.M9)
  vars_era5 <- vars[grepl('z', vars)]
  deg_list.lag[[station]] <- rep(0, times = length(vars_era5))
  names(deg_list.lag[[station]]) <- vars_era5
}

for (station in estaciones[20]) {
  for (var in rev(names(deg_list.lag[[station]]))){
    plot(gam(formula = as.formula(paste('Y ~ s(', var, ')')), data = X.MDO[[station]],
             family = binomial),
         main = paste(station, var, sep = '-'))
  }
  cat("Introduce los valores del vector de grados de la estación", station, 
      '\n', "(separados por espacios, y pulsa Enter al final):\n")
  vec <- scan(what = numeric(), quiet = TRUE)
  aux <- names(deg_list.lag[[station]])
  deg_list.lag[[station]] <- vec
  names(deg_list.lag[[station]]) <- aux
  dev.off()
}

qsave(deg_list.lag, 'deg_list.lag.MDO.qs' )

is.lag <- c()
for (station in estaciones){
  is.lag[station] <- paste0(station, '.p.lag') %in% names(MDO[[station]]$vars.M9)
  cat(paste0(station, '.p.lag') %in% names(MDO[[station]]$vars.M9), '\n')
}

for (station in rev(estaciones)){
  station.p <- paste0(station, '.p')
  plot(gam(formula = as.formula(paste0( 'Y ~ s(', station, '.p.lag)')), 
           data = X.MDO[[station]], family = binomial), main = station)
}

#-1 = log
deg_lag.new <- c(3, 3, 3, 3, 3, 3, 3, 3, -1, 3, 3, 3, 3, 3, 2, 3, 3, 3, -1, 3, -1, 
                 3, 3, 3, 3, 3, 3, -1 )
names(deg_lag.new) <- estaciones
deg_lag.new[!is.lag] <- 0

M10_list <- list()
for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  deg.list <- deg_list.lag[[station]]
  deg.lag <- deg_lag.new[station]
  #formula_null <- as.formula(paste0('Y ~ ', paste0('poly(', station, '.p.day, ',deg.day, ')'), '+',
  #                                  paste0('poly(', station, '.p.lag, ',deg.lag, ')')))
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = X.MDO[[station]])
  
  aux <- paste0('poly(', names(deg.list)[1], ', ', deg.list[1], ')')
  for (i in 2:length(deg.list)){
    aux <- c(aux,  paste0('poly(', names(deg.list)[i], ', ', deg.list[i], ')'))
  }
  if(deg.lag == 0){
    vars <- aux
  }else if(deg.lag == -1){
    vars <- c(aux, paste0('I(log(pmax(', paste0(station, '.p.lag,'), '1e-06)))'))
  }else{
    vars <- c(aux, 
              paste0('poly(', station, '.p.lag, ', deg.lag, ')'))
  }
  
  M10_list[[station]] <- step_rlog(mod_null, 
                                  data = X.MDO[[station]], 
                                  vars = vars, 
                                  harmonics.l)
  
}

for (station in estaciones){
  MDO[[station]][['M10']] <- M10_list[[station]]
  MDO[[station]][['vars.M10']] <- M10_list[[station]]$coefficients
}
rm('M10_list')
qsave(MDO, 'MDO.qs')


#----Comunalidades----
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

#plot all AUC and Roc curves: Comparison of models
auc.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 11))
colnames(auc.df) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8',
                      'M9', 'M10')
rownames(auc.df) <- estaciones
auc.df$station <- estaciones

auc.df.pc <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 11))
colnames(auc.df.pc) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8',
                         'M9', 'M10')
rownames(auc.df.pc) <- estaciones
auc.df.pc$station <- estaciones

AIC.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 11))
colnames(AIC.df) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8',
                      'M9', 'M10')
rownames(AIC.df) <- estaciones
AIC.df$station <- estaciones

for (station in estaciones){
  M1 <- MDO[[station]]$M1
  M2 <- MDO[[station]]$M2
  M3 <- MDO[[station]]$M3
  M4 <- MDO[[station]]$M4
  M5 <- MDO[[station]]$M5
  M6 <- MDO[[station]]$M6
  M7 <- MDO[[station]]$M7
  M8 <- MDO[[station]]$M8
  M9 <- MDO[[station]]$M9
  M10 <- MDO[[station]]$M10
  X <- X.MDO[[station]]
  
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = "%Y-%m-%d")
  ind <- which(X$date >= date.first & X$date <= date.last2)
  X_pc <- X[ind, ]
  
  #cat(dim(X_pc)[1], '\n')
  
  roc_M1 <- roc(X$Y, predict(M1, type = 'response'))
  roc_M2 <- roc(X$Y, predict(M2, type = 'response'))
  roc_M3 <- roc(X$Y, predict(M3, type = 'response'))
  roc_M4 <- roc(X$Y, predict(M4, type = 'response'))
  roc_M5 <- roc(X$Y, predict(M5, type = 'response'))
  roc_M6 <- roc(X$Y, predict(M6, type = 'response'))
  roc_M7 <- roc(X$Y, predict(M7, type = 'response'))
  roc_M8 <- roc(X$Y, predict(M8, type = 'response'))
  roc_M9 <- roc(X$Y, predict(M9, type = 'response'))
  roc_M10 <- roc(X$Y, predict(M10, type = 'response'))
  
  roc_M1.pc <- roc(X$Y[ind], M1$fitted.values[ind])
  roc_M2.pc <- roc(X$Y[ind], M2$fitted.values[ind])
  roc_M3.pc <- roc(X$Y[ind], M3$fitted.values[ind])
  roc_M4.pc <- roc(X$Y[ind], M4$fitted.values[ind])
  roc_M5.pc <- roc(X$Y[ind], M5$fitted.values[ind])
  roc_M6.pc <- roc(X$Y[ind], M6$fitted.values[ind])
  roc_M7.pc <- roc(X$Y[ind], M7$fitted.values[ind])
  roc_M8.pc <- roc(X$Y[ind], M8$fitted.values[ind])
  roc_M9.pc <- roc(X$Y[ind], M9$fitted.values[ind])
  roc_M10.pc <- roc(X$Y[ind], M10$fitted.values[ind])
  
  auc.df[station, 2:11] <- round(c(auc(roc_M1), auc(roc_M2), auc(roc_M3), 
                                  auc(roc_M4), auc(roc_M5), auc(roc_M6), 
                                  auc(roc_M7), auc(roc_M8),
                                  auc(roc_M9), auc(roc_M10)), 3)
  auc.df.pc[station, 2:11] <- round(c(auc(roc_M1.pc), auc(roc_M2.pc), auc(roc_M3.pc), 
                                     auc(roc_M4.pc), auc(roc_M5.pc), auc(roc_M6.pc), 
                                     auc(roc_M7.pc), auc(roc_M8.pc),
                                     auc(roc_M9.pc), auc(roc_M10.pc)), 3)
  AIC.df[station, 2:11] <- round(c(AIC(M1), AIC(M2), AIC(M3), 
                                  AIC(M4), AIC(M5), AIC(M6), 
                                  AIC(M7), AIC(M8),
                                  AIC(M9), AIC(M10)), 2)
  
}
# library(writexl)
# write_xlsx(AIC.df, "borrar.xlsx")

library(ggplot2)
library(reshape2)

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
    for (i in 1:nrow(df)) {
      df_rel[i, -1] <- df[i, -1] / df[i, 2]
    }
    
    df_long <- melt(df_rel, id.vars = "station")
    
    # 2️⃣ Calcular valores normalizados por fila solo para la escala de color
    df_long <- df_long %>%
      group_by(station) %>%
      mutate(color_val = (value - min(value)) / (max(value) - min(value))) %>%
      ungroup()
    
    # 3️⃣ Graficar usando color_val para el gradiente y value para el texto
    ggplot(df_long, aes(x = variable, y = station, fill = color_val)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), size = 5) +
      scale_fill_gradientn(
        colors = c("#277DF5", "white", "red"),
        name = "BIC relativo\npor fila",
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

mapa_calor(auc.df, tipo = 'AUC')
mapa_calor(auc.df.pc, tipo = 'AUC')
mapa_calor(AIC.df, tipo = 'AIC')



#extra 
X.MDO <- qread('X.MDO.qs')
station <- estaciones[1]
X <- X.MDO[[station]]
df.less.0 <- data.frame(matrix(NA, ncol = 12, nrow = 28))
rownames(df.less.0) <- estaciones
colnames(df.less.0) <- colnames(X)[grep('zt', colnames(X))]

for (i in 1:length(estaciones)){
  station <- estaciones[i]
  X <- X.MDO[[station]]
  cols <- grep('zt', colnames(X))
  X.sub <- X[, cols]
  cont <- apply(X.sub, 2, function(x) sum(x < 0))
  df.less.0[i, ] <- cont == nrow(X)
}

apply(df.less.0, 2, sum)

