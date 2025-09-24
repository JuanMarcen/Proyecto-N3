# MODELS OF DAIRLY DATA --> OCURRENCE

rm(list = ls())

load('data.RData')

library(lubridate)
# In the covariates matrix we add the climate variables
global_df <- readRDS('global_df.rds')
global_df$t <- year(global_df$date)
global_df <- global_df[, -which(colnames(global_df) %in% c('zg300.', 'zg500.', 'zg700.', 
                                                           'zt300.', 'zt500.', 'zt700.'))]

# cálculo de periodo comun
n <- length(estaciones)
aux <- apply(df_days[, paste0(estaciones, '.p')], 1, function(row) all(!is.na(row)))

first <- which(aux)[1]
last <- tail(which(aux), 1)
date.first <- as.Date(paste0(df_days[first, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
date.last <- as.Date(paste0(df_days[last, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')

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
  formula_M2 <- as.formula(paste('Y ~', paste(
    colnames(X_list[[station]])[15:ncol(X_list[[station]])], collapse = '+'
  )))
  
  formula_M3 <- as.formula(
    paste('Y ~', paste(
      colnames(X_list[[station]])[7:ncol(X_list[[station]])], collapse = '+')
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

M4_list <- list()
for (station in estaciones){
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

saveRDS(MDO, 'MDO.rds')
MDO <- readRDS('MDO.rds')

# M5 Y M6 ....
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

#plot all AUC and Roc curves: Comparison of models
auc.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 5))
colnames(auc.df) <- c('station', 'M1', 'M2', 'M3', 'M4')
rownames(auc.df) <- estaciones
auc.df$station <- estaciones

auc.df.pc <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 5))
colnames(auc.df.pc) <- c('station', 'M1', 'M2', 'M3', 'M4')
rownames(auc.df.pc) <- estaciones
auc.df.pc$station <- estaciones

AIC.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 5))
colnames(AIC.df) <- c('station', 'M1', 'M2', 'M3', 'M4')
rownames(AIC.df) <- estaciones
AIC.df$station <- estaciones

for (station in estaciones){
  M1 <- MDO[[station]]$M1
  M2 <- MDO[[station]]$M2
  M3 <- MDO[[station]]$M3
  M4 <- MDO[[station]]$M4
  # M5 <- MDO[[station]]$M5
  # M6 <- MDO[[station]]$M6
  X <- MDO[[station]]$X
  
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = "%Y-%m-%d")
  ind <- which(X$date >= date.first & X$date <= date.last)
  X_pc <- X[ind, ]
  
  #cat(dim(X_pc)[1], '\n')
  
  roc_M1 <- roc(X$Y, predict(M1, type = 'response'))
  roc_M2 <- roc(X$Y, predict(M2, type = 'response'))
  roc_M3 <- roc(X$Y, predict(M3, type = 'response'))
  roc_M4 <- roc(X$Y, predict(M4, type = 'response'))
  # roc_M5 <- roc(X$Y, predict(M5, type = 'response'))
  # roc_M6 <- roc(X$Y, predict(M6, type = 'response'))
  
  roc_M1.pc <- roc(X$Y[ind], M1$fitted.values[ind])
  roc_M2.pc <- roc(X$Y[ind], M2$fitted.values[ind])
  roc_M3.pc <- roc(X$Y[ind], M3$fitted.values[ind])
  roc_M4.pc <- roc(X$Y[ind], M4$fitted.values[ind])
  # roc_M5.pc <- roc(X$Y[ind], M5$fitted.values[ind])
  # roc_M6.pc <- roc(X$Y[ind], M6$fitted.values[ind])
  
  auc.df[station, 2:7] <- round(c(auc(roc_M1), auc(roc_M2), auc(roc_M3), 
                                  auc(roc_M4), auc(roc_M5), auc(roc_M6)), 3)
  auc.df.pc[station, 2:7] <- round(c(auc(roc_M1.pc), auc(roc_M2.pc), auc(roc_M3.pc), 
                                     auc(roc_M4.pc), auc(roc_M5.pc), auc(roc_M6.pc)), 3)
  AIC.df[station, 2:7] <- round(c(AIC(M1), AIC(M2), AIC(M3), 
                                  AIC(M4), AIC(M5), AIC(M6)), 2)
  # BIC.df[station, 2:7] <- round(c(BIC(M1), BIC(M2), BIC(M3), 
  #                                 BIC(M4), BIC(M5), BIC(M6)), 2)
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
mapa_calor(AIC.df[, -5], tipo = 'AIC')
mapa_calor(BIC.df[, -5], tipo = 'AIC')
