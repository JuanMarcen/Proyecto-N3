# model hourly ocurrence
# substract days in which it doesn't rain
rm(list = ls())

load('data.RData')

# cálculo de periodo comun
library(lubridate)
n <- length(estaciones)
aux <- apply(df_hours[, paste0(estaciones, '.p')], 1, function(row) all(!is.na(row)))

first <- which(aux)[1]
last <- tail(which(aux), 1)
date.first <- as.Date(paste0(df_hours[first, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')
date.last <- as.Date(paste0(df_hours[last, c('dia.mes', 'mes', 't')], collapse = '-'), format = '%d-%m-%Y')

# harmonics
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
h <- 0:23
harm_l <- data.frame(l = l,
                     cs(l, harmonics = 1:4, 365))
colnames(harm_l)[2:ncol(harm_l)] <- paste0(colnames(harm_l)[2:ncol(harm_l)], '.l')
harm_h <- data.frame(h = h,
                        cs(0:23, harmonics = 1:4, 24))
colnames(harm_h)[2:ncol(harm_h)] <- paste0(colnames(harm_h)[2:ncol(harm_h)], '.h')


library(dplyr)
df_hours <- df_hours %>%
  left_join(harm_l, by = 'l') %>%
  left_join(harm_h, by = 'h')


# creation of final X matrix -- Generalize
# for all stations
# Saving design matrix (X) and models M1, M2 and M3
# M1: only harmonics
# M2: only variables
# M3 : M1 + M2
M1_list <- list()
M2_list <- list()
M3_list <- list()
X_list <- list()
for (station in estaciones){
  cat('Ajuste modelos de la estación: ', station, '\n')
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
  
  formula_M1 <- as.formula(paste('Y ~', paste(c(colnames(harm_h)[2:ncol(harm_h)],
                                                colnames(harm_l)[2:ncol(harm_l)]),
                                              collapse = '+')))
  formula_M2 <- as.formula(paste('Y ~', paste0(station.p, '.day'), '+', paste0(station.p, '.lag')))
  
  formula_M3 <- as.formula(
    paste('Y ~', paste(colnames(X_list[[station]])[8:ncol(X_list[[station]])], collapse = '+'))
  )
  
  cat('Ajuste modelo M1: ', deparse(formula_M1), '\n')
  M1_list[[station]] <- glm(formula = formula_M1, family = binomial(logit), data = X_list[[station]])
  
  cat('Ajuste modelo M2: ', deparse(formula_M2), '\n')
  M2_list[[station]] <- glm(formula = formula_M2, family = binomial(logit), data = X_list[[station]])
  
  cat('Ajuste modelo M3: ', deparse(formula_M3), '\n')
  M3_list[[station]] <- glm(formula = formula_M3, family = binomial(logit), data = X_list[[station]])
  # assign(paste0('X.', station), X_final)
  # assign(paste0('mho.', station), mho)
  
  # rm(list = c('X', 'X_final', 'mho', 'formula', 'p_day', 'station.p', 'station'))
}

MHO <- list()
for (station in estaciones){
  MHO[[station]][['M1']] <- M1_list[[station]]
  MHO[[station]][['vars.M1']] <- M1_list[[station]]$coefficients
  MHO[[station]][['M2']] <- M2_list[[station]]
  MHO[[station]][['vars.M2']] <- M2_list[[station]]$coefficients
  MHO[[station]][['M3']] <- M3_list[[station]]
  MHO[[station]][['vars.M3']] <- M3_list[[station]]$coefficients
  MHO[[station]][['X']] <- X_list[[station]]
}

rm(list = c('M1_list', 'M2_list', 'M3_list', 'X_list', 'harm_h', 'harm_l', 'X_final', 'X', 'p_day'))
# estudio de covariables en los modelos
## ejemplo para EM71
# mod <- get(paste0('mho.', estaciones[1]))
# X <- get(paste0('X.', estaciones[1]))
# mod <- mho_list[[estaciones[15]]]
# X <- X_list[[estaciones[15]]]
# mod_rebuilt <- glm(formula = formula(mod), family = binomial(link = "logit"), data = X)
# step(mod_rebuilt, direction = "backward")
# 
# 
# #borrador
# basura <- update(mod, data = X, formula = .~. - P024.p.lag:P024.p.day)
# basura <- update(mho.EM71, data = X.EM71, formula = .~. + 
#                    poly(EM71.p.day, 2) + poly(EM71.p.lag, 2))

library(gam)
# aux.marcador <- is.element(X$mes, c(6,7,8)) #meses de mayor lluvia al parecer
# mod1 <- gam(formula = as.formula('Y ~ s(EM71.p.day)'), data = X[aux.marcador, ])
# summary(mod1)
# plot(gam(formula = as.formula('Y ~ s(EM71.p.day)'), data = X[aux.marcador, ]))
# plot(gam(formula = as.formula('Y ~ s(EM71.p.lag, df = 3)'), data = X[aux.marcador, ]))

#decisión del grado de polinomio (a mano)
for (station in estaciones[17:28]){
  X <- MHO[[station]]$X
  aux.marcador <- is.element(X$mes, c(6,7,8)) #meses de mayor lluvia al parecer
  #variable lag
  mod.lag <- gam(formula = as.formula(paste0('Y ~ s(', station, '.p.lag)')), data = X[aux.marcador, ])
  plot(mod.lag, main = paste(station))
  #variable del dia
  mod.day <- gam(formula = as.formula(paste0('Y ~ s(', station, '.p.day)')), data = X[aux.marcador, ])
  plot(mod.day, main = paste(station))
}

degrees.p.lag <- c(3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(degrees.p.lag) <- estaciones
degrees.p.day <- c(3, 1, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
names(degrees.p.day) <- estaciones

# M4: VARIABLES CON INTERACCIÓN Y POLY()
M4_list <- list()
for (station in estaciones){
  cat('Estación: ', station, '\n')
  X <- MHO[[station]]$X
  
  deg.lag <- degrees.p.lag[station]
  deg.day <- degrees.p.day[station]
  
  formula_M4 <- as.formula(paste('Y ~', paste(c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                                paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                                paste0(station,'.p.lag:',station,'.p.day')), 
                                              collapse = '+')))
  
  cat('Ajuste modelo M4: ', deparse(formula_M4), '\n')
  M4_list[[station]] <- glm(formula_M4, family = binomial(logit), data = X)
}

for (station in estaciones){
  MHO[[station]][['M4']] <- M4_list[[station]]
  MHO[[station]][['vars.M4']] <- M4_list[[station]]$coefficients
} 
rm('M4_list')
saveRDS(MHO, 'MHO.rds')

#----selección de variables --> M5 ----
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
                      harmonics.h, step = 'AIC'){ 
  # argumentos:
  # initial_model: modelo inicial (modelo nulo)
  # data: datos a utilizar (contiene explicada y explicativas)
  # vars: covariables no armónicas que se desean introducir en el modelo (char)
  # incluidas siempre? pues al modelo initial 
  # harmonics: lista de armónicos
  
  if (step == 'AIC'){
    k <- 2
  }else if (step == 'BIC'){
    n <- dim(data)[1]
    k <- log(n)
  }
  
  mod.aux <- initial_model
  cat('Initial model ', deparse(formula(mod.aux)), '\n')
  cat(step, ': ', AIC(mod.aux, k = k), '\n')
  
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
  
  #selección de variable por método step directamente
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
  
  cat('\nFinal model: ', deparse(formula(mod.aux)), '\n')
  cat(step, ': ', AIC(mod.aux, k = k), '\n\n')
  
  return(mod.aux)
}

#actualización de modelos
M5_list <- list()
M6_list <- list()
for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  deg.lag <- degrees.p.lag[station]
  deg.day <- degrees.p.day[station]
  
  #formula_null <- as.formula(paste0('Y ~ ', paste0('poly(', station, '.p.day, ',deg.day, ')'), '+',
  #                                  paste0('poly(', station, '.p.lag, ',deg.lag, ')')))
  
  mod_null <- glm(Y ~ 1, family = binomial(logit), data = MHO[[station]]$X)
  
  M5_list[[station]] <- step_rlog(mod_null, 
                                   data = MHO[[station]]$X, 
                                   vars = c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                            paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                            paste0(station,'.p.lag:',station,'.p.day')), 
                                   harmonics.l, harmonics.h)
  
  M6_list[[station]] <- step_rlog(mod_null, 
                                  data = MHO[[station]]$X, 
                                  vars = c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                           paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                           paste0(station,'.p.lag:',station,'.p.day')), 
                                  harmonics.l, harmonics.h,
                                  step = 'BIC')
}

# for (station in estaciones){
#   print(summary(mho_list[[station]]))
# }

# guardado de una lista definitiva por estacion
# estacion --> modelo, df, variables escogidas
for (station in estaciones){
  MHO[[station]][['M5']] <- M5_list[[station]]
  MHO[[station]][['vars.M5']] <- M5_list[[station]]$coefficients
  MHO[[station]][['M6']] <- M6_list[[station]]
  MHO[[station]][['vars.M6']] <- M6_list[[station]]$coefficients
} 

saveRDS(MHO, 'MHO.rds')
rm(list = c('M6_list', 'M5_list', 'X'))
MHO <- readRDS('MHO.rds')

# SELECCIÓN AUTOMÁTICA DE MODELO SEGÚN ESTUDIO DE MODELOS ANIDADOS (STAND BY)
data <- X_list[[15]]
initial_model <- glm(Y ~ 1, family = binomial(logit), data = data)

station <- estaciones[15]
vars <- c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
          paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
          paste0(station,'.p.lag:',station,'.p.day'))

#ALGORTIHM DID NOT CONVERGE?
for (i in 1:length(vars)){
  mod.aux <- update(initial_model, data = data, formula = as.formula(paste(". ~ . +", vars[i])),
                    control = glm.control(maxit = 50))
  
  aux <- anova(initial_model, mod.aux)
  pval <- aux$`Pr(>Chi)`[2]
  
  if (pval < 0.05){
    initial_model <- mod.aux
    cat('Variable añadida: ', vars[i], '\n')
  }
}


#----Estudio comunalidades (DIRIA QUE MEJOR PARA DIARIAS)----
# stand by 

#----evaluation #put in a Rmarkdow----
library(pROC)

for (station in estaciones){
  mho <- mho_list[[station]]
  X <- MHO[[station]]$X
  
  # summary(mho)
  # 
  # boxplot(mho$fitted.values ~ X$Y,
  #         xlab = 'Y real',
  #         ylab = 'P(Y = 1| X)',
  #         main = station)
  
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


#plot all AUC and Roc curves: Comparison of models
auc.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 7))
colnames(auc.df) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6')
rownames(auc.df) <- estaciones
auc.df$station <- estaciones

auc.df.pc <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 7))
colnames(auc.df.pc) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6')
rownames(auc.df.pc) <- estaciones
auc.df.pc$station <- estaciones

AIC.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 7))
colnames(AIC.df) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6')
rownames(AIC.df) <- estaciones
AIC.df$station <- estaciones

BIC.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 7))
colnames(BIC.df) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6')
rownames(BIC.df) <- estaciones
BIC.df$station <- estaciones

for (station in estaciones){
  M1 <- MHO[[station]]$M1
  M2 <- MHO[[station]]$M2
  M3 <- MHO[[station]]$M3
  M4 <- MHO[[station]]$M4
  M5 <- MHO[[station]]$M5
  M6 <- MHO[[station]]$M6
  X <- MHO[[station]]$X
  
  X$date <- as.Date(paste(X$t, X$mes, X$dia.mes, sep = "-"), format = "%Y-%m-%d")
  ind <- which(X$date >= date.first & X$date <= date.last)
  X_pc <- X[ind, ]
  
  #cat(dim(X_pc)[1], '\n')
  
  roc_M1 <- roc(X$Y, predict(M1, type = 'response'))
  roc_M2 <- roc(X$Y, predict(M2, type = 'response'))
  roc_M3 <- roc(X$Y, predict(M3, type = 'response'))
  roc_M4 <- roc(X$Y, predict(M4, type = 'response'))
  roc_M5 <- roc(X$Y, predict(M5, type = 'response'))
  roc_M6 <- roc(X$Y, predict(M6, type = 'response'))
  
  roc_M1.pc <- roc(X$Y[ind], M1$fitted.values[ind])
  roc_M2.pc <- roc(X$Y[ind], M2$fitted.values[ind])
  roc_M3.pc <- roc(X$Y[ind], M3$fitted.values[ind])
  roc_M4.pc <- roc(X$Y[ind], M4$fitted.values[ind])
  roc_M5.pc <- roc(X$Y[ind], M5$fitted.values[ind])
  roc_M6.pc <- roc(X$Y[ind], M6$fitted.values[ind])
  
  # plot(roc_M1,
  #      col = 'red',
  #      print.auc = T,
  #      print.auc.y = .5,
  #      main = paste('Comparación modelos', station))
  # plot(roc_M2,
  #      col = 'blue',
  #      add = TRUE,
  #      print.auc = T,
  #      print.auc.y = .45)
  # plot(roc_M3,
  #      col = 'green',
  #      add = TRUE,
  #      print.auc = T,
  #      print.auc.y = .4)
  # plot(roc_M4,
  #      col = 'purple',
  #      add = TRUE,
  #      print.auc = T,
  #      print.auc.y = .35)
  # plot(roc_M5,
  #      col = 'orange',
  #      add = TRUE,
  #      print.auc = T,
  #      print.auc.y = .3)
  # 
  # legend('topleft', lty = 1, lwd = 2,
  #        legend = c('M1', 'M2', 'M3', 'M4', 'M5'),
  #        col = c('red', 'blue', 'green', 'purple', 'orange'))
  
  auc.df[station, 2:7] <- round(c(auc(roc_M1), auc(roc_M2), auc(roc_M3), 
                               auc(roc_M4), auc(roc_M5), auc(roc_M6)), 3)
  auc.df.pc[station, 2:7] <- round(c(auc(roc_M1.pc), auc(roc_M2.pc), auc(roc_M3.pc), 
                                  auc(roc_M4.pc), auc(roc_M5.pc), auc(roc_M6.pc)), 3)
  AIC.df[station, 2:7] <- round(c(AIC(M1), AIC(M2), AIC(M3), 
                               AIC(M4), AIC(M5), AIC(M6)), 2)
  BIC.df[station, 2:7] <- round(c(BIC(M1), BIC(M2), BIC(M3), 
                               BIC(M4), BIC(M5), BIC(M6)), 2)
}
# library(writexl)
# write_xlsx(AIC.df, "borrar.xlsx")

library(ggplot2)
library(reshape2)

# mapa_calor <- function(df){
#   df_long <- melt(df, id.vars = "station")
#   df_long <- df_long %>%
#     group_by(station) %>%
#     mutate(rank = rank(-value, ties.method = "first")) %>%  # rank 1 = máximo
#     ungroup()
#   
#   df_long <- df_long %>%
#     mutate(color_cat = case_when(
#       rank == 1 ~ "Max",
#       rank == 2 ~ "Segundo",
#       rank == 3 ~ "Medio1",
#       rank == 4 ~ "Medio2",
#       rank == 5 ~ "Penúltimo",
#       rank == 6 ~ "Mínimo"
#     ))
#   
#   # Definir colores manualmente
#   colores <- c(
#     "Max" = "red",
#     "Segundo" = "tomato",
#     "Medio1" = "#F0C2BB",
#     "Medio2" = "#C7D8F0",
#     "Penúltimo" = "#7EACED",
#     "Mínimo" = "#277DF5"
#   )
#   
#   # Graficar
#   ggplot(df_long, aes(x = variable, y = station, fill = color_cat)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = round(value, 3)), size = 5) +
#     scale_fill_manual(values = colores) +
#     theme_minimal() +
#     theme(
#       axis.title = element_blank(),
#       legend.position = "none",
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid = element_blank()
#     )
#   
# }

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

plot(auc.df$M1, auc.df$M2, ylim = c(0.5, 1), xlim = c(0.5, 1))
points(auc.df$M1, auc.df$M3, col = 'blue')
points(auc.df$M3, auc.df$M5, col = 'purple')
abline(a = 0, b = 1, col = 'red')
