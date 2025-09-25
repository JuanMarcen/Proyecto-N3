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

#----Data building and M1, M2, m3, M4----
library(gamlss)
X_list <- list()
M1_list <- list()
M2_list <- list()
M3_list <- list()
M4_list <- list()
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
    filter(.data[[station.p]] > 0) %>%
    as.data.frame() %>%
    na.omit()
  
  formula_M1 <- as.formula(paste(station.p, '~', paste(c(colnames(harm_l)[2:ncol(harm_l)]),
                                                      collapse = '+')))
  formula_M2 <- as.formula(paste(station.p, '~', paste(
    colnames(X_list[[station]])[15:ncol(X_list[[station]])], collapse = '+'
  )))
  
  formula_M3 <- as.formula(
    paste(station.p, '~', paste(
      colnames(X_list[[station]])[7:ncol(X_list[[station]])], collapse = '+')
    ))
  
  cat('Ajuste modelo M1: ', deparse(formula_M1), '\n')
  M1_list[[station]] <- glm(formula = formula_M1, family = Gamma(link = 'log'), data = X_list[[station]])
  
  cat('Ajuste modelo M2: ', deparse(formula_M2), '\n')
  M2_list[[station]] <- glm(formula = formula_M2, family = Gamma(link = 'log'), data = X_list[[station]])
  
  cat('Ajuste modelo M3: ', deparse(formula_M3), '\n')
  M3_list[[station]] <- glm(formula = formula_M3, family = Gamma(link = 'log'), data = X_list[[station]])
  
  #modelos con CV variante segun armonicos
  
  formula_M4 <- formula_M3
  
  cat('Ajuste modelo M4: ', deparse(formula_M4), '\n')
  M4_list[[station]] <- gamlss(formula_M4, 
                               sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)), 
                               family = GA, 
                               data = X_list[[station]],
                               trace = F)
  
}

MDQ <- list()
for (station in estaciones){
  MDQ[[station]][['M1']] <- M1_list[[station]]
  MDQ[[station]][['vars.M1']] <- M1_list[[station]]$coefficients
  MDQ[[station]][['M2']] <- M2_list[[station]]
  MDQ[[station]][['vars.M2']] <- M2_list[[station]]$coefficients
  MDQ[[station]][['M3']] <- M3_list[[station]]
  MDQ[[station]][['vars.M3']] <- M3_list[[station]]$coefficients
  MDQ[[station]][['M4']] <- M4_list[[station]]
  MDQ[[station]][['vars.M4']] <- M4_list[[station]]$mu.coefficients
  MDQ[[station]][['X']] <- X_list[[station]]
} 
rm(list = c('M1_list', 'M2_list', 'M3_list', 'M4_list', 'X_list'))
saveRDS(MDQ, 'MDQ.rds')
MDQ <- readRDS('MDQ.rds')

#----model selection----
harmonics.l <- list(
  h1 = c('s.1.l', 'c.1.l'),
  h2 = c('s.2.l', 'c.2.l'),
  h3 = c('s.3.l', 'c.3.l'),
  h4 = c('s.4.l', 'c.4.l')
)

library(gamlss)
step_glm <- function(initial_model,
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
      mod.temp <- gamlss(formula.aux, 
                         sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)), 
                         family = GA, 
                         data = data,
                         trace = F)
      
      if (AIC(mod.temp) < AIC(mod.aux)){
        cat('Added: ', var, '\n')
        cat('AIC: ', AIC(mod.temp), '\n')
        mod.aux <- mod.temp
      }
    }
    
  }
  
  
  for (h in harmonics.l){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- gamlss(formula.aux, 
                       sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)), 
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
  
 
  cat('\nFinal model: ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n\n')
  
  return(mod.aux)
}


mdq_list <- list()

for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  mod_null <- gamlss(as.formula(paste(paste0(station, '.p'), '~ 1')), 
                     sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)), 
                     family = GA, 
                     data = X_list[[station]],
                     trace = F)
  
  mdq_list[[station]] <- step_glm(mod_null, 
                                  data = X_list[[station]], 
                                  vars = colnames(X_list[[station]])[15:ncol(X)], 
                                  harmonics.l)
  
}



MDQ <- list()
for (station in estaciones){
  MDQ[[station]][['mdo']] <- mdq_list[[station]]
  MDQ[[station]][['vars']] <- mdq_list[[station]]$coefficients
  MDQ[[station]][['X']] <- X_list[[station]]
} 

saveRDS(MDQ, 'MDQ.rds')


for (station in estaciones){
  station.p <- paste0(station, '.p')
  shape <- 1 / mdq_list[[station]]$sigma.fv^2
  rate <- shape / mdq_list[[station]]$mu.fv
  
  hist(X_list[[station]][, station.p], breaks = 50, prob = T, main = station)
  lines(density(X_list[[station]][, station.p]), col = 'blue')
  y_sim <- rgamma(length(shape), shape = shape, rate = rate)
  lines(density(y_sim), col = 'red', lwd = 2)
}


#----model comparison----
AIC.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 5))
colnames(AIC.df) <- c('station', 'M1', 'M2', 'M3', 'M4')
rownames(AIC.df) <- estaciones
AIC.df$station <- estaciones
for (station in estaciones){
  M1 <- MDQ[[station]]$M1
  M2 <- MDQ[[station]]$M2
  M3 <- MDQ[[station]]$M3
  M4 <- MDQ[[station]]$M4
  # M5 <- MDQ[[station]]$M5
  # M6 <- MDQ[[station]]$M6
  X <- MDQ[[station]]$X
  AIC.df[station, 2:5] <- round(c(AIC(M1), AIC(M2), AIC(M3), 
                                  AIC(M4)), 2)
  
}

library(ggplot2)
library(reshape2)

df_rel <- AIC.df
for (i in 1:nrow(AIC.df)) {
  df_rel[i, -1] <- AIC.df[i, -1] / AIC.df[i, 2]
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
  geom_text(aes(label = round(value, 4)), size = 5) +
  scale_fill_gradientn(
    colors = c("#277DF5", "white", "red"),
    name = "AIC relativo\npor fila",
    breaks = c(0, 0.5, 1),
    labels = c("min", "", "max")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

