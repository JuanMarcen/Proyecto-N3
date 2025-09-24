# MODEL HOURLY QUANTITY (GLM GAMMA)

rm(list = ls())

load('data.RData')

# harmonics (podría guardarlo ya con armónicos)
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

#----data.-----
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
    filter(.data[[station.p]] > 0) %>%
    na.omit() %>%
    as.data.frame
  
  X_list[[station]] <- X_final
}

#----models fitting----
# M1: harm + CV ct
# M2: rain var + CV ct
# M3: M1 + M2
# M4: M1 + M2 + CV varying
# M5: M1 + transf. var rain + CV varying

# grado polinomio (M5)
library(gam)
for (station in estaciones){
  X <- X_list[[station]]
  aux.marcador <- is.element(X$mes, c(6,7,8)) #meses de mayor lluvia al parecer
  #variable lag
  mod.lag <- gam(formula = as.formula(paste0(paste0(station, '.p'), '~ s(', station, '.p.lag)')), data = X[aux.marcador, ])
  plot(mod.lag, main = paste(station))
  #variable del dia
  mod.day <- gam(formula = as.formula(paste0(paste0(station, '.p'), '~ s(', station, '.p.day)')), data = X[aux.marcador, ])
  plot(mod.day, main = paste(station))
}

degrees.p.lag <- c(3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 2, 3,
                   3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3)
names(degrees.p.lag) <- estaciones
degrees.p.day <- c(1, 3, 3, 3, 3, 3, 3, 3, 1, 3, 1, 1, 1, 1, 1, 3,
                   3, 2, 2, 3, 1, 1, 3, 1, 1, 2, 1, 1)
names(degrees.p.day) <- estaciones

# models M1--M5 fitting
library(gamlss)
M1_list <- list()
M2_list <- list()
M3_list <- list()
M4_list <- list()
M5_list <- list()
for (station in estaciones){
  cat('Ajuste modelos de la estación: ', station, '\n')
  
  X <- X_list[[station]]
  station.p <- paste0(station, '.p')
  
  # modelos con CV constante
  
  formula_M1 <- as.formula(paste(station.p, '~', paste(c(colnames(harm_h)[2:ncol(harm_h)],
                                                colnames(harm_l)[2:ncol(harm_l)]),
                                              collapse = '+')))
  formula_M2 <- as.formula(paste(station.p, '~', paste0(station.p, '.day'), '+', paste0(station.p, '.lag')))
  
  formula_M3 <- as.formula(
    paste(station.p, '~', paste(colnames(X)[7:ncol(X)], collapse = '+'))
  )
  
  cat('Ajuste modelo M1: ', deparse(formula_M1), '\n')
  M1_list[[station]] <- glm(formula = formula_M1, family = Gamma(link = 'log'), data = X)
  
  cat('Ajuste modelo M2: ', deparse(formula_M2), '\n')
  M2_list[[station]] <- glm(formula = formula_M2, family = Gamma(link = 'log'), data = X)
  
  cat('Ajuste modelo M3: ', deparse(formula_M3), '\n')
  M3_list[[station]] <- glm(formula = formula_M3, family = Gamma(link = 'log'), data = X)
  
  #modelos con CV variante segun armonicos
  deg.lag <- degrees.p.lag[station]
  deg.day <- degrees.p.day[station]
  
  formula_M4 <- formula_M3
  
  cat('Ajuste modelo M4: ', deparse(formula_M4), '\n')
  M4_list[[station]] <- gamlss(formula_M4, 
                               sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                                 I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
                               family = GA, 
                               data = X,
                               trace = F)
  
  formula_M5 <- as.formula(paste(station.p, '~', 
                                 paste(c(colnames(harm_h)[2:ncol(harm_h)],
                                                colnames(harm_l)[2:ncol(harm_l)]),
                                              collapse = '+'), 
                                 '+' ,
                                 paste(c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                         paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                         paste0(station,'.p.lag:',station,'.p.day')), 
                                       collapse = '+')))
  
  cat('Ajuste modelo M5: ', deparse(formula_M5), '\n')
  M5_list[[station]] <- gamlss(formula_M5, 
                               sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                                 I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
                               family = GA, 
                               data = X,
                               trace = F)
  
  
}

MHQ <- list()
for (station in estaciones){
  MHQ[[station]][['M1']] <- M1_list[[station]]
  MHQ[[station]][['vars.M1']] <- M1_list[[station]]$coefficients
  MHQ[[station]][['M2']] <- M2_list[[station]]
  MHQ[[station]][['vars.M2']] <- M2_list[[station]]$coefficients
  MHQ[[station]][['M3']] <- M3_list[[station]]
  MHQ[[station]][['vars.M3']] <- M3_list[[station]]$coefficients
  MHQ[[station]][['M4']] <- M4_list[[station]]
  MHQ[[station]][['vars.M4']] <- M4_list[[station]]$mu.coefficients
  MHQ[[station]][['M5']] <- M5_list[[station]]
  MHQ[[station]][['vars.M5']] <- M5_list[[station]]$mu.coefficients
  MHQ[[station]][['X']] <- X_list[[station]]
} 
rm(list = c('M1_list', 'M2_list', 'M3_list', 'M4_list', 'M5_list', 'X_list',
            'p_day', 'X_final', 'mod.day', 'mod.lag'))
saveRDS(MHQ, 'MHQ.rds')

# model selection
# similar to MHQ but changing the function

step_glm <- function(initial_model,
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
  #     mod.temp <- gamlss(formula.aux, 
  #                        sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
  #                          I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
  #                        family = GA, 
  #                        data = data,
  #                        trace = F)
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
  
  
  mod.aux <- step(mod.aux, scope = scope.aux, direction = 'both', trace = FALSE)
  
  cat('Model after ', 'AIC ', 'step algorithm for variables: ', deparse(formula(mod.aux)), '\n')
  cat('AIC: ', AIC(mod.aux), '\n')
  
  
  for (h in harmonics.l){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- gamlss(formula.aux, 
                       sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                         I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
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
  
  for (h in harmonics.h){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- gamlss(formula.aux, 
                       sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                         I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
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

# armónicos
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



M6_list <- list()
for (station in estaciones[1]){
  cat('Estación ', station, '\n\n')
  
  deg.lag <- degrees.p.lag[station]
  deg.day <- degrees.p.day[station]
  
  mod_null <- gamlss(as.formula(paste(paste0(station, '.p'), '~ 1')), 
                     sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                       I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
                     family = GA, 
                     data = MHQ[[station]]$X,
                     trace = F)
  
  M6_list[[station]] <- step_glm(mod_null, 
                                   data = MHQ[[station]]$X, 
                                   vars = c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                            paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                            paste0(station,'.p.lag:',station,'.p.day')), 
                                   harmonics.l, harmonics.h)
  
}





for (station in estaciones[1]){
  MHQ[[station]][['M6']] <- M6_list[[station]]
  MHQ[[station]][['vars.M6']] <- M6_list[[station]]$mu.coefficients
} 

saveRDS(MHQ, 'MHQ.rds')
MHQ <- readRDS('MHQ.rds')

#----model comparison----
AIC.df <- data.frame(matrix(NA, nrow = length(estaciones), ncol = 7))
colnames(AIC.df) <- c('station', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6')
rownames(AIC.df) <- estaciones
AIC.df$station <- estaciones
for (station in estaciones){
  M1 <- MHQ[[station]]$M1
  M2 <- MHQ[[station]]$M2
  M3 <- MHQ[[station]]$M3
  M4 <- MHQ[[station]]$M4
  M5 <- MHQ[[station]]$M5
  M6 <- MHQ[[station]]$M6
  X <- MHQ[[station]]$X
  AIC.df[station, 2:7] <- round(c(AIC(M1), AIC(M2), AIC(M3), 
                                  AIC(M4), AIC(M5), AIC(M6)), 2)
  
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
  geom_text(aes(label = round(value, 3)), size = 5) +
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


#----model control----
M6 <- MHQ[[station]]$M6
M5 <- MHQ[[station]]$M5
X <- MHQ[[station]]$X

shape <- 1/M5$sigma.fv^2 
plot((X$R036.p - M5$mu.fv) / (sqrt(M5$mu.fv^2/shape)))
points(X$R036.p - M6$mu.fv, col = 'red')

for (station in estaciones[1]){
 station.p <- paste0(station, '.p')
  shape <- 1 / MHQ[[station]]$M6$sigma.fv^2
  rate <- shape / MHQ[[station]]$M6$mu.fv
  
  hist(MHQ[[station]]$X[, station.p], breaks = 50, prob = T, main = station)
  lines(density(MHQ[[station]]$X[, station.p]), col = 'blue')
  y_sim <- rgamma(length(shape), shape = shape, rate = rate)
  lines(density(y_sim), col = 'red', lwd = 2)
 }

station <- estaciones[12]
station.p <- paste0(station, '.p')
shape <- 1 / mhq_list[[station]]$sigma.fv^2
rate <- shape / mhq_list[[station]]$mu.fv
plot(shape, type = 'l')
hist(MHQ[[station]]$X[, station.p], breaks = 50, prob = T, main = station)
lines(density(MHQ[[station]]$X[, station.p]), col = 'blue')
y_sim <- rgamma(length(shape), shape = shape, rate = rate)
lines(density(y_sim), col = 'red', lwd = 2)


# comparison with uniform (0,1)
# example 1 station
station <- estaciones[1]

m <- MHQ[[station]]$M6
X <- MHQ[[station]]$X
p.obs <- X[[paste0(station, '.p')]]

mu <- m$mu.fv
shape <- 1 / m$sigma.fv
rate <- shape / mu

u <- pgamma(p.obs, shape = shape, rate = rate)

plot(density(u, from = 0, to = 1), col = 'blue', lwd = 2)
lines(density(runif(length(u), 0, 1), from = 0, to = 1), col = 'red', lwd = 2)

library(overlapping)
overlap(list(observed = u,
             theorical = runif(length(u), 0, 1)),
        type = '1',
        plot = T)

u_sorted <- sort(u)

# Cuantiles teóricos de una uniforme(0,1)
n <- length(u_sorted)
theoretical <- (1:n) / (n + 1)  # usar (i)/(n+1) es común para evitar 0 y 1 exactos

# QQ-plot
plot(theoretical, u_sorted, 
     main = "QQ-plot vs Uniforme(0,1)", 
     xlab = "Cuantiles teóricos U(0,1)", 
     ylab = "Cuantiles empíricos de u")
abline(0, 1, col = "red", lwd = 2) # línea de referencia


# 100 simulaciones y mirar quantiles
m <- MHQ[[station]]$M6
X <- MHQ[[station]]$X
p.obs <- X[[paste0(station, '.p')]]

mu <- m$mu.fv
shape <- 1 / m$sigma.fv
rate <- shape / mu

plot(density(p.obs), col = 'blue', lwd = 2)
for (i in 1:100){
  u <- rgamma(length(p.obs), shape = shape, rate = rate)
  lines(density(u), col = 'red')
}
lines(density(p.obs), col = 'blue', lwd = 2)

cuantiles <- c('q0.05', 'q0.50', 'q0.90', 'q0.95', 'q0.99')
q.obs <- quantile(p.obs, probs = c(0.05, 0.5, 0.90, 0.95, 0.99))
names(q.obs) <- cuantiles
q.obs

q.sim <- data.frame(matrix(NA, ncol = 5))
colnames(q.sim) <- cuantiles
for (i in 1:100){
  u <- rgamma(length(p.obs), shape = shape, rate = rate)
  q <- quantile(u, probs = c(0.05, 0.5, 0.90, 0.95, 0.99))
  names(q) <- cuantiles
  q.sim <- rbind(q.sim, q)
}
q.sim <- q.sim[-1, ]

plot(rep(q.obs[1], times = dim(q.sim)[1]), q.sim$q0.05, xlim = c(q.obs[1]-0.5, q.obs[5]+0.5), ylim = c(0, 10))
points(rep(q.obs[2], times = dim(q.sim)[1]), q.sim$q0.50)
points(rep(q.obs[3], times = dim(q.sim)[1]), q.sim$q0.90)
points(rep(q.obs[4], times = dim(q.sim)[1]), q.sim$q0.95)
points(rep(q.obs[5], times = dim(q.sim)[1]), q.sim$q0.99)

boxplot(q.sim,
        at = q.obs,                # 👈 coloca cada boxplot en la posición correspondiente
        names = round(q.obs, 2),   # (opcional) etiquetas en el eje x
        xlim = c(q.obs[1]-0.5, q.obs[5]+0.5),
        ylim = c(0, max(q.sim)),
        col = "lightblue",
        main = "Boxplots alineados con q.obs",
        ylab = "Valores simulados",
        xlab = "Cuantiles observados")
#------
