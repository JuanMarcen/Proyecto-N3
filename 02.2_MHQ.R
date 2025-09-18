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

# data
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

# model selection
# similar to MHO but changing the function
library(gamlss)
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
  if (!is.null(vars)) {
    for (var in vars){
      formula.aux <- update(formula(mod.aux), paste(". ~ . +", var))
      mod.temp <- gamlss(formula.aux, 
                         sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                           I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
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

# grado polinomio
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

degrees.p.lag <- c(3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 2, 3)
names(degrees.p.lag) <- estaciones
degrees.p.day <- c(1, 3, 3, 3, 3, 3, 3, 3, 1, 3, 1, 1, 1, 1, 1, 3)
names(degrees.p.day) <- estaciones

mhq_list <- list()
for (station in estaciones){
  cat('Estación ', station, '\n\n')
  
  deg.lag <- degrees.p.lag[station]
  deg.day <- degrees.p.day[station]
  
  mod_null <- gamlss(as.formula(paste(paste0(station, '.p'), '~ 1')), 
                     sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + 
                       I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
                     family = GA, 
                     data = X_list[[station]],
                     trace = F)
  
  mhq_list[[station]] <- step_glm(mod_null, 
                                   data = X_list[[station]], 
                                   vars = c(paste0('poly(', station, '.p.day, ',deg.day, ')'), 
                                            paste0('poly(', station, '.p.lag, ',deg.lag, ')'),
                                            paste0(station,'.p.lag:',station,'.p.day')), 
                                   harmonics.l, harmonics.h)
  
}




MHQ <- list()
for (station in estaciones){
  MHQ[[station]][['mho']] <- mhq_list[[station]]
  MHQ[[station]][['vars']] <- mhq_list[[station]]$coefficients
  MHQ[[station]][['X']] <- X_list[[station]]
} 

saveRDS(MHQ, 'MHQ.rds')

# prueba para una estación. Juguetear
station <- estaciones[1]
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

hist(X_final[, station.p], breaks = 50, prob = T)
lines(density(X_final[, station.p]), col = 'blue')

formula <- as.formula(paste(station.p, '~',  paste(colnames(X_final)[8:(ncol(X_final) - 2)], 
                                                   collapse = '+'), '+',
                            paste0('poly(', station, '.p.day, ',3, ')'), '+' ,
                            paste0('poly(', station, '.p.lag, ',3, ')'), '+',
                            paste0(station,'.p.lag:',station,'.p.day')))
formula

fit <- glm(formula,
           family = Gamma(link = 'log'),
           data = X_final)
summary(fit)

plot(fit)


fit$offset
sum <- summary(fit)
r <- 1 / sum$dispersion #estimation of shape parameter

mu <- exp(fit$fitted.values)
lambda <- r / mu

curve(dgamma(x, shape=r, rate = mean(lambda)), 
      col="red", lwd=2, add=TRUE)

plot(residuals(fit, type="deviance"))
qqnorm(residuals(fit, type="pearson")); qqline(residuals(fit, type="pearson"))


library(gamlss)
fit <- gamlss(formula, sigma.fo = ~ I(sin(2*pi*l/365)) + I(cos(2*pi*l/365)) + I(sin(2*pi*h/24)) + I(cos(2*pi*h/24)), 
              family = GA, data = X_final)

sum <- summary(fit)
plot(fit)



#rate
plot(1 / fit$sigma.fv^2, type = 'l')

coef.var <- tapply(fit$sigma.fv, INDEX = X_final$mes, mean)
plot(coef.var, type = 'l')

media_mensual <- tapply(1/fit$sigma.fv^2, INDEX = X_final$mes, mean)
plot(media_mensual, type = 'l')

for (station in estaciones){
  station.p <- paste0(station, '.p')
  shape <- 1 / mhq_list[[station]]$sigma.fv^2
  rate <- shape / mhq_list[[station]]$mu.fv
  
  hist(X_list[[station]][, station.p], breaks = 50, prob = T, main = station)
  lines(density(X_list[[station]][, station.p]), col = 'blue')
  y_sim <- rgamma(length(shape), shape = shape, rate = rate)
  lines(density(y_sim), col = 'red', lwd = 2)
}


#------

# chat ? 
# Ajuste con link log
fit_log <- glm(formula,
               family = Gamma(link = "log"),
               data = X_final)

# Ajuste con link inverse
fit_inv <- glm(formula,
               family = Gamma(link = "inverse"),
               data = X_final,
               start = rep(0.1, times = length(coef(fit_log))))

summary(fit_log)
summary(fit_inv)

pred_log <- predict(fit_log, type = "response")   # en escala de la media
pred_inv <- predict(fit_inv, type = "response")

# Comparar primeras filas
head(data.frame(pred_log, fit$mu.fv, real = X_final[[station.p]]), 10)
AIC(fit_log, fit_inv)


phi <- summary(fit_log)$dispersion  # estimación de la dispersión
alpha <- 1/phi                      # shape
mu <- predict(fit_log, type = "response")
scale <- mu / alpha

# Simulación de valores posibles de Y
y_sim <- rgamma(length(mu), shape = alpha, scale = scale)
plot(X_final[[station.p]], type = "p", pch = 16, col = "black",
    main = "Valores reales vs predichos y simulados",
    xlab = "Índice", ylab = "Respuesta")
aux <- cbind(X_final[[station.p]], y_sim)
head(aux)
lines(mu, col = "blue", lwd = 2)              # media predicha
points(y_sim, col = rgb(1, 0, 0, 0.5), pch = 16)  # simulaciones
legend("topright", legend = c("Observado", "Media predicha", "Simulado"),
       col = c("black", "blue", "red"), pch = c(16, NA, 16), lty = c(NA, 1, NA))

plot(density(X_final[[station.p]]), col = 'blue', lwd = 2)
lines(density(y_sim), col = 'red', lwd = 2)



