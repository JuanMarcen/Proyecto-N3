rm(list = ls())

load('data.RData')

# harmonics (podría guardarlo ya con armónicos)
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

formula <- as.formula(paste(station.p, '~',  paste(colnames(X_final)[8:ncol(X_final)], collapse = '+')))
formula

fit <- glm(formula,
           family = Gamma(link = 'log'),
           data = X_final)
summary(fit)
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
fit <- gamlss(formula, sigma.fo = ~ mes, family = GA, data = X_final)

sum <- summary(fit)

plot(1 / fit$sigma.fv, type = 'l')


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
lines(mu, col = "blue", lwd = 2)              # media predicha
points(y_sim, col = rgb(1, 0, 0, 0.5), pch = 16)  # simulaciones
legend("topright", legend = c("Observado", "Media predicha", "Simulado"),
       col = c("black", "blue", "red"), pch = c(16, NA, 16), lty = c(NA, 1, NA))

plot(density(X_final[[station.p]]), col = 'blue', lwd = 2)
lines(density(y_sim), col = 'red', lwd = 2)

qqplot(X_final[[station.p]], y_sim,
       main = "Q-Q Plot Observado vs Simulado",
       xlab = "Cuantiles observados", ylab = "Cuantiles simulados",
       pch = 16, col = "blue")

abline(0, 1, col = "red", lwd = 2)
