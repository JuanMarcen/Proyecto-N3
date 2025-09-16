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
