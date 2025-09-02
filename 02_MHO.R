# model hourly ocurrence
# substract days in which it doesn't rain
rm(list = ls())

load('data.RData')

# creation of design matrix (covariates)
# example for 1 station
station <- 'EM71.p'

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

# creation of final X matrix
X_example <- df_hours[, c('t', 'l', 'mes', 'dia.mes', 'h', station,
                           colnames(harm_l)[2:ncol(harm_l)],
                           colnames(harm_h)[2:ncol(harm_h)])]

# whole day
p_day <- df_days[,  c('t', 'l', 'mes', 'dia.mes', station)]
colnames(p_day)[colnames(p_day) == station] <- paste0(station, '.day')

X_example <- X_example %>%
  left_join(p_day, by = c('t', 'l', 'mes', 'dia.mes'))

# lags
X_example <- X_example %>% 
  mutate(!!paste0(station, '.lag') := lag(.data[[station]]))

colnames(X_example)

#eliminate days with no rain
X_final <- X_example[-which(X_example$EM71.p.day == 0), ]

# hourly indicator of rain
X_final <- X_final %>%
  mutate(Y = ifelse(.data[[station]] > 0, 1, 0)) %>% 
  relocate(Y, .after = !!station) %>%
  as.data.frame() %>% na.omit()

# MHO
formula <- as.formula(
  paste('Y ~', paste(colnames(X_final)[8:ncol(X_final)], collapse = '+'))
)

mho_example <- glm(formula = formula, family = binomial, data = X_final)
summary(mho_example)

boxplot(mho_example$fitted.values ~ X_final$Y)
