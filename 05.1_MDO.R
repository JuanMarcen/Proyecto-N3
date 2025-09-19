# MODELS OF DAIRLY DATA --> OCURRENCE

rm(list = ls())

load('data.RData')

# In the covariates matrix we add the climate variables
global_df <- readRDS('global_df.rds')
global_df$t <- year(global_df$date)

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

station <- estaciones[1]
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

X <- X %>% 
  mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0)) %>%
  relocate(Y, .after = !!station.p) %>%
  as.data.frame() %>%
  na.omit()



