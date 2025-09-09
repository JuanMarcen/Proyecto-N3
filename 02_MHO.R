# model hourly ocurrence
# substract days in which it doesn't rain
rm(list = ls())

load('data.RData')


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

# creation of final X matrix -- Generalize
# for all stations
# Saving design matrix (X) and model 
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
    mutate(Y = ifelse(.data[[station.p]] > 0, 1, 0)) %>% 
    relocate(Y, .after = !!station.p) %>%
    as.data.frame() %>% na.omit()
  
  #na omit necesario?
  
  # MHO
  formula <- as.formula(
    paste('Y ~', paste(colnames(X_final)[8:ncol(X_final)], collapse = '+'))
  )
  
  mho <- glm(formula = formula, family = binomial(logit), data = X_final)
  
  assign(paste0('X.', station), X_final)
  assign(paste0('mho.', station), mho)
  
  rm(list = c('X', 'X_final', 'mho', 'formula', 'p_day', 'station.p', 'station'))
}

# evaluation
library(pROC)

for (station in estaciones){
  mho <- get(paste0('mho.', station))
  X <- get(paste0('X.', station))
  
  summary(mho)
  
  boxplot(mho$fitted.values ~ X$Y,
          xlab = 'Y real',
          ylab = 'P(Y = 1| X)',
          main = station)
  
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

roc$thresholds
coords(roc, 'best')[1]

