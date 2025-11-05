# ocurrence of precipitacion simulators
rm(list = setdiff(ls(), c('MHQ', 'MHO', 'MDQ', 'MDO', 'common.models.final',
                          'per.comun.day', 'per.comun.h', 'estaciones')))

# previous day simulation only necessary in quantity models

# MHO
par(mfrow = c(2, 7))
set.seed(05052002)
nsim <- 100
for(station in estaciones){
  mho <- common.models.final[[station]][['MHO']]
  #mho <- MHO[[station]]$M9
  X <- mho$data
  p <- predict(mho, type = 'response')
  
  mean.sims <- replicate(nsim, mean(rbinom(length(p), size = 1, prob = p)))
  
  boxplot(mean.sims, 
          main = paste('MHO:', station),
          ylab = 'Proportion of hours with rain')
  points(mean(X$Y), pch = 19, col = 'red', cex = 1.5)
}



# MDO
set.seed(05052002)
n.sim <- 100
for (station in estaciones){
  mdo <- common.models.final[[station]][['MDO']]
  X <- mdo$data
  X <- X %>% filter(t == 2015)
  p <- predict(mdo, newdata = X, type = 'response')
  
  
  mean.sims <- replicate(n.sim, mean(rbinom(length(p), size = 1, prob = p)))
  
  boxplot(mean.sims, 
          main = paste('MDO:', station),
          ylab = 'Proportion of days with rain')
  points(mean(X$Y), pch = 19, col = 'red', cex = 1.5)
}



#MDQ
Xq <- X.MDQ[[station]]
Xq$date <- as.Date(paste(Xq$t, Xq$mes, Xq$dia.mes, sep = "-"),
                   format = '%Y-%m-%d')
Xq <- Xq %>% filter(date >= per.comun.day[1] & date <= per.comun.day[2])
Xq <- Xq[, -ncol(Xq)]
X.filt <- Xq %>% filter(t == 2015)
mq <- common.models.final[[station]][['MDQ']]

# formulas rewritten for new dataframe
aux.mu.formula <- mq$mu.formula
aux.mu.formula <- as.formula(
  paste(as.character(aux.mu.formula[2]), '~', 
        paste(labels(terms(aux.mu.formula)), collapse = '+'))
)

aux.sigma.formula <- mq$sigma.formula
aux.sigma.formula <- as.formula(
  paste('~', as.character(aux.sigma.formula[2]))
)


mu.form <- sanitize_formula(aux.mu.formula)
sigma.form <- sanitize_formula(aux.sigma.formula)

#re ajuste modelo cantidad (debido a las nuevas formulas)
mq2 <- gamlss(mu.form, sigma.fo = sigma.form,
              family = GA, data = Xq, trace = FALSE)


mu.fv <- predict(mq2, newdata = X.filt, what = 'mu', type = 'response',
                 data = Xq)
sigma.fv <- predict(mq2, newdata = X.filt, what = 'sigma', type = 'response',
                 data = Xq)

shape.fv <- 1 / sigma.fv ^2
rate.fv <- shape.fv / mu.fv

replicate(n.sim, rgamma(length(mu.fv), shape = shape.fv, rate = rate.fv))
