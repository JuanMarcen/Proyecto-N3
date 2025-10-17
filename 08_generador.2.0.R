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
nsim <- 1000
for (station in estaciones){
  mho <- common.models.final[[station]][['MDO']]
  X <- mho$data
  p <- predict(mho, type = 'response')
  
  mean.sims <- replicate(nsim, mean(rbinom(length(p), size = 1, prob = p)))
  
  boxplot(mean.sims, 
          main = paste('MDO:', station),
          ylab = 'Proportion of days with rain')
  points(mean(X$Y), pch = 19, col = 'red', cex = 1.5)
}

