# METRICS FOR COMPARISON OF MODELS WITH RESPECT TO THE OBSERVED
# SPATIAL COMPARISON, in order to see how well they behave globally

# SCC: Spatial Correlation Coefficient
SCC <- function(x.obs, y.sim){
  mean.obs <- mean(x.obs)
  mean.sim <- mean(y.sim)
  
  dif.obs <- x.obs - mean.obs
  dif.sim <- y.sim - mean.sim
  
  num <- sum(dif.obs * dif.sim)
  denom <- sqrt(sum(dif.obs^2)) * sqrt(sum(dif.sim^2))
  
  return(num / denom)
}


# RMSE: Root Mean Square Error
RMSE <- function(x.obs, y.sim){
  dif.sq <- sum((x.obs - y.sim)^2)
  
  return((sqrt(dif.sq / 28)))
}

# RB: Relative Bias
RB <- function(x.obs, y.sim){
  dif <- sum(x.obs - y.sim)
  
  return(dif / sum(x.obs) * 100)
}