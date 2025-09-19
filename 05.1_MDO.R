# MODELS OF DAIRLY DATA --> OCURRENCE

rm(list = ls())

load('data.RData')

# In the covariates matrix we add the climate variables
global_df <- readRDS('global_df.rds')
global_df$t <- year(global_df$date)
global_df <- global_df[, -which(colnames(global_df) %in% c('zg300.', 'zg500.', 'zg700.', 
                                                           'zt500.', 'zt700.'))]

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

formula <- as.formula(
  paste('Y ~', paste(colnames(X)[7:ncol(X)], collapse = '+'))
)

mdo <- glm(formula, family = binomial(logit), data = X)

harmonics.l <- list(
  h1 = c('s.1.l', 'c.1.l'),
  h2 = c('s.2.l', 'c.2.l'),
  h3 = c('s.3.l', 'c.3.l'),
  h4 = c('s.4.l', 'c.4.l')
)

mod_null <- glm(Y ~ 1, family = binomial(logit), data = X)

step_rlog <- function(initial_model,
                      data,
                      vars,
                      harmonics.l){ 
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
      mod.temp <- glm(formula.aux, data = data, family = binomial(link = "logit"))
      
      if (AIC(mod.temp) < AIC(mod.aux)){
        cat('Added: ', var, '\n')
        cat('AIC: ', AIC(mod.temp), '\n')
        mod.aux <- mod.temp
      }
    }
    
  }
  
  
  for (h in harmonics.l){
    
    formula.aux <- update(formula(mod.aux), paste(". ~ . +", paste(h, collapse = '+')))
    
    mod.temp <- glm(formula.aux, data = data, family = binomial(logit))
    
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


basura <- step_rlog(mod_null,
                    X,
                    colnames(X)[7:ncol(X)],
                    harmonics.l)

library(gam)
plot(gam(formula = as.formula(paste('Y ~ s(', colnames(X)[36], ')')), data = X))

summary(basura)
