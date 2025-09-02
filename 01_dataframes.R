# Datafrmae creation
rm(list = ls())
library(lubridate)
library(dplyr)

# Data importation
data_dir <- 'C:/Users/jumar/OneDrive/Escritorio/N3/datos confe'

files <- list.files(data_dir, full.names = TRUE)

estaciones <- c('EM71', 'P084', 'P018', 'E085', 'R026',
                'EM08', 'R037', 'R062',
                'EM70', 'P021', 'P087', 'P088',
                'A042', 'P023', 'P024', 'R036')

for (i in estaciones){
  
  file <- files[grepl(i, files) & grepl("PQUIN", files)]
  name <- paste0('df.',i)
  
  assign(name, read.csv(file, skip = 15, header = TRUE, sep = ';'))
  
  df <- get(name)
  
  df <- df %>%
    mutate(
      # date and hour
      fecha_hora = dmy_hms(paste(Fecha, Hora), tz = "UTC"),
      # value column as numeric. --- = NA
      value = as.numeric(ifelse(Valor..mm. == "---", NA, gsub(',', '.', Valor..mm.)))
    )
  
  # no 29 february
  df <- df %>%
    filter(!(month(fecha_hora) == 2 & day(fecha_hora) == 29))
  
  #rest of  desired columns
  df <- df %>%
    mutate(
      t = year(fecha_hora),
      mes = month(fecha_hora),
      dia.mes = day(fecha_hora),
      h = hour(fecha_hora),
      min = minute(fecha_hora),
      # days from 1-365
      l = {
        l0  <- yday(fecha_hora)
        adj <- if_else(leap_year(fecha_hora) & month(fecha_hora) > 2, 1L, 0L)
        as.integer(l0 - adj)
      }
    ) %>%
    select(t, l, mes, dia.mes, h, min, value)
  
  names(df)[names(df) == 'value'] <- paste0(i, '.p')
  
  assign(name, df)
  
}

# final data frame of minutes
df_minutes <- get(paste0('df.', estaciones[1]))

for(i in 2: length(estaciones)){
  df_minutes <- full_join(df_minutes, 
                      get(paste0('df.', estaciones[i])),
                      by = c('t', 'l', 'mes', 'dia.mes', 'h', 'min'))
}

df_minutes <- df_minutes %>%
  arrange(t, mes, dia.mes, h, min)
