#!/usr/bin/Rscript

# Clear workspace
rm(list = ls())

# Read stations information
stations <- readRDS("stations.rds")
stations$Grid <- paste0(round(stations$LAT),"N.",abs(round(stations$LON)),ifelse(stations$LON>-0.5,"E","W"))

# Target grid points
# Corner points, last 4 are the ones we will use
cpoints <- c("45N.10W","45N.5E","35N.10W","35N.5E", "42N.2W", "42N.1W", "41N.2W", "41N.1W")
cpoints.4 <- c("42N.2W", "42N.1W", "41N.2W", "41N.1W")
# Near-station points
tpoints <- unique(c(cpoints,unique(stations$Grid)))
tpoints <- cpoints.4
# Reference period
ref_period <- c("1997-10-01","2023-12-31")

# Read grid data
era5 <- cbind(readRDS("datos_elsa/g300.rds"), readRDS('datos_elsa/g500.rds'), readRDS('datos_elsa/g700.rds'))
# Convert ERA5 geopotentials into geopotential height
#era5 <- era5[-which(format(era5$Date, '%m-%d') == '02-29'), ]
era5_vars <- grep("g",names(era5))
era5[,era5_vars] <- era5[,era5_vars]/9.80665
names(era5)[era5_vars] <- paste0("z",names(era5)[era5_vars])

# Reshape data
era5$Date <- as.Date(era5$Date)
TT <- length(table(format(era5$Date,"%Y")))
LL <- nrow(era5)/TT
SS <- nrow(stations)

# Init global data.frame
gdf <- data.frame(STAID = rep(stations$STAID, each = nrow(era5)),
                  date = rep(era5$Date, SS),
                  t = rep(rep(1:TT, each = LL), SS),
                  l = rep(rep(1:LL, TT), SS),
                  LAT = rep(stations$LAT, each = nrow(era5)),
                  LON = rep(stations$LON, each = nrow(era5)),
                  grid = rep(stations$Grid, each = nrow(era5)))

# Variable names
eralevels <- paste0("zg",c("300", "500","700"))
# Corner points
cpoints <- c("45N.10W","45N.5E","35N.10W","35N.5E")
cpoints.4 <- c("42N.2W", "42N.1W", "41N.2W", "41N.1W")
# Define variable names
vnames <- paste(rep(eralevels, length(cpoints.4)+1),
                rep(c("",cpoints.4),each = length(eralevels)),
                sep = ".")

# Global variables
gmat <- matrix(data = as.numeric(NA),
               nrow = nrow(gdf),
               ncol = length(vnames))
colnames(gmat) <- vnames
gmat <- as.data.frame(gmat)

# Read closed point for each station and concatenate in the global matrix
for(ss in 1:nrow(stations)){
  
  cat(paste0("..",ss))
  
  # Concatenating indices
  cidx1 <- 1+nrow(era5)*(ss-1)
  cidx2 <- nrow(era5)*ss
  
  # Include variables of interest
  fivep <- c(stations$Grid[ss], cpoints.4)
  eranames <- paste(rep(eralevels, length(fivep)),
                    rep(fivep, each = length(eralevels)),
                    sep = ".")
  # Add to matrix
  gmat[cidx1:cidx2,] <- era5[,match(eranames, names(era5))]
  
}# for ss stations

# Combine global data.frame with ERA5 variables
gdf <- cbind(gdf,gmat)
gdf <- gdf[gdf$date >= ref_period[1] & gdf$date <= ref_period[2], ]

#now the TEMPERATURES
# Target grid points
# Corner points, last 4 are the ones we will use
cpoints.4 <- c("42N.2W", "42N.1W", "41N.2W", "41N.1W")
# Near-station points
tpoints <- cpoints.4
# Reference period
# Read grid data
era5 <- cbind(readRDS('datos_elsa/t500.rds'), readRDS('datos_elsa/t700.rds'))
# Convert ERA5 geopotentials into geopotential height
#era5 <- era5[-which(format(era5$Date, '%m-%d') == '02-29'), ]
era5_vars <- grep("t",names(era5))
era5_vars <- era5_vars[-1]
era5[,era5_vars] <- era5[,era5_vars] - 273.15 #temp in ºC
names(era5)[era5_vars] <- paste0("z",names(era5)[era5_vars])

# Reshape data
era5$Date <- as.Date(era5$Date)
TT <- length(table(format(era5$Date,"%Y")))
LL <- nrow(era5)/TT
SS <- nrow(stations)

# Init global data.frame
tdf <- data.frame(STAID = rep(stations$STAID, each = nrow(era5)),
                  date = rep(era5$Date, SS),
                  t = rep(rep(1:TT, each = LL), SS),
                  l = rep(rep(1:LL, TT), SS),
                  LAT = rep(stations$LAT, each = nrow(era5)),
                  LON = rep(stations$LON, each = nrow(era5)),
                  grid = rep(stations$Grid, each = nrow(era5)))

# Variable names
eralevels <- paste0("zt",c( "500","700"))
# Corner points
cpoints.4 <- c("42N.2W", "42N.1W", "41N.2W", "41N.1W")
# Define variable names
vnames <- paste(rep(eralevels, length(cpoints.4)+1),
                rep(c("",cpoints.4),each = length(eralevels)),
                sep = ".")

# Global variables
gmat <- matrix(data = as.numeric(NA),
               nrow = nrow(tdf),
               ncol = length(vnames))
colnames(gmat) <- vnames
gmat <- as.data.frame(gmat)

# Read closed point for each station and concatenate in the global matrix
for(ss in 1:nrow(stations)){
  
  cat(paste0("..",ss))
  
  # Concatenating indices
  cidx1 <- 1+nrow(era5)*(ss-1)
  cidx2 <- nrow(era5)*ss
  
  # Include variables of interest
  fivep <- c(stations$Grid[ss], cpoints.4)
  eranames <- paste(rep(eralevels, length(fivep)),
                    rep(fivep, each = length(eralevels)),
                    sep = ".")
  # Add to matrix
  gmat[cidx1:cidx2,] <- era5[,match(eranames, names(era5))]
  
}# for ss stations

# Combine global data.frame with ERA5 variables
tdf <- cbind(tdf,gmat)

tdf <- tdf[tdf$date >= ref_period[1] & tdf$date <= ref_period[2], ]

global_df <- cbind(gdf, tdf[, c(8:ncol(tdf))])
# Write data.frame
saveRDS(object = global_df,
        file = "global_df.rds")

