#!/usr/bin/Rscript

# Clear workspace
rm(list = ls())

# Read stations information
stations <- read.csv("Data/ECAD/stations.csv")
stations$Grid <- paste0(round(stations$LAT),"N.",abs(round(stations$LON)),ifelse(stations$LON>-0.5,"E","W"))

# Target grid points
# Corner points
cpoints <- c("45N.10W","45N.5E","35N.10W","35N.5E")
# Near-station points
tpoints <- c(cpoints,unique(stations$Grid))
# Reference period
ref_period <- c("1981-06-01","2010-08-31")

# Read grid data
era5 <- cbind(read.csv("Data/ERA5/g500_grid.csv"),
              read.csv("Data/ERA5/g700_grid.csv"))
# Convert ERA5 geopotentials into geopotential height
era5_vars <- grep("g",names(era5))
era5[,era5_vars] <- era5[,era5_vars]/9.80665
names(era5)[era5_vars] <- paste0("z",names(era5)[era5_vars])

# Reshape data
era5$date <- as.Date(era5$date)
TT <- length(table(format(era5$date,"%Y")))
LL <- nrow(era5)/TT
SS <- nrow(stations)

# Init global data.frame
gdf <- data.frame(STAID = rep(stations$STAID, each = nrow(era5)),
                  date = rep(era5$date, SS),
                  t = rep(rep(1:TT, each = LL), SS),
                  l = rep(rep(1:LL, TT), SS),
                  LAT = rep(stations$LAT, each = nrow(era5)),
                  LON = rep(stations$LON, each = nrow(era5)),
                  grid = rep(stations$Grid, each = nrow(era5)),
                  Altitude = rep(stations$HGHT, each = nrow(era5)),
                  CoastDist = rep(stations$CoastDist, each = nrow(era5)))

# Variable names
eralevels <- paste0("zg",c("500","700"))
# Corner points
cpoints <- c("45N.10W","45N.5E","35N.10W","35N.5E")
# Define variable names
vnames <- paste(rep(eralevels, length(cpoints)+1),
                rep(c("",cpoints),each = length(eralevels)),
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
  fivep <- c(stations$Grid[ss], cpoints)
  eranames <- paste(rep(eralevels, length(fivep)),
                    rep(fivep, each = length(eralevels)),
                    sep = ".")
  # Add to matrix
  gmat[cidx1:cidx2,] <- era5[,match(eranames, names(era5))]
  
}# for ss stations

# Combine global data.frame with ERA5 variables
gdf <- cbind(gdf,gmat)

# Write data.frame
saveRDS(object = gdf,
        file = "Data/ERA5/global_df.rds")