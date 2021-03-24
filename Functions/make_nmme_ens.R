library(ncdf4)
library(dplyr)
library(lubridate)


#' @param var variable to grab, also a dir, ex "tasmin"
#' @param start.month forecast to grab, also a dir, ex "20190201"
#' @param dec.lat decimal latitude from NEON data 
#' @param dec.lon decimal longitude from NEON data 
#' @return matrix for the given nmme variable

make_nmme_ens <- function(var, start.month, dec.lat, dec.lon){
  
  # convert to get correct ncdf index
  dec.lon <- dec.lon %% 360 
  dec.lon <- round(dec.lon) # nmme has one deg resolution
  
  # get correct lat index
  lat.index.nc <- -90:90 # ncdf values
  dec.lat <- which(round(dec.lat) == lat.index.nc) # match
  
  # list ensemble members for given variable and start month
  path <- file.path("Data/NMME", start.month, var)
  ens <- list.files(path)
  
  nmme.ens <- matrix(NA, 365, length(ens))
  for(e in seq_along(ens)){
    nc <- nc_open(file.path(path, ens[e])) # open file
    name <- names(nc$var) # variable name
    nmme <- ncvar_get(nc, name)
    
    # daily forecast; everyday from begining of start month
    nmme.ens[,e] <- nmme[dec.lon, dec.lat, ]
    
  }
  
  nmme.days <- ymd(start.month) + nc$dim$TIME$vals # every day in nmme forecast
  epiWeek <- epiweek(nmme.days) # convert to epiWeek
  
  colnames(nmme.ens) <- paste0("ens", seq_along(ens))
  
  nmme.return <- cbind(epiWeek, nmme.ens)
  
  return(nmme.return)
}


