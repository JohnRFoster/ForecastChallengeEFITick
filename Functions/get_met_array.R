library(tidyverse)

get_met_array <- function(csv, weeks.per.year, filter.week, year.weeks, met.var, met.uncertainty){
  
  if("cumGDD" %in% met.var){
    daily.bio.temp <- read_csv("Data/BioTemperatureDaily.csv")
    met <- daily.bio.temp %>% 
      group_by(siteID, Year) %>% 
      mutate(growingDegree = if_else(bioTempMaximum > 0, bioTempMaximum, 0)) %>% 
      mutate(cumGDD = cumsum(growingDegree)) %>% 
      mutate(yearWeek = as.numeric(paste0(Year, epiWeek))) %>% 
      rename(cumGDDVariance = bioTempMaximumVariance) %>% 
      group_by(siteID, Year, epiWeek) %>% 
      slice(which.min(cumGDD)) %>% 
      filter(Year >= 2016) %>%
      filter(yearWeek < filter.week) %>%
      select(c(Year,
               yearWeek,
               siteID,
               cumGDDVariance,
               cumGDD))
    
  } else {
    met <- read_csv(csv)
    met <- met %>% 
      filter(Year >= 2016) %>% 
      mutate(yearWeek = as.numeric(gsub("-", "", yearWeek))) %>% 
      filter(yearWeek < filter.week) 
  }
  
  met$yearWeek <- as.character(met$yearWeek)
  
  n.sites <- unique(met$siteID) %>% length()
  sites <- unique(met$siteID)
  
  year.col <- NA
  for(i in seq_along(weeks.per.year)) year.col <- c(year.col, rep(i, weeks.per.year[i]))
  year.col <- year.col[-1]
    
  met.var.array <- array(NA, dim = c(4, max(weeks.per.year), length(sites)))
  met.unc.array <- array(NA, dim = c(4, max(weeks.per.year), length(sites)))
  
  for(i in seq_along(sites)){
    met.site <- met %>% 
      filter(siteID == sites[i]) %>% 
      filter(yearWeek %in% year.weeks$yearWeek) %>% 
      select(c(yearWeek, Year,
               all_of(c(met.var, met.uncertainty))))
    
    met.join <- left_join(year.weeks, met.site)
    met.join.var <- met.join %>% select(all_of(met.var))
    met.join.unc <- met.join %>% select(all_of(met.uncertainty))
    
    met.join.var <- cbind(met.join.var, year.col)
    met.join.unc <- cbind(met.join.unc, year.col)
    
    for(yy in seq_along(weeks.per.year)){

      met.center <- met.join.var %>% 
        filter(year.col == yy) %>% 
        pull(met.var)
      met.var.array[yy, 1:weeks.per.year[yy], i] <- met.center
      
      met.unc.pull <- met.join.unc %>% 
        filter(year.col == yy) %>%  
        pull(met.uncertainty) 
      
      if("cumGDD" %in% met.var){
        met.unc.array[yy, 1:weeks.per.year[yy], i] <- cumsum(met.unc.pull)
      } else {
        met.center <- met.center - mean(met.center, na.rm = TRUE)
        met.unc.array[yy, 1:weeks.per.year[yy], i] <- met.unc.pull  
      }
      met.var.array[yy, 1:weeks.per.year[yy], i] <- met.center
    }
  }
  
  x.inits <- met.var.array
  for(yy in seq_along(weeks.per.year)){
    if(any(is.na(x.inits[yy, 1:weeks.per.year[yy], ]))){
      for(cc in 1:length(sites)){
        na.rows <- which(is.na(x.inits[yy,1:weeks.per.year[yy],cc]))
        x.inits[yy, na.rows, cc] <- approx(as.vector(x.inits[yy,1:weeks.per.year[yy],cc]), xout = na.rows)$y  
      }
    }
  }
  
  return(list(met.vals = met.var.array,
              met.uncertainty = met.unc.array,
              met.inits = x.inits,
              site.dim = sites))
}



# going to want a site vector to match to tick data
