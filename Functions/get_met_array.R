library(tidyverse)

get_met_array <- function(csv, weeks.per.year, year.weeks, met.var, met.uncertainty){
  
  # csv <- "Data/RelativeHumidityWeekly.csv"
  met <- read_csv(csv)
  met <- met %>% 
    filter(Year <= 2018 & Year >= 2016) %>% 
    mutate(yearWeek = as.numeric(gsub("-", "", yearWeek)))
  
  n.sites <- unique(met$siteID) %>% length()
  sites <- unique(met$siteID)
  
  met.var.array <- array(NA, dim = c(3, max(weeks.per.year), length(sites)))
  met.unc.array <- array(NA, dim = c(3, max(weeks.per.year), length(sites)))
  
  for(i in seq_along(sites)){
    met.site <- met %>% 
      filter(siteID == sites[i]) %>% 
      filter(yearWeek %in% year.weeks$yearWeek) %>% 
      select(c(yearWeek, 
               all_of(c(met.var, met.uncertainty))))
    
    met.join <- left_join(year.weeks, met.site)
    met.join.var <- met.join %>% select(all_of(met.var))
    met.join.unc <- met.join %>% select(all_of(met.uncertainty))
    
    for(yy in seq_along(weeks.per.year)){
      if(yy == 1){
        year.index <- 1:weeks.per.year[1]  
      } else {
        year.index <- (weeks.per.year[yy-1]+1):(weeks.per.year[yy-1]+weeks.per.year[yy])  
      }
      met.center <- met.join.var[year.index, ] %>% pull()
      met.center <- met.center - mean(met.center, na.rm = TRUE)
      met.var.array[yy, 1:weeks.per.year[yy], i] <- met.center
      met.unc.array[yy, 1:weeks.per.year[yy], i] <- as.matrix(met.join.unc[year.index, ])
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
