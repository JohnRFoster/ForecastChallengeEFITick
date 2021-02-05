# ====================================================== #
#              NEON driver data download                 #
# ====================================================== #

library(neonstore)
library(tidyverse)
library(lubridate)

# products we want
driver.products <- c(
  "DP1.00098.001", # relative humidity 
  "DP1.00003.001", # triple aspirated air temperature
  "DP1.00005.001", # IR biological temperature (ground temp)
  "DP1.00094.001", # soil water content and water salinity
  "DP1.00041.001"  # soil temperature
)

make.rh <- FALSE
make.taat <- FALSE
make.ir.temp <- FALSE
make.soil.temp <- FALSE
make.soil.water <- TRUE

# target sites
target.sites <- c("BLAN", "ORNL", "SCBI", "SERC", "KONZ", "TALL", "UKFS")

message("Updating NEON Data...")
for(i in seq_along(driver.products)){
  neon_download(product = driver.products[i],
                end_date = "2019-12-31",   
                site = target.sites,       
                .token = Sys.getenv("NEON_TOKEN"))  
}

# set up end date
day.run <- lubridate::today() # the day the script is called

# end dates for each month. seq of start dates then subtract 1
date.end.month <- seq(as.Date("2019-01-01"), length = 12, by = "months") - 1 

# anytime we run this script before the start of the challenge we want to exclude all of 2019
if(day.run < ymd("2021-03-31")){
  end.date <- "2018-12-31"
} else { 
  month.run <- month(day.run)
  end.date <- date.end.month[month(day.run)]
}

make_dataset <- function(table, end.date, QF.col, min.col, max.col, var.col, exp.unc.col, 
                         new.min.var, new.min.exp.unc, new.max.var, new.max.exp.unc, QF.flag = 0){
  
  cols.select <- c("Date", "Year", "epiWeek", "yearWeek", "domainID", "siteID",
                    min.col, max.col, var.col, exp.unc.col)
  
  cols.select <- cols.select[!is.na(cols.select)] # remove NAs
  
  df <- table %>% 
    filter(.data[[QF.col]] == QF.flag) %>% # remove observations that fail QF
    filter(endDateTime <= end.date) %>%  # work with correct data
    mutate(Date = date(endDateTime), # get date
           Year = year(Date),    # add year column
           epiWeek = sprintf("%02d", epiweek(Date)), # add week column, leading 0 for single digits
           yearWeek = as.character(paste(Year, epiWeek, sep = "-"))) %>% 
    select(all_of(cols.select))
  
  weekly.min <- df %>% 
    group_by(siteID, yearWeek) %>%
    slice(which.min(.data[[min.col]])) %>% # use slice to retain variance column
    select(-c(Date, .data[[max.col]])) %>% 
    rename_at(vars(var.col), ~ new.min.var) %>% 
    rename_at(vars(exp.unc.col), ~ new.min.exp.unc)
  
  weekly.max <- df %>% 
    group_by(siteID, yearWeek) %>%
    slice(which.max(.data[[max.col]])) %>% # use slice to retain variance column
    select(-c(Date, .data[[min.col]])) %>% 
    rename_at(vars(var.col), ~ new.max.var) %>% 
    rename_at(vars(exp.unc.col), ~ new.max.exp.unc)
  
  weekly.data <- left_join(weekly.max, weekly.min)
  
  daily.min <- df %>% 
    group_by(siteID, Date) %>%
    slice(which.min(.data[[min.col]])) %>% # use slice to retain variance column
    select(-c(yearWeek, .data[[max.col]])) %>% 
    rename_at(vars(var.col), ~ new.min.var) %>% 
    rename_at(vars(exp.unc.col), ~ new.min.exp.unc)
  
  daily.max <- df %>% 
    group_by(siteID, Date) %>%
    slice(which.max(.data[[max.col]])) %>% # use slice to retain variance column
    select(-c(yearWeek, .data[[min.col]])) %>% 
    rename_at(vars(var.col), ~ new.max.var) %>% 
    rename_at(vars(exp.unc.col), ~ new.max.exp.unc)
  
  daily.data <- left_join(daily.max, daily.min)
  
  return(list(weekly.data = weekly.data,
              daily.data = daily.data))
  
}

# =========================== #
#      relative humidity      #
# =========================== #
if(make.rh){
  message("Extracting relative humidity")
  rh.data <- neon_read("RH_30min-expanded")
  relative.humidity <- make_dataset(
    table = rh.data,
    end.date = end.date,
    QF.col = "RHFinalQF",
    min.col = "RHMinimum",
    max.col = "RHMaximum",
    var.col = "RHVariance",
    exp.unc.col = "RHExpUncert",
    new.min.var = "RHMinimumVariance",
    new.min.exp.unc = "RHMinimumExpUncert",
    new.max.var = "RHMaximumVariance",
    new.max.exp.unc = "RHMaximumExpUncert"
  )  
  
  message("Writing relative humidity CSVs")
  write_csv(relative.humidity$weekly.data, "Data/RelativeHumidityWeekly.csv")
  write_csv(relative.humidity$daily.data, "Data/RelativeHumidityDaily.csv")

}

# =========================== #
#       Air temperature       #
# =========================== #
if(make.taat){
  message("Extracting air temperature")
  air.temp.data <- neon_read("TAAT_30min-expanded")
  air.temperature <- make_dataset(
    table = air.temp.data,
    end.date = end.date,
    QF.col = "finalQF",
    min.col = "tempTripleMinimum",
    max.col = "tempTripleMaximum",
    var.col = "tempTripleVariance",
    exp.unc.col = "tempTripleExpUncert",
    new.min.var = "airTempMinimumVariance",
    new.min.exp.unc = "airTempMinimumExpUncert",
    new.max.var = "airTempMaximumVariance",
    new.max.exp.unc = "airTempMaximumExpUncert"
  )
  
  message("Writing air temperature CSVs")
  write_csv(air.temperature$weekly.data, "Data/AirTemperatureWeekly.csv")
  write_csv(air.temperature$daily.data, "Data/AirTemperatureDaily.csv")
  
}

# =========================== #
#  IR Biological Temperature  #
# =========================== #
if(make.ir.temp){
  message("Extracting IR biological temperature")
  ir.temp.data <- neon_read("IRBT_30_minute-expanded")
  bio.temperature <- make_dataset(
    table = ir.temp.data,
    end.date = end.date,
    QF.col = "finalQF",
    min.col = "bioTempMinimum",
    max.col = "bioTempMaximum",
    var.col = "bioTempVariance",
    exp.unc.col = "bioTempExpUncert",
    new.min.var = "bioTempMinimumVariance",
    new.min.exp.unc = "bioTempMinimumExpUncert",
    new.max.var = "bioTempMaximumVariance",
    new.max.exp.unc = "bioTempMaximumExpUncert"
  )
  
  message("Writing IR biological temperature CSVs")
  write_csv(bio.temperature$weekly.data, "Data/BioTemperatureWeekly.csv")
  write_csv(bio.temperature$daily.data, "Data/BioTemperatureDaily.csv")  

}

# =========================== #
#      Soil Temperature       #
# =========================== #
if(make.soil.temp){
  message("Extracting soil temperature")
  soil.temp.data <- neon_read("ST_30_minute-expanded")
  soil.temperature <- make_dataset(
    table = soil.temp.data,
    end.date = end.date,
    QF.col = "finalQF",
    min.col = "soilTempMinimum",
    max.col = "soilTempMaximum",
    var.col = "soilTempVariance",
    exp.unc.col = "soilTempExpUncert",
    new.min.var = "soilTempMinimumVariance",
    new.min.exp.unc = "soilTempMinimumExpUncert",
    new.max.var = "soilTempMaximumVariance",
    new.max.exp.unc = "soilTempMaximumExpUncert"
  )
  
  message("Writing soil temperature CSVs")
  write_csv(soil.temperature$weekly.data, "Data/SoilTemperatureWeekly.csv")
  write_csv(soil.temperature$daily.data, "Data/SoilTemperatureDaily.csv")
  
}

# =========================== #
#     Soil Water Content      #
# =========================== #
if(make.soil.water){
  message("Extracting soil water content")
  soil.water.data <- neon_read("SWS_30_minute-expanded")
  soil.water <- make_dataset(
    table = soil.water.data,
    end.date = end.date,
    QF.col = "VSWCFinalQF",
    min.col = "VSWCMinimum",
    max.col = "VSWCMaximum",
    var.col = "VSWCVariance",
    exp.unc.col = "VSWCExpUncert",
    new.min.var = "VSWCMinimumVariance",
    new.min.exp.unc = "VSWCMinimumExpUncert",
    new.max.var = "VSWCMaximumVariance",
    new.max.exp.unc = "VSWCMaximumExpUncert"
  )
  
  message("Writing soil water content CSVs")
  write_csv(soil.water$weekly.data, "Data/SoilWaterContentWeekly.csv")
  write_csv(soil.water$daily.data, "Data/SoilWaterContentDaily.csv")
  
}



message("--- DONE ---")



