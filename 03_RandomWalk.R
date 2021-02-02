# ====================================================== #
#        Random Walk for tick nymphs in parallel         #
# ====================================================== #

renv::restore()

library(nimble) # for building models and mcmc
library(coda) # for mcmc stats
library(tidyverse) # for data wrangling and piping (dplyr probably ok)
library(lubridate) # for finding year from dates
library(stringr) # for searching within character strings 
library(here) # for working within subdirectories
library(parallel) # for running models in parallel
library(ISOweek) # to convert year-weeks to date

# first load the target data set
data <- read.csv("ticks-targets.csv.gz")

# for the random walk all we need are the targets and yearWeek
data <- data %>% 
  select(all_of(c("yearWeek", 
                  "targetSpecies", 
                  "targetCount", 
                  "targetPlotID", 
                  "specificTarget")))

# next, we need to extract the targets into their respective groups
# so we need data frames for each species x plot combination, retain
# the NA rows (weeks without observations), and make sure that the 
# plots that have both species present are separated into differnt
# target data sets

specific.targets <- data %>% 
  pull(specificTarget) %>% 
  unique()

data.list.master <- list()
for(i in seq_along(specific.targets)){
  data.list.master[[i]] <- data %>% 
    filter(specificTarget == specific.targets[i])
}

# names list elments
data.list.master <- set_names(data.list.master, specific.targets)

# set up target dates
start.epi.weeks <- c(10, 14, 19, 23, 28, 32, 36, 41) # all possible start weeks
day.run <- lubridate::today() # the day the script is called

# anytime we run this script before the start of the challenge we want to forecast all 2019 target weeks
if(day.run < "2021-03-31"){ 
  start.week <- start.epi.weeks[1]
} else { # otherwise use the appropriate starting week (months are 2 ahead)
  start.week <- start.epi.weeks[month(day.run) - 2]
}

end.week <- 44 # does not change
target.weeks <- start.week:end.week
n.weeks <- length(target.weeks)       # the number of weeks forecasted, same across all targets
obs.dim <- 1                          # 1 = latent state; 2 = latent state + observation error
n.ens <- 500                          # how many samples do we want to save? 
n.targets <- length(data.list.master) # the number of unique species_plot combinations 
targetName <- names(data.list.master) # the names of specific targets
forecast <- rep(0, n.weeks)           # is this a hindcast or forecast? hindcast = 0
n.cores <- as.numeric(Sys.getenv("NSLOTS")) # number of cores to use for cluster 

# forecast identifiers used in EFI standards
forecast.issue.time <- paste0(2019,"-", start.week)    # start date of the forecast YYYYWW
Forecast.id <- day.run    # ISO datetime should make a valid Forecast_id
ForecastProject.id <- "RandomWalk"     # Some ID that applies to a set of forecasts

# random walk code
random_walk <- nimbleCode({
  
  # process error prior
  tau.process ~ dgamma(1, 1)
  
  # first state gets it's own prior
  x[1] ~ dpois(x.ic)
  
  # data model, counts so Poisson
  for(t in 1:n){
    y[t] ~ dpois(x[t])
  }
  
  # process model
  for(c in 2:n){
    x[c] ~ T(dnorm(x[c-1], tau = tau.process), 0, Inf)
  }
})


# initialize array and list for storage
fit.ls <- list()
fx.array <- array(data = NA, 
                  dim = c(n.ens,   # the number of ensembles we want to save      
                          n.weeks,   # the number of weeks forecasted
                          obs.dim,   # observation dimension
                          n.targets)) # the number of targets (species_site)

for(i in seq_along(data.list.master)){
  df <- pluck(data.list.master, i)
  
  # check we are only modeling one species at one plot
  spp <- df %>% pull(targetSpecies) %>% unique() 
  spp <- spp[!is.na(spp)]
  if(length(spp) > 1) stop("Both species in data! \n", call. = FALSE)
  
  plot <- df %>% pull(targetPlotID) %>% unique()
  if(length(plot) > 1) stop("More than one plot in data! \n", call. = FALSE)
  
  # n.weeks the number of weeks to forecast into the future
  # need to pad df to make the forecast in jags
  na.pad <- rep(NA, length(target.weeks))
  y <- c(df$targetCount, na.pad)
  
  # mean counts to set initial condition prior
  ic <- mean(df$targetCount, na.rm = TRUE)
  
  constants <- list(n = length(y),
                    x.ic = rpois(1, 5))
  data <- list(y = y)
  inits <- function(){list(x = y + abs(round(rnorm(length(y), 0, 4))))}
  
  source("run_nimble_parallel.R")
  cl <- makeCluster(n.cores)
  out.par <- run_nimble_parallel(cl, random_walk, constants, data, inits)
  stopCluster(cl)  
  samples <- as.matrix(out.par) 
  fit.ls[[i]] <- samples # save everything for later analysis
  
  # extract forecast period and save in list
  states <- samples[,-grep("tau.process", colnames(samples))] # remove process error column
  
  # columns that represent the forecast
  forecast.cols <- (nrow(df)+1):ncol(states) 
  
  # we are going to save a thinned matrix (raw mcmc objects can get big)
  thin.seq <- round(seq(1, nrow(states), length.out = n.ens)) # sequence of samples to keep

  fx.array[ , , obs.dim, i] <- states[thin.seq, forecast.cols] # store
}

# get vector of species names, first extract Ixodes_scapularis
species.name <- str_extract(targetName, "Ixodes_scapularis")

# NAs are Ambloyomma_americanum
species.name[is.na(species.name)] <- "Ambloyomma_americanum"

# extract plotIDs
plot.id <- str_extract(targetName, "[[:upper:]]{4}_\\d{3}")

# convert week to date mapping to the first day of the week
date.col <- ISOweek2date(paste0("2019-W", target.weeks, "-1")) %>% as.character()

plot.id.unique <- plot.id %>% unique()
fx.df <- tibble()
for(p in seq_along(plot.id.unique)){
  plot.subset <- plot.id.unique[p]
  fx.index <- grep(plot.subset, targetName)
  
  # if only one species present at the plot
  if(length(fx.index) == 1){
    fx <- fx.array[,,obs.dim,fx.index]
    colnames(fx) <- date.col
    fx <- fx %>% 
      as_tibble() %>% 
      pivot_longer(all_of(date.col), 
                   names_to = "time",
                   values_to = species.name[fx.index]) %>% 
      mutate(plot = plot.subset,
             ensemble = rep(1:n.ens, each = length(target.weeks)),
             data_assimilation = 0,
             forecast = 0,
             obs_flag = obs.dim) 
    
    # if both species are present at the plot
  } else if (length(fx.index == 2)){
    
    fx.1 <- fx.array[,,obs.dim,fx.index[1]] 
    fx.2 <- fx.array[,,obs.dim,fx.index[2]] 
    colnames(fx.1) <- colnames(fx.2) <- date.col
    
    fx.1 <- fx.1 %>% 
      as_tibble() %>% 
      pivot_longer(all_of(date.col), 
                   names_to = "time",
                   values_to = species.name[fx.index[1]]) %>% 
      mutate(plot = plot.subset,
             ensemble = rep(1:n.ens, each = length(target.weeks)),
             data_assimilation = 0,
             forecast = 0,
             obs_flag = obs.dim)
    
    fx.2 <- fx.2 %>% 
      as_tibble() %>% 
      pivot_longer(all_of(date.col), 
                   names_to = "time",
                   values_to = species.name[fx.index[2]]) %>% 
      mutate(plot = plot.subset,
             ensemble = rep(1:n.ens, each = length(target.weeks)),
             data_assimilation = 0,
             forecast = 0,
             obs_flag = obs.dim)
    
    fx <- left_join(fx.1,
                    fx.2,
                    by = c("time", "plot", "ensemble", "data_assimilation", "forecast", "obs_flag"))
    
  }
  fx.df <- bind_rows(fx.df, fx)
}

# each forecast will get its own dir 
dir.ncfname <- file.path("Random_Walk_Fits", as.character(forecast.issue.time)) 

if(!dir.exists(dir.ncfname)) dir.create(dir.ncfname, recursive = TRUE)

# Save file as CSV in the
# [theme_name]-[yearWeek]-[team_name].csv
fx.file.name <- paste0("ticks-", 
                       as.character(date.col[1]), 
                       "-", 
                       ForecastProject.id, 
                       ".csv.gz")

write.csv(fx.df,
          file = file.path(dir.ncfname, fx.file.name))

# jags.csv <- read.csv("../RCN_tick_population/Random_Walk_Fits/201910/random-walk-forecast-summary-jags.csv")
# jags.fx <- jags.csv %>% 
#   filter(target == data.test$specificTarget[1])
# tail(jags.fx)


