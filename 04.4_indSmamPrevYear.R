# ====================================================== #
#        for fitting models to sites indepently          #
#        Using a normal-normal state-space model         #
#        Using tick density per trapping (1600^m-2)      #
# ====================================================== #


script.start <- Sys.time()
renv::restore()

library(nimble) # for building models and mcmc
library(coda) # for mcmc stats
library(tidyverse) # for data wrangling and piping (dplyr probably ok)
library(lubridate) # for finding year from dates
library(stringr) # for searching within character strings 
library(here) # for working within subdirectories
library(parallel) # for running models in parallel
library(ISOweek) # to convert year-weeks to date

check.models <- TRUE

## load model
ForecastProject.id <- "SmamPrevMonth"     # Some ID that applies to a set of forecasts
out.dir <- here("ModelOut", ForecastProject.id)
if(!dir.exists(out.dir)) dir.create(out.dir)

# names of models already run - some may not have converged!
reduce.models <- FALSE
if(check.models){
  models.done <- list.files(out.dir)
  for(m in seq_along(models.done)){
    load(file.path(out.dir, models.done[m]))
    if(!samples$convergence) models.done[m] <- NA
  }
  completed.fits <- models.done[!is.na(models.done)] # models that have converged
  completed.fits <- gsub(".RData", "", completed.fits)
  reduce.models <- TRUE
}

completed.fits <- completed.fits[which(completed.fits!="Amblyomma_americanum_ORNL_007")]

smam.data <- read_csv("Data/smallMammalAnimalsByBout.csv")
smam.data <- smam.data %>% 
  mutate(Year = collectYear,
         month = collectMonth) %>% 
  select(c(month,
           Year,
           siteID,
           totalPeromyscus,
           totalAnimals))

# split data by species
source("Functions/split_species.R")
data.ixodes <- split_species("Ixodes")
data.amblyomma <- split_species("Amblyomma")

data.ixodes <- data.ixodes %>% 
  mutate(month = month(time))
data.amblyomma <- data.amblyomma %>% 
  mutate(month = month(time))

data.ixodes <- left_join(data.ixodes, smam.data,
                         by = c("siteID", "Year", "month"))
data.amblyomma <- left_join(data.amblyomma, smam.data,
                         by = c("siteID", "Year", "month"))

# get unique plots
plot.ixodes <- unique(data.ixodes$plotID)
plot.amblyomma <- unique(data.amblyomma$plotID)

# restructure data to have plots combined by site
ls.ixodes <- list()
for(i in seq_along(plot.ixodes)){
  subset.df <- data.ixodes %>% 
    filter(plotID == plot.ixodes[i]) %>% 
    filter(Year >= 2016 & Year <= 2018) %>% 
    mutate(density = Ixodes_scapularis / totalSampledArea * 1600,
           yearEffect = Year - 2015,
           individuals = Ixodes_scapularis) %>%
    arrange(yearWeek)
  
  ls.ixodes[[i]] <- subset.df
}

ls.amblyomma <- list()
for(i in seq_along(plot.amblyomma)){
  subset.df <- data.amblyomma %>% 
    filter(plotID == plot.amblyomma[i]) %>% 
    filter(Year >= 2016 & Year <= 2018) %>% 
    mutate(density = Ambloyomma_americanum / totalSampledArea * 1600,
           yearEffect = Year - 2015,
           individuals = Ambloyomma_americanum) %>% 
    arrange(yearWeek)
  
  if(plot.amblyomma[i] == "UKFS_003"){
    subset.df <- subset.df %>% 
      filter(!(yearWeek == 201635 & Ambloyomma_americanum == 0))
  }
  
  ls.amblyomma[[i]] <- subset.df
}

# name list elemetns by species_site
ls.ixodes <- set_names(ls.ixodes, 
                       paste("Ixodes_scapularis", plot.ixodes, sep = "_"))
ls.amblyomma <- set_names(ls.amblyomma, 
                          paste("Amblyomma_americanum", plot.amblyomma, sep = "_"))

# combine
data.list.master <- prepend(ls.ixodes, ls.amblyomma)

source("Functions/target_weeks.R")
day.run <- today()
target.weeks <- target_weeks(day.run)
n.weeks <- length(target.weeks)       # the number of weeks forecasted, same across all targets
obs.dim <- 1                          # 1 = latent state; 2 = latent state + observation error
n.ens <- 500                          # how many samples do we want to save? 
n.targets <- length(data.list.master) # the number of unique species_plot combinations 
targetName <- names(data.list.master) # the names of specific targets
forecast <- rep(0, n.weeks)           # is this a hindcast or forecast? hindcast = 0
n.cores <- as.numeric(Sys.getenv("NSLOTS")) # number of cores to use for cluster 

# forecast identifiers used in EFI standards
forecast.issue.time <- paste0(2019,"-", target.weeks[1])    # start date of the forecast YYYYWW
Forecast.id <- day.run    # ISO datetime should make a valid Forecast_id

# initialize array and list for storage
fit.ls <- group.names <- list()
converged <- rep(NA, length(data.list.master))
fx.array <- array(data = NA, 
                  dim = c(n.ens,   # the number of ensembles we want to save      
                          n.weeks,   # the number of weeks forecasted
                          obs.dim,   # observation dimension
                          n.targets)) # the number of targets (species_site)

model <- nimbleCode({
  
  # priors
  sd.process ~ dgamma(10, 1) # process error
  sd.obs ~ dgamma(10, 1) # observation error
  beta ~ dnorm(0, sd = 5) # intercept
  rho ~ dnorm(0, sd = 5) # slope
  sd.x ~ dgamma(10, 1) # met error
  
  # data model, density so using dnorm
  for(t in 1:n.weeks){
    y[t] ~ dnorm(z[t], sd = sd.obs) # tick data model
    
    # driver data model
    x[t] ~ dnorm(x.mu, sd = x.sd) # driver prior
    x.obs[t] ~ dnorm(x[t], sd = sd.x)
  }
  
  # process model
  # z[1] ~ dnorm(x.ic, sd = 10) # first state gets it's own prior
  for(t in 2:n.weeks){
    
    ex[t-1] <- beta*(x[t-1]) + rho*z[t-1] 
    m[t] ~ dnorm(ex[t-1], sd = sd.process)
    z[t] <- max(0, m[t])
  }
})

monitor <- c("sd.process", "sd.obs", "beta", "rho", "z", "sd.x")

inits <- function() {
  
  list(
    z = pmax(0, jitter(mu.inits)) + 1,
    y = jitter(mu.inits),
    ex = jitter(mu.inits),
    m = jitter(mu.inits),
    x = jitter(x.inits),
    x.obs = jitter(x.inits),
    sd.process = runif(1, 4, 30), 
    sd.obs = rgamma(1, 10, 1),
    sd.x = rgamma(1, 15, 1),
    beta = rnorm(1, -3, 3),
    rho = rnorm(1, 1, 0.1)
  )
}

# discard data sets already fit
if(reduce.models){
  data.list.master <- data.list.master %>% 
    discard(names(data.list.master) %in% completed.fits)
}

# want mice caught in the previous month
for(i in seq_along(data.list.master)){
  
  message("=================================================")
  message(paste("Fitting", names(data.list.master[i])))
  message(paste("Model", i, "of", length(data.list.master)))
  
  df <- pluck(data.list.master, i)
  
  y <- df$density
  
  months <- df %>% 
    select(c(Year, month, density, totalAnimals))
  
  prev.months <- months %>% 
    select(c(Year, month)) %>% 
    mutate(month = month - 1) %>% 
    mutate(Year = if_else(month == 0, Year - 1, Year))
  
  smam.months <- df %>% 
    group_by(Year, month) %>% 
    select(Year, month, totalAnimals) %>% 
    distinct()
  
  smam.x <- left_join(prev.months, smam.months, by = c("Year", "month"))
  x.obs <- smam.x$totalAnimals
  x.obs <- x.obs - mean(x.obs, na.rm = TRUE)
  
  # make inits and check for NAs
  x.inits <- x.obs
  if(any(is.na(x.obs))){
    x.inits[is.na(x.inits)] <- approx(x.inits, xout = which(is.na(x.inits)))$y
  }
  
  constants <- list(n.weeks = length(y))
    
  data <- list(
    y = y,
    x.obs = x.obs,
    x.mu = 0, # driver data is centered 
    x.sd = sd(x.obs, na.rm = TRUE)
  )

  mu.inits <- y
  mu.inits[is.na(mu.inits)] <- round(approx(mu.inits, xout = which(is.na(mu.inits)))$y)
  
  init <- inits()
  # print(init)
  model.rw <- nimbleModel(model,
                          constants = constants,
                          data = data,
                          inits = init)
  model.rw$initializeInfo()
  cModel.rw <- compileNimble(model.rw)
  mcmcConf <- configureMCMC(cModel.rw, monitors = monitor)
  mcmcBuild <- buildMCMC(mcmcConf)
  mcmcBuild$run(1)
  # compMCMC <- compileNimble(mcmcBuild)
  # out.1 <- runMCMC(compMCMC, niter = 10000)
  
  mcmc.start <- Sys.time()
  message(paste("Start time:", mcmc.start))
  
  source("Functions/run_nimble_parallel.R")
  cl <- makeCluster(n.cores)
  samples <- run_nimble_parallel(cl, model, constants, data, inits, monitor, n.iter = 100000,
                                 use.dzip = FALSE, check.params.only = TRUE)
  stopCluster(cl)  
  
  save(samples,
       file = here("ModelOut", 
                   ForecastProject.id, 
                   paste0(names(data.list.master[i]), ".RData")))
  
  
  out <- samples$samples
  converged[i] <- samples$convergence 
  
  print(apply(out[,1:5], 2, quantile, c(0.025, 0.5, 0.975)))
  
}
