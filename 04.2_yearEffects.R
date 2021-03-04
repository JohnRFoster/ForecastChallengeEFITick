# ====================================================== #
#        year effect for each species at each plot       #
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

# split data by species
source("Functions/split_species.R")
data.ixodes <- split_species("Ixodes")
data.ambloyoma <- split_species("Amblyomma")

# get unique plots
plot.ixodes <- unique(data.ixodes$plotID)
plot.amblyomma <- unique(data.ambloyoma$plotID)

# restructure data to have plots combined by site
ls.ixodes <- list()
for(i in seq_along(plot.ixodes)){
  subset.df <- data.ixodes %>% 
    filter(plotID == plot.ixodes[i]) %>% 
    filter(Year >= 2016 & Year <= 2018) %>% 
    select(c(Year,
             yearWeek,
             time,
             plotID,
             Ixodes_scapularis,
             totalSampledArea)) %>% 
    mutate(density = Ixodes_scapularis / totalSampledArea * 1600,
           yearEffect = Year - 2015,
           individuals = Ixodes_scapularis) %>%
    arrange(yearWeek)
  
  ls.ixodes[[i]] <- subset.df
}

ls.amblyomma <- list()
for(i in seq_along(plot.amblyomma)){
  subset.df <- data.ambloyoma %>% 
    filter(plotID == plot.amblyomma[i]) %>% 
    filter(Year >= 2016 & Year <= 2018) %>% 
    select(-c(decimalLatitude,
              decimalLongitude,
              elevation)) %>% 
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

model <- nimbleCode({

  # priors
  sd.process ~ dgamma(10, 1) # process error
  sd.obs ~ dgamma(10, 1) # observation error
  sd.year ~ dgamma(5, 1) # year-to-year error
  mu ~ dnorm(0, sd = 5) # global mean

  for(y in 1:n.years){ # year effect
    year.effect[y] ~ dnorm(0, sd = sd.year)
  }

  # data model, density so using dnorm
  for(t in 1:n.weeks){
      y[t] ~ dnorm(z[t], sd = sd.obs)
  }

  # process model
  # z[1] ~ dnorm(x.ic, sd = 10) # first state gets it's own prior
  for(t in 2:n.weeks){
    ex[t-1] <- mu + year.effect[yearEffect[t]]
    x[t] ~ dnorm(ex[t-1], sd = sd.process)
    z[t] <- max(0, x[t])
  }
})
monitor <- c("sd.process", "sd.obs", "sd.year", "z", "year.effect", "mu")

# model <- nimbleCode({
#   
#   # priors
#   sd.process ~ dgamma(10, 1) # process error
#   mu ~ dnorm(0, sd = 10) # overall mean
#   sd.year ~ dgamma(10, 1) # year-to-year error 
#   
#   for(y in 1:n.years){ # year effect
#     year.effect[y] ~ dnorm(0, sd = sd.year)
#   }
#   
#   # data model
#   for(t in 1:n.weeks){
#       y[t] ~ dpois(z[t]) #, zeroProb = pi)
#   }
#   
#   # process model
#   # z[1] ~ dpois(x.ic) # first state gets it's own prior
#   for(t in 2:n.weeks){
#     ex[t-1] <- mu + year.effect[yearEffect[t]]
#     x[t] ~ dnorm(ex[t-1], sd = sd.process)
#     z[t] <- max(0, x[t])
#   }
# })
# source("Functions/ZIP.R")
# monitor <- c("sd.process", "sd.year", "mu", "z", "year.effect")

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
ForecastProject.id <- "yearEffect"     # Some ID that applies to a set of forecasts

# initialize array and list for storage
fit.ls <- group.names <- list()
fx.array <- array(data = NA, 
                  dim = c(n.ens,   # the number of ensembles we want to save      
                          n.weeks,   # the number of weeks forecasted
                          obs.dim,   # observation dimension
                          n.targets)) # the number of targets (species_site)


for(i in seq_along(data.list.master)){
  df <- pluck(data.list.master, i)
  
  # y <- df$individuals
  # pi.init <- sum(is.na(y)) / length(y)
  y <- df$density
  yearEffect <- df$yearEffect
  
  # n.weeks the number of weeks to forecast into the future
  # need to pad df to make the forecast in jags
  # na.pad <- rep(NA, length(target.weeks))
  # y <- c(y, na.pad)
  # yearEffect <- c(yearEffect, rep(max(yearEffect)+1, length(na.pad)))
  
  constants <- list(n.weeks = length(y),
                    n.years = length(unique(yearEffect)),
                    yearEffect = yearEffect
                    # x.ic = rpois(1, 5)
                    )
  
  data <- list(y = y)
  
  mu.inits <- y
  mu.inits[is.na(mu.inits)] <- round(approx(mu.inits, xout = which(is.na(mu.inits)))$y)
  
  inits <- function(){list(mu = mean(y, na.rm = TRUE),
                           z = pmax(0, jitter(mu.inits)) + 1,
                           y = mu.inits,
                           sd.process = runif(1, 5, 10),
                           # pi = pi.init,
                           # tau.obs = rgamma(1, 1, 1),
                           sd.year = rgamma(1, 5, 1),
                           year.effect = rnorm(length(unique(yearEffect)), 0, 1)
                           )}
  
  init <- inits()
  # print(init)
  model.rw <- nimbleModel(model,
                          constants = constants,
                          data = data,
                          inits = init)
  cModel.rw <- compileNimble(model.rw)
  mcmcConf <- configureMCMC(cModel.rw, monitors = monitor)
  mcmcBuild <- buildMCMC(mcmcConf)
  mcmcBuild$run(1)
  # compMCMC <- compileNimble(mcmcBuild)
  # out.1 <- runMCMC(compMCMC, niter = 10000)
  
  mcmc.start <- Sys.time()
  message("=================================================")
  message(paste("Fitting", names(data.list.master[i])))
  message(paste("Model", i, "of", length(data.list.master)))
  message(paste("Start time:", mcmc.start))
  
  source("Functions/run_nimble_parallel.R")
  cl <- makeCluster(n.cores)
  samples <- run_nimble_parallel(cl, model, constants, data, inits, monitor,
                                 use.dzip = FALSE, check.params.only = TRUE)
  stopCluster(cl)  
  
  fit.ls[[i]] <- samples # save everything for later analysis
  
  print(apply(samples[,1:5], 2, quantile, c(0.025, 0.5, 0.975)))

}

fit.ls <- set_names(fit.ls, names(data.list.master))

# save entire model fit (params, calibration, fx)
save(fit.ls, group.names,
     file = here("ModelOut", paste0(ForecastProject.id, ".RData")))

