# ====================================================== #
#        plot effect for plots within sites              #
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

# get unique sites
sites.ixodes <- unique(data.ixodes$siteID)
sites.amblyomma <- unique(data.ambloyoma$siteID)

# restructure data to have plots combined by site
ls.ixodes <- list()
for(i in seq_along(sites.ixodes)){
  subset.df <- data.ixodes %>% 
    filter(siteID == sites.ixodes[i]) %>% 
    filter(Year >= 2016 & Year <= 2018) %>% 
    select(-c(decimalLatitude,
              decimalLongitude,
              elevation)) %>% 
    mutate(density = Ixodes_scapularis / totalSampledArea * 1600) %>%
    pivot_wider(id_cols = c(yearWeek, time), 
                names_from = plotID, 
                values_from = density) %>% 
    arrange(yearWeek)
  
  if(ncol(subset.df) == 3) next
  
  ls.ixodes[[i]] <- subset.df
}

ls.amblyomma <- list()
for(i in seq_along(sites.amblyomma)){
  
  subset.df <- data.ambloyoma %>% 
    filter(siteID == sites.amblyomma[i]) %>% 
    filter(Year >= 2016 & Year <= 2018) %>% 
    select(-c(decimalLatitude,
              decimalLongitude,
              elevation)) %>% 
    mutate(density = Ambloyomma_americanum / totalSampledArea * 1600) %>% 
    pivot_wider(id_cols = c(yearWeek, time), 
                names_from = plotID, 
                values_fn = {max},
                values_from = density) %>% 
    arrange(yearWeek)
  if(ncol(subset.df) == 3) next
  ls.amblyomma[[i]] <- subset.df
}

# name list elemetns by species_site
ls.ixodes <- set_names(ls.ixodes, 
                       paste("Ixodes_scapularis", sites.ixodes, sep = "_"))
ls.amblyomma <- set_names(ls.amblyomma, 
                          paste("Amblyomma_americanum", sites.amblyomma, sep = "_"))

# combine
data.list.master <- prepend(ls.ixodes, ls.amblyomma)
data.list.master <- purrr::compact(data.list.master)

# par(mfrow = c(2, 5))
# for(jj in seq_along(data.list.master)){
#   df <- pluck(data.list.master, jj)
#   y <- select(df, matches("\\d{3}"))
#   hist(log(y[!is.na(y)]))
# }

# source("Functions/ZIP.R") 
model <- nimbleCode({

  # priors
  sd.process ~ dgamma(10, 1) # process error
  sd.obs ~ dgamma(10, 1) # observation error
  sd.plot ~ dgamma(10, 1) # plot error
  mu ~ dnorm(0, sd = 5)

  # data model, density so using dnorm
  for(s in 1:n.sites){
    for(t in 1:n.weeks){
      y[t, s] ~ dnorm(z[t, s], sd = sd.obs)
    }
  }

  # process model
  for(s in 1:n.sites){

    # plot effect prior
    plot.effect[s] ~ dnorm(0, sd = sd.plot)

    # first state gets it's own prior
    z[1, s] ~ dnorm(x.ic[s], sd = 10)

    for(t in 2:n.weeks){
      ex[t-1, s] <- mu + plot.effect[s]
      x[t, s] ~ dnorm(ex[t-1, s], sd = sd.process)
      z[t, s] <- max(0, x[t, s])
    }
  }
})
monitor <- c("sd.process", "sd.obs", "sd.plot", "z", "plot.effect", "mu")

# model <- nimbleCode({
#   
#   # priors
#   sd.process ~ dgamma(10, 1) # process error
#   # pi ~ dunif(0, 1) # probability of structural zero 
#   sd.plot ~ dgamma(10, 1) # plot error 
#   
#   # data model, counts so poisson
#   for(s in 1:n.sites){
#     for(t in 1:n.weeks){
#       y[t, s] ~ dpois(z[t, s])
#     }  
#   }
#   
#   # process model
#   for(s in 1:n.sites){
#     
#     # plot effect prior
#     plot.effect[s] ~ dnorm(0, sd = sd.plot)
#     
#     # first state gets it's own prior
#     z[1, s] ~ dpois(x.ic[s])
#     
#     for(t in 2:n.weeks){
#       mu[t-1, s] <- z[t-1, s] + plot.effect[s]
#       ex[t, s] ~ dnorm(mu[t-1, s], sd = sd.process)
#       z[t, s] <- max(0, ex[t, s])
#     }
#   }
# })
# monitor <- c("sd.process", "sd.plot", "z", "plot.effect")

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
ForecastProject.id <- "plotEffect"     # Some ID that applies to a set of forecasts

# initialize array and list for storage
fit.ls <- group.names <- list()
fx.array <- array(data = NA, 
                  dim = c(n.ens,   # the number of ensembles we want to save      
                          n.weeks,   # the number of weeks forecasted
                          obs.dim,   # observation dimension
                          n.targets)) # the number of targets (species_site)


for(i in seq_along(data.list.master)){
  df <- pluck(data.list.master, i)
  
  y <- select(df, matches("\\d{3}"))
  pi.init <- sum(is.na(y)) / (nrow(y) * ncol(y))
  
  y.means <- colMeans(y, na.rm = TRUE)
  y.sd <- apply(y, 2, sd, na.rm = TRUE)
  mu.inits <- y
  for(cc in 1:ncol(y)){
    na.rows <- which(is.na(mu.inits[,cc]))
    mu.inits[na.rows, cc] <- round(approx(pull(mu.inits, cc), xout = na.rows)$y)
  }
  
  # n.weeks the number of weeks to forecast into the future
  # need to pad df to make the forecast in jags
  # na.pad <- matrix(NA, length(target.weeks), ncol(y))
  # colnames(na.pad) <- colnames(y)
  # group.names[[i]] <- colnames(y)
  # y <- rbind(y, na.pad)
  
  constants <- list(n.weeks = nrow(y),
                    n.sites = ncol(y),
                    x.ic = c(1 + rpois(ncol(y), 5), 0))
  
  data <- list(y = y)
  
  # lets say 2019 was like 2018 for inits
  # need to add more explicit matching for when challenge starts
  # na.fill <- mu.inits[(nrow(mu.inits)-length(target.weeks)+1):nrow(mu.inits),]
  # mu.inits <- rbind(mu.inits, na.fill)
  
  inits <- function(){list(y = mu.inits,
                           z = abs(apply(mu.inits, 2, jitter)) + 1,
                           # tau.process = rgamma(1, 1, 3),
                           # pi = pi.init,
                           # tau.plot = rgamma(1, 1, 3),
                           plot.effect = rnorm(ncol(data$y), 0, 1)
                           )}
  
  init <- inits()
  model.rw <- nimbleModel(model,
                          constants = constants,
                          data = data,
                          inits = init)
  cModel.rw <- compileNimble(model.rw)
  mcmcConf <- configureMCMC(cModel.rw, monitors = monitor)
  mcmcBuild <- buildMCMC(mcmcConf)
  mcmcBuild$run(1)
  # compMCMC <- compileNimble(mcmcBuild)
  # out.1 <- runMCMC(compMCMC, niter = 100000)
  
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
  print(apply(samples[,1:8], 2, quantile, c(0.025, 0.5, 0.975)))
  
  # extract forecast period 
  # forecast.index <- (nrow(df)+1):nrow(y) 
  # forecast.cols <- paste0("z[", forecast.index)
  # forecast.samps <- samples %>% 
  #   as_tibble() %>% 
  #   select(starts_with(forecast.cols))
  
  # columns that represent the forecast
  
  # fx.array[ , , obs.dim, i] <- states[thin.seq, forecast.cols] # store
  
  # mcmc.end <- Sys.time()
  # total.time <- mcmc.end - mcmc.start
  # message(paste("End time:", mcmc.end))
  # print(total.time)
}

fit.ls <- set_names(fit.ls, names(data.list.master))

# save entire model fit (params, calibration, fx)
save(fit.ls, group.names,
     file = here("ModelOut", paste0(ForecastProject.id, ".RData")))


