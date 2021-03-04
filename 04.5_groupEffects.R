# ====================================================== #
#        group effects (spp, site, year).                #
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

# where do we want random effects?
species.effect <- TRUE
site.effect <- FALSE
year.effect <- FALSE

# n.cores <- as.numeric(Sys.getenv("NSLOTS")) # number of cores to use for cluster 
n.cores <- 3

ForecastProject.id <- "groupEffects"     # Some ID that applies to a set of forecasts
out.dir <- here("ModelOut", ForecastProject.id)
if(!dir.exists(out.dir)) dir.create(out.dir)

# split data by species
source("Functions/split_species.R")
data.ixodes <- split_species("Ixodes")
data.ambloyoma <- split_species("Amblyomma")

# get unique sites
sites.ixodes <- unique(data.ixodes$siteID)
sites.amblyomma <- unique(data.ambloyoma$siteID)

# get unique plots
plots.ixodes <- unique(data.ixodes$plotID)
plots.amblyomma <- unique(data.ambloyoma$plotID)

wide.ixodes <- data.ixodes %>% 
  filter(plotID %in% plots.ixodes) %>% 
  filter(Year >= 2016 & Year <= 2018) %>%
  mutate(density = Ixodes_scapularis / totalSampledArea * 1600) %>%
  pivot_wider(id_cols = c(yearWeek, time), 
              names_from = plotID, 
              values_from = density,
              names_glue = "Ixodes_{plotID}") 

wide.amblyomma <- data.ambloyoma %>% 
  filter(plotID %in% plots.amblyomma) %>% 
  filter(Year >= 2016 & Year <= 2018) %>%
  mutate(density = Ambloyomma_americanum / totalSampledArea * 1600) %>% 
  pivot_wider(id_cols = c(yearWeek, time), 
              names_from = plotID, 
              values_fn = {max},
              values_from = density,
              names_glue = "Amblyomma_{plotID}")

wide.data <- full_join(wide.ixodes, wide.amblyomma,
                       by = c("yearWeek", "time")) %>% 
  arrange(yearWeek)

# extract groupings (species, plots, sites)
ids <- colnames(wide.data[-c(1, 2)])
species.vec <- c(rep("ixodes", length(plots.ixodes)), rep("amblyomma", length(plots.amblyomma)))
plot.vec <- stringr::str_extract(ids, "[[:upper:]]{4}_\\d{3}")
site.vec <- stringr::str_extract(ids, "[[:upper:]]{4}")

# convert grouings to numeriv vectors
year.vec <- year(wide.data$time) - min(year(wide.data$time)) + 1
species.vec <- as.factor(species.vec) %>% as.numeric()
plot.vec <- as.factor(plot.vec) %>% as.numeric()
site.vec <- as.factor(site.vec) %>% as.numeric()

model <- nimbleCode({
  
  # priors
  sd.process ~ dgamma(10, 1) # process error
  sd.obs ~ dgamma(10, 1) # observation error
  mu ~ dnorm(0, sd = 100) # global mean 
  
  if(species.effect){
    sd.species ~ dgamma(10, 1)
    for(spp in 1:n.species){
      alpha.species[spp] ~ dnorm(0, sd = sd.species)  
    }  
  } else if(site.effect){
    sd.site ~ dgamma(10, 1)
    for(se in 1:n.sites){
      alpha.site[se] ~ dnorm(0, sd = sd.site)
    }
  } else if(year.effect){
    sd.year ~ dgamma(10, 1)
    for(ye in 1:n.years){
      alpha.year[ye] ~ dnorm(0, sd = sd.year)
    }
  }
  
  
  # data model, density so using dnorm
  for(p in 1:n.plots){
    for(t in 1:n.weeks){
      y[t, p] ~ dnorm(z[t, p], sd = sd.obs)
    }
  }
  
  # process model
  for(p in 1:n.plots){
    
    # first state gets it's own prior
    x[1, p] ~ dnorm(x.ic[p], sd = 10)
    z[1, p] <- max(0, x[1, p])
    
    for(t in 2:n.weeks){
      if(species.effect){
        ex[t, p] <- mu + alpha.species[species[p]]
      } else if(site.effect){
        ex[t, p] <- mu + alpha.site[site[p]]
      } else if(year.effect){
        ex[t, p] <- mu + alpha.year[year[t]]
      }
      
      x[t, p] ~ dnorm(ex[t, p], sd = sd.process)
      z[t, p] <- max(0, x[t, p])
    
    } # time
  } # plots
})

y <- wide.data[,-c(1,2)]
data <- list(y = y)

n.species <- length(unique(species.vec))
n.sites <- length(unique(site.vec))
n.years <- length(unique(year.vec))

mu.inits <- y
for(cc in 1:ncol(y)){
  na.rows <- which(is.na(mu.inits[,cc]))
  mu.inits[na.rows, cc] <- round(approx(pull(mu.inits, cc), xout = na.rows)$y)
}
mu.inits[is.na(mu.inits)] <- 0

inits <- function(){list(
  mu = jitter(mean(as.matrix(mu.inits))),
  y = apply(mu.inits, 2, jitter, pmax(0)),
  z = apply(mu.inits, 2, jitter, pmax(0)),
  x = apply(mu.inits, 2, jitter, pmax(0)),
  ex = apply(mu.inits, 2, jitter, pmax(0)),
  sd.process = runif(1, 4, 30), 
  sd.obs = rgamma(1, 10, 1),
  alpha.species = rnorm(n.species, 0, 1),
  alpha.site = rnorm(n.sites, 0, 1),
  alpha.year = rnorm(n.years, 0, 1),
  sd.species = rgamma(1, 15, 1)
)}

monitor <- c("sd.process", "sd.obs", "mu", "z")
if(species.effect) monitor <- c(monitor, "sd.species", "alpha.species")
if(site.effect) monitor <- c(monitor, "sd.site", "alpha.site")
if(year.effect) monitor <- c(monitor, "sd.year", "alpha.year")


constants <- list(
  n.weeks = nrow(y),
  n.plots = ncol(y),
  x.ic = runif(ncol(y), 0, 1),
  n.species = n.species,
  n.sites = n.sites,
  n.years = n.years,
  species = species.vec,
  site = site.vec,
  year = year.vec,
  species.effect = species.effect,
  site.effect = site.effect,
  year.effect = year.effect
)

model.rw <- nimbleModel(model, 
                        constants = constants,
                        data = data,
                        inits = inits())
model.rw$initializeInfo()
cModel.rw <- compileNimble(model.rw)
mcmcConf <- configureMCMC(cModel.rw, monitors = monitor)
mcmcBuild <- buildMCMC(mcmcConf)
mcmcBuild$run(1)

mcmc.start <- Sys.time()
message(paste("Start time:", mcmc.start))

file.name <- ForecastProject.id
if(species.effect) file.name <- paste0(file.name, "Species")
if(site.effect) file.name <- paste0(file.name, "Site")
if(year.effect) file.name <- paste0(file.name, "Year")
rdata.name <- paste0(file.name, ".RData")

source("Functions/run_nimble_parallel.R")
cl <- makeCluster(n.cores)
samples <- run_nimble_parallel(
    cl,
    model,
    constants,
    data,
    inits,
    monitor,
    n.iter = 25000,
    max.iter = 2e6,
    check.params.only = TRUE,
    file.name = rdata.name
  )
stopCluster(cl)  




save(samples,
     file = here("ModelOut", 
                 ForecastProject.id, 
                 rdata.name))





