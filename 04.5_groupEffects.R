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
species.effect <- FALSE
site.effect <- FALSE
year.effect <- FALSE

# n.cores <- as.numeric(Sys.getenv("NSLOTS")) # number of cores to use for cluster 
n.cores <- 3

source("Models/yearTimeSeriesMetDiff.R")
group.dir <- "singleObsError"

# for getting met variable
csv <- "Data/BioTemperatureWeekly.csv"
met.var <- "bioTempMaximum"
met.uncertainty <- "bioTempMaximumVariance"

ForecastProject.id <- "SppSiteSurv"     # Some ID that applies to a set of forecasts
out.dir <- here("ModelOut", group.dir)
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
year.weeks <- tibble(yearWeek = wide.data$yearWeek)

# convert groupings to numeric vectors for indexing
year.index <- year(wide.data$time) - min(year(wide.data$time)) + 1
species.index <- as.factor(species.vec) %>% as.numeric()
plot.index <- as.factor(plot.vec) %>% as.numeric()
site.index <- as.factor(site.vec) %>% as.numeric()

# number of each group
n.species <- length(unique(species.index))
n.sites <- length(unique(site.index))
n.years <- length(unique(year.index))
years <- year(wide.data$time) %>% unique()

# put data in an array by year
y.table <- wide.data[,-c(1,2)]
weeks.per.year <- year.index %>% table()
max.weeks <-  max(weeks.per.year)
y <- array(NA, dim = c(3, # years
                       max.weeks, # weeks
                       ncol(y.table))) # site


source("Functions/get_met_array.R")
met.list <- get_met_array(csv, 
                          weeks.per.year,
                          year.weeks,
                          met.var,
                          met.uncertainty)
met.vals <- met.list$met.vals
x.var <- met.list$met.uncertainty
met.inits <- met.list$met.inits

# pull out driver variance, check for NAs, set to mean
if(any(is.na(x.var))){
  x.var[is.na(x.var)] <- mean(x.var, na.rm = TRUE)
}
x.var[which(x.var == 0)] <- 0.1

# scaling met by standard deviation by site
met.inits <- met.inits / sd(met.inits, na.rm = TRUE)
for(ss in 1:n.sites){
  site.sd <- sd(met.vals[,,ss], na.rm = TRUE)  
  met.vals[,,ss] <- met.vals[,,ss] / site.sd
  
  # variance is scaled by the square of the scale (sd)
  x.var[,,ss] <- x.var[,,ss] / site.sd^2
}



for(i in 1:3){
  y.sub <- wide.data %>% 
    filter(year(time) == years[i]) 
  
  y[i, 1:weeks.per.year[i], 1:ncol(y.table)] <- as.matrix(y.sub[,-c(1,2)]) 
}

data <- list(
  y = y,
  life.constraint = array(1, dim = dim(y)),
  x.obs = met.vals,
  x.tau = 1 / x.var
)

mu.inits <- y
for(i in 1:3){
  for(cc in 1:ncol(y.table)){
    na.rows <- which(is.na(mu.inits[i,,cc]))
    mu.inits[, na.rows, cc] <- approx(as.vector(mu.inits[i,,cc]), xout = na.rows)$y
  }  
}

mu.inits[is.na(mu.inits)] <- 0

total.mu.index <- dim(mu.inits)[1]*dim(mu.inits)[2]*dim(mu.inits)[3]

x.ic <- matrix(runif(n.years*ncol(y), 0, 1), n.years, ncol(y))

constants <- list(
  n.weeks = weeks.per.year,
  n.plots = ncol(y.table),
  x.ic = x.ic, #runif(ncol(y), 0, 1),
  met.mu = 0,
  tau.met = 1 / var(met.list$met.vals, na.rm = TRUE),
  n.sites = n.sites,
  n.species = n.species,
  species = species.index,
  n.years = n.years,
  site = site.index
)

# monitor originally defined in model script sourced at top
if(site.effect) monitor <- c(monitor, "tau.site", "alpha.site")
if(species.effect){
  monitor <- c(monitor, "tau.species", "alpha.species")
  constants$n.species <- n.species
  constants$species <- species.index
} 
if(year.effect){
  monitor <- c(monitor, "tau.year", "alpha.year")
  constants$year <- year.index
} 



model.rw <- nimbleModel(model, 
                        constants = constants,
                        data = data,
                        inits = inits())
model.rw$initializeInfo()
cModel.rw <- compileNimble(model.rw)
mcmcConf <- configureMCMC(cModel.rw, monitors = monitor)
# mcmcConf$addSampler(target = c("tau.obs[1]", "tau.obs[2]"), type = "RW_block")
mcmcBuild <- buildMCMC(mcmcConf)
mcmcBuild$run(1)

mcmc.start <- Sys.time()
message(paste("Start time:", mcmc.start))

file.name <- ForecastProject.id
if(species.effect) file.name <- paste0(file.name, "SpeciesEffect")
if(site.effect) file.name <- paste0(file.name, "SiteEffect")
if(year.effect) file.name <- paste0(file.name, "YearEffect")
rdata.name <- paste0(file.name, ".RData")

save.path <- here("ModelOut",
                  group.dir,
                  rdata.name)

message(paste0("Writing model output to ", save.path))

source("Functions/run_nimble_parallel.R")
cl <- makeCluster(n.cores)
samples <- run_nimble_parallel(
    cl,
    model,
    constants,
    data,
    inits,
    monitor,
    n.iter = 30000,
    max.iter = 5e6,
    thin = 10,
    check.interval = 5,
    check.params.only = TRUE,
    file.name = save.path
  )
stopCluster(cl)  




save(samples, file = save.path)





