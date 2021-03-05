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

source("Models/yearTimeSeriesBasic.R")

ForecastProject.id <- "yearTimeSeriesBasic"     # Some ID that applies to a set of forecasts
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

# convert groupings to numeric vectors
year.vec <- year(wide.data$time) - min(year(wide.data$time)) + 1
species.vec <- as.factor(species.vec) %>% as.numeric()
plot.vec <- as.factor(plot.vec) %>% as.numeric()
site.vec <- as.factor(site.vec) %>% as.numeric()

# number of each group
n.species <- length(unique(species.vec))
n.sites <- length(unique(site.vec))
n.years <- length(unique(year.vec))
years <- year(wide.data$time) %>% unique()

# put data in an array by year
y.table <- wide.data[,-c(1,2)]
weeks.per.year <- year.vec %>% table()
max.weeks <-  max(weeks.per.year)
y <- array(NA, dim = c(3, # years
                       max.weeks, # weeks
                       ncol(y.table))) # site

for(i in 1:3){
  y.sub <- wide.data %>% 
    filter(year(time) == years[i]) 
  
  y[i, 1:weeks.per.year[i], 1:ncol(y.table)] <- as.matrix(y.sub[,-c(1,2)]) 
}

data <- list(y = y)

mu.inits <- y
for(i in 1:3){
  for(cc in 1:ncol(y.table)){
    na.rows <- which(is.na(mu.inits[i,,cc]))
    mu.inits[, na.rows, cc] <- approx(as.vector(mu.inits[i,,cc]), xout = na.rows)$y
  }  
}

mu.inits[is.na(mu.inits)] <- 0

total.mu.index <- dim(mu.inits)[1]*dim(mu.inits)[2]*dim(mu.inits)[3]

inits <- function(){list(
  rho = rnorm(1, 0, 10),
  y = mu.inits + runif(total.mu.index, 0, 5),
  z = mu.inits + runif(total.mu.index, 0, 5),
  x = mu.inits + runif(total.mu.index, 0, 5),
  ex = mu.inits + runif(total.mu.index, 0, 5),
  tau.process = runif(1, 0, 0.1), 
  # tau.site = runif(1, 0, 0.1),
  tau.obs = runif(2, 0, 0.1)
  # alpha.species = rnorm(n.species, 0, 1),
  # alpha.site = rnorm(n.sites, 0, 1),
  # alpha.year = rnorm(n.years, 0, 1),
  # tau.species = runif(1, 0, 0.1)
)}

monitor <- c("tau.process", "tau.obs", "rho", "x")
if(species.effect) monitor <- c(monitor, "tau.species", "alpha.species")
if(site.effect) monitor <- c(monitor, "tau.site", "alpha.site")
if(year.effect) monitor <- c(monitor, "tau.year", "alpha.year")

x.ic <- matrix(runif(n.years*ncol(y), 0, 1), n.years, ncol(y))


constants <- list(
  n.weeks = weeks.per.year,
  n.plots = ncol(y.table),
  x.ic = x.ic, #runif(ncol(y), 0, 1),
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

save.path <- here("ModelOut",
                  ForecastProject.id,
                  rdata.name)

source("Functions/run_nimble_parallel.R")
cl <- makeCluster(n.cores)
samples <- run_nimble_parallel(
    cl,
    model,
    constants,
    data,
    inits,
    monitor,
    n.iter = 50000,
    max.iter = 2e6,
    thin = 10,
    check.interval = 5,
    check.params.only = TRUE,
    file.name = save.path
  )
stopCluster(cl)  




save(samples, file = save.path)





