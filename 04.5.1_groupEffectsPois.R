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
library(MMWRweek) # to convert year-weeks to date

start.epi.weeks <- c(10, 14, 19, 23, 28, 32, 36, 41) # all possible start weeks
met.weeks <- start.epi.weeks[2:length(start.epi.weeks)] - 1
day.run <- lubridate::today() # the day the script is called
day.run <- "2021-05-25"

# anytime we run this script before the start of the challenge we want to forecast all 2019 target weeks
if(day.run < "2021-03-31"){ 
  start.week <- start.epi.weeks[1]
  end.met.obs <- met.weeks[1]
} else { # otherwise use the appropriate starting week (months are 2 ahead)
  start.week <- start.epi.weeks[month(day.run) - 2]
  end.met.obs <- met.weeks[month(day.run) - 2]
  nmme.start <- month(day.run) + 1
}
filter.week <- paste0("2019", start.week) %>% as.integer()
filter.week.met <- paste0("2019", end.met.obs)


# where do we want random effects?
species.effect <- FALSE
site.effect <- FALSE
year.effect <- FALSE

# n.cores <- as.numeric(Sys.getenv("NSLOTS")) # number of cores to use for cluster 
n.cores <- 3


group.dir <- "poisDataArea"

# for getting met variable
csv <- "Data/BioTemperatureWeekly.csv"
met.var <- "cumGDD"
met.uncertainty <- "cumGDDVariance"

ForecastProject.id <- "BetaSiteCGDDLogitTempSurvivalTruncatedTheta"     # Some ID that applies to a set of forecasts
out.dir <- here("ModelOut", group.dir)
if(!dir.exists(out.dir)) dir.create(out.dir)

# split data by species
source("Functions/split_species.R")
data.ixodes <- split_species("Ixodes")
data.ambloyoma <- split_species("Amblyomma")

plots.ixodes <- data.ixodes %>% 
  pull(plotID) %>% 
  unique()
plots.amblyomma <- data.ambloyoma %>% 
  pull(plotID) %>% 
  unique()

source("Functions/make_wide.R")
wide.data <- make_wide(data.ixodes, data.ambloyoma)

### set forecast variables ###
end.week <- 44 # does not change
target.weeks <- paste0("2019", start.week:end.week)
na.mat <- data.frame(yearWeek = target.weeks,
                     time = as.character(MMWRweek2Date(rep(2019, length(start.week:end.week)), 
                                                       start.week:end.week)))

wide.data <- bind_rows(wide.data, na.mat)

# get lat.lon
ticks <- read_csv("Data/ticks-targets.csv.gz")
lat.lon <- ticks %>% 
  select(c(siteID, decimalLatitude, decimalLongitude)) %>% 
  mutate(decimalLatitude = round(decimalLatitude), # nmme has one deg resolution 
         decimalLongitude = round(decimalLongitude)) %>% 
  distinct()

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
y <- area <- array(NA, dim = c(n.years, # years
                       max.weeks, # weeks
                       ncol(y.table))) # site


source("Functions/get_met_array.R")
met.list <- get_met_array(csv, 
                          weeks.per.year,
                          filter.week.met,
                          year.weeks,
                          met.var,
                          met.uncertainty)
met.vals <- met.list$met.vals
x.var <- met.list$met.uncertainty
met.inits <- met.list$met.inits

survival.met <- get_met_array(csv, 
                              weeks.per.year,
                              filter.week.met,
                              year.weeks,
                              "bioTempMaximum",
                              "bioTempMaximumVariance")
survival.met.vals <- survival.met$met.vals
survival.x.var <- survival.met$met.uncertainty
survival.met.inits <- survival.met$met.inits

source("Functions/make_nmme_ens.R")
start.month <- paste0("20190", nmme.start, "01")
message("Getting NMME from ", start.month)
for(s in 1:n.sites){
  met.fx <- make_nmme_ens(var = "tasmax", 
                          start.month = start.month,
                          dec.lat = lat.lon$decimalLatitude[s], 
                          dec.lon = lat.lon$decimalLongitude[s],
                          gdd = TRUE) %>% as_tibble() %>% 
    filter(epiWeek > end.met.obs & epiWeek <= end.week)
  
  start.cumgdd <- met.vals[4,3,s]

  met.fx.epi <- select(met.fx, epiWeek)
  met.fx.gdd <- met.fx 
  met.fx.cumgdd <- matrix(NA, nrow(met.fx.gdd), ncol(met.fx.gdd))
  for(g in 1:ncol(met.fx.cumgdd)){
    cumgdd <- cumsum(c(start.cumgdd, pull(met.fx.gdd, g)))
    met.fx.cumgdd[,g] <- cumgdd[-1]
  }
  cumgdd.mu <- apply(met.fx.cumgdd, 1, mean)
  cumgdd.var <- apply(met.fx.cumgdd, 1, var)
  met.fx.cumgdd <- tibble(epiWeek = met.fx.epi,
                          Year = 2019,
                          mean = cumgdd.mu,
                          variance = cumgdd.var) %>% 
    group_by(epiWeek) %>% 
    slice(which.min(mean))
  
  place <- seq(4, length.out = nrow(met.fx.cumgdd), by = 1)
  met.vals[4, place, s] <- pull(met.fx.cumgdd, "mean")
  x.var[4, place, s] <- pull(met.fx.cumgdd, "variance")
  
  nmme.mu <- apply(met.fx, 1, mean)
  nmme.var <- apply(met.fx, 1, var)
  met.fx.temp <- tibble(epiWeek = met.fx.epi,
                        Year = 2019,
                        mean = nmme.mu,
                        variance = nmme.var) %>% 
    group_by(epiWeek) %>% 
    slice(which.min(mean))
  survival.met.vals[4, place, s] <- pull(met.fx.temp, "mean")
  survival.x.var[4, place, s] <- pull(met.fx.temp, "variance")
}

# pull out driver variance, check for NAs, set to mean
if(any(is.na(x.var))){
  x.var[is.na(x.var)] <- mean(x.var, na.rm = TRUE)
}
if(any(is.na(survival.x.var))){
  survival.x.var[is.na(survival.x.var)] <- mean(survival.x.var, na.rm = TRUE)
}
x.var[which(x.var == 0)] <- 0.1
survival.x.var[which(survival.x.var == 0)] <- 0.1

# scaling met by standard deviation by site
met.inits <- met.inits / sd(met.inits, na.rm = TRUE)
survival.met.inits <- survival.met.inits / sd(survival.met.inits, na.rm = TRUE)
for(ss in 1:n.sites){
  site.sd <- sd(met.vals[,,ss], na.rm = TRUE)
  temp.site.sd <- sd(survival.met.vals[,,ss], na.rm = TRUE)
  
  met.vals[,,ss] <- met.vals[,,ss] / site.sd
  survival.met.vals[,,ss] <- survival.met.vals[,,ss] / temp.site.sd

  # variance is scaled by the square of the scale (sd)
  x.var[,,ss] <- x.var[,,ss] / site.sd^2
  survival.x.var[,,ss] <- survival.x.var[,,ss] / temp.site.sd^2
}

for(i in 1:n.years){
  y.sub <- wide.data %>% 
    as_tibble() %>% 
    filter(year(time) == years[i]) 
  
  y[i, 1:weeks.per.year[i], 1:ncol(y.table)] <- as.matrix(y.sub[,-c(1,2)]) 
  
  # area.sub <- wide.data.area %>% 
  #   filter(year(time) == years[i]) 
  # 
  # area[i, 1:weeks.per.year[i], 1:ncol(y.table)] <- as.matrix(area.sub[,-c(1,2)]) 
}

area[is.na(area)] <- 0

# load previous forecast
load("ModelOut/poisDataArea/BetaSiteCGDDLogitTempSurvival_2019-14.RData")
out.mcmc <- samples$samples %>% as.matrix()
# rm(save.ls)
x.cols <- grep("x[", colnames(out.mcmc), fixed = TRUE)
params <- out.mcmc[,-x.cols]
states <- out.mcmc[,x.cols]
rm(out.mcmc)

states <- apply(states, 2, pmax, 0)
states.mu <- apply(states, 2, mean)
# states.var <- apply(states, 2, var)

x.init <- array(0, dim = dim(y))
n.plots <- ncol(y.table)
for(p in 1:n.plots){
  for(k in 1:n.years){
    for(t in 1:weeks.per.year[k]){
      x.name <- paste0("x[", k, ", ", t, ", ", p, "]")
      x.init[k, t, p] <- states.mu[x.name]
    }
  }
}

params.mu <- apply(params, 2, mean)
params.var <- apply(params, 2, var)

beta.mu <- params.mu["beta"]
beta.tau <- 1 / params.var["beta"]
beta.site.mu <- params.mu[grep("beta.site[", names(params.mu), fixed = TRUE)]

psi.mu <- params.mu["psi"]
psi.tau <- 1 / params.var["psi"]
psi.site.mu <- params.mu[grep("psi.site[", names(params.mu), fixed = TRUE)]

tau.beta.site.shape <- params.mu["tau.beta.site"]^2 / params.var["tau.beta.site"]
tau.beta.site.rate <- params.mu["tau.beta.site"] / params.var["tau.beta.site"]
tau.psi.site.shape <- params.mu["tau.psi.site"]^2 / params.var["tau.psi.site"]
tau.psi.site.rate <- params.mu["tau.psi.site"] / params.var["tau.psi.site"]
tau.process.shape <- params.mu["tau.process"]^2 / params.var["tau.process"]
tau.process.rate <- params.mu["tau.process"] / params.var["tau.process"]

data <- list(
  y = y,
  x.obs = met.vals,
  x.tau = 1 / x.var,
  x.obs.phi = survival.met.vals,
  x.phi.tau = 1 / survival.x.var,
  beta.mu = beta.mu,
  beta.tau = beta.tau,
  beta.site.mu = beta.site.mu,
  psi.mu = psi.mu,
  psi.tau = psi.tau,
  psi.site.mu = rep(0, n.sites),
  tau.beta.site.shape = tau.beta.site.shape,
  tau.beta.site.rate = tau.beta.site.rate,
  tau.psi.site.shape = tau.psi.site.shape,
  tau.psi.site.rate = tau.psi.site.rate,
  tau.process.shape = tau.process.shape,
  tau.process.rate = tau.process.rate
)

mu.inits <- y
mu.inits[is.na(mu.inits)] <- mean(y, na.rm = TRUE)

total.mu.index <- dim(mu.inits)[1]*dim(mu.inits)[2]*dim(mu.inits)[3]

x.ic <- matrix(runif(n.years*ncol(y), 500, 1000), n.years, ncol(y))

constants <- list(
  n.weeks = weeks.per.year,
  n.plots = n.plots,
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

inits <- function(){list(
  # phi = array(runif(total.mu.index, 0.95, 1), dim = dim(y)),
  # pi = array(runif(total.mu.index, 0.4, 0.6), dim = dim(y)),
  theta = runif(1, 0, 0.1),
  beta = runif(1, params.mu["beta"]-0.1, params.mu["beta"]+0.1),
  psi = runif(1, 1, 2),
  beta.site = runif(n.sites, 0.9, 1.1),
  psi.site = runif(n.sites, 0.9, 1.1),
  x = round(x.init),    # + runif(total.mu.index, 0, 5)),
  y = abs(round(x.init)),    # + jitter(x.init) + runif(total.mu.index, 0, 5)),
  ex = round(x.init),  # + runif(total.mu.index, 0, 5),
  ex.z = round(x.init), # + runif(total.mu.index, 6, 8),
  cgdd = met.inits,
  tau.beta.site = abs(rnorm(1, params.mu["tau.beta.site"], 0.001)),
  tau.psi.site = abs(rnorm(1, params.mu["tau.psi.site"], 0.001)),
  tau.process = abs(rnorm(1, params.mu["tau.process"], 0.001))
)}

source("Models/yearTimeSeriesCGDDPois.R")
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
if(species.effect) file.name <- paste0(file.name, "SpeciesEffect")
if(site.effect) file.name <- paste0(file.name, "SiteEffect")
if(year.effect) file.name <- paste0(file.name, "YearEffect")
rdata.name <- paste0(file.name, 
                     paste0("_2019-", start.week),
                     ".RData")

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
    n.iter = 50000,
    max.iter = 3e6,
    thin = 10,
    check.interval = 2,
    check.params.only = TRUE,
    file.name = save.path
  )
stopCluster(cl)  




save(samples, model,
     file = save.path)





