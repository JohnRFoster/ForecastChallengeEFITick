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

efi_server <- TRUE

start.epi.weeks <- c(10, 14, 19, 23, 28, 32, 36, 41) # all possible start weeks
met.weeks <- start.epi.weeks[2:length(start.epi.weeks)] - 1
day.run <- lubridate::today() # the day the script is called

# anytime we run this script before the start of the challenge we want to forecast all 2019 target weeks
if(day.run < "2021-03-31"){ 
  start.week <- start.epi.weeks[1]
  end.met.obs <- met.weeks[1]
} else { # otherwise use the appropriate starting week (months are 2 ahead)
  start.week <- start.epi.weeks[month(day.run) - 2]
  end.met.obs <- met.weeks[month(day.run) - 2]
}
filter.week <- paste0("2019", start.week) %>% as.integer()
filter.week.met <- paste0("2019", end.met.obs)


# where do we want random effects?
species.effect <- FALSE
site.effect <- FALSE
year.effect <- FALSE

# n.cores <- as.numeric(Sys.getenv("NSLOTS")) # number of cores to use for cluster 
n.cores <- 3

source("Models/yearTimeSeriesCGDDPois.R")
group.dir <- "poisDataArea"

# for getting met variable
csv <- "Data/BioTemperatureWeekly.csv"
met.var <- "cumGDD"
met.uncertainty <- "cumGDDVariance"

ForecastProject.id <- "BetaSiteCGDDLogit"     # Some ID that applies to a set of forecasts
out.dir <- here("ModelOut", group.dir)
if(!dir.exists(out.dir)) dir.create(out.dir)

# split data by species
source("Functions/split_species.R")
data.ixodes <- split_species("Ixodes")
data.ambloyoma <- split_species("Amblyomma")

# get lat.lon
ticks <- read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz")
lat.lon <- ticks %>% 
  select(c(siteID, decimalLatitude, decimalLongitude)) %>% 
  mutate(decimalLatitude = round(decimalLatitude), # nmme has one deg resolution 
         decimalLongitude = round(decimalLongitude)) %>% 
  distinct()

# get unique sites
sites.ixodes <- unique(data.ixodes$siteID)
sites.amblyomma <- unique(data.ambloyoma$siteID)

# get unique plots
plots.ixodes <- unique(data.ixodes$plotID)
plots.amblyomma <- unique(data.ambloyoma$plotID)

wide.ixodes <- data.ixodes %>% 
  filter(plotID %in% plots.ixodes) %>% 
  filter(Year >= 2016) %>%
  filter(yearWeek < filter.week) %>% 
  # mutate(density = Ixodes_scapularis / totalSampledArea * 1600) %>%
  pivot_wider(id_cols = c(yearWeek, time), 
              names_from = plotID, 
              values_from = ixodes_scapularis,
              names_glue = "Ixodes_{plotID}") 

wide.amblyomma <- data.ambloyoma %>% 
  filter(plotID %in% plots.amblyomma) %>% 
  filter(Year >= 2016) %>%
  filter(yearWeek < filter.week) %>% 
  # mutate(density = amblyomma_americanum / totalSampledArea * 1600) %>% 
  pivot_wider(id_cols = c(yearWeek, time), 
              names_from = plotID, 
              values_fn = {max},
              values_from = amblyomma_americanum,
              names_glue = "Amblyomma_{plotID}")


wide.data <- full_join(wide.ixodes, wide.amblyomma,
                       by = c("yearWeek", "time")) %>% 
  arrange(yearWeek) %>% 
  mutate(yearWeek = as.character(yearWeek),
         time = as.character(time))



### set forecast variables ###
end.week <- 44 # does not change
target.weeks <- paste0("2019", start.week:end.week)
na.mat <- data.frame(yearWeek = target.weeks,
                     time = as.character(ISOweek2date(paste0("2019-W", start.week:end.week, "-1"))))

wide.data <- bind_rows(wide.data, na.mat)

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
y <- array(NA, dim = c(n.years, # years
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

source("Functions/make_nmme_ens.R")
for(s in 1:n.sites){
  met.fx <- make_nmme_ens(var = "tasmax", 
                          start.month = "20190401",
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
}

# pull out driver variance, check for NAs, set to mean
if(any(is.na(x.var))){
  x.var[is.na(x.var)] <- mean(x.var, na.rm = TRUE)
}
x.var[which(x.var == 0)] <- 0.1

# scaling met by standard deviation by site
for(ss in 1:n.sites){
  site.sd <- sd(met.vals[,,ss], na.rm = TRUE)
  met.vals[,,ss] <- met.vals[,,ss] / site.sd
  
  # variance is scaled by the square of the scale (sd)
  x.var[,,ss] <- x.var[,,ss] / site.sd^2
}

for(i in 1:n.years){
  y.sub <- wide.data %>% 
    as_tibble() %>% 
    filter(year(time) == years[i]) 
  
  y[i, 1:weeks.per.year[i], 1:ncol(y.table)] <- as.matrix(y.sub[,-c(1,2)]) 
}

data <- list(
  y = y,
  x.obs = met.vals,
  x.tau = 1 / x.var
)



load("ModelOut/poisDataArea/BetaSiteCGDDLogit.RData")


out.mcmc <- as.matrix(samples$samples)
x.cols <- grep("x[", colnames(out.mcmc), fixed = TRUE)
params <- out.mcmc[,-x.cols] 
z <- apply(out.mcmc[,x.cols], 2, pmax, 0)

pred <- pi <- array(NA, dim = c(n.years, 
                                max(weeks.per.year), 
                                ncol(y.table), 
                                nrow(params)))
for(p in 1:ncol(y.table)){
  print(p)
  for(k in 1:n.years){
    cgdd.p.k <- met.vals[k, 1:weeks.per.year[k], site.index[p]]
    cgdd.sd <- sqrt(x.var[k, 1:weeks.per.year[k], site.index[p]])
    
    # na.rows <- which(is.na(cgdd.p.k))
    # cgdd.p.k[na.rows] <- approx(as.vector(cgdd.p.k), xout = na.rows)$y 
    
    cgdd.p.k[is.na(cgdd.p.k)] <- mean(cgdd.p.k, na.rm = TRUE)
    
    for(t in 1:weeks.per.year[k]){
      z.col <- paste0("x[", k, ", ", t, ", ", p, "]")
      z.latent <- z[,z.col]
      
      cgdd <- rnorm(nrow(params), cgdd.p.k[t], cgdd.sd[t])
      pi.x <- params[, "beta"] + params[,"beta.site[1]"] * cgdd
      pi[k, t, p, ] <- boot::inv.logit(pi.x)
      lambda <- pi[k, t, p, ]*z.latent
      pred[k, t, p, ] <- rpois(nrow(params), lambda)
    }
  }
}

pred.2019 <- 1:weeks.per.year[4]
target.index <- start.week:end.week
fx.start <- which(target.index == start.week)
fx.index <- fx.start:last(pred.2019)
species.fx <- c(rep("ixodes_scapularis", length(plots.ixodes)), 
                 rep("amblyomma_americanum", length(plots.amblyomma)))


obs.error <- 2
data.assimilation <- 0
forecast <- 0
date.col <- ISOweek2date(paste0("2019-W", target.index, "-1")) %>% as.character()

n.ens <- 2000
draws <- sample.int(nrow(params), n.ens)

check <- read.csv("Data/ticks-2019-03-04-tickGlobalNull_RandomWalk.csv.gz")

fx.df <- tibble()
for(p in 1:28){
  fx.ens <- t(pred[4,fx.index,p,draws])
  colnames(fx.ens) <- date.col
  
  fx.plot <- plot.vec[p]
  fx.site <- site.vec[p]
  fx.species <- species.fx[p]
  
  fx <- fx.ens %>% 
    as_tibble() %>% 
    pivot_longer(all_of(date.col), 
                 names_to = "time",
                 values_to = fx.species) %>%
    mutate(plotID = fx.plot,
           siteID = fx.site,
           ensemble = rep(1:nrow(fx.ens), each = length(target.index)),
           data_assimilation = data.assimilation,
           forecast = forecast,
           obs_flag = obs.error) %>% 
    pivot_longer(cols = all_of(fx.species), 
                 names_to = "species") # %>%  ## make species a column
    # group_by(time, plotID, siteID, obs_flag, species, forecast, data_assimilation) %>%
    # summarize(mean = mean(value),
    #           sd   = sd(value),
    #           Pred_interv_02.5 = quantile(value, 0.025),
    #           Pred_interv_97.5 = quantile(value, 0.975)) %>%
    # pivot_longer(cols = c(mean, sd, Pred_interv_02.5, Pred_interv_97.5),
    #              names_to = "statistic")
  
    fx.df <- bind_rows(fx.df, fx)
}

fx.submit <- fx.df %>% 
  pivot_wider(names_from = species, 
              values_from = value) ## go back to species wide

file.name <- paste0("ticks-", date.col[1], "-BU_Dem.csv")
file.dest <- file.path("ForecastSubmissionFiles", file.name)
write_csv(fx.submit, file.dest)

if(efi_server){
  library(aws.s3) 
  
  Sys.setenv("AWS_DEFAULT_REGION" = "data",
             "AWS_S3_ENDPOINT" = "ecoforecast.org")
  
  aws.s3::put_object(file = file.dest,
                     object = file.name, 
                     bucket = "submissions")  
}

