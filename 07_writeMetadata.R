#### create forecast metadata ####

library(EML)
library(tidyverse)
library(lubridate)

fx.dir <- "ForecastSubmissionFiles"
fx.file <- "ticks-2019-03-04-BU_Dem"
file.dest <- file.path(fx.dir, paste0(fx.file, ".csv"))
fx <- read.csv(file.dest)
time <- fx %>% 
  pull(time) %>% 
  ymd() %>% 
  unique() %>% 
  sort()

# get lat.lon
ticks <- read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz")
lat.lon <- ticks %>% 
  select(c(siteID, decimalLatitude, decimalLongitude)) %>% 
  mutate(decimalLatitude = decimalLatitude, # nmme has one deg resolution 
         decimalLongitude = decimalLongitude) %>% 
  distinct()

n.ens <- 2000

forecast_project_id <- "BU_Dem"
forecast_model_id <- "79feb004ce3aec1c825227a470a3f8bd456f0a67"
forecast_iteration_id <- time[1]
forecast_issue_time <- lubridate::today()

attributes <- tibble::tribble(
  ~attributeName,           ~attributeDefinition,                          ~unit,                  ~formatString, ~numberType, ~definition,
  "time",                   "[dimension]{time}",                          "year",                 "YYYY-MM-DD",  "numberType", NA,
  "ensemble",               "[dimension]{index of ensemble member}",      "dimensionless",         NA,           "integer",    NA,
  "obs_flag",               "[dimension]{observation error}",             "dimensionless",         NA,           "integer",    NA,
  "Ixodes_scapularis",      "[variable]{Pop. abundace of I. scap.}",      "numberPerMeterSquared", NA,           "integer",    "Number of ticks observed per total sampled area",
  "Amblyomma_americanum",   "[variable]{Pop. abundace of A. amer.}",      "numberPerMeterSquared", NA,           "integer",    "Number of ticks observed per total sampled area",
  "forecast",               "[flag]{whether time step assimilated data}", "dimensionless",         NA,           "integer",    NA,
  "data_assimilation",      "[flag]{whether time step assimilated data}", "dimensionless",         NA,           "integer",    NA
) 

## note: EML uses a different unit standard than UDUNITS. For now use EML. EFI needs to provide a custom unitList.
attributes
attrList <- set_attributes(attributes, 
                           col_classes = c("Date", "numeric", "numeric","numeric", 
                                           "numeric","numeric", "numeric"))

## sets metadata about the file itself (name, file type, size, MD5, etc)
physical <- set_physical(file.dest,
                         recordDelimiter='\n')

## set metadata for the file as a whole
dataTable <- eml$dataTable(
  entityName = "forecast",  ## this is a standard name to allow us to distinguish this entity from 
  entityDescription = "Forecast of population size using a basic demographic model",
  physical = physical,
  attributeList = attrList)

me <- list(individualName = list(givenName = "John", 
                                 surName = "Foster"),
           electronicMailAddress = "fosterj@bu.edu")

taxa <- tibble::tribble(
  ~Genus,      ~Species,
  "Ixodes",    "scapularis",
  "Amblyomma", "americanum")

coverage <- 
  set_coverage(begin = first(time), 
               end = last(time),
               sci_names = taxa,
               geographicDescription = "NEON sites BLAN, ORNL, SCBI, SERC, KONZ, TALL, UKFS",
               west = min(lat.lon$decimalLongitude), 
               east = max(lat.lon$decimalLongitude), 
               north = max(lat.lon$decimalLatitude), 
               south = min(lat.lon$decimalLatitude))

keywordSet <- list(
  list(
    keywordThesaurus = "EFI controlled vocabulary",
    keyword = list("forecast",
                   "population",
                   "timeseries",
                   "ticks")
  ))

abstract <- "Demography state-spave model. Process error is Gaussian and observation 
error is Poisson. Cumalative growing degree days are used it estimate a observation 
probability for each week for each site. Survival is constant. Driver data is from  
NEONs IR temperature (ground temp) and NOAAs NMME max temp."

dataset = eml$dataset(
  title = "State-space survival and capture probability model",
  creator = me,
  contact = list(individualName = list(givenName = "John", 
                                       surName = "Foster")),
  pubDate = forecast_issue_time,
  intellectualRights = "http://www.lternet.edu/data/netpolicy.html.",
  abstract = abstract,
  dataTable = dataTable,
  keywordSet = keywordSet,
  coverage = coverage
)

additionalMetadata <- eml$additionalMetadata(
  metadata = list(
    forecast = list(
      ## Basic elements
      timestep = "1 week", ## should be udunits parsable; already in coverage -> temporalCoverage?
      forecast_horizon = paste0(time, " weeks"),
      forecast_issue_time = forecast_issue_time,
      forecast_iteration_id = forecast_iteration_id,
      forecast_project_id = forecast_project_id,
      metadata_standard_version = "0.3",
      model_description = list(
        forecast_model_id = forecast_model_id,
        name = "discrete State-space survival and capture probability model",
        type = "process-based",
        repository = "https://github.com/JohnRFoster/ForecastChallengeEFITick.git"
      ),
      ## MODEL STRUCTURE & UNCERTAINTY CLASSES
      initial_conditions = list(
        # Possible values: absent, present, data_driven, propagates, assimilates
        status = "propagates",
        # Number of parameters / dimensionality
        complexity = 2  ## [species 1, species 2] per plot
      ),
      drivers = list(
        status = "propagates"
      ),
      parameters = list(
        status = "present",
        complexity = 3 # survival, capture prob (beta + beta.site) per site
      ),
      random_effects = list(
        status = "absent"
      ),
      process_error = list(
        status = "propagates",
        propagation = list(
          type = "ensemble", # ensemble vs analytic
          size = n.ens          # required if ensemble
        ),
        complexity = 2,   
        covariance = FALSE
      ),
      obs_error = list(
        status = "present",
        complexity = 2,   
        covariance = FALSE
      )
    ) # forecast
  ) # metadata
) # eml$additionalMetadata

my_eml <- eml$eml(dataset = dataset,
                  additionalMetadata = additionalMetadata,
                  packageId = forecast_iteration_id , 
                  system = "datetime"  ## system used to generate packageId
)

eml.file <- file.path(fx.dir, paste0(fx.file, "-eml.xml"))
write_eml(my_eml, eml.file)

## check base EML
eml_validate(my_eml)

## check that the EML is also a valid EFI forecast
source("Functions/forecast_validator.R")
forecast_validator(my_eml)



