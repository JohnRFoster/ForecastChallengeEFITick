# ====================================================== #
#            NEON small mammal data download             #
# ====================================================== #

library(neonstore)
library(tidyverse)
library(lubridate)

# target sites
target.sites <- c("BLAN", "ORNL", "SCBI", "SERC", "KONZ", "TALL", "UKFS")

message("Updating NEON Small Mammal Data...")
# neon_download(product = "DP1.10072.001",
#               end_date = "2019-12-31",   
#               site = target.sites,       
#               .token = Sys.getenv("NEON_TOKEN"))  

# neon_store("DP1.10072.001", db_dir = Sys.getenv("NEONSTORE_HOME"))

# tables - can join by "nightuid"
plot.night <- neon_read("mam_perplotnight") # per trapping grid data
trap.night <- neon_read("mam_pertrapnight") # per trap night data

smam.data <- trap.night %>% 
  select(all_of(c("uid", # unique id
                  "nightuid", # night id, maps to plot.night
                  "domainID",
                  "siteID",
                  "plotID",
                  "nlcdClass", # land cover classification
                  "decimalLatitude",
                  "decimalLongitude",
                  "trapStatus", 
                  "trapCoordinate", # trap location on grid - A1
                  "collectDate",
                  "endDate",
                  "tagID",
                  "taxonID", # four letter species/ID code??
                  "scientificName",
                  "taxonRank",
                  "identificationQualifier", # not sure
                  "recapture",
                  "fate",
                  "replacedTag",
                  "lifeStage",
                  "larvalTicksAttached",
                  "nymphalTicksAttached",
                  "adultTicksAttached")))

smam.data %>% pull(uid) %>% unique() %>% length()
smam.data %>% pull(nightuid) %>% unique() %>% length()
smam.data %>% pull(collectDate) %>% unique() %>% length()

# looks like collectDate and endDate match
which(smam.data$collectDate == smam.data$endDate) %>% length()
nrow(smam.data)

smam.data <- smam.data %>% select(!endDate)

# need to deal with trap status
# trap status codes:
# 1 - trap not set   
# 2 - trap disturbed/door closed but empty 
# 3 - trap door open or closed w/ spoor left        
# 4 - more than 1 capture in one trap 
# 5 - capture                     
# 6 - trap set and empty 

smam.data %>% pull(trapStatus) %>% table()

# we know that codes 5 and 6 represent good efforts and are easiest to deal with

# measures we may want to use
# minimum number of hosts alive each trapping bout
# density of all hosts alive at each trapping bout
# species richness at each trapping bout
# total number Peromyscus alive at each trapping bout
# Peromyscus density alive at each trapping bout 

# bout = consecutive nights trapped
# 3 nights = pathogen grid
# 1 night = diversity grid

# some plots seem to switch from being designated a pathogen grid to diversity
# will need to check for consecutive nights to determine bouts
# there should only be one bout per month, just depends on whether a month
# has one night of trapping or three consecutive nights




# 1 - trap not set  
# want the total number of traps per trapnight that are not set
trap.not.set <- smam.data %>% 
  filter(trapStatus == "1 - trap not set") %>% 
  group_by(nightuid) %>% 
  summarise(totalTrapsNotSet = n())





write_csv(smam.data, "Data/smallMammalTrapNightRaw.csv")
write_csv(smam.data, "Data/smallMammalTrapNightRaw.csv")
write_csv(smam.data, "Data/smallMammalTrapNightRaw.csv")
