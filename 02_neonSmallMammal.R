# ====================================================== #
#            NEON small mammal data download             #
# ====================================================== #

library(neonstore)
library(tidyverse)
library(lubridate)
library(uuid)

# target sites
target.sites <- c("BLAN", "ORNL", "SCBI", "SERC", "KONZ", "TALL", "UKFS")

message("Updating NEON Small Mammal Data...")
# neon_download(product = "DP1.10072.001",
#               end_date = "2019-12-31",   
#               site = target.sites,       
#               .token = Sys.getenv("NEON_TOKEN"))  

# neon_store("DP1.10072.001", db_dir = Sys.getenv("NEONSTORE_HOME"))

# tables - can join by "nightuid"
# plot.night <- neon_read("mam_perplotnight") # per trapping grid data
trap.night <- neon_read("mam_pertrapnight") # per trap night data

# bout = consecutive nights trapped
# 3 nights = pathogen grid
# 1 night = diversity grid

# some plots seem to switch from being designated a pathogen grid to diversity
# will need to check for consecutive nights to determine bouts
# there should only be one bout per month, just depends on whether a month
# has one night of trapping or three consecutive nights
# collectYearMonth with plotID identifies unique bouts 

# determine bouts and keep select columns
smam.data <- trap.night %>% 
  mutate(collectYear = year(collectDate),
         collectMonth = month(collectDate),
         collectYearMonth = paste(collectYear, collectMonth, sep = "-")) %>%  # for determining bouts
  select(all_of(c("uid", # unique id
                  "nightuid", # night id, maps to plot.night
                  "collectDate",
                  "collectYear",
                  "collectMonth",
                  "collectYearMonth",
                  "domainID",
                  "siteID",
                  "plotID",
                  "nlcdClass", # land cover classification
                  "decimalLatitude",
                  "decimalLongitude",
                  "trapStatus", 
                  "trapCoordinate", # trap location on grid - A1
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

# generate a uuid for bouts
smam.bout <- smam.data %>% 
  group_by(plotID, collectYearMonth) %>% 
  mutate(boutuid = UUIDgenerate()) %>% 
  ungroup()

# want the number of unique animals each bout
# will go by unique tag
# some tags have been replaced, deal with those first
# the majority of tags replaced are just left or right
# and the tag number is the same
# Some have a new tag number, need to demarcate those
# the tagID in the replacedTag column is the old tag

# arrange by collect date to make replacing easier
smam.bout <- smam.bout %>% 
  group_by(plotID) %>% 
  arrange(collectDate)

tag.pattern <- "[[:upper:]]\\d{4}" # end of tagID - what is replaced; old tag
new.tag.rows <- grep(tag.pattern, smam.bout$replacedTag) # rows when new tags used
old.tag.id <- smam.bout$replacedTag[new.tag.rows] # old tag Ids
new.tag.id <- smam.bout$tagID[new.tag.rows] # the tags that were replaced

for(i in seq_along(new.tag.id)){
  
  smam.subset <- smam.bout %>%
    filter(tagID == new.tag.id[i])
    
  if(nrow(smam.subset) > 0){ # revert back to old tag
    index <- which(smam.bout$tagID == new.tag.id[i])
    smam.bout$tagID[index] <- gsub(tag.pattern, old.tag.id[i], new.tag.id[i]) 
  }
}


smam.bout %>% 
  filter(tagID == old.tag.id[3]) %>% 
  view()

# need to create a number animals column first before making wider
smam.data.wide <- smam.bout %>% 
  mutate(animalInTrap = ifelse(trapStatus == "4 - more than 1 capture in one trap" | 
                                 trapStatus == "5 - capture",
                               1, 0)) %>% 
  pivot_wider(names_from = scientificName,
              values_from = animalInTrap,
              values_fill = 0)

smam.data.wide %>% 
  group_by(boutuid) %>% 
  summarise()
  


replaced.tag <- smam.bout %>% 
  filter(!is.na(replacedTag))

table(replaced.tag$replacedTag)





# filtering out traps that were found with a dead animal
# this also removes all the NA entries for the 'fate' column
# which are just traps without animals because they didn't 
# catch anything or the trap was not set
smam.alive <- smam.data %>% filter(fate != "dead")

# looks like all escaped data has an ID of some sort
smam.escaped <- smam.data %>% filter(fate == "escaped")
# smam.escaped %>% pull(scientificName) %>% table()


smam.data %>% pull(recapture) %>% table()
smam.alive %>% pull(recapture) %>% table()

smam.data %>% filter(recapture == "Y") %>% pull(fate) %>% table()

smam.data %>% pull(uid) %>% unique() %>% length()
smam.data %>% pull(nightuid) %>% unique() %>% length()
smam.data %>% pull(collectDate) %>% unique() %>% length()

# need to deal with trap status
# trap status codes:
smam.data %>% pull(trapStatus) %>% table()
# 1 - trap not set   
# 2 - trap disturbed/door closed but empty 
# 3 - trap door open or closed w/ spoor left        
# 4 - more than 1 capture in one trap 
# 5 - capture                     
# 6 - trap set and empty 

# 1 - trap not set  
# want the total number of traps per trapnight that are not set
trap.not.set <- smam.data %>% 
  filter(trapStatus == "1 - trap not set") %>% 
  group_by(nightuid) %>% 
  summarise(totalTrapsNotSet = n())

trap.not.set.per.bout <- smam.data %>% 
  filter(trapStatus == "1 - trap not set") %>% 
  group_by(plotID, collectYearMonth) %>% 
  summarise(totalTrapsNotSet = n())

# trapStatus 2 and 3 suggest an animal could be present, but we don't
# know if the animal was captured in another trap (3) or if the disturbance
# was from something other than an animal (2). So, not going to assume
# these are animals
trap.empty <- smam.data %>% 
  filter(trapStatus == "2 - trap disturbed/door closed but empty" |
           trapStatus == "3 - trap door open or closed w/ spoor left" |
           trapStatus == "6 - trap set and empty")

# these records map to same nightuid and trap with different tagIDs
# want total animals in traps
trap.two.or.more <- smam.data %>% 
  filter(trapStatus == "4 - more than 1 capture in one trap") %>% 
  group_by(nightuid) %>% 
  summarise(animalsAlive = n())

captures <- smam.data %>% 
  filter(trapStatus == "5 - capture")


# first thing I want is total animal caught
total.animals <- smam.data %>% 
  filter(trapStatus == "4 - more than 1 capture in one trap" |
           trapStatus == "5 - capture") %>%
  group_by(nightuid) %>% 
  summarise(totalAnimalsObserved = n())

total.animals.per.bout <- smam.data %>% 
  filter(trapStatus == "4 - more than 1 capture in one trap" |
           trapStatus == "5 - capture") %>%
  group_by(plotID, collectYearMonth) %>% 
  summarise(totalAnimalsObserved = n())





smam.data.restructure <- left_join(smam.data, total.animals, by = "nightuid")






# we know that codes 5 and 6 represent good efforts and are easiest to deal with

# measures we may want to use
# minimum number of hosts alive each trapping bout
# density of all hosts alive at each trapping bout
# species richness at each trapping bout
# total number Peromyscus alive at each trapping bout
# Peromyscus density alive at each trapping bout 









write_csv(smam.data, "Data/smallMammalTrapNightRaw.csv")
write_csv(smam.data, "Data/smallMammalTrapNightRaw.csv")
write_csv(smam.data, "Data/smallMammalTrapNightRaw.csv")
