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
plot.night <- neon_read("mam_perplotnight") # per trapping grid data
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
                  # "plotID",
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
  group_by(siteID, collectYearMonth) %>% 
  mutate(boutuid = UUIDgenerate()) %>% 
  ungroup()

# pull out site/plot info for joining later
smam.site.info <- smam.bout %>% 
  select(all_of(c("boutuid",
                  "collectYear",
                  "collectMonth",
                  "collectYearMonth",
                  "domainID",
                  "siteID",
                  # "plotID",
                  "nlcdClass", # land cover classification
                  "decimalLatitude",
                  "decimalLongitude"))) %>% 
  distinct(boutuid, .keep_all = TRUE)

# want the number of unique animals each bout
# will go by unique tag
# some tags have been replaced, deal with those first
# the majority of tags replaced are just left or right
# and the tag number is the same
# Some have a new tag number, need to demarcate those
# the tagID in the replacedTag column is the old tag

# arrange by collect date to make replacing easier
smam.bout <- smam.bout %>% 
  group_by(siteID) %>% 
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

# change NA in scientificName to unknown
smam.bout$scientificName[is.na(smam.bout$scientificName)] <- "unknown sp."

# get scientific names vector
animals.names <- smam.bout %>% 
  filter(fate != "dead") %>% # only want the animals found alive
  pull(scientificName) %>% 
  unique()

# need to create a number animals column first before making wider
smam.data.wide <- smam.bout %>% 
  ungroup() %>% 
  filter(is.na(fate) | fate != "dead") %>% # only want the animals found alive, keep NAs for 0
  # select(-trapCoordinate) %>% 
  mutate(animalInTrap = ifelse(trapStatus == "4 - more than 1 capture in one trap" | 
                                 trapStatus == "5 - capture",
                               1, 0)) %>% 
  group_by(boutuid) %>% 
  pivot_wider(names_from = scientificName,
              values_from = animalInTrap,
              values_fill = 0) 

# need to make sure that we only count individuals once per bout
smam.wide.individuals <- smam.data.wide %>% 
  group_by(boutuid) %>% 
  distinct(tagID, .keep_all = TRUE)

# change NA column to unknown sp.
# index <- which(colnames(smam.wide.individuals) == "NA")
# colnames(smam.wide.individuals)[index] <- "unknown sp."

# want various animal totals
animals.by.bout <-  smam.wide.individuals %>% 
  group_by(boutuid) %>% 
  summarise(vars(animal.names), across(all_of(animals.names), ~ sum(.))) %>% # each animal
  mutate(totalPeromyscus = rowSums(across(starts_with("Peromyscus")))) %>% # all Peromyscus
  mutate(totalAnimals = rowSums(across(all_of(animals.names)))) %>%  # total animals
  mutate(speciesRichness = rowSums(across(all_of(animals.names)) > 0))

# join
smam.final <- left_join(animals.by.bout, smam.site.info, by = "boutuid")

write_csv(smam.data, "Data/smallMammalTrapNightRaw.csv")
write_csv(smam.final, "Data/smallMammalAnimalsByBout.csv")

