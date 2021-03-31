# function to divide data into species specific data sets
library(dplyr)

split_species <- function(spp){
 
   ixodes.plots <- c(
    "BLAN_012",
    "BLAN_005",
    "SCBI_013",
    "SCBI_002",
    "SERC_001",
    "SERC_005",
    "SERC_006",
    "SERC_012",
    "ORNL_007"  
  )
  
  ambloyomma.plots <- c(
    "SCBI_013",
    "SERC_001",
    "SERC_005",
    "SERC_006",
    "SERC_002",
    "SERC_012",
    "KONZ_025",
    "UKFS_001",
    "UKFS_004",
    "UKFS_003",
    "ORNL_002",
    "ORNL_040",
    "ORNL_008",
    "ORNL_007",
    "ORNL_009",
    "ORNL_003",
    "TALL_001",
    "TALL_008",
    "TALL_002"
  )
  
  # first load the target data set
  # data <- read.csv("Data/ticks-targets.csv.gz")
  data <- read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz")
  data <- data %>% 
    select(-all_of(c("RHMin_precent",
                     "RHMin_variance",
                     "RHMax_precent",
                     "RHMax_variance",
                     "airTempMin_degC",
                     "airTempMin_variance",
                     "airTempMax_degC",
                     "airTempMax_variance")))
  
  if(spp == "Ixodes"){
    data.return <- data %>% 
      filter(plotID %in% ixodes.plots) %>% 
      select(-amblyomma_americanum) %>% 
      group_by(plotID) %>%
      distinct(yearWeek, ixodes_scapularis, .keep_all = TRUE) %>% 
      mutate(indicator = !yearWeek %in% yearWeek[duplicated(yearWeek)]) %>% 
      filter(ixodes_scapularis >= 0 | indicator) %>% 
      select(-indicator)
  } else if(spp == "Amblyomma"){
    data.return <- data %>% 
      filter(plotID %in% ambloyomma.plots) %>% 
      select(-ixodes_scapularis) %>% 
      group_by(plotID) %>%
      distinct(yearWeek, amblyomma_americanum, .keep_all = TRUE) %>% 
      mutate(indicator = !yearWeek %in% yearWeek[duplicated(yearWeek)]) %>% 
      filter(amblyomma_americanum >= 0 | indicator) %>% 
      select(-indicator)
  }
  
  return(data.return)
  
}
