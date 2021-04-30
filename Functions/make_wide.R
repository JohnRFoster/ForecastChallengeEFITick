library(tidyverse)

make_wide <- function(data.ixodes = NULL, data.amblyomma = NULL){
  
  # get unique sites
  if(!is.null(data.ixodes)){
    plots.ixodes <- unique(data.ixodes$plotID)
    
    wide.ixodes <- data.ixodes %>% 
      filter(plotID %in% plots.ixodes) %>% 
      filter(Year >= 2016) %>%
      filter(yearWeek < filter.week) %>% 
      # mutate(density = Ixodes_scapularis / totalSampledArea * 1600) %>%
      pivot_wider(id_cols = c(yearWeek, time), 
                  names_from = plotID, 
                  values_from = ixodes_scapularis,
                  names_glue = "Ixodes_{plotID}") %>% 
      arrange(yearWeek) %>% 
      mutate(yearWeek = as.character(yearWeek),
             time = as.character(time))  
    
    if(is.null(data.amblyomma)) return(wide.ixodes)
  } 
  
  if(!is.null(data.amblyomma)){
    plots.amblyomma <- unique(data.ambloyoma$plotID)  
    
    wide.amblyomma <- data.ambloyoma %>% 
      filter(plotID %in% plots.amblyomma) %>% 
      filter(Year >= 2016) %>%
      filter(yearWeek < filter.week) %>% 
      # mutate(density = amblyomma_americanum / totalSampledArea * 1600) %>% 
      pivot_wider(id_cols = c(yearWeek, time), 
                  names_from = plotID, 
                  values_fn = {max},
                  values_from = amblyomma_americanum,
                  names_glue = "Amblyomma_{plotID}") %>% 
      arrange(yearWeek) %>% 
      mutate(yearWeek = as.character(yearWeek),
             time = as.character(time))  
    
    if(is.null(data.ixodes)) return(wide.amblyomma)
  }
  
  if(!is.null(data.ixodes) & !is.null(data.amblyomma)){
    wide.data <- full_join(wide.ixodes, wide.amblyomma,
                           by = c("yearWeek", "time")) %>% 
      arrange(yearWeek) %>% 
      mutate(yearWeek = as.character(yearWeek),
             time = as.character(time))  
    
    return(wide.data)
  }
}
