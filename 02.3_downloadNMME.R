library(ncdf4)
library(utils)

get_nmme <- function(var, start.month, end.month,  
                     dest.dir = "Data/NMME"){
  
  for(i in 1:10){ # 10 files per var x month
    # construct file name
    nc.2.grab <- paste0(var,
                        "_day_ccsm4_",
                        start.month,
                        "_r",
                        i,
                        "i1p1_",
                        start.month,
                        "-",
                        end.month,
                        ".nc")
    
    
    url <- "https://www.ncei.noaa.gov/data/north-american-multi-model-ensemble/access/ccsm4/2019/"
    
    save.dir <- file.path(dest.dir, start.month, var)
    if(!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)
    
    destfile <- file.path(save.dir, nc.2.grab)
    if(file.exists(destfile)) next
    
    # download
    download.file(
      url = paste0(url, nc.2.grab),
      destfile = destfile
    )
    
  }
}

# variables to grab
vars <- c(
  "tasmin", # daily min temp
  "tasmax", # faily max temp
  "hus" # specific humidity
)

# start months to grab
# 20190301 does not exist
start.month <- c(
  "20190201", # feb 2019 - jan 2020
  "20190401",
  "20190501",
  "20190601",
  "20190701",
  "20190801",
  "20190901",
  "20191001"
)

# end months
end.month <- c(
  "20200131",
  "20200331",
  "20200430",
  "20200531",
  "20200630",
  "20200731",
  "20200831",
  "20200930"
)

for(v in seq_along(vars)){
  for(m in seq_along(start.month)){
    get_nmme(vars[v], start.month[m], end.month[m])
  }
}
