## get target weeks for the forecast

target_weeks <- function(day.run){
  # set up target dates
  start.epi.weeks <- c(10, 14, 19, 23, 28, 32, 36, 41) # all possible start weeks
  
  # anytime we run this script before the start of the challenge we want to forecast all 2019 target weeks
  if(day.run < "2021-03-31"){ 
    start.week <- start.epi.weeks[1]
  } else { # otherwise use the appropriate starting week (months are 2 ahead)
    start.week <- start.epi.weeks[month(day.run) - 2]
  }
  
  end.week <- 44 # does not change
  target.weeks <- start.week:end.week
  
  return(target.weeks)
}