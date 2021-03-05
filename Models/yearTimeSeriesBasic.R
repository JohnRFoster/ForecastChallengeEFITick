## basic code for fitting all data at once
# treating years as seperate time series 

library(nimble)


model <- nimbleCode({
  
  # priors
  rho ~ dnorm(0, tau = 1/100)
  tau.process ~ dgamma(0.5, 1) # process error
  for(n.tau.obs in 1:n.species){
    tau.obs[n.tau.obs] ~ dgamma(0.5, 1) # observation error  
  }
  
  # data model, density so using dnorm
  for(k in 1:n.years){
    for(p in 1:n.plots){
      for(t in 1:n.weeks[k]){
        y[k, t, p] ~ dnorm(z[year[t], t, p], tau = tau.obs[species[p]])
      }
    }  
  }
  
  # process model
  for(k in 1:n.years){
    for(p in 1:n.plots){
      
      # first latent state of each year gets it's own prior
      x[k, 1, p] ~ dnorm(x.ic[k, p], tau = 10)
      z[k, 1, p] <- max(0, x[1, 1, p])
      
      for(t in 2:n.weeks[k]){
        ex[k, t, p] <- rho*x[year[t], t-1, p]
        x[k, t, p] ~ dnorm(ex[year[t], t, p], tau = tau.process)
        z[k, t, p] <- max(0, x[year[t], t, p])  
        
      } # weeks
    } # plots
  } # years
})



# two states
# z = latent: the true number of ticks, instead of a year effect 
# y = observed: what we see






