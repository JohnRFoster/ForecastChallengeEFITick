## basic code for fitting all data at once
# treating years as seperate time series 

library(nimble)


model <- nimbleCode({
  
  # priors
  phi ~ dbeta(10, 1) # survival
  theta ~ dbeta(1, 10) # transition / death
  life.constraint ~ dconstraint(phi + theta < 1) # probabilities must sum to one
  tau.obs ~ dgamma(0.5, 1) # observation error
  tau.process ~ dgamma(0.5, 1) # process error
  # for(n.tau.obs in 1:n.sites){
  #   tau.obs[n.tau.obs] ~ dgamma(0.5, 1) # observation error
  # }
  
  # data model, density so using dnorm
  for(p in 1:n.plots){
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        y[k, t, p] ~ dnorm(z[k, t, p], tau = tau.obs[site[p]])
      }
    }  
  }
  
  # process model
  for(p in 1:n.plots){
    for(k in 1:n.years){
    
      # first latent state of each year gets it's own prior
      x[k, 1, p] ~ dnorm(x.ic[k, p], tau = 1/10)
      z[k, 1, p] <- max(0, x[k, 1, p])
      
      for(t in 2:n.weeks[k]){
        ex[k, t, p] <- phi*z[k, t-1, p] - theta*z[k, t-1, p]
        x[k, t, p] ~ dnorm(ex[k, t, p], tau = tau.process)
        z[k, t, p] <- max(0, x[k, t, p])  
        
      } # weeks
    } # years
  } # plots
})



# two states
# z = latent: the true number of ticks, instead of a year effect 
# y = observed: what we see






