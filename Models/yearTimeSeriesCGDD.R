## basic code for fitting all data at once
# treating years as seperate time series 

library(nimble)


model <- nimbleCode({
  
  # priors
  phi ~ dbeta(10, 1) # survival
  # theta ~ dbeta(1, 10) # transition / death
  # alpha ~ dnorm(0, tau = 0.1) # intercept; observation model
  tau.beta.species ~ dgamma(0.5, 1)  # across site variability in detection
  beta ~ dnorm(0, tau = 0.01) 
  # life.constraint ~ dconstraint(phi + theta < 1) # probabilities must sum to one
  tau.obs ~ dgamma(0.5, 1) # observation error
  tau.process ~ dgamma(0.5, 1) # process error

  
  # data model, density so using dnorm
  for(p in 1:n.plots){
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        logit(pi[k, t, p]) <- beta + beta.species[species[p]]*cgdd[k, t, site[p]]
        ex.z[k, t, p] <- pi[k, t, p]*z[k, t, p]
        # ex.z[k, t, p] <- z[k, t, p] + beta.site[site[p]]*cgdd[k, t, site[p]]
        y[k, t, p] ~ dnorm(ex.z[k, t, p], tau = tau.obs)
      }
    }  
  }
  
  # species effect priors
  for(spp in 1:n.species){
  # alpha.species[spp] ~ dnorm(0, tau = tau.alpha.species) # random intercept by site
    beta.species[spp] ~ dnorm(0, tau = tau.beta.species)
  }
  
  # driver data model
  for(s in 1:n.sites){
    # beta.site[s] ~ dnorm(0, tau = tau.beta.site)
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        cgdd[k, t, s] ~ dunif(0, 6) # driver prior
        x.obs[k, t, s] ~ dnorm(cgdd[k, t, s], tau = x.tau[k, t, s])
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
        ex[k, t, p] <- phi*z[k, t-1, p] #- theta*z[k, t-1, p]
        x[k, t, p] ~ dnorm(ex[k, t, p], tau = tau.process)
        z[k, t, p] <- max(0, x[k, t, p])  
        
      } # weeks
    } # years
  } # plots
})

monitor <- c("tau.process", 
             "tau.obs", 
             "tau.beta.species",
             # "theta",
             "phi",
             # "alpha",
             # "alpha.species",
             "beta", 
             "beta.species",
             # "cgdd",
             "x")

inits <- function(){list(
  phi = runif(1, 0.95, 1),
  # theta = runif(1, 0, 0.001),
  # pi = array(runif(total.mu.index, 0.4, 0.6), dim = dim(y)),
  # alpha = rnorm(1, 0, 0.01),
  # alpha.species = rnorm(n.species, 0, 0.01),
  beta = runif(1, 1, 2),
  beta.species = runif(n.species, 0.9, 1.1),
  y = mu.inits,
  z = mu.inits,# + runif(total.mu.index, 0, 5),
  x = mu.inits,# + runif(total.mu.index, 0, 5),
  ex = mu.inits,# + runif(total.mu.index, 0, 5),
  ex.z = mu.inits,# + runif(total.mu.index, 6, 8),
  cgdd = met.inits,
  tau.beta.site = runif(1, 0, 0.001),
  tau.process = runif(1, 0, 1/2000^2), 
  tau.obs = rnorm(1, 150, 5)
)}


