## basic code for fitting all data at once
# treating years as seperate time series 

library(nimble)


model <- nimbleCode({
  
  # priors
  tau.alpha.species ~ dgamma(0.5, 1)  # across species variability in survival
  # tau.beta.species ~ dgamma(0.5, 1)  # across species variability in survival slope
  alpha ~ dnorm(0, tau = 0.01) # across site mean survival
  beta ~ dnorm(0, tau = 0.01) # across site mean survival
  theta ~ dunif(0, 1) # across site mean removal 
  tau.process ~ dgamma(0.1, 50) # process error
  tau.obs ~ dgamma(0.5, 1) # observation error
  
  # data model, density so using dnorm
  for(p in 1:n.plots){
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        y[k, t, p] ~ dnorm(z[k, t, p], tau = tau.obs)
      }
    }  
  }
  
  # species effect priors
  for(spp in 1:n.species){
    alpha.species[spp] ~ dnorm(0, tau = tau.alpha.species) # random intercept by site
    beta.species[spp] ~ dnorm(0, tau = 0.01)
  }
  
  # driver data model
  for(s in 1:n.sites){
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        met.true[k, t, s] ~ dnorm(met.mu, tau = tau.met) # driver prior
        x.obs[k, t, s] ~ dnorm(met.true[k, t, s], tau = x.tau[k, t, s])
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
        
        logit(phi[k, t, p]) <- alpha +
          alpha.species[species[p]] +
          beta.species[species[p]]*(met.true[k, t, site[p]] - met.true[k, t-1, site[p]])
        
        # probabilities must sum to one
        life.constraint[k, t, p] ~ dconstraint(phi[k, t, p] + theta <= 1)
        
        ex[k, t, p] <- phi[k, t, p]*z[k, t-1, p] - theta*z[k, t-1, p]
        x[k, t, p] ~ dnorm(ex[k, t, p], tau = tau.process)
        z[k, t, p] <- max(0, x[k, t, p])  
        
      } # weeks
    } # years
  } # plots
})

monitor <- c("tau.process", 
             "tau.obs", 
             "tau.alpha.species",
             # "tau.beta.species",
             "theta",
             "alpha",
             # "beta", 
             "alpha.species",
             "beta.species",
             # "met.true",
             "x")


inits <- function(){list(
  phi = array(runif(total.mu.index, 0.99, 1), dim = dim(y)),
  theta = runif(1, 0, 0.001),
  alpha.species = rnorm(n.species, 0, 0.01),
  beta.species = rnorm(n.species, 0, 0.01),
  alpha = rnorm(1, 0, 0.01),
  # beta = rnorm(1, 0, 0.01),
  y = mu.inits + runif(total.mu.index, 0, 5),
  z = mu.inits + runif(total.mu.index, 0, 5),
  x = mu.inits + runif(total.mu.index, 0, 5),
  ex = mu.inits + runif(total.mu.index, 0, 5),
  met.true = met.inits,
  tau.process = runif(1, 0, 1/200^2), 
  tau.alpha.species = runif(1, 0, 0.001),
  # tau.beta.species = runif(1, 0, 0.001),
  tau.obs = rnorm(1, 150, 5)
)}





