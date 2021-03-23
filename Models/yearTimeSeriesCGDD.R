## basic code for fitting all data at once
# treating years as seperate time series 

library(nimble)


model <- nimbleCode({
  
  # priors
  phi ~ dbeta(10, 1) # survival
  theta ~ dbeta(1, 10) # transition / death
  alpha ~ dnorm(0, tau = 0.1) # intercept; observation model
  beta ~ dnorm(0, tau = 0.1) # slope; observation model
  life.constraint ~ dconstraint(phi + theta < 1) # probabilities must sum to one
  tau.obs ~ dgamma(0.5, 1) # observation error
  tau.process ~ dgamma(0.5, 1) # process error

  
  # data model, density so using dnorm
  for(p in 1:n.plots){
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        logit(pi[k, t, p]) <- alpha + beta*cgdd[k, t, site[p]]
        ex.z[k, t, p] <- pi*z[k, t, p]
        y[k, t, p] ~ dnorm(ex.z[k, t, p], tau = tau.obs)
      }
    }  
  }
  
  # driver data model
  for(s in 1:n.sites){
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        cgdd[k, t, s] ~ dnorm(met.mu, tau = tau.met) # driver prior
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
        ex[k, t, p] <- phi*z[k, t-1, p] - theta*z[k, t-1, p]
        x[k, t, p] ~ dnorm(ex[k, t, p], tau = tau.process)
        z[k, t, p] <- max(0, x[k, t, p])  
        
      } # weeks
    } # years
  } # plots
})

monitor <- c("tau.process", 
             "tau.obs", 
             "theta",
             "phi",
             "pi",
             "alpha",
             "beta", 
             # "cgdd",
             "x")

inits <- function(){list(
  phi = runif(1, 0.95, 1),
  theta = runif(1, 0, 0.001),
  pi = runif(1, 0.4, 0.6),
  alpha = rnorm(1, 0, 0.01),
  beta = rnorm(1, 0, 0.01),
  y = mu.inits + runif(total.mu.index, 0, 5),
  z = mu.inits + runif(total.mu.index, 0, 5),
  x = mu.inits + runif(total.mu.index, 0, 5),
  ex = mu.inits + runif(total.mu.index, 0, 5),
  ex.z = mu.inits + runif(total.mu.index, 6, 8),
  met.true = met.inits,
  tau.process = runif(1, 0, 1/200^2), 
  tau.obs = rnorm(1, 150, 5)
)}


