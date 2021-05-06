## basic code for fitting all data at once
# treating years as seperate time series 
# observation probability driven by cumulative gdd
# survival rate driven by temperature

library(nimble)


model <- nimbleCode({
  
  # priors
  # phi ~ dbeta(10, 1) # survival
  theta ~ dbeta(1, 10) # transition / death
  
  beta ~ dnorm(beta.mu, tau = beta.tau) # intercept; observation model
  psi ~ dnorm(0, tau = 0.1) # intercept; survival
  tau.beta.site ~ dgamma(shape = tau.beta.site.shape, rate = tau.beta.site.rate)  # across site variability in detection
  tau.psi.site ~ dgamma(0.5, 1)  # across site variability in survival
  tau.process ~ dgamma(shape = tau.process.shape, rate = tau.process.rate) # process error

  
  # data model, density so using dnorm
  for(p in 1:n.plots){
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        
        logit(phi[k, t, p]) <- psi + psi.site[site[p]] * x.phi[k, t, site[p]]
        
        logit(pi[k, t, p]) <- beta + beta.site[site[p]] * cgdd[k, t, site[p]]
        ex.z[k, t, p] <- x[k, t, p] * pi[k, t, p]
        
        y[k, t, p] ~ dpois(ex.z[k, t, p])
      }
    }  
  }
  
  # species effect priors
  # for(spp in 1:n.species){
  # # alpha.species[spp] ~ dnorm(0, tau = tau.alpha.species) # random intercept by site
  #   beta.species[spp] ~ dnorm(beta, tau = tau.beta.species)
  # }
  
  # driver data model
  for(s in 1:n.sites){
    beta.site[s] ~ dnorm(mean = beta.site.mu[s], tau = tau.beta.site)
    psi.site[s] ~ dnorm(0, tau = tau.psi.site)
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        cgdd[k, t, s] ~ dunif(0, 6) # driver prior
        x.obs[k, t, s] ~ dnorm(cgdd[k, t, s], tau = x.tau[k, t, s])
        
        x.phi[k, t, s] ~ dunif(-4, 8)
        x.obs.phi[k, t, s] ~ dnorm(x.phi[k, t, s], tau = x.phi.tau[k, t, s])
      }
    }
  }
  
  # process model
  for(p in 1:n.plots){
    for(k in 1:n.years){
      
      # first latent state of each year gets it's own prior
      # x[k, 1, p] ~ T(dnorm(x.ic[k, p], tau = 1/10), 0, Inf)
      # z[k, 1, p] <- max(0, x[k, 1, p])
      
      for(t in 2:n.weeks[k]){
        ex[k, t, p] <- x[k, t-1, p] + phi[k, t, site[p]]*x[k, t-1, p] - theta*x[k, t-1, p]
        x[k, t, p] ~ T(dnorm(ex[k, t, p], tau = tau.process), 0, Inf)
        # z[k, t, p] <- max(0, x[k, t, p])  
        
      } # weeks
    } # years
  } # plots
})

monitor <- c("tau.process", 
             "tau.psi.site",
             "tau.beta.site",
             "psi",
             "psi.site",
             "beta",
             "beta.site",
             "theta",
             "x")

inits <- function(){list(
  # phi = array(runif(total.mu.index, 0.95, 1), dim = dim(y)),
  # pi = array(runif(total.mu.index, 0.4, 0.6), dim = dim(y)),
  theta <- runif(1, 0, 0.1),
  beta = runif(1, 1, 2),
  psi = runif(1, 1, 2),
  beta.site = runif(n.sites, 0.9, 1.1),
  psi.site = runif(n.sites, 0.9, 1.1),
  y = round(mu.inits),
  z = round(mu.inits),# + runif(total.mu.index, 0, 5),
  x = round(mu.inits),# + runif(total.mu.index, 0, 5),
  ex = round(mu.inits),# + runif(total.mu.index, 0, 5),
  ex.z = round(mu.inits),# + runif(total.mu.index, 6, 8),
  cgdd = met.inits,
  tau.beta.site = runif(1, 0.3, 0.4),
  tau.psi.site = runif(1, 0, 0.001),
  tau.process = runif(1, 0, 1/2000^2) 
)}


