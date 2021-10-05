## basic code for fitting all data at once
# treating years as seperate time series 
# observation probability driven by cumulative gdd
# survival rate driven by temperature
# observation error: dpois
# process error: dnorm

library(nimble)


model <- nimbleCode({
  
  # priors
  # phi ~ dbeta(10, 1) # survival
  
  
  beta ~ dnorm(beta.mu, tau = beta.tau) # intercept; observation model
  psi ~ dnorm(psi.mu, tau = psi.tau) # intercept; survival
  tau.beta.site ~ dgamma(shape = tau.beta.site.shape, 
                         rate = tau.beta.site.rate)  # across site variability in detection
  tau.psi.site ~ dgamma(shape = tau.psi.site.shape, 
                        rate = tau.psi.site.rate)  # across site variability in survival
  tau.process ~ dgamma(shape = tau.process.shape, 
                       rate = tau.process.rate) # process error

  
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
  for(spp in 1:n.species){
  # alpha.species[spp] ~ dnorm(0, tau = tau.alpha.species) # random intercept by site
    # beta.species[spp] ~ dnorm(beta, tau = tau.beta.species)
    theta[spp] ~ dbeta(theta.alpha[spp], theta.beta[spp])
  }
  
  # driver data model
  for(s in 1:n.sites){
    beta.site[s] ~ dnorm(mean = beta.site.mu[s], tau = tau.beta.site)
    psi.site[s] ~ dnorm(mean = psi.site.mu[s], tau = tau.psi.site)
    # theta[s] ~ dbeta(1, 10)
    for(k in 1:n.years){
      for(t in 1:n.weeks[k]){
        cgdd[k, t, s] ~ dnorm(x.obs.mu, tau = x.obs.tau) # driver prior
        x.obs[k, t, s] ~ dnorm(cgdd[k, t, s], tau = x.tau[k, t, s])
        
        x.phi[k, t, s] ~ dnorm(proc.met.mu, tau = proc.met.tau)
        x.obs.phi[k, t, s] ~ dnorm(x.phi[k, t, s], tau = x.phi.tau[k, t, s])
      }
    }
  }
  
  # process model
  for(p in 1:n.plots){
    for(k in 1:n.years){
      
      # first latent state of each year gets it's own prior
      x[k, 1, p] ~ T(dnorm(10, tau = 1/10), 0, Inf)
      # z[k, 1, p] <- max(0, x[k, 1, p])
      
      for(t in 2:n.weeks[k]){
        ex[k, t, p] <- x[k, t-1, p] + 
          phi[k, t, site[p]]*x[k, t-1, p] - 
          theta[species[p]]*x[k, t-1, p]
        x[k, t, p] ~ T(dnorm(ex[k, t, p], tau = tau.process), 0, Inf) 
        
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




