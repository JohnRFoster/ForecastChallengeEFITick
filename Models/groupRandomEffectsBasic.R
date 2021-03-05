# model to fit all data at once
# control flow to decide on basic random effects
#       year, species, or site

# data model has a seperate observation error for each species



model <- nimbleCode({
  
  # priors
  tau.process ~ dgamma(0.5, 1) # process error
  for(n.tau.obs in 1:n.species){
    tau.obs[n.tau.obs] ~ dgamma(0.5, 1) # observation error  
  }
  
  mu ~ dnorm(0, tau = 1/100) # global mean 
  
  if(species.effect){
    tau.species ~ dgamma(0.5, 1)
    for(spp in 1:n.species){
      alpha.species[spp] ~ dnorm(0, tau = tau.species)  
    }  
  } else if(site.effect){
    tau.site ~ dgamma(0.5, 1)
    for(se in 1:n.sites){
      alpha.site[se] ~ dnorm(0, tau = tau.site)
    }
  } else if(year.effect){
    tau.year ~ dgamma(0.5, 1)
    for(ye in 1:n.years){
      alpha.year[ye] ~ dnorm(0, tau = tau.year)
    }
  }
  
  
  # data model, density so using dnorm
  for(p in 1:n.plots){
    for(t in 1:n.weeks){
      y[t, p] ~ dnorm(z[t, p], tau = tau.obs[species[p]])
    }
  }
  
  # process model
  for(p in 1:n.plots){
    
    # first state gets it's own prior
    x[1, p] ~ dnorm(x.ic[p], tau = 10)
    z[1, p] <- max(0, x[1, p])
    
    for(t in 2:n.weeks){
      if(species.effect){
        ex[t, p] <- mu + alpha.species[species[p]]
      } else if(site.effect){
        ex[t, p] <- mu + alpha.site[site[p]]
      } else if(year.effect){
        ex[t, p] <- mu + alpha.year[year[t]]
      }
      
      x[t, p] ~ dnorm(ex[t, p], tau = tau.process)
      z[t, p] <- max(0, x[t, p])
      
    } # time
  } # plots
})