# ====================================================== #
#        one covaraiate in process model                 #
#        driver used as beta*(x[t]-x[t-1])
#        Using a normal-normal state-space model         #
#        Using tick density per trapping (1600^m-2)      #
# ====================================================== #

model <- nimbleCode({
  
  # priors
  sd.process ~ dgamma(10, 1) # process error
  sd.obs ~ dgamma(10, 1) # observation error
  beta ~ dnorm(0, sd = 5) # intercept
  rho ~ dnorm(0, sd = 5) # slope
  # sd.x ~ gamma(10, 1) # met error
  
  # data model, density so using dnorm
  for(t in 1:n.weeks){
    y[t] ~ dnorm(z[t], sd = sd.obs) # tick data model
    
    # driver data model
    x[t] ~ dnorm(x.mu, sd = x.sd) # driver prior
    x.obs[t] ~ dnorm(x[t], var = x.var[t])
  }
  
  # process model
  # z[1] ~ dnorm(x.ic, sd = 10) # first state gets it's own prior
  for(t in 2:n.weeks){
    
    ex[t-1] <- beta*(x[t]-x[t-1]) + rho*z[t-1] 
    m[t] ~ dnorm(ex[t-1], sd = sd.process)
    z[t] <- max(0, m[t])
  }
})

monitor <- c("sd.process", "sd.obs", "beta", "rho", "z")

inits <- function() {
  
  # mu.inits <- y
  # mu.inits[is.na(mu.inits)] <- round(approx(mu.inits, xout = which(is.na(mu.inits)))$y)
  
  list(
    z = pmax(0, jitter(mu.inits)) + 1,
    y = jitter(mu.inits),
    ex = jitter(mu.inits),
    m = jitter(mu.inits),
    x = jitter(x.inits),
    x.obs = jitter(x.inits),
    sd.process = runif(1, 30, 60), 
    sd.obs = rgamma(1, 5, 1),
    beta = rnorm(1, 0, 3),
    rho = rnorm(1, 1, 0.01)
  )
}
