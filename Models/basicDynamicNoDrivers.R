
# ====================================================== #
#        basic dynamic model (b + rho*[x-1])             #
#        Using a normal-normal state-space model         #
#        Using tick density per trapping (1600^m-2)      #
# ====================================================== #

model <- nimbleCode({
  
  # priors
  sd.process ~ dgamma(10, 1) # process error
  sd.obs ~ dgamma(10, 1) # observation error
  beta ~ dnorm(0, sd = 5) # intercept
  rho ~ dnorm(0, sd = 5) # slope
  
  # data model, density so using dnorm
  for(t in 1:n.weeks){
    y[t] ~ dnorm(z[t], sd = sd.obs)
  }
  
  # process model
  # z[1] ~ dnorm(x.ic, sd = 10) # first state gets it's own prior
  for(t in 2:n.weeks){
    ex[t-1] <- beta + rho*z[t-1] 
    x[t] ~ dnorm(ex[t-1], sd = sd.process)
    z[t] <- max(0, x[t])
  }
})

monitor <- c("sd.process", "sd.obs", "beta", "rho", "z")

inits <- function() {
  
  # mu.inits <- y
  # mu.inits[is.na(mu.inits)] <- round(approx(mu.inits, xout = which(is.na(mu.inits)))$y)
  
  list(
    beta = mean(y, na.rm = TRUE),
    z = pmax(0, jitter(mu.inits)) + 1,
    y = mu.inits,
    x = pmax(0, jitter(mu.inits)),
    sd.process = abs(rnorm(1, mean(mu.inits, na.rm = TRUE), 10)),
    sd.obs = rgamma(1, 5, 1),
    beta = jitter(mean(y, na.rm = TRUE)),
    rho = rnorm(1, 0, 3)
  )
}
