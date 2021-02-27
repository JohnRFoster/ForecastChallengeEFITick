

run_nimble_parallel <- function(cl, model, constants, data, inits, use.dzip = TRUE, n.ens = 5000,
                                n.iter = 50000, n.burnin = 5000, psrf.max = 1.1, max.iter = 2e10){
  library(parallel)
  library(nimble)
  library(coda)
  
  # message("Building cluster")
  # cl <- makeCluster(n.cores)
  
  if(use.dzip){
    source("Functions/ZIP.R")
    assign("dZIP", dZIP, envir = .GlobalEnv)
    clusterExport(cl, 
                  c("model", "constants", "data", "inits", "n.iter", "n.burnin", "dZIP"),
                  envir = environment())  
  } else {
    clusterExport(cl, 
                  c("model", "constants", "data", "inits", "n.iter", "n.burnin"),
                  envir = environment()) 
  }
  
  
  for(j in seq_along(cl)){
    set.seed(j)
    init <- inits()
    clusterExport(cl[j], "init", envir = environment())
  }
  
  message("Running mcmc")
  out <- clusterEvalQ(cl, {
    library(nimble)
    library(coda)
    model.rw <- nimbleModel(model,
                            constants = constants,
                            data = data,
                            inits = init)
    cModel.rw <- compileNimble(model.rw)
    mcmcConf <- configureMCMC(cModel.rw)
    mcmcBuild <- buildMCMC(mcmcConf)
    compMCMC <- compileNimble(mcmcBuild)
    out.1 <- runMCMC(compMCMC, niter = n.iter, nburnin = n.burnin)
    return(as.mcmc(out.1))
  })
  out.mcmc <- as.mcmc.list(out)
  message("Checking convergence") 
  g.diag <- gelman.diag(out.mcmc, multivariate = FALSE)$psrf
  convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
  message(paste("Convergence:", convergence)) 
  
  counter <- 1
  total.iter <- counter * n.iter
  while(!convergence | total.iter > max.iter){
    counter <- counter + 1
    total.iter <- counter * n.iter
    
    message("Running mcmc")
    out2 <- clusterEvalQ(cl, {
      out1 <- runMCMC(compMCMC, niter = n.iter)
      return(as.mcmc(out1))
    })
    out.mcmc.update1 <- as.mcmc(out2)
    
    out.mcmc.bind <- list()
    for (i in seq_len(n.cores)) {
      out.mcmc.bind[[i]] <- mcmc(rbind(out.mcmc[[i]], out.mcmc.update1[[i]]))
    }  
    out.mcmc <- out.mcmc.bind 
    message("Checking convergence") 
    g.diag <- gelman.diag(out.mcmc, multivariate = FALSE)$psrf  
    convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
    message(paste("Convergence:", convergence))  
    message(paste("Total iterations:", total.iter))
    # print(tail(g.diag))
  }
  
  # need to at effective sample size calculation
  
  GBR <- gelman.plot(out.mcmc, ask = FALSE)
  burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
  message(paste("Burnin after", burnin, "iterations"))
  
  # sometimes burnin happens at the end of the chain, need to keep sampling if that happens
  if(total.iter - burnin < 1000){
    message("Additional iterations because burnin is too close to the end of the chain")
    out3 <- clusterEvalQ(cl, {
      out <- runMCMC(compMCMC, niter = n.iter/2)
      return(as.mcmc(out))
    })
    samples <- as.mcmc.list(out)
  } else{
    samples <- window(as.mcmc.list(out.mcmc), start = burnin)
  }  
  
  # return a thinned matrix (raw mcmc objects can get big)
  samples <- as.matrix(samples)
  thin.seq <- round(seq(1, nrow(samples), length.out = n.ens)) # sequence of samples to keep
  samples <- samples[thin.seq,]
  
  return(samples)
}
