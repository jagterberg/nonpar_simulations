if(!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}
if(!require(parallel)) {
  install.packages("parallel")
  library(parallel)
}
if(!require(doParallel)) {
  install.packages('doParallel')
  library(doParallel)
}
install.packages("nonparGraphTesting_0.1.0.tar.gz", repos = NULL, type="source")
library(nonparGraphTesting)
if (!require(irlba)) {
  install.packages("irlba")
  library(irlba)
}
if(!require(igraph)) {
  install.packages("igraph")
  library(igraph)
}
if(!require(Rcpp)) {
  install.packages("Rcpp")
  library(Rcpp)
}
if(!require(Matrix)) {
  install.packages("Matrix")
  library(Matrix)
}

numcores <- detectCores()
registerDoParallel(cores=numcores)
epsilons <- c(0,.1,.2)
ns <- c(300,600,900)
print(paste0("packages loaded, running SBM simulation on ",numcores," cores."))

results_sbm <- list()
results_sbm <- foreach(eps=epsilons,n=ns
                       ,.packages=c('nonparGraphTesting','irlba','igraph','Rcpp','Matrix')
                       ,.noexport = "generateAdjacencyMatrix" ) %dopar% {
  source("./balanced_sbm/sbm_hyp_test.R")
  print(paste("eps = ",eps,", n = ",n))
  run_simulation_sbm(eps = eps,ntimes = 100,n=n)
  
}


# #j <- 1
# #for (n in ns) {
# #  print(paste0("n = ",n))
#   results_sbm[[j]] <- list()
#   i <- 1
#   for (eps in epsilons) {
#     print(paste0("eps = ",eps))
#     results_sbm[[j]][[i]] <- run_simulation_sbm(eps = eps,ntimes = 100,n=n)
#     i <- i+ 1
#   }
#   names(results_sbm[[j]]) <- epsilons
#   j <- j+1
# }
# names(results_sbm) <- ns

save(results_sbm,file = "sbm_results_10-19.Rdata")

print("finished SBM simulations.")#, starting DCSBM simulations")

# source("./balanced_vs_dcsbm/sbm_vs_dcsbm.R")
# j <- 1
# results_dcsbm <- list()
# for (n in ns) {
#   print(paste0("n = ",n))
#   results_dcsbm[[j]] <- run_simulation_dcsbm(n=n,ntimes = 100)
#   j <- j+1
# }
# names(results_dcsbm) <- ns
# save(results_dcsbm,file = "dcsbm_results_10-19.Rdata")

