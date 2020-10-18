#source("./balanced_sbm/sbm_hyp_test.R")
#epsilons <- c(0,.1,.2)
ns <- c(300,600,900)
#results_sbm <- list()

#j <- 1
#for (n in ns) {
#  print(paste0("n = ",n))
#  results_sbm[[j]] <- list()
#  i <- 1
#  for (eps in epsilons) {
#    print(paste0("eps = ",eps))
#    results_sbm[[j]][[i]] <- run_simulation_sbm(eps = eps,ntimes = 500,n=n)
#    i <- i+ 1
#  }
#  names(results_sbm[[j]]) <- epsilons
#  j <- j+1
#}
#names(results_sbm) <- ns

#save(results_sbm,file = "sbm_results_10-14.Rdata")

#print("finished SBM simulations, starting DCSBM simulations")

source("./balanced_vs_dcsbm/sbm_vs_dcsbm.R")
j <- 1
results_dcsbm <- list()
for (n in ns) {
  print(paste0("n = ",n))
  results_dcsbm[[j]] <- run_simulation_dcsbm(n=n,ntimes = 500)
  j <- j+1
}
names(results_dcsbm) <- ns
save(results_dcsbm,file = "dcsbm_results_10-14.Rdata")

