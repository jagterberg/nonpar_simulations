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

# numcores <- detectCores()
# registerDoParallel(cores=3)
epsilons <- c(0,.1,.2)
ns <- c(100,200,300,400,500,600,700)#600,900)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


#cl <- makeCluster() #not to overload your computer
#registerDoParallel(cl)
print(paste0("packages loaded, running SBM simulation"))

results <- foreach(n = ns,.packages=c('nonparGraphTesting','irlba','igraph','Rcpp','Matrix')
        ,.noexport = "generateAdjacencyMatrix" )  %:% 

            foreach(eps=epsilons,.packages=c('nonparGraphTesting','irlba','igraph','Rcpp','Matrix')
                                     ,.noexport = "generateAdjacencyMatrix" )  %dopar% {
                                       source("./balanced_sbm/sbm_hyp_test.R")
                                       #print(paste("eps = ",eps,", n = ",n))
                                       run_simulation_sbm(eps = eps,ntimes = 100,n=n,nMC = 500,seed=13)
        }
          
        


stopCluster(cl)
save(results,file = "sbm_results_11-3.Rdata")
print("finished.")#,ns[1]))
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# results_sbm_500 <- foreach(eps=epsilons,.packages=c('nonparGraphTesting','irlba','igraph','Rcpp','Matrix')
#                             ,.noexport = "generateAdjacencyMatrix" )  %dopar% {
#                               source("./balanced_sbm/sbm_hyp_test.R")
#                               #print(paste("eps = ",eps,", n = ",n))
#                               run_simulation_sbm(eps = eps,ntimes = 100,n=ns[2],nMC = 500)
#                             }
# stopCluster(cl)
# save(results_sbm_500,file = "sbm_results_10-21_500.Rdata")
# print(paste("finished n =",ns[2]))
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# results_sbm_700 <- foreach(eps=epsilons,.packages=c('nonparGraphTesting','irlba','igraph','Rcpp','Matrix')
#                             ,.noexport = "generateAdjacencyMatrix" )  %dopar% {
#                               source("./balanced_sbm/sbm_hyp_test.R")
#                               #print(paste("eps = ",eps,", n = ",n))
#                               run_simulation_sbm(eps = eps,ntimes = 100,n=ns[3],nMC = 500)
#                             }
# 
# 
# save(results_sbm_700,file = "sbm_results_10-21_700.Rdata")
# stopCluster(cl)
# results_sbm <- list(results_sbm_300,results_sbm_500,results_sbm_700)
# save(results_sbm,file="sbm_results_10-21.Rdata")
# print("finished SBM simulations.")#, starting DCSBM simulations")
# 
