
#if(!require(devtools)) {
#  install.packages("devtools")
#  library(devtools)
#}

if(!require(nonparGraphTesting)) {
  install.packages("nonparGraphTesting_0.1.0.tar.gz", repos = NULL, type="source")
  library(nonparGraphTesting)
  
}
if(!require(hitandrun)) {
  install.packages("hitandrun")
  library(hitandrun)
}

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

# functions to generate the graph given an arbitrary P matrix
# generateAdjacency simulates according to a P matrix
# generateMatrix creates the P matrix given b, the n dimensional assignment vector
# n, and bmatrix the B matrix


Rcpp::cppFunction("
  NumericMatrix generateAdjacencyMatrix(NumericMatrix pMatrix) {
  
    int n = pMatrix.cols();
    NumericMatrix A(n,n);
    for(int i = 0; i < n; i ++) {
      for (int j = i + 1; j < n; j++) {
        A(i,j) = (int)(rand()%100 < (pMatrix(i,j)* 100));
        A(j,i) = A(i,j);
      }
    }
    return A;
  }
")
set.seed(1136)
MCs <- 50
epsilons <- c(0,.1,.2,.3,.4)
ns <- c(500,1000,1500,2000)
vals <- list()
i <- 1
Q1 <- diag(1,2)
Q2 <- diag(-1,2)
Q3 <- diag(c(1,-1),2)
Q4 <- diag(c(-1,1),2)
signs <- list(Q1,Q2,Q3,Q4)#,Q5,Q6,Q7,Q8)
signs2 <- list(as.matrix(1),as.matrix(-1))

for (eps in epsilons) {
  vals[[i]] <- list()
  print(paste0("epsilon = ",eps))
  B <- matrix(c(.3,.8,.8,.8,.3,.8,.8,.8,.3),3,3)
  B2 <- B + diag(eps,3)
  nus <- eigen(B)
  nus_true1 <- nus$vectors %*% diag(abs(nus$values)^(1/2),3,3)
  nus <- eigen(B2)
  nus_true2 <- nus$vectors %*% diag(abs(nus$values)^(1/2),3,3)
  Ipq <- diag(c(1,-1,-1),3,3)
  pis <- c(.35,.35,.3)
  j <- 1
  for (n in ns) {
    print(paste0("n = ",n))
    vals[[i]][[j]] <- 0
    for (mc in c(1:MCs)) {
      print(paste0("Run: ",mc," out of ", MCs, " for n = ",n, ", eps = ",eps))#Matching datasets")
      m <- n
     # assignmentvector1 <- rmultinom(n,1,pis)
     # assignmentvector2 <- rmultinom(m,1,pis)
     # Xtrue <- t(assignmentvector1) %*% nus_true1
     # Ytrue <- t(assignmentvector2) %*% nus_true2
      assignmentvector1 <- simplex.sample(3,n)$samples#rmultinom(n,1,pis)
      assignmentvector2 <- simplex.sample(3,n)$samples
      Xtrue <- (assignmentvector1) %*% nus_true1
      Ytrue <- (assignmentvector2) %*% nus_true2
      P1 <- Xtrue %*%Ipq %*% t(Xtrue)
      P2 <- Ytrue %*% Ipq %*% t(Ytrue)
      A <- generateAdjacencyMatrix(P1)
      C <- generateAdjacencyMatrix(P2)
      Xhat <- irlba(A,3)
      Xhat <- Xhat$u %*% diag(Xhat$d)^(1/2)
      #Xhat <- Xhat$vectors[,c(1,n-1,n)] %*% diag(abs(Xhat$values[c(1,n-1,n)])^(1/2))
      Yhat <- irlba(C,3)
      Yhat <- Yhat$u %*% diag(Yhat$d)^(1/2)
      #Yhat <- Yhat$vectors[,c(1,n-1,n)] %*% diag(abs(Yhat$values[c(1,n-1,n)])^(1/2))
      get_matched_1 <- list()
      get_matched_2 <- list()
      
      for (l in c(1:length(signs))) {
        get_matched_1[[l]] <- match_support(Xhat[,c(2,3)],Yhat[,c(2,3)]
                                            ,lambda_init = .2
                                            ,alpha = .2
                                            , Q = signs[[l]],numReps = 0)
      }
      for (eta in c(1:2)) {
        get_matched_2[[eta]] <-match_support(as.matrix(Xhat[,1]),as.matrix(Yhat[,1])
                                             ,lambda_init = .2
                                             ,alpha = .2
                                             , Q = signs2[[eta]]
                                             , numReps = 0)
      }
      cs <- sapply(get_matched_1,`[[`,3)
      Q_neg <-  get_matched_1[[which.min(cs)]]$Q
      cs2 <- sapply(get_matched_2,`[[`,3)
      Q_pos <-get_matched_2[[which.min(cs2)]]$Q
      Xnew <- Xhat %*% diag(c(Q_pos,Q_neg),3)
      test <- nonpar.test(Xnew,Yhat,1000)
      vals[[i]][[j]] <-  sum(test$permutation_results) + vals[[i]][[j]]
      
      
    }
    vals[[i]][[j]] <- vals[[i]][[j]]/MCs
    j <- j + 1
  }
  names(vals[[i]]) <- ns
  
  i <- i + 1
}
#n <- 100

results <- matrix(0,length(epsilons),length(ns))
colnames(results) <- ns
rownames(results) <- epsilons


results2 <- t(results)
# rm(results)
for (eps in epsilons) {
  for(n in ns) {
    results2[as.character(n),as.character(eps)] <- vals[[as.character(eps)]][[as.character(n)]]
  }
}

#names(vals) <- ns  

save(results2,file = "MC_results_10-2.Rdata")
