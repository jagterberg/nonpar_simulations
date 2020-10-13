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



# Simulate two SBMs and test
ntimes <- 100
results <- list()
set.seed(1133) #1111 and #1112 is okay
n <- 300

Q1 <- diag(1,2)
Q2 <- diag(-1,2)
Q3 <- diag(c(1,-1),2)
Q4 <- diag(c(-1,1),2)
signs <- list(Q1,Q2,Q3,Q4)#,Q5,Q6,Q7,Q8)
signs2 <- list(as.matrix(1),as.matrix(-1))
m <- n
B <- matrix(c(.5,.8,.8,.8,.5,.8,.8,.8,.5),3,3)
nus <- eigen(B)
nus_true <- nus$vectors %*% diag(abs(nus$values)^(1/2),3,3)
Ipq <- diag(c(1,-1,-1),3,3)
pis <- c(1/3,1/3,1/3)
sigma <- 1/2
#n <- 100


for (i in c(1:ntimes)) {
  assignmentvector1 <- rmultinom(n,1,pis)
  assignmentvector2 <- rmultinom(m,1,pis)
  Xtrue <-t(assignmentvector1) %*% nus_true
  Ytrue <-  t(assignmentvector2) %*% nus_true
  P1 <- Xtrue %*%Ipq %*% t(Xtrue)
  P2 <- Ytrue %*% Ipq %*% t(Ytrue)
  A <- generateAdjacencyMatrix(P1)
  C <- generateAdjacencyMatrix(P2)
  Xhat <- irlba(A,3)
  Xhat <- Xhat$u %*% diag(Xhat$d)^(1/2)
  Yhat <- irlba(C,3)
  Yhat <- Yhat$u %*% diag(Yhat$d)^(1/2)
  
  #find the alignment
  cs1 <- rep(0,length(signs))
  cs2 <- rep(0,length(signs))
  get_matched_1 <- list()
  get_matched_2 <- list()
  for (l in c(1:(length(signs)))) {
    get_matched_1[[l]] <- iterative_optimal_transport(Xhat,Yhat
                                                       #,lambda_init = .5
                                                       #,alpha = .5
                                                       #,lambda_final = .27
                                                       , Q = bdiag(1,signs[[l]]),numReps = 10
                                                       ,p=1,q=2)
    cs1[l] <- get_matched_1[[l]]$obj.value
    cs1[l] <- gmmase::nonpar(Xhat%*% get_matched_1[[l]]$Q,Yhat)
    
    get_matched_2[[l]] <- iterative_optimal_transport(Xhat,Yhat
                                                       # ,lambda_init = .5
                                                       #,alpha = .5
                                                       # ,lambda_final = .27
                                                       , Q = bdiag(-1,signs[[l]]),numReps = 10
                                                       ,p=1,q=2)
    cs2[l] <- get_matched_2[[l]]$obj.value
    cs2[l] <- gmmase::nonpar(Xhat%*% get_matched_2[[l]]$Q,Yhat)
  }
  
  minval1 <- cs1[which.min(cs1)]
  minval2 <- cs2[which.min(cs2)]
  if( minval1 < minval2) {
    final_Q <- get_matched_1[[which.min(cs1)]]$Q
    
  } else {
    final_Q <- get_matched_2[[which.min(cs2)]]$Q
  }
  
  Xnew <- Xhat%*% final_Q
  
  results[[i]] <- nonpar.test(Xnew,Yhat)
}

save(results,file = "results_10-13.Rdata")
