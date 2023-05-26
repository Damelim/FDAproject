library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(ggplot2)

cppFunction("arma::mat armaInv(const arma::mat & x) {return arma::inv(x);}", depends = "RcppArmadillo")

sourceCpp('G:/My Drive/research/C_functions/MatMult.cpp')

# cppFunction(
# "SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
#   Eigen::MatrixXd C = A * B;
#   
#   return Rcpp::wrap(C);
# }"
# )




periodic_kernel = function(x1, x2, k){
  exp(sin(k * pi * abs(x1 - x2))^2)
}

gaussian_kernel <- function(x1, x2, sigma) {
  exp(-0.5 * sum((x1 - x2)^2) / sigma^2)
} # same as rbf kernel
# rbf_kernel = function(x1, x2, gamma){
#   exp(-gamma * sum((x1 - x2)^2))
# }
laplacian_kernel = function(x1, x2, gamma){
  exp(-gamma * sum(abs(x1 - x2)))
}

polynomial_kernel <- function(x1, x2, degree = 2) {
  (1 + sum(x1 * x2))^degree
}

sigmoid_kernel = function(x1,x2, gamma){
  tanh(gamma * sum(x1 * x2) + 1)
}


generate_kernel_matrix = function(x, kernel_type, params){
  
  n = length(x)
  K = matrix(0, n, n)
  
  for(i in 2:n){
    for(j in 1:(i-1)){
      if(kernel_type == "gaussian"){
        K[i,j] = gaussian_kernel(x[i], x[j], params)
      }else if(kernel_type == "polynomial"){
        K[i,j] = polynomial_kernel(x[i], x[j], params)
      }else if(kernel_type == "sigmoid"){
        K[i,j] = sigmoid_kernel(x[i], x[j], params)
      }else if(kernel_type == "laplacian"){
        K[i,j] = laplacian_kernel(x[i], x[j], params)
      }else if(kernel_type == "periodic"){
        K[i,j] = periodic_kernel(x[i], x[j], params)
      }
    }
  }
  
  K = K + t(K)
  
  for(i in 1:n){
    if(kernel_type == "gaussian"){
      K[i,i] = gaussian_kernel(x[i], x[i], params)
    }else if(kernel_type == "polynomial"){
      K[i,i] = polynomial_kernel(x[i], x[i], params)
    }else if(kernel_type == "sigmoid"){
      K[i,i] = sigmoid_kernel(x[i], x[i], params)
    }else if(kernel_type == "laplacian"){
      K[i,i] = laplacian_kernel(x[i], x[i], params)
    }else if(kernel_type == "periodic"){
      K[i,i] = periodic_kernel(x[i], x[i], params)
    }
  }
  return(K)
}




generate_prediction_kernel_matrix = function(trainx, testx, kernel_type, params){
  
  n = length(trainx)
  m = length(testx)
  
  K_pred = matrix(0, m, n)
  
  for(i in 1:m){
    for(j in 1:n){
      if(kernel_type == "gaussian"){
        K_pred[i,j] = gaussian_kernel(testx[i], trainx[j], params)
      }else if(kernel_type == "polynomial"){
        K_pred[i,j] = polynomial_kernel(testx[i], trainx[j], params)
      }else if(kernel_type == "sigmoid"){
        K_pred[i,j] = sigmoid_kernel(testx[i], trainx[j], params)
      }else if(kernel_type == "laplacian"){
        K_pred[i,j] = laplacian_kernel(testx[i], trainx[j], params)
      }else if(kernel_type == "periodic"){
        K_pred[i,j] = periodic_kernel(testx[i], trainx[j], params)
      }
    }
  }
  
  return(K_pred)
}










# K = generate_kernel_matrix(x, "gaussian", 0.5)
# K = generate_kernel_matrix(x, "polynomial", 3)

# KernelRidge = function(x, y, kernel_type, params, lamvec, trueftn){
#   
#   if((kernel_type == "gaussian") || (kernel_type == "laplacian") || (kernel_type == "sigmoid")){
#     plot(x,y, main = paste("n = ",length(x),', ', kernel_type,' kernel with parameter ',params, sep = ""))
#   }else if(kernel_type == "polynomial"){
#     plot(x,y, main = paste("n = ",length(x),', ', kernel_type,' kernel with degree ',params, sep = ""))
#   }
# 
#   lines(x,trueftn(x), type="l", col = "red", lwd = 4)
#   
#   K = generate_kernel_matrix(x, kernel_type, params)
#   
#   I = diag(1, length(x))
#   
#   for(i in 1:length(lamvec)){
#     
#     lambda = lamvec[i]
#     
#     alpha = eigenMapMatMult(armaInv(K + lambda*I), y)
#     
#     fittedval = eigenMapMatMult(K, alpha)
#     
#     lines(x, fittedval, col = i, lwd = 2, lty = 2)
#   }
#   legend("topright", legend = c('true',paste('lam = ', lamvec)), col = c('red', 1:length(lamvec)), lwd = 3)
# }


KernelRidge = function(x, y, kernel_type, params, lamvec, trueftn, return = FALSE){
  
  colvec = c('red','red3','red4','rosybrown','rosybrown2','rosybrown3','rosybrown4','royalblue','royalblue2','royalblue3','royalblue4')
  
  if((kernel_type == "gaussian") || (kernel_type == "laplacian") || (kernel_type == "sigmoid")){
    plot(x,y, cex = 0.4, main = paste("n = ",length(x),', ', kernel_type,' kernel with parameter ',params, sep = ""))
  }else if(kernel_type == "polynomial"){
    plot(x,y, cex = 0.4, main = paste("n = ",length(x),', ', kernel_type,' kernel with degree ',params, sep = ""))
  }
  
  lines(x,trueftn(x), type="l", col = "black", lwd = 5)
  
  K = generate_kernel_matrix(x, kernel_type, params)
  
  I = diag(1, length(x))
  
  fittedval_df = matrix(NA, length(x), length(lamvec))
  
  mny = mean(y) ; y_centered = y - mny
  
  for(i in 1:length(lamvec)){
    
    lambda = lamvec[i]
    
    alpha = eigenMapMatMult(armaInv(K + lambda*I), y_centered)
    
    fittedval = mny + eigenMapMatMult(K, alpha)
    
    fittedval_df[,i] = fittedval
    
    lines(x, fittedval, col = colvec[i], lwd = 2.5, lty = 2)
  }
  legend("topright", legend = c('true',paste('lam = ', lamvec)), col = c('black', colvec[1:length(lamvec)]), lwd = 3)
  
  if(return == TRUE){
    return(fittedval_df)
  }
}





