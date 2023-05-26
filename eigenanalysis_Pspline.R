library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(ggplot2)

cppFunction("arma::mat armaInv(const arma::mat & x) {return arma::inv(x);}", depends = "RcppArmadillo")

sourceCpp('MatMult.cpp')

get_basis_matrix = function(x, degree = 3, num_interior_knots, NCS = FALSE){
  
  knots = seq(0, 1, length = num_interior_knots + 2)
  knots = knots[-1] ; knots = knots[-length(knots)]
  
  if(NCS == TRUE){
    return(ns(x, knots = knots, degree = degree, intercept = T, Boundary.knots = c(0,1)))
  }else{
    return(bs(x, knots = knots, degree = degree, intercept= T))
  }
}

get_basis_matrix_nointercept = function(x, degree = 3, num_interior_knots, NCS = FALSE){
  
  knots = seq(0, 1, length = num_interior_knots + 2)
  knots = knots[-1] ; knots = knots[-length(knots)]
  
  if(NCS == TRUE){
    return(ns(x, knots = knots, degree = degree, intercept = F, Boundary.knots = c(0,1)))
  }else{
    return(bs(x, knots = knots, degree = degree, intercept= F))
  }
}


load('Simulationdata.rdata')

f = f8

train_sigma = train_mediumsigma # 0.4

n = 800
set.seed(1) ; x = sort(runif(n, 0, 1)) ; x_pred = seq(0.001, 0.999, by = 0.002) ## function evaluation points
trainfunction = data.frame(cbind(x, f(x))) ; colnames(trainfunction) = c('x','z')
testfunction = data.frame(cbind(x_pred, f(x_pred))) ; colnames(testfunction) = c('x','z')
test_Y_signal = testfunction$z

set.seed(1)
traindata = trainfunction ; n = nrow(traindata)
traindata$z = trainfunction$z + rnorm(n, mean = 0, sd = train_sigma)
y = traindata$z ; yty = eigenMapMatMult(t(y), y)


knot = 20

sigma2_0 = 0.01 ; nu_0 = 0.01
g_vec = c(Inf, 10^(seq(-0,-3, length = 3)))


testpredictions_list = list()
testpredictions_list[[1]] = data.frame(x = x_pred, y = test_Y_signal)

B = get_basis_matrix(x = x, degree = 3, num_interior_knots = knot)
BtB = eigenMapMatMult(t(B), B)
B_pred = get_basis_matrix(x = x_pred, degree = 3, num_interior_knots = knot)


p = ncol(B) ; I_p = diag(1,p)

P = diff(diag(1, ncol(B)), diff = 2) ; P = as.matrix(P); P = t(P) %*% P
P = P + diag(1e-10, ncol(P))

Bty = eigenMapMatMult( t(B), y )
nu0_n = nu_0 + n
nu0_sigma20 = nu_0*sigma2_0


logMSE_mat = matrix(NA, length(g_vec), 3) ; logMSE_mat = data.frame(logMSE_mat)
names(logMSE_mat) = c('knot','g','logMSE')
logMSE_mat[,1] = rep(knot, each = length(g_vec))
logMSE_mat[,2] = rep(g_vec, 1)


for(j in 1:length(g_vec)){
  
  testpredictions = matrix(NA, nrow = nrow(testfunction), ncol = n_mcmc_sample) # HAVE TO DEFINE THIS INSIDE THE LOOP TO AVOID INDEXERROR FROM BURNIN!!!
  
  g = g_vec[j]
  
  BtB_invgP = BtB + 1/g * P  
  BtB_invgP_inv = armaInv(BtB_invgP) # have to be invertible
  mn = eigenMapMatMult(BtB_invgP_inv , Bty)
  SSR = yty - eigenMapMatMult(t(y),eigenMapMatMult(B, mn)) ; SSR = SSR[1,1]
  
  for(iter in 1:n_mcmc_sample){
    
    sigma2 = rinvgamma(n = 1, shape = nu0_n/2 , scale = (nu0_sigma20 + SSR)/2)
    beta = mvnfast::rmvn(n = 1, mu = mn , sigma = sigma2 * BtB_invgP_inv) #, checkSymmetry = F)
    
    if(iter > nburnin){
      testpredictions[,iter] = eigenMapMatMult(B_pred, t(beta))
    }
  }
  
  testpredictions = testpredictions[,(nburnin+1):n_mcmc_sample] # get rid of burn in
  model_averaging = rowMeans(testpredictions) # posterior mean
  logMSE_mat[length(g_vec)+j,3] = log(mean((model_averaging - test_Y_signal)^2))
  testpredictions_list[[j+1]] = data.frame(x = x_pred, y = model_averaging)
}

Smoothermat1 = eigenMapMatMult(eigenMapMatMult(B, armaInv(BtB)), t(B))
Smoothermat2 = eigenMapMatMult(eigenMapMatMult(B, armaInv(BtB + 1/g_vec[2] * P)), t(B))
Smoothermat3 = eigenMapMatMult(eigenMapMatMult(B, armaInv(BtB + 1/g_vec[3] * P)), t(B))
Smoothermat4 = eigenMapMatMult(eigenMapMatMult(B, armaInv(BtB + 1/g_vec[4] * P)), t(B))

dim(Smoothermat1)
dim(Smoothermat2)
dim(Smoothermat3)
dim(Smoothermat4)

eig1 = eigen(Smoothermat1, symmetric = TRUE)
eig2 = eigen(Smoothermat2, symmetric = TRUE)
eig3 = eigen(Smoothermat3, symmetric = TRUE)
eig4 = eigen(Smoothermat4, symmetric = TRUE)

eval1 = eig1$values
evec1 = eig1$vectors

eval2 = eig2$values
evec2 = eig2$vectors

eval3 = eig3$values
evec3 = eig3$vectors

eval4 = eig4$values
evec4 = eig4$vectors









png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Pspline_prediction.png', width = 15, height = 8, units = "in", res = 300)

par(mfrow = c(1,3))

plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Cubic-Pspline: prediction, rank(B)=",ncol(B),', n=',n,', sigma=',train_sigma), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, testpredictions_list[[2]]$y, col = "red", lwd = 2)
lines(x_pred, testpredictions_list[[3]]$y, col = "green", lwd = 2)
lines(x_pred, testpredictions_list[[4]]$y, col = "blue", lwd = 2)
lines(x_pred, testpredictions_list[[5]]$y, col = "purple", lwd = 2)
legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)







f = f17
n = 800
set.seed(1) ; x = sort(runif(n, 0, 1)) ; x_pred = seq(0.001, 0.999, by = 0.001) ## function evaluation points
trainfunction = data.frame(cbind(x, f(x))) ; colnames(trainfunction) = c('x','z')
testfunction = data.frame(cbind(x_pred, f(x_pred))) ; colnames(testfunction) = c('x','z')
test_Y_signal = testfunction$z

set.seed(1)
traindata = trainfunction ; n = nrow(traindata)
traindata$z = trainfunction$z + rnorm(n, mean = 0, sd = train_sigma)
y = traindata$z ; yty = eigenMapMatMult(t(y), y)

knot = 20

sigma2_0 = 0.01 ; nu_0 = 0.01
g_vec = c(Inf, 10^(seq(-0,-3, length = 3)))


testpredictions_list = list()
testpredictions_list[[1]] = data.frame(x = x_pred, y = test_Y_signal)

B = get_basis_matrix(x = x, degree = 3, num_interior_knots = knot)
BtB = eigenMapMatMult(t(B), B)
B_pred = get_basis_matrix(x = x_pred, degree = 3, num_interior_knots = knot)


p = ncol(B) ; I_p = diag(1,p)

P = diff(diag(1, ncol(B)), diff = 2) ; P = as.matrix(P); P = t(P) %*% P
P = P + diag(1e-10, ncol(P))

Bty = eigenMapMatMult( t(B), y )
nu0_n = nu_0 + n
nu0_sigma20 = nu_0*sigma2_0


logMSE_mat = matrix(NA, length(g_vec), 3) ; logMSE_mat = data.frame(logMSE_mat)
names(logMSE_mat) = c('knot','g','logMSE')
logMSE_mat[,1] = rep(knot, each = length(g_vec))
logMSE_mat[,2] = rep(g_vec, 1)


for(j in 1:length(g_vec)){
  
  testpredictions = matrix(NA, nrow = nrow(testfunction), ncol = n_mcmc_sample) # HAVE TO DEFINE THIS INSIDE THE LOOP TO AVOID INDEXERROR FROM BURNIN!!!
  
  g = g_vec[j]
  
  BtB_invgP = BtB + 1/g * P  
  BtB_invgP_inv = armaInv(BtB_invgP) # have to be invertible
  mn = eigenMapMatMult(BtB_invgP_inv , Bty)
  SSR = yty - eigenMapMatMult(t(y),eigenMapMatMult(B, mn)) ; SSR = SSR[1,1]
  
  for(iter in 1:n_mcmc_sample){
    
    sigma2 = rinvgamma(n = 1, shape = nu0_n/2 , scale = (nu0_sigma20 + SSR)/2)
    beta = mvnfast::rmvn(n = 1, mu = mn , sigma = sigma2 * BtB_invgP_inv) #, checkSymmetry = F)
    
    if(iter > nburnin){
      testpredictions[,iter] = eigenMapMatMult(B_pred, t(beta))
    }
  }
  
  testpredictions = testpredictions[,(nburnin+1):n_mcmc_sample] # get rid of burn in
  model_averaging = rowMeans(testpredictions) # posterior mean
  logMSE_mat[length(g_vec)+j,3] = log(mean((model_averaging - test_Y_signal)^2))
  testpredictions_list[[j+1]] = data.frame(x = x_pred, y = model_averaging)
}


plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Cubic-Pspline: prediction, rank(B)=",ncol(B),', n=',n,', sigma=',train_sigma), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, testpredictions_list[[2]]$y, col = "red", lwd = 2)
lines(x_pred, testpredictions_list[[3]]$y, col = "green", lwd = 2)
lines(x_pred, testpredictions_list[[4]]$y, col = "blue", lwd = 2)
lines(x_pred, testpredictions_list[[5]]$y, col = "purple", lwd = 2)
legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)
















f = f12
n = 800
set.seed(1) ; x = sort(runif(n, 0, 1)) ; x_pred = seq(0.001, 0.999, by = 0.001) ## function evaluation points
trainfunction = data.frame(cbind(x, f(x))) ; colnames(trainfunction) = c('x','z')
testfunction = data.frame(cbind(x_pred, f(x_pred))) ; colnames(testfunction) = c('x','z')
test_Y_signal = testfunction$z

set.seed(1)
traindata = trainfunction ; n = nrow(traindata)
traindata$z = trainfunction$z + rnorm(n, mean = 0, sd = 0.1)
y = traindata$z ; yty = eigenMapMatMult(t(y), y)

sigma2_0 = 0.01 ; nu_0 = 0.01
g_vec = c(Inf, 10^(seq(-0,-3, length = 3)))


testpredictions_list = list()
testpredictions_list[[1]] = data.frame(x = x_pred, y = test_Y_signal)

B = get_basis_matrix(x = x, degree = 3, num_interior_knots = knot)
BtB = eigenMapMatMult(t(B), B)
B_pred = get_basis_matrix(x = x_pred, degree = 3, num_interior_knots = knot)


p = ncol(B) ; I_p = diag(1,p)

P = diff(diag(1, ncol(B)), diff = 2) ; P = as.matrix(P); P = t(P) %*% P
P = P + diag(1e-10, ncol(P))

Bty = eigenMapMatMult( t(B), y )
nu0_n = nu_0 + n
nu0_sigma20 = nu_0*sigma2_0


logMSE_mat = matrix(NA, length(g_vec), 3) ; logMSE_mat = data.frame(logMSE_mat)
names(logMSE_mat) = c('knot','g','logMSE')
logMSE_mat[,1] = rep(knot, each = length(g_vec))
logMSE_mat[,2] = rep(g_vec, 1)


for(j in 1:length(g_vec)){
  
  testpredictions = matrix(NA, nrow = nrow(testfunction), ncol = n_mcmc_sample) # HAVE TO DEFINE THIS INSIDE THE LOOP TO AVOID INDEXERROR FROM BURNIN!!!
  
  g = g_vec[j]
  
  BtB_invgP = BtB + 1/g * P  
  BtB_invgP_inv = armaInv(BtB_invgP) # have to be invertible
  mn = eigenMapMatMult(BtB_invgP_inv , Bty)
  SSR = yty - eigenMapMatMult(t(y),eigenMapMatMult(B, mn)) ; SSR = SSR[1,1]
  
  for(iter in 1:n_mcmc_sample){
    
    sigma2 = rinvgamma(n = 1, shape = nu0_n/2 , scale = (nu0_sigma20 + SSR)/2)
    beta = mvnfast::rmvn(n = 1, mu = mn , sigma = sigma2 * BtB_invgP_inv) #, checkSymmetry = F)
    
    if(iter > nburnin){
      testpredictions[,iter] = eigenMapMatMult(B_pred, t(beta))
    }
  }
  
  testpredictions = testpredictions[,(nburnin+1):n_mcmc_sample] # get rid of burn in
  model_averaging = rowMeans(testpredictions) # posterior mean
  logMSE_mat[length(g_vec)+j,3] = log(mean((model_averaging - test_Y_signal)^2))
  testpredictions_list[[j+1]] = data.frame(x = x_pred, y = model_averaging)
}


plot(x, y, xlab = "", ylab ="", cex = 0.4, main = paste("Cubic-Pspline: prediction, rank(B)=",ncol(B),', n=',n,', sigma=',0.1), cex.main = 1.25)
lines(x_pred, test_Y_signal, col = "black", lwd = 4)
lines(x_pred, testpredictions_list[[2]]$y, col = "red", lwd = 2)
lines(x_pred, testpredictions_list[[3]]$y, col = "green", lwd = 2)
lines(x_pred, testpredictions_list[[4]]$y, col = "blue", lwd = 2)
lines(x_pred, testpredictions_list[[5]]$y, col = "purple", lwd = 2)
legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)
dev.off()

















png('G:/My Drive/FunctionalDataAnalysis/Project/kernelridge/Pspline_eigenanalysis.png', width = 15, height = 8, units = "in", res = 300)

m1 = matrix(c(
  1,1,2,3,4,
  1,1,5,6,7,
  1,1,8,9,10), nrow = 3, ncol = 5, byrow = T
)
layout(m1)

plot(1:25, seq(-0.2,1.2, length = 25), type = "n", xaxt = "n", xlab = "order", cex.axis = 1.5, cex.lab = 1.5, ylab = "eigenvalues", main = paste("Cubic-Pspline:First 25 eigenvalues of smoother matrix, rank(B) =", ncol(B)), cex.main = 1.4)
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)

points(1:25,eval1[1:25], type = "o", col = "red", pch = 16)
points(1:25,eval2[1:25], type = "o", col = "green", pch = 16)
points(1:25,eval3[1:25], type = "o", col = "blue", pch = 16)
points(1:25,eval4[1:25], type = "o", col = "purple", pch = 16)
legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)
axis(1, xaxp = c(1,25,12), las = 2)



orders = 1:9

plot(x, eval1[orders[1]] * evec1[,orders[1]], col = "red", type = "l", xlab = "", ylab = "", main = "Pspline-1st eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, eval2[orders[1]] * evec2[,orders[1]] + rep(-0.002, n), col = "green")
lines(x, -eval3[orders[1]] * evec3[,orders[1]] + rep(0.002, n), col = "blue")
lines(x, eval4[orders[1]] * evec4[,orders[1]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)



plot(x, eval1[orders[2]] * evec1[,orders[2]], col = "red", type = "l", xlab = "",ylab = "", main = "Pspline-2nd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, eval2[orders[2]] * evec2[,orders[2]] + rep(-0.002, n), col = "green")
lines(x, eval3[orders[2]] * evec3[,orders[2]] + rep(0.002, n), col = "blue")
lines(x, eval4[orders[2]] * evec4[,orders[2]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)



plot(x, eval1[orders[3]] * evec1[,orders[3]], col = "red", type = "l", xlab = "",ylab = "", main = "Pspline-3rd eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x,-eval2[orders[3]] * -evec2[,orders[3]], col = "green")
lines(x, eval3[orders[3]] * evec3[,orders[3]] + rep(0.002, n), col = "blue")
lines(x, eval4[orders[3]] * evec4[,orders[3]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)


plot(x, eval1[orders[4]] * evec1[,orders[4]], col = "red", type = "l", xlab = "",ylab = "", main = "Pspline-4th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, eval2[orders[4]] * evec2[,orders[4]], col = "green")
lines(x, -eval3[orders[4]] * evec3[,orders[4]], col = "blue")
lines(x, eval4[orders[4]] * evec4[,orders[4]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)


plot(x, eval1[orders[5]] * evec1[,orders[5]], col = "red", type = "l", xlab = "",ylab = "", main = "Pspline-5th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, eval2[orders[5]] * evec2[,orders[5]], col = "green")
lines(x, -eval3[orders[5]] * evec3[,orders[5]], col = "blue")
lines(x, -eval4[orders[5]] * evec4[,orders[5]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)


plot(x, eval1[orders[6]] * evec1[,orders[6]], col = "red", type = "l", xlab = "",ylab = "", main = "Pspline-6th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x,-eval2[orders[6]] * -evec2[,orders[6]], col = "green")
lines(x, eval3[orders[6]] * evec3[,orders[6]], col = "blue")
lines(x, -eval4[orders[6]] * evec4[,orders[6]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)


plot(x, eval1[orders[7]] * evec1[,orders[7]], col = "red", type = "l", xlab = "",ylab = "", main = "Pspline-7th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -eval2[orders[7]] * evec2[,orders[7]], col = "green")
lines(x, -eval3[orders[7]] * evec3[,orders[7]], col = "blue")
lines(x, eval4[orders[7]] * evec4[,orders[7]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)


plot(x, eval1[orders[8]] * evec1[,orders[8]], col = "red", type = "l", xlab = "", ylab = "",main = "Pspline-8th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x, -eval2[orders[8]] * evec2[,orders[8]], col = "green")
lines(x, -eval3[orders[8]] * evec3[,orders[8]], col = "blue")
lines(x, eval4[orders[8]] * evec4[,orders[8]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)


plot(x, eval1[orders[9]] * evec1[,orders[9]], col = "red", type = "l", xlab = "", ylab = "",main = "Pspline-9th eigenvalue*eigenvector", cex.main = 1.2, ylim = c(-0.15, 0.15))
abline(h = 0, lty = 2)
lines(x,eval2[orders[9]] * -evec2[,orders[9]], col = "green")
lines(x, eval3[orders[9]] * evec3[,orders[9]], col = "blue")
lines(x, eval4[orders[9]] * evec4[,orders[9]], col = "purple")
# legend('topright',paste('lam = ', round(1/g_vec,4)), col = c('red','green','blue','purple'), lwd = 2)

dev.off()

